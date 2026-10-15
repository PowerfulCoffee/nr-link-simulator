#!/usr/bin/env python3
"""
NR PDSCH link-level simulation web console.

Stdlib-only HTTP server (WSL/Linux). Serves the SPA in this directory and
exposes a REST API to save/load scenarios, launch the pdsch_bler_web runner,
stream its log, and serve result curves from run history.

API:
  GET  /api/modules                     -> modules.json registry
  GET  /api/scenarios                   -> list saved scenarios
  GET  /api/scenarios/<name>            -> load scenario JSON
  POST /api/scenarios                   -> save scenario JSON {name, ...}
  POST /api/run                         -> start run {scenario} (409 if busy)
  GET  /api/runs                        -> history list (meta[])
  GET  /api/runs/<id>/meta              -> run meta (status/progress)
  GET  /api/runs/<id>/log?offset=N      -> incremental log {offset,text,progress}
  GET  /api/runs/<id>/curves            -> curves JSON from CSV
  GET  /api/runs/<id>/scenario          -> scenario used by the run
  POST /api/runs/<id>/stop              -> terminate running run
"""
import json
import os
import re
import subprocess
import threading
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import urlparse, parse_qs

BASE = os.path.dirname(os.path.abspath(__file__))          # web/
ROOT = os.path.dirname(BASE)                               # repo root
RUNNER = os.environ.get(
    "RUNNER_BIN", os.path.join(ROOT, "build_wsl", "examples", "pdsch_bler_web"))
SCEN_DIR = os.path.join(BASE, "scenarios")
RUNS_DIR = os.path.join(BASE, "runs")

for d in (SCEN_DIR, RUNS_DIR):
    os.makedirs(d, exist_ok=True)

with open(os.path.join(BASE, "modules.json"), encoding="utf-8") as f:
    MODULES = json.load(f)

STATE = {"active_id": None, "lock": threading.Lock()}
PROGRESS_RE = re.compile(r"\[PROGRESS\]\s+(\d+)/(\d+)\s+(.*)")


def slot_seconds(scs_hz: float) -> float:
    """Slot duration from SCS: 1 ms at 15 kHz, halved per doubling."""
    mu = 0
    f = scs_hz / 15000.0
    while f >= 2.0 and mu < 4:
        f /= 2.0
        mu += 1
    return 0.001 / (2 ** mu)


def scenario_to_args(params: dict, out_csv: str) -> list:
    """Map a flat scenario params dict onto pdsch_bler_web CLI arguments."""
    args = []
    def add(key, val):
        if val is None or val == "":
            return
        if isinstance(val, list):
            val = ",".join(str(v) for v in val)
        args.extend([key, str(val)])

    g = params
    # aliases: modules.json keys -> runner CLI keys
    g = dict(params)
    if "snr_stop" in g and "snr_end" not in g:
        g["snr_end"] = g.pop("snr_stop")
    if "max_errors" in g and "target_errors" not in g:
        g["target_errors"] = g.pop("max_errors")
    if "n_ldpc_iterations" in g and "ldpc_iters" not in g:
        g["ldpc_iters"] = g.pop("n_ldpc_iterations")
    add("--mcs", g.get("mcs_list"))
    # modules.json select values -> runner numeric mcs_table (Types.h: 1=qam64, 2=qam256)
    mt = g.get("mcs_table")
    if isinstance(mt, str):
        mt = {"qam64": 1, "qam256": 2}.get(mt.lower(), None)
    add("--mcs_table", mt)
    add("--n_rb", g.get("n_rb"))
    add("--n_tx", g.get("n_tx"))
    add("--n_rx", g.get("n_rx"))
    add("--n_layers", g.get("n_layers"))
    if g.get("scs"):
        add("--scs", int(float(g["scs"]) * 1000))
    add("--csi", g.get("csi"))
    add("--channel", g.get("channel"))
    add("--delay_spread_ns", g.get("delay_spread_ns"))
    add("--doppler", g.get("doppler"))
    add("--los", g.get("los"))
    add("--max_blocks", g.get("max_blocks"))
    add("--target_errors", g.get("target_errors"))
    add("--ldpc_iters", g.get("ldpc_iters"))
    add("--seed", g.get("seed"))
    add("--dmrs_type", g.get("dmrs_type"))
    add("--dmrs_addpos", g.get("dmrs_addpos"))
    add("--dmrs_duration", g.get("dmrs_duration"))
    add("--rnti", g.get("n_rnti"))
    add("--n_id", g.get("n_id"))
    add("--n_id_dmrs", g.get("n_id_dmrs"))
    add("--precoder", g.get("precoder"))
    add("--csi_n1", g.get("csi_n1"))
    add("--csi_n2", g.get("csi_n2"))
    add("--csi_i1_phi", g.get("csi_i1_phi"))
    add("--csi_i11", g.get("csi_i11"))
    add("--csi_i12", g.get("csi_i12"))
    add("--csi_i2", g.get("csi_i2"))
    add("--snr_start", g.get("snr_start"))
    add("--snr_end", g.get("snr_end"))
    add("--snr_step", g.get("snr_step"))
    args.extend(["--out", out_csv])
    return args


class RunManager:
    """Owns at most one active simulation process; persists run history."""

    def __init__(self):
        self.lock = threading.Lock()
        self.active = None  # dict: id, proc, dir, scenario

    def start(self, scenario: dict) -> dict:
        with self.lock:
            if self.active and self.active["proc"].poll() is None:
                raise RuntimeError("a simulation is already running")
            rid = time.strftime("%Y%m%d_%H%M%S")
            rdir = os.path.join(RUNS_DIR, rid)
            os.makedirs(rdir, exist_ok=True)
            out_csv = os.path.join(rdir, "result.csv")
            args = scenario_to_args(scenario.get("params", {}), out_csv)
            # pdsch_bler_web prints [PROGRESS] lines per SNR point to stdout
            logf = open(os.path.join(rdir, "log.txt"), "ab", buffering=0)
            proc = subprocess.Popen(
                [RUNNER] + args, stdout=logf, stderr=subprocess.STDOUT,
                cwd=ROOT, start_new_session=True)
            meta = {
                "id": rid, "name": scenario.get("name", rid),
                "scenario": scenario, "csv": out_csv,
                "status": "running", "pid": proc.pid,
                "start_time": time.time(), "cmd": [RUNNER] + args,
            }
            with open(os.path.join(rdir, "meta.json"), "w", encoding="utf-8") as f:
                json.dump(meta, f, indent=2, ensure_ascii=False)
            with open(os.path.join(rdir, "scenario.json"), "w", encoding="utf-8") as f:
                json.dump(scenario, f, indent=2, ensure_ascii=False)
            self.active = {"id": rid, "proc": proc, "dir": rdir, "logf": logf}
            threading.Thread(target=self._reaper, args=(proc, rdir), daemon=True).start()
            return meta

    def _reaper(self, proc, rdir):
        code = proc.wait()
        try:
            self.active["logf"].close()
        except Exception:
            pass
        mpath = os.path.join(rdir, "meta.json")
        try:
            with open(mpath, encoding="utf-8") as f:
                meta = json.load(f)
            meta["status"] = "finished" if code == 0 else "failed"
            meta["exit_code"] = code
            meta["end_time"] = time.time()
            with open(mpath, "w", encoding="utf-8") as f:
                json.dump(meta, f, indent=2, ensure_ascii=False)
        except Exception:
            pass
        with self.lock:
            if self.active and self.active["proc"] is proc:
                self.active = None

    def stop(self, rid: str) -> bool:
        with self.lock:
            if not self.active or self.active["id"] != rid:
                return False
            self.active["proc"].terminate()
            return True

    @property
    def active_id(self):
        return self.active["id"] if self.active else None


MGR = RunManager()


def load_meta(rid: str) -> dict:
    p = os.path.join(RUNS_DIR, rid, "meta.json")
    if not os.path.isfile(p):
        return None
    with open(p, encoding="utf-8") as f:
        return json.load(f)


def list_runs() -> list:
    out = []
    if not os.path.isdir(RUNS_DIR):
        return out
    for rid in sorted(os.listdir(RUNS_DIR), reverse=True):
        m = load_meta(rid)
        if m:
            out.append({k: m.get(k) for k in
                        ("id", "name", "status", "start_time", "end_time")})
    return out


def read_log(rid: str, offset: int) -> dict:
    """Return log bytes from offset plus parsed progress, if any."""
    path = os.path.join(RUNS_DIR, rid, "log.txt")
    text = ""
    if os.path.isfile(path):
        with open(path, "rb") as f:
            f.seek(max(offset, 0))
            text = f.read().decode("utf-8", errors="replace")
            new_offset = f.tell()
    else:
        new_offset = offset
    prog = None
    for line in reversed(text.splitlines()):
        mo = PROGRESS_RE.search(line)
        if mo:
            prog = {"done": int(mo.group(1)), "total": int(mo.group(2)),
                    "label": mo.group(3)}
            break
    return {"offset": new_offset, "text": text, "progress": prog}


def parse_curves(rid: str) -> dict:
    """Convert result.csv -> frontend contract:
    {mcs: [{mcs, qm, rate, tbs, points: [{snr, bler, tput, blocks}]}], tbs: {...}}
    tput unit: Mbps  (tbs bits per 1 ms slot -> tbs/1000 Mbps).
    """
    path = os.path.join(RUNS_DIR, rid, "result.csv")
    empty = {"mcs": [], "tbs": {}}
    if not os.path.isfile(path):
        return empty
    by_mcs = {}
    tbs_map = {}
    meta_map = {}
    with open(path, encoding="utf-8") as f:
        header = f.readline().strip().split(",")
        idx = {name: i for i, name in enumerate(header)}
        for line in f:
            parts = line.strip().split(",")
            if len(parts) < len(header):
                continue
            try:
                mcs = int(float(parts[idx["mcs"]]))
                snr = float(parts[idx["snr_db"] if "snr_db" in idx else idx["snr"]])
                bler = float(parts[idx["bler"]])
            except (KeyError, ValueError, IndexError):
                continue
            tbs = int(float(parts[idx["tbs"]])) if "tbs" in idx else 0
            qm = int(float(parts[idx["qm"]])) if "qm" in idx else 0
            rate = float(parts[idx["rate"]]) if "rate" in idx else 0.0
            bcol = "n_blocks" if "n_blocks" in idx else ("blocks" if "blocks" in idx else None)
            blocks = int(float(parts[idx[bcol]])) if bcol else 0
            if tbs:
                tbs_map[mcs] = tbs
                meta_map[mcs] = (qm, rate, tbs)
            by_mcs.setdefault(mcs, []).append(
                {"snr": snr, "bler": bler, "blocks": blocks,
                 "tput": tbs * (1.0 - bler) / 1000.0})
    mcs_list = []
    for mcs in sorted(by_mcs):
        qm, rate, tbs = meta_map.get(mcs, (0, 0.0, tbs_map.get(mcs, 0)))
        mcs_list.append({"mcs": mcs, "qm": qm, "rate": rate, "tbs": tbs,
                         "points": sorted(by_mcs[mcs], key=lambda p: p["snr"])})
    return {"mcs": mcs_list, "tbs": tbs_map}


class Handler(BaseHTTPRequestHandler):
    protocol_version = "HTTP/1.1"

    def log_message(self, fmt, *args):  # quiet access log
        pass

    # ---- helpers ----
    def _json(self, obj, code=200):
        body = json.dumps(obj, ensure_ascii=False).encode("utf-8")
        self.send_response(code)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def _read_body(self):
        n = int(self.headers.get("Content-Length", 0) or 0)
        raw = self.rfile.read(n) if n else b"{}"
        return json.loads(raw.decode("utf-8"))

    def _file(self, path, ctype):
        if not os.path.isfile(path):
            self._json({"error": "not found"}, 404)
            return
        with open(path, "rb") as f:
            data = f.read()
        self.send_response(200)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    # ---- routing ----
    def do_GET(self):
        u = urlparse(self.path)
        p = u.path
        q = parse_qs(u.query)
        if p == "/api/modules":
            return self._json(MODULES)
        if p == "/api/scenarios":
            names = sorted(x[:-5] for x in os.listdir(SCEN_DIR)
                           if x.endswith(".json")) if os.path.isdir(SCEN_DIR) else []
            return self._json(names)
        if p.startswith("/api/scenarios/"):
            name = p.rsplit("/", 1)[1]
            fp = os.path.join(SCEN_DIR, name + ".json")
            if not os.path.isfile(fp):
                return self._json({"error": "not found"}, 404)
            with open(fp, encoding="utf-8") as f:
                return self._json(json.load(f))
        if p == "/api/runs":
            return self._json(list_runs())
        if p == "/api/status":
            m = load_meta(MGR.active_id) if MGR.active_id else None
            return self._json({"active": MGR.active_id,
                               "meta": m})
        m = re.match(r"^/api/runs/([^/]+)/(\w+)$", p)
        if m:
            rid, what = m.group(1), m.group(2)
            if what == "meta":
                meta = load_meta(rid)
                return self._json(meta if meta else {"error": "not found"}, 200 if meta else 404)
            if what == "log":
                offset = int(q.get("offset", [0])[0])
                if not load_meta(rid):
                    return self._json({"error": "not found"}, 404)
                return self._json(read_log(rid, offset))
            if what == "curves":
                if not load_meta(rid):
                    return self._json({"error": "not found"}, 404)
                return self._json(parse_curves(rid))
            if what == "scenario":
                fp = os.path.join(RUNS_DIR, rid, "scenario.json")
                if not os.path.isfile(fp):
                    return self._json({"error": "not found"}, 404)
                with open(fp, encoding="utf-8") as f:
                    return self._json(json.load(f))
        # static SPA
        if p in ("/", "/index.html"):
            return self._file(os.path.join(BASE, "index.html"),
                              "text/html; charset=utf-8")
        fp = os.path.normpath(os.path.join(BASE, p.lstrip("/")))
        if fp.startswith(BASE) and os.path.isfile(fp):
            ctype = "application/json; charset=utf-8" if fp.endswith(".json") else \
                    "text/css; charset=utf-8" if fp.endswith(".css") else \
                    "application/javascript; charset=utf-8" if fp.endswith(".js") else \
                    "text/plain; charset=utf-8"
            return self._file(fp, ctype)
        self._json({"error": "not found"}, 404)

    def do_POST(self):
        u = urlparse(self.path)
        p = u.path
        if p == "/api/scenarios":
            body = self._read_body()
            name = re.sub(r"[^A-Za-z0-9_\-]", "_", body.get("name", "scenario"))
            fp = os.path.join(SCEN_DIR, name + ".json")
            with open(fp, "w", encoding="utf-8") as f:
                json.dump(body, f, indent=2, ensure_ascii=False)
            return self._json({"ok": True, "name": name})
        if p == "/api/run":
            body = self._read_body()
            scenario = body.get("scenario", body)
            try:
                meta = MGR.start(scenario)
                return self._json(meta, 201)
            except RuntimeError as e:
                return self._json({"error": str(e)}, 409)
            except FileNotFoundError as e:
                return self._json({"error": f"runner missing: {e}"}, 500)
        m = re.match(r"^/api/runs/([^/]+)/stop$", p)
        if m:
            ok = MGR.stop(m.group(1))
            return self._json({"stopped": ok})
        self._json({"error": "not found"}, 404)


def main():
    import sys
    port = int(sys.argv[1]) if len(sys.argv) > 1 else 8777
    if not os.path.isfile(RUNNER):
        print(f"[WARN] runner not found: {RUNNER} (build it first)", flush=True)
    srv = ThreadingHTTPServer(("0.0.0.0", port), Handler)
    print(f"[INFO] NR PDSCH link-level web console: http://localhost:{port}", flush=True)
    print(f"[INFO] runner: {RUNNER}", flush=True)
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    main()

