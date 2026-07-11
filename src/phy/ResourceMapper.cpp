#include "phy/PhyInterfaces.h"
#include "common/NrTables.h"

namespace nr {
namespace phy {

class ResourceMapperImpl : public IResourceMapper {
public:
    void map_pdsch(ResourceGrid& grid, const ComplexMat& precoded_symbols,
                   int rb_start, int rb_len, const SymbolAllocation& sym_alloc) override {
        int n_ant = grid.n_ant;
        int n_sc_per_rb = 12;
        int sym_start = sym_alloc.start;
        int sym_len = sym_alloc.length;
        
        int sym_idx = 0;
        for (int sym = sym_start; sym < sym_start + sym_len; sym++) {
            if (sym >= grid.n_symbols) break;
            
            for (int rb = rb_start; rb < rb_start + rb_len; rb++) {
                for (int re = 0; re < n_sc_per_rb; re++) {
                    int sc = rb * n_sc_per_rb + re;
                    if (sc >= grid.n_subcarriers) continue;
                    
                    if (sym_idx >= precoded_symbols.n_rows) return;
                    
                    for (int ant = 0; ant < n_ant; ant++) {
                        if (ant < precoded_symbols.n_cols) {
                            grid.set_re(ant, sym, sc, precoded_symbols(sym_idx, ant));
                        }
                    }
                    sym_idx++;
                }
            }
        }
    }
    
    ComplexMat extract_pdsch(const ResourceGrid& grid, int rb_start, int rb_len,
                             const SymbolAllocation& sym_alloc) override {
        int n_ant = grid.n_ant;
        int n_sc_per_rb = 12;
        int sym_start = sym_alloc.start;
        int sym_len = sym_alloc.length;
        
        int total_res = 0;
        for (int sym = sym_start; sym < sym_start + sym_len; sym++) {
            if (sym >= grid.n_symbols) break;
            for (int rb = rb_start; rb < rb_start + rb_len; rb++) {
                for (int re = 0; re < n_sc_per_rb; re++) {
                    int sc = rb * n_sc_per_rb + re;
                    if (sc < grid.n_subcarriers) {
                        total_res++;
                    }
                }
            }
        }
        
        ComplexMat extracted(total_res, n_ant, arma::fill::zeros);
        
        int sym_idx = 0;
        for (int sym = sym_start; sym < sym_start + sym_len; sym++) {
            if (sym >= grid.n_symbols) break;
            
            for (int rb = rb_start; rb < rb_start + rb_len; rb++) {
                for (int re = 0; re < n_sc_per_rb; re++) {
                    int sc = rb * n_sc_per_rb + re;
                    if (sc >= grid.n_subcarriers) continue;
                    
                    if (sym_idx >= total_res) continue;
                    
                    for (int ant = 0; ant < n_ant; ant++) {
                        extracted(sym_idx, ant) = grid.get_re(ant, sym, sc);
                    }
                    sym_idx++;
                }
            }
        }
        
        return extracted;
    }
};

std::unique_ptr<IResourceMapper> create_resource_mapper() {
    return std::make_unique<ResourceMapperImpl>();
}

} // namespace phy
} // namespace nr
