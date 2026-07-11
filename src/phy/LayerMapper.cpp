#include "phy/PhyInterfaces.h"
#include "common/NrTables.h"
#include <cmath>

namespace nr {
namespace phy {

class LayerMapperImpl : public ILayerMapper {
public:
    ComplexMat map(const ComplexVec& symbols, int n_layers) override {
        int n_symbols = symbols.n_elem;
        int sym_per_layer = (n_symbols + n_layers - 1) / n_layers;
        
        ComplexMat layered(sym_per_layer, n_layers, arma::fill::zeros);
        
        for (int i = 0; i < n_symbols; i++) {
            int layer = i % n_layers;
            int idx = i / n_layers;
            layered(idx, layer) = symbols(i);
        }
        
        return layered;
    }
    
    ComplexVec demap(const ComplexMat& layered_symbols, int n_layers) override {
        int sym_per_layer = layered_symbols.n_rows;
        int total_symbols = sym_per_layer * n_layers;
        
        ComplexVec symbols(total_symbols);
        
        for (int l = 0; l < n_layers; l++) {
            for (int i = 0; i < sym_per_layer; i++) {
                symbols(i * n_layers + l) = layered_symbols(i, l);
            }
        }
        
        return symbols;
    }
};

std::unique_ptr<ILayerMapper> create_layer_mapper() {
    return std::make_unique<LayerMapperImpl>();
}

} // namespace phy
} // namespace nr
