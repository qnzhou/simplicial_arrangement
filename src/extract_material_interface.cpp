#include "extract_material_interface.h"

#include "MIComplex.h"

#include <absl/container/flat_hash_set.h>

namespace simplicial_arrangement {

MaterialInterface<2> extract_material_interface(MIComplex<2>&& mi_complex)
{
    MaterialInterface<2> mi;
    mi.vertices = std::move(mi_complex.vertices);

    auto& edges = mi_complex.edges;
    size_t num_edges = edges.size();
    mi.faces.resize(num_edges);

    for (size_t i=0; i<num_edges; i++) {
        auto& ce = edges[i];
        auto& e = mi.faces[i];
        e.vertices = {ce.vertices[0], ce.vertices[1]};
        e.positive_material_label = ce.positive_material_label;
        e.negative_material_label = ce.negative_material_label;
    }

    auto& faces = mi_complex.faces;
    size_t num_faces = faces.size();
    mi.cells.resize(num_faces);

    for (size_t i=0; i<num_faces; i++) {
        auto& cf = faces[i];
        auto& f = mi.cells[i];
        f.faces = std::move(cf.edges);
        f.material_label = cf.material_label;
    }

    return mi;
}

MaterialInterface<3> extract_material_interface(MIComplex<3>&& mi_complex)
{
    MaterialInterface<3> mi;
    mi.vertices = std::move(mi_complex.vertices);

    auto& edges = mi_complex.edges;
    auto& faces = mi_complex.faces;
    size_t num_faces = faces.size();
    mi.faces.resize(num_faces);

    absl::flat_hash_set<size_t> vertex_set;
    vertex_set.reserve(mi.vertices.size());
    for (size_t i=0; i<num_faces; i++) {
        auto& cf = faces[i];
        auto& f = mi.faces[i];
        const size_t num_bd_edges = cf.edges.size();
        assert(num_bd_edges >= 3);
        f.vertices.reserve(num_bd_edges);

        vertex_set.clear();
        for (size_t j=0; j<num_bd_edges; j++) {
            auto& e = edges[cf.edges[j]];
            vertex_set.insert(e.vertices[0]);
            vertex_set.insert(e.vertices[1]);
        }

        for (size_t vid : vertex_set) {
            f.vertices.push_back(vid);
        }

        f.positive_material_label = cf.positive_material_label;
        f.negative_material_label = cf.negative_material_label;
    }

    auto& cells = mi_complex.cells;
    size_t num_cells = cells.size();
    mi.cells.resize(num_cells);

    for (size_t i=0; i<num_cells; i++) {
        auto& cc = cells[i];
        auto& c = mi.cells[i];
        c.faces = std::move(cc.faces);
        c.material_label = cc.material_label;
    }

    return mi;
}

}
