#include "extract_arrangement.h"
#include "ARComplex.h"
#include "SimplicialArrangementBuilder.h"

#include <absl/container/flat_hash_map.h>

#include <algorithm>

namespace {
using namespace simplicial_arrangement;

/**
 * Determines if 2 coplanar planes are consistently oriented.  Note that we are
 * assuming p1 and p2 are coplanar, otherwise the result is incorrect.
 */
template <typename Scalar, int DIM>
bool is_plane_consistently_oriented(const Plane<Scalar, DIM>& p1, const Plane<Scalar, DIM>& p2)
{
    for (size_t i = 0; i < DIM + 1; i++) {
        const Scalar v1 = p1[i];
        const Scalar v2 = p2[i];
        if (v1 == 0 && v2 == 0) continue;
        return (v1 > 0 && v2 > 0) || (v1 < 0 && v2 < 0);
    }
    // plane p1 and p2 are consistently 0 over the cell.
    return true;
}

template <typename Scalar, int DIM>
void extract_unique_planes(Arrangement<DIM>& r, SimplicialArrangementBuilder<Scalar, DIM>& builder)
{
    auto [coplanar_planes, unique_plane_indices] = builder.extract_coplanar_planes();
    std::vector<bool> coplanar_plane_orientations(unique_plane_indices.size(), true);
    for (const auto& planes : coplanar_planes) {
        const size_t num_planes = planes.size();
        assert(num_planes > 0);
        coplanar_plane_orientations[planes[0]] = true;
        const auto& ref_plane = builder.get_plane(planes[0]);
        for (size_t i = 1; i < num_planes; i++) {
            coplanar_plane_orientations[planes[i]] = is_plane_consistently_oriented<Scalar, DIM>(
                ref_plane, builder.get_plane(planes[i]));
        }
    }

    r.unique_plane_indices = std::move(unique_plane_indices);
    r.unique_planes = std::move(coplanar_planes);
    r.unique_plane_orientations = std::move(coplanar_plane_orientations);
}

template <typename Scalar, int DIM>
void extract_plane_orientations(
    Arrangement<DIM>& r, SimplicialArrangementBuilder<Scalar, DIM>& builder)
{
    const size_t num_cells = r.cells.size();
    const size_t num_planes = builder.get_num_planes();
    if (num_cells == 0) return;

    std::vector<bool> visited(num_cells, false);
    std::vector<bool> consistent(num_planes, true);
    std::vector<size_t> Q;
    Q.reserve(num_cells);
    Q.push_back(0);
    size_t id = 0;
    visited[id] = true;
    r.cells.front().plane_orientations = std::vector<bool>(num_planes, true);

    while (id < Q.size()) {
        size_t cell_id = Q[id];
        id++;
        auto& cell = r.cells[cell_id];
        size_t num_faces = cell.faces.size();

        for (size_t i = 0; i < num_faces; i++) {
            const auto& f = r.faces[cell.faces[i]];
            const bool ori = cell.face_orientations[i];
            size_t next_cell_id = ori ? f.negative_cell : f.positive_cell;
            assert(next_cell_id != cell_id);
            if (next_cell_id == Arrangement<3>::None) continue;
            if (visited[next_cell_id]) continue;

            auto& next_cell = r.cells[next_cell_id];
            next_cell.plane_orientations = cell.plane_orientations; // Copied intentionally.
            for (auto pid : r.unique_planes[r.unique_plane_indices[f.supporting_plane]]) {
                next_cell.plane_orientations[pid] = !next_cell.plane_orientations[pid];
                consistent[pid] =
                    (cell.plane_orientations[f.supporting_plane] == cell.face_orientations[i]) ==
                    (r.unique_plane_orientations[pid] ==
                        r.unique_plane_orientations[f.supporting_plane]);
            }
            visited[next_cell_id] = true;
            Q.push_back(next_cell_id);
        }
    }

    for (auto& cell : r.cells) {
        for (size_t i = 0; i < num_planes; i++) {
            cell.plane_orientations[i] = (cell.plane_orientations[i] == consistent[i]);
        }
    }
}

template <typename T>
typename T::iterator cyclic_next(T& container, const typename T::iterator& itr)
{
    assert(itr != container.end());
    typename T::iterator next_itr = std::next(itr);
    if (next_itr == container.end()) {
        next_itr = container.begin();
    }
    return next_itr;
};

template <typename T>
typename T::iterator cyclic_prev(T& container, const typename T::iterator& itr)
{
    assert(itr != container.end());
    typename T::iterator prev_itr;
    if (itr == container.begin()) {
        prev_itr = std::prev(container.end());
    } else {
        prev_itr = std::prev(itr);
    }
    return prev_itr;
};


template <typename Scalar>
Arrangement<2> extract_arrangement_2D(SimplicialArrangementBuilder<Scalar, 2>& builder)
{
    Arrangement<2> r;
    const size_t num_vertices = builder.get_vertex_count();
    r.vertices = builder.extract_vertices();

    // (vi, vj) |-> edge_index
    absl::flat_hash_map<std::array<size_t, 2>, size_t> edge_map;
    edge_map.reserve(num_vertices * 3); // Just a guess.

    auto add_face = [&](size_t v0,
                        size_t v1,
                        size_t supporting_plane_index,
                        bool cell_on_positive_side_of_supporting_plane) -> size_t {
        assert(v0 != v1);

        // Compute canonical ordering of the face.
        bool reversed = v0 > v1;
        if (reversed) std::swap(v0, v1);

        auto itr = edge_map.find({v0, v1});
        if (itr == edge_map.end()) {
            size_t fid = r.faces.size();

            // Generate face entry.
            Arrangement<2>::Face out_face;
            out_face.supporting_plane = supporting_plane_index;
            // By construction, a cell is always on the left side of the face in
            // its original ordering.  The following converts the ordering such
            // that positive side of the supporting plane is on the right side
            // of the edge.
            if (cell_on_positive_side_of_supporting_plane) {
                if (reversed)
                    out_face.vertices = {v0, v1};
                else
                    out_face.vertices = {v1, v0};
            } else {
                if (reversed)
                    out_face.vertices = {v1, v0};
                else
                    out_face.vertices = {v0, v1};
            }
            r.faces.push_back(std::move(out_face));

            // Note: Always use canonical ordering of the face as the key.
            edge_map.insert({{v0, v1}, fid});

            return fid;
        } else {
            return itr->second;
        }
    };

    auto add_cell = [&](const Cell<2>& cell, const std::vector<bool>& plane_orientations) {
        size_t num_cell_edges = cell.edges.size();
        const size_t cell_id = r.cells.size();

        Arrangement<2>::Cell out_cell;
        out_cell.faces.reserve(num_cell_edges);
        out_cell.face_orientations.reserve(num_cell_edges);

        std::vector<size_t> vertex_indices;
        vertex_indices.reserve(num_cell_edges);
        for (size_t i = 0; i < num_cell_edges; i++) {
            const auto& curr_plane = cell.edges[i];
            const auto& prev_plane = cell.edges[(i + num_cell_edges - 1) % num_cell_edges];

            vertex_indices.push_back(builder.get_vertex_index({prev_plane, curr_plane}));
        }

        for (size_t i = 0; i < num_cell_edges; i++) {
            size_t j = (i + 1) % num_cell_edges;

            bool oriented = plane_orientations[cell.edges[i]];
            auto fid = add_face(vertex_indices[i], vertex_indices[j], cell.edges[i], oriented);

            out_cell.faces.push_back(fid);
            out_cell.face_orientations.push_back(oriented);

            auto& out_face = r.faces[fid];
            if (oriented) {
                out_face.positive_cell = cell_id;
            } else {
                out_face.negative_cell = cell_id;
            }
        }
        r.cells.push_back(std::move(out_cell));
    };

    // Traverse the BSP tree and gather cells and faces.
    std::vector<bool> bounding_plane_orientations(builder.get_num_planes(), true);
    std::function<void(const BSPNode<2>&)> depth_first_traversal;
    depth_first_traversal = [&](const BSPNode<2>& node) {
        if (node.positive != nullptr && node.negative != nullptr) {
            bounding_plane_orientations[node.separating_plane] = true;
            depth_first_traversal(*node.positive);
            bounding_plane_orientations[node.separating_plane] = false;
            depth_first_traversal(*node.negative);
        } else {
            assert(node.positive == nullptr);
            assert(node.negative == nullptr);
            add_cell(node.cell, bounding_plane_orientations);
        }
    };

    extract_unique_planes(r, builder);
    depth_first_traversal(builder.get_root());
    extract_plane_orientations(r, builder);

    return r;
}

template <typename Scalar>
Arrangement<3> extract_arrangement_3D(SimplicialArrangementBuilder<Scalar, 3>& builder)
{
    Arrangement<3> r;
    r.vertices = builder.extract_vertices();

    absl::flat_hash_map<std::vector<size_t>, size_t> face_map;
    face_map.reserve(r.vertices.size() * 3); // TODO: need more accurate guess.

    auto add_face = [&](const Face<3>& face,
                        bool cell_on_positive_side_of_supporting_plane) -> size_t {
        const size_t num_edges = face.edge_planes.size();
        assert(num_edges >= 3);

        // Step 1: convert boundary edge plane loop into vertex index loop.
        std::vector<size_t> vertex_indices;
        vertex_indices.reserve(num_edges);

        for (size_t i = 0; i < num_edges; i++) {
            size_t prev_i = (i + num_edges - 1) % num_edges;
            vertex_indices.push_back(builder.get_vertex_index(
                {face.edge_planes[i], face.edge_planes[prev_i], face.supporting_plane}));
            assert(vertex_indices.back() < builder.get_vertex_count());
        }

        // Step 2: reorder vertex indices such that a face always starts
        // from the vertex with the smallest index, and the next vertex
        // is the smaller of the two neighboring vertices.  Together, they
        // uniquely determine the starting point and ordering of the boundary
        // vertex loop.

        // cell_on_positive_side_of_face == true if the cell is on the positive
        // side of the face after reordering.
        bool cell_on_positive_side_of_face = false;
        {
            auto itr = std::min_element(vertex_indices.begin(), vertex_indices.end());
            auto next_itr = cyclic_next(vertex_indices, itr);
            auto prev_itr = cyclic_prev(vertex_indices, itr);

            if (*next_itr < *prev_itr) {
                std::rotate(vertex_indices.begin(), itr, vertex_indices.end());
                cell_on_positive_side_of_face = false;
            } else {
                std::rotate(vertex_indices.begin(), next_itr, vertex_indices.end());
                std::reverse(vertex_indices.begin(), vertex_indices.end());
                cell_on_positive_side_of_face = true;
            }
        }

        size_t fid = simplicial_arrangement::INVALID;
        auto itr = face_map.find(vertex_indices);
        if (itr == face_map.end()) {
            fid = face_map.size();
            face_map.insert({vertex_indices, fid});

            // Step 3: generate new face entry.
            Arrangement<3>::Face out_face;
            out_face.supporting_plane = face.supporting_plane;
            if (cell_on_positive_side_of_face == cell_on_positive_side_of_supporting_plane) {
                out_face.vertices = std::move(vertex_indices);
            } else {
                std::reverse(vertex_indices.begin(), vertex_indices.end());
                out_face.vertices = std::move(vertex_indices);
            }

            r.faces.push_back(std::move(out_face));
        } else {
            fid = itr->second;
        }
        return fid;
    };

    auto add_cell = [&](const Cell<3>& cell, const std::vector<bool>& plane_orientations) {
        const size_t cell_id = r.cells.size();
        const size_t num_faces = cell.faces.size();

        Arrangement<3>::Cell out_cell;
        out_cell.faces.reserve(num_faces);
        out_cell.face_orientations.reserve(num_faces);

        for (const auto& face : cell.faces) {
            bool on_positive_side = plane_orientations[face.supporting_plane];
            auto fid = add_face(face, on_positive_side);
            out_cell.faces.push_back(fid);
            out_cell.face_orientations.push_back(on_positive_side);

            if (on_positive_side) {
                r.faces[fid].positive_cell = cell_id;
            } else {
                r.faces[fid].negative_cell = cell_id;
            }
        }

        r.cells.push_back(std::move(out_cell));
    };

    std::vector<bool> bounding_plane_orientations(builder.get_num_planes(), true);
    std::function<void(const BSPNode<3>&)> depth_first_traversal;
    depth_first_traversal = [&](const BSPNode<3>& node) {
        if (node.positive != nullptr && node.negative != nullptr) {
            bounding_plane_orientations[node.separating_plane] = true;
            depth_first_traversal(*node.positive);
            bounding_plane_orientations[node.separating_plane] = false;
            depth_first_traversal(*node.negative);
        } else {
            assert(node.positive == nullptr);
            assert(node.negative == nullptr);
            add_cell(node.cell, bounding_plane_orientations);
        }
    };

    extract_unique_planes(r, builder);
    depth_first_traversal(builder.get_root());
    extract_plane_orientations(r, builder);

    return r;
}

} // namespace

namespace simplicial_arrangement {

Arrangement<2> extract_arrangement(SimplicialArrangementBuilder<double, 2>& builder)
{
    return ::extract_arrangement_2D(builder);
}

Arrangement<2> extract_arrangement(SimplicialArrangementBuilder<Int, 2>& builder)
{
    return ::extract_arrangement_2D(builder);
}

Arrangement<3> extract_arrangement(SimplicialArrangementBuilder<double, 3>& builder)
{
    return ::extract_arrangement_3D(builder);
}

Arrangement<3> extract_arrangement(SimplicialArrangementBuilder<Int, 3>& builder)
{
    return ::extract_arrangement_3D(builder);
}

Arrangement<2> extract_arrangement(ARComplex<2>&& ar_complex)
{
    Arrangement<2> ar;
    ar.vertices = std::move(ar_complex.vertices);

    auto& edges = ar_complex.edges;
    size_t num_edges = edges.size();
    ar.faces.resize(num_edges);

    for (size_t i = 0; i < num_edges; i++) {
        auto& ce = edges[i];
        auto& e = ar.faces[i];
        e.vertices = {ce.vertices[0], ce.vertices[1]};
        e.positive_cell = INVALID; // TODO
        e.negative_cell = INVALID; // TODO
    }

    auto& faces = ar_complex.faces;
    size_t num_faces = faces.size();
    ar.cells.resize(num_faces);

    for (size_t i = 0; i < num_faces; i++) {
        auto& cf = faces[i];
        auto& f = ar.cells[i];
        f.faces = std::move(cf.edges);
        f.plane_orientations = std::move(cf.signs);
    }

    return ar;
}

Arrangement<3> extract_arrangement(ARComplex<3>&& ar_complex)
{
    Arrangement<3> ar;
    ar.vertices = std::move(ar_complex.vertices);

    auto& edges = ar_complex.edges;
    auto& faces = ar_complex.faces;
    size_t num_faces = faces.size();
    ar.faces.resize(num_faces);

    for (size_t i = 0; i < num_faces; i++) {
        auto& cf = faces[i];
        auto& f = ar.faces[i];
        const size_t num_bd_edges = cf.edges.size();
        assert(num_bd_edges >= 3);
        f.vertices.reserve(num_bd_edges);

        for (size_t j = 0; j < num_bd_edges; j++) {
            auto& curr_e = edges[cf.edges[j]];
            auto& next_e = edges[cf.edges[(j + 1) % num_bd_edges]];
            if (curr_e.vertices[0] == next_e.vertices[0] ||
                curr_e.vertices[0] == next_e.vertices[1]) {
                f.vertices.push_back(curr_e.vertices[0]);
            } else {
                assert(curr_e.vertices[1] == next_e.vertices[0] ||
                       curr_e.vertices[1] == next_e.vertices[1]);
                f.vertices.push_back(curr_e.vertices[1]);
            }
        }

        f.positive_cell = INVALID; // TODO
        f.negative_cell = INVALID; // TODO
    }

    auto& cells = ar_complex.cells;
    size_t num_cells = cells.size();
    ar.cells.resize(num_cells);

    for (size_t i = 0; i < num_cells; i++) {
        auto& cc = cells[i];
        auto& c = ar.cells[i];
        c.faces = std::move(cc.faces);
        c.plane_orientations = std::move(cc.signs);
    }

    return ar;
}

} // namespace simplicial_arrangement
