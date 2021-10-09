#include "extract_arrangement.h"
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

    auto extract_unique_planes = [&]() {
        auto [coplanar_planes, unique_plane_indices] = builder.extract_coplanar_planes();
        size_t num_unique_planes = coplanar_planes.size();
        std::vector<std::vector<bool>> coplanar_plane_orientations;
        coplanar_plane_orientations.reserve(num_unique_planes);
        for (const auto& planes : coplanar_planes) {
            if (planes.size() == 1)
                coplanar_plane_orientations.push_back({true});
            else {
                size_t num_planes = planes.size();
                std::vector<bool> orientations(num_planes);
                for (size_t i = 0; i < num_planes; i++) {
                    orientations[i] = is_plane_consistently_oriented<Scalar, 2>(
                        builder.get_plane(planes[0]), builder.get_plane(planes[i]));
                }
                coplanar_plane_orientations.push_back(std::move(orientations));
            }
        }

        r.unique_plane_indices = std::move(unique_plane_indices);
        r.unique_planes = std::move(coplanar_planes);
        r.unique_plane_orientations = std::move(coplanar_plane_orientations);
    };

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

    extract_unique_planes();
    depth_first_traversal(builder.get_root());

    return r;
}

template <typename Scalar>
Arrangement<3> extract_arrangement_3D(SimplicialArrangementBuilder<Scalar, 3>& builder)
{
    Arrangement<3> r;
    r.vertices = builder.extract_vertices();

    absl::flat_hash_map<std::vector<size_t>, size_t> face_map;
    face_map.reserve(r.vertices.size() * 3); // TODO: need more accurate guess.

    auto extract_unique_planes = [&]() {
        auto [coplanar_planes, unique_plane_indices] = builder.extract_coplanar_planes();
        size_t num_unique_planes = coplanar_planes.size();
        std::vector<std::vector<bool>> coplanar_plane_orientations;
        coplanar_plane_orientations.reserve(num_unique_planes);
        for (const auto& planes : coplanar_planes) {
            if (planes.size() == 1)
                coplanar_plane_orientations.push_back({true});
            else {
                size_t num_planes = planes.size();
                std::vector<bool> orientations(num_planes);
                for (size_t i = 0; i < num_planes; i++) {
                    orientations[i] = is_plane_consistently_oriented<Scalar, 3>(
                        builder.get_plane(planes[0]), builder.get_plane(planes[i]));
                }
                coplanar_plane_orientations.push_back(std::move(orientations));
            }
        }

        r.unique_plane_indices = std::move(unique_plane_indices);
        r.unique_planes = std::move(coplanar_planes);
        r.unique_plane_orientations = std::move(coplanar_plane_orientations);
    };


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

    extract_unique_planes();
    depth_first_traversal(builder.get_root());

    return r;
}

} // namespace

namespace simplicial_arrangement::internal {

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

} // namespace simplicial_arrangement::internal
