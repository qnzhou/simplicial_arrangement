#include <simplicial_arrangement/SimplicialArrangement.h>
#include <simplicial_arrangement/extract_arrangement.h>

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
    // plane p1 and p2 are constently 0 over the cell.
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
Arrangement<2> extract_arrangement_2D(SimplicialArrangement<Scalar, 2>& arrangement)
{
    Arrangement<2> r;
    const size_t num_vertices = arrangement.get_vertex_count();
    r.vertices = arrangement.extract_vertices();

    // (vi, vj) |-> edge_index
    absl::flat_hash_map<std::array<size_t, 2>, size_t> edge_map;
    edge_map.reserve(num_vertices * 3); // Just a guess.

    size_t num_unique_planes;
    std::vector<size_t> unique_plane_indices;
    std::tie(unique_plane_indices, num_unique_planes) = arrangement.get_unique_plane_indices();
    std::vector<std::vector<size_t>> coplanar_planes = arrangement.extract_coplanar_planes();

    auto add_face = [&](size_t v0, size_t v1) -> std::tuple<size_t, bool> {
        assert(v0 != v1);
        bool oriented = v0 > v1;
        if (oriented) std::swap(v0, v1);
        // oriented means the cell is on the positive side of the face.
        // Typically, if the face is oriented to point outwards, the cell is on
        // the negative side.

        auto itr = edge_map.find({v0, v1});
        if (itr == edge_map.end()) {
            Arrangement<2>::Face out_face;
            out_face.vertices = {v0, v1};

            size_t fid = r.faces.size();
            edge_map.insert({{v0, v1}, fid});
            r.faces.push_back(std::move(out_face));
            return {fid, oriented};
        } else {
            return {itr->second, oriented};
        }
    };

    auto add_cell = [&](const Cell<2>& cell, const std::vector<bool>& plane_orientations) {
        size_t num_cell_edges = cell.edges.size();
        const size_t cell_id = r.cells.size();
        const auto& planes = arrangement.get_planes();

        Arrangement<2>::Cell out_cell;
        out_cell.faces.reserve(num_cell_edges);
        out_cell.face_orientations.reserve(num_cell_edges);

        std::vector<size_t> vertex_indices;
        vertex_indices.reserve(num_cell_edges);
        for (size_t i = 0; i < num_cell_edges; i++) {
            const auto& curr_plane = cell.edges[i];
            const auto& prev_plane = cell.edges[(i + num_cell_edges - 1) % num_cell_edges];

            vertex_indices.push_back(arrangement.get_vertex_index({prev_plane, curr_plane}));
        }

        for (size_t i = 0; i < num_cell_edges; i++) {
            size_t j = (i + 1) % num_cell_edges;

            auto [fid, oriented] = add_face(vertex_indices[i], vertex_indices[j]);
            // oriented == true means cell is on the positive side of the face.

            out_cell.faces.push_back(fid);
            out_cell.face_orientations.push_back(oriented);

            auto& out_face = r.faces[fid];
            if (oriented) {
                out_face.positive_cell = cell_id;
            } else {
                out_face.negative_cell = cell_id;
            }

            const size_t seed_supporting_plane = cell.edges[i];
            bool seed_plane_orientation = plane_orientations[seed_supporting_plane];
            // seed_plane_orientation == true means cell is on the positive side of
            // the plane.
            size_t unique_plane_id = unique_plane_indices[seed_supporting_plane];
            out_face.supporting_planes.reserve(coplanar_planes[unique_plane_id].size());

            out_face.supporting_planes.push_back(seed_supporting_plane);
            out_face.supporting_plane_orientations.push_back(seed_plane_orientation == oriented);

            for (auto pid : coplanar_planes[unique_plane_id]) {
                if (pid == seed_supporting_plane) continue;
                out_face.supporting_planes.push_back(pid);
                if (is_plane_consistently_oriented<Scalar, 2>(
                        planes[seed_supporting_plane], planes[pid])) {
                    out_face.supporting_plane_orientations.push_back(
                        out_face.supporting_plane_orientations[0]);
                } else {
                    out_face.supporting_plane_orientations.push_back(
                        !out_face.supporting_plane_orientations[0]);
                }
            }
        }
        r.cells.push_back(std::move(out_cell));
    };

    // Traverse the BSP tree and gather cells and faces.
    std::vector<bool> bounding_plane_orientations(arrangement.get_num_planes(), true);
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
    depth_first_traversal(arrangement.get_root());

    return r;
}

template <typename Scalar>
Arrangement<3> extract_arrangement_3D(SimplicialArrangement<Scalar, 3>& arrangement)
{
    Arrangement<3> r;
    r.vertices = arrangement.extract_vertices();
    const auto& planes = arrangement.get_planes();

    absl::flat_hash_map<std::vector<size_t>, size_t> face_map;
    face_map.reserve(r.vertices.size() * 3); // TODO: need more accuate guess.

    size_t num_unique_planes;
    std::vector<size_t> unique_plane_indices;
    std::tie(unique_plane_indices, num_unique_planes) = arrangement.get_unique_plane_indices();
    std::vector<std::vector<size_t>> coplanar_planes = arrangement.extract_coplanar_planes();

    auto add_face = [&](const Face<3>& face,
                        const std::vector<bool>& plane_orientations) -> std::tuple<size_t, bool> {
        const size_t num_edges = face.edge_planes.size();
        assert(num_edges >= 3);

        // Step 1: convert boundary edge plane loop into vertex index loop.
        std::vector<size_t> vertex_indices;
        vertex_indices.reserve(num_edges);

        for (size_t i = 0; i < num_edges; i++) {
            size_t prev_i = (i + num_edges - 1) % num_edges;
            vertex_indices.push_back(arrangement.get_vertex_index(
                {face.edge_planes[i], face.edge_planes[prev_i], face.supporting_plane}));
            assert(vertex_indices.back() < arrangement.get_vertex_count());
        }

        // Step 2: reorder vertex indices such that a face alwasys starts
        // from the vertex with the smallest index, and the next vertex
        // is the smaller of the two neighboring vertices.  Together, they
        // uniquely determine the starting point and ordering of the boundary
        // vertex loop.

        // on_positive_side == true if the cell is on the positive side of the
        // face after reordering.
        bool on_positive_side = false;
        {
            auto itr = std::min_element(vertex_indices.begin(), vertex_indices.end());
            auto next_itr = cyclic_next(vertex_indices, itr);
            auto prev_itr = cyclic_prev(vertex_indices, itr);

            if (*next_itr < *prev_itr) {
                std::rotate(vertex_indices.begin(), itr, vertex_indices.end());
                on_positive_side = false;
            } else {
                std::rotate(vertex_indices.begin(), next_itr, vertex_indices.end());
                std::reverse(vertex_indices.begin(), vertex_indices.end());
                on_positive_side = true;
            }
        }

        size_t fid = simplicial_arrangement::INVALID;
        auto itr = face_map.find(vertex_indices);
        if (itr == face_map.end()) {
            fid = face_map.size();
            face_map.insert({vertex_indices, fid});

            // Step 3: generate new face entry.
            Arrangement<3>::Face out_face;
            out_face.vertices = std::move(vertex_indices);

            const size_t seed_plane_index = face.supporting_plane;
            const size_t unique_plane_id = unique_plane_indices[seed_plane_index];
            const size_t num_coplanar_planes = coplanar_planes[unique_plane_id].size();

            out_face.supporting_planes.reserve(num_coplanar_planes);
            out_face.supporting_plane_orientations.reserve(num_coplanar_planes);

            out_face.supporting_planes.push_back(seed_plane_index);
            out_face.supporting_plane_orientations.push_back(
                plane_orientations[seed_plane_index] == on_positive_side);

            for (auto plane_index : coplanar_planes[unique_plane_id]) {
                if (plane_index == seed_plane_index) continue;
                out_face.supporting_planes.push_back(plane_index);
                if (is_plane_consistently_oriented<Scalar, 3>(
                        planes[seed_plane_index], planes[plane_index])) {
                    out_face.supporting_plane_orientations.push_back(
                        out_face.supporting_plane_orientations[0]);
                } else {
                    out_face.supporting_plane_orientations.push_back(
                        !out_face.supporting_plane_orientations[0]);
                }
            }
            r.faces.push_back(std::move(out_face));
        } else {
            fid = itr->second;
        }
        return {fid, on_positive_side};
    };

    auto add_cell = [&](const Cell<3>& cell, const std::vector<bool>& plane_orientations) {
        const size_t cell_id = r.cells.size();
        const size_t num_faces = cell.faces.size();

        Arrangement<3>::Cell out_cell;
        out_cell.faces.reserve(num_faces);
        out_cell.face_orientations.reserve(num_faces);

        for (const auto& face : cell.faces) {
            auto [fid, on_positive_side] = add_face(face, plane_orientations);
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

    std::vector<bool> bounding_plane_orientations(arrangement.get_num_planes(), true);
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

    depth_first_traversal(arrangement.get_root());

    return r;
}

} // namespace

namespace simplicial_arrangement::internal {

Arrangement<2> extract_arrangement(SimplicialArrangement<double, 2>& arrangement)
{
    return ::extract_arrangement_2D(arrangement);
}

Arrangement<2> extract_arrangement(SimplicialArrangement<Int, 2>& arrangement)
{
    return ::extract_arrangement_2D(arrangement);
}

Arrangement<3> extract_arrangement(SimplicialArrangement<double, 3>& arrangement)
{
    return ::extract_arrangement_3D(arrangement);
}

Arrangement<3> extract_arrangement(SimplicialArrangement<Int, 3>& arrangement)
{
    return ::extract_arrangement_3D(arrangement);
}

} // namespace simplicial_arrangement::internal
