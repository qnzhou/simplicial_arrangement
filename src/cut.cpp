#include "BSPNode.h"
#include "SimplicialArrangement.h"
#include "cut.h"

#include <implicit_predicates/implicit_predicates.h>

#include <absl/container/flat_hash_map.h>

#include <algorithm>
#include <optional>

namespace {
using namespace simplicial_arrangement;
using namespace implicit_predicates;

template <int DIM>
using OptionalPoint = std::optional<Point<DIM>>;

template <int DIM>
using OptionalEdge = std::optional<Edge<DIM>>;

template <int DIM>
using OptionalFace = std::optional<Face<DIM>>;

template <int DIM>
using OptionalCell = std::optional<Cell<DIM>>;

using EdgeIndexMap = absl::flat_hash_map<std::array<size_t, 2>, size_t>;

/**
 * Cut a 3D face using a cut plane.
 *
 * @param[in]  arrangement  Arrangement context.
 * @param[in]  cut_plane_index
 * @param[in]  face  The face to be cut.
 *
 * TODO: Should cut edge be returned in this form (planes) or in index form?
 * @returns three optional values:
 *   * Cut edge: The intersection of face and cut plane, oriented ccw on
 *               the positive side of the cut plane.
 *   * Negative subface if any.
 *   * Positive subface if any.
 */
template <typename Scalar>
std::tuple<OptionalEdge<3>, OptionalFace<3>, OptionalFace<3>> cut_face(
    SimplicialArrangement<Scalar, 3>& arrangement,
    EdgeIndexMap& edge_map,
    size_t cut_plane_index,
    const Face<3>& face)
{
    const size_t num_edges = face.edge_planes.size();
    const auto& planes = arrangement.get_planes();
    const auto& cut_plane = planes[cut_plane_index];
    const auto& supporting_plane = planes[face.supporting_plane];

    std::vector<implicit_predicates::Orientation> orientations;
    orientations.reserve(num_edges);

    bool non_negative = true, non_positive = true;
    std::vector<size_t> vertices_on_cut_plane;
    vertices_on_cut_plane.reserve(2);

    for (size_t i = 0; i < num_edges; i++) {
        const size_t prev_edge = (i + num_edges - 1) % num_edges;
        orientations.push_back(
            implicit_predicates::orient3d(planes[face.edge_planes[prev_edge]].data(),
                supporting_plane.data(),
                planes[face.edge_planes[i]].data(),
                cut_plane.data()));
        assert(orientations.back() != implicit_predicates::INVALID);
        non_negative = non_negative && orientations.back() >= 0;
        non_positive = non_positive && orientations.back() <= 0;
        if (orientations.back() == implicit_predicates::ZERO) {
            vertices_on_cut_plane.push_back(i);
        }
    }

    auto add_vertex = [&](size_t e_prev, size_t e_curr, size_t e_next) {
        size_t v0 = arrangement.get_vertex_index({face.supporting_plane, e_prev, e_curr});
        size_t v1 = arrangement.get_vertex_index({face.supporting_plane, e_curr, e_next});
        assert(v0 != simplicial_arrangement::INVALID);
        assert(v1 != simplicial_arrangement::INVALID);

        if (v0 > v1) std::swap(v0, v1);
        std::array<size_t, 2> e({v0, v1});

        auto itr = edge_map.find(e);
        if (itr == edge_map.end()) {
            size_t id = arrangement.get_vertex_count();
            edge_map[e] = id;
            assert(edge_map.contains(e));
            arrangement.register_vertex({face.supporting_plane, e_curr, cut_plane_index}, id);
            arrangement.bump_vertex_count();
            logger().debug("Register ({}, {}, {}) as {} (new)",
                face.supporting_plane,
                e_curr,
                cut_plane_index,
                id);
        } else {
#ifndef NDEBUG
            // Just a sanity check.
            assert(arrangement.has_vertex({face.supporting_plane, e_curr, cut_plane_index}));
            assert(itr->second ==
                   arrangement.get_vertex_index({face.supporting_plane, e_curr, cut_plane_index}));
#endif
        }
    };

    auto add_touching_vertex = [&](size_t e0, size_t e1) {
        size_t id = arrangement.get_vertex_index({face.supporting_plane, e0, e1});
        assert(id != simplicial_arrangement::INVALID);
        arrangement.register_vertex({face.supporting_plane, e0, cut_plane_index}, id);
        arrangement.register_vertex({face.supporting_plane, e1, cut_plane_index}, id);
    };

    if (non_negative && vertices_on_cut_plane.size() < 2) {
        // Case 1: Face is on the positive side or tangent at a vertex.
        return {std::nullopt, std::nullopt, face};
    } else if (non_positive && vertices_on_cut_plane.size() < 2) {
        // Case 2: Face is on the negative side or tangent at a vertex.
        return {std::nullopt, face, std::nullopt};
    } else if (!non_positive && !non_negative) {
        // Case 3: Face is cut into 2 halves.
        Face<3> negative, positive;
        negative.supporting_plane = face.supporting_plane;
        positive.supporting_plane = face.supporting_plane;
        negative.edge_planes.reserve(num_edges + 1);
        positive.edge_planes.reserve(num_edges + 1);

        Edge<3> cut_edge(cut_plane_index,
            simplicial_arrangement::INVALID,
            face.supporting_plane,
            simplicial_arrangement::INVALID);

        for (size_t i = 0; i < num_edges; i++) {
            const size_t h = (i + num_edges - 1) % num_edges;
            const size_t j = (i + 1) % num_edges;
            const size_t prev_plane = face.edge_planes[h];
            const size_t curr_plane = face.edge_planes[i];
            const size_t next_plane = face.edge_planes[j];
            switch (orientations[i]) {
            case implicit_predicates::NEGATIVE:
                switch (orientations[j]) {
                case implicit_predicates::NEGATIVE:
                    // Negative - Negative: miss.
                    negative.edge_planes.push_back(curr_plane);
                    break;
                case implicit_predicates::POSITIVE:
                    // Negative - Positive: crossing.
                    negative.edge_planes.push_back(curr_plane);
                    negative.edge_planes.push_back(cut_plane_index);
                    positive.edge_planes.push_back(curr_plane);
                    cut_edge.next_plane = curr_plane;
                    add_vertex(prev_plane, curr_plane, next_plane);
                    break;
                case implicit_predicates::ZERO:
                    // Negative - Zero: touching.
                    negative.edge_planes.push_back(curr_plane);
                    negative.edge_planes.push_back(cut_plane_index);
                    cut_edge.next_plane = curr_plane;
                    add_touching_vertex(curr_plane, next_plane);
                    break;
                default: assert(false);
                }
                break;
            case implicit_predicates::POSITIVE:
                switch (orientations[j]) {
                case implicit_predicates::NEGATIVE:
                    // Positive - Negative: crossing.
                    positive.edge_planes.push_back(curr_plane);
                    positive.edge_planes.push_back(cut_plane_index);
                    negative.edge_planes.push_back(curr_plane);
                    cut_edge.prev_plane = curr_plane;
                    add_vertex(prev_plane, curr_plane, next_plane);
                    break;
                case implicit_predicates::POSITIVE:
                    // Positive - Positive: miss.
                    positive.edge_planes.push_back(curr_plane);
                    break;
                case implicit_predicates::ZERO:
                    // Positive - Zero: touching.
                    positive.edge_planes.push_back(curr_plane);
                    positive.edge_planes.push_back(cut_plane_index);
                    cut_edge.prev_plane = curr_plane;
                    add_touching_vertex(curr_plane, next_plane);
                    break;
                default: assert(false);
                }
                break;
            case implicit_predicates::ZERO:
                add_touching_vertex(prev_plane, curr_plane);
                switch (orientations[j]) {
                case implicit_predicates::NEGATIVE:
                    // Zero - Negative: touching.
                    negative.edge_planes.push_back(curr_plane);
                    break;
                case implicit_predicates::POSITIVE:
                    // Zero - Positive: touching.
                    positive.edge_planes.push_back(curr_plane);
                    break;
                case implicit_predicates::ZERO:
                    // Zero - Zero: impossible.
                    // This case should be captured by case 4 below, not here.
                    assert(false);
                    break;
                default: assert(false);
                }
                break;
            default: assert(false);
            }
        }

        assert(cut_edge.prev_plane != simplicial_arrangement::INVALID);
        assert(cut_edge.next_plane != simplicial_arrangement::INVALID);
        assert(cut_edge.prev_plane != cut_edge.next_plane);
        return {std::move(cut_edge), std::move(negative), std::move(positive)};
    } else if (!non_negative || !non_positive) {
        // Case 4: Face is tangent to the cut plane at an edge.
        assert(vertices_on_cut_plane.size() == 2);
        size_t e_curr = vertices_on_cut_plane[0];
        size_t e_next = vertices_on_cut_plane[1];
        if (e_curr == 0 && e_next == num_edges - 1) {
            std::swap(e_curr, e_next);
        }
        const size_t e_prev = (e_curr + num_edges - 1) % num_edges;
        assert(e_next == (e_curr + 1) % num_edges);
        assert(e_prev != e_next);

        // Basically, e_curr, supporting_plane and cut_plane intersect at a line.
        // Need to update vertex index map.
        size_t v0 = arrangement.get_vertex_index(
            {face.supporting_plane, face.edge_planes[e_prev], face.edge_planes[e_curr]});
        size_t v1 = arrangement.get_vertex_index(
            {face.supporting_plane, face.edge_planes[e_curr], face.edge_planes[e_next]});
        assert(v0 != simplicial_arrangement::INVALID);
        assert(v1 != simplicial_arrangement::INVALID);
        arrangement.register_vertex(
            {face.supporting_plane, face.edge_planes[e_prev], cut_plane_index}, v0);
        arrangement.register_vertex(
            {face.supporting_plane, face.edge_planes[e_next], cut_plane_index}, v1);

        if (non_negative) {
            return {Edge<3>({cut_plane_index,
                        face.edge_planes[e_prev],
                        face.supporting_plane,
                        face.edge_planes[e_next]}),
                std::nullopt,
                face};
        } else {
            assert(non_positive);
            return {Edge<3>({cut_plane_index,
                        face.edge_planes[e_next],
                        face.supporting_plane,
                        face.edge_planes[e_prev]}),
                face,
                std::nullopt};
        }
    } else {
        // Case 5: Face is coplanar with cut plane.
        assert(vertices_on_cut_plane.size() == num_edges);
        arrangement.register_coplanar_planes(face.supporting_plane, cut_plane_index);
        return {std::nullopt, std::nullopt, std::nullopt};
    }
}

template <typename Scalar>
std::optional<std::array<Cell<3>, 2>> cut_cell(SimplicialArrangement<Scalar, 3>& arrangement,
    EdgeIndexMap& edge_map,
    size_t cut_plane_index,
    const Cell<3>& cell)
{
    logger().debug("Cutting cell with cut plane {}", cut_plane_index);
    const size_t num_faces = cell.faces.size();
    Cell<3> negative, positive;
    negative.faces.reserve(num_faces + 1);
    positive.faces.reserve(num_faces + 1);
    std::vector<Edge<3>> cut_edges;
    cut_edges.reserve(num_faces);

    for (size_t i = 0; i < num_faces; i++) {
        auto r = cut_face(arrangement, edge_map, cut_plane_index, cell.faces[i]);
        if (std::get<0>(r)) {
            const auto& e = std::get<0>(r).value();
            logger().debug("supporting: {}  prev: {}  curr: {}  next: {}",
                e.supporting_plane,
                e.prev_plane,
                e.curr_plane,
                e.next_plane);
        }
        if (!std::get<1>(r) && !std::get<2>(r)) {
            // A face is coplanar with the cell, thus no cut is made.
            return std::nullopt;
        } else {
            if (std::get<1>(r)) {
                negative.faces.push_back(std::move(std::get<1>(r).value()));
            }
            if (std::get<2>(r)) {
                positive.faces.push_back(std::move(std::get<2>(r).value()));
                if (std::get<0>(r)) {
                    cut_edges.push_back(std::move(std::get<0>(r).value()));
                }
            }
        }
    }

    // Cut plane did not split the cell in two.
    if (negative.faces.empty() || positive.faces.empty()) {
        return std::nullopt;
    }

    // Chain cut edges into faces.
    const size_t num_cut_edges = cut_edges.size();
    absl::flat_hash_map<size_t, std::pair<size_t, size_t>> next_map;
    next_map.reserve(num_cut_edges);
    size_t start_vertex = simplicial_arrangement::INVALID;
    for (const auto& e : cut_edges) {
        const size_t v0 =
            arrangement.get_vertex_index({e.supporting_plane, e.curr_plane, e.prev_plane});
        const size_t v1 =
            arrangement.get_vertex_index({e.supporting_plane, e.curr_plane, e.next_plane});
        assert(v0 != simplicial_arrangement::INVALID);
        assert(v1 != simplicial_arrangement::INVALID);
        next_map.insert({v0, {v1, e.curr_plane}});
        if (start_vertex == simplicial_arrangement::INVALID) {
            start_vertex = v0;
        }
    }

    Face<3> cut_face;
    cut_face.supporting_plane = cut_plane_index;
    cut_face.edge_planes.reserve(num_cut_edges);
    std::vector<size_t> vertex_loop;
    vertex_loop.reserve(num_cut_edges);
    size_t curr_vertex = start_vertex;
    size_t curr_edge = 0;
    do {
        auto itr = next_map.find(curr_vertex);
        assert(itr != next_map.end());
        vertex_loop.push_back(curr_vertex);
        std::tie(curr_vertex, curr_edge) = itr->second;
        cut_face.edge_planes.push_back(curr_edge);
        assert(cut_face.edge_planes.size() <= num_cut_edges);
    } while (curr_vertex != start_vertex);
    assert(cut_face.edge_planes.size() == num_cut_edges);

    // Ensure all vertices around the face are registered.
    for (size_t i = 0; i < num_cut_edges; i++) {
        const size_t curr_edge = cut_face.edge_planes[i];
        const size_t prev_edge = cut_face.edge_planes[(i + num_cut_edges - 1) % num_cut_edges];
        arrangement.register_vertex(
            {curr_edge, prev_edge, cut_face.supporting_plane}, vertex_loop[i]);
    }

    Face<3> cut_face_reversed = cut_face;
    std::reverse(cut_face_reversed.edge_planes.begin(), cut_face_reversed.edge_planes.end());

    positive.faces.push_back(std::move(cut_face_reversed));
    negative.faces.push_back(std::move(cut_face));

    return std::array<Cell<3>, 2>({std::move(negative), std::move(positive)});
}

/**
 * Cut a 2D cell using the cut plane.
 *
 * @param[in]  arrangement      Arrangement object storing all planes.
 * @param[in]  cut_plane_index  Index of the cut plane.
 * @param[in]  cell             The cell to be cut.
 *
 * @returns A tuple of negative and positive sub-cells.  If a sub-cell is empty,
 * the cutting does not generate that cell.
 */
template <typename Scalar>
std::optional<std::array<Cell<2>, 2>> cut_cell(SimplicialArrangement<Scalar, 2>& arrangement,
    EdgeIndexMap& edge_map,
    size_t cut_plane_index,
    const Cell<2>& cell)
{
    std::vector<implicit_predicates::Orientation> orientations;
    const size_t num_edges = cell.edges.size();
    orientations.reserve(num_edges);
    const auto& planes = arrangement.get_planes();
    const auto& cut_plane = planes[cut_plane_index];

    for (size_t i = 0; i < num_edges; i++) {
        const size_t prev_i = (i + num_edges - 1) % num_edges;
        const auto& prev_plane = planes[cell.edges[prev_i]];
        const auto& curr_plane = planes[cell.edges[i]];
        orientations.push_back(
            implicit_predicates::orient2d(prev_plane.data(), curr_plane.data(), cut_plane.data()));
        assert(orientations.back() != implicit_predicates::INVALID);
    }

    bool non_negative = true, non_positive = true;
    for (size_t i = 0; i < num_edges; i++) {
        const size_t j = (i + 1) % num_edges;
        const auto oi = orientations[i];
        const auto oj = orientations[j];
        if (oi == implicit_predicates::ZERO && oj == implicit_predicates::ZERO) {
            // Cell contains an edge that is collinear with the cut plane.
            arrangement.register_coplanar_planes(cell.edges[i], cut_plane_index);
        }
        non_negative &= (oi == implicit_predicates::POSITIVE || oi == implicit_predicates::ZERO);
        non_positive &= (oi == implicit_predicates::NEGATIVE || oi == implicit_predicates::ZERO);
    }

    // Case 1: cell is on the positive side of or tangent to the cutting plane.
    if (non_negative) {
        return std::nullopt;
    }

    // Case 2: cell is on the negative side of or tangent to the cutting plane.
    if (non_positive) {
        return std::nullopt;
    }

    // Case 3: cell is cut into 2 halves.
    Cell<2> positive, negative;
    positive.edges.reserve(num_edges + 1);
    negative.edges.reserve(num_edges + 1);

    auto add_vertex = [&](size_t e_prev, size_t e_curr, size_t e_next) {
        size_t v0 = arrangement.get_vertex_index({e_prev, e_curr});
        size_t v1 = arrangement.get_vertex_index({e_curr, e_next});
        assert(v0 != simplicial_arrangement::INVALID);
        assert(v1 != simplicial_arrangement::INVALID);

        if (v0 > v1) std::swap(v0, v1);
        std::array<size_t, 2> e({v0, v1});

        auto itr = edge_map.find(e);
        if (itr == edge_map.end()) {
            size_t id = arrangement.get_vertex_count();
            edge_map[e] = id;
            assert(edge_map.contains(e));
            arrangement.register_vertex({e_curr, cut_plane_index}, id);
            arrangement.bump_vertex_count();
            logger().debug("Register ({}, {}) as {} (new)", e_curr, cut_plane_index, id);
        } else {
#ifndef NDEBUG
            // Just a sanity check.
            assert(arrangement.has_vertex({e_curr, cut_plane_index}));
            assert(itr->second == arrangement.get_vertex_index({e_curr, cut_plane_index}));
#endif
        }
    };

    auto add_touching_vertex = [&](size_t e0, size_t e1) {
        size_t id = arrangement.get_vertex_index({e0, e1});
        assert(id != simplicial_arrangement::INVALID);
        arrangement.register_vertex({e0, cut_plane_index}, id);
        arrangement.register_vertex({e1, cut_plane_index}, id);
    };

    for (size_t i = 0; i < num_edges; i++) {
        const size_t h = (i + num_edges - 1) % num_edges;
        const size_t j = (i + 1) % num_edges;

        switch (orientations[i]) {
        case implicit_predicates::NEGATIVE:
            switch (orientations[j]) {
            case implicit_predicates::NEGATIVE:
                // Negative - Negative: miss.
                negative.edges.push_back(cell.edges[i]);
                break;
            case implicit_predicates::POSITIVE:
                // Negative - Positive: crossing.
                negative.edges.push_back(cell.edges[i]);
                negative.edges.push_back(cut_plane_index);
                positive.edges.push_back(cell.edges[i]);
                add_vertex(cell.edges[h], cell.edges[i], cell.edges[j]);
                break;
            case implicit_predicates::ZERO:
                // Negative - Zero: touching.
                negative.edges.push_back(cell.edges[i]);
                negative.edges.push_back(cut_plane_index);
                add_touching_vertex(cell.edges[i], cell.edges[j]);
                break;
            default: assert(false);
            }
            break;
        case implicit_predicates::POSITIVE:
            switch (orientations[j]) {
            case implicit_predicates::NEGATIVE:
                // Positive - Negative: crossing.
                positive.edges.push_back(cell.edges[i]);
                positive.edges.push_back(cut_plane_index);
                negative.edges.push_back(cell.edges[i]);
                add_vertex(cell.edges[h], cell.edges[i], cell.edges[j]);
                break;
            case implicit_predicates::POSITIVE:
                // Positive - Positive : miss.
                positive.edges.push_back(cell.edges[i]);
                break;
            case implicit_predicates::ZERO:
                // Positive - Zero : touching.
                positive.edges.push_back(cell.edges[i]);
                positive.edges.push_back(cut_plane_index);
                add_touching_vertex(cell.edges[i], cell.edges[j]);
                break;
            default: assert(false);
            }
            break;
        case implicit_predicates::ZERO:
            add_touching_vertex(cell.edges[h], cell.edges[i]);
            switch (orientations[j]) {
            case implicit_predicates::NEGATIVE:
                // Zero - Negative: touching.
                negative.edges.push_back(cell.edges[i]);
                break;
            case implicit_predicates::POSITIVE:
                // Zero - positive: touching.
                positive.edges.push_back(cell.edges[i]);
                break;
            case implicit_predicates::ZERO:
                // Zero - Zero: impossible.
                // Cell is tangent to the cut plane, should be handled by case 1 or 2.
                assert(false);
                break;
            default: assert(false);
            }
            break;
        default: assert(false);
        }
    }
    return std::array<Cell<2>, 2>({std::move(negative), std::move(positive)});
}

template <typename Scalar, int DIM>
void cut(SimplicialArrangement<Scalar, DIM>& arrangement,
    EdgeIndexMap& edge_map,
    BSPNode<DIM>& root,
    size_t cut_plane_index)
{
    if (root.negative == nullptr && root.positive == nullptr) {
        auto r = cut_cell(arrangement, edge_map, cut_plane_index, root.cell);
        if (r) {
            root.separating_plane = cut_plane_index;
            root.negative = std::make_unique<BSPNode<DIM>>();
            root.positive = std::make_unique<BSPNode<DIM>>();
            root.negative->cell = std::move((*r)[0]);
            root.positive->cell = std::move((*r)[1]);
        }
    } else {
        assert(root.negative != nullptr);
        assert(root.positive != nullptr);
        cut(arrangement, edge_map, *root.negative, cut_plane_index);
        cut(arrangement, edge_map, *root.positive, cut_plane_index);
    }
}

} // namespace

namespace simplicial_arrangement::internal {

void cut(SimplicialArrangement<double, 2>& arrangement, BSPNode<2>& root, size_t cut_plane_index)
{
    EdgeIndexMap edge_map;
    edge_map.reserve(arrangement.get_num_planes());
    ::cut(arrangement, edge_map, root, cut_plane_index);
}

void cut(SimplicialArrangement<Int, 2>& arrangement, BSPNode<2>& root, size_t cut_plane_index)
{
    EdgeIndexMap edge_map;
    edge_map.reserve(arrangement.get_num_planes());
    ::cut(arrangement, edge_map, root, cut_plane_index);
}

void cut(SimplicialArrangement<double, 3>& arrangement, BSPNode<3>& root, size_t cut_plane_index)
{
    EdgeIndexMap edge_map;
    edge_map.reserve(arrangement.get_num_planes() * arrangement.get_num_planes());
    ::cut(arrangement, edge_map, root, cut_plane_index);
}

void cut(SimplicialArrangement<Int, 3>& arrangement, BSPNode<3>& root, size_t cut_plane_index)
{
    EdgeIndexMap edge_map;
    edge_map.reserve(arrangement.get_num_planes() * arrangement.get_num_planes());
    ::cut(arrangement, edge_map, root, cut_plane_index);
}

} // namespace simplicial_arrangement::internal
