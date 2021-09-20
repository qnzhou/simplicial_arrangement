#include <simplicial_arrangement/BSPNode.h>
#include <simplicial_arrangement/SimplicialArrangement.h>
#include <simplicial_arrangement/cut.h>

#include <implicit_predicates/implicit_predicates.h>

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

// template <typename Scalar>
// std::optional<std::tuple<Point<3>, Edge<3>, Edge<3>>> cut_edge(
//    SimplicialArrangement<Scalar, 3>& arrangement, size_t cutting_plane, const Edge<3>& edge)
//{
//    // TODO
//    return {};
//}

/**
 * Cut a 3D face using a cut plane.
 *
 * @param[in]  arrangement  Arrangement context.
 * @param[in]  cut_plane_index
 * @param[in]  face  The face to be cut.
 *
 * @returns three optional values:
 *   * Cut edge: The intersection of face and cut plane, oriented ccw on
 *               the positive side of the cut plane.
 *   * Negative subface if any.
 *   * Positive subface if any.
 */
template <typename Scalar>
std::tuple<OptionalEdge<3>, OptionalFace<3>, OptionalFace<3>> cut_face(
    const SimplicialArrangement<Scalar, 3>& arrangement,
    size_t cut_plane_index,
    const Face<3>& face)
{
    const size_t num_edges = face.edge_planes.size();
    const auto& planes = arrangement.get_planes();
    const auto& cut_plane = planes[cut_plane_index];
    const auto& supporting_plane = planes[face.supporting_plane];

    std::vector<implicit_predicates::Orientation> orientations;
    orientations.reserve(num_edges);

    bool all_positive = true, all_negative = true;
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
        all_positive = all_positive && orientations.back() >= 0;
        all_negative = all_negative && orientations.back() <= 0;
        if (orientations.back() == implicit_predicates::ZERO) {
            vertices_on_cut_plane.push_back(i);
        }
    }

    if (all_positive && vertices_on_cut_plane.size() < 2) {
        // Case 1: Face is on the positive side or tangent at a vertex.
        return {std::nullopt, std::nullopt, face};
    } else if (all_negative && vertices_on_cut_plane.size() < 2) {
        // Case 2: Face is on the negative side or tangent at a vertex.
        return {std::nullopt, face, std::nullopt};
    } else if (!all_negative && !all_positive) {
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
            const size_t j = (i + 1) % num_edges;
            const size_t curr_plane = face.edge_planes[i];
            switch (orientations[i]) {
            case implicit_predicates::NEGATIVE:
                switch (orientations[j]) {
                case implicit_predicates::NEGATIVE:
                    negative.edge_planes.push_back(curr_plane);
                    break;
                case implicit_predicates::POSITIVE:
                    negative.edge_planes.push_back(curr_plane);
                    negative.edge_planes.push_back(cut_plane_index);
                    positive.edge_planes.push_back(curr_plane);
                    cut_edge.next_plane = curr_plane;
                    break;
                case implicit_predicates::ZERO:
                    negative.edge_planes.push_back(curr_plane);
                    negative.edge_planes.push_back(cut_plane_index);
                    cut_edge.next_plane = curr_plane;
                    break;
                default: assert(false);
                }
                break;
            case implicit_predicates::POSITIVE:
                switch (orientations[j]) {
                case implicit_predicates::NEGATIVE:
                    positive.edge_planes.push_back(curr_plane);
                    positive.edge_planes.push_back(cut_plane_index);
                    negative.edge_planes.push_back(curr_plane);
                    cut_edge.prev_plane = curr_plane;
                    break;
                case implicit_predicates::POSITIVE:
                    positive.edge_planes.push_back(curr_plane);
                    break;
                case implicit_predicates::ZERO:
                    positive.edge_planes.push_back(curr_plane);
                    positive.edge_planes.push_back(cut_plane_index);
                    cut_edge.prev_plane = curr_plane;
                    break;
                default: assert(false);
                }
                break;
            case implicit_predicates::ZERO:
                switch (orientations[j]) {
                case implicit_predicates::NEGATIVE:
                    negative.edge_planes.push_back(curr_plane);
                    break;
                case implicit_predicates::POSITIVE:
                    positive.edge_planes.push_back(curr_plane);
                    break;
                case implicit_predicates::ZERO: assert(false); break;
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
    } else if (!all_positive || !all_negative) {
        // Case 4: Face is tangent to the cut plane at an edge.
        assert(vertices_on_cut_plane.size() == 2);
        size_t e_curr = vertices_on_cut_plane[0];
        size_t e_next = vertices_on_cut_plane[1];
        if (e_curr == 0 && e_next == num_edges-1) {
            std::swap(e_curr, e_next);
        }
        const size_t e_prev = (e_curr + num_edges - 1) % num_edges;
        assert(e_next == (e_curr + 1) % num_edges);
        assert(e_prev != e_next);
        if (all_positive) {
            return {Edge<3>({cut_plane_index,
                        face.edge_planes[e_prev],
                        face.supporting_plane,
                        face.edge_planes[e_next]}),
                std::nullopt,
                face};
        } else {
            assert(all_negative);
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
        return {std::nullopt, std::nullopt, std::nullopt};
    }
}

template <typename Scalar>
std::optional<std::array<Cell<3>, 2>> cut_cell(const SimplicialArrangement<Scalar, 3>& arrangement,
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
        auto r = cut_face(arrangement, cut_plane_index, cell.faces[i]);
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

    const auto& planes = arrangement.get_planes();
    auto is_same_point = [&](const Point<3>& p, const Point<3>& q) {
        // TODO: could be faster if Scalar is Int.
        using namespace implicit_predicates;
        return (orient3d(planes[p[0]].data(),
                    planes[p[1]].data(),
                    planes[p[2]].data(),
                    planes[q[0]].data()) == ZERO &&
                orient3d(planes[p[0]].data(),
                    planes[p[1]].data(),
                    planes[p[2]].data(),
                    planes[q[1]].data()) == ZERO &&
                orient3d(planes[p[0]].data(),
                    planes[p[1]].data(),
                    planes[p[2]].data(),
                    planes[q[2]].data()) == ZERO);
    };

    // Chain cut edges into faces.
    // TODO: This is a O(n^2) algorithm.
    const size_t num_cut_edges = cut_edges.size();
    assert(num_cut_edges >= 3);
    std::vector<bool> visited(num_cut_edges, false);
    Face<3> cut_face;
    cut_face.supporting_plane = cut_plane_index;
    cut_face.edge_planes.reserve(num_cut_edges);
    size_t curr_edge = 0;
    do {
        const Point<3> curr_end{cut_edges[curr_edge].supporting_plane,
            cut_edges[curr_edge].curr_plane,
            cut_edges[curr_edge].next_plane};
        Point<3> next_start;
        for (size_t i = 0; i < num_cut_edges; i++) {
            if (visited[i]) continue;
            next_start = {
                cut_edges[i].supporting_plane, cut_edges[i].curr_plane, cut_edges[i].prev_plane};
            if (is_same_point(curr_end, next_start)) {
                visited[i] = true;
                curr_edge = i;
                cut_face.edge_planes.push_back(cut_edges[i].curr_plane);
                break;
            }
        }
    } while (curr_edge != 0);
    assert(cut_face.edge_planes.size() == num_cut_edges);

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
std::optional<std::array<Cell<2>, 2>> cut_cell(const SimplicialArrangement<Scalar, 2>& arrangement,
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

    // Case 1: cell is on the positive side of or tangent to the cutting plane.
    bool all_positive = std::all_of(
        orientations.begin(), orientations.end(), [](implicit_predicates::Orientation o) {
            return o == implicit_predicates::POSITIVE || o == implicit_predicates::ZERO;
        });
    if (all_positive) {
        return std::nullopt;
    }

    // Case 2: cell is on the negative side of or tangent to the cutting plane.
    bool all_negative = std::all_of(
        orientations.begin(), orientations.end(), [](implicit_predicates::Orientation o) {
            return o == implicit_predicates::NEGATIVE || o == implicit_predicates::ZERO;
        });
    if (all_negative) {
        return std::nullopt;
    }

    // Case 3: cell is cut into 2 halves.
    Cell<2> positive, negative;
    positive.edges.reserve(num_edges + 1);
    negative.edges.reserve(num_edges + 1);

    for (size_t i = 0; i < num_edges; i++) {
        const size_t j = (i + 1) % num_edges;

        switch (orientations[i]) {
        case implicit_predicates::NEGATIVE:
            switch (orientations[j]) {
            case implicit_predicates::NEGATIVE: negative.edges.push_back(cell.edges[i]); break;
            case implicit_predicates::POSITIVE:
                negative.edges.push_back(cell.edges[i]);
                negative.edges.push_back(cut_plane_index);
                positive.edges.push_back(cell.edges[i]);
                break;
            case implicit_predicates::ZERO:
                negative.edges.push_back(cell.edges[i]);
                negative.edges.push_back(cut_plane_index);
                break;
            default: assert(false);
            }
            break;
        case implicit_predicates::POSITIVE:
            switch (orientations[j]) {
            case implicit_predicates::NEGATIVE:
                positive.edges.push_back(cell.edges[i]);
                positive.edges.push_back(cut_plane_index);
                negative.edges.push_back(cell.edges[i]);
                break;
            case implicit_predicates::POSITIVE: positive.edges.push_back(cell.edges[i]); break;
            case implicit_predicates::ZERO:
                positive.edges.push_back(cell.edges[i]);
                positive.edges.push_back(cut_plane_index);
                break;
            default: assert(false);
            }
            break;
        case implicit_predicates::ZERO:
            switch (orientations[j]) {
            case implicit_predicates::NEGATIVE: negative.edges.push_back(cell.edges[i]); break;
            case implicit_predicates::POSITIVE: positive.edges.push_back(cell.edges[i]); break;
            case implicit_predicates::ZERO:
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
void cut(const SimplicialArrangement<Scalar, DIM>& arrangement,
    BSPNode<DIM>& root,
    size_t cut_plane_index)
{
    if (root.negative == nullptr && root.positive == nullptr) {
        auto r = cut_cell(arrangement, cut_plane_index, root.cell);
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
        cut(arrangement, *root.negative, cut_plane_index);
        cut(arrangement, *root.positive, cut_plane_index);
    }
}

} // namespace

namespace simplicial_arrangement::internal {

void cut(
    const SimplicialArrangement<double, 2>& arrangement, BSPNode<2>& root, size_t cut_plane_index)
{
    ::cut(arrangement, root, cut_plane_index);
}

void cut(const SimplicialArrangement<Int, 2>& arrangement, BSPNode<2>& root, size_t cut_plane_index)
{
    ::cut(arrangement, root, cut_plane_index);
}

void cut(
    const SimplicialArrangement<double, 3>& arrangement, BSPNode<3>& root, size_t cut_plane_index)
{
    ::cut(arrangement, root, cut_plane_index);
}

void cut(const SimplicialArrangement<Int, 3>& arrangement, BSPNode<3>& root, size_t cut_plane_index)
{
    ::cut(arrangement, root, cut_plane_index);
}

} // namespace simplicial_arrangement::internal
