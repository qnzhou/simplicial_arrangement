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

template <typename Scalar>
std::optional<std::tuple<Point<3>, Edge<3>, Edge<3>>> cut_edge(
    SimplicialArrangement<Scalar, 3>& arrangement, size_t cutting_plane, const Edge<3>& edge)
{
    // TODO
    return {};
}

template <typename Scalar>
std::optional<std::tuple<Edge<3>, Face<3>, Face<3>>> cut_face(
    SimplicialArrangement<Scalar, 3>& arrangement, size_t cutting_plane, const Face<3>& face)
{
    // TODO
    return {};
}

template <typename Scalar>
std::optional<std::tuple<Face<3>, Cell<3>, Cell<3>>> cut_cell(
    SimplicialArrangement<Scalar, 3>& arrangement, size_t cutting_plane, const Cell<3>& cell)
{
    // TODO
    return {};
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

} // namespace

namespace simplicial_arrangement::internal {

void cut(
    const SimplicialArrangement<double, 2>& arrangement, BSPNode<2>& root, size_t cut_plane_index)
{
    if (root.negative == nullptr && root.positive == nullptr) {
        auto r = cut_cell(arrangement, cut_plane_index, root.cell);
        if (r) {
            root.separating_plane = cut_plane_index;
            root.negative = std::make_unique<BSPNode<2>>();
            root.positive = std::make_unique<BSPNode<2>>();
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

void cut(const SimplicialArrangement<Int, 2>& arrangement, BSPNode<2>& root, size_t cut_plane_index)
{
    if (root.negative == nullptr && root.positive == nullptr) {
        auto r = cut_cell(arrangement, cut_plane_index, root.cell);
        if (r) {
            root.separating_plane = cut_plane_index;
            root.negative = std::make_unique<BSPNode<2>>();
            root.positive = std::make_unique<BSPNode<2>>();
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

void cut(const SimplicialArrangement<double, 3>& arrangement, BSPNode<3>& root, size_t cut_plane)
{
    // TODO
}

void cut(const SimplicialArrangement<Int, 3>& arrangement, BSPNode<3>& root, size_t cut_plane)
{
    // TODO
}

} // namespace simplicial_arrangement::internal
