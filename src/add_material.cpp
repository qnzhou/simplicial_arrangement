#include "add_material.h"
#include "MaterialInterfaceBuilder.h"

#include <implicit_predicates/implicit_predicates.h>

#include <optional>

namespace {
using EdgeIndexMap = absl::flat_hash_map<std::array<size_t, 2>, size_t>;
using namespace simplicial_arrangement;

template <typename Scalar>
implicit_predicates::Orientation compute_orientation(
    const MaterialInterfaceBuilder<Scalar, 2>& builder,
    const Joint<2>& p,
    const Material<Scalar, 2>& material)
{
    short vertex_type = 0;
    if (p[0] > 2) vertex_type |= 1;
    if (p[1] > 2) vertex_type |= 2;
    if (p[2] > 2) vertex_type |= 4;

    auto get_corner_id = [](size_t i, size_t j) -> size_t {
        assert(i <= 2 && j <= 2);
        if (i != 0 && j != 0) return 0;
        if (i != 1 && j != 1) return 1;
        if (i != 2 && j != 2) return 2;
        logger().error("Invalid simplex boundary inputs: {}, {}", i, j);
        return simplicial_arrangement::INVALID;
    };

    auto compute_orientation_0d =
        [&](size_t i, size_t j, const Material<Scalar, 2>& m) -> implicit_predicates::Orientation {
        const size_t corner_id = get_corner_id(i, j);
        if (m[corner_id] > material[corner_id])
            return implicit_predicates::POSITIVE;
        else if (m[corner_id] < material[corner_id])
            return implicit_predicates::NEGATIVE;
        else
            return implicit_predicates::ZERO;
    };

    auto compute_orientation_1d =
        [&](size_t i,
            const Material<Scalar, 2>& m0,
            const Material<Scalar, 2>& m1) -> implicit_predicates::Orientation {
        assert(i <= 2);
        switch (i) {
        case 0: {
            const Scalar mm0[]{m0[1], m0[2]};
            const Scalar mm1[]{m1[1], m1[2]};
            const Scalar mm[]{material[1], material[2]};
            return implicit_predicates::mi_orient1d(mm0, mm1, mm);
        }
        case 1: {
            const Scalar mm0[]{m0[0], m0[2]};
            const Scalar mm1[]{m1[0], m1[2]};
            const Scalar mm[]{material[0], material[2]};
            return implicit_predicates::mi_orient1d(mm0, mm1, mm);
        }
        case 2: {
            const Scalar mm0[]{m0[0], m0[1]};
            const Scalar mm1[]{m1[0], m1[1]};
            const Scalar mm[]{material[0], material[1]};
            return implicit_predicates::mi_orient1d(mm0, mm1, mm);
        }
        default:
            logger().error("Invalid simplex boundary: {}", i);
            return implicit_predicates::INVALID;
        }
    };

    switch (vertex_type) {
    case 1: {
        const auto& m0 = builder.get_material(p[0]);
        return compute_orientation_0d(p[1], p[2], m0);
    }
    case 2: {
        const auto& m1 = builder.get_material(p[1]);
        return compute_orientation_0d(p[0], p[2], m1);
    }
    case 3: {
        const auto& m0 = builder.get_material(p[0]);
        const auto& m1 = builder.get_material(p[1]);
        return compute_orientation_1d(p[2], m0, m1);
    }
    case 4: {
        const auto& m2 = builder.get_material(p[2]);
        return compute_orientation_0d(p[0], p[1], m2);
    }
    case 5: {
        const auto& m0 = builder.get_material(p[0]);
        const auto& m2 = builder.get_material(p[2]);
        return compute_orientation_1d(p[1], m0, m2);
    }
    case 6: {
        const auto& m1 = builder.get_material(p[1]);
        const auto& m2 = builder.get_material(p[2]);
        return compute_orientation_1d(p[0], m1, m2);
    }
    case 7: {
        const auto& m0 = builder.get_material(p[0]);
        const auto& m1 = builder.get_material(p[1]);
        const auto& m2 = builder.get_material(p[2]);
        return implicit_predicates::mi_orient2d(m0.data(), m1.data(), m2.data(), material.data());
    }
    default:
        logger().error("Impossible vertex type case detected: {}", vertex_type);
        return implicit_predicates::INVALID;
    }
}

/**
 * Remove unused verices and faces.
 */
void consolidate(MaterialInterface<2>& material_interface)
{
    auto& vertices = material_interface.vertices;
    auto& faces = material_interface.faces;
    auto& cells = material_interface.cells;

    const size_t num_vertices = vertices.size();
    const size_t num_faces = faces.size();

    std::vector<bool> active_vertex(num_vertices, false);
    std::vector<bool> active_face(num_faces, false);

    size_t vertex_count = 0;
    size_t face_count = 0;

    for (const auto& c : cells) {
        for (auto fid : c.faces) {
            if (!active_face[fid]) {
                active_face[fid] = true;
                face_count++;
            } else {
                continue;
            }

            const auto& f = faces[fid];
            for (auto vid : f.vertices) {
                if (!active_vertex[vid]) {
                    active_vertex[vid] = true;
                    vertex_count++;
                }
            }
        }
    }

    size_t active_vertex_count = 0;
    std::vector<size_t> vertex_map(num_vertices, INVALID);
    for (size_t i = 0; i < num_vertices; i++) {
        if (!active_vertex[i]) continue;

        vertex_map[i] = active_vertex_count;
        std::swap(vertices[active_vertex_count], vertices[i]);
        active_vertex_count++;
    }
    vertices.resize(vertex_count);

    size_t active_face_count = 0;
    std::vector<size_t> face_map(num_faces, INVALID);
    for (size_t i = 0; i < num_faces; i++) {
        if (!active_face[i]) continue;

        face_map[i] = active_face_count;
        std::swap(faces[active_face_count], faces[i]);
        active_face_count++;
    }
    faces.resize(face_count);

    for (auto& f : faces) {
        std::transform(f.vertices.begin(), f.vertices.end(), f.vertices.begin(), [&](size_t i) {
            return vertex_map[i];
        });
    }
    for (auto& c : cells) {
        std::transform(
            c.faces.begin(), c.faces.end(), c.faces.begin(), [&](size_t i) { return face_map[i]; });
    }
}

template <typename Scalar>
void add_material(MaterialInterfaceBuilder<Scalar, 2>& builder, size_t material_index)
{
    using Face = MaterialInterface<2>::Face;
    using Cell = MaterialInterface<2>::Cell;

    const auto& material = builder.get_material(material_index);
    auto& material_interface = builder.get_material_interface();
    const size_t num_vertices = material_interface.vertices.size();
    const size_t num_faces = material_interface.faces.size();
    const size_t num_cells = material_interface.cells.size();
    logger().info("adding material {}", material_index);
    logger().debug("Before: {} {} {}", num_vertices, num_faces, num_cells);

    auto& vertices = material_interface.vertices;
    auto& faces = material_interface.faces;
    auto& cells = material_interface.cells;

    // Reserve capacity.
    vertices.reserve(num_vertices + num_faces);
    faces.reserve(num_faces * 3);
    cells.reserve(num_cells * 3);

    // Step 1: handle 0-faces.
    std::vector<implicit_predicates::Orientation> orientations;
    orientations.reserve(num_vertices);
    for (size_t i = 0; i < num_vertices; i++) {
        const auto& p = vertices[i];
        orientations.push_back(compute_orientation(builder, p, material));
        assert(orientations.back() != implicit_predicates::INVALID);
        logger().debug("vertex {}: {}, {}, {}", i, p[0], p[1], p[2]);
        logger().debug("orientation: {}", orientations.back());
    }

    // Step 2: handle 1-faces.
    std::vector<size_t> intersection_points(num_faces, INVALID);
    std::vector<size_t> positive_face_map(num_faces, INVALID);
    std::vector<size_t> negative_face_map(num_faces, INVALID);

    auto set_subface_no_cut = [&](size_t i, bool positive) {
        if (positive) {
            positive_face_map[i] = i;
        } else {
            negative_face_map[i] = i;
        }
    };

    auto set_subfaces = [&](size_t i, Face positive_subface, Face negative_subface) {
        faces.push_back(std::move(positive_subface));
        positive_face_map[i] = faces.size() - 1;

        faces.push_back(std::move(negative_subface));
        negative_face_map[i] = faces.size() - 1;
    };

    auto cut_face = [&](size_t i) {
        const auto& f = faces[i];
        const auto& end_points = f.vertices;
        const auto o0 = orientations[end_points[0]];
        const auto o1 = orientations[end_points[1]];

        if (o0 == 0) intersection_points[i] = end_points[0];
        if (o1 == 0) intersection_points[i] = end_points[1];

        if (o0 == 0 && o1 == 0) {
            // set_subface_no_cut(i, true);
            // set_subface_no_cut(i, false);
        } else if (o0 >= 0 && o1 >= 0) {
            set_subface_no_cut(i, true);
        } else if (o0 <= 0 && o1 <= 0) {
            set_subface_no_cut(i, false);
        } else {
            size_t m0 = f.positive_material_label;
            size_t m1 = f.negative_material_label;
            vertices.push_back({m0, m1, material_index});
            intersection_points[i] = vertices.size() - 1;
            Face positive_subface, negative_subface;

            if (o0 > 0 && o1 < 0) {
                positive_subface.vertices = {end_points[0], intersection_points[i]};
                negative_subface.vertices = {intersection_points[i], end_points[1]};
            } else {
                assert(o0 < 0);
                assert(o1 > 0);
                negative_subface.vertices = {end_points[0], intersection_points[i]};
                positive_subface.vertices = {intersection_points[i], end_points[1]};
            }
            positive_subface.positive_material_label = m0;
            positive_subface.negative_material_label = m1;
            negative_subface.positive_material_label = m0;
            negative_subface.negative_material_label = m1;
            set_subfaces(i, std::move(positive_subface), std::move(negative_subface));
        }
        logger().debug("face {}: {} {}",
            i,
            (positive_face_map[i] == INVALID ? "INVALID" : fmt::format("{}", positive_face_map[i])),
            (negative_face_map[i] == INVALID ? "INVALID"
                                             : fmt::format("{}", negative_face_map[i])));
    };

    for (size_t i = 0; i < num_faces; i++) {
        cut_face(i);
    }

    // Step 3: handle 2-faces.
    auto get_ordered_edge_end_points = [&](size_t cell_id,
                                           size_t edge_id) -> std::array<size_t, 2> {
        const auto& c = cells[cell_id];
        const auto& f = faces[edge_id];
        if (f.negative_material_label == c.material_label) {
            return {f.vertices[0], f.vertices[1]};
        } else {
            assert(f.positive_material_label == c.material_label);
            return {f.vertices[1], f.vertices[0]};
        }
    };

    auto cut_cell = [&](size_t i) -> std::array<std::optional<Cell>, 2> {
        const auto& c = cells[i];
        const size_t num_bd_faces = c.faces.size();
        logger().debug("cell {} has {} bd faces", i, num_bd_faces);

        // Find edges that straddle the material boundary.
        size_t first_positive_idx = INVALID;
        size_t first_negative_idx = INVALID;
        std::array<size_t, 2> cut_points{INVALID, INVALID};

        for (size_t j = 0; j < num_bd_faces; j++) {
            const size_t fid = c.faces[j];
            size_t intersection_point = intersection_points[fid];
            logger().debug("{}: intersection point {}", fid, intersection_point);

            if (intersection_point == INVALID) continue;

            size_t positive_subface = positive_face_map[fid];
            size_t negative_subface = negative_face_map[fid];
            size_t end_vertex = get_ordered_edge_end_points(i, fid)[1];
            auto end_orientation = orientations[end_vertex];

            if (positive_subface != INVALID && end_orientation > 0) {
                assert(first_positive_idx == INVALID);
                first_positive_idx = j;
                cut_points[0] = intersection_point;
            }
            if (negative_subface != INVALID && end_orientation < 0) {
                assert(first_negative_idx == INVALID);
                first_negative_idx = j;
                cut_points[1] = intersection_point;
            }
        }
        logger().debug("first positive face: {}", first_positive_idx);
        logger().debug("first negative face: {}", first_negative_idx);
        logger().debug("cut point: {}, {}", cut_points[0], cut_points[1]);

        if (first_positive_idx == INVALID || first_negative_idx == INVALID) {
            // No cut.
            logger().debug("No cut!");
            bool on_positive_side = false;
            for (size_t j = 0; j < num_bd_faces; j++) {
                const size_t fid = c.faces[j];
                const auto& f = faces[fid];
                if (orientations[f.vertices[0]] > 0) {
                    on_positive_side = true;
                    break;
                } else if (orientations[f.vertices[1]] < 0) {
                    on_positive_side = false;
                    break;
                }
            }
            if (on_positive_side) {
                return {c, std::nullopt};
            } else {
                return {std::nullopt, c};
            }
        }

        Face cut_face;
        cut_face.vertices = {cut_points[0], cut_points[1]};
        cut_face.positive_material_label = c.material_label;
        cut_face.negative_material_label = material_index;
        size_t cut_face_id = faces.size();
        faces.push_back(std::move(cut_face));

        Cell positive_cell, negative_cell;
        positive_cell.material_label = c.material_label;
        negative_cell.material_label = material_index;
        positive_cell.faces.reserve(num_bd_faces + 1);
        negative_cell.faces.reserve(num_bd_faces + 1);

        positive_cell.faces.push_back(cut_face_id);
        for (size_t j = 0; j < num_bd_faces; j++) {
            const size_t fid = c.faces[(j + first_positive_idx) % num_bd_faces];
            size_t positive_subface = positive_face_map[fid];
            if (positive_subface == INVALID) break;
            positive_cell.faces.push_back(positive_subface);
        }

        negative_cell.faces.push_back(cut_face_id);
        for (size_t j = 0; j < num_bd_faces; j++) {
            const size_t fid = c.faces[(j + first_negative_idx) % num_bd_faces];
            size_t negative_subface = negative_face_map[fid];
            if (negative_subface == INVALID) break;
            negative_cell.faces.push_back(negative_subface);

            // Update cell id.
            auto& f = faces[negative_subface];
            if (f.positive_material_label == c.material_label) {
                f.positive_material_label = material_index;
            } else {
                assert(f.negative_material_label == c.material_label);
                f.negative_material_label = material_index;
            }
        }

        assert(positive_cell.faces.size() > 2);
        assert(negative_cell.faces.size() > 2);
        return {positive_cell, negative_cell};
    };

    std::vector<Cell> positive_cells, negative_cells;
    positive_cells.reserve(num_cells + 1); // +1 for combined cell later.
    negative_cells.reserve(num_cells);
    for (size_t i = 0; i < num_cells; i++) {
        auto [positive_subcell, negative_subcell] = cut_cell(i);
        if (positive_subcell) {
            positive_cells.push_back(std::move(positive_subcell.value()));
        }
        if (negative_subcell) {
            negative_cells.push_back(std::move(negative_subcell.value()));
        }
    }
    logger().debug("{} positive subcells", positive_cells.size());
    logger().debug("{} negative subcells", negative_cells.size());

    // Step 4: combine negative cells into a single cell.
    if (!negative_cells.empty()) {
        logger().debug("combining {} negative cells.", negative_cells.size());
        Cell combined_cell;
        combined_cell.material_label = material_index;
        const size_t total_num_faces = faces.size();
        std::vector<bool> on_border(total_num_faces, false);
        std::vector<size_t> next_face(vertices.size(), INVALID);

        for (const auto& c : negative_cells) {
            for (const auto& f : c.faces) {
                on_border[f] = !on_border[f];
            }
        }

        auto get_ordered_vertex = [&](size_t fid, size_t i) {
            assert(i < 2);
            const auto& f = faces[fid];
            if (f.positive_material_label == material_index) {
                return f.vertices[(i + 1) % 2];
            } else {
                assert(f.negative_material_label == material_index);
                return f.vertices[i];
            }
        };

        size_t num_bd_faces = 0;
        size_t start_face = INVALID;
        for (size_t i = 0; i < total_num_faces; i++) {
            if (!on_border[i]) continue;
            start_face = i;
            next_face[get_ordered_vertex(i, 0)] = i;
            logger().debug("face {}: {} {}", i, get_ordered_vertex(i, 0), get_ordered_vertex(i, 1));
            num_bd_faces++;
        }

        combined_cell.faces.reserve(num_bd_faces);
        size_t fid = start_face;
        for (size_t i = 0; i < num_bd_faces && fid != INVALID; i++) {
            combined_cell.faces.push_back(fid);
            fid = next_face[get_ordered_vertex(fid, 1)];
            assert(fid != INVALID);
        }
        assert(fid == start_face);

        // Update corner boundary vertices.
        for (size_t i = 0; i < num_bd_faces; i++) {
            const auto fid = combined_cell.faces[i];
            const auto vid = get_ordered_vertex(fid, 0);

            auto& v = vertices[vid];
            if (v[0] <= 2 && v[1] <= 2) {
                v[2] = material_index;
            } else if (v[0] <= 2 && v[2] <= 2) {
                v[1] = material_index;
            } else if (v[1] <= 2 && v[2] <= 2) {
                v[0] = material_index;
            }
        }

        // Remove collinear boundary edges.
        auto vertex_on_boundary = [&](const Joint<2>& v) {
            size_t num_bd_material = 0;
            if (v[0] <= 2) num_bd_material++;
            if (v[1] <= 2) num_bd_material++;
            if (v[2] <= 2) num_bd_material++;
            return num_bd_material == 1;
        };

        bool skip_next = false;
        size_t output_face_count = 0;
        for (size_t i = 0; i < num_bd_faces; i++) {
            if (skip_next) {
                skip_next = false;
                continue;
            }
            const auto fid = combined_cell.faces[i];
            const auto v1 = get_ordered_vertex(fid, 1);

            const auto& v = vertices[v1];
            if (v[0] != material_index && v[1] != material_index && v[2] != material_index) {
                assert(vertex_on_boundary(v));
                logger().debug("remove boundary vertex {}: {} {} {}", v1, v[0], v[1] , v[2]);

                const auto next_fid = combined_cell.faces[(i + 1) % num_bd_faces];
                const auto v0 = get_ordered_vertex(fid, 0);
                const auto v2 = get_ordered_vertex(next_fid, 1);

                auto& f = faces[fid];
                f.vertices = {v0, v2};
                if (f.negative_material_label != material_index) {
                    std::swap(f.positive_material_label, f.negative_material_label);
                }
                skip_next = true;
            }
            if (i != output_face_count) {
                combined_cell.faces[output_face_count] = combined_cell.faces[i];
            }
            output_face_count++;
        }
        combined_cell.faces.resize(output_face_count);

        positive_cells.push_back(std::move(combined_cell));
    }
    cells = std::move(positive_cells);

    // Step 5: consolidate.
    consolidate(material_interface);
    logger().debug("After: {} {} {}",
        material_interface.vertices.size(),
        material_interface.faces.size(),
        material_interface.cells.size());
}

template <typename Scalar>
void add_material(MaterialInterfaceBuilder<Scalar, 3>& builder, size_t material_index)
{
    auto& material_interface = builder.get_material_interface();
    const size_t num_vertices = material_interface.vertices.size();

    std::vector<implicit_predicates::Orientation> orientations;
    orientations.reserve(num_vertices);
    // TODO.
}

} // namespace

namespace simplicial_arrangement::internal {

void add_material(MaterialInterfaceBuilder<double, 2>& builder, size_t material_index)
{
    ::add_material(builder, material_index);
}

void add_material(MaterialInterfaceBuilder<Int, 2>& builder, size_t material_index)
{
    ::add_material(builder, material_index);
}

void add_material(MaterialInterfaceBuilder<double, 3>& builder, size_t material_index)
{
    ::add_material(builder, material_index);
}

void add_material(MaterialInterfaceBuilder<Int, 3>& builder, size_t material_index)
{
    ::add_material(builder, material_index);
}

} // namespace simplicial_arrangement::internal
