#include "add_material.h"
#include "MaterialInterfaceBuilder.h"

#include <implicit_predicates/implicit_predicates.h>

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
        if (material[corner_id] > m[corner_id])
            return implicit_predicates::POSITIVE;
        else if (material[corner_id] < m[corner_id])
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

template <typename Scalar>
void add_material(MaterialInterfaceBuilder<Scalar, 2>& builder, size_t material_index)
{
    const auto& material = builder.get_material(material_index);
    auto& material_interface = builder.get_material_interface();
    const size_t num_vertices = material_interface.vertices.size();
    const size_t num_faces = material_interface.faces.size();
    const size_t num_cells = material_interface.cells.size();
    logger().info("Before: {} {} {}", num_vertices, num_faces, num_cells);

    auto& vertices = material_interface.vertices;
    auto& faces = material_interface.faces;
    auto& cells = material_interface.cells;

    std::vector<implicit_predicates::Orientation> orientations;
    orientations.reserve(num_vertices);
    for (size_t i = 0; i < num_vertices; i++) {
        const auto& p = vertices[i];
        orientations.push_back(compute_orientation(builder, p, material));
        assert(orienttions.back() != implicit_predicates::INVALID);
    }

    std::vector<size_t> material_points(num_faces, INVALID);
    std::vector<size_t> face_map(num_faces, INVALID);
    std::vector<MaterialInterface<2>::Face> new_faces;
    new_faces.reserve(num_faces * 2);

    for (size_t i = 0; i < num_faces; i++) {
        const auto& f = faces[i];
        const auto& end_points = f.vertices;

        if (orientations[end_points[0]] == implicit_predicates::ZERO) {
            // Material interface pass through existing joint.
            material_points[i] = end_points[0];
            if (orientations[end_points[1]] == implicit_predicates::NEGATIVE) {
                new_faces.push_back(f);
                face_map.push_back(new_faces.size() - 1);
            }
        } else if (orientations[end_points[1]] == implicit_predicates::ZERO) {
            // Material interface pass through existing joint.
            material_points[i] = end_points[1];
            if (orientations[end_points[0]] == implicit_predicates::NEGATIVE) {
                new_faces.push_back(f);
                face_map.push_back(new_faces.size() - 1);
            }
        } else if (orientations[end_points[0]] != orientations[end_points[1]]) {
            // New joint!
            size_t m0 = cells[f.positive_cell].material_label;
            size_t m1 = cells[f.negative_cell].material_label;
            vertices.push_back({m0, m1, material_index});
            material_points[i] = vertices.size() - 1;
            if (orientations[end_points[0]] == implicit_predicates::NEGATIVE) {
                MaterialInterface<2>::Face face;
                face.vertices = {end_points[0], material_points[i]};
                face.positive_cell = f.positive_cell;
                face.negative_cell = f.negative_cell;
                new_faces.push_back(std::move(face));
            } else {
                assert(orientations[end_points[1]] == implicit_predicates::NEGATIVE);
                MaterialInterface<2>::Face face;
                face.vertices = {material_points[i], end_points[1]};
                face.positive_cell = f.positive_cell;
                face.negative_cell = f.negative_cell;
                new_faces.push_back(std::move(face));
            }
            face_map.push_back(new_faces.size() - 1);
        } else if (orientations[end_points[0]] == implicit_predicates::NEGATIVE) {
            new_faces.push_back(f);
            face_map.push_back(new_faces.size() - 1);
        }
    }

    std::vector<MaterialInterface<2>::Cell> new_cells;
    new_cells.reserve(num_cells * 2);
    MaterialInterface<2>::Cell cut_cell;
    cut_cell.material_label = material_index;
    cut_cell.faces.reserve(num_cells * 2);

    for (size_t i = 0; i < num_cells; i++) {
        const auto& cell = cells[i];
        const size_t num_cell_faces = cell.faces.size();

        // Identify cut face and its location.
        MaterialInterface<2>::Face cut_face;
        cut_face.vertices.reserve(2);
        size_t cut_face_index = INVALID;

        for (size_t j = 0; j < num_cell_faces && cut_face_index == INVALID; j++) {
            const auto fid = cell.faces[j];
            const size_t vid = material_points[fid];
            if (vid != INVALID) {
                for (size_t k = 1; k < num_cell_faces; k++) {
                    const size_t idx = (j + k) % num_cell_faces;
                    const auto next_fid = cell.faces[idx];
                    const size_t next_vid = material_points[next_fid];
                    if (next_vid != INVALID && next_vid != vid) {
                        cut_face.vertices = {vid, next_vid};
                        cut_face.positive_cell = INVALID;
                        cut_face.negative_cell = i;
                        cut_face_index = idx;
                        cut_cell.faces.push_back(new_faces.size());
                        new_faces.push_back(std::move(cut_face));
                        break;
                    }
                }
            }
        }

        // Based on the location of cut face, add new faces to new cell.
        MaterialInterface<2>::Cell new_cell;
        new_cell.material_label = cell.material_label;
        new_cell.faces.reserve(num_cell_faces + 1);
        if (cut_face_index == INVALID) {
            // No cut face, starting face does not matter.
            cut_face_index = 0;
        } else {
            // Always start from cut face.
            new_cell.faces.push_back(new_faces.size() - 1);
        }

        for (size_t j = 0; j < num_cell_faces; j++) {
            const size_t fid = (cut_face_index + j) % num_cell_faces;
            const size_t new_fid = face_map[fid];

            if (new_fid != INVALID) {
                new_cell.faces.push_back(new_fid);
                assert(new_faces[new_cell.faces.back()].positive_cell == material_index ||
                       new_faces[new_cell.faces.back()].negative_cell == material_index ||);
            }
        }
        new_cells.push_back(std::move(new_cell));
    }

    // Lastly, reorder the faces in cut cell and we are done.
    const size_t num_cut_faces = cut_cell.faces.size();
    if (num_cut_faces > 0) {
        assert(num_cut_faces >= 3);
        const size_t cut_cell_id = new_cells.size();
        absl::flat_hash_map<size_t, size_t> next_face;
        next_face.reserve(num_cut_faces);
        for (size_t i = 0; i < num_cut_faces; i++) {
            auto& f = new_faces[cut_cell.faces[i]];
            assert(f.positive_cell == INVALID);
            assert(f.negative_cell != INVALID);
            f.positive_cell = cut_cell_id;
            next_face[f.vertices[1]] = i;
        }

        std::vector<size_t> ordered_faces;
        ordered_faces.reserve(num_cut_faces);

        size_t curr_index = 0;
        for (size_t i = 0; i < num_cut_faces; i++) {
            ordered_faces.push_back(cut_cell.faces[curr_index]);
            size_t curr_vid = new_faces[ordered_faces.back()].vertices[0];
            const auto itr = next_face.find(curr_vid);
            assert(itr != next_face.end());
            curr_index = itr->second;
        }
        assert(curr_index == 0);
        cut_cell.faces = std::move(ordered_faces);
        new_cells.push_back(std::move(cut_cell));
    }
    material_interface.faces = std::move(new_faces);
    material_interface.cells = std::move(new_cells);
    logger().info("After: {} {} {}",
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
