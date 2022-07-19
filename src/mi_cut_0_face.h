#pragma once

#include "MIComplex.h"
#include "MaterialRepo.h"
#include "robust_assert.h"
#include "utils.h"

#include <implicit_predicates/implicit_predicates.h>

namespace simplicial_arrangement {

template <typename Scalar>
int8_t mi_cut_0_face(const MaterialRepo<Scalar, 2>& materials,
    const MIComplex<2>& mi_complex,
    size_t vid,
    size_t material_index)
{
    ROBUST_ASSERT(material_index > 2);
    const auto& vertices = mi_complex.vertices;
    const auto& p = vertices[vid];
    const auto& material = materials.get_material(material_index);

    short vertex_type = 0;
    if (p[0] > 2) vertex_type |= 1;
    if (p[1] > 2) vertex_type |= 2;
    if (p[2] > 2) vertex_type |= 4;

    auto get_corner_id = [](size_t i, size_t j) -> size_t {
        ROBUST_ASSERT(i <= 2 && j <= 2);
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
        ROBUST_ASSERT(i <= 2);
        switch (i) {
        case 0: {
            const Scalar mm0[]{m0[1], m0[2]};
            const Scalar mm1[]{m1[1], m1[2]};
            const Scalar mm[]{material[1], material[2]};
#ifdef SIMPLICIAL_ARRANGEMENT_NON_ROBUST
            return implicit_predicates::mi_orient1d_nonrobust(mm0, mm1, mm);
#else
            return implicit_predicates::mi_orient1d(mm0, mm1, mm);
#endif
        }
        case 1: {
            const Scalar mm0[]{m0[0], m0[2]};
            const Scalar mm1[]{m1[0], m1[2]};
            const Scalar mm[]{material[0], material[2]};
#ifdef SIMPLICIAL_ARRANGEMENT_NON_ROBUST
            return implicit_predicates::mi_orient1d_nonrobust(mm0, mm1, mm);
#else
            return implicit_predicates::mi_orient1d(mm0, mm1, mm);
#endif
        }
        case 2: {
            const Scalar mm0[]{m0[0], m0[1]};
            const Scalar mm1[]{m1[0], m1[1]};
            const Scalar mm[]{material[0], material[1]};
#ifdef SIMPLICIAL_ARRANGEMENT_NON_ROBUST
            return implicit_predicates::mi_orient1d_nonrobust(mm0, mm1, mm);
#else
            return implicit_predicates::mi_orient1d(mm0, mm1, mm);
#endif
        }
        default:
            logger().error("Invalid simplex boundary: {}", i);
            return implicit_predicates::INVALID;
        }
    };

    switch (vertex_type) {
    case 1: {
        const auto& m0 = materials.get_material(p[0]);
        return utils::signof(compute_orientation_0d(p[1], p[2], m0));
    }
    case 2: {
        const auto& m1 = materials.get_material(p[1]);
        return utils::signof(compute_orientation_0d(p[0], p[2], m1));
    }
    case 3: {
        const auto& m0 = materials.get_material(p[0]);
        const auto& m1 = materials.get_material(p[1]);
        return utils::signof(compute_orientation_1d(p[2], m0, m1));
    }
    case 4: {
        const auto& m2 = materials.get_material(p[2]);
        return utils::signof(compute_orientation_0d(p[0], p[1], m2));
    }
    case 5: {
        const auto& m0 = materials.get_material(p[0]);
        const auto& m2 = materials.get_material(p[2]);
        return utils::signof(compute_orientation_1d(p[1], m0, m2));
    }
    case 6: {
        const auto& m1 = materials.get_material(p[1]);
        const auto& m2 = materials.get_material(p[2]);
        return utils::signof(compute_orientation_1d(p[0], m1, m2));
    }
    case 7: {
        const auto& m0 = materials.get_material(p[0]);
        const auto& m1 = materials.get_material(p[1]);
        const auto& m2 = materials.get_material(p[2]);
#ifdef SIMPLICIAL_ARRANGEMENT_NON_ROBUST
        return utils::signof(implicit_predicates::mi_orient2d_nonrobust(
            m0.data(), m1.data(), m2.data(), material.data()));
#else
        return utils::signof(
            implicit_predicates::mi_orient2d(m0.data(), m1.data(), m2.data(), material.data()));
#endif
    }
    default:
        logger().error("Impossible vertex type case detected: {}", vertex_type);
        return utils::signof(implicit_predicates::INVALID);
    }
}

template <typename Scalar>
int8_t mi_cut_0_face(const MaterialRepo<Scalar, 3>& materials,
    const MIComplex<3>& mi_complex,
    size_t vid,
    size_t material_index)
{
    ROBUST_ASSERT(material_index > 3);
    const auto& vertices = mi_complex.vertices;
    const auto& p = vertices[vid];
    const auto& material = materials.get_material(material_index);

    short vertex_type = 0;
    if (p[0] > 3) vertex_type |= 1;
    if (p[1] > 3) vertex_type |= 2;
    if (p[2] > 3) vertex_type |= 4;
    if (p[3] > 3) vertex_type |= 8;

    auto get_corner_id = [](size_t i, size_t j, size_t k) -> size_t {
        ROBUST_ASSERT(i <= 3 && j <= 3 && k <= 3);
        if (i != 0 && j != 0 && k != 0) return 0;
        if (i != 1 && j != 1 && k != 1) return 1;
        if (i != 2 && j != 2 && k != 2) return 2;
        if (i != 3 && j != 3 && k != 3) return 3;
        logger().error("Invalid simplex boundary inputs: {}, {}, {}", i, j, k);
        return INVALID;
    };

    auto compute_orientation_0d = [&](size_t i, size_t j, size_t k, const Material<Scalar, 3>& m)
        -> implicit_predicates::Orientation {
        const size_t corner_id = get_corner_id(i, j, k);
        if (m[corner_id] > material[corner_id])
            return implicit_predicates::POSITIVE;
        else if (m[corner_id] < material[corner_id])
            return implicit_predicates::NEGATIVE;
        else
            return implicit_predicates::ZERO;
    };
    auto compute_orientation_1d =
        [&](size_t i,
            size_t j,
            const Material<Scalar, 3>& m0,
            const Material<Scalar, 3>& m1) -> implicit_predicates::Orientation {
        ROBUST_ASSERT(i <= 3 && j <= 3);
        ROBUST_ASSERT(i != j);
        short t = 0;
        if (i == 0 || j == 0) t |= 1;
        if (i == 1 || j == 1) t |= 2;
        if (i == 2 || j == 2) t |= 4;
        if (i == 3 || j == 3) t |= 8;
        switch (t) {
        case 3: {
            const Scalar mm0[]{m0[2], m0[3]};
            const Scalar mm1[]{m1[2], m1[3]};
            const Scalar mm[]{material[2], material[3]};
#ifdef SIMPLICIAL_ARRANGEMENT_NON_ROBUST
            return implicit_predicates::mi_orient1d_nonrobust(mm0, mm1, mm);
#else
            return implicit_predicates::mi_orient1d(mm0, mm1, mm);
#endif
        }
        case 5: {
            const Scalar mm0[]{m0[1], m0[3]};
            const Scalar mm1[]{m1[1], m1[3]};
            const Scalar mm[]{material[1], material[3]};
#ifdef SIMPLICIAL_ARRANGEMENT_NON_ROBUST
            return implicit_predicates::mi_orient1d_nonrobust(mm0, mm1, mm);
#else
            return implicit_predicates::mi_orient1d(mm0, mm1, mm);
#endif
        }
        case 6: {
            const Scalar mm0[]{m0[0], m0[3]};
            const Scalar mm1[]{m1[0], m1[3]};
            const Scalar mm[]{material[0], material[3]};
#ifdef SIMPLICIAL_ARRANGEMENT_NON_ROBUST
            return implicit_predicates::mi_orient1d_nonrobust(mm0, mm1, mm);
#else
            return implicit_predicates::mi_orient1d(mm0, mm1, mm);
#endif
        }
        case 9: {
            const Scalar mm0[]{m0[1], m0[2]};
            const Scalar mm1[]{m1[1], m1[2]};
            const Scalar mm[]{material[1], material[2]};
#ifdef SIMPLICIAL_ARRANGEMENT_NON_ROBUST
            return implicit_predicates::mi_orient1d_nonrobust(mm0, mm1, mm);
#else
            return implicit_predicates::mi_orient1d(mm0, mm1, mm);
#endif
        }
        case 10: {
            const Scalar mm0[]{m0[0], m0[2]};
            const Scalar mm1[]{m1[0], m1[2]};
            const Scalar mm[]{material[0], material[2]};
#ifdef SIMPLICIAL_ARRANGEMENT_NON_ROBUST
            return implicit_predicates::mi_orient1d_nonrobust(mm0, mm1, mm);
#else
            return implicit_predicates::mi_orient1d(mm0, mm1, mm);
#endif
        }
        case 12: {
            const Scalar mm0[]{m0[0], m0[1]};
            const Scalar mm1[]{m1[0], m1[1]};
            const Scalar mm[]{material[0], material[1]};
#ifdef SIMPLICIAL_ARRANGEMENT_NON_ROBUST
            return implicit_predicates::mi_orient1d_nonrobust(mm0, mm1, mm);
#else
            return implicit_predicates::mi_orient1d(mm0, mm1, mm);
#endif
        }
        default:
            logger().error("Invalid 3D simplex edge: {}, {} with key {}", i, j, t);
            return implicit_predicates::INVALID;
        }
    };
    auto compute_orientation_2d =
        [&](size_t i,
            const Material<Scalar, 3>& m0,
            const Material<Scalar, 3>& m1,
            const Material<Scalar, 3>& m2) -> implicit_predicates::Orientation {
        switch (i) {
        case 0: {
            const Scalar mm0[]{m0[1], m0[2], m0[3]};
            const Scalar mm1[]{m1[1], m1[2], m1[3]};
            const Scalar mm2[]{m2[1], m2[2], m2[3]};
            const Scalar mm[]{material[1], material[2], material[3]};
#ifdef SIMPLICIAL_ARRANGEMENT_NON_ROBUST
            return implicit_predicates::mi_orient2d_nonrobust(mm0, mm1, mm2, mm);
#else
            return implicit_predicates::mi_orient2d(mm0, mm1, mm2, mm);
#endif
        }
        case 1: {
            const Scalar mm0[]{m0[0], m0[2], m0[3]};
            const Scalar mm1[]{m1[0], m1[2], m1[3]};
            const Scalar mm2[]{m2[0], m2[2], m2[3]};
            const Scalar mm[]{material[0], material[2], material[3]};
#ifdef SIMPLICIAL_ARRANGEMENT_NON_ROBUST
            return implicit_predicates::mi_orient2d_nonrobust(mm0, mm1, mm2, mm);
#else
            return implicit_predicates::mi_orient2d(mm0, mm1, mm2, mm);
#endif
        }
        case 2: {
            const Scalar mm0[]{m0[0], m0[1], m0[3]};
            const Scalar mm1[]{m1[0], m1[1], m1[3]};
            const Scalar mm2[]{m2[0], m2[1], m2[3]};
            const Scalar mm[]{material[0], material[1], material[3]};
#ifdef SIMPLICIAL_ARRANGEMENT_NON_ROBUST
            return implicit_predicates::mi_orient2d_nonrobust(mm0, mm1, mm2, mm);
#else
            return implicit_predicates::mi_orient2d(mm0, mm1, mm2, mm);
#endif
        }
        case 3: {
            const Scalar mm0[]{m0[0], m0[1], m0[2]};
            const Scalar mm1[]{m1[0], m1[1], m1[2]};
            const Scalar mm2[]{m2[0], m2[1], m2[2]};
            const Scalar mm[]{material[0], material[1], material[2]};
#ifdef SIMPLICIAL_ARRANGEMENT_NON_ROBUST
            return implicit_predicates::mi_orient2d_nonrobust(mm0, mm1, mm2, mm);
#else
            return implicit_predicates::mi_orient2d(mm0, mm1, mm2, mm);
#endif
        }
        default:
            logger().error("Invalid 3D simplex face: {}", i);
            return implicit_predicates::INVALID;
        }
    };

    switch (vertex_type) {
    case 1: {
        const auto& m0 = materials.get_material(p[0]);
        return utils::signof(compute_orientation_0d(p[1], p[2], p[3], m0));
    }
    case 2: {
        const auto& m1 = materials.get_material(p[1]);
        return utils::signof(compute_orientation_0d(p[0], p[2], p[3], m1));
    }
    case 3: {
        const auto& m0 = materials.get_material(p[0]);
        const auto& m1 = materials.get_material(p[1]);
        return utils::signof(compute_orientation_1d(p[2], p[3], m0, m1));
    }
    case 4: {
        const auto& m2 = materials.get_material(p[2]);
        return utils::signof(compute_orientation_0d(p[0], p[1], p[3], m2));
    }
    case 5: {
        const auto& m0 = materials.get_material(p[0]);
        const auto& m2 = materials.get_material(p[2]);
        return utils::signof(compute_orientation_1d(p[1], p[3], m0, m2));
    }
    case 6: {
        const auto& m1 = materials.get_material(p[1]);
        const auto& m2 = materials.get_material(p[2]);
        return utils::signof(compute_orientation_1d(p[0], p[3], m1, m2));
    }
    case 7: {
        const auto& m0 = materials.get_material(p[0]);
        const auto& m1 = materials.get_material(p[1]);
        const auto& m2 = materials.get_material(p[2]);
        return utils::signof(compute_orientation_2d(p[3], m0, m1, m2));
    }
    case 8: {
        const auto& m3 = materials.get_material(p[3]);
        return utils::signof(compute_orientation_0d(p[0], p[1], p[2], m3));
    }
    case 9: {
        const auto& m0 = materials.get_material(p[0]);
        const auto& m3 = materials.get_material(p[3]);
        return utils::signof(compute_orientation_1d(p[1], p[2], m0, m3));
    }
    case 10: {
        const auto& m1 = materials.get_material(p[1]);
        const auto& m3 = materials.get_material(p[3]);
        return utils::signof(compute_orientation_1d(p[0], p[2], m1, m3));
    }
    case 11: {
        const auto& m0 = materials.get_material(p[0]);
        const auto& m1 = materials.get_material(p[1]);
        const auto& m3 = materials.get_material(p[3]);
        return utils::signof(compute_orientation_2d(p[2], m0, m1, m3));
    }
    case 12: {
        const auto& m2 = materials.get_material(p[2]);
        const auto& m3 = materials.get_material(p[3]);
        return utils::signof(compute_orientation_1d(p[0], p[1], m2, m3));
    }
    case 13: {
        const auto& m0 = materials.get_material(p[0]);
        const auto& m2 = materials.get_material(p[2]);
        const auto& m3 = materials.get_material(p[3]);
        return utils::signof(compute_orientation_2d(p[1], m0, m2, m3));
    }
    case 14: {
        const auto& m1 = materials.get_material(p[1]);
        const auto& m2 = materials.get_material(p[2]);
        const auto& m3 = materials.get_material(p[3]);
        return utils::signof(compute_orientation_2d(p[0], m1, m2, m3));
    }
    case 15: {
        const auto& m0 = materials.get_material(p[0]);
        const auto& m1 = materials.get_material(p[1]);
        const auto& m2 = materials.get_material(p[2]);
        const auto& m3 = materials.get_material(p[3]);
#ifdef SIMPLICIAL_ARRANGEMENT_NON_ROBUST
        return utils::signof(implicit_predicates::mi_orient3d_nonrobust(
            m0.data(), m1.data(), m2.data(), m3.data(), material.data()));
#else
        return utils::signof(implicit_predicates::mi_orient3d(
            m0.data(), m1.data(), m2.data(), m3.data(), material.data()));
#endif
    }
    default:
        logger().error("Impossible vertex type case detected: {}", vertex_type);
        return utils::signof(implicit_predicates::INVALID);
    }
}

} // namespace simplicial_arrangement
