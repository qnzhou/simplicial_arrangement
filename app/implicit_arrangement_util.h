#pragma once
#include <simplicial_arrangement/simplicial_arrangement.h>
#include <simplicial_arrangement/material_interface.h>

#include <string>

#include <Eigen/Core>

using namespace simplicial_arrangement;



// a polygon in a polygonal mesh
struct PolygonFace
{
    // a list of polygon's vertex indices (index into some global list of vertices)
    std::vector<size_t> vert_indices;
    // the local index of this polygon in all the tets that contains it
    // each pair is (tet_Id, tet_face_Id)
    std::vector<std::pair<size_t, size_t>> tet_face_indices;
};

// vertex of isosurface
struct IsoVert
{
    // the tet containing the IsoVert
    size_t tet_index;
    // the index of IsoVert in tet.vertices
    size_t tet_vert_index;
    // minimal simplex that contains the IsoVert
    size_t simplex_size; // 1: point, 2: edge, 3: triangle, 4: tetrahedron
    // index into a list of tet vertices
    std::array<size_t, 4> simplex_vert_indices;
    // list of implicit functions whose isosurfaces pass IsoVert (indexed into a global list of
    // implicit functions)
    std::array<size_t, 3> func_indices;
};


// vertex of material interface
struct MI_Vert
{
    // the tet containing the MI_Vert
    size_t tet_index;
    // the index of MI_Vert in tet.vertices
    size_t tet_vert_index;
    // minimal simplex that contains the MI_Vert
    size_t simplex_size; // 1: point, 2: edge, 3: triangle, 4: tetrahedron
    // index into a list of tet vertices
    std::array<size_t, 4> simplex_vert_indices;
    // list of materials whose values are equal at MI_Vert (indexed into a global list of
    // material functions)
    std::array<size_t, 4> material_indices;
};


struct Edge
{
    size_t v1;
    size_t v2;
    // each pair is (face_Id, edge_face_Id)
    // face_Id: face index in the global list of faces
    // edge_face_Id: edge index in the face
    std::vector<std::pair<size_t, size_t>> face_edge_indices;
};


// Sphere: (center, radius)
typedef std::pair<std::array<double,3>, double> Sphere;

bool parse_config_file(const std::string &filename,
    std::string& tet_mesh_file,
    std::string& sphere_file,
    std::string& output_dir,
    bool& use_lookup,
    bool& use_2func_lookup,
    bool& use_topo_ray_shooting,
    bool& use_bbox,
    std::array<double,3> &bbox_min,
    std::array<double,3> &bbox_max);

bool parse_config_file_MI(const std::string &filename,
    std::string& tet_mesh_file,
    std::string& material_file,
    std::string& output_dir,
    bool& use_lookup,
    bool& use_3func_lookup);

bool load_tet_mesh(const std::string &filename,
    std::vector<std::array<double, 3>> &pts,
    std::vector<std::array<size_t, 4>> &tets);

bool load_spheres(const std::string &filename,
    std::vector<Sphere> &spheres);

bool load_seeds(const std::string& filename,
    std::vector<std::array<double,3>> &seeds);

inline double compute_Euclidean_distance(const std::array<double,3> &p, const std::array<double,3>& q)
{
    return sqrt((p[0]-q[0])*(p[0]-q[0]) + (p[1]-q[1])*(p[1]-q[1]) + (p[2]-q[2])*(p[2]-q[2]));
}

inline double compute_signed_sphere_distance(const std::array<double, 3>& center, double r, const std::array<double, 3>& p)
{
    return r - compute_Euclidean_distance(center, p);
}

inline double compute_unsigned_sphere_distance(const std::array<double, 3>& center, double r, const std::array<double, 3>& p)
{
    return -abs(r - compute_Euclidean_distance(center, p));
}

inline double compute_sphere_distance(const std::array<double, 3>& center, double r, const std::array<double, 3>& p)
{
    if (r >= 0) {
        return compute_signed_sphere_distance(center, r, p);
    } else {
        return compute_unsigned_sphere_distance(center, -r, p);
    }
}

bool save_result(const std::string& filename,
    const std::vector<std::array<double, 3>>& iso_pts,
    const std::vector<PolygonFace>& iso_faces,
    const std::vector<std::vector<size_t>>& patches,
    const std::vector<Edge>& edges,
    const std::vector<std::vector<size_t>>& chains,
    const std::vector<std::vector<size_t>>& non_manifold_edges_of_vert,
    const std::vector<std::vector<std::pair<size_t, int>>>& half_patch_list,
    const std::vector<std::vector<size_t>>& shells,
    const std::vector<std::vector<size_t>>& components,
    const std::vector<std::vector<size_t>>& arrangement_cells);

bool save_result_MI(const std::string& filename,
    const std::vector<std::array<double, 3>>& MI_pts,
    const std::vector<PolygonFace>& MI_faces,
    const std::vector<std::vector<size_t>>& patches,
    const std::vector<Edge>& edges,
    const std::vector<std::vector<size_t>>& chains,
    const std::vector<std::vector<size_t>>& non_manifold_edges_of_vert,
    const std::vector<std::vector<std::pair<size_t, int>>>& half_patch_list,
    const std::vector<std::vector<size_t>>& shells,
    const std::vector<std::vector<size_t>>& components,
    const std::vector<std::vector<size_t>>& material_cells);


bool save_result_msh(const std::string& filename,
    const std::vector<std::array<double, 3>>& iso_pts,
    const std::vector<PolygonFace>& iso_faces,
    const std::vector<std::vector<size_t>>& patches,
    const std::vector<Edge>& iso_edges,
    const std::vector<std::vector<size_t>>& chains,
    const std::vector<std::vector<size_t>>& non_manifold_edges_of_vert,
    const std::vector<std::vector<std::pair<size_t, int>>>& half_patch_list,
    const std::vector<std::vector<size_t>>& shells,
    const std::vector<std::vector<size_t>>& components,
    const std::vector<std::vector<size_t>>& arrangement_cells);

bool save_nesting_data(const std::string& filename,
    const std::vector<size_t> &next_vert,
    const std::vector<size_t> &extremal_edge_of_component);


bool save_result_mini(const std::string& filename,
    const std::vector<std::array<double, 3>>& iso_pts,
    const std::vector<PolygonFace>& iso_faces,
    const std::vector<std::vector<size_t>>& patches,
    const std::vector<std::vector<size_t>>& arrangement_cells);

// save a list of iso-meshes
bool save_iso_mesh_list(const std::string& filename,
    const std::vector<std::vector<std::array<double, 3>>>& iso_pts_list,
    const std::vector<std::vector<PolygonFace>>& iso_faces_list);


// save a triangle mesh
bool save_tri_mesh(const std::string& filename,
    const std::vector<std::array<double, 3>> &verts,
    const std::vector<std::array<size_t, 3>> &tris);

// extract the boundary triangle mesh of a tet mesh
// assume: the tet mesh represents a simply-connected 3D volume
// assume: the triangles of the boundary mesh don't need to be consistently oriented
// input:
//  verts: num_vert * 3, tet mesh vertices
//  tets: num_tet * 4, tet's corner vertices' indices
// output:
//  boundary_verts: num_boundary_vert * 3, boundary vertices
//  boundary_faces: num_boundary_face * 3, boundary triangle's vertex indices
void extract_tet_boundary_mesh(
    const std::vector<std::array<double, 3>> &verts,
    const std::vector<std::array<size_t, 4>> &tets,
    std::vector<std::array<double, 3>> &boundary_verts, std::vector<std::array<size_t,3>> &boundary_faces);

// save a triangle mesh
bool save_tri_mesh(const std::string& filename,
    const std::vector<std::array<double, 3>> &verts,
    const std::vector<std::array<size_t, 3>> &tris);

// save a list of triangle meshes
bool save_tri_mesh_list(const std::string& filename,
    const std::vector<std::vector<std::array<double, 3>>> &verts_list,
    const std::vector<std::vector<std::array<size_t, 3>>> &tris_list);

bool save_timings(const std::string& filename,
    const std::vector<std::string> &timing_labels,
    const std::vector<double> &timings);

// point (x,y,z): dictionary order
inline bool point_xyz_less(const std::array<double, 3>& p, const std::array<double, 3>& q)
{
    if (p[0] == q[0]) {
        if (p[1] == q[1]) {
            return p[2] < q[2];
        } else {
            return p[1] < q[1];
        }
    } else {
        return p[0] < q[0];
    }
}

// extract the boundary triangle mesh of a tet mesh
// assume: the tet mesh represents a simply-connected 3D volume
// assume: the triangles of the boundary mesh don't need to be consistently oriented
// input:
//  verts: num_vert * 3, tet mesh vertices
//  tets: num_tet * 4, tet's corner vertices' indices
// output:
//  boundary_verts: num_boundary_vert * 3, boundary vertices
//  boundary_faces: num_boundary_face * 3, boundary triangle's vertex indices
void extract_tet_boundary_mesh(
    const std::vector<std::array<double, 3>> &verts,
    const std::vector<std::array<size_t, 4>> &tets,
    std::vector<std::array<double, 3>> &boundary_verts, std::vector<std::array<size_t,3>> &boundary_faces);

// given the list of vertex indices of a face, return the unique key of the face: (the smallest vert Id,
// second-smallest vert Id, the largest vert Id) assume: face_verts is a list of non-duplicate natural
// numbers, with at least three elements.
template <typename IndexType>
void compute_iso_face_key(const std::vector<IndexType>& face_verts, std::array<IndexType, 3>& key)
{
    IndexType min_vert = face_verts[0];
    size_t min_pos = 0;
    IndexType max_vert = face_verts[0];
    for (size_t i = 1; i < face_verts.size(); i++) {
        if (face_verts[i] < min_vert) {
            min_vert = face_verts[i];
            min_pos = i;
        } else if (face_verts[i] > max_vert) {
            max_vert = face_verts[i];
        }
    }
    IndexType second_min_vert = max_vert + 1;
    for (size_t i = 0; i < face_verts.size(); i++) {
        if (i != min_pos && face_verts[i] < second_min_vert) {
            second_min_vert = face_verts[i];
        }
    }
    //
    key[0] = min_vert;
    key[1] = second_min_vert;
    key[2] = max_vert;
}

// extract iso-mesh (topology only) and create map: local index --> global index
//void extract_iso_mesh(const std::vector<bool>& has_isosurface,
//    const std::vector<Arrangement<3>>& cut_results,
//    const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> &func_in_tet,
//    const Eigen::VectorXi &num_func_in_tet,
//    const std::vector<std::array<size_t, 4>>& tets,
//    std::vector<IsoVert>& iso_verts,
//    std::vector<PolygonFace>& iso_faces,
//    std::vector<std::vector<size_t>>& global_vId_of_tet_vert,
//    std::vector<std::vector<size_t>>& iso_fId_of_tet_face);
void extract_iso_mesh(
    size_t num_1_func, size_t num_2_func, size_t num_more_func,
    const std::vector<Arrangement<3>>& cut_results,
    const std::vector<size_t>& cut_result_index,
    const std::vector<size_t>& func_in_tet,
    const std::vector<size_t>& start_index_of_tet,
    const std::vector<std::array<size_t, 4>>& tets,
    std::vector<IsoVert>& iso_verts,
    std::vector<PolygonFace>& iso_faces,
    std::vector<long long>& global_vId_of_tet_vert,
    std::vector<size_t>& global_vId_start_index_of_tet,
    std::vector<size_t>& iso_fId_of_tet_face,
    std::vector<size_t>& iso_fId_start_index_of_tet);

// extract iso-mesh (topology only)
void extract_iso_mesh_pure(
    size_t num_1_func, size_t num_2_func, size_t num_more_func,
    const std::vector<Arrangement<3>>& cut_results,
    const std::vector<size_t>& cut_result_index,
    const std::vector<size_t>& func_in_tet,
    const std::vector<size_t>& start_index_of_tet,
    const std::vector<std::array<size_t, 4>>& tets,
    std::vector<IsoVert>& iso_verts,
    std::vector<PolygonFace>& iso_faces);

// extract material interface mesh (topology only)
void extract_MI_mesh_pure(
    size_t num_2_func, size_t num_3_func, size_t num_more_func,
    const std::vector<MaterialInterface<3>>& cut_results,
    const std::vector<size_t>& cut_result_index,
    const std::vector<size_t>& material_in_tet,
    const std::vector<size_t>& start_index_of_tet,
    const std::vector<std::array<size_t, 4>>& tets,
    std::vector<MI_Vert>& MI_verts,
    std::vector<PolygonFace>& MI_faces);

// extract iso-mesh from marching tet (topology only)
void extract_iso_mesh_marching_tet(const std::vector<bool>& has_isosurface,
    const std::vector<Arrangement<3>>& cut_results,
    const std::vector<std::array<size_t, 4>>& tets,
    std::vector<IsoVert>& iso_verts,
    std::vector<PolygonFace>& iso_faces);

// compute xyz coordinates of iso-vertices
void compute_iso_vert_xyz(
    const std::vector<IsoVert> &iso_verts, 
    const Eigen::MatrixXd &funcVals,
    const std::vector<std::array<double, 3>> &pts,
    std::vector<std::array<double, 3>>& iso_pts);

// compute xyz coordinates of material interface vertices
void compute_MI_vert_xyz(
    const std::vector<MI_Vert> &MI_verts,
    const Eigen::MatrixXd &funcVals,
    const std::vector<std::array<double,3>> &pts,
    std::vector<std::array<double,3>>& MI_pts);

// compute xyz coordinates of iso-vertices from marching tet
void compute_iso_vert_xyz_marching_tet(const std::vector<IsoVert>& iso_verts,
    const std::vector<double>& funcVals,
    const std::vector<std::array<double, 3>>& pts,
    std::vector<std::array<double, 3>>& iso_pts);

// compute iso-edges and edge-face connectivity
//void compute_mesh_edges(std::vector<PolygonFace>& iso_faces, std::vector<Edge>& iso_edges);
//void compute_iso_edges_r(std::vector<PolygonFace>& iso_faces, std::vector<Edge>& iso_edges);

void compute_mesh_edges(const std::vector<PolygonFace>& mesh_faces,
    std::vector<std::vector<size_t>> & edges_of_face,
    std::vector<Edge>& mesh_edges);


// group iso-faces into patches
void compute_patches(const std::vector<std::vector<size_t>> & edges_of_face,
    const std::vector<Edge>& mesh_edges,
    std::vector<std::vector<size_t>>& patches);


// group non-manifold iso-edges into chains
void compute_chains(const std::vector<Edge>& mesh_edges,
    const std::vector<std::vector<size_t>>& non_manifold_edges_of_vert,
    std::vector<std::vector<size_t>>& chains);


// compute neighboring pair of half-faces around an iso-edge in a tetrahedron
// pair<size_t, int> : pair (iso-face index, iso-face orientation)
void compute_face_order_in_one_tet(const Arrangement<3>& tet_cut_result,
    const std::vector<PolygonFace>& iso_faces,
    const Edge& iso_edge,
    std::vector<std::pair<size_t, int>>& ordered_faces);

// compute neighboring pair of half-faces around an edge in a tetrahedron
// pair<size_t, int> : pair (MI-face index, MI-face orientation)
void compute_face_order_in_one_tet_MI(const MaterialInterface<3>& tet_cut_result,
    const std::vector<PolygonFace>& MI_faces,
    const Edge& edge,
    std::vector<std::pair<size_t, int>>& ordered_faces);

// Given tet mesh,
// build the map: v-->v_next, where v_next has lower order than v
void build_next_vert(const std::vector<std::array<double,3>> &pts,
    const std::vector<std::array<size_t,4>> &tets,
    std::vector<size_t> &next_vert);

// compute the order of iso-vertices on a tet edge v->u, v,u in {0,1,2,3}
// return a list of sorted vertex indices {v_id, i1, i2, ..., u_id}
void compute_edge_intersection_order(const Arrangement<3>& tet_cut_result, size_t v, size_t u,
    std::vector<size_t> &vert_indices);

// compute the order of MI-vertices on a tet edge v->u, v,u in {0,1,2,3}
// return a list of sorted vertex indices {v_id, i1, i2, ..., u_id}
void compute_edge_intersection_order_MI(const MaterialInterface<3>& tet_cut_result, size_t v, size_t u,
    std::vector<size_t> &vert_indices);

// find the face passing v, v->u is part of a tet edge, and u is a tet vertex
void compute_passing_face(const Arrangement<3>& tet_cut_result,
    size_t v, size_t u,
    std::pair<size_t,int> &face_orient);

// find the face passing v, v->u is part of a tet edge, and u is a tet vertex
void compute_passing_face_MI(const MaterialInterface<3>& tet_cut_result,
    size_t v, size_t u,
    std::pair<size_t,int> &face_orient);

// find the two faces passing v1 and v2, v1->v2 is part of a tet edge
void compute_passing_face_pair(const Arrangement<3>& tet_cut_result,
    size_t v1, size_t v2,
    std::pair<size_t,int> &face_orient1,
    std::pair<size_t,int> &face_orient2);

// find the two faces passing v1 and v2, v1->v2 is part of a tet edge
void compute_passing_face_pair_MI(const MaterialInterface<3>& tet_cut_result,
    size_t v1, size_t v2,
    std::pair<size_t,int> &face_orient1,
    std::pair<size_t,int> &face_orient2);

//void compute_arrangement_cells(size_t num_patch,
//    const std::vector<std::vector<std::pair<size_t, int>>>& half_patch_list,
//    std::vector<std::vector<size_t>>& arrangement_cells);

// group shells into arrangement cells
void compute_arrangement_cells(size_t num_shell,
    const std::vector<std::pair<size_t,size_t>> &shell_links,
    std::vector<std::vector<size_t>>& arrangement_cells);

// compute shells and connected components of isosurfaces
// each shell is a list of half-patches
// each component is a list of patches
// we also build maps: half-patch --> shell,  patch --> component
void compute_shells_and_components(size_t num_patch,
    const std::vector<std::vector<std::pair<size_t, int>>>& half_patch_list,
    std::vector<std::vector<size_t>>& shells,
    std::vector<size_t>& shell_of_half_patch,
    std::vector<std::vector<size_t>>& components,
    std::vector<size_t>& component_of_patch
);





// compute barycentric coordinate of Point (intersection of three planes)
// Point in tet cell
// template <typename Scalar>
inline void compute_barycentric_coords(const std::array<double, 4>& plane1,
    const std::array<double, 4>& plane2,
    const std::array<double, 4>& plane3,
    std::array<double, 4>& bary_coords)
{
    double n1 = plane1[3] * (plane2[2] * plane3[1] - plane2[1] * plane3[2]) +
                plane1[2] * (plane2[1] * plane3[3] - plane2[3] * plane3[1]) +
                plane1[1] * (plane2[3] * plane3[2] - plane2[2] * plane3[3]);
    double n2 = plane1[3] * (plane2[0] * plane3[2] - plane2[2] * plane3[0]) +
                plane1[2] * (plane2[3] * plane3[0] - plane2[0] * plane3[3]) +
                plane1[0] * (plane2[2] * plane3[3] - plane2[3] * plane3[2]);
    double n3 = plane1[3] * (plane2[1] * plane3[0] - plane2[0] * plane3[1]) +
                plane1[1] * (plane2[0] * plane3[3] - plane2[3] * plane3[0]) +
                plane1[0] * (plane2[3] * plane3[1] - plane2[1] * plane3[3]);
    double n4 = plane1[2] * (plane2[0] * plane3[1] - plane2[1] * plane3[0]) +
                plane1[1] * (plane2[2] * plane3[0] - plane2[0] * plane3[2]) +
                plane1[0] * (plane2[1] * plane3[2] - plane2[2] * plane3[1]);
    double d = n1 + n2 + n3 + n4;
    //
    bary_coords[0] = n1 / d;
    bary_coords[1] = n2 / d;
    bary_coords[2] = n3 / d;
    bary_coords[3] = n4 / d;
}

// Point on tet face
inline void compute_barycentric_coords(const std::array<double, 3>& plane1,
    const std::array<double, 3>& plane2,
    std::array<double, 3>& bary_coords)
{
    double n1 = plane1[2] * plane2[1] - plane1[1] * plane2[2];
    double n2 = plane1[0] * plane2[2] - plane1[2] * plane2[0];
    double n3 = plane1[1] * plane2[0] - plane1[0] * plane2[1];
    double d = n1 + n2 + n3;
    //
    bary_coords[0] = n1 / d;
    bary_coords[1] = n2 / d;
    bary_coords[2] = n3 / d;
}

// Point on tet edge
inline void compute_barycentric_coords(double f1, double f2, std::array<double, 2>& bary_coords)
{
    bary_coords[0] = f2 / (f2 - f1);
    bary_coords[1] = 1 - bary_coords[0];
}


// implicit functions
inline double sphere_function(
    const std::array<double, 3>& center, double r, const std::array<double, 3>& p)
{
    return (p[0] - center[0]) * (p[0] - center[0]) + (p[1] - center[1]) * (p[1] - center[1]) +
           (p[2] - center[2]) * (p[2] - center[2]) - r * r;
}

inline int sign(double x)
{
    return (x > 0) ? 1 : ((x < 0) ? -1 : 0);
}

// tetrahedron dual contouring
// input:
// pts & tets: tet mesh
// funcVals: |pts| * |num_func| matrix, function values at tet vertices
// funcSigns: |pts| * |num_func| matrix, true if function j is positive at vertex i
// output:
// mesh_verts
// mesh_tris
void tet_dual_contouring(const std::vector<std::array<double,3>> &pts,
    const std::vector<std::array<size_t,4>> &tets,
    const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> &funcVals,
    const Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> &funcSigns,
    std::vector<std::array<double,3>> &mesh_verts,
    std::vector<std::array<size_t,3>> &mesh_tris);

// tetrahedron dual contouring
// input:
// pts & tets: tet mesh
// funcVals: |pts| * |num_func| matrix, function values at tet vertices
// funcSigns: |pts| * |num_func| matrix, true if function j is positive at vertex i
// func_in_tet and start_index_of_tet: function indices in each tet, stored as CRS vector.
// output:
// mesh_verts
// mesh_tris
void tet_dual_contouring(size_t num_1func, size_t num_2func, size_t num_more_func,
    const std::vector<std::array<double,3>> &pts,
    const std::vector<std::array<size_t,4>> &tets,
    const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> &funcVals,
    const Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> &funcSigns,
    const std::vector<size_t> &func_in_tet,
    const std::vector<size_t> &start_index_of_tet,
    std::vector<std::array<double,3>> &mesh_verts,
    std::vector<std::array<size_t,3>> &mesh_tris);