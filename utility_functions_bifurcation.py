import numpy as np
import open3d as o3d

def sliceplane(mesh, axis, value, direction):
    # axis can be 0,1,2 (which corresponds to x,y,z)
    # value where the plane is on that axis
    # direction can be True or False (True means remove everything that is
    # greater, False means less
    # than)

    vertices = np.asarray(mesh.vertices)
    triangles = np.asarray(mesh.triangles)
    new_vertices = list(vertices)
    new_triangles = []

    # (a, b) -> c
    # c refers to index of new vertex that sits at the intersection between a,b
    # and the boundingbox edge
    # a is always inside and b is always outside
    intersection_edges = dict()

    # find axes to compute
    axes_compute = [0,1,2]
    # remove axis that the plane is on
    axes_compute.remove(axis)

    def compute_intersection(vertex_in_index, vertex_out_index):
        vertex_in = vertices[vertex_in_index]
        vertex_out = vertices[vertex_out_index]
        if (vertex_in_index, vertex_out_index) in intersection_edges:
            intersection_index = intersection_edges[(vertex_in_index, vertex_out_index)]
            intersection = new_vertices[intersection_index]
        else:
            intersection = [None, None, None]
            intersection[axis] = value
            const_1 = (value - vertex_in[axis])/(vertex_out[axis] - vertex_in[axis])
            c = axes_compute[0]
            intersection[c] = (const_1 * (vertex_out[c] - vertex_in[c])) + vertex_in[c]
            c = axes_compute[1]
            intersection[c] = (const_1 * (vertex_out[c] - vertex_in[c])) + vertex_in[c]
            assert not (None in intersection)
            # save new vertice and remember that this intersection already added an edge
            new_vertices.append(intersection)
            intersection_index = len(new_vertices) - 1
            intersection_edges[(vertex_in_index, vertex_out_index)] = intersection_index

        return intersection_index

    for t in triangles:
        v1, v2, v3 = t
        if direction:
            v1_out = vertices[v1][axis] > value
            v2_out = vertices[v2][axis] > value
            v3_out = vertices[v3][axis] > value
        else: 
            v1_out = vertices[v1][axis] < value
            v2_out = vertices[v2][axis] < value
            v3_out = vertices[v3][axis] < value

        bool_sum = sum([v1_out, v2_out, v3_out])
        # print(f"{v1_out=}, {v2_out=}, {v3_out=}, {bool_sum=}")

        if bool_sum == 0:
            # triangle completely inside --> add and continue
            new_triangles.append(t)
        elif bool_sum == 3:
            # triangle completely outside --> skip
            continue
        elif bool_sum == 2:
            # two vertices outside 
            # add triangle using both intersections
            vertex_in_index = v1 if (not v1_out) else (v2 if (not v2_out) else v3)
            vertex_out_1_index = v1 if v1_out else (v2 if v2_out else v3)
            vertex_out_2_index = v3 if v3_out else (v2 if v2_out else v1)
            # print(f"{vertex_in_index=}, {vertex_out_1_index=}, {vertex_out_2_index=}")
            # small sanity check if indices sum matches
            assert sum([vertex_in_index, vertex_out_1_index, vertex_out_2_index]) == sum([v1,v2,v3])

            # add new triangle 
            new_triangles.append([vertex_in_index, compute_intersection(vertex_in_index, vertex_out_1_index), 
                compute_intersection(vertex_in_index, vertex_out_2_index)])

        elif bool_sum == 1:
            # one vertice outside
            # add three triangles
            vertex_out_index = v1 if v1_out else (v2 if v2_out else v3)
            vertex_in_1_index = v1 if (not v1_out) else (v2 if (not v2_out) else v3)
            vertex_in_2_index = v3 if (not v3_out) else (v2 if (not v2_out) else v1)
            # print(f"{vertex_out_index=}, {vertex_in_1_index=}, {vertex_in_2_index=}")
            # small sanity check if outdices sum matches
            assert sum([vertex_out_index, vertex_in_1_index, vertex_in_2_index]) == sum([v1,v2,v3])

            new_triangles.append([vertex_in_1_index, compute_intersection(vertex_in_1_index, vertex_out_index), vertex_in_2_index])
            new_triangles.append([compute_intersection(vertex_in_1_index, vertex_out_index), 
                compute_intersection(vertex_in_2_index, vertex_out_index), vertex_in_2_index])

        else:
            assert False

    # remap indices and remove unused 

    mesh = o3d.geometry.TriangleMesh()
    mesh.vertices = o3d.utility.Vector3dVector(np.array(new_vertices))
    mesh.triangles = o3d.utility.Vector3iVector(np.array(new_triangles))
    return mesh

def clean_crop_x(mesh, min_x = None, max_x = None, min_slice = False, max_slice = False):
    #min_x = min(min_corner[0], max_corner[0])
    #min_y = min(min_corner[1], max_corner[1])
    #max_x = max(min_corner[0], max_corner[0])
    #max_y = max(min_corner[1], max_corner[1])
    mesh_sliced = mesh
    # mesh = sliceplane(mesh, 0, min_x, False)
    if max_slice:
        mesh_sliced = sliceplane(mesh_sliced, 0, max_x, True)
    if min_slice:    
        mesh_sliced = sliceplane(mesh_sliced, 0, min_x, False)
    #mesh_sliced = sliceplane(mesh_sliced, 1, max_y, True)
    #mesh_sliced = sliceplane(mesh_sliced, 1, min_y, False)
    # mesh_sliced = mesh_sliced.paint_uniform_color([0,0,1])
    return mesh_sliced


def rotation_matrix_about_y(theta):
    #If anyone ever reads this source code, I accidentally had a fucking 1 
    #instead of a zero in the last row of this 
    R_y = np.array([[np.cos(theta),      0, np.sin(theta)],
                    [0,                  1, 0            ],
                    [-1* np.sin(theta), 0, np.cos(theta)]])
    return R_y

def unit_vector_match_rotation(a, b):
    #ROTATES UNIT VECTOR A ONTO UNIT VECTOR B
    u = np.cross(a, b)
    s = np.linalg.norm(u)
    c = np.dot(a, b)
    if np.linalg.norm(b - a) == 0:
        R = np.identity(3)
    else:
        U_mat = np.array([[0,      -1*u[2], u[1]   ],
                          [u[2],    0,      -1*u[0]],
                          [-1*u[1], u[0],   0      ]])
        R = np.identity(3) + U_mat + (U_mat @ U_mat) * ((1-c)/s**2)
    return R

def rotation_about_vector_u(u, theta):
    U_mat = np.array([[0,      -1*u[2], u[1]   ],
                      [u[2],    0,      -1*u[0]],
                      [-1*u[1], u[0],   0      ]])
    R = np.cos(theta)*np.identity(3) + np.sin(theta)*U_mat + (1-np.cos(theta)) * np.outer(u, u)    
    return R