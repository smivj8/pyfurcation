import numpy as np
from stylianou_geometry_equations import *
from bifurcation_vectors_class import bifurcation_vector


def generate_pointclouds(parameters, n_circ = 100, n_streamline = 120, 
                           n_car_range = 80, n_car_rad = 8, scalar = 1.1):
    R_p = parameters[0]; L_p = parameters[1]; R_d = parameters[2]
    L_d = parameters[3]; R_o = parameters[4]; R_i = parameters[5]
    iota_b = parameters[6]; iota_gamma = parameters[8]
    
    L_p = 0.05*L_d
    L_d = 1.05*L_d
    
    carina_parameters = [R_p, R_d, R_o, R_i, iota_b]

    phi_s_pos_min = np.arctan(-1*L_p/R_o)  
    #phi_s_neg_min = 0

    eta = np.tan(iota_gamma - iota_b)
    phi_s_pos_max = iota_b + np.arctan(eta + (L_d/R_o))
    #phi_s_neg_max = phi_s_pos_max
    
    phi_s_pos = np.linspace(phi_s_pos_min, phi_s_pos_max, n_streamline)
    # phi_s_neg = np.linspace(phi_s_neg_min, phi_s_neg_max, n_streamline)
    phi_c = np.linspace(-1*np.pi/2, 3*np.pi/2, n_circ)
    
    phi_s_carina_max = iota_gamma
    phi_s_carina = np.linspace(0, phi_s_carina_max, n_car_range)
    
    
    xyz_pos = np.zeros([n_circ*n_streamline, 3])
    xyz_pos_normals = np.zeros([n_circ*n_streamline, 3])
    
    car_pos = np.zeros([n_car_range*n_car_rad, 3])
    car_pos_normals =  np.zeros([n_car_range*n_car_rad, 3])
    
    car_neg = np.zeros([n_car_range*n_car_rad, 3])
    car_neg_normals = np.zeros([n_car_range*n_car_rad, 3])


    # for j in range(len(phi_c)):
    #     for i in range(len(phi_s_pos)):
    #         xyz_pos[j*n_streamline+i,:], xyz_pos_normals[j*n_streamline+i,:] = r_positive(phi_s_pos[i], phi_c[j], 
    #                                          phi_c_pos_min, phi_c_pos_max, 
    #                                          parameters, scalar)
            
    # Radial Test
    iteration = 0
    for i in range(len(phi_s_pos)):
        for j in range(len(phi_c)):
            xyz_pos[iteration,:], xyz_pos_normals[iteration,:] = r_positive(phi_s_pos[i], phi_c[j], 
                                             phi_c_pos_min, phi_c_pos_max, 
                                             parameters, scalar)
            iteration += 1
    
    for m in range(len(phi_s_carina)):
        _, omega, _, _ = rho_omega_theta_parameters(phi_s_carina[m], parameters, carina_parameters, scalar)
        psi_range = np.linspace(-1*omega, omega, n_car_rad)
        for n in range(len(psi_range)):
            car_pos[n*n_car_range + m, :], car_pos_normals[n*n_car_range + m, :] = carina_positive(
                phi_s_carina[m], psi_range[n], parameters, carina_parameters, scalar)
            car_neg[n*n_car_range + m, :], car_neg_normals[n*n_car_range + m, :] = carina_negative(
                phi_s_carina[m], psi_range[n], parameters, carina_parameters, scalar)
    
    
    xyz_pc_positive = np.concatenate((xyz_pos, xyz_pos_normals), axis = 1)
    
    carina_pc_positive = np.concatenate((car_pos, car_pos_normals), axis = 1)
    carina_pc_negative = np.concatenate((car_neg, car_neg_normals), axis = 1)
    
    xyz_pos_trimmed = xyz_pc_positive[~np.isnan(xyz_pc_positive).any(axis=1),:] 
    xyz_neg_trimmed = np.zeros(np.shape(xyz_pos_trimmed))

    xyz_neg_trimmed[:,0] = xyz_pos_trimmed[:,0]
    xyz_neg_trimmed[:,1] = xyz_pos_trimmed[:,1]
    xyz_neg_trimmed[:,2] = -1*xyz_pos_trimmed[:,2]
    xyz_neg_trimmed[:,3] = xyz_pos_trimmed[:,3]
    xyz_neg_trimmed[:,4] = xyz_pos_trimmed[:,4]
    xyz_neg_trimmed[:,5] = -1*xyz_pos_trimmed[:,5]
    
    car_pos_trimmed = carina_pc_positive[~np.isnan(carina_pc_positive).any(axis=1),:]
    car_neg_trimmed = carina_pc_negative[~np.isnan(carina_pc_negative).any(axis=1),:]
    
    # bifurcation_pointcloud = np.concatenate((xyz_pos_trimmed, xyz_neg_trimmed))
    # carina_pointcloud = np.concatenate((car_pos_trimmed, car_neg_trimmed))
    
    
    pointcloud = np.concatenate((xyz_pos_trimmed, xyz_neg_trimmed, car_pos_trimmed, car_neg_trimmed))
    
    #Shift pointcloud to x: 0
    x_min = np.min(pointcloud[:,0])
    
    pointcloud[:,0] = pointcloud[:,0] - x_min - L_p
    

    
    return pointcloud



def pointcloud_to_open3d(pointcloud):
    pcd = o3d.geometry.PointCloud()
    if len(pointcloud[0]) == 6:
        pcd.points = o3d.utility.Vector3dVector(pointcloud[:,:3])
        pcd.normals = o3d.utility.Vector3dVector(pointcloud[:,3:])
    else:
        pcd.points = o3d.utility.Vector3dVector(pointcloud)
    return pcd


def poisson_algorithm(pointcloud, depth = 7, width = 0, scale = 1.1, 
                      linear_fit = False, visualization = False):
    #Also requires normals
    pcd = pointcloud_to_open3d(pointcloud)
    #computing the mesh
    poisson_mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_poisson(pcd, depth = depth, width = width,
                                                                             scale = scale, linear_fit = linear_fit)[0]
    #visualization???
    #poisson_mesh.compute_triangle_normals()
    #decimate to 100,000 triangles
    #poisson_mesh.simplify_quadric_decimation(100000)
    return poisson_mesh, pcd

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

def crop_bifurcation_poisson_mesh(mesh, parameters):
    iota_b = parameters[6]
    #crop inlet end
    mesh = clean_crop_x(mesh, min_x = 0, min_slice = True)
    #Crop Daughter ends
        # z pos crop
        #find x_max
    
    position_vec = bifurcation_vector(parameters, False)
    r_D_pos, _ = position_vec.calculate_outlet_positions()
    normal_vec = bifurcation_vector(parameters, True)
    n_D_pos, n_D_neg = normal_vec.calculate_outlet_normals()
    
    cos_theta = np.dot(n_D_pos, r_D_pos)/np.linalg.norm(r_D_pos)
    hypot_max = np.linalg.norm(r_D_pos)*cos_theta
    #rotate poisson cloud and do the first daughter mesh crop

    mesh.rotate((rotation_matrix_about_y(iota_b)), center = np.zeros(3))
    mesh = clean_crop_x(mesh, max_x = hypot_max, max_slice = True)
    
    mesh.rotate((rotation_matrix_about_y(-1*iota_b)), center = np.zeros(3))
    mesh.rotate((rotation_matrix_about_y(-1*iota_b)), center = np.zeros(3))
    mesh = clean_crop_x(mesh, max_x = hypot_max, max_slice = True)
    
    mesh.rotate((rotation_matrix_about_y(iota_b)), center = np.zeros(3))

    return mesh


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


D_0 = 100; R_0 = D_0/2; R_1 = D_0/2.48888; L_0 = 0; L_1 = 4.08888*R_1
R_o_0 = 5*R_1; iota_b_test = 35*np.pi/180; iota_gamma_test = 39.83 * np.pi/180
delta_alpha_test = 3.6*np.pi/180

D_1 = R_1*2; R_2 = D_1/2.48888; L_0 = 0; L_2 = 4.08888*R_2
R_o_1 = 5*R_2

parameters_0 = [R_0, 0, R_1, L_1, R_o_0, R_i, iota_b_test, delta_alpha_test, iota_gamma_test, 1]
parameters_1 = [R_1, 0, R_2, L_2, R_o_1, R_i, iota_b_test, delta_alpha_test, iota_gamma_test, 1]

xyz_0 = generate_pointclouds(parameters_0)
mesh_0, pcd_0 = poisson_algorithm(xyz_0)
mesh_0 = crop_bifurcation_poisson_mesh(mesh_0, parameters_0)

o3d.visualization.draw_geometries([mesh_0], window_name = "Visualization",
                                   mesh_show_wireframe = True, mesh_show_back_face = True)