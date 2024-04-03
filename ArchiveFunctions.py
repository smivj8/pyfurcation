# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 16:11:24 2024

@author: SMIVJ8
"""
def fuse_two_meshes(mesh_A, mesh_B, mesh_B_parameters, fuse_location, fuse_normal, fuse_range):
    #  fuse distance might be wrong?
    # is the correct fuse distance the x value of the rotated fuse vector? 
    # wait duh yes that makes sense. leaving this comment here in case I was indeed 
    # wrong and need to redo it
    
    #Rotate meshes from original to crop
    
    R1 = unit_vector_match_rotation(fuse_normal, np.array([1,0,0]))
    fuse_distance = (R1 @ fuse_location)[0]
    mesh_A.rotate(R1, center = np.zeros(3))
    mesh_B.rotate(R1, center = np.zeros(3))
    
    #Initialize Variables from meshes
    mesh_A_vertices = np.asarray(mesh_A.vertices)
    mesh_A_triangles = np.asarray(mesh_A.triangles)
    
    mesh_B_vertices = np.asarray(mesh_B.vertices)
    mesh_B_triangles = np.asarray(mesh_B.triangles)
    

    
    new_A_triangles = []
    new_B_triangles = []
    
    fuse_A_vertices = []
    fuse_B_vertices = []
    
    # TODO add bounds to crop plane (of mesh A only, mesh B is a single bifurcation unit
    # so it should be fine. I don't want to crop earlier generations on accident!!!
    # Need some sort of metric to measure the yz coordinates of the triangle in relation
    # to the fuse location. Radius of mesh B would be perfect
    mesh_radius = mesh_B_parameters[0]
    fuse_y = (R1 @ fuse_location)[1]
    fuse_z = (R1 @ fuse_location)[2]
    

    
    # find axes to compute
    axis = 0    
    
    # Evaluate mesh A triangles
    value_A = fuse_distance - (fuse_range/2)
    for t_A in mesh_A_triangles:
        #NOTE: t is a 1D array [v1_index, v2_index, v3_index]! Not an integer
        #If someone else is debugging this, that is mainly a note for me, I did
        #understand how iterables worked lol.
        v1_index, v2_index, v3_index =  t_A
        #Set out boolian to TRUE if vertex is outside of range
        # Set to FALSE if in range
        v1_out_bool = mesh_A_vertices[v1_index][axis] > value_A
        v2_out_bool = mesh_A_vertices[v2_index][axis] > value_A
        v3_out_bool = mesh_A_vertices[v3_index][axis] > value_A
    
        if v1_out_bool:
            v1_y = mesh_A_vertices[v1_index][1]
            v1_z = mesh_A_vertices[v1_index][2]
            v1_out_bool = (np.abs(v1_y - fuse_y) <= 2 * mesh_radius) and \
                np.abs(v1_z - fuse_z) <= 2 * mesh_radius
        else:
            v1_out_bool = mesh_A_vertices[v1_index][axis] > value_A
        if v2_out_bool:
            v2_y = mesh_A_vertices[v2_index][1]
            v2_z = mesh_A_vertices[v2_index][2]
            v2_out_bool = np.abs(v2_y - fuse_y) <= 2 * mesh_radius and \
                np.abs(v2_z - fuse_z) <= 2 * mesh_radius
        else:
            v2_out_bool = mesh_A_vertices[v2_index][axis] > value_A
        if v3_out_bool:
            v3_y = mesh_A_vertices[v3_index][1]
            v3_z = mesh_A_vertices[v3_index][2]
            v3_out_bool = np.abs(v3_y - fuse_y) <= 2 * mesh_radius and \
                np.abs(v3_z - fuse_z) <= 2 * mesh_radius
        else:
            v3_out_bool = mesh_A_vertices[v3_index][axis] > value_A

        
    
        bool_sum = sum([v1_out_bool, v2_out_bool, v3_out_bool])
        
        if bool_sum == 0:
            #Triangle is completely within range, add it to new list of triangles
            #t is a triangle, not an integer!
            new_A_triangles.append(t_A)
        elif bool_sum == 3:
            #Triangle is completely outside of range, ignore and move on
            #I could leave this alse the else condition I guess, but idk
            #I'm copying the guys code above me
            continue
        elif bool_sum == 2:
            # 2 vertices are outside of range, squish the two vertices!
            if not v1_out_bool:
                #vertex 1 is in, vertex 2 and vertex 3 are out
                mesh_A_vertices[v2_index, 0] = value_A
                mesh_A_vertices[v3_index, 0] = value_A
                fuse_A_vertices.append(mesh_A_vertices[v2_index])
                fuse_A_vertices.append(mesh_A_vertices[v3_index])
            elif not v2_out_bool:
                #vertex 2 is in, vertex 1 and vertex 3 are out
                mesh_A_vertices[v1_index, 0] = value_A
                mesh_A_vertices[v3_index, 0] = value_A
                fuse_A_vertices.append(mesh_A_vertices[v1_index])
                fuse_A_vertices.append(mesh_A_vertices[v3_index])
            elif not v3_out_bool:
                #vertex 3 is in, vertex 1 and vertex 2 are out
                mesh_A_vertices[v1_index, 0] = value_A
                mesh_A_vertices[v2_index, 0] = value_A
                fuse_A_vertices.append(mesh_A_vertices[v1_index])
                fuse_A_vertices.append(mesh_A_vertices[v2_index])
            else:
                print("bool sum incorrectly calculated somehow - 2 verts out")
                assert False
            new_A_triangles.append(t_A)
        elif bool_sum == 1:
            #1 vertex is outside of range, squish vertex into place
            if v1_out_bool:
                #vertex 1 is out
                mesh_A_vertices[v1_index, 0] = value_A
                fuse_A_vertices.append(mesh_A_vertices[v1_index])
            elif v2_out_bool:
                #vertex 2 is outside of range
                mesh_A_vertices[v2_index, 0] = value_A
                fuse_A_vertices.append(mesh_A_vertices[v2_index])
            elif v3_out_bool:
                #vertex 3 out of range
                mesh_A_vertices[v3_index, 0] = value_A
                fuse_A_vertices.append(mesh_A_vertices[v3_index])
            else:
                print("bool sum incorrectly calculated somehow - 1 vert out")
                assert False
            new_A_triangles.append(t_A)

        else:
            print("Bool sum not in range")
            assert False
        
    value_B = fuse_distance + (fuse_range/2)
    for t_B in mesh_B_triangles:
        #NOTE: t is a 1D array [v1_index, v2_index, v3_index]! Not an integer
        #If someone else is debugging this, that is mainly a note for me, I did
        #not understand how iterables worked lol.
        v1_index, v2_index, v3_index =  t_B
        #Set out boolian to TRUE if vertex is outside of range
        # Set to FALSE if in range
        v1_out_bool = mesh_B_vertices[v1_index][axis] < value_B
        v2_out_bool = mesh_B_vertices[v2_index][axis] < value_B
        v3_out_bool = mesh_B_vertices[v3_index][axis] < value_B
    
        bool_sum = sum([v1_out_bool, v2_out_bool, v3_out_bool])
        
        if bool_sum == 0:
            #Triangle is completely within range, add it to new list of triangles
            #t is a triangle, not an integer!
            new_B_triangles.append(t_B)
        elif bool_sum == 3:
            #Triangle is completely outside of range, ignore and move on
            #I could leave this alse the else condition I guess, but idk
            #I'm copying the guys code above me
            
            #Okay but this does not remove the free vertices????
            #might not be an issue, could be auto resolved with o3d
            continue
        elif bool_sum == 2:
            # 2 vertices are outside of range, squish the two vertices!
            if not v1_out_bool:
                #vertex 1 is in, vertex 2 and vertex 3 are out
                mesh_B_vertices[v2_index, 0] = value_B
                mesh_B_vertices[v3_index, 0] = value_B
                fuse_B_vertices.append(mesh_B_vertices[v2_index])
                fuse_B_vertices.append(mesh_B_vertices[v3_index])
            elif not v2_out_bool:
                #vertex 2 is in, vertex 1 and vertex 3 are out
                mesh_B_vertices[v1_index, 0] = value_B
                mesh_B_vertices[v3_index, 0] = value_B
                fuse_B_vertices.append(mesh_B_vertices[v1_index])
                fuse_B_vertices.append(mesh_B_vertices[v3_index])
            elif not v3_out_bool:
                #vertex 3 is in, vertex 1 and vertex 2 are out
                mesh_B_vertices[v1_index, 0] = value_B
                mesh_B_vertices[v2_index, 0] = value_B
                fuse_B_vertices.append(mesh_B_vertices[v1_index])
                fuse_B_vertices.append(mesh_B_vertices[v2_index])
            else:
                print("bool sum incorrectly calculated somehow - 2 verts out")
                assert False
            new_B_triangles.append(t_B)
        elif bool_sum == 1:
            #1 vertex is outside of range, squish vertex into place
            if v1_out_bool:
                #vertex 1 is out
                mesh_B_vertices[v1_index, 0] = value_B
                fuse_B_vertices.append(mesh_B_vertices[v1_index])
            elif v2_out_bool:
                #vertex 2 is outside of range
                mesh_B_vertices[v2_index, 0] = value_B
                fuse_B_vertices.append(mesh_B_vertices[v2_index])
            elif v3_out_bool:
                #vertex 3 out of range
                mesh_B_vertices[v3_index, 0] = value_B
                fuse_B_vertices.append(mesh_B_vertices[v3_index])
            else:
                print("bool sum incorrectly calculated somehow - 1 vert out")
                assert False
            new_B_triangles.append(t_B)

        else:
            print("Bool sum not in range")
            assert False
    
    mesh_A = o3d.geometry.TriangleMesh()
    mesh_A.vertices = o3d.utility.Vector3dVector(mesh_A_vertices)
    mesh_A.triangles = o3d.utility.Vector3iVector(np.array(new_A_triangles))
    
    mesh_B = o3d.geometry.TriangleMesh()
    mesh_B.vertices = o3d.utility.Vector3dVector(mesh_B_vertices)
    mesh_B.triangles = o3d.utility.Vector3iVector(np.array(new_B_triangles))
    
    # TODO actually fuse the two meshes together !!!!!!
    #convert list of vertices to numpy array
    fuse_A_vertices = np.array(fuse_A_vertices)
    fuse_B_vertices = np.array(fuse_B_vertices)
    
    
    #Shift y,z to be about origin by subtracting from the list of vertices 
    #the shifted fuse location. This should be right, and should be independent 
    #of mesh radius? Think cartesian to polar using arcos (range (0, 2pi))
    shifted_A_vertices = fuse_A_vertices - (R1 @ fuse_location)
    shifted_B_vertices = fuse_B_vertices - (R1 @ fuse_location)

    n_A = np.shape(fuse_A_vertices)[0]
    n_B = np.shape(fuse_B_vertices)[0]
    
    theta_A = np.zeros(n_A)
    theta_B = np.zeros(n_B)
    
    #There should be EXACTLY this many triangles I think!!
    fuse_triangles = np.zeros((n_A + n_B, 3))
    

    #Parameterize the yz coordinate by theta about (x, 0, 0). This gets a 1D 
    #"ordering" of the points
    for j in range(n_A):
        z = shifted_A_vertices[j,2]
        y = shifted_A_vertices[j,1]
        H = np.sqrt(shifted_A_vertices[j,1]**2 + shifted_A_vertices[j,2]**2)
        assert H != 0.0
        theta_A[j] = np.arccos(z/H) + np.pi * bool(y < 0)
    for i in range(n_B):
        z = shifted_B_vertices[i,2]
        y = shifted_B_vertices[i,1]
        H = np.sqrt(shifted_B_vertices[i,1]**2 + shifted_B_vertices[i,2]**2)
        assert H != 0.0
        theta_B[i] = np.arccos(z/H) + np.pi * bool(y < 0)
    
    #Sort the fuseA and fuseB vertices by thetaA and thetaB respectively
    #also sort thetaA and theta B of course
    sort_index_A = theta_A.argsort()
    sort_index_B = theta_B.argsort()
    
    sorted_A_vertices = fuse_A_vertices[sort_index_A]
    theta_A.sort()
    sorted_B_vertices = fuse_B_vertices[sort_index_B]
    theta_B.sort()
    fuse_vertices = np.concatenate((sorted_A_vertices, sorted_B_vertices))
    index_A = 0
    index_B = 0
    downward_indices = np.zeros(n_A)
    
    # print(theta_A)
    # print(theta_B)
    
    #Create DOWNWARD TRIANGLES (there will be nA total, create the first nA-1)
    for d in range(n_A - 1):
        v1_index = d
        v2_index = d+1
        #print(d); print(index_B); print("\n")
        index_B = 0
        while theta_B[index_B] <= theta_A[d]:
            index_B += 1
        v3_index = n_A + index_B
        downward_indices[d] = index_B
        fuse_triangles[d] = np.array([v1_index, v2_index, v3_index])
    
    #create the final downward triangle
    v1_index = 0; v2_index = (n_A - 1); v3_index = (n_A + n_B - 1)
    downward_indices[n_A - 1] = n_B - 1
    fuse_triangles[n_A - 1] = np.array([v1_index, v2_index, v3_index])
    
    
    
    #Okay, I have them sorted. Now, I need to finish up and create the triangles
    #I am tired and don't feel like thinking through this logic right now tbh
    #so problem for tomorrow!
    
    mesh_fuse = o3d.geometry.TriangleMesh()
    mesh_fuse.vertices = o3d.utility.Vector3dVector(fuse_vertices)
    mesh_fuse.triangles = o3d.utility.Vector3iVector(fuse_triangles)
    
    mesh = mesh_A + mesh_B + mesh_fuse
    
    R2 = unit_vector_match_rotation(np.array([1,0,0]), fuse_normal)
    mesh.rotate(R2, center = np.zeros(3))
   

    return mesh



def combine_parent_daughter_mesh(parent_mesh, daughter_mesh, r_P, n_P, 
                                 daughter_parameters, z_positive = True):
    if z_positive:
        r_D, _ = outlet_positions(daughter_parameters)
        n_D, _ = outlet_normals(daughter_parameters) 
    else:
        _, r_D = outlet_positions(daughter_parameters)
        _, n_D = outlet_normals(daughter_parameters)
    
        # r_D_pos, r_D_neg = outlet_positions(parameters_test)

        # n_D_pos, n_D_neg = outlet_normals(parameters_test)

    
    #Rotate daughter mesh [1,0,0] to match the normal face of the parent
    R_1 = unit_vector_match_rotation(np.array([1,0,0]), n_P)
    daughter_mesh.rotate(R_1, center = np.array([0,0,0]))
    #rotate random angle about PARENT normal axis
    # TODO make this a custom input not random...
    theta = np.deg2rad(180*np.random.rand())
    R_2 = rotation_about_vector_u(n_P, theta)
    daughter_mesh.rotate(R_2, center = np.array([0,0,0]))
    #translate [0,0,0] of daughter to r_P
    daughter_mesh.translate(r_P)
    #COMBINE MESHES
    
    
    #merge_clos_vertices
    
    #update r and n
    r_new = (R_2 @ (R_1 @ r_D)) + r_P
    n_new = R_2 @ (R_1 @ n_D)
    
    #there's no way it's this fucking easy...
    #combined_mesh = parent_mesh + daughter_mesh
    #it was LMAO
    #BUT I changed my mind, I'm not going to combine them
    
    #NOW the harder part: cleaning up the two geometries
    #combined_mesh.merge_close_vertices(eps = merge_eps)
    #combined_mesh.simplify_vertex_clustering(voxel_size = merge_eps)
    
    parent_mesh, daughter_mesh = remesh_junction(parent_mesh, daughter_mesh, 
                                                                r_P, n_P, daughter_parameters, 300)
    #TODO replace junciton
    junction_mesh = None
    
    return parent_mesh, daughter_mesh, junction_mesh, r_new, n_new

def construct_meshes(tree_parameters):
    n_gen = len(tree_parameters)
    meshes = [0]*n_gen
    
    for k in range(n_gen):
        pointcloud = generate_pointclouds(tree_parameters[k])
        mesh, _ = poisson_algorithm(pointcloud)
        print(k)
        mesh = crop_bifurcation_poisson_mesh(mesh, tree_parameters[k])
        meshes[k] = mesh
        
        del(mesh)
        del(pointcloud)
        
    return meshes

def define_branching_structure(tree_parameters, axial_rotations, r_init = np.zeros(3),
                               n_init = np.array([1,0,0])):
    assert len(tree_parameters) == len(axial_rotations)
    
    #Initialize Variable
    n_gen = len(tree_parameters)
    r_k_pos = np.zeros((n_gen, 3))
    r_k_neg = np.zeros((n_gen, 3))
    n_k_pos = np.zeros((n_gen, 3))
    n_k_neg = np.zeros((n_gen, 3))
    
    r_mat_init = unit_vector_match_rotation(np.array([1,0,0]), n_init)
    r_axial_init = rotation_about_vector_u(n_init, axial_rotations[0])
    r_k_pos[0] = (r_axial_init @ (r_mat_init @ outlet_positions(tree_parameters[0])[0])) + r_init
    r_k_neg[0] = (r_axial_init @ (r_mat_init @ outlet_positions(tree_parameters[0])[1])) + r_init

    n_k_pos[0] = r_mat_init @ (outlet_normals(tree_parameters[0])[0])
    n_k_neg[0] = r_mat_init @ (outlet_normals(tree_parameters[0])[1])
    
    
    # #update r and n
    # r_new = (R_2 @ (R_1 @ r_D)) + r_P
    # n_new = R_2 @ (R_1 @ n_D)
    
    
    #Define Positions
    for k in range(1, n_gen):
        r_pos_D, r_neg_D = outlet_positions(tree_parameters[k])
        n_pos_D, n_neg_D =  outlet_normals(tree_parameters[k])
        R1 = unit_vector_match_rotation(np.array([1,0,0]), n_pos_D)
        R2 = rotation_about_vector_u(n_k_pos[k-1], axial_rotations[k])
        r_k_pos[k] = (R2 @ (R1 @ r_pos_D)) + r_k_pos[k-1]
        r_k_neg[k] = (R2 @ (R1 @ r_neg_D)) + r_k_pos[k-1]
        n_k_pos[k] = (R2 @ (R1 @ n_pos_D))
        n_k_neg[k] = (R2 @ (R1 @ n_neg_D))
        
    return r_k_pos, r_k_neg, n_k_pos, n_k_neg


def construct_bronchi_tree(tree_parameters, r_init = np.zeros(3),
                           n_init = np.array([1,0,0])):
    # Check that the right number of parameters was input
    n_gen = len(tree_parameters)
    #initialize lists and arrays
    r_k_pos = np.zeros((n_gen, 3))
    r_k_neg = np.zeros((n_gen, 3))
    n_k_pos = np.zeros((n_gen, 3))
    n_k_neg = np.zeros((n_gen, 3))
    meshes = [0]*n_gen
    junctions = [0] * (n_gen-1)
    
    # TODO figure out shift if we are not starting from r_init, n_init
    
    r_mat_init = unit_vector_match_rotation(np.array([1,0,0]), n_init)
    
    #gen 0 pointcloud and meshes
    r_k_pos[0], r_k_neg[0] = outlet_positions(tree_parameters[0]) + r_init
    n_k_pos[0] = r_mat_init @ (outlet_normals(tree_parameters[0])[0])
    n_k_neg[0] = r_mat_init @ (outlet_normals(tree_parameters[0])[1])
    
    pcloud_k = generate_pointclouds(tree_parameters[0])
    meshes[0], _ = poisson_algorithm(pcloud_k)
    meshes[0] = crop_bifurcation_poisson_mesh(meshes[0], tree_parameters[0])
    del(pcloud_k)
    
    for k in range(1, n_gen):
        pcloud_k = generate_pointclouds(tree_parameters[k])
        meshes[k], _ = poisson_algorithm(pcloud_k)
        del(pcloud_k)
        meshes[k] = crop_bifurcation_poisson_mesh(meshes[k], tree_parameters[k])
        
        meshes[k-1], meshes[k], junctions[k-1], r_k_pos[k], n_k_pos[k] = \
            combine_parent_daughter_mesh(meshes[k-1], meshes[k], r_k_pos[k-1],
                                         n_k_pos[k-1], tree_parameters[k])
        print(k)
        
    total_mesh = meshes[0]
    
    for j in range(1, n_gen):
        total_mesh += meshes[j]
        #total_mesh += junctions[j - 1]
        
    
    return total_mesh, r_k_pos, r_k_neg, n_k_pos, n_k_neg


def bbox_transformation(pcd, bbox, pos_norms):
    #Get vertices
    bbox_vertices = bbox.get_box_points()
    bbox_vertices_np = np.asarray(bbox_vertices)
    print(type(bbox_vertices))
    bbox_sorted_index = np.lexsort((-1*bbox_vertices_np[:,1], bbox_vertices_np[:,2], bbox_vertices_np[:,0]))

    bbox_sorted = bbox_vertices_np[bbox_sorted_index]
    pos_norms[0][1] = 0
    normal = pos_norms[0]
    d = np.dot(pos_norms[0], pos_norms[1])

    bbox_sorted[0] = np.array([0, bbox_sorted[0,1], (d/normal[2]) ])
    bbox_sorted[1] = np.array([0, bbox_sorted[1,1], (d/normal[2]) ])

    bbox_sorted[2] = np.array([0, bbox_sorted[2,1], 0])
    bbox_sorted[3] = np.array([0, bbox_sorted[3,1], 0])

    bbox_sorted[4] = np.array([(d/normal[0]), bbox_sorted[4,1], 0 ])
    bbox_sorted[5] = np.array([(d/normal[0]), bbox_sorted[5,1], 0 ])

    bbox_sorted[6] = np.array([0, bbox_sorted[6,1], -1*(d/normal[2]) ])
    bbox_sorted[7] = np.array([0, bbox_sorted[7,1], -1*(d/normal[2]) ])

    bbox_shifted_vertices_np = np.zeros(np.shape(bbox_vertices_np))

    for j in range(len(bbox_sorted_index)):
        bbox_shifted_vertices_np[j] = bbox_sorted[bbox_sorted_index[j]]
    
    #convert points BACK to utility  vector
    #bbox_shifted_o3d = o3d.utility.Vector3dVector(bbox_shifted_vertices_np)
    
    return bbox_shifted_vertices_np

# def delaunay_surface_triangulation(parameters,  n_circ = 80, n_streamline = 80, 
#                            n_car_range = 48, n_car_rad = 16):
    
#     xyz_pos, xyz_neg, car_pos, car_neg = generate_pointclouds(parameters, 
#                                                                n_circ = n_circ, 
#                                                                n_streamline = n_streamline, 
#                                                                n_car_range = n_car_range, 
#                                                                n_car_rad = n_car_rad)
    
#     len_xyz = np.shape(xyz_pos)[0]
    
#     y_pos = np.concatenate((xyz_pos[0:int(len_xyz/2), :], xyz_neg[0:int(len_xyz/2),:], car_pos))
#     y_neg = np.concatenate((xyz_pos[int(len_xyz/2):, :], xyz_neg[int(len_xyz/2):,:], car_neg))

#     pos_sorted_index = np.lexsort((y_pos[:,0], y_pos[:,2], y_pos[:,1]))
#     neg_sorted_index = np.lexsort((y_neg[:,0], y_neg[:,2], -1*y_neg[:,1]))

#     y_pos = y_pos[pos_sorted_index]
#     y_neg = y_neg[neg_sorted_index]
    
#     #Don't forget to duplicate some of the vertices HERE, but skipping for now
#     #to try to simply get the surface triangulation
    
    
#     xz_pos_projection = np.zeros((np.shape(y_pos)[0], 2))
#     xz_neg_projection = np.zeros((np.shape(y_neg)[0], 2))
#     xz_pos_projection[:,0] = y_pos[:,0]
#     xz_pos_projection[:,1] = y_pos[:,2]
#     xz_neg_projection[:,0] = y_neg[:,0]
#     xz_neg_projection[:,1] = y_neg[:,2]
    
#     # fig = plt.figure(figsize =(10,10))
#     # ax = fig.add_subplot(111, projection='3d')
    
#     tri_pos = spatial.Delaunay(xz_pos_projection)
#     tri_neg = spatial.Delaunay(xz_neg_projection)
    
#     # ax.plot_trisurf(
#     #     y_pos[:,0], y_pos[:,1], y_pos[:,2],
#     #     triangles=tri_pos.simplices, cmap="magma"
#     # )
#     # plt.show()
    
#     # ax.plot_trisurf(
#     #     y_neg[:,0], y_neg[:,1], y_neg[:,2],
#     #     triangles=tri_neg.simplices, cmap="magma"
#     # )
#     # plt.show()
    
#     return tri_pos, tri_neg, y_pos, y_neg

    def compute_intersection_mesh_A(vertex_in_index, vertex_out_index):
        vertex_in = mesh_A_vertices[vertex_in_index]
        vertex_out = mesh_A_vertices[vertex_out_index]
        if (vertex_in_index, vertex_out_index) in intersection_A_edges:
            intersection_A_index = intersection_A_edges[(vertex_in_index, vertex_out_index)]
            intersection = new_A_vertices[intersection_A_index]
        else:
            value_A = fuse_distance - (fuse_range/2)
            intersection = [None, None, None]
            intersection[0] = value_A
            const_1 = (value_A - vertex_in[axis])/(vertex_out[axis] - vertex_in[axis])
            c = axes_compute[0]
            intersection[c] = (const_1 * (vertex_out[c] - vertex_in[c])) + vertex_in[c]
            c = axes_compute[1]
            intersection[c] = (const_1 * (vertex_out[c] - vertex_in[c])) + vertex_in[c]
            assert not (None in intersection)
            # save new vertice and remember that this intersection already added an edge
            new_A_vertices.append(intersection)
            intersection_A_index = len(new_A_vertices) - 1
            intersection_A_edges[(vertex_in_index, vertex_out_index)] = intersection_A_index

        return intersection_A_index
    
    def compute_intersection_mesh_B(vertex_in_index, vertex_out_index):
        vertex_in = mesh_B_vertices[vertex_in_index]
        vertex_out = mesh_B_vertices[vertex_out_index]
        if (vertex_in_index, vertex_out_index) in intersection_B_edges:
            intersection_B_index = intersection_B_edges[(vertex_in_index, vertex_out_index)]
            intersection = new_B_vertices[intersection_B_index]
        else:
            value_B = fuse_distance + (fuse_range/2)
            intersection = [None, None, None]
            intersection[axis] = value_B
            const_1 = (value_B - vertex_in[axis])/(vertex_out[axis] - vertex_in[axis])
            c = axes_compute[0]
            intersection[c] = (const_1 * (vertex_out[c] - vertex_in[c])) + vertex_in[c]
            c = axes_compute[1]
            intersection[c] = (const_1 * (vertex_out[c] - vertex_in[c])) + vertex_in[c]
            assert not (None in intersection)
            # save new vertice and remember that this intersection already added an edge
            new_B_vertices.append(intersection)
            intersection_index = len(new_B_vertices) - 1
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