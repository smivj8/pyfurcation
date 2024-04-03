# -*- coding: utf-8 -*-
"""
Samuel Smith
MATH6019 - Nonlinear Programming
Date: Thu Feb  1 17:22:48 2024

Assignment
Problem

Description:

Author: smit3sl@mail.uc.edu
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.spatial as spatial
from stl import mesh


L_p = 0.75
R_p = 0.25
n = 200

zeta = np.linspace(0, L_p, n)
phi_c = np.linspace(-1*np.pi/2, 3*np.pi/2, n)


def coord(zeta, phi_c):
    x = zeta
    y = R_p * np.cos(phi_c)
    z = R_p * np.sin(phi_c)
    return np.array([x,y,z])

xyz = np.zeros([n*n, 3])



for j in range(n):
    for i in range(n):
        xyz[j*n+i,:] = coord(zeta[j], phi_c[i])



fig = plt.figure(figsize =(10,10))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(xyz[:,0],xyz[:,1],xyz[:,2])
# for k in range(6):
#     ax.text(xn[k], yn[k], zn[k], str(k))
#ax.plot_trisurf(xn,yn,zn, triangles = tri.simplices)
plt.title('Test')

plt.show()

convexhull = spatial.ConvexHull(xyz)
faces = convexhull.simplices

myramid_mesh = mesh.Mesh(
  np.zeros(faces.shape[0], dtype=mesh.Mesh.dtype)
)
for i, f in enumerate(faces):
  for j in range(3):
    myramid_mesh.vectors[i][j] = xyz[f[j],:]
    
myramid_mesh.save('numpy_stl_example_02.stl')

#triangulation = spatial.Delaunay(xyz)

