# -*- coding: utf-8 -*-
"""
Created on Sat Mar  2 22:15:43 2024

@author: SMIVJ8
"""

import numpy as np

# D_p = 10; R_p = D_p/2; R_d = D_p/2.48888; L_p = 3.92857*R_p; L_d = 4.08888*R_d
# R_o = 5*R_d; iota_b = 35*np.pi/180; iota_gamma = 39.83 * np.pi/180
# delta_alpha = 3.6*np.pi/180


# g0_diameter = 2

arr1 = np.array([[1,2],
                 [3,4],
                 [5,6],
                 [6,7],
                 [8,9],
                 [10,11],
                 [12,13]])

arr2 = np.array([3,9,4,7,5,1,9])

sort_index = arr2.argsort()

sorted_arr1 = arr1[sort_index]
arr2.sort()
