import numpy as np

def my_func(x,y,z):
    return np.array([2*x, 2*y, 2*z])

[new_x, _, new_z] = my_func(1,2,3)
new_y = 'lol'
print(new_x)
print(new_y)
print(new_z)