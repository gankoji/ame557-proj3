import numpy as np
import math

d2r = math.pi/180.0

# Inputs are from the 30 day propagation of the orbit from TLE
r0 = np.array( [ 1895874.38838166,  4497007.59724938, -4507886.60518927])
v0 = np.array( [-7120.17746759,   3051.88529362,    73.34216649])
oe0 = np.array( [6.64582965e6, 2.06999581e-3, d2r*42.7320639,
                 d2r*156.211622, d2r*8.19853023, d2r*80.7678293])

# km, NA, deg, deg, deg, deg
oef = np.array([6657391.879, 0.002595, d2r*42.748082,
                d2r*345.325883, d2r*124.412552, d2r*287.718564])

print(oef - oe0)

