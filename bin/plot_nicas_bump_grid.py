# plot some bump data
# 
# D. J. Lea
# Oct 2022

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np

testdata="/data/users/sophia.moreton/jedi/mo-bundle-dev/orca-jedi/src/tests/Data"
# testfile="orca2_nicas_grids_local_000001-000001.nc"
# testfile="orca2_nicas_grids_local_000002-000001.nc"
testfile="orca2_nicas_grids_local_000002-000002.nc"

fname=testdata+"/"+testfile
print("Opening ",fname)

f = Dataset(fname, "r", format="NETCDF4")

print(f.groups)

for group in f.groups:
   print("group ",group)
   for group2 in f.groups[group].groups:
       print("group2 ",group2)
   lon = f.groups[group].groups[group2]["lon_sc"][:]
   lat = f.groups[group].groups[group2]["lat_sc"][:]
   
   plt.scatter(lon,lat)
   plt.savefig("bump_grid_"+group+"_"+group2+"_2MPI_2_2.png")

f.close()
    

