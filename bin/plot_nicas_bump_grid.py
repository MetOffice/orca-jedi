# plot some bump data
# 
# D. J. Lea
# Oct 2022

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np

# original 1 MPI
# testdata="/data/users/sophia.moreton/jedi/mo-bundle-dev/orca-jedi/src/tests/Data"
# # testfile="orca2_nicas_grids_local_000001-000001.nc"
# global_flag = False

# old:
# testdata="/data/users/sophia.moreton/jedi/mo-bundle-dev/orca-jedi/src/tests/Data/100res"
# testfile="orca2_nicas.nc"
# global_flag = True

testdata="/data/users/sophia.moreton/jedi/mo-bundle-dev/orca-jedi/src/tests/Data"
# testfile="orca2_my_bump_files_nicas_1e6_hlength.nc" #8 resolution #looks ok
# testfile="orca2_my_bump_files_nicas_4e7_hlength.nc" #8 resolution #not good
testfile="orca2_my_bump_files_nicas_4e6_hlength.nc" #8 resolution
global_flag = True

fname=testdata+"/"+testfile
print("Opening ",fname)

def find_var(group, names):
    for name in names:
        if name in group.variables:
            return group.variables[name][:]
    for subgroup in group.groups.values():
        value = find_var(subgroup, names)
        if value is not None:
            return value
    return None


if global_flag:
    f = Dataset(fname, "r", format="NETCDF4")
    print(f.groups)

    for group in f.groups:
        print("group ", group)
        for group2 in f.groups[group].groups:
            print("group2 ", group2)

            subgrp = f.groups[group].groups[group2]
            lon = find_var(subgrp, ["lon_c1", "lon_sc", "lon"])
            lat = find_var(subgrp, ["lat_c1", "lat_sc", "lat"])

            if lon is None or lat is None:
                print("  no lon/lat found in", group, group2)
                continue

            plt.figure()
            plt.scatter(lon, lat, s=1)
            plt.savefig(f"bump_grid_{group}_{group2}_2MPI_global_4e6_hlength.png")
            plt.close() 

    f.close()

else:
    f = Dataset(fname, "r", format="NETCDF4")
    print(f.groups)

    for group in f.groups:
        print("group ",group)
        for group2 in f.groups[group].groups:
            print("group2 ",group2)
        lon = f.groups[group].groups[group2]["lon_sc"][:]
        lat = f.groups[group].groups[group2]["lat_sc"][:]
        
        plt.scatter(lon,lat)
        plt.savefig("bump_grid_"+group+"_"+group2+"_2MPI_2_1.png")

    f.close()
    

