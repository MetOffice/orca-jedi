# (C) British Crown Copyright 2026 Met Office
# This software is licensed under the terms of the Apache Licence Version 2.0 which can be obtained at
# http://www.apache.org/licenses/LICENSE-2.0.
#

"""JOPA development container with LLVM Clang

This container provides a full LLVM/Clang toolchain. All dependencies are built
with Clang for ABI consistency, avoiding mixed GCC/Clang builds.

Usage:
hpccm --recipe THIS-FILE [options] > DEST
"""


def github_url(repo, vn):
    return f'https://github.com/{repo}/archive/refs/tags/{vn}.tar.gz'


def gitlab_url(repo, vn):
    name = repo.rsplit('/', 1)[1]
    return f'https://gitlab.com/{repo}/-/archive/{vn}/{name}-{vn}.tar.gz'


# get versions via --userarg options
# build hdf5/netcdf with zstd (libzstd from epel)
atlas_orca_vn = USERARG.get('atlas_orca_vn', '0.4.3')
atlas_vn = USERARG.get('atlas_vn', '0.44.1')
blitz_vn = USERARG.get('blitz_vn', '1.0.2')
boost_vn = USERARG.get('boost_vn', '1.88.0')
bufr_query_vn = USERARG.get('bufr_query_vn', '0.0.4')
cmake_vn = USERARG.get('cmake_vn', '3.26.5')
ecbuild_vn = USERARG.get('ecbuild_vn', '3.11.0')
eccodes_vn = USERARG.get('eccodes_vn', '2.30.1')  # requires AEC (libaec-devel) by default
eckit_vn = USERARG.get('eckit_vn', '1.29.3')
ectrans_vn = USERARG.get('ectrans_vn', '1.6.2')
fckit_vn = USERARG.get('fckit_vn', '0.13.4')
fiat_vn = USERARG.get('fiat_vn', '1.5.1')
fparser_vn = USERARG.get('fparser_vn', '0.2.0')
gsl_lite_vn = USERARG.get('gsl_lite_vn', '0.42.0')
gsw_fortran_vn = USERARG.get('gsw_fortran_vn', '3.08')
hdf5_vn = USERARG.get('hdf5_vn', '1.14.6')
json_schema_validator_vn = USERARG.get('json_schema_validator_vn', '2.3.0')
json_vn = USERARG.get('json_vn', '3.11.3')
lapack_vn = USERARG.get('lapack_vn', '3.12.1')
nccmp_vn = USERARG.get('nccmp_vn', '1.9.1.0')
nceplibs_bufr_vn = USERARG.get('nceplibs_bufr_vn', '12.2.0')
netcdf_vn = USERARG.get('netcdf_vn', '4.9.2')
netcdfcxx_vn = USERARG.get('netcdfcxx_vn', '4.3.1')
netcdfftn_vn = USERARG.get('netcdfftn_vn', '4.6.1')
netcdf4python_vn = USERARG.get('netcdf4python_vn', '1.6.5')
numpy_vn = USERARG.get('numpy_vn', '1.26.4')
odc_vn = USERARG.get('odc_vn', '1.6.1')
openmpi_vn = USERARG.get('openmpi_vn', '4.1.5')
pycodestyle_vn = USERARG.get('pycodestyle_vn', '2.10')
qhull_vn = USERARG.get('qhull_vn', '8.0.2')  # no qhull-devel rpm
udunits_vn = USERARG.get('udunits_vn', '2.2.28')
yaxt_vn = USERARG.get('yaxt_vn', '528-0.10.0')  # URL has a number and a version

COMMON_PACKAGES = [
    'bison',
    'binutils',  # System linker and binary tools
    'bzip2',
    'clang-tools-extra',
    'eigen3-devel',
    'expat-devel',
    'file',
    'flex',
    'fftw-devel',
    'gcc',  # System GCC for linker and runtime support
    'gcc-c++',  # System G++ for complete C++ toolchain
    'gcc-gfortran',  # For Fortran support (no flang available)
    'git',
    'git-lfs',
    'gmp-devel',
    'gnupg2',
    'graphviz',
    'jq',
    'lcov',
    'less',
    'libaec-devel',
    'libcurl-devel',
    'libX11-devel',
    'libxml2-devel',
    'libzstd-devel',
    'lz4-devel',
    'mpfr-devel',
    'ncurses-devel',
    'ninja-build',
    'openssh-server',
    'openssl-devel',
    'patch',
    'pkgconfig',
    'pybind11-devel',
    'python3-devel',
    'python3-pip',
    'python3-pytest',
    'python3-pyyaml',
    'python3-scipy',
    'rsync',
    'time',
    'unzip',
    'vim-minimal',
    'wget',
    'xz',
    'zlib-devel',
    'zstd',
]

Stage0 += baseimage(image='almalinux:9', _as='build', _distro='rhel')
Stage0 += shell(commands=[
    'dnf install -y \'dnf-command(config-manager)\'',
    'dnf config-manager -y --set-enabled crb',
])

# Full LLVM/Clang stack for AlmaLinux 9
# Using libstdc++ (GCC's standard library) since libc++ packages are not readily available
# Note: System GCC (gcc, gcc-c++, gfortran) provides linker and runtime needed by Clang
# Key: ALL dependencies are built with Clang for ABI consistency
LLVM_PACKAGES = [
    'clang',
    'clang-tools-extra',
    'libomp-devel',
    'llvm',
    'compiler-rt',  # LLVM compiler runtime
    'libstdc++-devel',  # GCC standard library for C++
]

Stage0 += packages(epel=True, ospackages=COMMON_PACKAGES + LLVM_PACKAGES)

# Set up compiler paths BEFORE building anything
# System LLVM provides: clang, clang++
# System GCC provides: gcc, g++, gfortran (gcc-c++ package needed for complete toolchain)
# Set compilers to clang for the build stage
# Using libstdc++ (default on AlmaLinux) - all dependencies built with Clang for consistency
# Note: gcc-c++ provides linker (ld) and standard library support that Clang requires
# CRITICAL: Force Clang to use system linker with -fuse-ld=/usr/bin/ld (not gcc-toolset paths)
# Only set for C/C++ - gfortran doesn't understand this flag
# Explicitly include /usr/local/lib for libraries built by this container
Stage0 += environment(variables={
    'CC': 'clang',
    'CXX': 'clang++',
    'FC': 'gfortran',
    'CFLAGS': '-fuse-ld=/usr/bin/ld',
    'CXXFLAGS': '-fuse-ld=/usr/bin/ld',
    'PATH': '/usr/local/bin:/usr/bin:$PATH',
    'LD_LIBRARY_PATH': '/usr/local/lib:/usr/lib64:$LD_LIBRARY_PATH',
})

Stage0 += cmake(eula=True, version=cmake_vn)

Stage0 += generic_cmake(
    prefix='/usr/local',
    url=github_url('Reference-LAPACK/lapack', f'v{lapack_vn}'),
    directory=f'lapack-{lapack_vn}',
    cmake_opts=['-DCMAKE_BUILD_TYPE=Release', '-DBUILD_SHARED_LIBS=ON'],
)

Stage0 += packages(ospackages=['bzip2', 'bzip2-devel', 'tar', 'wget', 'which'])
# Build Boost with GCC (simpler than fighting with Clang's linker issues in bootstrap)
# This is fine because GCC and Clang both use libstdc++ - ABI compatible
Stage0 += shell(commands=[
    f'mkdir -p /var/tmp && wget -q -nc --no-check-certificate -P /var/tmp https://archives.boost.io/release/{boost_vn}/source/boost_1_88_0.tar.bz2',
    'mkdir -p /var/tmp && tar -x -f /var/tmp/boost_1_88_0.tar.bz2 -C /var/tmp -j',
    'cd /var/tmp/boost_1_88_0 && ./bootstrap.sh --prefix=/usr/local --with-libraries=chrono,date_time,filesystem,program_options,regex,serialization,system,thread --with-toolset=gcc',
    'cd /var/tmp/boost_1_88_0 && ./b2 toolset=gcc cxxflags="-std=c++17" -j$(nproc) -q install',
    'rm -rf /var/tmp/boost_1_88_0.tar.bz2 /var/tmp/boost_1_88_0',
])

mpi = openmpi(
    prefix='/usr/local',
    version=openmpi_vn,
    cuda=False,
    infiniband=False,
    configure_opts=['--enable-mpi-fortran', '--enable-mpi-cxx'],
)
Stage0 += mpi

Stage0 += generic_cmake(
    prefix='/usr/local',
    url=github_url('HDFGroup/hdf5', f'hdf5-{hdf5_vn}'),
    directory=f'hdf5-hdf5-{hdf5_vn}',
    cmake_opts=[
        '-DCMAKE_BUILD_TYPE=Release',
        '-DBUILD_SHARED_LIBS=ON',
        '-DHDF5_ENABLE_PARALLEL=ON',
        '-DHDF5_BUILD_FORTRAN=ON',
        '-DHDF5_ENABLE_ZLIB_SUPPORT=ON',
        '-DHDF5_ENABLE_SZIP_SUPPORT=ON',
    ],
    toolchain=mpi.toolchain,
)

Stage0 += environment(variables={'H5DIR': '/usr/local', 'LIBS': '-ldl'})
Stage0 += netcdf(
    version=netcdf_vn,
    version_cxx=netcdfcxx_vn,
    version_fortran=netcdfftn_vn,
    prefix='/usr/local',
    cxx=True,
    fortran=True,
    enable_netcdf_4=True,
    enable_shared=True,
    disable_zstandard_plugin=True,
    toolchain=mpi.toolchain,
)
Stage0 += environment(variables={
    'NETCDF_DIR': '/usr/local',
    'NetCDF_ROOT': '/usr/local'})
Stage0 += generic_cmake(
    prefix='/usr/local',
    url=gitlab_url('remikz/nccmp', nccmp_vn),
    cmake_opts=['-DCMAKE_BUILD_TYPE=Release', '-DBUILD_SHARED_LIBS=ON'],
    toolchain=mpi.toolchain,
)

Stage0 += generic_autotools(
    prefix='/usr/local',
    url=f'https://downloads.unidata.ucar.edu/udunits/{udunits_vn}/udunits-{udunits_vn}.tar.gz',
    configure_opts=['--enable-shared=yes'],
)

Stage0 += generic_cmake(
    prefix='/usr/local',
    url=github_url('gsl-lite/gsl-lite', f'v{gsl_lite_vn}'),
    directory=f'gsl-lite-{gsl_lite_vn}',
    cmake_opts=['-DCMAKE_BUILD_TYPE=Release'],
)

Stage0 += generic_cmake(
    prefix='/usr/local',
    url=github_url('blitzpp/blitz', blitz_vn),
    directory=f'blitz-{blitz_vn}',
    cmake_opts=['-DCMAKE_BUILD_TYPE=Release', '-DBUILD_SHARED_LIBS=ON'],
)

Stage0 += generic_cmake(
    prefix='/usr/local',
    url=github_url('nlohmann/json', f'v{json_vn}'),
    directory=f'json-{json_vn}',
    cmake_opts=['-DCMAKE_BUILD_TYPE=Release', '-DJSON_BuildTests=OFF'],
)

Stage0 += generic_cmake(
    prefix='/usr/local',
    url=github_url('qhull/qhull', f'v{qhull_vn}'),
    directory=f'qhull-{qhull_vn}',
    cmake_opts=['-DCMAKE_BUILD_TYPE=Release'],
)

Stage0 += generic_cmake(
    prefix='/usr/local',
    url=github_url('NOAA-EMC/NCEPLIBS-bufr', f'v{nceplibs_bufr_vn}'),
    directory=f'NCEPLIBS-bufr-{nceplibs_bufr_vn}',
    cmake_opts=['-DCMAKE_BUILD_TYPE=Release', '-DBUILD_TESTS=OFF'],
)

Stage0 += generic_cmake(
    prefix='/usr/local',
    url=github_url('pboettch/json-schema-validator', json_schema_validator_vn),
    directory=f'json-schema-validator-{json_schema_validator_vn}',
    cmake_opts=[
        '-DCMAKE_BUILD_TYPE=Release',
        '-DCMAKE_POLICY_DEFAULT_CMP0074=NEW',
        '-DBUILD_SHARED_LIBS=ON',
        '-DBUILD_TESTS=OFF',
        '-DBUILD_EXAMPLES=OFF',
    ],
)

for org, name, vn in (
    ('ecmwf', 'ecbuild', ecbuild_vn),
    ('ecmwf', 'eckit', eckit_vn),
    ('ecmwf', 'fckit', fckit_vn),
    ('ecmwf', 'odc', odc_vn),
    ('ecmwf-ifs', 'fiat', fiat_vn),
    ('ecmwf-ifs', 'ectrans', ectrans_vn),
    ('ecmwf', 'atlas', atlas_vn),
    ('ecmwf', 'atlas-orca', atlas_orca_vn),
    ('ecmwf', 'eccodes', eccodes_vn),
):
    Stage0 += generic_cmake(
        prefix='/usr/local',
        url=github_url(f'{org}/{name}', vn),
        directory=f'{name}-{vn}',
        cmake_opts=['-DCMAKE_BUILD_TYPE=Release', '-DMPI=ON', '-DOMP=ON'],
    )

# -- bufr-query
Stage0 += generic_cmake(
    prefix='/usr/local',
    url=github_url('NOAA-EMC/bufr-query', f'v{bufr_query_vn}'),
    directory=f'bufr-query-{bufr_query_vn}',
    cmake_opts=['-DCMAKE_BUILD_TYPE=Release', '-DMPI=ON', '-DOMP=ON'],
    build_environment={'LDFLAGS': '-lnetcdf'},
)

# -- GSW-Fortran
Stage0 += generic_cmake(
    prefix='/usr/local',
    url=github_url('JCSDA-internal/GSW-Fortran', f'v{gsw_fortran_vn}'),
    directory=f'GSW-Fortran-{gsw_fortran_vn}',
    cmake_opts=['-DCMAKE_BUILD_TYPE=Release', '-DMPI=ON', '-DOMP=ON'],
)

yaxt_vns = yaxt_vn.split('-', 1)
Stage0 += generic_autotools(
    prefix='/usr/local',
    url=(
        'https://swprojects.dkrz.de/redmine/attachments/download/' +
        f'{yaxt_vns[0]}/yaxt-{yaxt_vns[1]}.tar.xz'
    ),
    configure_opts=['--with-idxtype=long', '--without-regard-for-quality'],
)

Stage0 += pip(pip='pip3', packages=[
    f"pycodestyle=={pycodestyle_vn}",
    f"numpy=={numpy_vn}",
    f"netcdf4=={netcdf4python_vn}",
])
Stage1 += baseimage(image='almalinux:9', _distro='rhel')
Stage1 += comment('JEDI development image with LLVM Clang and OpenMPI')
Stage1 += label(metadata={
    'Maintainer': 'darth@metoffice.gov.uk',
    'Species': 'JOPA',
    'Version': 'v0.2'})
Stage1 += shell(commands=[
    'dnf install -y \'dnf-command(config-manager)\'',
    'dnf config-manager -y --set-enabled crb',
])
Stage1 += packages(epel=True, ospackages=COMMON_PACKAGES + LLVM_PACKAGES)
Stage1 += pip(pip='pip3', packages=[
    'cpplint',
])
Stage1 += copy(_from='build', src='/usr/local', dest='/usr/local')
Stage1 += shell(commands=['ln -sfT python3 /usr/bin/python'])
Stage1 += environment(variables={
    'PATH': '/usr/local/bin:/usr/bin:$PATH',
    'LD_LIBRARY_PATH': '/usr/local/lib:/usr/lib64:$LD_LIBRARY_PATH',
    'VALIDATE_PARAMETERS': '1',
})
Stage1 += workdir(directory='/var/tmp')
