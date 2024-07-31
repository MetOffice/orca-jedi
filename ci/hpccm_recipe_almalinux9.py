# © Copyright 2021 Met Office
# This software is licensed under the terms of the Apache Licence Version 2.0 which can be obtained at
# http://www.apache.org/licenses/LICENSE-2.0.
#

"""JOPA development container (simplified)

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
atlas_orca_vn = USERARG.get('atlas_orca_vn', '0.3.1')
atlas_vn = USERARG.get('atlas_vn', '0.37.0')
blitz_vn = USERARG.get('blitz_vn', '1.0.2')
boost_vn = USERARG.get('boost_vn', '1.85.0')
bufr_query_vn = USERARG.get('bufr_query_vn', 'v0.0.1')
cmake_vn = USERARG.get('cmake_vn', '3.30.0')
ecbuild_vn = USERARG.get('ecbuild_vn', '3.8.5')
eccodes_vn = USERARG.get('eccodes_vn', '2.30.1')  # requires AEC (libaec-devel) by default
eckit_vn = USERARG.get('eckit_vn', '1.26.2')
ectrans_vn = USERARG.get('ectrans_vn', '1.2.0')
fckit_vn = USERARG.get('fckit_vn', '0.13.0')
fiat_vn = USERARG.get('fiat_vn', '1.4.1')
fparser_vn = USERARG.get('fparser_vn', '0.1.4')
gsl_lite_vn = USERARG.get('gsl_lite_vn', '0.41.0')
gsw_fortran_vn = USERARG.get('gsw_fortran_vn', '3.07')
hdf5_vn = USERARG.get('hdf5_vn', 'hdf5-1_14_2')
json_schema_validator_vn = USERARG.get('json_schema_validator_vn', '2.3.0')
json_vn = USERARG.get('json_vn', '3.9.1')
lapack_vn = USERARG.get('lapack_vn', '3.11.0')
nccmp_vn = USERARG.get('nccmp_vn', '1.9.1.0')
nceplibs_bufr_vn = USERARG.get('nceplibs_bufr_vn', '12.0.1')
netcdf_vn = USERARG.get('netcdf_vn', '4.9.2')
netcdfcxx_vn = USERARG.get('netcdfcxx_vn', '4.3.1')
netcdfftn_vn = USERARG.get('netcdfftn_vn', '4.6.1')
odc_vn = USERARG.get('odc_vn', '1.5.2')
openmpi_vn = USERARG.get('openmpi_vn', '4.1.5')
pycodestyle_vn = USERARG.get('pycodestyle_vn', '2.10')
udunits_vn = USERARG.get('udunits_vn', '2.2.28')
yaxt_vn = USERARG.get('yaxt_vn', '528-0.10.0')  # URL has a number and a version

COMMON_PACKAGES = [
    'bison',
    'bzip2',
    'clang-tools-extra',
    'eigen3-devel',
    'expat-devel',
    'flex',
    'gcc',
    'gcc-c++',
    'gcc-gfortran',
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
Stage0 += packages(epel=True, ospackages=COMMON_PACKAGES)

Stage0 += cmake(eula=True, version=cmake_vn)

Stage0 += generic_cmake(
    prefix='/usr/local',
    url=github_url('Reference-LAPACK/lapack', f'v{lapack_vn}'),
    directory=f'lapack-{lapack_vn}',
    cmake_opts=['-DCMAKE_BUILD_TYPE=Release', '-DBUILD_SHARED_LIBS=ON'],
)

Stage0 += boost(
    prefix='/usr/local',
    version=boost_vn,
    b2_opts=['toolset=gcc', 'cxxflags="-std=c++17"'],
    bootstrap_opts=[
        '--with-libraries=chrono,date_time,filesystem,program_options,regex,serialization,system,thread',
    ],
)

mpi = openmpi(
    prefix='/usr/local',
    version=openmpi_vn,
    cuda=False,
    infiniband=False,
    configure_opts=['--enable-mpi-fortran', '--enable-mpi-cxx'],
)
Stage0 += mpi

Stage0 += generic_autotools(
    prefix='/usr/local',
    repository='https://github.com/HDFGroup/hdf5',
    commit=hdf5_vn,
    enable_build_mode='production',
    enable_cxx=True,
    enable_fortran=True,
    enable_parallel=True,
    enable_threadsafe=True,
    enable_unsupported=True,
    toolchain=mpi.toolchain,
    with_szlib='/usr/local',
    with_zlib='/usr/local',
)
# - Post-process h5pcc to avoid f951: Warning:
# Nonexistent include directory '/var/tmp/hdf5-1.14.0/src/H5FDsubfiling'
# sed -i "s|-I${TMPDIR}/hdf5-${hdf5_vn}/src/H5FDsubfiling||g" "/usr/local/bin/h5pcc"
Stage0 += shell(commands=[
    f'sed -i "s|-I/var/tmp/hdf5/src/H5FDsubfiling||g" "/usr/local/bin/h5pcc"'
])
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
    repository='https://github.com/NOAA-EMC/bufr-query',
    commit=bufr_query_vn,
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

Stage1 += baseimage(image='almalinux:9', _distro='rhel')
Stage1 += comment('JEDI development image with GNU and OpenMPI')
Stage1 += label(metadata={
    'Maintainer': 'darth@metoffice.gov.uk',
    'Species': 'NextGen',
    'Version': 'v0.1'})
Stage1 += shell(commands=[
    'dnf install -y \'dnf-command(config-manager)\'',
    'dnf config-manager -y --set-enabled crb',
])
Stage1 += packages(epel=True, ospackages=COMMON_PACKAGES)
Stage1 += pip(pip='pip3', packages=['cpplint', 'yamlprocessor'])
Stage1 += copy(_from='build', src='/usr/local', dest='/usr/local')
Stage1 += shell(commands=['ln -sfT python3 /usr/bin/python'])
Stage1 += environment(variables={
    'PATH': '/usr/local/bin:$PATH',
    'LD_LIBRARY_PATH': '/usr/local/lib64:/usr/local/lib:/usr/lib64:/usr/lib:$LD_LIBRARY_PATH',
})
Stage1 += workdir(directory='/var/tmp')
