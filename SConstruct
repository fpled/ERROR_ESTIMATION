from LMT import *

arch = os.uname()[4]

opt = 'opt' in ARGUMENTS and int(ARGUMENTS['opt'])
debug = 'debug' in ARGUMENTS and int(ARGUMENTS['debug'])
timdavis = 'timdavis' in ARGUMENTS and int(ARGUMENTS['timdavis'])
arch = ( 'arch' in ARGUMENTS and ARGUMENTS['arch'] ) or os.uname()[4]

cppfl = ""
#cppfl += ' -DWITH_CHOLMOD '
cppfl += ' -DWITH_UMFPACK -DWITH_LDL '
#cppfl += ' -DWITH_MUMPS -DWITH_MKL '
cppfl += ' -DPRINT_ALLOC ' * 0 # si on veut surveiller la memoire geree par LMTpp
cppfl += ' -pipe ' * 0 # si on veut accelerer le processus de compilation, a condition qu'il y ait assez de memoire
cppfl += ' -Wall ' * 0 # tous les warning standards que le compilateur sait generer
cppfl += ' -pthread '
cppfl += ' -g3 ' # pour avoir la correspondance binaire (assembleur) -> code C++ (ligne fichier...) pour pouvoir utiliser gdb ou valgrind
cppfl += ' -O1 ' * (opt==1) + ' -O2 ' * (opt==2) + ' -O3 ' * (opt==3) + ' -ffast-math -fexpensive-optimizations ' * (opt!=0) # -ffast-math pour autoriser le compilateur a utiliser la commutation et l'associativite ds les calc arithmetiques
cppfl += ' -fpermissive '
cppfl += ' -DDEBUG ' * debug # flag pour LMTpp pour des tests supplementaires du genre depassement d'indice -> ralentit tout
cppfl += ' -DH5Acreate_vers=1 -DH5Gcreate_vers=1 -DH5Dcreate_vers=1 -DH5Dopen_vers=1 -DH5Gopen_vers=1 ';

env = Environment(
    CC = 'gcc',
    CXX = 'g++',
    CPPPATH = [ '#LMT/include', '/usr/local/include', '/opt/local/include', '/usr/local/opt/lapack/include', '/opt/intel/mkl/include' ],
    CPPFLAGS = cppflags( ['xml2-config'] ) + cppfl,
    LINKFLAGS = linkflags( ['xml2-config'] ) + ' -L/usr/local/lib -L/opt/local/lib -L/usr/local/opt/lapack/lib -L/usr/local/gfortran/lib -L/opt/intel/lib -L/opt/intel/mkl/lib ',
    AS = 'yasm -f Elf32 -g Dwarf2',
    LIBS = [ 'cholmod', 'umfpack', 'blas', 'lapack', 'pthread', 'mkl_core', 'mkl_intel_lp64', 'mkl_sequential', 'iomp5', 'dmumps', 'mumps_common', 'mpi', 'mpi_cxx', 'hdf5' ]
)

make_dep_py(env) # sert a generer le code Metil en C++

env.VariantDir( 'build/LMT', 'LMT/include', duplicate=0 )

libs = SConscript( 'LMT/include/SConscript', exports='env', variant_dir='build/LMT' )

env.Command( ["build/problem_error_estimation/all_in_one.h"], ["formulation.met","LMT/LmtppFormulation.met"], "export METILPATH=../METIL/MET; ../METIL-install/bin/metil formulation.met export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/opt/intel/lib:/opt/intel/mkl/lib" )

env.Program( "main", ["main.cpp"] + libs )
