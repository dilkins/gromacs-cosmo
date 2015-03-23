#!/bin/bash --login

mpi="false"


if [ ${mpi} == "true" ] 
then

    export MPICCDIR=/usr/bin/

    export CXX=${MPICCDIR}/mpicxx
    export CC=${MPICCDIR}/mpicc

    cmake ../ \
         -DGMX_BUILD_OWN_FFTW=ON \
         -DREGRESSIONTEST_DOWNLOAD=OFF \
         -DREGRESSIONTEST_PATH=/home/tocci/workspace/codes/gromacs-5.0.4/build/regressiontests-5.0.4 \
         -DCMAKE_INSTALL_PREFIX=/home/tocci/workspace/codes/gromacs-5.0.4/installation \
         -DGMX_MPI=ON \
         -DCMAKE_CXX_COMPILER=${CXX} \
         -DCMAKE_C_COMPILER=${CC} \
         -DBUILD_SHARED_LIBS=OFF 


#    make -j 4 mdrun
    make
    sudo make install-mdrun_mpi
#    make install-mdrun_mpi




else

      export MPICCDIR=/usr/bin/

    export CXX=${MPICCDIR}/mpicxx
    export CC=${MPICCDIR}/mpicc

    cmake ../ \
         -DGMX_BUILD_OWN_FFTW=ON \
         -DREGRESSIONTEST_DOWNLOAD=OFF \
         -DREGRESSIONTEST_PATH=/home/tocci/workspace/codes/gromacs-5.0.4/build/regressiontests-5.0.4 \
         -DCMAKE_INSTALL_PREFIX=/home/tocci/workspace/codes/gromacs-5.0.4/installation \
         -DGMX_MPI=OFF \
         -DCMAKE_CXX_COMPILER=${CXX} \
         -DCMAKE_C_COMPILER=${CC} \
         -DBUILD_SHARED_LIBS=OFF


#    make -j 4 mdrun
    make 
  
    make install 

fi
