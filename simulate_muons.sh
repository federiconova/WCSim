
en=$1

source ../../hk-hyperk/setup.sh
export WCSIM_BASE_DIR=/opt/ppd/t2k/users/federico_nova/HyperK/WCSim-WLSplate-original/WCSim
export G4WORKDIR=${WCSIM_BASE_DIR}/exe


    vecfile="muminus_${en}MeV_random_4pi_000.kin"
    macfile="muminus_${en}MeV_random_4pi_000.mac"
    rootfile="wcsim_muminus_${en}MeV_random_4pi_000.root"
    cp allmacs/${vecfile} .
    cp allmacs/${macfile} .
    ./exe/bin/Linux-g++/WCSim ${macfile}
    rm ${vecfile}
    rm ${macfile}

    mv ${rootfile} sample-root-scripts/outputs/muminus_random_x16TriangularTiles/

