
n_particles=100000
type_particle="mu-"
script="MakeKin_RADH2O.py"

energy_min=150

energy_max=5150

energy_step=100

for((en=${energy_min}; en<${energy_max}; en+=${energy_step})); do

    python ${script} -n ${n_particles} -t ${type_particle} -e ${en} -v random
    kinfile="muminus_${en}MeV_random_4pi_000.kin"
    macfile="muminus_${en}MeV_random_4pi_000.mac"
    rootfile="wcsim_muminus_${en}MeV_random_4pi_000.root"
    cp WCSim_x16TriangularTiles_base.mac ${macfile}
    sed -i "s/.*vecfile.*/\/mygen\/vecfile ${kinfile}/g" ${macfile}
    sed -i "s/.*RootFile.*/\/WCSimIO\/RootFile ${rootfile}/g" ${macfile}
    sed -i "s/.*beamOn.*/\/run\/beamOn ${n_particles}/g" ${macfile}


done
