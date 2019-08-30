
COMMAND="bqsub"
TIME="10:00"
EXECUTABLE="./exe/bin/Linux-g++/WCSim"

#MACFILE="WCSim_600MeVmuon.mac"
#MACFILE="WCSim_x1TriangularTiles_600MeVmuon.mac"
#MACFILE="WCSim_600MeVmuon_darknoise.mac"
#MACFILE="WCSim_x1TriangularTiles_600MeVmuon_darknoise.mac"
#MACFILE="WCSim_FastStars_600MeVmuon.mac"
#MACFILE="WCSim_x16TriangularTiles_600MeVmuon.mac"
#MACFILE="WCSim_x16TriangularTiles_600MeVmuon_darknoise.mac"

#MACFILE="WCSim_x1SquarePlate_600MeVmuon.mac"
#MACFILE="WCSim_x1TriangularTiles_600MeVmuon.mac"

MACFILE="WCSim_OD_600MeVmuon.mac"
#MACFILE="WCSim_OD_SquarePlate_600MeVmuon.mac"
#MACFILE="WCSim_OD_Triangle_600MeVmuon.mac"
#MACFILE="WCSim_OD_ID_Triangle_600MeVmuon.mac"
#MACFILE="WCSim_OD_ID_SquarePlate_600MeVmuon.mac"

echo ${COMMAND} -c ${TIME} ${EXECUTABLE} ${MACFILE}
${COMMAND} -c ${TIME} ${EXECUTABLE} ${MACFILE} > log_${MACFILE}.log


