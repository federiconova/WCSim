
COMMAND="bqsub"
TIME="10:00"
EXECUTABLE="./exe/bin/Linux-g++/WCSim"

#MACFILE="WCSim.mac"
MACFILE="WCSim_x16TriangularTiles_600MeVmuon.mac"

echo ${COMMAND} -c ${TIME} ${EXECUTABLE} ${MACFILE}
${COMMAND} -c ${TIME} ${EXECUTABLE} ${MACFILE} > log_${MACFILE}.log


