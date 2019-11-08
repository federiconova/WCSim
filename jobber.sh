
COMMAND="bqsub"
TIME="12:00"
EXECUTABLE="./simulate_muons.sh"

energy_min=150

energy_max=5150

energy_step=100

for((en=${energy_min}; en<${energy_max}; en+=${energy_step})); do

    echo ${COMMAND} -c ${TIME} ${EXECUTABLE} ${en} 
    ${COMMAND} -c ${TIME} ${EXECUTABLE} ${en}  > log_${en}.log

done


