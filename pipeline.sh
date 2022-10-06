#!/bin/bash

output=$1
ages=$2
path_to_relate=$3
map=$4

#Simulate
echo Simulating
python simulate.py -out ${output} -ages ${ages}

num_seqs=10
echo $num_seqs
echo Running Relate
#Run Relate
for ((i=0; i<num_seqs; i++)); do
	echo $i
	./relate.sh $output $i ${path_to_relate} ${map} &
done
wait

#Train model 
echo Training model
python training.py -ts ${output} -out ${output} -nn 5 -poplab ${output}.poplabels
#Testing
echo Testing model
python testing_bygroup.py -ts ${output} -out ${output} -model model_${output}.h5 -nn 5 -poplab ${output}.poplabels

