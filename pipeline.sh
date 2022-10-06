#!/bin/bash

output=$1

#Simulate
echo Simulating
python simulate.py -out ${output} -ages mesoneo_ages_pop.txt

num_seqs=10
echo $num_seqs
echo Running Relate
#Run Relate
for ((i=0; i<num_seqs; i++)); do
	echo $i
	./relate.sh $output $i &
done
wait

#Train model 
echo Training model
python training.py -ts ${output} -out ${output} -nn 5 -poplab ${output}.poplabels
#Testing
echo Testing model
python testing_bygroup.py -ts ${output} -out ${output} -model model_${output}.h5 -nn 5 -poplab ${output}.poplabels

