#!/bin/bash

output=$1
ages=$2
path_to_relate=$3
mkdir "${output}"
cd ${output}
mkdir relate_files
mkdir simulated_files
cp ../genetic_map_simulated.txt genetic_map_simulated.txt
cd ../

#Simulate
echo Simulating
python simulate.py -out ${output}/${output} -ages ${ages} > simulate.log

num_seqs=10
echo $num_seqs
echo Running Relate
cd $output
#Run Relate
for ((i=0; i<num_seqs; i++)); do
	echo $i
	.././relate.sh ${output} $i ${path_to_relate} > relate.log &
done
wait

#Train model 
echo Training model
python ../training.py -ts ${output} -out ${output} -nn 5 -poplab ${output}.poplabels
#Testing
echo Testing model
python ../testing_bygroup.py -ts ${output} -out ${output} -model model_${output}.h5 -nn 5 -poplab ${output}.poplabels

### Clean up
mv *.sample relate_files
mv *.haps relate_files
mv *.annot relate_files
mv ${output}_relate* relate_files
mv ${output}*.vcf simulated_files
mv ${output}*.trees simulated_files
mv ${output}.ages simulated_files
mv ${output}.poplabels simulated_files
mv genetic_map* simulated_files
mv relate.log relate_files
mv simulate.log simulated_files

