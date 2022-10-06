output=$1
i=$2
path_to_relate=$3

${path_to_relate}/bin/RelateFileFormats --mode ConvertFromVcf --haps ${output}_${i}.haps --sample ${output}_${i}.sample -i ${output}_${i}
${path_to_relate}/bin/RelateFileFormats --mode GenerateSNPAnnotations --haps ${output}_${i}.haps --sample ${output}_${i}.sample --poplabels ${output}.poplabels -o ${output}_${i}
${path_to_relate}/scripts/RelateParallel/RelateParallel.sh -m 1.25e-8 --threads 8 -N 30000 --sample_ages ${output}.ages --haps ${output}_${i}.haps --sample ${output}_${i}.sample --map genetic_map_simulated.txt -o ${output}_relate_${i} --annot ${output}_${i}.annot
${path_to_relate}/scripts/EstimatePopulationSize/EstimatePopulationSize.sh -i ${output}_relate_${i} -m 1.25e-8 -o ${output}_relate_popsize_${i} --threads 8 --num_iter 5 --poplabels ${output}.poplabels
${path_to_relate}/bin/RelateFileFormats --mode ConvertToTreeSequence -i ${output}_relate_popsize_${i} -o ${output}_relate_popsize_${i}

