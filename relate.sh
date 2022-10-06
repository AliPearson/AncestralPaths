output=$1
i=$2

/home/ap885/rds/rds-durbin-group-8b3VcZwY7rY/projects/human/ancient/relate/relate/bin/RelateFileFormats --mode ConvertFromVcf --haps ${output}_${i}.haps --sample ${output}_${i}.sample -i ${output}_${i}
/home/ap885/rds/rds-durbin-group-8b3VcZwY7rY/projects/human/ancient/relate/relate/bin/RelateFileFormats --mode GenerateSNPAnnotations --haps ${output}_${i}.haps --sample ${output}_${i}.sample --poplabels ${output}.poplabels -o ${output}_${i}
/home/ap885/rds/rds-durbin-group-8b3VcZwY7rY/projects/human/ancient/relate/relate/scripts/RelateParallel/RelateParallel.sh -m 1.25e-8 --threads 8 -N 30000 --sample_ages ${output}.ages --haps ${output}_${i}.haps --sample ${output}_${i}.sample --map genetic_map_simulated.txt -o ${output}_relate_${i} --annot ${output}_${i}.annot
/home/ap885/rds/rds-durbin-group-8b3VcZwY7rY/projects/human/ancient/relate/relate/scripts/EstimatePopulationSize/EstimatePopulationSize.sh -i ${output}_relate_${i} -m 1.25e-8 -o ${output}_relate_popsize_${i} --threads 8 --num_iter 5 --poplabels ${output}.poplabels
/home/ap885/rds/rds-durbin-group-8b3VcZwY7rY/projects/human/ancient/relate/relate/bin/RelateFileFormats --mode ConvertToTreeSequence -i ${output}_relate_popsize_${i} -o ${output}_relate_popsize_${i}

