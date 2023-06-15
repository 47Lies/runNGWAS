j=0
while [ $j -le 99 ]
do	
       	python3 3_LMM.py --genotype data/genotype --phenotype data/y_50.pheno --nbs results/settings/permutations/neighborhoods/nbs_ --kernel lin --j $j --odir results/llr/permuted --ofile llr_
	((j++))
done
