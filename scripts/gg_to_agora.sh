#!/bin/bash

#Map gg OTUs to AGORA 16S seqs thoroughly

tail -6904 gg_13_8_99_otus.fasta > gg_13_8_99_otus_41.fasta

for num in {1..41}; do
	#split=$(($num * 10000))
	#echo $split
	#echo gg_13_8_99_otus_$num.fasta
	#head -$split gg_13_8_99_otus.fasta | tail -10000 > gg_13_8_99_otus_$num.fasta
	echo "~/bin/vsearch --usearch_global gg_13_8_99_otus_$num.fasta --db ../blastDB/agora_NCBI_16S.udb --id 0.97 --strand both --maxrejects 3636 --blast6out gg_13_8_99_to_agora_97id_$num.txt" | qsub -l mfree=30g -l h_rt=96:0:0
done


