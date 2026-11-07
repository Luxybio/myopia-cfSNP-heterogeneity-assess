#!/bin/sh
# Replace <LOCAL_DATA_ROOT>/<SERVER_DATA_ROOT> with your own paths.

#PBS -N intron1-GWAS-left
#PBS -l walltime=24:00:00   # TODO: set your cluster walltime
#PBS -q workq   # TODO: set your cluster queue name
#PBS -l nodes=1
#PBS -l ncpus=2
#PBS -j oe

cd $PBS_O_WORKDIR

date
bowtie2 -x <SERVER_DATA_ROOT>/data/EAS-SNP-DNAelement-relocation/intron-phase1-base/intron-phase1 -f <SERVER_DATA_ROOT>/data/dbSNPGRCh38_myopia_flanking_30bp/EAS-SNP-flanking/GWAS-EAS-SNP-left.fa --n-ceil C,3 --np 0 --end-to-end -a --score-min C,0 --quiet --no-unal --no-head --no-sq -S <SERVER_DATA_ROOT>/data/Alignment/intron-align/EAS-GWAS/EAS-intron-phase1-GWAS-left.sam
date
