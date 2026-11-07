#!/bin/sh
# Replace <LOCAL_DATA_ROOT>/<SERVER_DATA_ROOT> with your own paths.

#PBS -N intronphase3
#PBS -l walltime=24:00:00   # TODO: set your cluster walltime
#PBS -q workq   # TODO: set your cluster queue name
#PBS -l nodes=1
#PBS -l ncpus=2
#PBS -j oe

cd $PBS_O_WORKDIR

date
bowtie2-build <SERVER_DATA_ROOT>/data/DNAelement-seq/intron_phase_3.fa <SERVER_DATA_ROOT>/data/EAS-SNP-DNAelement-relocation/intron-phase3-base/intron-phase3
date
