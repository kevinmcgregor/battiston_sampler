#!/bin/bash
#PBS -N bsamp
#PBS -o log/
#PBS -e log/
#PBS -t 1-100
#PBS -l mem=4G
#PBS -l vmem=4G
#PBS -l walltime=6:00:00
#PBS -l nodes=1:ppn=4
#PBS -m ae
#PBS -M kevin.mcgregor@mail.mcgill.ca

cd $PBS_O_WORKDIR
mkdir -p log/

# Argument filename
ARG=/mnt/GREENWOOD_BACKUP/home/kevin.mcgregor/research/pitman_yor/battiston_sampler/sim/arg_files/af1.txt

echo Started: `date`
Rscript sim_b_sampler.R $PBS_ARRAYID $ARG
echo Ended: `date`
