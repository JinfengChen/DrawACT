#!/bin/sh
#PBS -l nodes=1:ppn=1
#PBS -l mem=2gb
#PBS -l walltime=100:00:00

cd $PBS_O_WORKDIR

#perl drawRegionPipe.pl --table test.inf --flanking 0 --ACT --project HEG4_SV_test
#perl drawRegionPipe.pl --table Meerkat.filtered.test.inf  --flanking 50000 --ACT --reorder --project HEG4_SV_Meerkat
#perl drawRegionPipe.pl --table Meerkat.filtered.TEass.draw.inf  --flanking 200000 --ACT --reorder --project HEG4_SV_Meerkat_TEass_anno
perl drawRegionPipe.pl --table Meerkat.100kb.draw.inf  --flanking 200000 --ACT --reorder --project HEG4_SV_Meerkat_100kb
perl drawRegionPipe.pl --table Bd.100kb.draw.inf  --flanking 200000 --ACT --reorder --project HEG4_SV_Bd_100kb
