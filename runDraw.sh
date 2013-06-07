#!/bin/sh
#PBS -l nodes=1:ppn=1
#PBS -l mem=2gb
#PBS -l walltime=100:00:00

cd $PBS_O_WORKDIR

#perl drawRegionPipe.pl --table test.inf --flanking 0 --ACT --project HEG4_SV_test
#perl drawRegionPipe.pl --table Meerkat.filtered.test.inf  --flanking 50000 --ACT --reorder --project HEG4_SV_Meerkat
#perl drawRegionPipe.pl --table Meerkat.filtered.TEass.draw.inf  --flanking 200000 --ACT --reorder --project HEG4_SV_Meerkat_TEass_anno
#perl drawRegionPipe.pl --table Meerkat.100kb.draw.inf  --flanking 200000 --ACT --reorder --project HEG4_SV_Meerkat_100kb
#perl drawRegionPipe.pl --table Bd.100kb.draw.inf  --flanking 200000 --ACT --reorder --project HEG4_SV_Bd_100kb

#10_100kb.merge.inf
#perl drawRegionPipe.pl --table 100kb.merge.draw.inf  --flanking 200000 --ACT --reorder --project HEG4_SV_100kb.merge
#perl drawRegionPipe.pl --table 10_100kb.merge.draw.inf  --flanking 200000 --ACT --reorder --project HEG4_SV_10_100kb.merge
#perl drawRegionPipe.pl --table 1_10kb.merge.draw.inf  --flanking 200000 --ACT --reorder --project HEG4_SV_1_10kb.merge
#perl drawRegionPipe.pl --table 1kb.merge.draw.inf  --flanking 200000 --ACT --reorder --project HEG4_SV_1kb.merge
perl drawRegionPipe.pl --table 1kb.merge.draw.inf1  --flanking 200000 --ACT --reorder --project HEG4_SV_1kb.merge
