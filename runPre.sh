#!/bin/sh
#PBS -l nodes=1:ppn=1
#PBS -l mem=2gb
#PBS -l walltime=100:00:00

cd $PBS_O_WORKDIR

#perl table2inf.pl --table ../input/meerkat.filtered.table --project Meerkat.filtered
#perl table2inf.pl --table ../input/meerkat.filtered.TEass.table --project Meerkat.filtered.TEass

#perl table2inf.pl --table ../input/SVlargerTH100kb/HEG4.meerkat.table --project Meerkat.100kb
#perl table2inf.pl --table ../input/SVlargerTH100kb/HEG4.bd.table --project Bd.100kb
#perl table2inf.pl --table ../input/SVlargerTH100kb/HEG4.pindel.filter.table --project Pindel.100kb

perl table2inf.pl --table ../input/SVlargerTH100kb/100kb.merge.table --project 100kb.merge
perl table2inf.pl --table ../input/SVlargerTH100kb/10_100kb.merge.table --project 10_100kb.merge
perl table2inf.pl --table ../input/SVlargerTH100kb/1_10kb.merge.table --project 1_10kb.merge
perl table2inf.pl --table ../input/SVlargerTH100kb/1kb.merge.table --project 1kb.merge

echo "Done"


