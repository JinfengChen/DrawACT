echo "Draw ACT compare for N regions"
perl drawRegionPipe.pl --table test.inf --flanking 0 --ACT --project HEG4_SV_test
qsub runDraw.sh

echo "collinear gene list"
cd ../input
cat colinear.lst | grep -v "^#" | awk '{print "LOC_"$1}' | sort | uniq > colinear.Os.lst
perl /rhome/cjinfeng/software/bin/getidgff.pl -l colinear.Os.lst -g MSU7.gene.gff -o MSU7.collinear.gene.gff

echo "Prepare inf"
qsub runPre.sh
awk '/OG/ && /HEG4/' Meerkat.filtered.inf > Meerkat.filtered.draw.inf
awk '/OG/ && /HEG4/' Meerkat.filtered.TEass.inf > Meerkat.filtered.TEass.draw.inf
echo "Draw"
qsub runDraw.sh

echo "larger SV > 100kb"
qsub runPre.sh
awk '/OG/ || /HEG4/' Meerkat.100kb.inf > Meerkat.100kb.draw.inf
awk '/OG/ || /HEG4/' Bd.100kb.inf > Bd.100kb.draw.inf

