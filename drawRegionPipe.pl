#!/usr/bin/perl
use Getopt::Long;
use FindBin qw($Bin);

GetOptions (\%opt,"table:s","flanking:s","ACT","reorder","project:s","help");


my $help=<<USAGE;
perl drawRegionPipe.pl --table test.inf --flanking 0 --reorder --ACT --project HEG4_SV_test/
--table: information table used to draw
OS	Chr1	2575898	2576036	OG	1	1905399	1917790	HEG4	Superscaffold1	1979832	1999309
OS	Chr1	2577382	2577394	OG	1	1905399	1917790	HEG4	Superscaffold1	1979832	1999309	OS	Chr5	6721766	6721996
--flanking: if get N of flanking sequence from both side
--reorder: if reorder the sequence order in the table. We put OG at first and do not move others.
--ACT: if draw ACT figure, if not set just do file convert
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

$opt{flanking} ||=0;
$opt{project} ||="HEG4_SV_Meerkat";
`mkdir $opt{project}`;
readtable($opt{table});


#############################################
#OS	chr01   9408106 9415703 OG	chr02   3731742 3739474
sub readtable
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/ or $_=~/^#/);
    my @unit=split("\t",$_);
    my $title;
    my @header;
    for(my $i=0;$i<@unit;$i+=4){
        my $prefix=$unit[$i];
        my $chr   =$unit[$i+1];
        $chr=~s/\_//;
        my ($s,$e)=sort {$a <=> $b} ($unit[$i+2],$unit[$i+3]);
        $s = $s-$opt{flanking} > 0 ? $s-$opt{flanking} : 0;
        $e = $e+$opt{flanking};
        #$title=$prefix."_".$chr."_".$s."_".$e if ($i == 0);
        my $head=$prefix."_".$chr."_".$s."_".$e;
        push @header, $head;
        $title=$head if ($i == 0);
        `perl $Bin/scripts/getsubdata2genome_sd.pl --genegff ../input/MSU7.gene.anno.gff --tegff ../input/MSU_r7.fa.RepeatMasker.out.gff --fasta ../input/MSU_r7.fa --refhead $head` if ($prefix=~/OS/i);
        `perl $Bin/scripts/getsubdata2genome_sd.pl --genegff ../input/OGL.fgenesh.clean.gff --tegff ../input/OGL.fa.RepeatMasker.out.gff --fasta ../input/OGL.fa --refhead $head` if ($prefix=~/OG/i);
        `perl $Bin/scripts/getsubdata2genome_sd.pl --genegff ../input/HEG4_REF.fgenesh.clean.gff --tegff ../input/HEG4_RAW.refassist.fa.RepeatMasker.out.gff --fasta ../input/HEG4_RAW.refassist.fa --refhead $head` if ($prefix=~/HEG4/i);
        if ($opt{ACT}){
           `perl $Bin/act/GFF2embl_anno.pl -gff $head.gene.gff -embl $head.gene.embl -fasta $head.fasta`;
           `perl $Bin/act/gffrepeat2embl.pl -repeat $head.te.gff -embl $head.gene.embl -title $head` if (-f "$head.te.gff");
           if (-f "$head.te.gff"){
                  `mv $head.merge $head.embl`;
           }else{
                  `mv $head.gene.embl $head.embl`;
           }
        }
    }
    if ($opt{reorder}){
        for (my $i=0;$i<@header;$i++){
            if ($header[$i]=~/OG/){
               my $tmp=$header[$i];
               $header[$i]=$header[0];
               $header[0] =$tmp;
            }
        }   
    }
    if ($opt{ACT}){
    my $headers=join(",",@header);
    print "Headers: $headers\n";
    `perl $Bin/act/runblast2seq.pl`;
    `perl $Bin/act/run2act.pl`;
    `perl $Bin/scripts/drawRegionNway.pl --headers $headers --project $title`;
    `rm *.fasta.n* *.blast *.temp *.gene.embl formatdb.log`; 
    }
    print "TEST\t$opt{project}\t$title\n";
    `mkdir $opt{project}/$title`;
    `mv *.gff $opt{project}/$title`;
    `mv *.svg $opt{project}/$title`;
    `mv *.fasta $opt{project}/$title`;
    `mv *.pdf $opt{project}/$title`;
    `mv *.embl $opt{project}/$title`;
    `mv *4ACT $opt{project}/$title`;
    `mv *.shell $opt{project}/$title`;
}
close IN;
return \%hash;
}

