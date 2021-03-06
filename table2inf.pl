#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"table:s","project:s","help");


my $help=<<USAGE;
perl $0 --table
--table: position of SV, the second position is the donor of insertion if exists.
Chr1	2575898	2576036
Chr1	2577382	2577394	Chr5	6721766	6721996
project.inf is the file used for draw ACT by drawRegionPipe.pl, the order will be used to order the lines in figure.
OS	Chr1	2575898 2576036	HEG4	scafold1	2575898 2576036	OG	1	195102	417841
OS      Chr1    2575898 2576036 HEG4    scafold1        2575898 2576036	OG      1       195102  417841	OS	Chr5    6721766 6721996
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

my $os="/rhome/cjinfeng/HEG4_cjinfeng/Variations/SV2ACT/input/MSU_r7.fa.RepeatMasker.masked";
my $osgff="/rhome/cjinfeng/HEG4_cjinfeng/Variations/SV2ACT/input/MSU7.collinear.gene.gff";
my $og="/rhome/cjinfeng/HEG4_cjinfeng/Variations/SV2ACT/input/OGL.fa.RepeatMasker.masked";
my $heg4="/rhome/cjinfeng/HEG4_cjinfeng/Variations/SV2ACT/input/HEG4_RAW.refassist.fa.RepeatMasker.masked";


$opt{project} ||= "Meerkat";

inf($opt{table},$os,$osgff,$og,$heg4);

#my $refseq=getfastaseq($os);

sub inf
{
my ($file,$os,$osgff,$og,$heg4)=@_;
my @table;
my $refseq=getfastaseq($os);
my $refgff=parseGFF($osgff);
open INF, ">$opt{project}.inf" or die "$!";
open IN1, "$file" or die "$!";
while(<IN1>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    ## use collinear genes position as blast seed
    if (exists $refseq->{$unit[0]}){
       my $h1="upstream";
       my ($s1,$len1)=getposition($unit[1],$refgff->{$unit[0]},0);##getposition(start of SV, @gff->[gene start,gene end],0 to get upstream gene)
       my $seq1=substr($refseq->{$unit[0]},$s1,$len1);
       my $h2="downstream"; 
       my ($s2,$len2)=getposition($unit[2],$refgff->{$unit[0]},1);
       my $seq2=substr($refseq->{$unit[0]},$s2,$len2);
       open OUT, ">$opt{project}.query.fa" or die "$!";
       print OUT ">$h1\n$seq1\n>$h2\n$seq2\n";
       close OUT;
       `blastall -p blastn -U -i $opt{project}.query.fa -d $og -o OGL.blast -e 1e-5 -m 8`;
       `blastall -p blastn -U -i $opt{project}.query.fa -d $heg4 -o HEG4.blast -e 1e-5 -m 8`;
       my $oginf=target("OGL.blast","OG");
       my $heg4inf=target("HEG4.blast","HEG4");
       if ($heg4inf eq "NA"){
          $seq1=substr($refseq->{$unit[0]},$unit[1]-20000,20000);
          $seq2=substr($refseq->{$unit[0]},$unit[2],20000);
          open OUT, ">$opt{project}.query.fa" or die "$!";
               print OUT ">$h1\n$seq1\n>$h2\n$seq2\n";
          close OUT;
          `blastall -p blastn -U -i $opt{project}.query.fa -d $heg4 -o HEG4.blast -e 1e-5 -m 8`;
          $heg4inf=target("HEG4.blast","HEG4");
       }
       if ($og4inf eq "NA"){
          $seq1=substr($refseq->{$unit[0]},$unit[1]-20000,20000);
          $seq2=substr($refseq->{$unit[0]},$unit[2],20000);
          open OUT, ">$opt{project}.query.fa" or die "$!";
               print OUT ">$h1\n$seq1\n>$h2\n$seq2\n";
          close OUT;
          `blastall -p blastn -U -i $opt{project}.query.fa -d $og -o OGL.blast -e 1e-5 -m 8`;
          $heg4inf=target("OGL.blast","OG");
       }
       print "$oginf\t$heg4inf\n";
       print INF "OS\t$unit[0]\t$unit[1]\t$unit[2]";
       print INF "\t$oginf" if ($oginf ne "NA");
       print INF "\t$heg4inf" if ($heg4inf ne "NA");
       print INF "\tOS\t$unit[3]\t$unit[4]\t$unit[5]" if (@unit > 3);
       print INF "\n";
       #`rm $opt{project}.query.fa OGL.blast HEG4.blast`;
    }
}
close IN1;
close INF;
}

#upstream	1	99.61	1020	4	0	1	1020	41366	42385	0.0	1990
#downstream	1	99.17	1330	6	2	1	1325	57997	59326	0.0	2456
sub target
{
my ($blast,$title)=@_;
my $target;
my %hash;
print "Find position in $title:\n";
open IN, "$blast" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    my $start = $unit[8] > $unit[9] ? $unit[9] : $unit[8];
    my $end   = $unit[8] > $unit[9] ? $unit[8] : $unit[9];
    if ($unit[0] eq "upstream" and ! exists $hash{"upstream"}){
       print "Upstream\t$unit[1]\t$start\t$end\n";
       @{$hash{"upstream"}}=($unit[1],$start,$end);
    }
    if ($unit[0] eq "downstream" and ! exists $hash{"downstream"}){
       print "Downstream\t$unit[1]\t$start\t$end\n";
       @{$hash{"downstream"}}=($unit[1],$start,$end);
    }
}
close IN;
   print "start and end chromosome: $hash{\"upstream\"}->[0]\t$hash{\"downstream\"}->[0]\n";
if ($hash{"upstream"}->[0] eq $hash{"downstream"}->[0]){
   $target="$title\t$hash{\"upstream\"}->[0]\t$hash{\"upstream\"}->[1]\t$hash{\"downstream\"}->[2]"; ## OS	Chr1	start	end
}else{
   print "Hit on different chromosome: $hash{\"upstream\"}->[0]\t$hash{\"upstream\"}->[1]\t$hash{\"upstream\"}->[2]\t$hash{\"downstream\"}->[0]\t$hash{\"downstream\"}->[1]\t$hash{\"downstream\"}->[2]\n";
   
   $target="NA";
}
return $target;
}


sub getposition
{
my ($pos,$gff,$flag)=@_;
my @sort=$flag == 1 ? sort {$a->[1] <=> $b->[1]} @$gff : sort {$b->[1] <=> $a->[1]} @$gff;
print "Position: $pos\n";
for (my $i=0;$i<@sort;$i++){
    #print "$sort[$i]->[0]\t$sort[$i]->[1]\t$sort[$i]->[2]\n";
    if ($flag == 0){
       if ($sort[$i]->[0] < $pos and $sort[$i]->[1] < $pos){ ## gene start and end 
          print "Find: $sort[$i]->[0]\t$sort[$i]->[1]\t$sort[$i]->[2]\n";
          return ($sort[$i]->[0],$sort[$i]->[1]-$sort[$i]->[0]+1);
       }
    }else{
       if ($sort[$i]->[0] > $pos and $sort[$i]->[1] > $pos){ ## gene start and end
          print "Find: $sort[$i]->[0]\t$sort[$i]->[1]\t$sort[$i]->[2]\n";
          return ($sort[$i]->[0],$sort[$i]->[1]-$sort[$i]->[0]+1);
       }
    }
}
} 



sub getfastaseq
{
$/=">";
my %hash;
my ($file)=@_;
open IN,"$file" or die "$!";
while (<IN>){
    next if (length $_ < 2);
    my @unit=split("\n",$_);
    my $temp=shift @unit;
    my @temp1=split(" ",$temp);
    my $head=$temp1[0];
    my $seq=join("\n",@unit);
    $seq=~s/\>//g;
    $seq=~s/\n//g;
    $seq=~s/\s//g;
    #print "$head\n";
    $hash{$head}=$seq;
}
$/="\n";
return \%hash;
}


sub parseGFF
{
my ($gff)=@_;
my %hash;  ##hash to store every record by key of Seq_id
my $seq;   ##Scaffold
my $id;    ##ID for element
open IN, "$gff" or die "$!";
while (<IN>){
    chomp $_;
    next if ($_=~/^#/);
    my @unit=split("\t",$_);
    if ($unit[2]=~/mRNA/){
        $seq=$unit[0];
        if ($unit[8]=~/ID=(.*?);/ or $unit[8] =~/ID=(.*)/){
            $id=$1;
        }
        push @{$hash{$seq}},[$unit[3],$unit[4],$id];
    }

}
close IN;
return \%hash;
}




