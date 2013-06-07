#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"genegff:s","tegff:s","fasta:s","qrygenegff:s","qrytegff:s","qryfasta:s","refhead:s","qryhead:s","draw","ACT","project:s","qryname:s","help");


my $help=<<USAGE;
Can get subfregment of gene gff and te gff for one gene or two genes. And draw dotplot of two regions with annotation.
perl $0 --genegff --tegff --fasta --qrygenegff --qrytegff --qryfasta --refhead --qryhead --draw --project Brd

USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

$opt{project} ||="Brd";
$opt{qryname} ||=$opt{project};

my $refhead=$opt{refhead};
print "Refhead:$refhead\n";
& getsubfasta($opt{fasta},$refhead,"+");
& getsubgff3($opt{genegff},$refhead,"+");
& getsubrepeat($opt{tegff},$refhead,"+") if $opt{tegff};

my ($qrychr,$qrystart,$qryend,$qryhead);
if (defined $opt{qryhead}){
$qryhead=$opt{qryhead};
print "Qryhead:$qryhead\n";
& getsubfasta($opt{qryfasta},$qryhead,"+");
& getsubgff3($opt{qrygenegff},$qryhead,"+");
& getsubrepeat($opt{qrytegff},$qryhead,"+") if $opt{qrytegff};
}

if ($opt{ACT}){
`perl ./act/GFF2embl.pl -gff $refhead.gene.gff -embl $refhead.gene.embl -fasta $refhead.fasta`;
`perl ./act/gffrepeat2embl.pl -repeat $refhead.te.gff -embl $refhead.gene.embl -title $refhead` if (-f "$refhead.te.gff");
#`mv $refhead.merge $refhead.embl` if (-f "$refhead.te.gff");
if (-f "$refhead.te.gff"){
   `mv $refhead.merge $refhead.embl`;
}else{
   `mv $refhead.gene.embl $refhead.embl`;
}
`perl ./act/GFF2embl.pl -gff $qryhead.gene.gff -embl $qryhead.gene.embl -fasta $qryhead.fasta`;
`perl ./act/gffrepeat2embl.pl -repeat $qryhead.te.gff -embl $qryhead.gene.embl -title $qryhead` if (-f "$qryhead.te.gff");
#`mv $qryhead.merge $qryhead.embl` if (-f "$qryhead.te.gff");
if (-f "$qryhead.te.gff"){
   `mv $qryhead.merge $qryhead.embl`;
}else{
   `mv $qryhead.gene.embl $qryhead.embl`;
}
`perl ./act/runblast2seq.pl`;
`perl ./act/run2act.pl`;
`rm *.fasta.n* *.blast *.temp *.gene.embl`;
}


if ($opt{draw}){
#`perl drawRegion4heterchromatin.pl --refseq $refhead.fasta --qryseq $qryhead.fasta --refgenegff $refhead.gene.gff --reftegff $refhead.te.gff --qrygenegff $qryhead.gene.gff --qrytegff $qryhead.te.gff --project $opt{refhead}_$opt{qryhead}`;
`perl drawRegion.pl --refseq $refhead.fasta --qryseq $qryhead.fasta --refgenegff $refhead.gene.gff --qrygenegff $qryhead.gene.gff --project $opt{refhead}_$opt{qryhead}`;
my $cmd="perl drawRegion.pl --refseq $refhead.fasta --qryseq $qryhead.fasta --refgenegff $refhead.gene.gff --qrygenegff $qryhead.gene.gff --project $opt{refhead}_$opt{qryhead}";
my $cmdfile="$opt{refhead}_$opt{qryhead}".".shell";
open CMD, ">$cmdfile" or die "$!";
     print CMD "$cmd\n";
close CMD;
#`perl drawRegion.pl --refseq $refhead.fasta --qryseq $qryhead.fasta --refgenegff $refhead.gene.gff --reftegff $refhead.te.gff --qrygenegff $qryhead.gene.gff --qrytegff $qryhead.te.gff --project $opt{refhead}_$opt{qryhead}`;
`rm compare.blast`;
#`perl dotplotRegion.pl --refseq $refhead.fasta --qryseq $qryhead.fasta --refgenegff $refhead.gene.gff --reftegff $refhead.te.gff --qrygenegff $qryhead.gene.gff --qrytegff $qryhead.te.gff --project $opt{refgene}_$opt{qrygene}`;
}

####################
sub parseGFF
{
my ($gff)=@_;
my %hash;  ##hash to store every record by key of Seq_id
my $seq;
my $id;    ##ID for element
my $record;##all line for this record
open IN, "$gff" or die "$!";
while (<IN>){
    chomp $_;
    next if ($_=~/^#/);
    my @unit=split("\t",$_);
    if ($unit[2]=~/mRNA/){
        $seq=$unit[0];
        if ($unit[8]=~/ID=(.*?);/){
            $id=$1;
            $id=$1 if ($id=~/LOC_(.*)/);
        }
        $hash{$id}=[$unit[0],$unit[3],$unit[4],$unit[6]];
    }

}
close IN;
return \%hash;
}

#########
sub getsubrepeat
{
my ($gff,$head,$strand)=@_;
my @temp=split("_",$head);
my $chr=$temp[1];
my $start=$temp[2];
my $end  =$temp[3];
my $len=$end-$start+1;
#print "$chr\t$start\t$end\n";
if ($strand eq "+"){
   open IN1, "$gff" or die "$!";
   open OUT1, ">>$head.te.gff" or die "$!";
   while (<IN1>){
     my @unit=split("\t",$_);
     $unit[0]=~s/\_//;
     if ($unit[0] eq $chr and $unit[3] >= $start and $unit[4] <= $end){
       $unit[3]=$unit[3]-$start;
       $unit[4]=$unit[4]-$start;
       my $line=join("\t",@unit);
       print OUT1 "$line";
     }

   }
   close OUT1;
   close IN1;
}else{
   open IN1, "$gff" or die "$!";
   open OUT1, ">>$head.te.gff" or die "$!";
   while (<IN1>){
     my @unit=split("\t",$_);
     $unit[0]=~s/\_//;
     if ($unit[0] eq $chr and $unit[3] >= $start and $unit[4] <= $end){
        my $tempend   =$len-($unit[3]-$start);
        my $tempstart =$len-($unit[4]-$start);
        $unit[3]=$tempstart;
        $unit[4]=$tempend;
        if ($unit[6] eq "+"){
            $unit[6] = "-";
        }else{
            $unit[6] = "+";
        }
        my $line=join("\t",@unit);
        print OUT1 "$line";
     }
   }
   close OUT1;
   close IN1;
}
}

sub getsubgff3
{
my ($gff,$head,$strand)=@_;
my @temp=split("_",$head);
my $chr=$temp[1];
my $start=$temp[2];
my $end  =$temp[3];
my $len=$end-$start+1;
print "$gff\t$chr\t$start\t$end\n";
if ($strand eq "+"){
   open IN1, "$gff" or die "$!";
   open OUT1, ">>$head.gene.gff" or die "$!";
   while (<IN1>){
     my @unit=split("\t",$_);
     $unit[0]=~s/\_//;
     if ($unit[0] eq $chr and $unit[3] >= $start and $unit[4] <= $end){ 
       $unit[3]=$unit[3]-$start;
       $unit[4]=$unit[4]-$start;
       my $line=join("\t",@unit);
       print OUT1 "$line";
     } 

   }
   close OUT1;
   close IN1;
}else{
   open IN1, "$gff" or die "$!";
   open OUT1, ">>$head.gene.gff" or die "$!";
   while (<IN1>){
     my @unit=split("\t",$_);
     $unit[0]=~s/\_//;
     if ($unit[0] eq $chr and $unit[3] >= $start and $unit[4] <= $end){ 
        my $tempend   =$len-($unit[3]-$start);
        my $tempstart =$len-($unit[4]-$start); 
        $unit[3]=$tempstart;
        $unit[4]=$tempend;
        if ($unit[6] eq "+"){
            $unit[6] = "-";
        }else{
            $unit[6] = "+";
        }
        my $line=join("\t",@unit);
        print OUT1 "$line";
     }
   }
   close OUT1;
   close IN1;
}
}


################
sub getsubfasta
{
my ($fasta,$head,$strand)=@_;
my @temp=split("_",$head);
my $chr=$temp[1];
my $start=$temp[2];
my $end  =$temp[3];
my $len=$end-$start+1;
my $refseq=getfastaseq($fasta);
if ($strand eq "+"){
   if (exists $refseq->{$chr}){
       my $subseq=substr($refseq->{$chr},$start,$len);
       open FAS, ">$head.fasta" or die "$!"; 
             print FAS ">$head\n$subseq\n";
       close FAS;  
   }else{
       print "$chr can not found in $fasta\n";
   }
}else{
   if (exists $refseq->{$chr}){
       my $subseq=substr($refseq->{$chr},$start,$len);
       $subseqrec=revcom($subseq);
       open FAS, ">$head.fasta" or die "$!";
             print FAS ">$head rec\n$subseqrec\n";
       close FAS;
   }else{
       print "$chr can not found in $fasta\n";
   }   
}
}

################
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
    $head=~s/\_//;
    my $seq=join("\n",@unit);
    $seq=~s/\>//g;
    $seq=~s/\n//g;
    $seq=~s/\s//g;
    #print "$head\n";
    $hash{$head}=$seq;
}
close IN;
$/="\n";
return \%hash;
}
 
