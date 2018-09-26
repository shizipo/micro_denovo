#!/usr/bin/perl -w
use strict;
use warnings;
use FindBin qw($Bin);
use Getopt::Long;
use Cwd;
my $qsub="/home/fanyucai/software/qsub/qsub-pbs.pl";
my $diamond="/home/fanyucai/software/diamond/diamond_v0.9.7/diamond";
my $index="/public/land/database/nr_diamond/nr_2017.dmnd";
my $idba_ud="/home/fanyucai/software/idba/idba-master/bin";
my $blast2lca="/home/fanyucai/software/MEGAN6/tools/blast2lca";
my $g2t="/public/land/database/MEGAN6/prot_acc2tax-May2017.abin";
my $R_env="export PATH=/home/fanyucai/software/gcc/gcc-v6.1.0/bin/:/home/fanyucai/software/R/R-v3.4.0/bin/:\$PATH && export LD_LIBRARY_PATH=/home/fanyucai/software/gcc/gcc-v6.1.0/lib64/:\$LD_LIBRARY_PATH";
my $xvf_run="/home/fanyucai/bin/xvfb-run";
my($outdir,@prefix,@pe1,@pe2,@fa,$num,$type,$db,$queue);
$outdir||=getcwd;
GetOptions(
     "pe1:s{1,}"=>\@pe1,
     "pe2:s{1,}"=>\@pe2,
     "fa:s{1,}"=>\@fa,
     "p:s{1,}"=>\@prefix,
     "o:s"=>\$outdir,
     "queue:s"=>\$queue,
           );
$num ||=15;
$queue||="all";

sub usage{
    print qq{
This script will map fasta(fastq) data to nr database by diamond and parse species abudance.
usage:
perl $0 -pe1 sample1_1.fq sample2_1.fq -pe2 sample1_2.fq sample2_2.fq  -o $outdir -p sample1 sample2
            or
perl $0 -fa sample1.fa sample2.fa -o $outdir -p sample1 sample2

options:
-pe1                    5' read
-pe2                    3' read
-p                      the prefix of output(force)
-o                      output directory(default:$outdir)
-fa                     input fasta file
-queue			which queue you will run(default:all)

Buchfink B, Xie C, Huson D H. Fast and sensitive protein alignment using DIAMOND[J]. Nature methods, 2015, 12(1): 59-60.

Email:fanyucai1\@126.com
2017.9.21
};
    exit;
}
if(!@prefix||!$outdir)
{
    &usage();
}
system "mkdir -p $outdir/";
###############################run in fastq mode
#if(@pe1 && @pe2)
#{
#    open(D1,">$outdir/diamond1.sh");
#    for(my $k=0;$k<=$#prefix;$k++)
#    {
#        print D1 "$R_env && $idba_ud/fq2fa --merge $pe1[$k] $pe2[$k] $outdir/$prefix[$k].fa && $diamond blastx -q $outdir/$prefix[$k].fa -d $index -p 20 -k 250 --min-score 40 --outfmt 100 -e 0.00001 -o $outdir/$prefix[$k].daa\n";
#    }
#    system "perl $qsub --queue $queue --ppn 5 $outdir/diamond1.sh";
#}
################################run in fasta mode
#if(@fa)
#{
#    open(D1,">$outdir/diamond1.sh");
#    for(my $k=0;$k<=$#prefix;$k++)
#    {
#        print D1 "$diamond blastx -q $fa[$k] -d $index  -p 20 -k 250 --outfmt 100 --min-score 40 -e 0.00001 -o $outdir/$prefix[$k].daa\n";
#    }
#   system "perl $qsub --queue $queue --ppn 5 $outdir/diamond1.sh";
#}
################################parse the diamond
#open(D2,">$outdir/diamond2.sh");
#for(my $k=0;$k<=$#prefix;$k++)
#{
#    print D2 "$blast2lca -a2t $g2t -i $outdir/$prefix[$k].daa -top 1 -sr false -o $outdir/$prefix[$k].tax.out\n";
#}
#system "perl $qsub --queue $queue $outdir/diamond2.sh";
###############################statistics and plot
my $tmp;
for(my $k=0;$k<=$#prefix;$k++)
{
    my %hash;
    open(TAX,"$outdir/$prefix[$k].tax.out");
    while(<TAX>)
    {
        chomp;
        $tmp++;
        my @array=split(/;/,$_);
        if(defined $array[2] && $array[2]=~/[0-9A-Z]/i)
        {
            $array[$#array-1]=~s/^\s+//;
            $array[$#array-1]=~s/\s+/_/g;
            $hash{$array[$#array-1]}++;
        }
    }
    open(TAX2,">$outdir/$prefix[$k].tax.txt");
    my $num=0;
    my $other=0;
    print TAX2 "Species\tCounts\n";
    foreach my $key(sort {$hash{$b}<=>$hash{$a}} keys %hash)
    {
        $num++;
        if($num<=10)
        {
            $tmp-=$hash{$key};
            print TAX2 $key,"\t",$hash{$key},"\n";
        }
        else
        {
            $other+=$hash{$key};
        }
    }
    $tmp-=$other;
    #print TAX2 "other","\t",$other,"\n";
    #print TAX2 "unclassfied","\t",$tmp,"\n";
    system "echo '#!/usr/bin/env Rscript
    library(ggplot2)
    a=read.table(\"$outdir/$prefix[$k].tax.txt\",header=T,sep=\"\\t\",quote=\"\",check.names=F)
    a\$Species=factor(a\$Species,levels=a\$Species);
    lbls=a[,1]
    slices=a[,2]
    pct <- round(slices/sum(slices)*100,digits=2)
    lbls <- paste(lbls,\" (\",slices,\", \", pct,\"%)\",sep=\"\")
    p=ggplot(a, aes(x = \"\", y = Counts, fill = Species)) +
    geom_bar(stat = \"identity\", width = 1) +
    coord_polar(theta = \"y\") +
    theme_bw() +
    theme(axis.ticks = element_blank()) +
    theme(axis.text.x = element_blank()) +
    labs(x = \"\", y = \"\", title = \"Top 10 species distribution of $prefix[$k]\") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.title = element_blank(), legend.position = \"right\") + 
    scale_fill_discrete(breaks = a\$Species, labels = lbls) +
    theme(panel.grid=element_blank()) + 
    theme(panel.border=element_blank())
    ggsave(\"$outdir/$prefix[$k].tax.pdf\",width=8,height=6,plot=p)
    ggsave(\"$outdir/$prefix[$k].tax.png\",type=\"cairo-png\",width=8,height=6,plot=p)
    '>$outdir/$prefix[$k].pie.Rscript";
    system "$R_env && Rscript $outdir/$prefix[$k].pie.Rscript";
}
