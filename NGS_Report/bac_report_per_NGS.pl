#!/usr/bin/perl
use strict;
use warnings;
use Cwd;
use FindBin qw($Bin $Script);
use File::Basename;
use File::Spec;
use IO::File;
use Getopt::Long;
use File::Find;
use XML::Writer;
use Config::IniFiles;


my $xml2html="/home/fanyucai/script/xml2html/xml2HtmlConverter.py";
my $python="/home/fanyucai/software/python/Python-v2.7.13/bin/python";
my $topng="sh /home/yangxinchao/micro/bac_pipeline/bin/topng_v1.0.sh";

my ($idir, $odir, $sample, $config, $step, $help);
GetOptions(
	"idir:s"    =>\$idir,
	"odir:s"    =>\$odir,
	"sample:s"  =>\$sample,
	"config:s"  =>\$config,
	"step:s"    =>\$step,
	"help!"     =>\$help,
);

my $usage =<<"USAGE";
Program: $0
Description: generate a html format annotation report
Options:
	-idir      <indir>     input analysis dir or step1 outdir     [Required]
	-sample    <str>       sample name,eg:AA,BB,CC                [Required]
	-config    <file>      config file                            [Required for step2]
	                       if run:perl bac_report.pl -config help, get file format.
	-step      <step>      1:copy, 2:report. default:1,2          [Optional]
	-odir      <outdir>    output dir, default: ./                [Optional]
	-help                  print help information
Example
	1. only run step1:   perl bac_report.pl -idir OExxxxx -sample AA,BB -step 1 -odir result
	2. only run step2:   perl bac_report.pl -idir result -sample AA,BB -config config.ini -step 2 -odir result
	3. run step1,2   :   perl bac_report.pl -idir OExxxxx -sample AA,BB -config config.ini -step 1,2 -odir result

USAGE

$step ||= "1,2";
$odir ||= getcwd;

if(defined $config && $config=~/help/){
	print "[report]\ntitle=二代细菌框架图结题报告\n\n[par]\n#客户姓名\ncustomer=XXX\n#项目号\ncompact_NUM=OE2017HXXXX\n#样本数目\nsample_num=6\n#执行编号\nJob_number=OEXXXX\n#组装策略\nassembly=\"以三代测序为主，使用二代测序数据作为辅助，使用SPAdes组装软件，对测序数据进行denovo组装。\"\n";
	exit;
}elsif(defined $config && $config!~/help/){
	(-s $config) || die "Error: don't open file $config!\n";
	$config=File::Spec->rel2abs($config);
}

die $usage if($help || !$idir);
die "Error: -sample is required !\n" if(!$sample);
die "Error: -config is required for step2!\n" if($step=~/2/ && !$config);

(-d $idir) || die "Error: not find DIR: $idir!\n";
$idir = File::Spec->rel2abs($idir);$idir=~s/\/$//g;
(-d $odir) || mkdir $odir;
$odir = File::Spec->rel2abs($odir);$odir=~s/\/$//g;

if($step=~/1/){
	print "########## step1 start ######################\n";
	my @sample=split /,/,$sample;
	open CP,">$odir/copy_tmp.sh" || die $!;
	my $dirname=0;
#	if(-d "$idir/Assembly/"){
#		print CP "##Assembly\n";
#		$dirname++;
#		for my $s(@sample){
#			my $assembly_outdir;
#			if(@sample==1){
#				$assembly_outdir="$odir/$dirname.Assembly";
#			}else{
#				$assembly_outdir="$odir/$s/$dirname.Assembly";
#			}
#			print CP "[ -d $assembly_outdir ] || mkdir -p $assembly_outdir\n";
#			print CP "cp $idir/Assembly/{$s.fasta,$s.assembly_summary.xls} $assembly_outdir\n";
#		}
#		print CP "\n";
#	}
    $dirname++;
    $dirname++;
	if(-d "$idir/Genome_structure"){
		print CP "###Genome_structure\n" if(-d "$idir/Genome_structure/");
		my %combind_gff;
		$dirname++;
		if(-d "$idir/Genome_structure/Gene_predict"){
			print CP "###Gene_predict\n";
			for my $s(@sample){
				my $gene_outdir;
				if(@sample==1){
					$gene_outdir="$odir/$dirname.Genome_structure/Gene_predict";
				}else{
					$gene_outdir="$odir/$s/$dirname.Genome_structure/Gene_predict";
				}
				print CP "[ -d $gene_outdir ] || mkdir -p $gene_outdir\n";
				print CP "cp $idir/Genome_structure/Gene_predict/{$s.ffn,$s.faa,$s.gene.gff,$s.gene.summary.xls} $gene_outdir\n";
				$combind_gff{$s}.="$gene_outdir/$s.gene.gff ";
			}
			print CP "\n";
		}
		if(-d "$idir/Genome_structure/ncRNA"){
			print CP "###ncRNA\n";
			for my $s(@sample){
				my $ncRNA_outdir;
				if(@sample==1){
					$ncRNA_outdir="$odir/$dirname.Genome_structure/ncRNA";
				}else{
					$ncRNA_outdir="$odir/$s/$dirname.Genome_structure/ncRNA";
				}
				print CP "[ -d $ncRNA_outdir ] || mkdir -p $ncRNA_outdir\n";
				print CP "cp $idir/Genome_structure/ncRNA/{$s.ncRNA.gff,$s.rRNA.fna,$s.tRNA.structure,$s.sRNA.fna} $ncRNA_outdir\n";
				$combind_gff{$s}.="$ncRNA_outdir/$s.ncRNA.gff ";
			}
			print CP "\n";
		}
		if(-d "$idir/Genome_structure/Repeat_sequence"){
			print CP "###Repeat_sequence\n";
			for my $s(@sample){
				my $repeat_outdir=(@sample==1)?("$odir/$dirname.Genome_structure/Repeat_sequence"):("$odir/$s/$dirname.Genome_structure/Repeat_sequence");
				print CP "[ -d $repeat_outdir ] || mkdir -p $repeat_outdir\n";
				print CP "cp $idir/Genome_structure/Repeat_sequence/{$s.masked.fna,$s.repeat.gff,$s.repeat.xls} $repeat_outdir\n";
				my $delete="yes";
				open IN,"<$idir/Genome_structure/Repeat_sequence/$s.repeat.xls" || die $!;
				while(<IN>){
					chomp;
					next if(/^#/ || /^\s/ || /^$/ || /^--/);
					my @line=split /\t/;
					$delete="no" if($line[2]!=0);
				}
				close IN;
				if($delete eq "yes"){
					print CP "###$s: Repeat_sequence result is empty, directory[Repeat_sequence] will be removed\n";
					print CP "rm -rf $repeat_outdir\n";
					print "###$s: Repeat_sequence result is empty, directory[Repeat_sequence] will be removed\n";
				}else{
					$combind_gff{$s}.="$repeat_outdir/$s.repeat.gff ";
				}
			}
			print CP "\n";
		}
		if((keys %combind_gff)>1){
			for my $tmp(keys %combind_gff){
				print CP "cat $combind_gff{$tmp} > $odir/$tmp/$dirname.Genome_structure/$tmp.gff\nrm -rf $combind_gff{$tmp}\n";
			}
		}elsif((keys %combind_gff)==1){
			for my $tmp(keys %combind_gff){
				print CP "cat $combind_gff{$tmp} > $odir/$dirname.Genome_structure/$tmp.gff\nrm -rf $combind_gff{$tmp}\n";
			}
		}
		print CP "\n";
		if(-d "$idir/Genome_structure/CRISPR"){
			print CP "###CRISPR\n";
			for my $s(@sample){
				my $crispr_outdir=(@sample==1)?("$odir/$dirname.Genome_structure/CRISPR"):("$odir/$s/$dirname.Genome_structure/CRISPR");
				print CP "[ -d $crispr_outdir ] || mkdir -p $crispr_outdir\n";
				print CP "cp $idir/Genome_structure/CRISPR/{$s.crt.xls,$s.pilercr.xls,$s.CRISPR.summary.xls} $crispr_outdir\n";
				my $delete="yes";
				open IN,"<$idir/Genome_structure/CRISPR/$s.CRISPR.summary.xls" || die $!;
				while(<IN>){
					chomp;
					next if(/^Sample\tPILER-CR\(\#\)/);
					my @line=split /\t/;
					$delete="no" if($line[1]!=0);
					$delete="no" if($line[2]!=0);
				}
				close IN;
				if($delete eq "yes"){
					print CP "###$s: CRISPR result is empty, directory[CRISPR] will be removed\n";
					print CP "rm -rf $crispr_outdir\n";
					print "###$s: CRISPR result is empty, directory[CRISPR] will be removed\n";
				}
			}
			print CP "\n";
		}
		if(-d "$idir/Genome_structure/Prophage"){
			print CP "###Prophage\n";
			for my $s(@sample){
				my $prophage_outdir=(@sample==1)?("$odir/$dirname.Genome_structure/Prophage"):("$odir/$s/$dirname.Genome_structure/Prophage");
				print CP "[ -d $prophage_outdir ] || mkdir -p $prophage_outdir\n";
				print CP "cp $idir/Genome_structure/Prophage/$s.prophage.xls $prophage_outdir\n";
				my $delete="yes";my $phage_num=0;
				open IN,"<$idir/Genome_structure/Prophage/$s.prophage.xls" || die $!;
				while(<IN>){
					chomp;
					next if(/^Prophage ID\tChr\tStart/);
					next if(/Not found prophage/i);
					$phage_num++;
				}
				close IN;
				$delete="no" if($phage_num>0);
				if($delete eq "yes"){
					print CP "###$s: Prophage result is empty, directory[Prophage] will be removed\n";
					print CP "rm -rf $prophage_outdir\n";
					print "###$s: Prophage result is empty, directory[Prophage] will be removed\n";
				}
			}
			print CP "\n";
		}
	}
	print CP "\n############################################################\n";
	
	if(-d "$idir/Annotation/Basic_function"){
		print CP "###Basic_function\n";
		$dirname++;
		for my $s(@sample){
			my $Basic_function=(@sample==1)?("$odir/$dirname.Basic_function"):("$odir/$s/$dirname.Basic_function");
			print CP "[ -d $Basic_function ] || mkdir -p $Basic_function\n";
			print CP "cp $idir/Annotation/Basic_function/$s/result/{$s.generalDB.anno.xls,$s.generalDB.stat.xls} $Basic_function/\n";
			if(-d "$idir/Annotation/Basic_function/$s/result/NR"){
				print CP "[ -d $Basic_function/NR ] || mkdir -p $Basic_function/NR\n";
				print CP "cp $idir/Annotation/Basic_function/$s/result/NR/* $Basic_function/NR\n";
				print CP "$topng $Basic_function/NR/*pdf\n";
			}
			if(-d "$idir/Annotation/Basic_function/$s/result/Swissprot"){
				print CP "[ -d $Basic_function/Swissprot ] || mkdir -p $Basic_function/Swissprot\n";
				print CP "cp $idir/Annotation/Basic_function/$s/result/Swissprot/* $Basic_function/Swissprot\n";
			}
			if(-d "$idir/Annotation/Basic_function/$s/result/eggNOG"){
				print CP "[ -d $Basic_function/eggNOG ] || mkdir -p $Basic_function/eggNOG\n";
				print CP "cp $idir/Annotation/Basic_function/$s/result/eggNOG/* $Basic_function/eggNOG\n";
				print CP "$topng $Basic_function/eggNOG/*pdf\n";
			}
			if(-d "$idir/Annotation/Basic_function/$s/result/GO"){
				print CP "[ -d $Basic_function/GO ] || mkdir -p $Basic_function/GO\n";
				print CP "cp $idir/Annotation/Basic_function/$s/result/GO/* $Basic_function/GO\n";
				print CP "$topng $Basic_function/GO/*pdf\n";
			}
			if(-d "$idir/Annotation/Basic_function/$s/result/KEGG"){
				print CP "[ -d $Basic_function/KEGG ] || mkdir -p $Basic_function/KEGG\n";
				print CP "cp -a $idir/Annotation/Basic_function/$s/result/KEGG/* $Basic_function/KEGG\n";
				print CP "$topng $Basic_function/KEGG/*pdf\n";
			}
			if(-d "$idir/Annotation/Basic_function/$s/result/KOG"){
				print CP "[ -d $Basic_function/KOG ] || mkdir -p $Basic_function/KOG\n";
				print CP "cp $idir/Annotation/Basic_function/$s/result/KOG/* $Basic_function/KOG\n";
				print CP "$topng $Basic_function/KOG/*pdf\n";
			}
			if(-d "$idir/Annotation/Basic_function/$s/result/COG"){
				print CP "[ -d $Basic_function/COG ] || mkdir -p $Basic_function/COG\n";
				print CP "cp $idir/Annotation/Basic_function/$s/result/COG/* $Basic_function/COG\n";
				print CP "$topng $Basic_function/COG/*pdf\n";
			}
			if(-d "$idir/Annotation/Basic_function/$s/result/Pfam"){
				print CP "[ -d $Basic_function/Pfam ] || mkdir -p $Basic_function/Pfam\n";
				print CP "cp $idir/Annotation/Basic_function/$s/result/Pfam/* $Basic_function/Pfam\n";
				my $delete="yes";my $pfam_num=0;
				open IN,"<$idir/Annotation/Basic_function/$s/result/Pfam/$s.Pfam.align.anno.xls" || die $!;
				while(<IN>){
					chomp;
					next if(/^#Gene_ID/);
					$pfam_num++;
				}
				close IN;
				$delete="no" if($pfam_num>0);
				if($delete eq "yes"){
					print CP "###$s: Pfam result is empty, directory[Pfam] will be removed\n";
					print CP "rm -rf $Basic_function/Pfam\n";
					print "###$s: Pfam result is empty, directory[Pfam] will be removed\n";
				}
			}
		}
		print CP "\n";
	}
	if(-d "$idir/Annotation/Specialized_function"){
		print CP "###Specialized_function\n";
		$dirname++;
		if(-d "$idir/Annotation/Specialized_function/CARD"){
			print CP "###CARD\n";
			for my $s(@sample){
				my $card_outdir=(@sample==1)?("$odir/$dirname.Specialized_function/CARD"):("$odir/$s/$dirname.Specialized_function/CARD");
				print CP "[ -d $card_outdir ] || mkdir -p $card_outdir\n";
				print CP "cp $idir/Annotation/Specialized_function/CARD/{$s.card.anno.xls,$s.card.pdf,$s.statics.txt} $card_outdir\n";
				print CP "$topng $card_outdir/*pdf\n";
				my $delete="yes";my $line_num=0;
				open IN,"<$idir/Annotation/Specialized_function/CARD/$s.card.anno.xls" || die $!;
				while(<IN>){
					chomp;
					next if(/^#Query_seq/);
					$line_num++;
				}
				close IN;
				$delete="no" if($line_num>0);
				if($delete eq "yes"){
					print CP "###$s: CARD result is empty, directory[CARD] will be removed\n";
					print CP "rm -rf $card_outdir\n";
					print "###$s: CARD result is empty, directory[CARD] will be removed\n";
				}
			}
			print CP "\n";
		}
		if(-d "$idir/Annotation/Specialized_function/CAZy"){
			print CP "###CAZy\n";
			for my $s(@sample){
				my $cazy_outdir=(@sample==1)?("$odir/$dirname.Specialized_function/CAZy"):("$odir/$s/$dirname.Specialized_function/CAZy");
				print CP "[ -d $cazy_outdir ] || mkdir -p $cazy_outdir\n";
				print CP "cp $idir/Annotation/Specialized_function/CAZy/{$s.CAZy.anno.xls,$s.CAZy.class.xls,$s.CAZy.family.xls,$s.CAZy.class.pdf,$s.CAZy.class.png} $cazy_outdir\n";
				print CP "$topng $cazy_outdir/*pdf\n";
				my $delete="yes";my $line_num=0;
				open IN,"<$idir/Annotation/Specialized_function/CAZy/$s.CAZy.anno.xls" || die $!;
				while(<IN>){
					chomp;
					next if(/^Gene\tFamily\tEvalue/);
					$line_num++;
				}
				close IN;
				$delete="no" if($line_num>0);
				if($delete eq "yes"){
					print CP "###$s: CAZy result is empty, directory[CAZy] will be removed\n";
					print CP "rm -rf $cazy_outdir\n";
					print "###$s: CAZy result is empty, directory[CAZy] will be removed\n";
				}
			}
			print CP "\n";
		}
		if(-d "$idir/Annotation/Specialized_function/PHI"){
			print CP "###PHI\n";
			for my $s(@sample){
				my $phi_outdir=(@sample==1)?("$odir/$dirname.Specialized_function/PHI"):("$odir/$s/$dirname.Specialized_function/PHI");
				print CP "[ -d $phi_outdir ] || mkdir -p $phi_outdir\n";
				print CP "cp $idir/Annotation/Specialized_function/PHI/$s.PHI.anno.xls $phi_outdir\n";
				my $delete="yes";my $line_num=0;
				open IN,"<$idir/Annotation/Specialized_function/PHI/$s.PHI.anno.xls" || die $!;
				while(<IN>){
					chomp;
					next if(/^#Gene_ID/ || /^\s/);
					$line_num++;
				}
				close IN;
				$delete="no" if($line_num>0);
				if($delete eq "yes"){
					print CP "###$s: PHI result is empty, directory[PHI] will be removed\n";
					print CP "rm -rf $phi_outdir\n";
					print "###$s: PHI result is empty, directory[PHI] will be removed\n";
				}
			}
			print CP "\n";
		}
		if(-d "$idir/Annotation/Specialized_function/VFDB"){
			print CP "###VFDB\n";
			for my $s(@sample){
				my $vfdb_outdir=(@sample==1)?("$odir/$dirname.Specialized_function/VFDB"):("$odir/$s/$dirname.Specialized_function/VFDB");
				print CP "[ -d $vfdb_outdir ] || mkdir -p $vfdb_outdir\n";
				print CP "cp $idir/Annotation/Specialized_function/VFDB/{$s.VFDB.anno.xls,VFs_reference.xls} $vfdb_outdir\n";
				my $delete="yes";my $line_num=0;
				open IN,"<$idir/Annotation/Specialized_function/VFDB/$s.VFDB.anno.xls" || die $!;
				while(<IN>){
					chomp;
					next if(/^#Gene_ID/);
					$line_num++;
				}
				close IN;
				$delete="no" if($line_num>0);
				if($delete eq "yes"){
					print CP "###$s: VFDB result is empty, directory[VFDB] will be removed\n";
					print CP "rm -rf $vfdb_outdir\n";
					print "###$s: VFDB result is empty, directory[VFDB] will be removed\n";
				}
			}
			print CP "\n";
		}
		if(-d "$idir/Annotation/Specialized_function/SignalP"){
			print CP "###SignalP\n";
			$dirname++;
			for my $s(@sample){
				my $signalp_outdir=(@sample==1)?("$odir/$dirname.SignalP"):("$odir/$s/$dirname.SignalP");
				print CP "[ -d $signalp_outdir ] || mkdir -p $signalp_outdir\n";
				print CP "cp $idir/Annotation/Specialized_function/SignalP/$s.signalp.xls $signalp_outdir\n";
			}
			print CP "\n";
		}
		if(-d "$idir/Annotation/Specialized_function/EffectiveT3"){
			print CP "###EffectiveT3\n";
			$dirname++;
			for my $s(@sample){
				my $t3ss_outdir=(@sample==1)?("$odir/$dirname.EffectiveT3"):("$odir/$s/$dirname.EffectiveT3");
				print CP "[ -d $t3ss_outdir ] || mkdir -p $t3ss_outdir\n";
				print CP "cp $idir/Annotation/Specialized_function/EffectiveT3/$s.T3SS.xls $t3ss_outdir\n";
			}
			print CP "\n";
		}
		if(-d "$idir/Annotation/Specialized_function/TMHMM"){
			print CP "###TMHMM\n";
			$dirname++;
			for my $s(@sample){
				my $tmhmm_outdir=(@sample==1)?("$odir/$dirname.TMHMM"):("$odir/$s/$dirname.TMHMM");
				print CP "[ -d $tmhmm_outdir ] || mkdir -p $tmhmm_outdir\n";
				print CP "cp $idir/Annotation/Specialized_function/TMHMM/$s.tmhmm.xls $tmhmm_outdir\n";
				print CP "[ -d $tmhmm_outdir/$s\_TMHMM_plot ] || mkdir -p $tmhmm_outdir/$s\_TMHMM_plot\n";
				print CP "cp -a $idir/Annotation/Specialized_function/TMHMM/$s\_TMHMM_plot/*png $tmhmm_outdir/$s\_TMHMM_plot\n";
			}
			print CP "\n";
		}
	}
	if(-d "$idir/Circle_map"){
		print CP "############################################################\n";
		print CP "###Circle_map\n";
		$dirname++;
		for my $s(@sample){
			my $Circle_map=(@sample==1)?("$odir/$dirname.Circle_map"):("$odir/$s/$dirname.Circle_map");
			print CP "[ -d $Circle_map ] || mkdir -p $Circle_map\n";
			print CP "cp $idir/Circle_map/$s/{$s.circle.png,$s.circle.svg,COG_legend.png,COG_legend.svg} $Circle_map\n";
		}
		print CP "\n";
	}
	close CP;
	system("sh $odir/copy_tmp.sh");
	$idir=$odir;
}

if($step=~/2/){
	print "\n########## step2 start ######################";
	system("rm -rf $odir/copy_tmp.sh");
	my @sample=split /,/,$sample;
	my %ini;
	tie %ini, 'Config::IniFiles', (-file =>$config);
	my $initial_idir=$idir;
	my $initial_odir=$odir;
	for my $s(@sample){
		print "\n##$s web report start ...\n";
		$odir=(@sample==1)?("$initial_odir"):("$initial_odir/$s");
		$idir=(@sample==1)?("$initial_idir"):("$initial_idir/$s");
        my $sample=$ini{par}{sample};
        my $species=$ini{par}{species};
		if(-d "$odir/pic"){system("rm -rf $odir/pic");};
		system("mkdir -p $odir/pic");
		system("cp $Bin/pic/* $odir/pic");
		open PRO,">$odir/pic/project.txt" || die $!;
		print PRO "客户姓名\t$ini{par}{customer}\n";
		print PRO "合同编号\t$ini{par}{compact_NUM}\n";
		print PRO "样本数量\t$ini{par}{sample_num}\n";
        print PRO "物种名\t$ini{par}{species}\n";
		print PRO "样品名称\t$s\n";
		print PRO "执行编号\t$ini{par}{Job_number}\n";
		close PRO;
		my $output = IO::File->new(">$odir/$s.xml");
		my $writer = XML::Writer->new(OUTPUT => $output,NAMESPACES => 1);
		$writer->xmlDecl("UTF-8");
		
		#first report title
		print $output "\n";$writer->startTag("report");
		print $output "\n";$writer->emptyTag("report_name", "value"=> "$ini{report}{title}");
		#abstract
		print $output "\n";$writer->emptyTag("report_abstract", "value"=>"");
		#客户信息表格
		print $output "\n";$writer->emptyTag('h1', 'name'=>"项目信息",'type'=>"一级标题显示样式",'desc'=>"一级标题描述");
		print $output "\n";$writer->emptyTag('table','name'=>"",'type'=>"type1|full",'desc'=>"",'path'=>"pic/project.txt");
		
		print $output "\n";$writer->emptyTag('h1', 'name'=>"生物信息分析流程图",'type'=>"一级标题显示样式",'desc'=>"一级标题描述");
		my $title1Num=0;my $title2Num=0;my $title3Num=0;my $table=0;my $pic=0;
		print $output "\n";$writer->emptyTag('pic','name'=>"",'desc'=>"",'path'=>"pic/microbial_denovo.png");
        if(-d "$idir/1.QC/"){

            #建库测序
            print $output "\n";$writer->emptyTag('h1', 'name'=>++$title1Num."建库测序",'type'=>"一级标题显示样式",'desc'=>"一级标题描述");
            print $output "\n";$writer->emptyTag('p','type'=>'正文段落显示样式','desc'=>"提取样本中的基因组DNA，电泳检测合格后进行建库。");
            print $output "\n";$writer->emptyTag('p','type'=>'正文段落显示样式','desc'=>"检测合格的DNA样品先经过Covaris随机打断成350bp的片段，采用TruSeq DNA LT Sample Prep kit试剂盒进行建库，DNA片段经过末端修复、加A尾、加测序接头、纯化、PCR扩增等步骤，最终完成文库构建。文库检验合格后利用测序仪进行双端测序。");
            print $output "\n";$writer->emptyTag('pic',name=>"图 ".++$pic."  建库流程",type=>'img-width-max',desc=>'',path=>"pic/method.png");


            ##测序数据质控及统计
            print $output "\n";$writer->emptyTag('h1', 'name'=>++$title1Num." 原始数据处理以及质控",'type'=>"一级标题显示样式",'desc'=>"一级标题描述");
            #原始测序数据
            print $output "\n";$writer->emptyTag('h2', 'name'=>$title1Num.".".++$title2Num." 原始测序数据",'type'=>"二级标题显示样式",'desc'=>"二级标题描述");
            print $output "\n";$writer->emptyTag('p','type'=>'正文段落显示样式','desc'=>"高通量测序得到的原始图像数据文件经碱基识别 (Base Calling) 分析转化为原始测序序列，我们称之为 Raw Data 或 Raw Reads，结果以 FASTQ（简称为 fq）文件格式存储，其中包含测序序列 (reads) 的序列信息以及其对应的测序质量信息。FASTQ 格式文件中每个 read 由四行描述，如下图所示：");
            print $output "\n";$writer->emptyTag('pic',name=>"图 ".++$pic." 原始数据示例",type=>'img-width-max',desc=>'',path=>"pic/FASTQ.png");
            
            print $output "\n";$writer->emptyTag('p','type'=>'正文段落显示样式','desc'=>"每个序列共有4行，第1行和第3行是序列名称（有的 fq 文件为了节省存储空间会省略第三行“＋”后面的序列名称），由测序仪产生；第2行是序列；第4行是序列的测序质量，每个字符对应第2行每个碱基。第四行中每个字符对应的 ASCII 值减去33，即为对应第二行碱基的测序质量值。如果测序错误率用 E 表示，Illumina 的碱基质量值用 Qphred 表示，则有下列关系：");
            print $output "\n";$writer->emptyTag('pic',name=>"",type=>'img-width-max',desc=>'',path=>"pic/Qphred.png");
            print $output "\n";$writer->emptyTag('table','name'=>"表" .++$table. " Illumina 测序错误率与测序质量值简明对应关系表",'type'=>"type1|full",'desc'=>"",'path'=>"pic/filter_reads.xls");
            
            #有效数据
            #数据过滤的方法
            print $output "\n";$writer->emptyTag('h2', 'name'=>$title1Num.".".++$title2Num." 数据过滤的方法",'type'=>"二级标题显示样式",'desc'=>"二级标题描述");
            print $output "\n";$writer->emptyTag('p','type'=>'正文段落显示样式','desc'=>"原始测序数据中会包含接头信息，低质量碱基，未测出的碱基（以N表示），这些信息会对后续的信息分析造成很大的干扰，通过精细的过滤方法将这些干扰信息去除掉，最终得到的数据即为有效数据，我们称之为clean data 或 clean reads，该文件的数据格式与Raw data完全一样。");
            print $output "\n";$writer->emptyTag('p','type'=>'正文段落显示样式','desc'=>"通过 Illumina 平台，得到了大量的样本双端测序数据。鉴于数据错误率对结果的影响，我们采用 Trimmomatic<sup>[<a href='#ref1'>1</a>]</sup> 软件对原始数据进行质量预处理，并对整个质控过程中的 reads 数进行统计汇总。");
            print $output "\n";$writer->emptyTag('p','type'=>'正文段落显示样式','desc'=>"质量预处理步骤：");
            print $output "\n";$writer->emptyTag('p','type'=>'正文段落显示样式','desc'=>"(1) 去接头 (Adaptor)；");
            print $output "\n";$writer->emptyTag('p','type'=>'正文段落显示样式','desc'=>"(2) 去除低质量 Reads；");
            print $output "\n";$writer->emptyTag('p','type'=>'正文段落显示样式','desc'=>"(3) 从3’端及5‘端以不同方式去除低质量碱基；");
            print $output "\n";$writer->emptyTag('p','type'=>'正文段落显示样式','desc'=>"(4) 统计原始测序量，有效测序量，Q30，GC 含量，并进行综合评价。");
            
            #质控结果
            print $output "\n";$writer->emptyTag('h2', 'name'=>$title1Num.".".++$title2Num." 质控结果",'type'=>"二级标题显示样式",'desc'=>"二级标题描述");
            print $output "\n";$writer->emptyTag('p','type'=>'正文段落显示样式','desc'=>"经过对测序数据的严格过滤，得到高质量的clean data。对产出数据进行统计（见表2），包括测序read数量，数据产量，Q30含量，GC含量等。");
            print $output "\n";$writer->emptyTag('table','name'=>"表" .++$table. " 测序数据产出统计 ",'type'=>"type1|full",'desc'=>"
                备注：<br/>
                    (1) Sample：样本名;<br/>
                    (2) Clean_Reads ：质控后的有效reads数，以M（106）为单位; <br/>
                    (3) Clean_Bases ：质控后的有效数据量，以Gp为单位;<br/>
                    (4) Q30 (%）：大于Q30的碱基占Raw_Bases的比例;<br/>
                    (5) GC (%) ：原始测序数据G和C碱基数占总碱基数的比例;<br/>
                    (6) Clean_Reads ：质控后的有效reads数，以M（106）为单位; <br/>
                    (7) Clean_Bases ：质控后的有效数据量，以Gp为单位;<br/>
                    (8) GC (%) ：质控后的有效数据G和C碱基数占总碱基数的比例;<br/>
                    (9) Q30 (%）：大于Q30的碱基占Clean_Bases的比例; <br/>
                    (10) Valid_Bases (%）：Clean_Reads占Raw_Reads的比例;<br/>
                    ",'path'=>"1.QC/report_data_stat.xls");

            #K-mer分析
            $title2Num=0;
            print $output "\n";$writer->emptyTag('h1', 'name'=>++$title1Num." K-mer分析",'type'=>"一级标题显示样式",'desc'=>"一级标题描述");
            #K-mer分析原理
            #K-mer估计基因组大小
            print $output "\n";$writer->emptyTag('h2', 'name'=>$title1Num.".".++$title2Num." K-mer估计基因组大小",'type'=>"二级标题显示样式",'desc'=>"二级标题描述");
            print $output "\n";$writer->emptyTag('p','type'=>'正文段落显示样式','desc'=>"根据测序对基因组大小的估计，可以 ①判断测序得到的基因组与样本的基因组大小预计值是否一致； ②对后续基因组序列的拼接组装所消耗的资源进行预估； ③评估测序深度及测序数据量是否足够<sup>[<a href='#ref2'>2</a>]</sup>。");
            print $output "\n";$writer->emptyTag('p','type'=>'正文段落显示样式','desc'=>"K-mer 是指一段长度为K bp 的序列，对于一个长度为L 的序列，可以得到(L-K+1)条长度为K 的K-mer 片段。");
            print $output "\n";$writer->emptyTag('p','type'=>'正文段落显示样式','desc'=>"对测序得到的reads 分别取K-mer，在read 上每移动一个碱基取一次K-mer，然后统计每个K-mer 出现的频数，K-mer 出现的频数（即深度）理论上是服从泊松分布的。");
            print $output "\n";$writer->emptyTag('p','type'=>'正文段落显示样式','desc'=>"从K-mer 的深度分布图，可以得到深度的峰值，此即为该泊松分布的平均值和方差。基因组估计值可以用下面公式得到：");
            print $output "\n";$writer->emptyTag('pic',name=>"",type=>'img-width-max',desc=>'',path=>"pic/caculated_genesize.png");
            
            #K-mer分析结果
            print $output "\n";$writer->emptyTag('h2', 'name'=>$title1Num.".".++$title2Num." K-mer分析结果",'type'=>"二级标题显示样式",'desc'=>"二级标题描述");
            print $output "\n";$writer->emptyTag('p','type'=>'正文段落显示样式','desc'=>"使用K-mer 分析基因组大小统计情况：");
            print $output "\n";$writer->emptyTag('table','name'=>"表".++$table." K-mer 分析得到的基因组特征统计情况",'type'=>"type1|full",'desc'=>"
            备注：<br/>
            (1)K-mer：分析选取的K-mer的大小。<br/>
            (2)K-mer Depth：K-mer深度，即为泊松分布对应的期望值。<br/>
            (3)K-mer number：采取Soapdenovo软件得到的K-mer总数。<br/>
            (4)Genome size：计算得到的基因组大小，K-mer number/K-mer Depth，以M（兆）为单位。<br/>
            ",'path'=>"2.assembly/Result.xls");
            
            #基因组组装
            $title2Num=0;
            print $output "\n";$writer->emptyTag('h1', 'name'=>++$title1Num." 基因组组装",'type'=>"一级标题显示样式",'desc'=>"一级标题描述");
            #组装原理
            #Contig组装
            print $output "\n";$writer->emptyTag('h2', 'name'=>$title1Num.".".++$title2Num." Contig组装",'type'=>"二级标题显示样式",'desc'=>"二级标题描述");
            print $output "\n";$writer->emptyTag('p','type'=>'正文段落显示样式','desc'=>"将所有小片段库测序得到的reads截断成更小的序列片段，通过它们之间的重叠关系构建de Brujin图；de Brujin图非常复杂，为了简化de Brujin图，需要去掉无法继续连接的分支、低覆盖度的分支（一般这两种情况由于测序错误造成），并且利用reads信息化简重复序列在de Brujin图的分叉通路，对于少量的杂合位点，采用随机选择策略，合并杂合位点；得到一个简化后的de Brujin图，这样的图仍然会有很多分叉位点无法确定真正的连接关系，因此在每个分叉位点将序列截断，得到了最初的contigs。");
            
            #Scaffold组装
            print $output "\n";$writer->emptyTag('h2', 'name'=>$title1Num.".".++$title2Num." Scaffold组装",'type'=>"二级标题显示样式",'desc'=>"二级标题描述");
            print $output "\n";$writer->emptyTag('p','type'=>'正文段落显示样式','desc'=>"将所有文库测序得到的reads比对回初步得到的contigs，利用reads之间的连接关系和插入片段大小信息，将contigs组装成scaffolds。");
            print $output "\n";$writer->emptyTag('pic',name=>"图 ".++$pic." 组装过程示意图",type=>'img-width-max',desc=>'',path=>"pic/assembly_flow.png");
            #组装结果
            #组装结果统计
            print $output "\n";$writer->emptyTag('h2', 'name'=>$title1Num.".".++$title2Num." 组装结果统计",'type'=>"二级标题显示样式",'desc'=>"二级标题描述");
            print $output "\n";$writer->emptyTag('table','name'=>"表".++$table." 组装结果数据统计",'type'=>"type1|full",'desc'=>"注：contigs_num (>=x bp)表示大于等于x bp的contigs数量；Longest contig表示最长contigs长度；Total_len表示组装基因组大小；GC (%)表示组装基因组GC含量；N50表示将各个序列按长度大小排序，从大至小逐一扫描各个序列的长度值，进行累加，当该累加值第一次超过所有序列总长的50%时，此时扫描到的序列，其长度值即为N50值，N75值亦同理；L50表示当累计长度达到N50时contigs的数量，L75值亦同理。",'path'=>"2.assembly/Assembly_report.xls");





            print $output "\n";$writer->emptyTag('file','name'=>"组装scaffold结果：2.assembly/$s.fasta",'type'=>"文件显示样式",'desc'=>"",'path'=>"2.assembly/$s.fasta",'action'=>"文件类型");
}
#        if(-d "$idir/1.Assembly/"){
#            print $output "\n";$writer->emptyTag('h1', 'name'=>++$title1Num." 基因组拼接",'type'=>"一级标题显示样式",'desc'=>"一级标题描述");
#            print $output "\n";$writer->emptyTag('file','name'=>"结果目录：1.Assembly",'type'=>"文件显示样式",'desc'=>"",'path'=>"1.Assembly",'action'=>"文件类型");
#            print $output "\n";$writer->emptyTag('p','type'=>'正文段落显示样式','desc'=>"基于三代数据进行基因组的从头组装，组装原理如下：测序获得的三代测序数据(Subreads)首先使用falcon<sup>[<a href='#ref1'>1</a>]</sup>进行原始数据的自我矫正及基因组初步组装，基于Overlap-Layout-Consensus算法初步组装得到一致性序列(Consensus)，之后采用GenomicConsensus 软件基于其中的quiver算法以原始数据(Subreads)为辅助再次对一致性序列进行矫正，再使用sprai（single pass read accuracy improver）矫正过的原始序列(Corrected Subreads)作为辅助数据将组装出来的一致性序列(Corrected Consensus)进行circulator<sup>[<a href='#ref2'>2</a>]</sup>环化处理，最终得到环化的细菌基因组(Genome)。组装流程如下：");
#            print $output "\n";$writer->emptyTag('pic','name'=>"图 ".++$pic."基因组组装流程",'path'=>"pic/Denovo_script.png");
#            system("cp $idir/1.Assembly/$s.assembly_summary.xls $odir/pic/");
#            print $output "\n";$writer->emptyTag('file','name'=>"组装结果统计：$s.assembly_summary.xls",'type'=>"文件显示样式",'desc'=>"",'path'=>"1.Assembly/$s.assembly_summary.xls",'action'=>"文件类型");
#            print $output "\n";$writer->emptyTag('table','name'=>"表 ".++$table." 组装结果统计",'type'=>"type1|full",'desc'=>"注：contigs_num (>=x bp)表示大于等于x bp的contigs数量；Longest contig表示最长contigs长度；Total_len表示组装基因组大小；GC (%)表示组装基因组GC含量；N50表示将各个序列按长度大小排序，从大至小逐一扫描各个序列的长度值，进行累加，当该累加值第一次超过所有序列总长的50%时，此时扫描到的序列，其长度值即为N50值，N75值亦同理；L50表示当累计长度达到N50时contigs的数量，L75值亦同理。",'path'=>"pic/$s.assembly_summary.xls");
#            print $output "\n";$writer->emptyTag('file','name'=>"组装基因组序列：$s.fasta",'type'=>"文件显示样式",'desc'=>"",'path'=>"1.Assembly/$s.fasta",'action'=>"文件类型");
#        }        
#

		if(-d "$idir/3.Genome_structure/"){
            $title2Num=0;
			print $output "\n";$writer->emptyTag('h1', 'name'=>++$title1Num." 基因组组分分析",'type'=>"一级标题显示样式",'desc'=>"一级标题描述");
			print $output "\n";$writer->emptyTag('p','type'=>'正文段落显示样式','desc'=>"概述：微生物基因组不大，但包含的功能区域非常丰富，比例甚至占基因组大小的90%以上。除编码基因区域，更有非编码区域实现转录调控、转录后调控、翻译调控、表观遗传调控等功能，部分功能区域还与物种进化的多样性存在关系。");
			print $output "\n";$writer->emptyTag('p','type'=>'正文段落显示样式','desc'=>"通过多种方法，对编码基因、重复序列、非编码RNA 等进行预测，获取测序菌株基因组的组成情况。");
			my $totalgff="$idir/3.Genome_structure/$s.gff";
			print $output "\n";$writer->emptyTag('file','name'=>"基因/ncRNA/repeat总gff文件：3.Genome_structure/$s.gff",'type'=>"文件显示样式",'desc'=>"",'path'=>"3.Genome_structure/$s.gff",'action'=>"文件类型");
			if(-d "$idir/3.Genome_structure/Gene_predict"){
				print $output "\n";$writer->emptyTag('h2', 'name'=>$title1Num.".".++$title2Num." 基因预测",'type'=>"二级标题显示样式",'desc'=>"二级标题描述");
				print $output "\n";$writer->emptyTag('file','name'=>"结果目录：3.Genome_structure/Gene_predict",'type'=>"文件显示样式",'desc'=>"",'path'=>"3.Genome_structure/Gene_predict",'action'=>"文件类型");
				print $output "\n";$writer->emptyTag('p','type'=>'正文段落显示样式','desc'=>"在原核生物中，由于基因组结构比真核生物简单，基因紧密连续的分布在基因组上，在启动子区以及转录翻译起始位点附近可以找到有较高保守性的序列模式，因此原核生物基因预测比真核生物相对来说更为简单。");
				print $output "\n";$writer->emptyTag('p','type'=>'正文段落显示样式','desc'=>"使用软件prodigal(v2.6.3)<sup>[<a href='#ref3'>3</a>]</sup>来进行编码基因预测，prodigal的全称是Prokaryotic Dynamic Programming Genefinding Algorithm。");
				#print $output "\n";$writer->emptyTag('p','type'=>'正文段落显示样式','desc'=>"其中$s.ffn为预测基因的碱基序列文件；$s.faa为预测基因的蛋白序列。");
				print $output "\n";$writer->emptyTag('file','name'=>"基因预测统计结果：$s.gene.summary.xls",'type'=>"文件显示样式",'desc'=>"",'path'=>"3.Genome_structure/Gene_predict/$s.gene.summary.xls",'action'=>"文件类型");
				system("cp $idir/3.Genome_structure/Gene_predict/$s.gene.summary.xls $odir/pic/");
				print $output "\n";$writer->emptyTag('table','name'=>"表 ".++$table." 基因预测统计结果",'type'=>"type1|full",'desc'=>"注：Sample样品名；Genome_size(bp)基因组大小；Gene_num(#)基因数量；Total_len(bp)基因总长；Average_len(bp)基因平均长度；GC% (gene region)基因区GC含量；Gene_len/Genome(%)基因总长与基因组长度比值。",'path'=>"pic/$s.gene.summary.xls");
				print $output "\n";$writer->emptyTag('file','name'=>"基因的碱基序列文件：$s.ffn",'type'=>"文件显示样式",'desc'=>"",'path'=>"3.Genome_structure/Gene_predict/$s.ffn",'action'=>"文件类型");
				print $output "\n";$writer->emptyTag('file','name'=>"基因的蛋白序列文件：$s.faa",'type'=>"文件显示样式",'desc'=>"",'path'=>"3.Genome_structure/Gene_predict/$s.faa",'action'=>"文件类型");
			}
			if(-d "$idir/3.Genome_structure/ncRNA"){
				print $output "\n";$writer->emptyTag('h2', 'name'=>$title1Num.".".++$title2Num." 非编码RNA预测",'type'=>"二级标题显示样式",'desc'=>"二级标题描述");
				print $output "\n";$writer->emptyTag('file','name'=>"结果目录：3.Genome_structure/ncRNA",'type'=>"文件显示样式",'desc'=>"",'path'=>"3.Genome_structure/ncRNA",'action'=>"文件类型");
				print $output "\n";$writer->emptyTag('p','type'=>'正文段落显示样式','desc'=>"非编码RNA（ncRNA）是一类广泛存在于细菌、古菌和真核生物生物体内，执行多种生物学功能的RNA分子，但其本身并不携带翻译为蛋白质的信息，终产物是RNA。主要类型包括sRNA、rRNA、tRNA、snRNA以及miRNA等。对于细菌，ncRNA的类型主要指tRNA、rRNA以及sRNA三种，其中：");
				print $output "\n";$writer->emptyTag('p','type'=>'正文段落显示样式','desc'=>"tRNA：具有携带并转运氨基酸功能的一类小分子核糖核酸，由70-90个核苷酸组成，预测软件：tRNAscan-SE(v1.3.1)<sup>[<a href='#ref4'>4</a>]</sup>；");
				print $output "\n";$writer->emptyTag('p','type'=>'正文段落显示样式','desc'=>"rRNA：核糖体RNA。原核生物的rRNA分类为三类： 5S rRNA、16S rRNA和23S rRNA，预测软件：RNAmmer(v1.2)<sup>[<a href='#ref5'>5</a>]</sup>；");
				print $output "\n";$writer->emptyTag('p','type'=>'正文段落显示样式','desc'=>"sRNA：细菌中，长度在50~500nt的ncRNA通常定义为小RNA（sRNA）,包括一些Cis-reg（顺式作用元件）、snoRNA等等。预测软件：Rfam(v10.0)<sup>[<a href='#ref6'>6</a>]</sup>。");
				#print $output "\n";$writer->emptyTag('p','type'=>'正文段落显示样式','desc'=>"其中sample.rRNA.fna为rRNA序列；sample.sRNA.fna为sRNA序列；sample.tRNA.structure为tRNA序列和结构文件");
				print $output "\n";$writer->emptyTag('file','name'=>"tRNA结构和序列：$s.tRNA.structure",'type'=>"文件显示样式",'desc'=>"",'path'=>"3.Genome_structure/ncRNA/$s.tRNA.structure",'action'=>"文件类型");
				print $output "\n";$writer->emptyTag('file','name'=>"rRNA序列：$s.rRNA.fna",'type'=>"文件显示样式",'desc'=>"",'path'=>"3.Genome_structure/ncRNA/$s.rRNA.fna",'action'=>"文件类型");
				print $output "\n";$writer->emptyTag('file','name'=>"sRNA序列：$s.sRNA.fna",'type'=>"文件显示样式",'desc'=>"",'path'=>"3.Genome_structure/ncRNA/$s.sRNA.fna",'action'=>"文件类型");
			}
			if(-d "$idir/3.Genome_structure/Repeat_sequence"){
				print $output "\n";$writer->emptyTag('h2', 'name'=>$title1Num.".".++$title2Num." 重复序列预测",'type'=>"二级标题显示样式",'desc'=>"二级标题描述");
				print $output "\n";$writer->emptyTag('file','name'=>"结果目录：3.Genome_structure/Repeat_sequence",'type'=>"文件显示样式",'desc'=>"",'path'=>"3.Genome_structure/Repeat_sequence",'action'=>"文件类型");
				print $output "\n";$writer->emptyTag('p','type'=>'正文段落显示样式','desc'=>"重复序列可分为串联重复序列（Tandem repeat）和散在重复序列（Interspersed repeat）两大类。其中串联重复序列包括有微卫星序列、小卫星序列和卫星DNA等；散在重复序列又称转座子元件，包括以DNA-DNA方式转座的DNA转座子和反转录转座子。常见的反转录转座子类别有LTR、LINE和SINE等。");
				print $output "\n";$writer->emptyTag('p','type'=>'正文段落显示样式','desc'=>"通过RepeatMasker(v4.0.7)<sup>[<a href='#ref7'>7</a>]</sup>软件进行重复序列预测。");
				#print $output "\n";$writer->emptyTag('p','type'=>'正文段落显示样式','desc'=>"其中sample.masked.fna为屏蔽掉重复序列的基因组；sample.repeat.gff为重复序列预测结果；sample.repeat.xls为重复序列统计信息。");
				print $output "\n";$writer->emptyTag('file','name'=>"重复序列统计文件：$s.repeat.xls",'type'=>"文件显示样式",'desc'=>"",'path'=>"3.Genome_structure/Repeat_sequence/$s.repeat.xls",'action'=>"文件类型");
				system("awk -F \"\\t\" \'{if(\$0~/^##/ || \$0~/^\\t/ || \$0~/^--/ || \$0==\"\"){}else{print \$1\"\\t\"\$3\"\\t\"\$4\"\\t\"\$5}}\' $idir/3.Genome_structure/Repeat_sequence/$s.repeat.xls|sed \'s/://g\' > $odir/pic/$s.repeat.xls");
				print $output "\n";$writer->emptyTag('table','name'=>"表 ".++$table." 重复序列统计文件",'type'=>"type1|full",'desc'=>"注：sample name样品名；sequences num序列数量；total length基因组总长；GC%表示GC含量百分比；bases masked屏蔽重复序列长度；class重复序列类别；family重复序列家族；number对应重复类型数量；length(bp)对应重复类型总长度；percentage(%)对应重复类型百分比。重复序列类型：SINE：短散在重复序列，LINE：长散在重复序列，LTR elements：长末端重复序列，DNA elements：DNA转座子，Unclassified：未知类别散在重复序列，Total interspersed repeats：全部散在重复序列，Small RNA：小RNA，Satellites：卫星，Simple repeats：简单重复，Low complexity：低复杂度序列。",'path'=>"pic/$s.repeat.xls");
				print $output "\n";$writer->emptyTag('file','name'=>"屏蔽掉重复序列的基因组：$s.masked.fna",'type'=>"文件显示样式",'desc'=>"",'path'=>"3.Genome_structure/Repeat_sequence/$s.masked.fna",'action'=>"文件类型");
			}
			if(-d "$idir/3.Genome_structure/CRISPR"){
				print $output "\n";$writer->emptyTag('h2', 'name'=>$title1Num.".".++$title2Num." CRISPR预测",'type'=>"二级标题显示样式",'desc'=>"二级标题描述");
				print $output "\n";$writer->emptyTag('file','name'=>"结果目录：3.Genome_structure/CRISPR",'type'=>"文件显示样式",'desc'=>"",'path'=>"3.Genome_structure/CRISPR",'action'=>"文件类型");
				print $output "\n";$writer->emptyTag('p','type'=>'正文段落显示样式','desc'=>"CRISPR是一串包含多个短而重复的序列的碱基序列，重复序列之间是一些约30bp的\"spacer DNA\"。在原核生物中，CRISPR起到免疫系统的作用，对外来的质粒和噬菌体序列具有抵抗作用。CRISPR能识别并使入侵的功能元件沉默。在这里分别使用PILER-CR v1.06<sup>[<a href='#ref8'>8</a>]</sup>和CRT1.2-CLI<sup>[<a href='#ref9'>9</a>]</sup>对基因组进行CRISPR预测。");
				print $output "\n";$writer->emptyTag('file','name'=>"CRISPR预测统计结果：$s.CRISPR.summary.xls",'type'=>"文件显示样式",'desc'=>"",'path'=>"3.Genome_structure/CRISPR/$s.CRISPR.summary.xls",'action'=>"文件类型");
				system("cp $idir/3.Genome_structure/CRISPR/$s.CRISPR.summary.xls $odir/pic/");
				print $output "\n";$writer->emptyTag('table','name'=>"表 ".++$table." CRISPR预测统计结果",'type'=>"type1|full",'desc'=>"注：Sample样品名；PILER-CR(#)预测到CRISPR的数量；CRT(#)预测到CRISPR的数量。",'path'=>"pic/$s.CRISPR.summary.xls");
				print $output "\n";$writer->emptyTag('table','name'=>"表 ".++$table." PILER-CR v1.06预测结果说明",'type'=>"type1|full",'desc'=>"",'path'=>"pic/header.sample.pilercr.xls");
				print $output "\n";$writer->emptyTag('table','name'=>"表 ".++$table." CRT1.2-CLI预测结果说明",'type'=>"type1|full",'desc'=>"",'path'=>"pic/header.sample.crt.xls");
			}
			if(-d "$idir/3.Genome_structure/Prophage"){
				print $output "\n";$writer->emptyTag('h2', 'name'=>$title1Num.".".++$title2Num." Prophage预测",'type'=>"二级标题显示样式",'desc'=>"二级标题描述");
				print $output "\n";$writer->emptyTag('file','name'=>"结果目录：3.Genome_structure/Prophage",'type'=>"文件显示样式",'desc'=>"",'path'=>"3.Genome_structure/Prophage",'action'=>"文件类型");
				print $output "\n";$writer->emptyTag('p','type'=>'正文段落显示样式','desc'=>"整合在宿主基因组上的温和噬菌体的核酸称之为前噬菌体。带有前噬菌体的菌称为溶源菌，它们具有无需由外部感染而可产生噬菌体的遗传能力，并且这种能力可传递给后代。前噬菌体序列的存在可能也会允许一些细菌获取抗生素抗性，增强对环境的适应性，提高粘附力或使细菌成为致病菌。同时，通过前噬菌体的研究可能找到特异的抗生素甚至是先进的癌症治疗方法。");
				print $output "\n";$writer->emptyTag('p','type'=>'正文段落显示样式','desc'=>"通过软件PhiSpy(v2.3)<sup>[<a href='#ref10'>10</a>]</sup>预测前噬菌体。");
				print $output "\n";$writer->emptyTag('file','name'=>"Prophage预测结果：$s.prophage.xls",'type'=>"文件显示样式",'desc'=>"",'path'=>"3.Genome_structure/Prophage/$s.prophage.xls",'action'=>"文件类型");
				system("cp $idir/3.Genome_structure/Prophage/$s.prophage.xls $odir/pic/");
				print $output "\n";$writer->emptyTag('table','name'=>"表 ".++$table." Prophage预测结果",'type'=>"type1|full",'desc'=>"注：Prophage ID前噬菌体ID；Chr所在染色体；Start(bp)起始位置；End(bp)结束位置。",'path'=>"pic/$s.prophage.xls");
			}
		}
		if(-d "$idir/4.Basic_function"){
			$title2Num=0;
			print $output "\n";$writer->emptyTag('h1', 'name'=>++$title1Num." 常见功能数据库注释",'type'=>"一级标题显示样式",'desc'=>"一级标题描述");
			print $output "\n";$writer->emptyTag('p','type'=>'正文段落显示样式','desc'=>"常见数据库注释包括NR注释、COG/KOG功能注释、GO分类、Swissprot、eggNOG、KEGG和Pfam。对于NR、COG/KOG、GO、Swissprot、eggNOG、KEGG数据库注释，我们采用diamond<sup>[<a href='#ref10'>10</a>]</sup>软件进行比对，取e<1e-5的注释，筛选具有最高序列相似性的蛋白，从而得到功能注释信息。Pfam数据库是一系列蛋白质家族的集合，其中每一个蛋白家族都以多序列比对和隐马尔科夫模型的形式来表示。我们采用HMMER<sup>[<a href='#ref12'>12</a>]</sup>软件和蛋白家族模型比对，从而筛选出得分最高的家族。");
			print $output "\n";$writer->emptyTag('file','name'=>"各数据库注释比例统计表：$s.generalDB.stat.xls",'type'=>"文件显示样式",'desc'=>"",'path'=>"4.Basic_function/$s.generalDB.stat.xls",'action'=>"文件类型");
			system("cp $idir/4.Basic_function/$s.generalDB.stat.xls $odir/pic/");
			print $output "\n";$writer->emptyTag('table','name'=>"表 ".++$table." 各数据库注释统计表",'type'=>"type1|full",'desc'=>"注：Anno_Database：数据库；Annotated_Number：注释到对应数据库的基因数量（占总基因百分比）；300<=length<1000：注释到对应数据库的满足基因长度>=300bp且<1000bp的基因数量（占总基因百分比）；length>=1000：注释到对应数据库的长度>=1000bp的基因数量（占总基因百分比）。",'path'=>"pic/$s.generalDB.stat.xls");
			chomp(my $db_num=`wc -l $idir/4.Basic_function/$s.generalDB.stat.xls | awk '{print $1}'`); $db_num=$db_num-1;
			my $picdes=$db_num <= 5 ? "各个不同颜色的环代表各数据库，图上的数字代表基因数量。" : "上方条形图上的数字代表下面矩阵中对应加黑点的数据库交集的结果，左边的柱子代表各数据库全部注释到基因数量。";
			#system("cp $idir/4.Basic_function/$s.VennGraph.png $odir/pic/");
			#print $output "\n";$writer->emptyTag('file','name'=>"各数据库注释维恩图: $s.VennGraph.pdf",'type'=>"文件显示样式",'desc'=>"",'path'=>"4.Basic_function/$s.VennGraph.pdf",'action'=>"文件类型");
			#print $output "\n";$writer->emptyTag('pic','name'=>"图 ".++$pic." 各数据库注释维恩图",'desc'=>"注：$picdes",'path'=>"pic/$s.VennGraph.png");
			print $output "\n";$writer->emptyTag('file','name'=>"各数据库注释总表：$s.generalDB.anno.xls",'type'=>"文件显示样式",'desc'=>"",'path'=>"4.Basic_function/$s.generalDB.anno.xls",'action'=>"文件类型");
			print $output "\n";$writer->emptyTag('table','name'=>"表 ".++$table." 注释总表表头说明",'type'=>"type1|full",'desc'=>"",'path'=>"pic/annotationHeader.txt");
			if(-d "$idir/4.Basic_function/NR/"){
				print $output "\n";$writer->emptyTag('h2', 'name'=>$title1Num.".".++$title2Num." NR数据库注释",'type'=>"二级标题显示样式",'desc'=>"二级标题描述");
				print $output "\n";$writer->emptyTag('file','name'=>"结果目录：4.Basic_function/NR",'type'=>"文件显示样式",'desc'=>"",'path'=>"4.Basic_function/NR",'action'=>"文件类型");
				print $output "\n";$writer->emptyTag('p','type'=>'正文段落显示样式','desc'=>"NR是non-redundant的缩写，即非冗余数据库。它剔除了冗余序列最主要的来源（EST、STS、GSS、HTGS）。NR数据库含有整个NCBI GenBank保存的DNA序列所翻译的蛋白质序列，该数据库信息比较全面，但是大多数据未经过人工校验。数据库链接：https://www.ncbi.nlm.nih.gov/。");
					
				print $output "\n";$writer->emptyTag('file','name'=>"NR库Top10物种分布图：$s.NR.species.top10.pdf",'type'=>"文件显示样式",'desc'=>"",'path'=>"4.Basic_function/NR/$s.NR.species.top10.pdf",'action'=>"文件类型");
				system("cp $idir/4.Basic_function/NR/$s.NR.species.top10.png $odir/pic/");
				print $output "\n";$writer->emptyTag('pic','name'=>"图 ".++$pic." NR库Top10物种分布图",'desc'=>"注：不同颜色的扇形代表注释到的各物种占比。",'path'=>"pic/$s.NR.species.top10.png");
				print $output "\n";$writer->emptyTag('file','name'=>"NR库注释结果：$s.NR.blast.best.anno.xls",'type'=>"文件显示样式",'desc'=>"",'path'=>"4.Basic_function/NR/$s.NR.blast.best.anno.xls",'action'=>"文件类型");
				#print $output "\n";$writer->emptyTag('table','name'=>"表 ".++$table." NR库注释表头说明(sample.NR.blast.best.anno.xls)",'type'=>"type1|full",'desc'=>"",'path'=>"pic/nr_swissprot_header.txt");
			}
			if(-d "$idir/4.Basic_function/Swissprot/"){
				print $output "\n";$writer->emptyTag('h2', 'name'=>$title1Num.".".++$title2Num." Swissprot数据库注释",'type'=>"二级标题显示样式",'desc'=>"二级标题描述");
				print $output "\n";$writer->emptyTag('file','name'=>"结果目录：4.Basic_function/Swissprot",'type'=>"文件显示样式",'desc'=>"",'path'=>"4.Basic_function/Swissprot",'action'=>"文件类型");
				print $output "\n";$writer->emptyTag('p','type'=>'正文段落显示样式','desc'=>"Swiss-Prot，是2002 年由UniProt consortium 建立的基因数据库，其特点在注释结果经过实验验证，可靠性较高，可用作其他数据的参考。数据库链接：http://www.uniprot.org/。");
				print $output "\n";$writer->emptyTag('file','name'=>"Swissprot库注释结果：$s.Swissprot.blast.best.anno.xls",'type'=>"文件显示样式",'desc'=>"",'path'=>"4.Basic_function/Swissprot/$s.Swissprot.blast.best.anno.xls",'action'=>"文件类型");
				#print $output "\n";$writer->emptyTag('table','name'=>"表 ".++$table." Swissprot库注释表头说明(sample.Swissprot.blast.best.anno.xls)",'type'=>"type1|full",'desc'=>"",'path'=>"pic/nr_swissprot_header.txt");
			}
			if(-d "$idir/4.Basic_function/KEGG/"){
				print $output "\n";$writer->emptyTag('h2', 'name'=>$title1Num.".".++$title2Num." KEGG数据库注释",'type'=>"二级标题显示样式",'desc'=>"二级标题描述");
				print $output "\n";$writer->emptyTag('file','name'=>"结果目录：4.Basic_function/KEGG",'type'=>"文件显示样式",'desc'=>"",'path'=>"4.Basic_function/KEGG",'action'=>"文件类型");
				print $output "\n";$writer->emptyTag('p','type'=>'正文段落显示样式','desc'=>"KEGG 京都基因与基因组百科全书，核心数据库为KEGG PATHWAY 数据库。该数据库将生物通路划分为八大类，每一大类下还有细分，以通路图的方式展示出来。通过该数据库注释，可以快速地查询行使某一功能相关的所有注释基因。数据库链接：http://www.genome.jp/kegg/pathway.html。");
				print $output "\n";$writer->emptyTag('p','type'=>'正文段落显示样式','desc'=>"首先和KEGG数据库比对，然后利用kobas<sup>[<a href='#ref13'>13</a>]</sup>软件（http://kobas.cbi.pku.edu.cn/kobas-2.1.1/kobas-3.0.3.tar.gz）将比对结果转化为注释上的KO和通路。");
				
				print $output "\n";$writer->emptyTag('file','name'=>"Level2水平下的KEGG注释统计图：$s.KEGG.classification.pdf",'type'=>"文件显示样式",'desc'=>"",'path'=>"4.Basic_function/KEGG/$s.KEGG.classification.pdf",'action'=>"文件类型");
				system("cp $idir/4.Basic_function/KEGG/$s.KEGG.classification.png $odir/pic/");
				print $output "\n";$writer->emptyTag('pic','name'=>"图 ".++$pic." Level2水平下的KEGG注释统计图",'desc'=>"注：横轴表示gene数量，纵轴表示Level2 pathway通路名，柱子右边数字代表注释到该Level2 pathway下的gene数量。",'path'=>"pic/$s.KEGG.classification.png");
				print $output "\n";$writer->emptyTag('file','name'=>"KEGG库注释表格：$s.KEGG.gene.anno.xls",'type'=>"文件显示样式",'desc'=>"",'path'=>"4.Basic_function/KEGG/$s.KEGG.gene.anno.xls",'action'=>"文件类型");
				#print $output "\n";$writer->emptyTag('table','name'=>"表 ".++$table." sample.KEGG.gene.anno.xls表头说明",'type'=>"type1|full",'desc'=>"",'path'=>"pic/KEGG.gene.anno.header.txt");
				print $output "\n";$writer->emptyTag('file','name'=>"KEGG库注释通路统计表格：$s.KEGG.pathway.xls",'type'=>"文件显示样式",'desc'=>"",'path'=>"4.Basic_function/KEGG/$s.KEGG.pathway.xls",'action'=>"文件类型");
				#print $output "\n";$writer->emptyTag('table','name'=>"表 ".++$table." sample.KEGG.pathway.xls表头说明",'type'=>"type1|full",'desc'=>"",'path'=>"pic/KEGG.pathway.header.txt");
				
				my @png=glob("$idir/4.Basic_function/KEGG/$s\_pathway/*png");
				system("cp $png[0] $odir/pic/");
				my $png=$png[0];$png=basename($png);$png=~s/\.png$//g;
				print $output "\n";$writer->emptyTag('file','name'=>"KEGG通路图目录：$s\_pathway",'type'=>"文件显示样式",'desc'=>"",'path'=>"4.Basic_function/KEGG/$s\_pathway",'action'=>"文件类型");
				print $output "\n";$writer->emptyTag('pic','name'=>"图 ".++$pic." KEGG通路图（以$png为例）",'desc'=>"注：注释上的酶标记为浅绿色。详细说明请参见：http://www.genome.jp/kegg/document/help_pathway.html",'path'=>"pic/$png\.png");
			}
			if(-d "$idir/4.Basic_function/COG/"){
				print $output "\n";$writer->emptyTag('h2', 'name'=>$title1Num.".".++$title2Num." COG数据库注释",'type'=>"二级标题显示样式",'desc'=>"二级标题描述");
				print $output "\n";$writer->emptyTag('file','name'=>"结果目录：4.Basic_function/COG",'type'=>"文件显示样式",'desc'=>"",'path'=>"4.Basic_function/COG",'action'=>"文件类型");
				print $output "\n";$writer->emptyTag('p','type'=>'正文段落显示样式','desc'=>"Cluster of Orthologous Groups of proteins，即同源蛋白簇。该数据库分为COG（原核生物）和KOG（真核生物），是由NCBI创建并维护的蛋白数据库，根据细菌、藻类和真核生物完整基因组的编码蛋白系统进化关系分类构建而成。通过比对可以将某个蛋白序列注释到某一个COG中，每一簇COG由直系同源序列构成，从而可以推测该序列的功能，数据库链接：https://www.ncbi.nlm.nih.gov/COG/。");
				print $output "\n";$writer->emptyTag('file','name'=>"COG功能分类图：$s.COG.cluster.pdf",'type'=>"文件显示样式",'desc'=>"",'path'=>"4.Basic_function/COG/$s.COG.cluster.pdf",'action'=>"文件类型");
				system("cp $idir/4.Basic_function/COG/$s.COG.cluster.png $odir/pic/");
				print $output "\n";$writer->emptyTag('pic','name'=>"图 ".++$pic." COG功能分类图",'desc'=>"注：横轴表示COG功能分类，纵轴表示gene数量。",'path'=>"pic/$s.COG.cluster.png");
				print $output "\n";$writer->emptyTag('file','name'=>"COG注释结果：$s.COG.gene.anno.xls",'type'=>"文件显示样式",'desc'=>"",'path'=>"4.Basic_function/COG/$s.COG.gene.anno.xls",'action'=>"文件类型");
				#print $output "\n";$writer->emptyTag('table','name'=>"表 ".++$table." $s.COG.gene.anno.xls表头说明",'type'=>"type1|full",'desc'=>"",'path'=>"pic/COG.gene.anno.header.txt");
				print $output "\n";$writer->emptyTag('file','name'=>"COG分类水平统计结果：$s.COG.class.xls",'type'=>"文件显示样式",'desc'=>"",'path'=>"4.Basic_function/COG/$s.COG.class.xls",'action'=>"文件类型");
				#print $output "\n";$writer->emptyTag('table','name'=>"表 ".++$table." $s.COG.class.xls表头说明",'type'=>"type1|full",'desc'=>"",'path'=>"pic/COG.class.header.txt");
			}
			if(-d "$idir/4.Basic_function/KOG/"){
				print $output "\n";$writer->emptyTag('h2', 'name'=>$title1Num.".".++$title2Num." KOG数据库注释",'type'=>"二级标题显示样式",'desc'=>"二级标题描述");
				print $output "\n";$writer->emptyTag('file','name'=>"结果目录：4.Basic_function/KOG",'type'=>"文件显示样式",'desc'=>"",'path'=>"4.Basic_function/KOG",'action'=>"文件类型");
				print $output "\n";$writer->emptyTag('p','type'=>'正文段落显示样式','desc'=>"Cluster of Orthologous Groups of proteins，即同源蛋白簇。该数据库分为COG（原核生物）和KOG（真核生物），是由NCBI创建并维护的蛋白数据库，根据细菌、藻类和真核生物完整基因组的编码蛋白系统进化关系分类构建而成。通过比对可以将某个蛋白序列注释到某一个KOG中，每一簇KOG由直系同源序列构成，从而可以推测该序列的功能，数据库链接：https://www.ncbi.nlm.nih.gov/COG/。");
				print $output "\n";$writer->emptyTag('file','name'=>"KOG功能分类图：$s.KOG.cluster.pdf",'type'=>"文件显示样式",'desc'=>"",'path'=>"4.Basic_function/KOG/$s.KOG.cluster.pdf",'action'=>"文件类型");
				system("cp $idir/4.Basic_function/KOG/$s.KOG.cluster.png $odir/pic/");
				print $output "\n";$writer->emptyTag('pic','name'=>"图 ".++$pic." KOG功能分类图",'desc'=>"注：横轴表示KOG功能分类，纵轴表示gene数量。",'path'=>"pic/$s.KOG.cluster.png");
				print $output "\n";$writer->emptyTag('file','name'=>"KOG注释结果：$s.KOG.gene.anno.xls",'type'=>"文件显示样式",'desc'=>"",'path'=>"4.Basic_function/KOG/$s.KOG.gene.anno.xls",'action'=>"文件类型");
				#print $output "\n";$writer->emptyTag('p','type'=>'正文段落显示样式','desc'=>"KOG库注释表格说明：");
				#print $output "\n";$writer->emptyTag('table','name'=>"表 ".++$table." $s.KOG.gene.anno.xls表头说明",'type'=>"type1|full",'desc'=>"",'path'=>"pic/KOG.gene.anno.header.txt");
				print $output "\n";$writer->emptyTag('file','name'=>"KOG分类水平统计结果：$s.KOG.class.xls",'type'=>"文件显示样式",'desc'=>"",'path'=>"4.Basic_function/KOG/$s.KOG.class.xls",'action'=>"文件类型");
				#print $output "\n";$writer->emptyTag('table','name'=>"表 ".++$table." $s.KOG.class.xls表头说明",'type'=>"type1|full",'desc'=>"",'path'=>"pic/KOG.class.header.txt");
			}
			if(-d "$idir/4.Basic_function/eggNOG/"){
				print $output "\n";$writer->emptyTag('h2', 'name'=>$title1Num.".".++$title2Num." eggNOG数据库注释",'type'=>"二级标题显示样式",'desc'=>"二级标题描述");
				print $output "\n";$writer->emptyTag('file','name'=>"结果目录：4.Basic_function/eggNOG",'type'=>"文件显示样式",'desc'=>"",'path'=>"4.Basic_function/eggNOG",'action'=>"文件类型");
				print $output "\n";$writer->emptyTag('p','type'=>'正文段落显示样式','desc'=>"eggNOG（evolutionary genealogy of genes: Non-supervised Orthologous Groups）数据库是国际上普遍认可的同源聚类基因群的专业注释数据库，对NCBI COG/KOG数据库进行了扩展，囊括了真核和原核生物信息。数据库链接：http://eggnog.embl.de/");
				print $output "\n";$writer->emptyTag('file','name'=>"eggNOG功能分类图：$s.eggNOG.barplot.pdf",'type'=>"文件显示样式",'desc'=>"",'path'=>"4.Basic_function/eggNOG/$s.eggNOG.barplot.pdf",'action'=>"文件类型");
				system("cp $idir/4.Basic_function/eggNOG/$s.eggNOG.barplot.png $odir/pic/");
				print $output "\n";$writer->emptyTag('pic','name'=>"图 ".++$pic." eggNOG功能分类图",'desc'=>"注：横轴表示eggNOG功能分类，纵轴表示gene数量。",'path'=>"pic/$s.eggNOG.barplot.png");
				print $output "\n";$writer->emptyTag('file','name'=>"eggNOG注释结果：$s.eggNOG.anno.xls",'type'=>"文件显示样式",'desc'=>"",'path'=>"4.Basic_function/eggNOG/$s.eggNOG.anno.xls",'action'=>"文件类型");
				#print $output "\n";$writer->emptyTag('table','name'=>"表 ".++$table." $s.eggNOG.anno.xls表头说明",'type'=>"type1|full",'desc'=>"",'path'=>"pic/eggNOG.anno.header.txt");
			}
			if(-d "$idir/4.Basic_function/GO/"){
				print $output "\n";$writer->emptyTag('h2', 'name'=>$title1Num.".".++$title2Num." GO数据库注释",'type'=>"二级标题显示样式",'desc'=>"二级标题描述");
				print $output "\n";$writer->emptyTag('file','name'=>"结果目录：4.Basic_function/GO",'type'=>"文件显示样式",'desc'=>"",'path'=>"4.Basic_function/GO",'action'=>"文件类型");
				print $output "\n";$writer->emptyTag('p','type'=>'正文段落显示样式','desc'=>"GO全称是Gene Ontology，GO总共有三个本体，分别为细胞学组件（Cellular Component）、分子功能（Molecular Function）、生物学途径（Biological Process）。三个本体下又细分为约64个小类。GO的最小单位是term（词条、节点），每个term都对应一个或多个类。此分析是基于Swissprot的结果，将蛋白ID映射到GO term。");
				
				print $output "\n";$writer->emptyTag('file','name'=>"GO功能分类图：$s.GO.classification.stat.pdf",'type'=>"文件显示样式",'desc'=>"",'path'=>"4.Basic_function/GO/$s.GO.classification.stat.pdf",'action'=>"文件类型");
				system("cp $idir/4.Basic_function/GO/$s.GO.classification.stat.png $odir/pic/");
				print $output "\n";$writer->emptyTag('pic','name'=>"图 ".++$pic." GO功能分类图",'desc'=>"注：横轴表示GO功能分类，左边纵轴表示注释到该类gene数量占比，右边纵轴表示注释到该类gene数量。",'path'=>"pic/$s.GO.classification.stat.png");
				print $output "\n";$writer->emptyTag('file','name'=>"GO注释结果：$s.GO.gene.anno.xls",'type'=>"文件显示样式",'desc'=>"",'path'=>"4.Basic_function/GO/$s.GO.gene.anno.xls",'action'=>"文件类型");
				#print $output "\n";$writer->emptyTag('p','type'=>'正文段落显示样式','desc'=>"GO库注释表格说明：");
				#print $output "\n";$writer->emptyTag('table','name'=>"表 ".++$table." $s.GO.gene.anno.xls表头说明",'type'=>"type1|full",'desc'=>"",'path'=>"pic/GO.gene.anno.header.txt");
				print $output "\n";$writer->emptyTag('file','name'=>"GO分类水平统计表：$s.GO.class.stat.xls",'type'=>"文件显示样式",'desc'=>"",'path'=>"4.Basic_function/GO/$s.GO.class.stat.xls",'action'=>"文件类型");
				#print $output "\n";$writer->emptyTag('table','name'=>"表 ".++$table." $s.GO.class.stat.xls表头说明",'type'=>"type1|full",'desc'=>"",'path'=>"pic/GO.class.stat.header.txt");
			}
			if(-d "$idir/4.Basic_function/Pfam/"){
				print $output "\n";$writer->emptyTag('h2', 'name'=>$title1Num.".".++$title2Num." Pfam数据库注释",'type'=>"二级标题显示样式",'desc'=>"二级标题描述");
				print $output "\n";$writer->emptyTag('file','name'=>"结果目录：4.Basic_function/Pfam",'type'=>"文件显示样式",'desc'=>"",'path'=>"4.Basic_function/Pfam",'action'=>"文件类型");
				print $output "\n";$writer->emptyTag('p','type'=>'正文段落显示样式','desc'=>"Pfam数据库是一系列蛋白质家族的集合，其中每一个蛋白家族都以多序列比对和隐马尔科夫模型的形式来表示。该分析利用HMMER3<sup>[<a href='#ref14'>14</a>]</sup>软件与蛋白家族模型比对完成注释。数据库链接：http://pfam.xfam.org/。");
				print $output "\n";$writer->emptyTag('file','name'=>"Pfam库注释结果：$s.Pfam.align.anno.xls",'type'=>"文件显示样式",'desc'=>"",'path'=>"4.Basic_function/Pfam/$s.Pfam.align.anno.xls",'action'=>"文件类型");
				#print $output "\n";$writer->emptyTag('table','name'=>"表 ".++$table." $s.Pfam.align.anno.xls表头说明",'type'=>"type1|full",'desc'=>"",'path'=>"pic/pfam_anno_header.txt");
			}
		}
		
		if(-d "$idir/5.Specialized_function"){
			print $output "\n";$writer->emptyTag('h1', 'name'=>++$title1Num." 专业功能数据库注释",'type'=>"一级标题显示样式",'desc'=>"一级标题描述");
			$title2Num=0;
			if(-d "$idir/5.Specialized_function/CARD/"){
				print $output "\n";$writer->emptyTag('h2', 'name'=>$title1Num.".".++$title2Num." 耐药基因（CARD）数据库注释",'type'=>"二级标题显示样式",'desc'=>"二级标题描述");
				print $output "\n";$writer->emptyTag('file','name'=>"结果目录：5.Specialized_function/CARD",'type'=>"文件显示样式",'desc'=>"",'path'=>"5.Specialized_function/CARD",'action'=>"文件类型");
				print $output "\n";$writer->emptyTag('p','type'=>'正文段落显示样式','desc'=>"CARD 抗性基因数据库<sup>[<a href='#ref14'>14</a>]</sup>，ARDB自从09年上线以来鲜有更新，CARD整合了ARDB的全部数据，搭建了一个基于志愿者贡献的数据共享平台，以ARO（Antibiotic Resistance Ontology）为核心对抗性数据进行组织，以达到数据实时更新的效果。利用blast<sup>[<a href='#ref16'>16</a>]</sup>与数据库比对，取e<1e-10的注释。数据库链接：https://card.mcmaster.ca/");
				print $output "\n";$writer->emptyTag('file','name'=>"CARD数据库注释top10 term展示图：$s.card.pdf",'type'=>"文件显示样式",'desc'=>"",'path'=>"5.Specialized_function/CARD/$s.card.pdf",'action'=>"文件类型");
				system("cp $idir/5.Specialized_function/CARD/$s.card.png $odir/pic/");
				print $output "\n";$writer->emptyTag('pic','name'=>"图 ".++$pic." CARD数据库注释top10 term展示图",'desc'=>"注：不同颜色的扇形大小代表注释到该term的基因数量占比，左上角图例后面的数字是注释到该term的gene数量。",'path'=>"pic/$s.card.png");
				print $output "\n";$writer->emptyTag('file','name'=>"CARD数据库注释结果：$s.card.anno.xls",'type'=>"文件显示样式",'desc'=>"",'path'=>"5.Specialized_function/CARD/$s.card.anno.xls",'action'=>"文件类型");
				print $output "\n";$writer->emptyTag('table','name'=>"表 ".++$table." $s.card.anno.xls表头说明",'type'=>"type1|full",'desc'=>"",'path'=>"pic/header.card.anno.xls");
			}
			if(-d "$idir/5.Specialized_function/CAZy/"){
				print $output "\n";$writer->emptyTag('h2', 'name'=>$title1Num.".".++$title2Num." 碳水化合物相关酶（CAZy）数据库注释",'type'=>"二级标题显示样式",'desc'=>"二级标题描述");
				print $output "\n";$writer->emptyTag('file','name'=>"结果目录：5.Specialized_function/CAZy",'type'=>"文件显示样式",'desc'=>"",'path'=>"5.Specialized_function/CAZy",'action'=>"文件类型");
				print $output "\n";$writer->emptyTag('p','type'=>'正文段落显示样式','desc'=>"CAZy 全称为Carbohydrate-Active enZYmes Database<sup>[<a href='#ref16'>16</a>]</sup>，碳水化合物酶相关的专业数据库，内容包括能催化碳水化合物降解、修饰、以及生物合成的相关酶系家族。其包含四个主要分类：糖苷水解酶（Glycoside Hydrolases, GHs）、糖基转移酶（Glycosyl Transferases, GTs）、多糖裂解酶（Polysaccharide Lyases, PLs）和糖类酯解酶（Carbohydrate Esterases, CEs）。此外，还包含与碳水化合物结合结构域（Carbohydrate-Binding Modules, CBMs）。其下又细分为若干小的家族。利用HMMER<sup>[<a href='#ref11'>11</a>]</sup>软件与数据库比对，取e<1e-10的注释。数据库链接：http://www.cazy.org/");
				print $output "\n";$writer->emptyTag('file','name'=>"CAZy注释class分布图：$s.CAZy.class.pdf",'type'=>"文件显示样式",'desc'=>"",'path'=>"5.Specialized_function/CAZy/$s.CAZy.class.pdf",'action'=>"文件类型");
				system("cp $idir/5.Specialized_function/CAZy/$s.CAZy.class.png $odir/pic/");
				print $output "\n";$writer->emptyTag('pic','name'=>"图 ".++$pic." CAZy注释class分布图",'desc'=>"注：横坐标是CAZy分类，纵坐标是注释到对应分类的基因数量。",'path'=>"pic/$s.CAZy.class.png");
				print $output "\n";$writer->emptyTag('file','name'=>"CAZy注释结果：$s.CAZy.anno.xls",'type'=>"文件显示样式",'desc'=>"",'path'=>"5.Specialized_function/CAZy/$s.CAZy.anno.xls",'action'=>"文件类型");
				print $output "\n";$writer->emptyTag('table','name'=>"表 ".++$table." $s.CAZy.anno.xls表头说明",'type'=>"type1|full",'desc'=>"",'path'=>"pic/header.sample.CAZy.anno.xls");
				print $output "\n";$writer->emptyTag('file','name'=>"CAZy注释class水平统计表：$s.CAZy.class.xls",'type'=>"文件显示样式",'desc'=>"",'path'=>"5.Specialized_function/CAZy/$s.CAZy.anno.xls",'action'=>"文件类型");
				print $output "\n";$writer->emptyTag('table','name'=>"表 ".++$table." $s.CAZy.class.xls表头说明",'type'=>"type1|full",'desc'=>"",'path'=>"pic/header.sample.CAZy.class.xls");
				print $output "\n";$writer->emptyTag('file','name'=>"CAZy注释family水平统计表：$s.CAZy.family.xls",'type'=>"文件显示样式",'desc'=>"",'path'=>"5.Specialized_function/CAZy/$s.CAZy.family.xls",'action'=>"文件类型");
				print $output "\n";$writer->emptyTag('table','name'=>"表 ".++$table." $s.CAZy.family.xls表头说明",'type'=>"type1|full",'desc'=>"",'path'=>"pic/header.sample.CAZy.family.xls");
			}
			if(-d "$idir/5.Specialized_function/VFDB/"){
				print $output "\n";$writer->emptyTag('h2', 'name'=>$title1Num.".".++$title2Num." 细菌致病菌毒力因子（VFDB）注释",'type'=>"二级标题显示样式",'desc'=>"二级标题描述");
				print $output "\n";$writer->emptyTag('file','name'=>"结果目录：5.Specialized_function/VFDB",'type'=>"文件显示样式",'desc'=>"",'path'=>"5.Specialized_function/VFDB",'action'=>"文件类型");
				print $output "\n";$writer->emptyTag('p','type'=>'正文段落显示样式','desc'=>"VFDB 数据库全称为Virulence Factors of Pathogenic Bacteria<sup>[<a href='#ref18'>18</a>]</sup>，用于专门研究致病细菌、衣原体和支原体致病因子的数据库。利用diamond<sup>[<a href='#ref10'>10</a>]</sup>与数据库比对，取e<1e-10的注释。数据库链接：http://www.mgc.ac.cn/VFs/main.htm");
				print $output "\n";$writer->emptyTag('file','name'=>"VFDB库注释结果：$s.VFDB.anno.xls",'type'=>"文件显示样式",'desc'=>"",'path'=>"5.Specialized_function/VFDB/$s.VFDB.anno.xls",'action'=>"文件类型");
				print $output "\n";$writer->emptyTag('table','name'=>"表 ".++$table." $s.VFDB.anno.xls表头说明",'type'=>"type1|full",'desc'=>"注：某些致病因子的详细信息，可通过注释结果中的VF_Name这一列，在结果目录下VFs_reference.xls文件中对应的第一列进行对应检索，检索更为详细的信息。",'path'=>"pic/header.sample.VFDB.anno.xls");
			}
			if(-d "$idir/5.Specialized_function/PHI/"){
				print $output "\n";$writer->emptyTag('h2', 'name'=>$title1Num.".".++$title2Num." 病原与宿主互作数据库（PHI）注释",'type'=>"二级标题显示样式",'desc'=>"二级标题描述");
				print $output "\n";$writer->emptyTag('file','name'=>"结果目录：5.Specialized_function/PHI",'type'=>"文件显示样式",'desc'=>"",'path'=>"5.Specialized_function/PHI",'action'=>"文件类型");
				print $output "\n";$writer->emptyTag('p','type'=>'正文段落显示样式','desc'=>"PHI全称为Pathogen Host Interactions Database<sup>[<a href='#ref19'>19</a>]</sup>，病原与宿主互作数据库，该数据库经过了实验验证。PHI是目前医疗、农业真菌和卵菌候选靶向位点的重要在线资源。利用diamond<sup>[<a href='#ref10'>10</a>]</sup>与数据库比对，取e<1e-10的注释。数据库链接：http://www.phi-base.org/");
				print $output "\n";$writer->emptyTag('file','name'=>"PHI库注释结果：$s.PHI.anno.xls",'type'=>"文件显示样式",'desc'=>"",'path'=>"5.Specialized_function/PHI/$s.PHI.anno.xls",'action'=>"文件类型");
				print $output "\n";$writer->emptyTag('table','name'=>"表 ".++$table." $s.PHI.anno.xls表头说明",'type'=>"type1|full",'desc'=>"",'path'=>"pic/header.sample.PHI.anno.xls");
			}
		}
		my @dir=glob "$idir/*\.SignalP/";
		if(@dir!=0 && -d $dir[0]){
			my $nametmp=$1 if($dir[0]=~/\/(\d\.SignalP)/);
			print $output "\n";$writer->emptyTag('h1', 'name'=>++$title1Num." 分泌蛋白预测",'type'=>"一级标题显示样式",'desc'=>"一级标题描述");
			print $output "\n";$writer->emptyTag('file','name'=>"结果目录：$nametmp",'type'=>"文件显示样式",'desc'=>"",'path'=>"$nametmp",'action'=>"文件类型");
			print $output "\n";$writer->emptyTag('p','type'=>'正文段落显示样式','desc'=>"分泌蛋白是指在细胞内合成后，分泌到细胞外起作用的蛋白质。分泌蛋白的N端有一般由5～30个氨基酸组成的信号肽。使用信号肽预测工具SignalP(v4.1)<sup>[<a href='#ref20'>20</a>]</sup>注释蛋白序列是否是分泌蛋白。");
			print $output "\n";$writer->emptyTag('file','name'=>"分泌蛋白预测结果：$s.signalp.xls",'type'=>"文件显示样式",'desc'=>"",'path'=>"$nametmp/$s.signalp.xls",'action'=>"文件类型");
			print $output "\n";$writer->emptyTag('table','name'=>"表 ".++$table." $s.signalp.xls表头说明",'type'=>"type1|full",'desc'=>"",'path'=>"pic/header.sample.signalp.xls");
		}
		@dir=glob "$idir/*\.EffectiveT3/";
		if(@dir!=0 && -d $dir[0]){
			my $nametmp=$1 if($dir[0]=~/\/(\d\.EffectiveT3)/);
			print $output "\n";$writer->emptyTag('h1', 'name'=>++$title1Num." III型分泌系统效应蛋白",'type'=>"一级标题显示样式",'desc'=>"一级标题描述");
			print $output "\n";$writer->emptyTag('file','name'=>"结果目录：$nametmp",'type'=>"文件显示样式",'desc'=>"",'path'=>"$nametmp",'action'=>"文件类型");
			print $output "\n";$writer->emptyTag('p','type'=>'正文段落显示样式','desc'=>"目前在革兰氏阴性菌中已经报道了6种分泌系统，III型分泌系统（T3SS）在大肠杆菌、沙门氏菌、耶尔森氏菌和霍乱弧菌等重要病原菌中均有发现，和细菌的致病性密切相关。其主要特点为效应蛋白可以直接通过分泌系统从细菌胞质进入宿主细胞，发挥毒力作用。因此，对效应蛋白的研究有重要的应用价值，效应蛋白的信号肽可以用于携带蛋白疫苗、携带药物等多种小分子进入细胞，进行细胞水平的免疫或治疗。");
			print $output "\n";$writer->emptyTag('p','type'=>'正文段落显示样式','desc'=>"根据预测基因的氨基酸序列文件，使用软件EffectiveT3(v1.0)<sup>[<a href='#ref21'>21</a>]</sup>（http://effectivedb.org/）来进行III型分泌系统效应蛋白的预测。");
			print $output "\n";$writer->emptyTag('file','name'=>" III型分泌系统效应蛋白预测结果：$s.T3SS.xls",'type'=>"文件显示样式",'desc'=>"",'path'=>"$nametmp/$s.T3SS.xls",'action'=>"文件类型");
			print $output "\n";$writer->emptyTag('table','name'=>"表 ".++$table." $s.T3SS.xls表头说明",'type'=>"type1|full",'desc'=>"",'path'=>"pic/header.sample.T3SS.xls");
		}
		@dir=glob "$idir/*\.TMHMM/";
		if(@dir!=0 && -d $dir[0]){
			my $nametmp=$1 if($dir[0]=~/\/(\d\.TMHMM)/);
			print $output "\n";$writer->emptyTag('h1', 'name'=>++$title1Num." 蛋白跨膜结构预测",'type'=>"一级标题显示样式",'desc'=>"一级标题描述");
			print $output "\n";$writer->emptyTag('file','name'=>"结果目录：$nametmp",'type'=>"文件显示样式",'desc'=>"",'path'=>"$nametmp",'action'=>"文件类型");
			#print $output "\n";$writer->emptyTag('p','type'=>'正文段落显示样式','desc'=>"由双层脂类组成的生物膜结构是细胞的基本构成单元，然而生物膜的大部分功能是由镶嵌在生物膜中的蛋白质来完成的，膜蛋白决定了生物膜的功能特性，因此，不同类型的生物膜，其蛋白构成比例有很大的差别，以质量来算，从25%到75%不等。另一个角度来说，生物体内的蛋白质中，20-30%都属于膜蛋白，因此膜蛋白在细胞的功能中占据重要地位是显而易见的。");
			print $output "\n";$writer->emptyTag('p','type'=>'正文段落显示样式','desc'=>"蛋白质跨越细胞膜的区域通常为α-螺旋结构，由约20~25个氨基酸残基组成，且大部分是疏水性氨基酸。根据预测出来的蛋白序列，使用TMHMM(v2.0)<sup>[<a href='#ref22'>22</a>]</sup>（http://www.cbs.dtu.dk/services/TMHMM/）软件来预测蛋白序列中的跨膜区域结构，TMHMM(v2.0)通过隐马尔科夫模型来预测蛋白质的跨膜螺旋   ，预测准确度为97-98%。虽然当存在信号肽时预测准确度会下降，但软件对可溶性蛋白和疏水跨膜蛋白预测的特异性和敏感性依然高于99%。");
			print $output "\n";$writer->emptyTag('file','name'=>" 蛋白跨膜结构预测结果：$s.tmhmm.xls",'type'=>"文件显示样式",'desc'=>"",'path'=>"$nametmp/$s.tmhmm.xls",'action'=>"文件类型");
			print $output "\n";$writer->emptyTag('table','name'=>"表 ".++$table." $s.tmhmm.xls表头说明",'type'=>"type1|full",'desc'=>"",'path'=>"pic/header.sample.tmhmm.xls");
			print $output "\n";$writer->emptyTag('p','type'=>'正文段落显示样式','desc'=>"同时下图展示了α螺旋在细胞膜inside/outside的后验概率，通过作图可以判断预测的跨膜α螺旋的确定性，较弱的跨膜α螺旋没有被预测。在图形的上方(1~1.2之间)展示了N-best算法的预测结果。");
			print $output "\n";$writer->emptyTag('file','name'=>" TMHMM后验概率图目录：$s\_TMHMM_plot",'type'=>"文件显示样式",'desc'=>"",'path'=>"$nametmp/$s\_TMHMM_plot",'action'=>"文件类型");
			print $output "\n";$writer->emptyTag('pic','name'=>"图 ".++$pic." TMHMM后验概率示意图",'desc'=>"注：该图是通过计算氨基酸残基在细胞膜内/外的概率得到的。有时图形和N-best的预测结果相矛盾，这是因为图形展示的是每个残基在螺旋内外的概率，而N-best是基于对结构的预测。因此图形应该作为一种互补的信息展示方式。",'path'=>"pic/tmhmm_plot.png");
		}
		#文件格式说明
		#print $output "\n";$writer->emptyTag('h1', 'name'=>++$title1Num." 文件格式说明",'type'=>"一级标题显示样式",'desc'=>"一级标题描述");
		#print $output "\n";$writer->emptyTag('table','name'=>"表 ".++$table." gff文件说明",'type'=>"type1|full",'desc'=>"",'path'=>"pic/header.sample.gff");
		#print $output "\n";$writer->emptyTag('table','name'=>"表 ".++$table." 注释总表表头说明",'type'=>"type1|full",'desc'=>"",'path'=>"pic/annotationHeader.txt");
		#print $output "\n";$writer->emptyTag('table','name'=>"表 ".++$table." NR/Swissprot anno文件说明",'type'=>"type1|full",'desc'=>"",'path'=>"pic/nr_swissprot_header.txt");
		
		#参考文献
		print $output "\n";$writer->startTag('ref_list','name'=>"参考文献",'type'=>"参考文献列表显示样式",'desc'=>"参考文献列表描述");
		print $output "\n";$writer->emptyTag('ref','id'=>"1",'name'=>"Bolger A M, Lohse M, Usadel B. Trimmomatic: a flexible trimmer for Illumina sequence data[J]. Bioinformatics, 2014, 30(15): 2114-2120.",'link'=>"https://academic.oup.com/bioinformatics/article/30/15/2114/2390096");
		print $output "\n";$writer->emptyTag('ref','id'=>"2",'name'=>"Liu B, Shi Y, Yuan J, et al. Estimation of genomic characteristics by analyzing k-mer frequency in de novo genome projects[J]. arXiv preprint arXiv:1308.2012, 2013.",'link'=>"https://arxiv.org/abs/1308.2012");
		print $output "\n";$writer->emptyTag('ref','id'=>"3",'name'=>"Hyatt D, Chen G L, LoCascio P F, <em>et al</em>. Prodigal: prokaryotic gene recognition and translation initiation site identification[J]. BMC bioinformatics, 2010, 11(1): 119.",'link'=>"https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-119");
		print $output "\n";$writer->emptyTag('ref','id'=>"4",'name'=>"Lowe T M, Eddy S R. tRNAscan-SE: a program for improved detection of transfer RNA genes in genomic sequence[J]. Nucleic acids research, 1997, 25(5): 955-964.",'link'=>"http://nar.oxfordjournals.org/content/25/5/955.abstract");
		print $output "\n";$writer->emptyTag('ref','id'=>"5",'name'=>"Lagesen K, Hallin P, Rødland E A, <em>et al</em>. RNAmmer: consistent and rapid annotation of ribosomal RNA genes[J]. Nucleic acids research, 2007, 35(9): 3100-3108.",'link'=>"https://academic.oup.com/nar/article-abstract/35/9/3100/2401119");
		print $output "\n";$writer->emptyTag('ref','id'=>"6",'name'=>"Griffiths-Jones S, Bateman A, Marshall M, <em>et al</em>. Rfam: an RNA family database[J]. Nucleic acids research, 2003, 31(1): 439-441.",'link'=>"https://academic.oup.com/nar/article-abstract/31/1/439/2401159");
		print $output "\n";$writer->emptyTag('ref','id'=>"7",'name'=>"Tarailo‐Graovac M, Chen N. Using RepeatMasker to identify repetitive elements in genomic sequences[J]. Current protocols in bioinformatics, 2009: 4.10. 1-4.10. 14.",'link'=>"http://onlinelibrary.wiley.com/doi/10.1002/0471250953.bi0410s25/full");
		print $output "\n";$writer->emptyTag('ref','id'=>"8",'name'=>"Edgar R C. PILER-CR: fast and accurate identification of CRISPR repeats[J]. BMC bioinformatics, 2007, 8(1): 18.",'link'=>"https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-8-18");
		print $output "\n";$writer->emptyTag('ref','id'=>"9",'name'=>"Bland C, Ramsey T L, Sabree F, <em>et al</em>. CRISPR recognition tool (CRT): a tool for automatic detection of clustered regularly interspaced palindromic repeats[J]. BMC bioinformatics, 2007, 8(1): 209.",'link'=>"https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-8-209");
		print $output "\n";$writer->emptyTag('ref','id'=>"10",'name'=>"Akhter S, Aziz R K, Edwards R A. PhiSpy: a novel algorithm for finding prophages in bacterial genomes that combines similarity-and composition-based strategies[J]. Nucleic acids research, 2012, 40(16): e126-e126.",'link'=>"https://academic.oup.com/nar/article-abstract/40/16/e126/1027055");
		print $output "\n";$writer->emptyTag('ref','id'=>"11",'name'=>"Benjamin Buchfink, Chao Xie, Daniel H Huson. Fast and sensitive protein alignment using DIAMOND[J]. Nature Methods, 2015, 12: 59-60.",'link'=>"https://www.nature.com/articles/nmeth.3176");
		print $output "\n";$writer->emptyTag('ref','id'=>"12",'name'=>"SEAN R. EDDY. A new generation of homology search tools based on probabilistic inference[J]. Genome Informatics, 2009, 205-211.",'link'=>"http://www.worldscientific.com/doi/abs/10.1142/9781848165632_0019");
		print $output "\n";$writer->emptyTag('ref','id'=>"13",'name'=>"Xie C, Mao X, Huang J, <em>et al</em>. KOBAS 2.0: a web server for annotation and identification of enriched pathways and diseases[J]. Nucleic Acids Research, 2011, 39, W316–W322.",'link'=>"https://academic.oup.com/nar/article-lookup/doi/10.1093/nar/gkr483");
		print $output "\n";$writer->emptyTag('ref','id'=>"14",'name'=>"Jaina Mistry, Robert D. Finn, Sean R. Eddy, <em>et al</em>. Challenges in homology search: HMMER3 and convergent evolution of coiled-coil regions[J]. Nucleic Acids Research, 2013, 41(12): e121.",'link'=>"https://academic.oup.com/nar/article/41/12/e121/1025950");
		print $output "\n";$writer->emptyTag('ref','id'=>"15",'name'=>"Jia B, Raphenya A R, Alcock B, <em>et al</em>. CARD 2017: expansion and model-centric curation of the comprehensive antibiotic resistance database[J]. Nucleic acids research, 2017, 45(D1): D566-D573.",'link'=>"https://academic.oup.com/nar/article-abstract/45/D1/D566/2333912");
		print $output "\n";$writer->emptyTag('ref','id'=>"16",'name'=>"Altschul, S.F., Gish, W., Miller, <em>et al</em>. Basic local alignment search tool[J]. Journal of Molecular Biology, 1990, 215(3): 403-410.",'link'=>"http://linkinghub.elsevier.com/retrieve/pii/S0022-2836(05)80360-2");
		print $output "\n";$writer->emptyTag('ref','id'=>"17",'name'=>"Cantarel B L, Coutinho P M, Rancurel C, <em>et al</em>. The Carbohydrate-Active EnZymes database (CAZy): an expert resource for glycogenomics[J]. Nucleic acids research, 2008, 37(suppl_1): D233-D238.",'link'=>"https://academic.oup.com/nar/article/37/suppl_1/D233/1003505");
		print $output "\n";$writer->emptyTag('ref','id'=>"18",'name'=>"Chen L, Yang J, Yu J, <em>et al</em>. VFDB: a reference database for bacterial virulence factors[J]. Nucleic acids research, 2005, 33(suppl_1): D325-D328.",'link'=>"https://academic.oup.com/nar/article-abstract/33/suppl_1/D325/2505203");
		print $output "\n";$writer->emptyTag('ref','id'=>"19",'name'=>"Winnenburg R, Baldwin T K, Urban M, <em>et al</em>. PHI-base: a new database for pathogen host interactions[J]. Nucleic acids research, 2006, 34(suppl_1): D459-D464.",'link'=>"https://academic.oup.com/nar/article-abstract/34/suppl_1/D459/1132790");
		print $output "\n";$writer->emptyTag('ref','id'=>"20",'name'=>"Petersen T N, Brunak S, von Heijne G, <em>et al</em>. SignalP 4.0: discriminating signal peptides from transmembrane regions[J]. Nature methods, 2011, 8(10): 785-786.",'link'=>"https://www.nature.com/articles/nmeth.1701");
		print $output "\n";$writer->emptyTag('ref','id'=>"21",'name'=>"Jehl M A, Arnold R, Rattei T. Effective—a database of predicted secreted bacterial proteins[J]. Nucleic acids research, 2010, 39(suppl_1): D591-D595.",'link'=>"https://academic.oup.com/nar/article-abstract/39/suppl_1/D591/2507337");
		print $output "\n";$writer->emptyTag('ref','id'=>"22",'name'=>"Sonnhammer E L L, Von Heijne G, Krogh A. A hidden Markov model for predicting transmembrane helices in protein sequences. Proceedings of the 6th International Conference on Intelligent Systems for Molecular Biology[C]. 1998: 175-182.",'link'=>"http://www.aaai.org/Papers/ISMB/1998/ISMB98-021.pdf");
		print $output "\n";$writer->endTag("ref_list");
		
		print $output "\n";$writer->endTag("report");
		print $output "\n";$writer->end();
		print $output "\n";$output->close();
		system("cd $odir/ && $python $xml2html -i $s.xml -n 样品$s\项目报告 && cd -");
#		system("rm -rf $odir/pic $odir/$s.xml");
	}
}
