#!/usr/bin/perl
use strict;
use Getopt::Long;
Getopt::Long::Configure qw(no_ignore_case);
use File::Basename;
use FindBin qw($Bin);
use lib "$Bin/.Modules";
use Parameter::BinList;
use Sort::ChrPos;
use BamRelated::BaseAndQual;
use VcfRelated::VcfParse;
use SeqRelated::Seq;
use warnings;

my ($HelpFlag,$BinList,$BeginTime);
my ($BwaBin,$SamtoolsBin,$tmpDir,$ThreadNum,$ReadLen,$StrictFlag);
my ($SortedBam,$Ref,$QsubNode);
my @Fq;
my $HelpInfo = <<USAGE;

 Fq2SortedBam_v2.pl
 Auther: zhangdong_xie\@foxmail.com

  This script was used to map fq.gz files to ref.

 -i      ( Required ) Fq.gz files;
 -b      ( Required ) Bam logging file;
 -ref    ( Optional ) Genome reference which should be indexed (default: ucsc.hg19.fasta);
 -qsub   ( Optional ) Qsub node if qsub needed (default: off, like -qsub c0035);
 -t      ( Optional ) Thread number (default: 5);
 -d      ( Optional ) Temporary directory;
 -r      ( Optional ) Read length;
 -s      ( Optional ) Strict mapping when its a snv for short length mapping (with -s only);

 -bin    List for searching of related bin or scripts; 
 -h      Help infomation;

USAGE

GetOptions(
	'i=s' => \@Fq,
	'b=s' => \$SortedBam,
	'ref:s' => \$Ref,
	'qsub:s' => \$QsubNode,
	't:i' => \$ThreadNum,
	'd:s' => \$tmpDir,
	'r:i' => \$ReadLen,
	's!' => \$StrictFlag,
	'bin:s' => \$BinList,
	'h!' => \$HelpFlag
) or die $HelpInfo;

if($HelpFlag || !@Fq || !$SortedBam)
{
	die $HelpInfo;
}
else
{
	$BeginTime = ScriptBegin();
	
	for my $i (0 .. $#Fq)
	{
		die "[ Error ] Fq not exist ($Fq[$i]).\n" unless(-s $Fq[$i]);
	}
	my $Dir = dirname $SortedBam;
	die "[ Error ] Directory not exist ($Dir).\n" unless(-d $Dir);
	$ThreadNum = 5 if(!$ThreadNum);
	$BinList = BinListGet() unless($BinList);
	$Ref = BinSearch("hg19",$BinList) unless($Ref);
	$tmpDir = BinSearch("tmpOnly",$BinList) unless($tmpDir);
	$BwaBin = BinSearch("bwa",$BinList);
	$SamtoolsBin = BinSearch("samtools",$BinList);
}

if(@Fq && $SortedBam)
{
	my $tmpId = RandString(10);
	my $tmpShell = $tmpDir . "/tmpFq2Bam" . $tmpId . ".sh";
	my $tmpLog = $tmpDir . "/tmpFq2Bam" . $tmpId . ".log";
	
	my $tFq = join(" ",@Fq);
	my $Line;
	$Line = "$BwaBin mem -t $ThreadNum -M -R \'\@RG\\tID:1\\tLB:lib1\\tPL:illumina\\tPU:Unit1\\tSM:ASHD753\\tCN:Annoroad\' $Ref $tFq | ";
	$Line .= "$SamtoolsBin view -bh | $SamtoolsBin sort -o $SortedBam";
	
	# 单端且reads读长短于70bp;
	if($#Fq == 0 && $ReadLen > 0 && $ReadLen < 70 && !$StrictFlag)
	{
		my $Prefix = $SortedBam;
		$Prefix =~ s/bam$//;
		my $Sai = $Prefix . "sai";
		my $Sam = $Prefix . "sam";
		
		$Line = "$BwaBin aln -n 5 -t 5 -o 1 -m 2000000000 -l 20 -i 20 -d 20 -q 10  $Ref $Fq[0] > $Sai\n";
		$Line .= "$BwaBin samse $Ref $Sai $Fq[0] -f $Sam -r \'\@RG\\tID:1\\tLB:lib1\\tPL:illumina\\tPU:Unit1\\tSM:ASHD753\\tCN:Annoroad\'\n";
		$Line .= "$SamtoolsBin view -bh $Sam | $SamtoolsBin sort -o $SortedBam\n";
		$Line .= "rm $Sai $Sam";
	}
	
	open(SH,"> $tmpShell") or die $!;
	print SH $Line;
	close SH;
	`chmod 770 $tmpShell`;
	
	if($QsubNode)
	{
		`qsub -wd $tmpDir -q cancer.q\@$QsubNode.local -l vf=4G $tmpShell`;
	}
	else
	{
		`sh $tmpShell > $tmpLog`;
		`chmod 770 $tmpLog`;
	}
	
	`rm $tmpShell` if($tmpShell);
	`rm $tmpLog` if($tmpLog);
}
printf "[ %.2fmin ] The end.\n",(time - $BeginTime)/60;


######### Sub functions ##########
