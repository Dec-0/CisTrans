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
my ($Bam,$LogBam,$ThreadNum,$SclipFlag,$ExtendFlag,$RMFlag,$UnMapSaveFlag,$StrictFlag);
my ($Reference,$SamtoolsBin,$BwaBin,$Fq2BamScript,$BedtoolsBin,$ExtendNum,$MaxDiffNum);
my (@Var,@Fq);
my $HelpInfo = <<USAGE;

 IndelConfirm_v4.pl
 Auther: zhangdong_xie\@foxmail.com

 This script was used to confirm specified indel(s) from fq.gz or bam one at a time.
 和v3版本相比，本次修改了输入变异的格式，不再要求前期的人工转换。

 -v       ( Required ) Variant info, only support like 'chr1,12345,12345,A,T' or 'chr1,12345,12345,-,ATG' or 'chr1,12345,12347,ATG,-' [ multi times ];
 -i       ( Required ) fq.gz files (can be specified 1 or 2 times);
 -b       ( Required ) Bam file (.bai indexed);
                       '-i' and '-b' can be specified only one at a time;
 -o       ( Required ) Result logging bam file (in bam format);
 
 -t       ( Optional ) The threads\' number which bwa will adapt (default: 5);
 -sclip   ( Optional ) If the soft clip will be checked (with -sclip only);
 -extend  ( Optional ) For the checking of flanking sequences which can match partially;
          Only support two numbers seperated by ',', like '4,2' which means extending 4bp at each side and only tolerates
          2bp mismatch in total;
 -move    ( Optional ) If deleting the raw re-mapping bam;
 -save    ( Optional ) If saving the unmapped reads after mapping to the new ref;
 -s       ( Optional ) Strict mapping when its a snv for short length mapping (with -s only);

 -bin     List for searching of related bin or scripts; 
 -h       Help infomation;

USAGE

GetOptions(
	'v=s' => \@Var,
	'i=s' => \@Fq,
	'b=s' => \$Bam,
	'o=s' => \$LogBam,
	't:i' => \$ThreadNum,
	'sclip!' => \$SclipFlag,
	'extend:s' => \$ExtendFlag,
	'move!' => \$RMFlag,
	'save!' => \$UnMapSaveFlag,
	's!' => \$StrictFlag,
	'bin:s' => \$BinList,
	'h!' => \$HelpFlag
) or die $HelpInfo;

if($HelpFlag || !@Var || (!@Fq && !$Bam) || (@Fq && $Bam) || !$LogBam)
{
	die $HelpInfo;
}
else
{
	$BeginTime = ScriptBegin();
	
	if(@Fq)
	{
		die "[ Error ] Too much fq.gz.\n" if($#Fq > 1);
		for my $i (0 .. $#Fq)
		{
			die "[ Error ] Fq not exist: $Fq[$i].\n" unless(-e $Fq[$i]);
		}
	}
	elsif($Bam)
	{
		die "[ Error ] Bam not exist: $Bam.\n" unless(-e $Bam);
		my $BamIndex = $Bam . ".bai";
		die "[ Error ] Index of bam not exist: $BamIndex.\n" unless(-e $BamIndex);
	}
	# 默认值;
	$SclipFlag = 0 unless($SclipFlag);
	$ThreadNum = 5 if(!$ThreadNum);
	$StrictFlag = 0 unless($StrictFlag);
	$BinList = BinListGet() if(!$BinList);
	$Reference = BinSearch("hg19",$BinList);
	$SamtoolsBin = BinSearch("samtools",$BinList);
	$BwaBin = BinSearch("bwa",$BinList);
	$Fq2BamScript = BinSearch("Fq2SortedBam",$BinList);
	$BedtoolsBin = BinSearch("BedTools-2.26",$BinList);
}

if(@Var)
{
	# initial;
	my $LogDir = dirname $LogBam;
	`mkdir -p $LogDir` unless(-d $LogDir || !$LogDir);
	
	# 输入检查以及变异规整;
	my (@Chr,@LBoundary,@RBoundary,@Ref,@Alt);
	for my $i (0 .. $#Var)
	{
		my ($tChr,$tLBoundary,$tRBoundary,$tRef,$tAlt) = &VarCheck($Var[$i]);
		push @Chr, $tChr;
		push @LBoundary, $tLBoundary;
		push @RBoundary, $tRBoundary;
		push @Ref, $tRef;
		push @Alt, $tAlt;
	}
	# 变异排序;
	if(1)
	{
		my @OtherInfo = ();
		for my $i (0 .. $#Chr)
		{
			$OtherInfo[$i] = $Ref[$i] . "\t" . $Alt[$i];
		}
		
		my @tRef = ChrPosAndOther(\@Chr,\@LBoundary,\@RBoundary,\@OtherInfo);
		@Chr = @{$tRef[0]};
		@LBoundary = @{$tRef[1]};
		@RBoundary = @{$tRef[2]};
		@OtherInfo = @{$tRef[3]};
		for my $i (0 .. $#OtherInfo)
		{
			my @Cols = split /\t/, $OtherInfo[$i];
			$Cols[1] = "" unless($Cols[1]);
			
			$Ref[$i] = $Cols[0];
			$Alt[$i] = $Cols[1];
			$Var[$i] = $Chr[$i] . "," . $LBoundary[$i] . "," . $Ref[$i] . "," . $Alt[$i];
		}
	}
	# 检查各种参数、条件是否合理;
	my $ReadLen = 0;
	$ReadLen = ReadLenConfirmFromFq($Fq[0]) if(@Fq);
	$ReadLen = ReadLenConfirmFromBam($Bam) if($Bam && !$ReadLen);
	printf "[ %.2fmin Info ] Read length: %d.\n",(time - $BeginTime)/60,$ReadLen;
	printf "[ %.2fmin Info ] Variant: %s,%s,%s,%s,%s.\n",(time - $BeginTime)/60,$Chr[0],$LBoundary[0],$RBoundary[0],$Ref[0],$Alt[0];
	for my $i (1 .. $#Chr)
	{
		die "[ Error ] Chr not consistent ($Chr[$i] vs $Chr[0]).\n" if($Chr[$i] ne $Chr[0]);
		die "[ Error ] Var pos crossing ($Chr[$i],$LBoundary[$i],$RBoundary[$i] vs $Chr[$i - 1],$LBoundary[$i - 1],$RBoundary[$i - 1]).\n" if($LBoundary[$i] + 1 < $RBoundary[$i - 1]);
		die "[ Error ] Out the limit of one read ($Chr[$i],$LBoundary[$i] vs $Chr[0],$LBoundary[0]).\n" if(abs($LBoundary[$i] - $LBoundary[0]) >= $ReadLen);
		printf "[ %.2fmin Info ] Variant: %s,%s,%s,%s,%s.\n",(time - $BeginTime)/60,$Chr[$i],$LBoundary[$i],$RBoundary[$i],$Ref[$i],$Alt[$i];
	}
	
	
	# 构建新的参考基因序列;
	my $MinPos = $LBoundary[0] - $ReadLen + 1;
	my $MaxPos = $RBoundary[-1] + $ReadLen - 1;
	printf "[ %.2fmin New reference ]\t",(time - $BeginTime)/60;
	my $NewRef = RefGet($Reference,$Chr[0],$MinPos,$LBoundary[0],$BedtoolsBin) . $Alt[0];
	print RefGet($Reference,$Chr[0],$MinPos,$LBoundary[0],$BedtoolsBin),"-(",RefGetWithOutWarning($Reference,$Chr[0],$LBoundary[0] + 1,$RBoundary[0] - 1,$BedtoolsBin),")",$Alt[0];
	# Boundary的坐标是1-based的;
	$LBoundary[0] = $ReadLen;
	for my $i (1 .. $#Chr)
	{
		my $tLen = length($NewRef);
		$NewRef .= RefGetWithOutWarning($Reference,$Chr[0],$RBoundary[$i - 1],$LBoundary[$i],$BedtoolsBin) . $Alt[$i];
		print "-",RefGetWithOutWarning($Reference,$Chr[0],$RBoundary[$i - 1],$LBoundary[$i],$BedtoolsBin),"-(",RefGetWithOutWarning($Reference,$Chr[0],$LBoundary[$i] + 1,$RBoundary[$i] - 1,$BedtoolsBin),")",$Alt[$i];
		$RBoundary[$i - 1] = $tLen + 1;
		$LBoundary[$i] = length($NewRef) - length($Alt[$i]);
	}
	$NewRef .= RefGet($Reference,$Chr[0],$RBoundary[-1],$MaxPos,$BedtoolsBin);
	print "-",RefGet($Reference,$Chr[0],$RBoundary[-1],$MaxPos,$BedtoolsBin),"\n";
	$RBoundary[-1] = length($NewRef) - $ReadLen + 1;
	
	
	# 完全匹配是为了确定的确是发生了这种变异所代表的变化;
	my (@ExactSeq,@ExactFrom,@ExactTo) = ();
	my (@ExcludSeq,@ExcludFrom,@ExcludTo) = ();
	for my $i (0 .. $#Var)
	{
		my $LSeq = substr($NewRef,0,$LBoundary[$i]);
		my $RSeq = substr($NewRef,$RBoundary[$i] - 1);
		# 1-based;
		($ExactSeq[$i],$ExactFrom[$i],$ExactTo[$i]) = &ExactConfirm($Ref[$i],$Alt[$i],$LSeq,$RSeq);
		printf "[ %.2fmin Info ] Seq has to match exactly for %s:\t%s(%d - %d)[1-based].\n",(time - $BeginTime)/60,$Var[$i],$ExactSeq[$i],$ExactFrom[$i],$ExactTo[$i];
		
		# 假如ref有序列，则需要排除这个段序列的存在，以免引起误判;
		if($Ref[$i])
		{
			($ExcludSeq[$i],$ExcludFrom[$i],$ExcludTo[$i]) = &ExcludConfirm($Ref[$i],$LSeq,$RSeq);
			printf "[ %.2fmin Info ] Seq has to excluding for %s:\t%s(%d - %d)[1-based].\n",(time - $BeginTime)/60,$Var[$i],$ExcludSeq[$i],$ExcludFrom[$i],$ExcludTo[$i];
		}
	}
	# 部分匹配是为了过滤周围序列杂乱无章的情况，此时很可能是由于错配导致的;
	my @PartSeq;
	if($ExtendFlag)
	{
		($ExtendNum,$MaxDiffNum) = split /,/, $ExtendFlag;
		die "[ Error ] Empty ExtendFlag.\n" if(!$ExtendNum || !$MaxDiffNum);
		die "[ Error ] ExtendFlag not right.\n" if($ExtendNum =~ /\D/ || $MaxDiffNum =~ /\D/);
		
		for my $i (0 .. $#Var)
		{
			my $tLPartSeq = substr($NewRef,$ExactFrom[$i] - $ExtendNum - 1,$ExtendNum);
			my $tRPartSeq = substr($NewRef,$ExactTo[$i],$ExtendNum);
			$PartSeq[$i] = $tLPartSeq . $tRPartSeq;
		}
	}
	
	
	# new ref seq logging and index;
	my $tId = RandString(10);
	my $NewRefFa = $LogDir . "/Ref4IndelConform." . $tId . ".fa";
	&RefMake($NewRefFa,$Chr[0],$NewRef);
	
	
	# bam to fq and fq to bam;
	if($Bam)
	{
		# bam to fq;
		# reads on specific area and those un-mapped; One fq in case paired not found;
		@Fq = ();
		$Fq[0] = $LogDir . "/Fq4IndelConfirm." . $tId . ".fq.gz";
		&Bam2Fq($Bam,$Fq[0],$Chr[0],$MinPos + $ReadLen - 2,$MaxPos - $ReadLen + 1);
	}
	# fq re-mapping to new ref;
	my $SortedBam = $LogBam;
	$SortedBam =~ s/bam$//;
	$SortedBam .= $tId . ".original.bam";
	&Fq2Bam(\@Fq,$NewRefFa,$SortedBam,$ThreadNum,$LogDir,$ReadLen,$StrictFlag);
	if($Bam)
	{
		`rm $Fq[0]` unless(!$Fq[0] || $Fq[0] =~ /\*/);
	}
	`rm $NewRefFa\*` if($NewRefFa && -e $NewRefFa);
	
	
	# bam filtering for all variants;
	&UnMapFilter($SortedBam);
	# counting;
	my %CombineCount = ();
	my @SoloCount = ();
	my ($FullDp,$AltDp,$FNum,$RNum) = (0,0,0,0);
	`$SamtoolsBin index $SortedBam`;
	my $BamIndex = $SortedBam . ".bai";
	`rm $BamIndex` if($BamIndex && -e $BamIndex);
	open(TBAM,"$SamtoolsBin view -h $SortedBam |") or die $!;
	open(LOG,"| $SamtoolsBin view -bh > $LogBam") or die $!;
	while(my $Line = <TBAM>)
	{
		if($Line =~ /^@/)
		{
			print LOG $Line;
			next;
		}
		last if(&IfUnMap($Line));
		
		my ($RangeChr,$RangeFrom,$RangeTo) = ReadCoverRange($Line,$SclipFlag);
		next unless($RangeFrom <= $ExactFrom[0] && $RangeTo >= $ExactTo[-1]);
		next if($ExcludSeq[0] && $ExcludFrom[0] < $RangeFrom);
		next if($ExcludSeq[-1] && $ExcludTo[-1] > $RangeTo);
		$FullDp ++;
		
		my $MatchFlag = 1;
		my $CombineVar = "";
		# exact matching check;
		for my $i (0 .. $#Var)
		{
			my $tFlag = 1;
			my $RExactSeq = AltGet($Chr[0],$ExactFrom[$i],$ExactTo[$i],$Line,$SclipFlag);
			if($RExactSeq ne $ExactSeq[$i])
			{
				$MatchFlag = 0;
				$tFlag = 0;
			}
			
			if($ExcludSeq[$i])
			{
				my $RExcludSeq = AltGet($Chr[0],$ExcludFrom[$i],$ExcludTo[$i],$Line,$SclipFlag);
				if($RExcludSeq eq $ExcludSeq[$i])
				{
					$MatchFlag = 0;
					$tFlag = 0;
				}
			}
			
			if($ExtendFlag && $tFlag)
			{
				my $RLPartSeq = AltGet($Chr[0],$ExactFrom[$i] - $ExtendNum,$ExactFrom[$i] - 1,$Line,$SclipFlag);
				my $RRPartSeq = AltGet($Chr[0],$ExactTo[$i] + 1,$ExactTo[$i] + $ExtendNum,$Line,$SclipFlag);
				my $RPartSeq = $RLPartSeq . $RRPartSeq;
				
				if(&BaseCompare($PartSeq[$i],$RPartSeq) > $MaxDiffNum)
				{
					$MatchFlag = 0;
					$tFlag = 0;
				}
			}
			
			$CombineVar .= ";" . $Var[$i] if($tFlag);
			$SoloCount[$i] ++ if($tFlag);
		}
		
		if($CombineVar =~ s/^;//)
		{
			if($CombineCount{$CombineVar})
			{
				$CombineCount{$CombineVar} ++;
			}
			else
			{
				$CombineCount{$CombineVar} = 1;
			}
		}
		
		if($MatchFlag)
		{
			my @Cols = split /\t/, $Line;
			if(($Cols[1] / 16) % 2)
			{
				$RNum ++;
			}
			else
			{
				$FNum ++;
			}
			print LOG $Line;
		}
	}
	close TBAM;
	close LOG;
	$AltDp = $FNum + $RNum;
	if($RMFlag && $SortedBam && -e $SortedBam)
	{
		`rm $SortedBam`;
	}
	elsif(-e $SortedBam)
	{
		my $tSortedBam = $SortedBam;
		$tSortedBam =~ s/\.$tId//;
		`rename $SortedBam $tSortedBam $SortedBam`;
		print "[ Info ] Original bam is:\t$tSortedBam.\n";
	}
	
	
	print "[ Info ] Final bam is:\t$LogBam.\n";
	print "# 请注意，本脚本只分析了覆盖所有变异的reads，而且对于紧挨着变异的判断可能会有问题。\n";
	for my $i (0 .. $#Var)
	{
		$SoloCount[$i] = 0 unless($SoloCount[$i]);
		print "[ Info ] Solo count for $Var[$i] is:\t$SoloCount[$i]\n";
	}
	foreach my $Key (sort keys %CombineCount)
	{
		print "[ Info ] Total reads that match these var (",$Key,") are:\t$CombineCount{$Key}\n";
	}
	my $Freq = "-";
	$Freq = sprintf("%.4f",$AltDp / $FullDp) if($FullDp);
	print "[ Info ] Total reads that match all these vars are (Freq,Forward,Reverse,FullDepth,AltDepth):\t$Freq,$FNum,$RNum,$FullDp,$AltDp\n";
}
printf "[ %.2fmin ] The end.\n",(time - $BeginTime)/60;


######### Sub functions ##########
sub VarCheck
{
	my $Var = $_[0];
	
	my ($Chr,$From,$To,$Ref,$Alt) = split /,/, $Var;
	my $tmp = VarTrans($Reference,$Chr,$From,$To,$Ref,$Alt,$BedtoolsBin);
	my ($LB,$RB) = ($From,$To);
	($Chr,$LB,$Ref,$Alt) = split /,/, $tmp;
	$Ref = "" if($Ref eq "-");
	$Alt = "" if($Alt eq "-");
	$RB = $LB + length($Ref) + 1;
	
	return $Chr,$LB,$RB,$Ref,$Alt;
}

sub ExactConfirm
{
	my ($tRef,$tAlt,$tLSeq,$tRSeq) = @_;
	my $tSeq = "";
	my ($tFrom,$tTo) = (0,0);
	
	# 这一步主要是防止多个拷贝中一个的插入或者删除，此时需要囊括所有的拷贝;
	# 0-based;
	if(!$tRef)
	{
		# insertion;
		$tFrom = &SameTruncLocate(0,$tLSeq,$tAlt);
		# 需要左侧的碱基和ref匹配以判断indel;
		$tFrom --;
		$tSeq = substr($tLSeq,$tFrom);
		
		my $tLen = &SameTruncLocate(1,$tRSeq,$tAlt);
		$tTo = length($tLSeq) + length($tAlt) + $tLen;
		$tSeq .= $tAlt . substr($tRSeq,0,$tLen + 1);
	}
	elsif(!$tAlt)
	{
		# deletion;
		$tFrom = &SameTruncLocate(0,$tLSeq,$tRef);
		# 需要左侧的碱基和ref匹配以判断indel;
		$tFrom --;
		$tSeq = substr($tLSeq,$tFrom);
		
		my $tLen = &SameTruncLocate(1,$tRSeq,$tRef);
		$tTo = length($tLSeq) + length($tAlt) + $tLen;
		$tSeq .= $tAlt . substr($tRSeq,0,$tLen + 1);
	}
	elsif(length($tRef) == 1 && length($tAlt) == 1)
	{
		$tSeq = $tAlt;
		$tFrom = length($tLSeq);
		$tTo = length($tLSeq);
	}
	else
	{
		# replacing;
		$tFrom = &SameTruncLocate(0,$tLSeq,$tAlt);
		$tSeq = substr($tLSeq,$tFrom);
		
		my $tLen = &SameTruncLocate(1,$tRSeq,$tAlt);
		$tTo = length($tLSeq) + length($tAlt) + $tLen;
		$tSeq .= $tAlt . substr($tRSeq,0,$tLen + 1);
	}
	
	# when the seq was empty, tTo will be smaller than tFrom;
	# 1-based;
	$tFrom ++;
	$tTo ++;
	
	return $tSeq,$tFrom,$tTo;
}

sub ExcludConfirm
{
	my ($tRef,$tLSeq,$tRSeq) = @_;
	my $tSeq = "";
	my ($tFrom,$tTo) = (0,0);
	
	# 这一步主要是防止多个拷贝中一个的插入或者删除，此时需要囊括所有的拷贝;
	# 0-based;
	$tFrom = &SameTruncLocate(0,$tLSeq,$tRef);
	$tSeq = substr($tLSeq,$tFrom);
	
	my $tLen = &SameTruncLocate(1,$tRSeq,$tRef);
	$tTo = length($tLSeq) + length($tRef) + $tLen;
	$tSeq .= $tRef . substr($tRSeq,0,$tLen + 1);
	
	# when the seq was empty, tTo will be smaller than tFrom;
	# 1-based;
	$tFrom ++;
	$tTo ++;
	
	return $tSeq,$tFrom,$tTo;
}

# 完全匹配，并非部分匹配;
sub SameTruncLocate
{
	my ($Flag,$Full,$M) = @_;
	my $tPos = length($Full);
	
	# 0-based && exact match;
	my $MFlag = 1;
	if($Flag)
	{
		# from left to right;
		while($MFlag)
		{
			$MFlag = 0 unless($Full =~ s/^$M//);
		}
		$tPos = $tPos - length($Full) - 1;
	}
	else
	{
		# from right to left;
		while($MFlag)
		{
			$MFlag = 0 unless($Full =~ s/$M$//);
		}
		$tPos = length($Full);
	}
	
	return $tPos;
}

sub RefMake
{
	my ($Fafile,$Chr,$Seq) = @_;
	
	open(REF,"> $Fafile") or die $!;
	print REF ">$Chr\n$Seq\n";
	close REF;
	
	# index;
	my $Return = `$BwaBin index $Fafile`;
	print "$Return\n" if($Return);
	
	$Fafile =~ s/fa$//;
	$Return = `chmod 660 $Fafile\*`;
	print "$Return\n" if($Return);
	
	return 1;
}

sub Bam2Fq
{
	my ($OriBam,$Fq,$Chr,$From,$To) = @_;
	
	open(FQ,"| gzip > $Fq") or die $!;
	my $ReadFlag = 0;
	my $Return = `$SamtoolsBin view $OriBam | head -n 1 | cut -f 1`;
	chomp $Return;
	my @Items = split /:/, $Return;
	if($#Items == 7 && ($Items[-1] eq "1" || $Items[-1] eq "2"))
	{
		$ReadFlag = 1;
		open(BAM,"$SamtoolsBin view -F 0x500 $OriBam | cut -f 1,2,6,10,11 |") or die $!;
	}
	else
	{
		open(BAM,"$SamtoolsBin view -F 0x500 $OriBam $Chr:$From\-$To | cut -f 1,2,6,10,11 |") or die $!;
	}
	while(my $Line = <BAM>)
	{
		chomp $Line;
		my @Cols = split /\t/, $Line;
		next if($Cols[2] =~ /H/);
		
		# read1 or read2;
		unless($ReadFlag)
		{
			if(($Cols[1] / 64) % 2)
			{
				$Cols[0] .= ":1";
			}
			elsif(($Cols[1] / 128) % 2)
			{
				$Cols[0] .= ":2";
			}
			else
			{
				# SE时前两个条件都不能满足;
				$Cols[0] .= ":0";;
			}
		}
		$Cols[0] = "@" . $Cols[0];
		
		# forward or reverse;
		if(($Cols[1] / 16) % 2)
		{
			$Cols[3] =~ tr/atcgATCG/tagcTAGC/;
			$Cols[3] = reverse $Cols[3];
			
			$Cols[4] = reverse $Cols[4];
		}
		
		print FQ join("\n",$Cols[0],$Cols[3],"+",$Cols[4]),"\n";
	}
	close BAM;
	open(BAM,"$SamtoolsBin view -f 0x4 $OriBam| cut -f 1,2,6,10,11 |") or die $!;
	while(my $Line = <BAM>)
	{
		chomp $Line;
		my @Cols = split /\t/, $Line;
		
		# read1 or read2;
		unless($ReadFlag)
		{
			if(($Cols[1] / 64) % 2)
			{
				$Cols[0] .= ":1";
			}
			elsif(($Cols[1] / 128) % 2)
			{
				$Cols[0] .= ":2";
			}
			else
			{
				# SE时前两个条件都不能满足;
				$Cols[0] .= ":0";;
			}
		}
		$Cols[0] = "@" . $Cols[0];
		
		# forward or reverse;
		if(($Cols[1] / 16) % 2)
		{
			$Cols[3] =~ tr/atcgATCG/tagcTAGC/;
			$Cols[3] = reverse $Cols[3];
			
			$Cols[4] = reverse $Cols[4];
		}
		
		print FQ join("\n",$Cols[0],$Cols[3],"+",$Cols[4]),"\n";
	}
	close BAM;
	close FQ;
	
	return 1;
}

sub Fq2Bam
{
	my @Fq = @{$_[0]};
	my $Fa = $_[1];
	my $SortedBam = $_[2];
	my $TNum = $_[3];
	my $LogDir = $_[4];
	my $Len = $_[5];
	my $StrictFlag = $_[6];
	
	my $tFq = join(" -i ",@Fq);
	`perl $Fq2BamScript -i $tFq -b $SortedBam -ref $Fa -t $TNum -d $LogDir -r $Len` unless($StrictFlag);
	`perl $Fq2BamScript -i $tFq -b $SortedBam -ref $Fa -t $TNum -d $LogDir -r $Len -s` if($StrictFlag);
	
	return 1;
}

sub UnMapFilter
{
	my $OriBam = $_[0];
	
	my $tBam = $OriBam;
	$tBam =~ s/bam//;
	$tBam .= "tmp.bam";
	open(TBAM,"$SamtoolsBin view -h $OriBam |") or die $!;
	open(LOG,"| $SamtoolsBin view -bh > $tBam") or die $!;
	while(my $Line = <TBAM>)
	{
		if($Line =~ /^@/)
		{
			print LOG $Line;
			next;
		}
		if(&IfUnMap($Line) && !$UnMapSaveFlag)
		{
			last;
		}
		my @Cols = split /\t/, $Line;
		next if($Cols[4] == 0 && $Cols[5] eq "*" && !$UnMapSaveFlag);
		
		print LOG $Line;
	}
	close TBAM;
	close LOG;
	`rm $OriBam` unless(!$OriBam || $OriBam =~ /\*/);;
	`rename $tBam $OriBam $tBam`;
	
	return 1;
}

sub IfUnMap
{
	my $tLine = $_[0];
	
	chomp $tLine;
	my @Cols = split /\t/, $tLine;
	my $Flag = 0;
	if($Cols[2] eq "*" && $Cols[3] == 0)
	{
		$Flag = 1;
	}
	
	return $Flag;
}

sub BaseCompare
{
	my ($SeqA,$SeqB) = @_;
	my @ColA = split //, $SeqA;
	my @ColB = split //, $SeqB;
	my $Diff = 0;
	
	my ($LeftDiff,$RightDiff) = (0,0);
	if($#ColA >= $#ColB)
	{
		$LeftDiff = $#ColA - $#ColB;
		for my $i (0 .. $#ColB)
		{
			if($ColA[$i] ne $ColB[$i])
			{
				$LeftDiff ++;
			}
		}
	}
	else
	{
		$LeftDiff = $#ColB - $#ColA;
		for my $i (0 .. $#ColA)
		{
			if($ColA[$i] ne $ColB[$i])
			{
				$LeftDiff ++;
			}
		}
	}
	
	if($#ColA >= $#ColB)
	{
		$RightDiff = $#ColA - $#ColB;
		my $tId = $RightDiff;
		for my $i (0 .. $#ColB)
		{
			if($ColA[$i + $tId] ne $ColB[$i])
			{
				$RightDiff ++;
			}
		}
	}
	else
	{
		$RightDiff = $#ColB - $#ColA;
		my $tId = $RightDiff;
		for my $i (0 .. $#ColA)
		{
			if($ColA[$i] ne $ColB[$i + $tId])
			{
				$RightDiff ++;
			}
		}
	}
	
	$Diff = $LeftDiff;
	$Diff = $RightDiff if($RightDiff < $LeftDiff);
	
	return $Diff;
}