#!/usr/bin/perl
use strict;
use Getopt::Long;
Getopt::Long::Configure qw(no_ignore_case);
use File::Basename;
use FindBin qw($Bin);
use lib "$Bin/.Modules";
use Parameter::BinList;
use Sort::ChrPos;
use Sort::PureNum;
use BamRelated::BaseAndQual;
use VcfRelated::VcfParse;
use warnings;

my ($HelpFlag,$BinList,$BeginTime);
my ($Bam,$Dir,$LogFile,$AllFlag,$SubFlag,$SingleFlag,$BlastFlag,$RefGen,$IndelConfirmScript,$BlastScript,$Samtools);
my @PanVar;
my $HelpInfo = <<USAGE;

 CisTransOfAnnoFile_v2.pl
 Auther: zhangdong_xie\@foxmail.com

  This script was used to check the cis or trans status of potential variants.
  本脚本用于对PanCancer文件夹下snv和indel的检测结果进行潜在顺反式的自动判断。
  和v1版本相比，v2会自动核实所有可能是顺式的最长集合（配合-sub参数），同时新增了blast功能。

 -var    ( Required ) Variants (snv or indel) logging file in PanCancer directory [ multi times ];
                      PanCancer文件夹下(或者血液病文件夹下.RESULT.xls文件)的snv及indel结果文件，可以指定一次或多次；
                      支持3种格式：实体瘤流程下的pancancer文件，血液病下的RESULT文件以及vcf文件；
 -bam    ( Required ) The corresponding bam file;
                      Alignment文件夹下的比对文件，最好是.final.bam;
 -dir    ( Required ) The Directory for result logging;

 -name   ( Optional ) File name for result logging [ default:CisTrans.Result ];
 -all    ( Optional ) If there was need to confirm the truth of a single variant (with -all only);
                      在检查两个或多个变异的顺反式时，也要同时检查单个变异的真实与否。
 -sub    ( Optional ) If there was need to analyse all the subset rather than the longest one (with -sub only);
                      在多个变异组合不属于顺式时会检查其子集，直到检查到顺式为止，用于定位最小的顺式集合。
 -single ( Optional ) When checking single variants only (with '-single');
                      假如只是逐个处理每个变异，则采用‘-all’ + ‘-single’组合，它只会判断单个变异。
 -blast  ( Optional ) Blast when '-single';
                      假如需要对搜索出来的相关reads进行种属判断，则可以加‘-blast’。
 -bin    List for searching of related bin or scripts; 
 -h      Help infomation;

USAGE

GetOptions(
	'var=s' => \@PanVar,
	'bam=s' => \$Bam,
	'dir=s' => \$Dir,
	'name:s' => \$LogFile,
	'all!' => \$AllFlag,
	'sub!' => \$SubFlag,
	'single!' => \$SingleFlag,
	'blast!' => \$BlastFlag,
	'bin:s' => \$BinList,
	'h!' => \$HelpFlag
) or die $HelpInfo;

if($HelpFlag || !@PanVar || !$Bam || !$Dir)
{
	die $HelpInfo;
}
else
{
	$BeginTime = ScriptBegin();
	
	$Dir =~ s/\/$//;
	die "[ Error ] Directory not exist.\n" unless(-d $Dir);
	die "[ Error ] Bam not exist or empty ($Bam).\n" unless(-s $Bam);
	for my $i (0 .. $#PanVar)
	{
		die "[ Error ] File not exist or empty ($PanVar[$i]).\n" unless(-s $PanVar[$i]);
	}
	$LogFile = "CisTrans.Result" unless($LogFile);
	$AllFlag = "" unless($AllFlag);
	$BinList = BinListGet() if(!$BinList);
	$IndelConfirmScript = BinSearch("IndelConfirm",$BinList);
	$RefGen = BinSearch("hg19",$BinList);
	$BlastScript = BinSearch("BlastSuit",$BinList) if($BlastFlag);
	$Samtools = BinSearch("samtools",$BinList) if($BlastFlag);
}

if(1)
{
	# 准备工作;
	my $SHAndLogDir = $Dir . "/ShellAndLog";
	`mkdir -p $SHAndLogDir` unless(-d $SHAndLogDir);
	my $TmpFileDir = $Dir . "/tmpFile";
	`mkdir -p $TmpFileDir` unless(-d $TmpFileDir);
	my $BlastDir = $Dir . "/BlastResult";
	if($BlastFlag)
	{
		`mkdir -p $BlastDir` unless(-d $BlastDir);
	}
	$LogFile = $Dir . "/" . $LogFile;
	`rm $LogFile` if($LogFile && -e $LogFile);
	my $BaseName = basename $Dir;
	`echo '本次分析的Key为: $BaseName' >> $LogFile`;
	$BaseName = basename $Bam;
	my @Cols = split /\./, $BaseName;
	my $SP = $Cols[0];
	
	
	# 判断变异信息的存储格式;
	my $FileFlag = "Pan";
	if(1)
	{
		my $MatchFlag = 0;
		
		unless($MatchFlag)
		{
			my $Return = `cat $PanVar[0] | grep -v ^# | grep -v 'Chr'\$'\\t' | cut -f 1-3 | head -n 1`;
			chomp $Return;
			my ($Chr,$From,$To) = split /\t/, $Return;
			unless($Chr !~ /^chr/ && $From =~ /\D/)
			{
				$FileFlag = "Vcf";
				$MatchFlag = 1;
			}
		}
		
		unless($MatchFlag)
		{
			my $Return = `cat $PanVar[0] | grep -v ^# | grep -v 'Chr'\$'\\t' | cut -f 11-13,17-18 | head -n 1`;
			chomp $Return;
			my ($Chr,$From,$To,$Ref,$Alt) = split /\t/, $Return;
			unless($From =~ /\D/ || $To =~ /\D/ || $Ref =~ /\D/ || $Alt =~ /\D/)
			{
				$FileFlag = "Pan";
				$MatchFlag = 1;
			}
		}
		
		unless($MatchFlag)
		{
			$FileFlag = "HB";
			$MatchFlag = 1;
		}
	}
	print "[ Info ] File format: $FileFlag\n";
	
	
	# 获取变异信息;
	my (@Chr,@From,@To,@Ref,@Alt,@OtherInfo) = ();
	for my $i (0 .. $#PanVar)
	{
		if($FileFlag eq "HB")
		{
			my $Return = `cat $PanVar[$i] | grep -E ^'\\.|NM' | cut -f 2,11-15 | awk '{if(\$2 ~ /^chr/){print \$0}}'`;
			my @Items = split /\n/, $Return;
			for my $j (0 .. $#Items)
			{
				my @Cols = split /\t/, $Items[$j];
				push @Chr, $Cols[1];
				push @From, $Cols[2];
				push @To, $Cols[3];
				push @OtherInfo, join("\t",@Cols[4 .. $#Cols],$Cols[0]);
			}
		}
		elsif($FileFlag eq "Pan")
		{
			my $Return = `cut -f 11-18 $PanVar[$i] | grep -v ^# | grep -v ^'Chr'\$'\\t'`;
			my @Items = split /\n/, $Return;
			for my $j (0 .. $#Items)
			{
				my @Cols = split /\t/, $Items[$j];
				push @Chr, $Cols[0];
				push @From, $Cols[1];
				push @To, $Cols[2];
				push @OtherInfo, join("\t",@Cols[3 .. $#Cols]);
			}
		}
		elsif($FileFlag eq "Vcf")
		{
			my $Return = `cut -f 1-6 $PanVar[$i] | grep -v ^# | grep -v ^'Chr'\$'\\t'`;
			my @Items = split /\n/, $Return;
			for my $j (0 .. $#Items)
			{
				my @Cols = split /\t/, $Items[$j];
				$Cols[2] = $Cols[1] + length($Cols[3]) - 1 if($Cols[2] eq ".");
				$Cols[5] = "-" unless($Cols[5]);
				
				my @tAlt = split /,/, $Cols[4];
				for my $k (0 .. $#tAlt)
				{
					push @Chr, $Cols[0];
					push @From, $Cols[1];
					push @To, $Cols[2];
					push @OtherInfo, join("\t",$Cols[3],$tAlt[$k],$Cols[5]);
				}
			}
		}
		else
		{
			die "[ Error ] Unknow file type ($FileFlag).\n";
		}
	}
	# 排序;
	my @tRef = ChrPosAndOther(\@Chr,\@From,\@To,\@OtherInfo);
	@Chr = @{$tRef[0]};
	@From = @{$tRef[1]};
	@To = @{$tRef[2]};
	@OtherInfo = @{$tRef[3]};
	for my $i (0 .. $#OtherInfo)
	{
		my @Cols = split /\t/, $OtherInfo[$i];
		push @Ref, $Cols[0];
		push @Alt, $Cols[1];
		$OtherInfo[$i] = join("\t",@Cols[2 .. $#Cols]);
	}
	
	
	# 循环处理所有可能需要合并的变异集合;
	my $ReadLen = ReadLenConfirmFromBam($Bam);
	print "[ Info ] Read length is $ReadLen.\n";
	my $PreId = 0;
	# 防止重复检测;
	my %DupH = ();
	while($PreId <= $#Chr)
	{
		# 确定需要判断顺反式的坐标区间;
		my $CurrId = $PreId;
		unless($SingleFlag)
		{
			# 假如不需要判断顺反式则逐个变异处理;
			for my $i ($PreId + 1 .. $#Chr)
			{
				last unless($Chr[$i] eq $Chr[$PreId] && $From[$i] - $From[$PreId] <= $ReadLen);
				$CurrId ++;
			}
		}
		
		# 当该区间内变异数量大于1时才需要处理;
		if($CurrId > $PreId || $AllFlag)
		{
			# 得到互斥的最长变异组合;
			my @VarString = ();
			$VarString[0] = "," . $PreId . ",";
			for my $i ($PreId + 1 .. $CurrId)
			{
				# 相关坐标有交叉的变异是不可能共存或者顺式存在的;
				my @ConflictId = ();
				for my $j ($PreId .. $i - 1)
				{
					# 向前检查，哪些有交叉;
					push @ConflictId, $j if($From[$i] <= $To[$j]);
				}
				
				if(@ConflictId)
				{
					# 按一定规则组合成新的变异列表;
					my @tString = ();
					for my $j (0 .. $#VarString)
					{
						my $tmp = $VarString[$j];
						my $ConflictFlag = 0;
						for my $k (0 .. $#ConflictId)
						{
							$ConflictFlag = 1 if($tmp =~ s/,$ConflictId[$k],//);
						}
						# 有冲突时替换掉冲突项;
						if($tmp && $ConflictFlag)
						{
							$tmp .= "," . $i . ",";
							push @tString, $tmp;
						}
						# 没有冲突时自然延伸;
						$VarString[$j] .= "," . $i . "," unless($ConflictFlag);
					}
					# 实在没有就从它自己开始;
					$tString[0] = "," . $i . "," unless(@tString);
					
					# 合并多余的，留长不留短;
					for my $j (0 .. $#tString - 1)
					{
						for my $k ($j + 1 .. $#tString)
						{
							if($tString[$j] =~ /$tString[$k]/)
							{
								$tString[$k] = "";
							}
							elsif($tString[$k] =~ /$tString[$j]/)
							{
								$tString[$j] = "";
								last;
							}
						}
					}
					
					# 补充;
					for my $j (0 .. $#tString)
					{
						next unless($tString[$j]);
						# 再次去短留长;
						my $tmp = $tString[$j];
						$tmp =~ tr/,,/\t/;
						$tmp =~ tr/,//;
						my @Var = split /\t/, $tmp;
						my $ConflictFlag = 0;
						for my $k (0 .. $#VarString)
						{
							$ConflictFlag = 1;
							for my $z (0 .. $#Var)
							{
								$tmp = "," . $Var[$z] . ",";
								unless($VarString[$k] =~ /$tmp/)
								{
									$ConflictFlag = 0;
									last;
								}
							}
							
							last if($ConflictFlag);
						}
						next if($ConflictFlag);
						
						push @VarString, $tString[$j];
					}
				}
				else
				{
					# 没有冲突时需要打包一起判断;
					for my $j (0 .. $#VarString)
					{
						$VarString[$j] .= "," . $i . ",";
					}
				}
			}
			
			# 假如没有特殊需求则跳过只有一条变异的组合（因为一个变异没有顺反式之分）;
			unless($AllFlag)
			{
				my @tVarString = ();
				for my $i (0 .. $#VarString)
				{
					my $tString = $VarString[$i];
					$tString =~ s/,,/\t/g;
					$tString =~ s/,//g;
					my @tVar = split /\t/, $tString;
					next if($#tVar <= 0);
					
					push @tVarString, $VarString[$i];
				}
				
				@VarString = @tVarString;
			}
			
			# 对潜在的组合进行检查并记录;
			for my $i (0 .. $#VarString)
			{
				$VarString[$i] =~ s/,,/\t/g;
				$VarString[$i] =~ s/,//g;
				my @VarId = split /\t/, $VarString[$i];
				
				my @VarIdIndex = ();
				if($SubFlag)
				{
					# 得到当前组合所有可能的子集，比如5个变异不是顺式，但其中3个可能是顺式;
					my @tString = ();
					for my $j (0 .. $#VarId)
					{
						push @tString, $VarId[$j];
					}
					@VarIdIndex = @{AllSubSet(\@tString,0)};
				}
				else
				{
					# 默认情况下只分析当前这个最长的集合;
					for my $j (0 .. $#VarId)
					{
						push @{$VarIdIndex[0]}, $VarId[$j];
					}
				}
				
				
				# 从最长的子集开始检查;
				for my $k (0 .. $#VarIdIndex)
				{
					# 当检出顺式之后，剩下子集就没有检测的必要了;
					my $JumpFlag = 0;
					foreach my $Key (keys %DupH)
					{
						# 顺式会被标记;
						next unless($DupH{$Key} > 1);
						
						$JumpFlag = 1;
						for my $j (0 .. $#{$VarIdIndex[$k]})
						{
							my $tmp = "," . $VarIdIndex[$k][$j] . ",";
							if($Key !~ /$tmp/)
							{
								# 有冲突则肯定不是子集;
								$JumpFlag = 0;
								last;
							}
						}
						# 是子集则下一个循环;
						last if($JumpFlag);
					}
					next if($JumpFlag);
					
					my @TransVar = ();
					for my $j (0 .. $#{$VarIdIndex[$k]})
					{
						push @TransVar, join(",",$Chr[$VarIdIndex[$k][$j]],$From[$VarIdIndex[$k][$j]],$To[$VarIdIndex[$k][$j]],$Ref[$VarIdIndex[$k][$j]],$Alt[$VarIdIndex[$k][$j]]);
					}
					# 确认待检测的变异集合是否分析过，避免重复;
					my $DupHKey = "," . join(",",@{$VarIdIndex[$k]}) . ",";
					next if($DupH{$DupHKey});
					
					my $Prefix = "CisTransConfirm." . $SP . "." . $Chr[$PreId] . "." . $From[$PreId] . "." . $From[$CurrId] . "." . $i . "." . $k;
					$Prefix = "CisTransConfirm." . $SP . "." . $Chr[$CurrId] . "." . $From[$CurrId] . "." . $To[$CurrId] . "." . $Alt[$CurrId] if($SingleFlag);
					my $ResultBam = $TmpFileDir . "/" . $Prefix . ".bam";
					my $Shell = $SHAndLogDir . "/" . $Prefix . ".sh";
					my $Log = $SHAndLogDir . "/" . $Prefix . ".sh.log";
					
					open(SH,"> $Shell") or die $!;
					print SH "perl $IndelConfirmScript \\\n";
					print SH "\t-b $Bam \\\n";
					print SH "\t-o $ResultBam \\\n";
					print SH "\t-extend 1,1 \\\n";
					print SH "\t-v '",join("' \\\n\t-v '",@TransVar),"'";
					close SH;
					`chmod 750 $Shell`;
					
					my $ErrorInfo;
					$ErrorInfo = `sh $Shell > $Log`;
					chomp $ErrorInfo;
					die "[ Error ] Something went wrong when running $Shell.\n" if($ErrorInfo);
					
					# blast;
					if($BlastFlag)
					{
						my $tBltDir = $BlastDir . "/" . $SP . "." . $Chr[$CurrId] . "." . $From[$CurrId] . "." . $To[$CurrId] . "." . $Alt[$CurrId];
						`mkdir -p $tBltDir` unless(-d $tBltDir);
						my $BltFa = $tBltDir . "/" . $SP . "." . $Chr[$CurrId] . "." . $From[$CurrId] . "." . $To[$CurrId] . "." . $Alt[$CurrId] . ".fa";
						`$Samtools view $ResultBam | awk '{print ">"\$1"\\n"\$10}' > $BltFa`;
						`perl $BlastScript -fq $BltFa -dir $tBltDir`;
					}
					
					# 对分析过的变异集合做好标记;
					$DupH{$DupHKey} = 1;
					
					# 记录检测结果;
					my $Return = '';
					$Return = `grep 'Total reads that match all these vars are' $Log` if($Log && -e $Log);
					chomp $Return;
					if($Return)
					{
						# 顺式集合做好标记;
						my @Items = split /,/, $Return;
						die "[ Error ] The number of alt was not pure number ($Items[-1]).\n" if($Items[-1] =~ /\D/);
						$DupH{$DupHKey} = 10 if($Items[-1] >= 20);
						
						my $tLogInfo = "\n# 可能需要分析或合并的变异信息:\n";
						$tLogInfo .= join("\t","# [ . ]","Chr","From","To","Ref","Alt","Freq","RefDepth","AltDepth") . "\n" if($FileFlag eq "Pan");
						$tLogInfo .= join("\t","# [ . ]","Chr","From","To","Ref","Alt","Freq") . "\n" if($FileFlag eq "HB");
						$tLogInfo .= join("\t","# [ . ]","Chr","From","To","Ref","Alt","Other") . "\n" if($FileFlag eq "Vcf");
						for my $j (0 .. $#{$VarIdIndex[$k]})
						{
							my $tId = $j + 1;
							$tLogInfo .= "# [ $tId ]\t" . join("\t",$Chr[$VarIdIndex[$k][$j]],$From[$VarIdIndex[$k][$j]],$To[$VarIdIndex[$k][$j]],$Ref[$VarIdIndex[$k][$j]],$Alt[$VarIdIndex[$k][$j]],$OtherInfo[$VarIdIndex[$k][$j]]) . "\n";
						}
						$tLogInfo .= "# 结果(这里只分析了覆盖所有变异的reads，所以总深度可能较低):\n";
						$tLogInfo .= $Return . "\n";
						`echo '$tLogInfo' >> $LogFile`;
					}
				}
			}
			
			$CurrId -- if($CurrId > $PreId);
		}
		
		# 准备下一个循环;
		$PreId = $CurrId + 1;
	}
	
	# 记录结束标记;
	`echo 'All Done' >> $LogFile`;
}
printf "[ %.2fmin ] All Done.\n",(time - $BeginTime)/60;


######### Sub functions ##########
