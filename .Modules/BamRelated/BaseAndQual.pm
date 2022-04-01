# Package Name
package BamRelated::BaseAndQual;

# Exported name
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(ReadLenConfirmFromBam ReadLenConfirmFromFq CigarSplit CigarLenCheck ReadCoverRange BamItemGet MDSplit AltGet);

sub ReadLenConfirmFromBam
{
	my $Bam = $_[0];
	my $ReadLen = 0;
	
	die "[ Error ] Bam file not exist ($Bam).\n" unless($Bam && -s $Bam);
	my $Return = `samtools view -F 0x900 $Bam | head -n 100 | cut -f 10`;
	my @Read = split /\n/, $Return;
	for my $i (0 .. $#Read)
	{
		my $tLen = length($Read[$i]);
		if($tLen > $ReadLen)
		{
			$ReadLen = $tLen;
		}
	}
	
	return $ReadLen;
}

sub ReadLenConfirmFromFq
{
	my $Fq = $_[0];
	my $ReadLen = 0;
	
	die "[ Error ] Fastq not exist ($Fq).\n" unless($Fq && -s $Fq);
	my $Return;
	if($Fq =~ /.gz$/)
	{
		$Return = `zcat $Fq | head -n 100 | awk '{if(NR < 40 && NR % 4 == 2){print \$0}}'`;
	}
	else
	{
		$Return = `cat $Fq | head -n 100 | awk '{if(NR < 40 && NR % 4 == 2){print \$0}}'`;
	}
	my @Read = split /\n/, $Return;
	for my $i (0 .. $#Read)
	{
		my $tLen = length($Read[$i]);
		if($tLen > $ReadLen)
		{
			$ReadLen = $tLen;
		}
	}
	
	return $ReadLen;
}

sub CigarSplit
{
	my $Cigar = $_[0];
	
	my @Items = ();
	my $ItemId = 0;
	while($Cigar =~ s/^(\d+)(\D)//)
	{
		$Items[$ItemId][0] = $2;
		$Items[$ItemId][1] = $1;
		
		die "[ Error ] Non-single characer in cigar string ($Cigar,$Items[$ItemId][0]).\n" unless(length($Items[$ItemId][0]) == 1);
		die "[ Error ] Non Capital in Cigar ($Cigar,$Items[$ItemId][0]).\n" unless(ord($Items[$ItemId][0]) >= 65 && ord($Items[$ItemId][0]) <= 90);
		
		$ItemId ++;
	}
	
	return @Items;
}

sub CigarLenCheck
{
	my ($Cigar,$SeqLen) = @_;
	
	my @Items = &CigarSplit($Cigar);
	my $CigarLen = 0;
	for my $i (0 .. $#Items)
	{
		$CigarLen += $Items[$i][1] unless($Items[$i][0] eq "H" || $Items[$i][0] eq "D" || $Items[$i][0] eq "N" || $Items[$i][0] eq "P");
	}
	my $Diff = $CigarLen - $SeqLen;
	
	return $Diff;
}

# 用于判断bam中的一行reads所覆盖的范围，返回基因组坐标;
sub ReadCoverRange
{
	my ($Line,$SclipFlag) = @_;
	$SclipFlag = 0 unless($SclipFlag);
	
	chomp $Line;
	$Line =~ s/\s+/\t/g unless($Line =~ /\t/);
	my @Cols = split /\t/, $Line;
	my ($Chr,$Pos,$Cigar,$Seq) = ($Cols[2],$Cols[3],$Cols[5],$Cols[9]);
	die "[ Error ] Cigar not match with the length of the read's sequence ($Cigar, $Seq).\n" unless(&CigarLenCheck($Cigar,length($Seq)) == 0);
	my ($From,$To) = ($Pos,$Pos - 1);
	my @Items = &CigarSplit($Cigar);
	for my $i (0 .. $#Items)
	{
		next if($Items[$i][0] eq "H" || $Items[$i][0] eq "I");
		next if(!$SclipFlag && $Items[$i][0] eq "S");
		
		if($Items[$i][0] eq "S" && $i == 0)
		{
			$From -= $Items[$i][1];
		}
		else
		{
			$To += $Items[$i][1];
		}
	}
	# coordinate should not minus 1;
	# $From = 1 if($From < 1);
	
	return $Chr,$From,$To;
}

sub BamItemGet
{
	my ($Line,$Item) = @_;
	
	my $Value = "";
	chomp $Line;
	if($Line =~ /\t($Item[^\t\s]+)/)
	{
		$Value = $1;
	}
	elsif($Line =~ /\s($Item[^\t\s]+)/)
	{
		$Value = $1;
	}
	elsif($Line =~ /^($Item[^\t\s]+)/)
	{
		$Value = $1;
	}
	else
	{
		die "[ Error ] No item named $Item in $Line.\n";
	}
	my @Cols = split /:/, $Value;
	$Value = $Cols[-1];
	
	return $Value;
}

sub MDSplit
{
	my $String = $_[0];
	
	my @Value = ();
	while($String)
	{
		if($String=~ s/^(\d+)//)
		{
			push @Value, $1;
		}
		elsif($String=~ s/^(\D+)//)
		{
			push @Value, $1;
		}
	}
	
	return @Value;
}

sub AltGet
{
	my ($Chr,$From,$To,$Line,$SclipFlag) = @_;
	$SclipFlag = 0 unless($SclipFlag);
	
	# make sure it was covered by this read;
	my ($cChr,$cFrom,$cTo) = &ReadCoverRange($Line,$SclipFlag);
	my $AltSeq = "";
	return $AltSeq if($cChr ne $Chr || $From > $cTo || $To < $cFrom);
	
	# initial and cigar assignment;
	chomp $Line;
	$Line =~ s/\s+/\t/g unless($Line =~ /\t/);
	my @Cols = split /\t/, $Line;
	my ($Chr,$Pos,$Cigar,$Seq) = ($Cols[2],$Cols[3],$Cols[5],$Cols[9]);
	my @Base = split //, $Seq;
	my @Items = &CigarSplit($Cigar);
	my @Alt = ();
	my ($AltId,$BaseId,$SclipLen) = (0,0,0);
	for my $i (0 .. $#Items)
	{
		next if($Items[$i][0] eq "H");
		
		if($Items[$i][0] eq "M" || $Items[$i][0] eq "=" || $Items[$i][0] eq "X")
		{
			# match or mis-match;
			for my $j (1 .. $Items[$i][1])
			{
				$Alt[$AltId][0] = $Items[$i][0];
				$Alt[$AltId][1] = $Base[$BaseId];
				$Alt[$AltId][2] = "";
				$AltId ++;
				$BaseId ++;
			}
		}
		elsif($Items[$i][0] eq "D" || $Items[$i][0] eq "N" || $Items[$i][0] eq "P")
		{
			# deletion;
			for my $j (1 .. $Items[$i][1])
			{
				$Alt[$AltId][0] = $Items[$i][0];
				$Alt[$AltId][1] = "";
				$Alt[$AltId][2] = "";
				$AltId ++;
			}
		}
		elsif($Items[$i][0] eq "I")
		{
			# insertion;
			$Alt[$AltId - 1][2] = $Items[$i][0];
			for my $j (1 .. $Items[$i][1])
			{
				$Alt[$AltId - 1][1] .= $Base[$BaseId];
				$BaseId ++;
			}
		}
		elsif($Items[$i][0] eq "S")
		{
			# soft-clip;
			for my $j (1 .. $Items[$i][1])
			{
				$Alt[$AltId][0] = $Items[$i][0];
				$Alt[$AltId][1] = $Base[$BaseId] if($SclipFlag);
				$Alt[$AltId][1] = "" unless($SclipFlag);
				$Alt[$AltId][2] = "";
				$AltId ++;
				$BaseId ++;
				$SclipLen ++ unless($SclipFlag);
			}
		}
		else
		{
			die "[ Error ] Unknown Cigar character $Items[$i][0] in $Cigar ($Line)\n";
		}
	}
	die "[ Error ] Cigar spliting not correct ($Line).\n" unless($AltId - $SclipLen == $cTo - $cFrom + 1 && length($Seq) == $BaseId);
	
	# ----------------------------------------------------------------------------------------
	# revising when indel was followed by a snv (assign snv with the help of MD value in bam);
	my $MDString = &BamItemGet($Line,"MD:Z");
	my @MD = &MDSplit($MDString);
	my $MDId = 0;
	# jump off sclip;
	for my $i (0 .. $AltId - 1)
	{
		last unless($Alt[$i][0] eq "S");
		$MDId ++;
	}
	# 调整因slcipflag导致的坐标移位;
	unless($SclipFlag)
	{
		$cFrom = $cFrom - $MDId;
		$cTo = $cTo + $SclipLen - $MDId;
	}
	die "[ Error ] Cigar ($cFrom -> $cTo) not match with the length of read cover ($Line).\n" unless($cTo - $cFrom + 1 == $AltId);
	for my $i (0 .. $#MD)
	{
		# no insertion, ignore ...;
		if($MD[$i] =~ /^\^/)
		{
			my $String = substr($MD[$i],1);
			my @RefBase = split //, $String;
			for my $j (0 .. $#RefBase)
			{
				die "[ Error ] Not deletion, MD not matching with Cigar ($MD[$i] in $MDString vs $Alt[$MDId][0] in $Cigar).\n" unless($Alt[$MDId][0] eq "D");
				
				$Alt[$MDId][3] = $RefBase[$j];
				$MDId ++;
			}
		}
		elsif($MD[$i] =~ /\D/)
		{
			# mis-match;
			my @RefBase = split //, $MD[$i];
			for my $j (0 .. $#RefBase)
			{
				die "[ Error ] Not non-matching, MD not matching with Cigar ($MD[$i] in $MDString vs $Alt[$MDId][0] in $Cigar).\n" unless($Alt[$MDId][0] eq "M" || $Alt[$MDId][0] eq "=" || $Alt[$MDId][0] eq "X");
				
				$Alt[$MDId][0] = "X";
				$Alt[$MDId][3] = $RefBase[$j];
				$MDId ++;
			}
		}
		else
		{
			# match;
			for my $j (1 .. $MD[$i])
			{
				die "[ Error ] Not matching, MD not matching with Cigar ($MD[$i] in $MDString vs $Alt[$MDId][0] in $Cigar).\n" unless($Alt[$MDId][0] eq "M" || $Alt[$MDId][0] eq "=" || $Alt[$MDId][0] eq "X");
				
				$MDId ++;
			}
		}
	}
	die "[ Error ] MD string and cigar are not consistent ($MDId vs $AltId).\n" unless($MDId == $AltId || $Alt[$MDId][0] eq "S");
	# 将上下游最靠近该区域的1个indel尽量往该区域移动;
	# 根据bwa的特点，考虑5'就好;
	my $SpaceFlag = 1;
	while($SpaceFlag)
	{
		$SpaceFlag = 0;
		for(my $i = $From - $cFrom;$i > 0;$i --)
		{
			# indel + match;
			my $MatchFlag = 0;
			if($Alt[$i][0] eq "M" || $Alt[$i][0] eq "=")
			{
				if($Alt[$i - 1][0] eq "D")
				{
					my $DelBegin = $i - 1;
					for(my $j = $i - 2;$j > 0;$j --)
					{
						last unless($Alt[$j][0] eq "D");
						$DelBegin --;
					}
					
					if($Alt[$DelBegin][3] eq $Alt[$i][1])
					{
						$Alt[$DelBegin][0] = $Alt[$i][0];
						$Alt[$DelBegin][1] = $Alt[$DelBegin][3];
						
						$Alt[$i][0] = "D";
						$Alt[$i][3] = $Alt[$i][1];
						$Alt[$i][1] = "";
						
						$SpaceFlag = 1;
					}
					
					$MatchFlag = 1;
				}
				elsif($Alt[$i - 1][2] eq "I")
				{
					my $tBase = substr($Alt[$i - 1][1],1,1);
					if($tBase eq $Alt[$i][1])
					{
						$Alt[$i][1] = substr($Alt[$i - 1][1],1) . $Alt[$i][1];
						$Alt[$i][2] = "I";
						
						$Alt[$i - 1][1] = substr($Alt[$i - 1][1],0,1);
						$Alt[$i - 1][2] = "";
						
						$SpaceFlag = 1;
					}
					
					$MatchFlag = 1;
				}
			}
			
			last if($MatchFlag);
		}
	}
	
	# final;
	my $tFrom = $From - $cFrom;
	$tFrom = 0 if($tFrom < 0);
	my $tTo = $To - $cFrom;
	$tTo = $cTo if($tTo >= $AltId);
	#print "($tFrom - $tTo)\n";
	for my $i ($tFrom .. $tTo)
	{
		$AltSeq .= $Alt[$i][1];
	}
	
	# in case normal seq treated as insert + del;
	if(length($Alt[$tFrom - 1][1]) > 1)
	{
		$AltSeq = substr($Alt[$tFrom - 1][1],1) . $AltSeq;
	}
	
	return $AltSeq;
}

1;