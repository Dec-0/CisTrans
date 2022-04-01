# Package Name
package Sort::PureNum;

# Exported name
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(PureNumSort AllSubSet);

sub PureNumSort
{
	my @Ori = @{$_[0]};
	
	my @Index = ();
	my ($Id1,$Id2) = (0,1);
	for my $i (0 .. $#Ori)
	{
		push @{$Index[$Id1]}, $i;
	}
	my $tmp = @Ori;
	my $To = $#Ori;
	my $MaxDulSpan = $tmp * 2;
	for(my $DulSpan = 2;$DulSpan < $MaxDulSpan;$DulSpan = $DulSpan * 2)
	{
		my $MinSpan = $DulSpan / 2;
		my $tmpId = 0;
		
		for(my $i = 0;$i <= $To;$i += $DulSpan)
		{
			my $LeftBegin = $i;
			my $RightBegin = $LeftBegin + $MinSpan;
			if($RightBegin <= $To)
			{
				my $LeftEnd = $RightBegin - 1;
				my $RightEnd = $LeftEnd + $MinSpan;
				if($RightEnd > $To)
				{
					$RightEnd = $To;
				}
				
				while($LeftBegin <= $LeftEnd || $RightBegin <= $RightEnd)
				{
					if($LeftBegin > $LeftEnd)
					{
						$Index[$Id2][$tmpId] = $Index[$Id1][$RightBegin];
						$RightBegin ++;
					}
					elsif($RightBegin > $RightEnd)
					{
						$Index[$Id2][$tmpId] = $Index[$Id1][$LeftBegin];
						$LeftBegin ++;
					}
					else
					{
						if($Ori[$Index[$Id1][$LeftBegin]] > $Ori[$Index[$Id1][$RightBegin]])
						{
							$Index[$Id2][$tmpId] = $Index[$Id1][$RightBegin];
							$RightBegin ++;
						}
						else
						{
							$Index[$Id2][$tmpId] = $Index[$Id1][$LeftBegin];
							$LeftBegin ++;
						}
					}
					$tmpId ++;
				}
			}
			else
			{
				for(my $j = $LeftBegin;$j <= $To;$j ++)
				{
					$Index[$Id2][$tmpId] = $Index[$Id1][$j];
					$tmpId ++;
				}
			}
		}
		$tmp = $Id2;
		$Id2 = $Id1;
		$Id1 = $tmp;
	}
	
	my @tOri = ();
	for my $i (0 .. $#{$Index[$Id1]})
	{
		push @tOri, $Ori[$Index[$Id1][$i]];
	}
	
	return \@tOri;
}

sub AllSubSet
{
	my @OriString = @{$_[0]};
	my $LenFlag = $_[1];
	$LenFlag = 0 unless($LenFlag);
	my @SubSet = ();
	my $MaxLen = @OriString;
	die "[ Error ] Exceed the maximal number of items $LenFlag vs. $MaxLen .\n" if($LenFlag > $MaxLen);
	
	for my $i (0 .. $#OriString)
	{
		my @tSubSet = @SubSet;
		for my $j (0 .. $#tSubSet)
		{
			$tSubSet[$j] .= "\t" . $OriString[$i];
		}
		push @SubSet, @tSubSet;
		push @SubSet, $OriString[$i];
	}
	
	if($LenFlag)
	{
		my @tSubSet = ();
		for my $i (0 .. $#SubSet)
		{
			my @Items = split /\t/, $SubSet[$i];
			my $Len = @Items;
			if($Len == $LenFlag)
			{
				push @tSubSet, $SubSet[$i];
			}
		}
		@SubSet = @tSubSet;
	}
	else
	{
		my @tSubSet = ();
		for(my $i = $MaxLen; $i >= 1; $i --)
		{
			for my $j (0 .. $#SubSet)
			{
				my @Items = split /\t/, $SubSet[$j];
				my $Len = @Items;
				if($Len == $i)
				{
					push @tSubSet, $SubSet[$j];
				}
			}
		}
		@SubSet = @tSubSet;
	}
	
	if(1)
	{
		my @tSubSet = ();
		for my $i (0 .. $#SubSet)
		{
			my @Items = split /\t/, $SubSet[$i];
			@{$tSubSet[$i]} = @Items;
		}
		@SubSet = @tSubSet;
	}
	
	return \@SubSet;
}

1;