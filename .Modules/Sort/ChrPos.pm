# Package Name
package Sort::ChrPos;

# Exported name
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(ChrPosAndOther);

sub ChrPosAndOther
{
	my @tRef = @_;
	die "[ Error ] ChrPosAndOther requires 3 parameters.\n" unless($#tRef == 3);
	my @Chr = @{$tRef[0]};
	my @From = @{$tRef[1]};
	my @To = @{$tRef[2]};
	my @OtherInfo = @{@tRef[3]};
	
	my @Index = ();
	my ($Id1,$Id2) = (0,1);
	for my $i (0 .. $#Chr)
	{
		push @{$Index[$Id1]}, $i;
	}
	my $tmp = @Chr;
	my $To = $#Chr;
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
						if($Chr[$Index[$Id1][$LeftBegin]] gt $Chr[$Index[$Id1][$RightBegin]] || ($Chr[$Index[$Id1][$LeftBegin]] eq $Chr[$Index[$Id1][$RightBegin]] && $From[$Index[$Id1][$LeftBegin]] > $From[$Index[$Id1][$RightBegin]]) || ($Chr[$Index[$Id1][$LeftBegin]] eq $Chr[$Index[$Id1][$RightBegin]] && $From[$Index[$Id1][$LeftBegin]] == $From[$Index[$Id1][$RightBegin]] && $To[$Index[$Id1][$LeftBegin]] > $To[$Index[$Id1][$RightBegin]]))
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
	
	my (@FinalChr,@FinalFrom,@FinalTo,@FinalOther) = ();
	for my $i (0 .. $#{$Index[$Id1]})
	{
		push @FinalChr, $Chr[$Index[$Id1][$i]];
		push @FinalFrom, $From[$Index[$Id1][$i]];
		push @FinalTo, $To[$Index[$Id1][$i]];
		push @FinalOther, $OtherInfo[$Index[$Id1][$i]];
	}
	
	return \@FinalChr,\@FinalFrom,\@FinalTo,\@FinalOther;
}

1;