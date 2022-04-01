# Package Name
package Parameter::BinList;

# Exported name
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(BinListGet BinSearch DayStamp TimeStamp DayTimeStamp ScriptBegin RandString);

use FindBin qw($Bin);

sub BinListGet
{
	my $BinList = $Bin . "/.BinList.xls";
	return $BinList if(-s $BinList);
	
	$BinList = "./Scripts/.BinList.xls";
	unless(-s $BinList)
	{
		die "[ Error ] File not exist ($BinList).\n";
	}
	
	return $BinList;
}

sub BinSearch
{
	my ($BinName,$BinList) = @_;
	
	die "[ Error ] File not exist ($BinList).\n" unless(-s $BinList);
	my $Return = `grep '^$BinName'\$'\\t' $BinList`;
	chomp $Return;
	
	die "Name ($BinName) not unique.\n($Return)\n" if($Return =~ /\n/);
	#print "[ Info ] $Return\n";
	my @Cols = split /\t/, $Return;
	# 替换成Bin目录
	if($Cols[1] =~ s/^BINDIR//)
	{
		$Cols[1] = $Bin . $Cols[1];
	}
	# 替换成Scripts目录
	elsif($Cols[1] =~ s/^DEFAULTDIR//)
	{
		my @Items = split /\//, $Bin;
		my $LastId = -1;
		for($i = $#Items;$i >= 0;$i --)
		{
			if($Items[$i] eq "Scripts")
			{
				$LastId = $i - 1;
				last;
			}
			elsif($Items[$i] eq "xiezhangdong")
			{
				$LastId = $i;
				last;
			}
		}
		
		if($LastId >= 0)
		{
			my $RDir = join("/",@Items[0 .. $LastId]);
			$Cols[1] = $RDir . $Cols[1];
		}
	}
	die "[ Error ] $Cols[1] for $BinName not exist!\n" unless($Cols[1] && -e $Cols[1]);
	
	return $Cols[1];
}

sub DayStamp
{
	my @temp_time = localtime();
	my $localtime_year = $temp_time[5] + 1900;
	my $localtime_month = $temp_time[4] + 1;
	my $DayStamp = $localtime_year . "/" . $localtime_month . "/" . $temp_time[3];
	
	return $DayStamp;
}

sub TimeStamp
{
	my @temp_time = localtime();
	$TimeStamp = $temp_time[2] . ":" . $temp_time[1] . ":" . $temp_time[0];
	
	return $TimeStamp;
}

sub DayTimeStamp
{
	my $DayStamp = &DayStamp();
	my $TimeStamp = &TimeStamp();
	my $Stamp = $DayStamp . " " . $TimeStamp;
	
	return $Stamp;
}

sub ScriptBegin
{
	$BeginTime = time;
	if($BeginTime)
	{
		my $Stamp = &DayTimeStamp();
		print "[ $Stamp ] This script begins.\n";
	}
	else
	{
		die "[ Error ] Fail to get system time with 'time'.\n";
	}
	
	return $BeginTime;
}

sub RandString
{
	my $Len = $_[0];
	
	die "[ Error ] The length of rand string is not pure number ($Len).\n" if($Len =~ /\D/);
	my @Items = ();
	for my $i (0 .. 9)
	{
		push @Items, $i;
	}
	for my $i (65 .. 90)
	{
		push @Items, chr($i);
	}
	for my $i (97 .. 122)
	{
		push @Items, chr($i);
	}
	my $Full = @Items;
	my $String = "";
	for my $i (1 .. $Len)
	{
		$String .= $Items[int(rand($Full))];
	}
	
	return $String;
}

1;