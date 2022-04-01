# Package Name
package SeqRelated::Seq;

# Exported name
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(RefGet RefGetWithOutWarning);

use FindBin qw($Bin);

sub RefGet
{
	my ($RefGen,$Chr,$From,$To,$DefaultBedtools) = @_;
	
	die "[ Error ] End position can not be smaller than start position when getting seq ($Chr,$From,$To) from ref($RefGen).\n" if($From > $To);
	die "[ Error ] Reference not exist ($RefGen).\n" unless(-e $RefGen);
	
	$From --;
	$DefaultBedtools = "bedtools" unless($DefaultBedtools);
	my $Ref = `echo -e '$Chr\\t$From\\t$To' | $DefaultBedtools getfasta -fi $RefGen -bed - | tail -n 1`;
	chomp $Ref;
	$Ref = uc($Ref);
	
	return $Ref;
}

sub RefGetWithOutWarning
{
	my ($RefGen,$Chr,$From,$To,$DefaultBedtools) = @_;
	
	my $Ref = "";
	return $Ref if($From > $To);
	die "[ Error ] Reference not exist ($RefGen).\n" unless(-e $RefGen);
	
	$From --;
	$DefaultBedtools = "bedtools" unless($DefaultBedtools);
	$Ref = `echo -e '$Chr\\t$From\\t$To' | $DefaultBedtools getfasta -fi $RefGen -bed - | tail -n 1`;
	chomp $Ref;
	$Ref = uc($Ref);
	
	return $Ref;
}

1;