# Package Name
package VcfRelated::VcfParse;

# Exported name
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(VarSimplify VarUniform RefConfirm VarTrans);
use SeqRelated::Seq;

# 去除多余的碱基;
sub VarSimplify
{
	my ($Chr,$From,$To,$Ref,$Alt) = @_;
	
	if($Ref =~ /^$Alt/)
	{
		# deletion;
		$From += length($Alt);
		$Ref =~ s/^$Alt//;
		$Alt = "-";
	}
	elsif($Alt =~ /^$Ref/)
	{
		# insertion;
		$From = $To;
		$Alt =~ s/^$Ref//;
		$Ref = "-";
	}
	
	return $Chr,$From,$To,$Ref,$Alt;
}

# 将所有的变异都尽量往5'移动;
sub VarUniform
{
	my ($RefGen,$Chr,$From,$To,$Ref,$Alt,$DefaultBedtools) = @_;
	
	($Chr,$From,$To,$Ref,$Alt) = &VarSimplify($Chr,$From,$To,$Ref,$Alt);
	if($Ref && $Alt eq "-")
	{
		# deletion;
		my $Len = length($Ref);
		
		my $LeftBase = RefGet($RefGen,$Chr,$From - 1,$From - 1,$DefaultBedtools);
		my $RightBase = RefGet($RefGen,$Chr,$To,$To,$DefaultBedtools);
		while($LeftBase eq $RightBase)
		{
			$From --;
			$To --;
			$LeftBase = RefGet($RefGen,$Chr,$From - 1,$From - 1,$DefaultBedtools);
			$RightBase = RefGet($RefGen,$Chr,$To,$To,$DefaultBedtools);
		}
		$Ref = RefGet($RefGen,$Chr,$From,$To,$DefaultBedtools);
	}
	elsif($Ref eq "-" && $tAlt)
	{
		# insertion;
		my $AltLen = length($Alt);
		my @AltBase = split //, $Alt;
		
		my $tPos = $From;
		my $tBase = RefGet($RefGen,$Chr,$tPos,$tPos,$DefaultBedtools);
		my $tNum = 0;
		my $tId = $#AltBase - ($tNum % $AltLen);
		while($tBase eq $AltBase[$tId])
		{
			$tPos --;
			$tBase = RefGet($RefGen,$Chr,$tPos,$tPos,$DefaultBedtools);
			$tNum ++;
			$tId = $#AltBase - ($tNum % $AltLen);
		}
		
		if($tPos < $From)
		{
			$tBase = RefGet($RefGen,$Chr,$tPos + 1,$From,$DefaultBedtools);
			$tBase .= $Alt;
			$Alt = substr($tBase,0,$AltLen);
			$From = $tPos;
			$To = $From;
		}
	}
	
	return $Chr,$From,$To,$Ref,$Alt;
}

# 确认ref序列是否正确;
sub RefConfirm
{
	my ($RefGen,$Chr,$From,$To,$Ref,$DefaultBedtools) = @_;
	my $Flag = 1;
	
	die "[ Error in VcfRelated::VcfParse ] To smaller than From ($Chr,$From,$To,$Ref).\n" if($To =~ /\D/ || $To < $From);
	if($Ref ne "-")
	{
		$Flag = 0 unless(length($Ref) == $To - $From + 1);
		
		if($Flag)
		{
			my $tSeq = RefGet($RefGen,$Chr,$From,$To,$DefaultBedtools);
			$Flag = 0 if($tSeq ne $Ref);
		}
	}
	
	return $Flag;
}

# 转换成顺反式判断统一的格式;
sub VarTrans
{
	my ($RefGen,$Chr,$From,$To,$Ref,$Alt,$DefaultBedtools) = @_;
	
	die "[ Error ] Ref seq not correct for $Chr,$From,$To,$Ref (Reference: $RefGen).\n" unless(&RefConfirm($RefGen,$Chr,$From,$To,$Ref,$DefaultBedtools));
	($Chr,$From,$To,$Ref,$Alt) = &VarUniform($RefGen,$Chr,$From,$To,$Ref,$Alt,$DefaultBedtools);
	
	$From -- unless($Ref eq "-");
	my $TransVar = join(",",$Chr,$From,$Ref,$Alt);
	
	return $TransVar;
}

1;