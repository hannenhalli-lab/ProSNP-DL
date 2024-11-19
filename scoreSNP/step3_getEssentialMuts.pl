#!/usr/bin/perl -w
#use List::Utils qw(min max sum);
use strict;

my $file = "slidingWindowOverlapHighFstSNP.withDLscore.example";###slidingWindowOverlapHighFstSNP.withDLscore.example is the example snp rs10095018 showing the format of the input

my $threshold = 0.58;#### enhancer model FPR = 0.01
my %snpWindow2refScore = ();
my %snpWindow2mutScore = ();
my @lines = `cat $file`;
foreach my $line(@lines){
	chomp $line;
	my @spline = split(/\t/, $line);
	$spline[0]=~/Dis\:(\S+)$/ or die("windowDis: $spline[0]\n");
	my $window = $1 + 1 - 1;
	next unless($spline[1]=~/rs\d+/);
	$spline[1]=~/(rs\d+)\-(\S+)Ale\:(\S+)$/ or die("2nd element: $spline[1]\n");
	my $snpID = $1;

	$spline[-1]=~/\[(\S+)\]/ or die("last element: $spline[-1]\n");
	my $score = $1;

	if($spline[1]=~/refAle/){
		$snpWindow2refScore{$snpID}{$window} = $score;
	}else{
		$spline[1]=~/Ale\:(\S+)$/ or die("aleInfo: $spline[2]\n");
		my $mut = $1;
		my $mutScore = $mut."|".$score;
		push(@{$snpWindow2mutScore{$snpID}{$window}}, $mutScore);
		 
	}
	
}

my %gainMuts = ();### 
my %lossMuts = ();#### 
my %essentialGainMuts = ();### from zero to one
my %essentialLossMuts = ();#### from one to zero
my %nonEssentialGainMuts = ();### from zero to zero
my %nonEssentialLossMuts = ();#### from one to one

foreach my $snpID(keys(%snpWindow2refScore)){
	my %window2refScore = ();
	my %window2mutScore = ();
	
	foreach my $window(keys(%{$snpWindow2refScore{$snpID}})){
		my $refScore = $snpWindow2refScore{$snpID}{$window};
		unless(exists($snpWindow2mutScore{$snpID})){
			print ("SNP $snpID does not exist in snpWindow2mutScore hash\n");
			next;
		}
		unless(exists($snpWindow2mutScore{$snpID}{$window})){
			print ("window: $window does not exist in snpWindow2mutScore:$snpID hash\n");
			next;
		}
		my @mutScores = @{$snpWindow2mutScore{$snpID}{$window}};
		foreach my $mutScore(@mutScores){
			my @tmp = split(/\|/, $mutScore);
			my $mut = $tmp[0];
			my $score = $tmp[1];
			$window2refScore{$window} = $refScore;
			$window2mutScore{$mut}{$window} = $score;
		}
	}
	
	foreach my $mut(keys(%window2mutScore)){
		my $sum = 0;
		my @weights = ();
		my $weightSum = 0;
		my $positiveCount = 0;### positive delta
		my $negativeCount = 0;
		my $essentialPositiveCount = 0;## one mut chagen 0 to 1
		my $essentialNegativeCount = 0;## one mut change 1 to 0
		my $windowCount = 0;
		foreach my $window(keys(%window2refScore)){
			$windowCount ++;
			my $absDis = abs($window);
			my $weight;
			if($absDis <= 200){
				$weight = 1;
			}else{
				$weight = exp(-$absDis);
			}
			push(@weights, $weight);
			$weightSum += $weight;
			my $refScore = $window2refScore{$window};
			my $mutScore = $window2mutScore{$mut}{$window};
			my $delta = $mutScore - $refScore;
			$sum += $weight*$delta;
			
			if($absDis <= 200){
				if(($refScore < $threshold)&&($mutScore >= $threshold)){
					$essentialPositiveCount ++;
	#				print "$snpID\tmut:$mut\tfrom 0 to 1\n";
				}
				if(($refScore >= $threshold)&&($mutScore < $threshold)){
					$essentialNegativeCount ++;
	#				print "$snpID\tmut:$mut\tfrom 1 to 0\n";
				}
				$positiveCount ++ if($delta > 0);
				$negativeCount ++ if($delta < 0);
			}
		}
		my $average = $sum/$weightSum;
		
		my $positiveFrac = $positiveCount/$windowCount;
		my $negativeFrac = $negativeCount/$windowCount;
	#	print "$snpID\tmut:$mut\tpositiveCount:$positiveCount\tnegativeCount:$negativeCount\n";
		if($positiveCount >= 18){
			$gainMuts{$snpID}{$mut} = "";
			if($essentialPositiveCount >= 1){
				$essentialGainMuts{$snpID}{$mut} = $essentialPositiveCount;
			}
		}else{
				$nonEssentialGainMuts{$snpID}{$mut} = 0 if(($positiveCount >= 3)&&($essentialPositiveCount == 0));
		}
		if($negativeCount >= 18){
			$lossMuts{$snpID}{$mut} = "";
			if($essentialNegativeCount >= 1){
				$essentialLossMuts{$snpID}{$mut} = $essentialNegativeCount;
			}
		}else{
				$nonEssentialLossMuts{$snpID}{$mut} = 0 if(($negativeCount >= 3)&&($essentialNegativeCount == 0));
		}
	}
}


#}
my $out1 = $file.".essentialGainedMuts.list";
my $out2 = $file.".essentialLostMuts.list";
open(OUT1,">$out1");
open(OUT2,">$out2");
print OUT1 "#SNP\tmutation\tessentialGainWindowNumber\n";
foreach my $snp(keys(%essentialGainMuts)){
	my @emuts = keys(%{$essentialGainMuts{$snp}});
	my @semuts = sort{$essentialGainMuts{$snp}{$b} <=> $essentialGainMuts{$snp}{$a}} keys(%{$essentialGainMuts{$snp}});
	foreach my $mut(@semuts){
		my $essentialGwindowN = $essentialGainMuts{$snp}{$mut};
		next if($essentialGwindowN < 5);
		print OUT1 "$snp\t$mut\t$essentialGwindowN\n";
	}
	
}
close OUT1;

print OUT2 "#SNP\tmutation\tessentiaLosslWindowNumber\n";
foreach my $snp(keys(%essentialLossMuts)){
	my @emuts = keys(%{$essentialLossMuts{$snp}});
	my @semuts = sort{$essentialLossMuts{$snp}{$b} <=> $essentialLossMuts{$snp}{$a}} keys(%{$essentialLossMuts{$snp}});
	foreach my $mut(@semuts){
		my $essentialLwindowN = $essentialLossMuts{$snp}{$mut};
		next if($essentialLwindowN < 5);
		print OUT2 "$snp\t$mut\t$essentialLwindowN\n";
	}

}

close OUT2;



