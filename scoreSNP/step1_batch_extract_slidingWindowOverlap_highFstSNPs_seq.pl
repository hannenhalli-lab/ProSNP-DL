#!/usr/bin/perl -w
use List::Util qw(min max shuffle);
use strict;

### first partition the SNP file 1Kgenomes_phase3_commonAFR_or_commonEUR_snps.commonChrs.top5percentFst.bed to X bins (100..1000)
### one example of a bin file called "1Kgenomes_phase3_commonAFR_or_commonEUR_snps.commonChrs.top5percentFst.bed.bin_1" is listed here.
### In each bin file, extract 20 1KB sliding windows of fasta sequence for each SNP, using nibFrag


my $nfDir = "/data/lis11/tools/nibFrag";
my $NIB_DIR = "/data/lis11/tools/nibFrag/hg19/NIB";


my $binID = shift;
my $ef = "1Kgenomes_phase3_commonAFR_or_commonEUR_snps.commonChrs.top5percentFst.bed.bin_".$binID;
my @elines = `cat $ef`;
$ef=~/^(\S+)\.bed/ or die("$ef ??\n");
my $prefix = $1;
my $out = $prefix.".slidingWindowOverlapSNP.DLinput.bin_".$binID;
open(OUT,">$out");

for my $i(0..$#elines){

	my $line = $elines[$i];
	chomp $line;
	my @spline = split(/\t/, $line);
	my $chrom = $spline[0];
	my $snpCoord = $spline[2];
	my $aleStr = $spline[3];
	my @ales = split(/\-/, $aleStr);
	my $rsID = shift @ales;
	shift @ales;### chrom
	shift @ales;### snp position
	my $refAle = shift @ales;
	my @refNTs = split(//, $refAle);
	next if(scalar(@refNTs)>1);	
	my $leftBorder = $snpCoord - 1000;### 0-based
	my $i = $leftBorder;
	my $from;
	my $to;
	my $windowCenter;
	my $windowDis;
	for($i = $leftBorder; $i <= $snpCoord; $i = $i + 20){
		$from = $i - 1;
		$to = $from + 1 + 1000;
		$windowCenter = $from + 1 + 500;
		$windowDis = $windowCenter - $snpCoord;
		my $seq = `$nfDir\/nibFrag -upper $NIB_DIR\/$chrom.nib $from $to + stdout | grep -v "\>"`;
		$seq =~ s/\n//g;
		my $snpRelaPos = $snpCoord - ($from + 1);
		if($snpRelaPos < 0){
			die("snpCoord:$snpCoord\ti:$i\tfrom:$from\n");
		}
		my $ont = substr($seq, $snpRelaPos, 1);
		if($ont ne $refAle){
			die("nt in the seq: $ont is not the same as refAle: $refAle\n");
		}
		my $windowName = "windowDis:".$windowDis;
		print OUT "$seq\t$windowName\t$rsID\-refAle:$refAle\n";
		my @nts = split(//, $seq);
		my $seqL = scalar(@nts);
		foreach my $mutAleStr(@ales){
			my @mutAles = split(/\,/, $mutAleStr);
			foreach my $mutAle(@mutAles){
				my @mutAleNTs = split(//, $mutAle);
				next if(scalar(@mutAleNTs) > 1);
				my @mutSeqNTs = @nts;
				splice(@mutSeqNTs, $snpRelaPos, 1, $mutAle);
				my $mutSeqL = scalar(@mutSeqNTs);
				die("mutSeqL:$mutSeqL\toriginalSeqL:$seqL\n") if($mutSeqL != $seqL);
				my $mutSeq = join("", @mutSeqNTs);
				print OUT "$mutSeq\t$windowName\t$rsID\-mutAle:$mutAle\n";
			}
		}
	}
	if($i < $snpCoord){
		$from = $snpCoord - 1;
		$to = $snpCoord + 1000;
		$windowCenter = $snpCoord + 500;
		$windowDis = $windowCenter - $snpCoord;
		my $seq = `$nfDir\/nibFrag -upper $NIB_DIR\/$chrom.nib $from $to + stdout | grep -v "\>"`;
		$seq =~ s/\n//g;
		my $snpRelaPos = $snpCoord - ($from + 1);
		my $ont = substr($seq, $snpRelaPos, 1);
		if($ont ne $refAle){
			die("nt in the seq: $ont is not the same as refAle: $refAle\n");
		}
		my $windowName = "windowDis:".$windowDis;
		print OUT "$seq\t$windowName\t$rsID\-refAle:$refAle\n";
		my @nts = split(//, $seq);
		foreach my $mutAleStr(@ales){
			my @mutAles = split(/\,/, $mutAleStr);
			foreach my $mutAle(@mutAles){
				my @mutAleNTs = split(//, $mutAle);
				next if(scalar(@mutAleNTs) > 1);
				my @mutSeqNTs = @nts;
				splice(@mutSeqNTs, $snpRelaPos, 1, $mutAle);
				my $mutSeq = join("", @mutSeqNTs);
				print OUT "$mutSeq\t$windowName\t$rsID\-mutAle:$mutAle\n";
			}
		}
	}
	
}


close OUT;

