#!/usr/local/bin/perl

#batch_normalize.pl input_file.txt batch_no.txt output_file.txt

use strict;
my $usage = "$0 input output\n";
my $arg=0;
my $input_file = $ARGV[$arg++] or die $usage; 
#Input file is a log-transformed file that combines several batches of data
#1st column=probeID, then - log-transformed (Cy5-normalized) data for each array.

my $batchno_file = $ARGV[$arg++] or die $usage;
#This file specifies batch number on each line for each array. Each line has one number
#Number in line i specifies the batch number (1,2,3...) for array i.

my $output_file = $ARGV[$arg++] or die $usage;

my @batchno;
my $Nbatch = 0;
open(INFO, $batchno_file) or die $!;
while(my $line = <INFO>){
	$line =~ s/\n$//;
	push(@batchno,$line);
	if($Nbatch<$line){ $Nbatch=$line; }
}
close INFO;

my $log10 = log(10);
open(INFO, $input_file) or die $!;
open (OUTPUT, ">$output_file") or die $!;
my $line=<INFO>;
print OUTPUT $line;
my $count;
while(my $line = <INFO>){
	my ($probeID,@data)=split(/\t/, $line);
	my @batches;
	my @dataNorm=();
	for(my $i=0; $i<@data; $i++){
		my $x = $data[$i];
		push(@dataNorm, $x);
		my $ib = $batchno[$i];
		if(!$ib){ print "ERROR: batch not specified for array #$i\n"; exit(0); }
		push(@{$batches[$ib]},$x);
	}
	my $median = median(\@dataNorm);
	my @med=();
	for(my $ib=1; $ib<=$Nbatch; ++$ib){
		$med[$ib] = median($batches[$ib]);
	}
	for(my $i=0; $i<@data; $i++){
		my $x = $dataNorm[$i];
		my $ib = $batchno[$i];
		$x = $x - $med[$ib] + $median;
		$dataNorm[$i] = int(10000*$x+0.5)/10000;
	}
	print OUTPUT "$probeID\t".join("\t",@dataNorm)."\n";
	$count++;
}
close OUTPUT;
close INFO;
exit(0);

#**************************
sub   median
#**************************
{
my $xr = shift;
my $n = shift;

if(!$xr){ return 0; }
if(!$n){ $n = @$xr; }
if(!$n){ return 0; }
my $median = 0;
my @sorted = sort {$a<=>$b} @$xr;
my $i = $n/2;
if($i > int($i)){
	$i = int($i);
	$median = $sorted[$i];
}else{
	$median = ($sorted[$i] + $sorted[$i-1])/2;
}
return $median;
}

