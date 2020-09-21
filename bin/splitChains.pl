#!/usr/bin/perl -w
use strict;

my $usage = "It splits the protein either in files as many as the number of chains or it keeps only a single chain\n./splitChains.pl -in <PDB FILE> <-hetero/-homo> -chain (if homo)\n\n";

if($#ARGV < 0){ die $usage; }

my $infile = "";
my $chain="A";
my $homo = 1;

while(my $args = shift @ARGV){
    if($args =~ /^-in$/){ $infile = shift @ARGV; next; }
    if($args =~ /^-hetero$/){ $homo = 0; next;}
    if($args =~/^-homo$/){
	if($homo == 0){ die "Cannot give both -homo and -hetero\n"; }
	$homo = 1;
	next;
    }
    if($args =~ /^-chain$/){ $chain = shift @ARGV; next; }
    die "Argument $args is not valid. $usage\n";
}

open(IN, $infile) or die "Couldn't open $infile for input\n";

my $definedCh="";
my @pream = ();
while(defined(my $ln = <IN>)){
    #print $ln, "\n";
    chomp($ln);
    if($ln =~ /^\s*$/){next;}
    if($ln !~/^\s*ATOM/){
#	push @pream, $ln."\n";
	print $ln, "\n";
	next;
    }
    my @l = split(//, $ln);
    my $ch = $l[21];
    if($ch ne $definedCh && (($homo == 1 && $ch eq $chain)||($homo == 0 ) ) ){
	if($definedCh ne ""){ close OUT; }
	my $outputFile = $infile;
	$outputFile =~ s/\.pdb$//;
	$outputFile .= $ch.".pdb";
	open(OUT, ">$outputFile") or die "Couldn't open $outputFile for output\n";
	$definedCh = $ch;
    }

    if($ch eq $definedCh){
	print OUT $ln, "\n";
    }
}










