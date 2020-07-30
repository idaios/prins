#!/usr/bin/perl -w

use strict;

my $usage = "
It annotates a protein 3D file based on the scoring of each aminoacid
proteing3Dannotation.pl score=<file> pdb=<file>

";

if($#ARGV < 0)
{
    die $usage; 
}

my $scoreFile = "";
my $pdbFile = "";
my $schemeString = "5,5,10";
my $thresholdString = "0,0";
my $tf = 0;
my $scoreIndex = 5;
my $avScoreIndex = 6;
my $maxScoreIndex = 7;
my $chainIndex = 9;
my $posIndex = 2;
my $exclude = "";

my $temperature = 0;

while(my $args = shift @ARGV)
{
    if($args =~ /^score=(.*)/i){ $scoreFile = $1; next; }
    elsif( $args =~ /^scoreIndex=(.*)/i){ $scoreIndex = $1; next; }
    elsif($args =~ /^pdb=(.*)/i){ $pdbFile = $1; next; }
    elsif($args =~/^scheme=(.*)/){ $schemeString = $1; next;}
    elsif($args =~ /^threshold=(.*)/i){ $thresholdString = $1; next; }
    elsif($args =~ /^exclude=(.*)/){ $exclude= $1; next; }
    elsif($args =~ /^-temp/i){ $temperature = 1; next; }
    else{die "Argument $args is invalid\n$usage\n"; }
}

my @threshold = ();
if($thresholdString ne "0,0")
{
    @threshold = split(/,/, $thresholdString);
    push(@threshold, -9999999); # a hypothetical minimum value
    $tf = 1;
}


my @scheme = split(/,/, $schemeString);
my $remaining = 0;

my $sum = 0.;
for(my $i = 0; $i < @scheme; ++$i)
{
    $sum += $scheme[$i];
}

if($sum > 100)
{
    print STDERR "SUM cannot be greater than 100 in schemeString\n";
}
else
{
    $remaining = 100 - $sum;
}

if($remaining > 0)
{
    push @scheme, $remaining;
}

for(my $i = 1; $i < @scheme; ++$i)
{
    $scheme[$i] += $scheme[$i-1];
}

open(PDB, $pdbFile) or die "Couldn't open PDB file for input\n";
open(SCORE, $scoreFile) or die "Coudln't open SCORE file for input\n";



my %at = ();
# my %index = ();
my @w = ();
my %sc = ();
my %chain = ();
my $identifier = "";

while(defined(my $ln = <SCORE> ) )
{
    if($ln =~ /^\s*$/){ next; }
    if($ln =~ /score/i){ next; }
    chomp($ln);

    @w = split(/\s+/, $ln);

    $identifier = $w[ $posIndex ].$w[ $chainIndex ];

    $at{ $identifier } = $w[0];
    $sc{ $identifier } = $w[ $scoreIndex ];
    $chain{ $identifier } = $w[ $chainIndex ];

}

my @sortedKeys = sort{ $sc{$b} <=> $sc{$a} } keys %sc;
my $l = scalar @sortedKeys;
my %coloredAtom = ();
my @orderedCols = ("S", "P", "O", "N", "C", "H");

my $j = 0;

if($tf == 0)
{

    for(my $i = 0; $i < $l; ++$i)
    {
	
	if( $i < $scheme[$j] * $l / 100 )
	{
	    
	    $coloredAtom{ $sortedKeys[$i] } = $orderedCols[$j];
	    if($j == $#scheme)
	    {
		$coloredAtom{ $sortedKeys[$i] } = "H";
	    }
	}
	else
	{
	    $j++;
	    $i--;
	    next;
	}
	
    }
}
elsif($tf == 1)
{
    
    for(my $i = 0; $i < $l; ++$i)
    {
	#print STDERR $j, "\t", $threshold[$j], "\n";
	
	if( $sc{ $sortedKeys[$i] } > $threshold[$j] ) 
	{
	    
	    $coloredAtom{ $sortedKeys[$i] } = $orderedCols[$j];

	    if($j == $#threshold)
	    {
		$coloredAtom{ $sortedKeys[$i] } = "H";
	    }
	}
	else
	{
	    $j++;
	    $i--;
	    next;
	}
	
    }
}

 # foreach my $k (keys %coloredAtom){
 #     print $k, "\t", $coloredAtom{$k}, "\t$sc{$k}\n";
 # }
my $atomID;
my $chain;
my $position;
my $aa;

my $minTemp = 9999999;

my $maxTemp = -9999999;


while( defined(my $ln = <PDB> ) )
{
    if($ln =~ /^\s*$/){ next; }
    chomp($ln);

    $atomID = substr($ln, 0, 6);
    $atomID =~ s/ //g;

    $chain = substr($ln, 21, 1);
    $chain =~ s/ //g;

    $position = substr($ln, 22, 4);
    $position =~ s/ //g;

    $aa = substr($ln, 17, 4);
    $aa =~ s/ //g;

    if( $ln =~ /^\s*END/ )
    {
	print "END\n";
	last;
    }
	
    
    #print "$atomID\t$chain\t$position\t$aa\n";
    

     if($atomID ne "ATOM")
     {

	
	if(($exclude ne "" ) && ( $ln =~ /$exclude/i) )
	{
	    next;
	}
	print $ln, "\n";

     }
    else{
	
	
	if(($exclude ne "" ) && ( $ln =~ /$exclude/i) )
	{
	    next;
	}

	
     	$identifier = $position.$chain; 

	#print $identifier, "\n";

	my $tempScore = 0.0; 

     	if( defined($at{ $identifier }) && $at{$identifier} ne $aa )
     	{
     	    print STDERR "Error!: aa $at{ $identifier } ( $identifier ) does not correspond to $aa ($identifier)\n";
     	    exit;
    	}



    # 	#my $letter = substr($ln, 77, 1);
    # 	#print $letter, "\t$w[5]\t", $coloredAtom{$w[5]}, "\n";
	##print STDERR $temperature, "\n";
	if(defined ($coloredAtom{ $identifier} ) ){
	    $tempScore = sprintf("%5.2f", $sc{ $identifier } );
	    if($temperature == 0)
	    {
		substr($ln, 77, 1) = $coloredAtom{ $identifier };
	    }
	}
	else
	{
	    $tempScore = "000.00";
	    print STDERR "Warning (UNDEFINED) at line $ln (key: $identifier)\n";
	    next;
	    #substr($ln, 77, 1) = "H";
	}

	substr($ln, 60, 5) = $tempScore;
	
	if($tempScore > $maxTemp)
	{
	    $maxTemp = $tempScore;
	}

	if($tempScore < $minTemp )
	{
	    $minTemp = $tempScore;
	}

     	print $ln, "\n";
    }
    
}

if( $temperature == 1)
{

    print STDERR "minTempScore: $minTemp, maxTempScore: $maxTemp\n";

}
