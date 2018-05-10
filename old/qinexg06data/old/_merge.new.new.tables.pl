#!/usr/bin/perl

use strict;	use warnings; use Cwd;
use lib '/home/hqin/lib/perl/';     use BeginPerlBioinfo;   use Util;

my $debug = 1;

my $longline = "101S BY4743 M1-2 M13 M14 M22 M2-8 M32 M34 M5 M8 RM112N YPS128 YPS163";
my @strains = split( /\s+/, $longline );
my @out = ();
my $header = ();

foreach my $stain ( @strains ) {
 my $in = "$stain/new.new.$stain.out";
 open (IN, "<$in" ); my @lines = <IN>; close (IN);
 $header = $lines[0];
 shift @lines;;
 push ( @out, @lines);
}

open (OUT, ">042305.new.new.all.for.predictG.by.max.out");
print OUT $header;
foreach my $line (@out){
 print OUT $line;
}
close(OUT);

exit; 
