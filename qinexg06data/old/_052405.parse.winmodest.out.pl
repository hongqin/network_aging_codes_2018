#!/usr/bin/perl
use strict; use Cwd; use warnings;

#my $longline = "101S BY4716 BY4742 BY4743 M1-2 M13 M14 M22 M2-8 M32 M34 M5 M8 PSY316  
#RM112N RXB SGU57 sir2Da sir2Dalpha sir2DSIR2 SK1 srp1ts W303 YPS128 YPS163"; 

my $longline = "101S M5 SK1 YPS163";

my @strains = split( /[\s\n]+/, $longline );

my $rootdir = getcwd;

foreach my $strain ( @strains ) {
 my $dir = "$rootdir/$strain";

 chdir $dir ; 
 my $wd = getcwd;
 print "now in ..$wd/\n";

 my $infile = "$strain.out";
 print "\$infile = [$infile]\n";

 open (IN, "<$infile"); my @lines=<IN>; close (IN); 
 shift @lines;
 shift @lines;
 shift @lines;
 print "\@lines=[$lines[0]]\n";

 my $outfile = "$strain.out.tab";
 open (OUT, ">$outfile"); print OUT @lines; close (OUT);

 system ("cp $outfile wm.out.tab" );

}#strain


