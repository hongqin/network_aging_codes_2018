#!/usr/bin/perl
use strict; use Cwd; use warnings;

#my $longline = "101S M1-2 M13 M14 M22 M2-8 M32 M34 M5 M8 RM112N SGU57 W303 BY4743 YPS128 YPS163 SK1";
my $longline = "101S M5 M8  W303 BY4743 YPS163 SK1";

my @strains = split( /\s+/, $longline );

my $rootdir = getcwd;

foreach my $strain ( @strains ) {
 my $dir = "$rootdir/$strain";
 my $file = "$strain.merged";
 print "\$file = [$file]\n";

 chdir $dir ; 
 my $wd = getcwd;
 print "now in ..$wd/\n";

 system( "cp $file rls.tab" );
 system( "cp ../bootstrap.00.r . ");
 system( "cp ../r.sh ." );
 system( "sh r.sh bootstrap.00.r");

 system( "parse_tab_to_winmodest.00.pl -i boot.tab -o $strain.boot.dat" );
 
 #system ( " rm rls.tab ");
}#strain


