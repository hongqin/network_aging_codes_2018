#!/usr/bin/perl
use strict; use Cwd; use warnings;

my $longline = "101S M1-2 M13 M14 M22 M2-8 M32 M34 M5 M8 RM112N SGU57 W303 BY4743 YPS128 YPS163";
my @strains = split( /\s+/, $longline );

my $rootdir = getcwd;
my $outfile = "052305.strain.a.sd.b.sd.tab";

open ( OUT, ">$outfile" );
print OUT "strain\ta\ta.sd\tb\tb.sd\n";

foreach my $strain ( @strains ) {
 my $dir = "$rootdir/$strain";

 chdir $dir ; 
 my $wd = getcwd;
 print "now in ..$wd/\n";

 system( " cp ../_052405b.get.mean.sd.from.WMout.r .");
 system( "sh r.sh _052405b.get.mean.sd.from.WMout.r ");

 my $infile = "summary.wm.out.tab";

 open (IN, "<$infile"); my @lines=<IN>; close (IN); 
 print OUT "$strain\t$lines[1]";

}#strain

close ( OUT);
exit;

