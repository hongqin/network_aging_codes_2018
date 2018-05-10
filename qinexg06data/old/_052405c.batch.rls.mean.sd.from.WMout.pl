#!/usr/bin/perl
use strict; use Cwd; use warnings;

#my $longline = "101S M1-2 M13 M14 M22 M2-8 M32 M34 M5 M8 RM112N SGU57 W303 BY4743 YPS128 YPS163";
#my $longline = "101S BY4716 BY4742 BY4743 M1-2 M13 M14 M22 M2-8 M32 M34 M5 M8 PSY316 
#RM112N RXB SGU57 sir2Da sir2Dalpha sir2DSIR2 SK1 srp1ts W303 YPS128 YPS163";

my $longline = "101S M1-2 M13 M14 M22 M2-8 M32 M34 M5 M8 SGU57 
RM112N 
YPS128 YPS163
BY4743 
SK1 
W303 
PSY316 
";


my @strains = split( /[\s\n]+/, $longline );

my $rootdir = getcwd;
#my $outfile = "052305.strain.a.sd.b.sd.tab";
my $outfile = "061405.rls.for.publication.tab";

open ( OUT, ">$outfile" );
print OUT "strain\tn\trls\trls.sd\ta\ta.sd\tb\tb.sd\n";

foreach my $strain ( @strains ) {
 my $dir = "$rootdir/$strain";

 chdir $dir ; 
 my $wd = getcwd;
 print "now in ..$wd/\n";

 my $fl = "$strain.merged";
 open (IN, "<$fl"); my @lines2 = <IN>; close (IN);
 my $samplesize = $#lines2 + 1;
 print "***\$samplesize= $samplesize \n";

 system(" cp ../_052405b.get.rls.a.b.with.sd.r .");
 system(" sh r.sh _052405b.get.rls.a.b.with.sd.r ");

 my $infile = "summary.wm.out.tab";

 open (IN, "<$infile"); my @lines=<IN>; close (IN); 
 print OUT "$strain\t$samplesize\t$lines[1]";

}#strain

close ( OUT);
exit;

