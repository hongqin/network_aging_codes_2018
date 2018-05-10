#!/usr/bin/perl

my $debug = 3;

BEGIN { unshift(@INC,"/home/hqin/lib/perl/");   }

# output format
# experiment mean se median n 

use strict; use warnings; use Getopt::Long;
use Util; use Cwd; 

#my $longline = "101S BY4743 M1-2 M13 M14 M22 M2-8 M32 M34 M5 M8 RM112N YPS128 YPS163 M5-1abcd ";
#my @strains = split( /\s+/, $longline );

open (IN, "<042305.survfit.out.merged.txt"); my @lines=<IN>; close(IN);

my %reports = ();

my ( @starts ) = ();
my ( @pos_results ) = ();

# find the starting positions
for( my $i = 0; $i<= $#lines; $i++) {
 if ( $lines[$i] =~ /^> ######/ ) { push (@starts, $i); }
 if ( $lines[$i] =~ /^\s+n\s+events\s+rmean/ ) { push (@pos_results, $i+1); }
 }
# now @starts stores the position of each record
print "\@start=[@starts]\n";
print "\@pos_starts=[@pos_results]\n";

  # experiment mean se median n 
# get the result section
for ( my $j=0; $j<= $#starts; $j++ ) {
 # get experiment name 
 my @els = split ( /\"/, $lines[ $starts[$j]  + 4]);
 my $name = $els[1];  # print $name."\n";
 $name =~ s/^X//o;
 
 # results 
 my $result = $lines[ $pos_results[$j] ] ;
 $result =~ s/^\s+//o;
 chomp $result;
 my ($n, $events, $rmean, $se, $median, @rest ) = split (/\s+/, $result);
 
 $reports{ $name } = "$rmean\t$se\t$median\t$n";
}

if ($debug) { showHashTable( \%reports ); }

open ( OUT, ">042305.average.rls.tab");
 print OUT "expt\tavg\tavg.se\tmedian\tn\n";
 foreach my $name ( sort( keys %reports ) )  {
   print OUT "$name\t$reports{$name}\n";
 }
close (OUT);

exit;

