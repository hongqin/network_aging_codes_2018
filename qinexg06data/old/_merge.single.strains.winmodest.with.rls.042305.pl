#!/usr/bin/perl

use strict;	use warnings; use Cwd;
use lib '/home/hqin/lib/perl/';     use BeginPerlBioinfo;   use Util;

#note 101S is X101S in R-output

my $debug = 1;
my ( 	%rls_avg,  	%rls_sd, 
	%I, %I_low,	%I_up,
	%G, %G_low, 	%G_up,
	%buffer,	%buffer2,
	) = ();

###### parse winmodest results
 #expt    Line       Start(a)       Start(b)    a        b        LCI(a)    LCI(b)     UCI(a)     UCI(b)    Like      Inform
open( IN, "<natural.winmodest.out"); 
my @lines = <IN>;
shift @lines;
foreach my $line (@lines) {
 chomp $line; 
 my ($expt, $e2, $e3, $e4, $a, $b, $a_low, $b_low, $a_up, $b_up, @rest) = split ( /[\s\t]+/, $line);
 $buffer{ $expt } = "$a\t$b\t$a_low\t$a_up\t$b_low\t$b_up";
}
close(IN);

if ($debug) { showHashTable( \%buffer ); }


##### parse rls mean and sd values
#my $longline = "101S BY4743 M1-2 M13 M14 M22 M2-8 M32 M34 M5 M8 RM112N YPS128 YPS163";
 #my @strains = split( /\s+/, $longline );

open ( IN, "<042305.average.rls.tab"); my @lines2 = <IN>; close (IN);
shift @lines2;
foreach my $line (@lines2 ) {  
  chomp $line;
  my ($expt, $avg, $avg_se, $median, $n, @rest) = split (/[\s\t]+/, $line);
  if ($expt =~ m/M1\.2/) { $expt =~ s /M1\.2/M1\-2/o; }
  if ($expt =~ m/M2\.8/) { $expt =~ s /M2\.8/M2\-8/o; }
  
  $buffer2{ $expt } = "$avg\t$avg_se\t$median\t$n";
}

if ($debug) { showHashTable( \%buffer2 ); 
  print "\n";
#  print "expts from winmdoest:[". (sort (keys %buffer)  ). "]\n\n";
#  print "expts from survfit  :[". (sort (keys %buffer2) )	. "]\n\n";
}


#### now generate the merged tab file for R
open (OUT, ">042405.rls.I.G.summary.tab");
print OUT "expt\ta\tb\ta_low\ta_up\tb_low\tb_up";  #header from winmodest
print OUT "\tavg\tavg_se\tmedian\tn\n"; #header from survfit

my $i=0;
foreach my $expt ( sort(keys %buffer) ) {
  print OUT "$expt\t$buffer{ $expt }";
  if ( exists $buffer2{$expt} ) {
      print OUT "\t$buffer2{ $expt }\n";
  } else { 
      print OUT "\tERROR\n";
  }
  
  $i++;
  print "$i ... $expt\n"; 
}
close(OUT);

exit; 
