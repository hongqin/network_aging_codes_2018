#!/usr/bin/perl
use strict;

use Cwd;
my $dir = getcwd;

opendir(DIR, $dir );
my @tokens = split( /\//, $dir );
my $strain = pop (@tokens);

print "current dir is: [$dir]\n";
print "strain is [$strain]\n";
my $outfile = "$strain.rls.tab";

my @lines = readdir(DIR);
print "\@lines=[@lines]\n";
my @files = grep ( /\d{6}\.$strain\.rls\.tab/, @lines );
print "\@files=[@files]\n";

my @new_fls = ();

my $merged_file  = "$strain.merged";
open(OUT, ">$merged_file"); close(OUT);

foreach my $fl (@files) { 
 my @tokens = split ( /\./, $fl );
 print "\@tokens = [@tokens]\n";
 my $new_fl = "$strain.$tokens[0]";
 print "old [$fl]  new[$new_fl]\n"; 
 push ( @new_fls, $new_fl );
 system( "cp $fl $new_fl" );
 system( "cat $new_fl >> $merged_file");
}

close(DIR);

open( TMP, ">_tmp");
foreach my $fl ( @new_fls ) { print TMP "$fl\n"; }
print TMP "$merged_file\n";
close( TMP );

system( "gtrls01 -i _tmp -o $outfile" );
system( "ln -sf $outfile rls.tab" );

system( "rm _tmp" );
foreach my $fl ( @new_fls ) { system( "rm $fl " ); }

system( "parse_tab_to_winmodest.00.pl -i $outfile");
