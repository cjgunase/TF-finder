#!/usr/bin/perl
use strict;
use Getopt::Std;
use vars qw($opt_i $opt_h);
getopts('hi:');

my (@fields, @gnams, $num, $prevf1, $prevf2, $prevf3, $presg1, $presf1, $presf2, $presf3, $avg_coef, $abs_coef, $coef, $n);

$n = 1;
$num = 1;
$coef = 0;

open(Ref,"$opt_i") || die "Can't open $opt_i";

# open(REF,"wheat_gid_oligo_AG_wheat_chl_TAGI.100604_P70.txt");
    while (<Ref>) {
       chomp;
       @fields = split(/\t/, $_);

       if ($n == 1) {
	    $prevf1 = $fields[0];
            $prevf2 = $fields[1];
            $num = 1;
            $coef = $fields[2];
            $abs_coef = abs($fields[2]);
	    $n =  0;
       } else  {
            $presf1 = $fields[0];
            $presf2 = $fields[1];
            if (($presf1 eq $prevf1) && ($presf2 eq $prevf2) ){
                  $num++; 
                  $coef = $coef + $fields[2];
                  $abs_coef = $abs_coef + abs($fields[2]);
#                  print "num-sum:$num\t$coef\n";
            } else {
		 my $avg_coef= $coef/$num;
                 my $abs_avg = $abs_coef/$num;
                 print "$prevf1\t$prevf2\t$num\t$avg_coef\t$abs_coef\n";
		 $num  =  1 ;
                 $coef = $fields[2];
                 $abs_coef = abs($fields[2]);
	         $prevf1 = $presf1;
	         $prevf2 = $presf2;
	     }
	}
   }
my $avg_coef= $coef/$num;
my $abs_avg = $abs_coef/$num;
print "$prevf1\t$prevf2\t$num\t$avg_coef\t$abs_coef\n";
close(Ref);
