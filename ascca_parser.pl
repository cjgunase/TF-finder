#! /usr/bin/perl -w

use Getopt::Std;

my(@fields, @columns, $prevClu, $str, $str2, $i, $len, $num);

my %ghash = (); # force %ghash empty

#"31" -0.267898777020056 "AT1G01290" -0.132906335072070 "AT5G40890" -0.0501630091762539 "AT5G62730"
#"1" 0.158839605819669 "AT1G03230" 0.191547000684979 "AT5G43060" -0.199247980037148 "AT3G28220"

%opt=();
getopts("hi:",\%opt);
 
open(Ref,"$opt{i}") or die "Cannot open the database file: $!";     

while (<Ref>) {
    s/\"//g;

    if (/Transcription/) {
       $type = "TF";
   } elsif (/Target Genes/) {
       $type = "Target";
   }

    if (/^\d+\W-*\d+\.\d+\W\w+\W-*\d+\.\d+\W\w+\W(-*\d+\.\d+)\W(\w+)/){
        print "$type\t$2\t$1\n";
    }
}
close(Ref);
