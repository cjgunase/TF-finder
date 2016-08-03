#! /usr/bin/perl -w

use Getopt::Std;

my(@fields, @columns, $prevClu, $str, $row,  $str2, $len, $num);
$row=1;

%opt=();
getopts("hk:",\%opt);
 
open(PAR,"$opt{k}") or die "Cannot open the database file: $!";     

while (<PAR>) {

    next if ((/GROUP/) || (/^\s*$/));
    # chomp();
     s/\r\n?/\n/;  # remove return
     chomp();
     @fields = split(/\t/, $_);
     if  ($row == 1) {
	$str = $fields[0];
        $prevClu=$fields[1];
        $row++; 
     } else {
         if ($fields[1] eq $prevClu) {        
             $str = $str . "\t" . $fields[0];
          }  else {
             print "$prevClu\t$str\n";         
#             $ghash{$prevClu} = $str;
             $prevClu = $fields[1];
             $str = $fields[0];
	  }
    }
}
print "$prevClu\t$str\n";
close(PAR);
