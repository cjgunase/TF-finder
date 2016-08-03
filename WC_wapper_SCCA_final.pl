#! c:/Perl/perl.exe -w

use strict;
use Sort::Fields;

my (@fields, @carray, @cluster, $temp, @temp1, @temp2, @temp3, $i, $num, $gnum, $str, $len, $clust, $item, @ascca, @sorted_ascca, $fname, $mod1, $mod2, $max, $min, $TF_num);
my ($N_ptf, $n, $m, $N_tf, $N_ptar);
my %ghash=();
my %phash=();
my $output_cluster="output_cluster.txt";
my $tmp_cluster = "tmp_cluster.txt";
my $output_ASCCA_freq="output_ASCCA_freq.txt";   # file with final result: the frequency of each gene being caught
my $output_ASCCA="output_ASCCA.txt";             # file contain each gene caught, the 3rd col is the correlation cofficient 

my $TGfile = shift;    # this is target input file
my $TFfile =shift;      # this is transcription factor file
my $PTF=shift;         # gene of positive TF
my $EF= 3;            # Enrichment factor, you can reset this

unless($TGfile and  $TFfile and $PTF){
    print "Usage:                                                                                                                              
    $0 target_file TF_file postiveTF_list_file\n";
    exit;
}

# push expression data of target (positive) genes into hash
$N_ptar = 0;
open(TAR,"$TGfile") or die "Cannot open the database file: $!";
while (<TAR>) {
    chomp();
    unless (/^\s*$/) {
       $N_ptar = $N_ptar + 1; 
       @fields = split(/\t/, $_);    
       $num = @fields;
       $str=$fields[1];
       for ($i =2; $i < $num; $i++) {
           $str = $str . "\t". $fields[$i];
       }
       $ghash{$fields[0]} = $str;
       undef($str);
   }
}
close(TAR);

$N_ptf=0;
open(PTF,"$PTF") or die "Cannot open the database file: $!";
while (<PTF>) {
    chomp();   
    unless (/^\s*$/) {                                                                          
       $N_ptf=$N_ptf + 1 ;
       @fields = split(/\t/, $_); 
       $phash{$fields[0]} = $fields[0];
   }
}
close(PTF);

$N_tf=0;
open(TF,"$TFfile") or die "Cannot open the database file: $!";
while (<TF>) {
    unless (/^\s*$/) {
       chomp();
       $N_tf=$N_tf + 1 ;
   }
}
close(TF);


open(OUT,">$output_cluster") or die "Cannot open the database file: $!";

# perform mutiple cluster analysis

$mod1 = $N_ptar%4;  # A cluster on average allow to have minmal 4 genes
$mod2 =$N_ptar%20;  # A cluster on average allow to have maximal 20 genes

if ($mod1 == 0) {
     $max=($N_ptar-$mod1)/4;
} else {
    $max=($N_ptar-$mod1)/4+1;
}
if ($mod2 == 0){
    $min=($N_ptar-$mod2)/20;
} else {
    $min=($N_ptar-$mod2)/20+1;
}


# repeatedly cluster the target genes with k-mean and save all results

for ($i=$min;$i <=$max; $i++){
     qx(cluster -f $TGfile -g 7 -k $i);
     $fname = $TGfile;
     $fname =~ s/.txt//;  
     $fname = $fname  . "_K_G" . $i .".cdt" ; 
     qx (rm $fname);
     $fname =~ s/.cdt//;
     $fname = $fname .".kgg" ;
     @temp1= `perl cluster_parse.pl -k $fname`;
     qx(rm $fname);
     foreach $item (@temp1) {
         push (@carray, $item);
     }         
 }
print OUT  @carray;
close(OUT);

open(SCCA,">$output_ASCCA") or die "Cannot open the database file: $!";
open(FREQ,">$output_ASCCA_freq") or die "Cannot open the database file: $!";

foreach $clust (@carray) {
     open(INCL,">$tmp_cluster") or die "Cannot open the database file: $!";
     chomp($clust);
     @cluster = split(/\t/, $clust);
     $len=@cluster;
     if ($len > 3) {   #only cluster with more than 2 genes are used for ASCCA analysis
         for ($i=1; $i < $len; $i++) {
            print INCL "$cluster[$i]\t$ghash{$cluster[$i]}\n";
         }
     }
      
     close(INCL); 

     qx(R --no-restore --no-save --args "$TFfile" "$tmp_cluster" < ascca_CV_gamma_command_line.R);
     @temp2=qx(perl ascca_parser.pl -i out_ascca_CV_gamma.txt);
     @temp3 = grep /^TF/, @temp2;     # number of TF hooked
     $n=@temp3; 
     $m=0;
     foreach $item (@temp3) {
         @fields = split(/\t/, $item);
         if ($phash{$fields[1]}) {
	     $m=$m + 1;   # number of positive TF hooked 
         }
     }

     $TF_num = $EF * ($N_ptf/$N_tf) * $n;  # enriched theory number of positive TFs in hooked TF list
     if ($TF_num < 1) {$TF_num=1;}  
     if ($m >= $TF_num) {
         foreach $item (@temp2) {
            push (@ascca, $item);
         }
     }
 }

     @sorted_ascca= fieldsort [1, 2, '3n'], @ascca;
     print SCCA @sorted_ascca;
     
     @temp3 = `perl cal_scca_freq.pl -i $output_ASCCA`;
     print  FREQ @temp3;
     qx(rm out_ascca_CV_gamma.txt);    
close (FREQ);
close(SCCA);

  

