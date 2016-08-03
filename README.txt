
Instructions:



1. Download everything to a directory in Linux/Unix


2. Install Eisen's cluster software (http://bonsai.hgc.jp/~mdehoon/software/cluster/software.htm)


3. Make sure R and Perl are pre-installed


4. Run it like this 



% Perl  WC_wapper_SCCA_final.pl positve_target.txt  all_transcription_factor.txt positive_guide_genes.txt  



Since it usually takes 5~20 h to run, it is better to run your ananlysis in background like this 



%nohup  WC_wapper_SCCA_final.pl positve_target.txt  all_transcription_factor.txt  positive_guide_genes.txt  &



Alternatively you can use screen to detach your terminal from server after the job is submitted.
 



We provided two sets of sample data sets for testing.


Codes:

A. WC_wapper_SCCA_final.pl    main Perl script
B. ascca_parser.pl            parser 
C. cluster_parse.pl           parser
D. cal_scca_freq.pl           calculate frequency for each TF discovered

E.  ascca_CV_gamma_command_line.R     main ASCCA code
F.  scca_command_line.R               SCCA code that not used for the package, a gift.  You can use it for standalone analysis
G.  func.R                            required by ASCCA and SCCA


  
A set of test data


Salt_stress_allTF1640.txt
Salt_stress_positive_TF13_guide.txt
Salt_stress_target_Gene157_bait.txt


How to run?
%nohup  WC_wapper_SCCA_final.pl  Salt_stress_target_Gene157_bait.txt  Salt_stress_allTF1640.txt  Salt_stress_positive_TF13_guide.txt  &


Check the output next day.
