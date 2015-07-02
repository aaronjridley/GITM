#!/usr/bin/perl -w
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf
#
# Script: replace_amrinit.pl
#
# Updates "PARAM" file defined in the optional input parameter (e.g. PARAM.in) to correspond to the current version of SWMF
# E.g., at this moment replaces "#AMRINIT" section in "PARAM" file by #USERINPUTBEGIN 
# and produces two files:
#   if incoming file does not have a version/distribution name in it (e.g. PARAM.in or PARAM.in_init): 
#      initial {PARAM} file is saved as {PARAM}_old_version, e.g. PARAM.in_old_version
#
#   otherwise (if the file name contains a version/distribution component, e.g. PARAM.in_v8.01_init)
#      initial {PARAM} file is left unchanged
#      updated file is saved as PARAM.in_v{DISTRIBUTION_NAME}_...
#
# Three input parameters (optional, not positional):
#     1) Name of PARAM file - if not specified,then PARAM.in will be used
#     2) Grid resolution - if not specified, will use default value:0.03125
#     3) Distribution version - if not specified than use current date YYYYMMDD
#    
$file_to_be_updated='';
$new_distribution = '';
$gridresolution='';
   
for ($argnum = 0; $argnum <= $#ARGV; $argnum++) {
   if ($ARGV[$argnum] =~ m/PARAM/){
       $file_to_be_updated=$ARGV[$argnum];
   } elsif (($ARGV[$argnum] =~ m/200/)||($ARGV[$argnum] =~ m/8.01/)) {
       $new_distribution=$ARGV[$argnum];
       $new_distribution = "_".$new_distribution."_";
   } else {
       $gridresolution=$ARGV[$argnum];
   }
}

# Default cases:
if ($file_to_be_updated eq '' ){
     $file_to_be_updated="PARAM.in";  # set default
     print "Name of file to be updated is set to default value: $file_to_be_updated\n";
} else {
     print "File to be updated: $file_to_be_updated\n";
}
if ($gridresolution eq '' ){
     $gridresolution=0.03125;  # set default
     print "Grid resolution set to default value: $gridresolution\n";
} else {
     print "Grid resolution requested: $gridresolution\n";
}
if ($new_distribution eq '' ){
#    Distribution date is calculated based on the current local date
     my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
     
     $new_distribution = sprintf ("%4d", $year+1900).sprintf ("%02d", $mon+1).sprintf ("%02d",$mday);
     $new_distribution = "_".$new_distribution."_";
     print "Distribution set to default value: $new_distribution\n";
} else {
     print "Distribution selected: $new_distribution\n";
}


if ( ! -e "$file_to_be_updated" ){
     print "Please run this script in the directory containing $file_to_be_updated file\n";
     print  "exiting script\n";    exit;
}


# Set the names for the updated file and the original file:
if ($file_to_be_updated =~ m/_v/){
#    input file name contains "Version/Distribution" string 
     $old_file_name=$file_to_be_updated;
     $new_file_name = $file_to_be_updated;
     $file_to_be_updated =~ m/(_v.*?_)/;
     $old_version = $1;
     print "Old version: $old_version\n";
     $new_file_name =~ s/$old_version/$new_distribution/;

} else {
     $old_file_name=$file_to_be_updated."_old_version"; 
     $new_file_name=$file_to_be_updated;
}
print "Old file name: $old_file_name\n";
print "New file name: $new_file_name\n";

# Determine which line contains AMRINIT identifier:
$line_with_amrinit=`grep -n "#AMRINIT" $file_to_be_updated |head -1 |cut -f1 -d:`;

if ( $line_with_amrinit lt 0){
     print "$file_to_be_updated file does not contain #AMRINIT block\n";
     print  "exiting script\n"; exit;
}

$line_with_amrinit= $line_with_amrinit - 1;
$line_with_grid_name= $line_with_amrinit+1;
$line_with_grid_level= $line_with_amrinit+2;

$first_line_after_amrinit_block=$line_with_amrinit + 3;
$last_line_in_file=`wc -l $file_to_be_updated | awk '{print \$1}'`;
$final_lines = $last_line_in_file - $first_line_after_amrinit_block;

open (PARAMFILE,"<$file_to_be_updated");
@logmes= <PARAMFILE>;
close(PARAMFILE);
@PARAMFILE=@logmes;

chomp($grid_name = `echo "$PARAMFILE[$line_with_grid_name]" | awk '{print \$1}'`);
chomp($grid_name = $grid_name);
chomp($grid_level = $PARAMFILE[$line_with_grid_level]);

$replacement_lines="#USERINPUTBEGIN

#CCMCGRID
$grid_name    InitialRefinementType (character)

#USERINPUTEND

#GRIDLEVEL
$grid_level
init

#GRIDRESOLUTION
$gridresolution
user
";
if ($file_to_be_updated ne $old_file_name) {
   system("cp -p $file_to_be_updated $old_file_name");
}
system("head -$line_with_amrinit $old_file_name > $new_file_name");
open (PARAMFILE,">>$new_file_name");
print PARAMFILE $replacement_lines;
close PARAMFILE;

system("tail -$final_lines $old_file_name >> $new_file_name");

exit(0);
