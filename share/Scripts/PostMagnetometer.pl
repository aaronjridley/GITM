#!/usr/bin/perl -s
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

my $help = ($h or $H or $help);

use strict;

$help = 1 if not @ARGV; 

&print_help if $help;

my $ERROR = "ERROR in PostMangetometer.pl";

my %stations;  # hash of station names (all capitals)
my %data;      # hash of magnetometer data indexed by station and date
my %read;      # hash to check which file types have been read
my %check;     # hash to check multiple occurance of dates in one type

# list of variables in the output files
my $variablesout  = "year mo dy hr mn sc ms X Y Z sumBn sumBe sumBd";

my $variablesgm   = " dBn dBe dBd facdBn facdBe facdBd";
my $variablesie   = " JhdBn JhdBe JhdBd JpBn JpBe JpBd";
my $variablesdata = " dataBn dataBe dataBd";

# list of original data items in the output files
my $ndataout = 3;     # X Y Z are alway there

# Error collected for all stations
my %rms;

my $file; 
my $type;
foreach $file (sort @ARGV){
    print "reading $file\n";
    open FILE, $file or die "$ERROR: could not open file $file\n";
    if($file =~ /^(GM|IE)/){
	$type = $1;
	&process_swmf_file;
    }else{
	$type = "data";
	&process_data_file;
    }
    # Remember all file types read
    $read{$type} = 1;
    close FILE;
}

&write_output;

if(%rms){
    my @avgrms;
    my $nstation;
    print "-" x 40, "\n";
    print "RMS errors per station:\n";
    print "-" x 40, "\n";
    foreach my $station (sort keys %rms){
	$nstation++;
	my @rms = @{$rms{$station}};
	printf "$station: %8.4f %8.4f %8.4f\n", @rms;

	$avgrms[0] += $rms[0];
	$avgrms[1] += $rms[1];
	$avgrms[2] += $rms[2];
    }
    $avgrms[0] = $avgrms[0]/$nstation;
    $avgrms[1] = $avgrms[1]/$nstation;
    $avgrms[2] = $avgrms[2]/$nstation;

    print "-" x 40, "\n";
    printf "all: %6.2f %6.2f %6.2f\n", @avgrms;
}

exit 0;

##########################################################################
sub process_swmf_file{

# GM header
#       12 magnetometers:  YKC MEA NEW FRN IQA PBQ OTT FRD HRN ABK WNG FUR
#nstep year mo dy hr mn sc msc station X Y Z dBn dBe dBd facdBn facdBe facdBd
#
# IE header
#   12 magnetometers:  YKC MEA NEW FRN IQA PBQ OTT FRD HRN ABK WNG FUR
#nsolve year mo dy hr mn sc msc station X Y Z JhdBn JhdBe JhdBd JpBn JpBe JpBd

    # Read header and extract station names
    my $header = <FILE>;

    $header =~ /magnetometers:/ 
	or die "$ERROR for $file: invalid header=$header\n";

    my $stations = $';
    my @stations = split(' ',$stations);

    my $station;
    foreach $station (@stations){
	$stations{$station} .= $type;
    }

    my $variables = <FILE>;
    if($type eq "GM"){
	die "$ERROR for $file: invalid variables=$variables\n" unless
	    $variables =~ /nstep year mo dy hr mn sc msc station X Y Z$variablesgm/;

	if($variablesout !~ /$variablesgm/){
	    $variablesout .= $variablesgm;
	    $ndataout     += 6;
	}
    }else{
	die "$ERROR for $file: invalid variables=$variables\n" unless
	    $variables =~ /nsolve year mo dy hr mn sc msc station X Y Z$variablesie/;

	if($variablesout !~ /$variablesie/){
	    $variablesout .= $variablesie;
	    $ndataout     += 6;
	}
    }

    # Read data
    while(<FILE>){
	# split line into 18 items (columns)
	my @items = split(' ',$_);

	die "$ERROR in $file at line $.: incorrect number of items in $_"
	    unless $#items == 17;

	# Extract date and station number/name
        $items[6] =~ s/\d$/0/;  # replace second digit of seconds with a 0
	my $date = join(' ',@items[1..6]);
	my $imag = $items[8]-1;
	$station = $stations[$imag];

	# Only record the first occurance
	next if $check{$station}{$date}{$type}++;

	# Store coordinates if it has not been done already
	$data{$station}{$date} .= join(' ',@items[9..11])
	    unless $data{$station}{$date};

	# Store magnetic field
	$data{$station}{$date} .= ' ' . join(' ',@items[12..17]);

	# There is some data for this station
	$read{$station}{$type} = 1;
    }
}
##########################################################################
sub process_data_file{

    # Data files contain the following lines
    #31029  0 0 aae XYZ    12    -4     2     0     0     0
    # The first number has to do with YearMonthDay (" 3" --> 2003 !!)
    # Then hour, minute, station name, Dbn, Dbe, Dbd, and 3 zeros
    # Measurement errors are signaled as -99999
    # Note that the columns are NOT separated by spaces.

    my $nodata = -99999;

    if($variablesout !~ /$variablesdata/){
	$variablesout .= $variablesdata;
	$ndataout += 3;
    }

    while(<FILE>){
	/^(..)(..)(..) (..)(..) (...) XYZ(......)(......)(......)/
	    or die "$ERROR in $file at line $.: could not parse line=$_";

	my $station = uc($6);

	# Skip stations that do not occur in the GM / IE files
	next unless $stations{$station};

	# Extract data
	my $dBn = $7 + 0;
	my $dBe = $8 + 0;
	my $dBd = $9 + 0;

	# skip bad data
	next if $dBn == $nodata or $dBe == $nodata or $dBd == $nodata;

	# Extract date and data
	my $year  = $1; 
	my $month = $2;
	my $day   = $3;
	my $hour  = $4;
	my $minute= $5;

	# Fix 2 digit year
	$year  += 1900; $year += 100 if $year < 1950;

	# Spaces to 0-s.
	$month =~ s/ /0/g;
	$day   =~ s/ /0/g;
	$hour  =~ s/ /0/g;
	$minute=~ s/ /0/g;

	# Construct date in SWMF format
	my $date = "$year $month $day $hour $minute";

	# ignore multiple occurances
	next if $check{$station}{$date}{$type}++;

	# Check if there is GM/IE data for this date
	next unless $data{$station}{$date};

	# Add measured data
	$data{$station}{$date} .= " $dBn $dBe $dBd";

	# There is some data for this station
	$read{$station}{$type} = 1;
    }

}
#############################################################################
sub write_output{

    # write out perturbations for each station separately
    # calculate total (GM+IE) perturbation
    # calculate RMS error and correlation with
    # respect to data if present

    my $station;
    foreach $station (sort keys %stations){

	my $ndata = $ndataout;
	my $variables = $variablesout;

	if($read{"GM"} and not $read{$station}{"GM"}){
	    $variables =~ s/$variablesgm//;
	    $ndata -= 6;
	    print "No GM data for station $station\n";
	}
	if($read{"IE"} and not $read{$station}{"IE"}){
	    $variables =~ s/$variablesie//;
	    $ndata -= 6;
	    print "No IE data for station $station\n";
	}
	if($read{"data"} and not $read{$station}{"data"}){
	    $variables =~ s/$variablesdata//;
	    $ndata -= 3;
	    print "No measured data for station $station\n";
	}

	# No point writing output for stations with positions only
	next if $ndata < 4;

	# Construct header line
	my $header = "Magnetic perturbations from ";
	$header .= "GM and " if $read{$station}{"GM"};
	$header .= "IE and " if $read{$station}{"IE"};
	$header .= "measurements and " if $read{$station}{"data"};
	$header =~ s/and $//;

	# Process data and collect output to be written into file
	my $output; # string storing output

	# Variables for the error calculations
	my $ndate;
	my @rms_swmf;
	my @rms_data;
	my @rms_diff;

	my $date;
	foreach $date (sort keys %{ $data{$station} }){
	    my $data = $data{$station}{$date};
	    my @data = split(' ',$data);

	    # skip lines with missing data
	    next if $ndata != scalar @data;

	    # Add up the first two perturbations
	    my $sumBn = $data[3] + $data[6];
	    my $sumBe = $data[4] + $data[7];
	    my $sumBd = $data[5] + $data[8];

	    if($ndata >= 15){
		# Add up GM and IE components
		$sumBn += $data[9]  + $data[12];
		$sumBe += $data[10] + $data[13];
		$sumBd += $data[11] + $data[14];
	    }

	    # Print date, add a "msc" column with 0s, sums, and data
	    $output .= "$date 0 @data[0..2] ".
		"$sumBn $sumBe $sumBd @data[3..($ndata-1)]\n";

	    # no error calculation if there is no data
	    next unless $read{$station}{"data"};

	    # Count number of dates
	    $ndate++;

	    my $dataBn = $data[$ndata-3];
	    my $dataBe = $data[$ndata-2];
	    my $dataBd = $data[$ndata-1];

	    # Calculate averages and errors
	    $rms_swmf[0] += $sumBn**2;
	    $rms_swmf[1] += $sumBe**2;
	    $rms_swmf[2] += $sumBd**2;

	    $rms_data[0] += $dataBn**2;
	    $rms_data[1] += $dataBe**2;
	    $rms_data[2] += $dataBd**2;

	    $rms_diff[0] += ($sumBn - $dataBn)**2;
	    $rms_diff[1] += ($sumBe - $dataBe)**2;
	    $rms_diff[2] += ($sumBd - $dataBd)**2;

	}

	if($ndate > 0){
	    $header .= ", RMS:";
	    for my $i (0..2){

		$rms_swmf[$i] = sqrt($rms_swmf[$i] / $ndate);
		$rms_data[$i] = sqrt($rms_data[$i] / $ndate);
		$rms_diff[$i] = sqrt($rms_diff[$i] / $ndate) / $rms_data[$i];

		$header .= sprintf("%6.2f", $rms_diff[$i]);

		$rms{$station}[$i] = $rms_diff[$i];
	    }
	}

	my $fileout = "SWMF_mag_$station.log";
	print "writing $fileout\n";

	open FILE, ">$fileout";

	print FILE "$header\n";
	print FILE "$variables\n";
	print FILE $output;

	close FILE;
    }

}
##############################################################################
sub print_help{
    
    print "
Purpose: combine GM, IE and measured magnetic perturbations. Calculate errors.
The input files are parsed and the output is split into separate files for
each magnetometer. Only data points that occur in all 3 files are used. 

Usage:

PostMagnetometer.pl GM*.mag IE*.mag *.final

";
    exit 0;
}
