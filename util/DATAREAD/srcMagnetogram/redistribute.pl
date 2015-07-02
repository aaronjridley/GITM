#!/usr/bin/perl -s
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

my $Help        = $h or $H or $help;
my $Debug       = $D;
my $Snapshots   = $s;
my $FullGridOut = ($f or $g);
my $GhostCell   = ($g or 2);

use strict;

&help if $Help;

&help if @ARGV != 2;

# Range of snapshots to be read
my $firstread = 1;
my $lastread  = 100000;
my $dread     = 1;
if($Snapshots){
    if($Snapshots =~ /(^\d+)(:(\d+))?(:(\d+))?$/){
	$firstread=$1;
	$lastread= ($3 or $firstread);
	$dread   = ($5 or 1);
    }else{
	print "ERROR: inncorrect format -s=$Snapshots\n";
	&help;
    }
    print "firstread=$firstread lastread=$lastread stride=$dread\n" if $Debug;
}

my $filein  = $ARGV[0];
my $fileout = $ARGV[1];

# Check whether output files would overwrite input files
my $tmpfilein=$filein;   $tmpfilein  =~ s/(_np\d+)_\d\d\d\./$1./;
my $tmpfileout=$fileout; $tmpfileout =~ s/(_np\d+)_\d\d\d\./$1./;
die "Input and output file names are in conflict: $tmpfilein\n"
    if ($tmpfilein eq $tmpfileout);

my @npeDin  = &get_npeD($filein);
my @npeDout = &get_npeD($fileout);

my $ndimin = @npeDin;
my $ndimout= @npeDout;

die "Neither $filein nor $fileout contain processor number _np...\n"
    if $ndimin == 0 and $ndimout == 0;

# If one of the input and output files is not distributed, 
# take the dimensionality form the other one

$ndimin  or $ndimin = $ndimout;
$ndimout or $ndimout = $ndimin;

die "Dimensionality of input and output files differ: $ndimin != $ndimout\n"
    if $ndimin != $ndimout;

my $ndim = $ndimin;  # number of dimensions

print "ndim=$ndim\n" if $Debug;

# Calculate number of files for input and output

my $npein = 1; foreach my $npeD (@npeDin)  {$npein  *= $npeD};
my $npeout= 1; foreach my $npeD (@npeDout) {$npeout *= $npeD};

print "npein=$npein npeout=$npeout\n" if $Debug;

# Some global variables
my $filetype;   # 'ascii' or 'binary'
my $precision;  # 'real4' or 'real8'
my @ipeD;       # directional processor indexes
my $idim;       # dimension index

# Extract local grid sizes and global index ranges of input files.
# Calculate the size of the headers and the data records and the
# number of snapshots saved in the input files.

my @nxin;         # grid sizes for all input files
my @ixminin;      # minimum global indexes for all input files
my @ixmaxin;      # maximum global indexes for all input files

my $nw;           # number of flow variables
my $headersize;   # size of the header lines in bytes
my $headertail;   # constant string following the nx-s in the headers
my $recordsize;   # size of the variables in one line
my $nsnapshot;    # number of snapshots in the first input file

# Fill up @ipeD with ndim 0-s
for $idim (0..$ndim-1){$ipeD[$idim]=0};

for my $ipe (0..$npein-1){
    my $file = &set_filename($filein,$ipe);
    print "reading file=$file\n" if $Debug;

    my $type = (-T $file)? 'ascii' : 'binary';
    if(not $filetype){
	$filetype=$type;
    }else{
	die "file types differ for input files\n" if $type ne $filetype;
    }
    # Read nx from the file. Also check ndim and set nw and nsnapshot
    @{$nxin[$ipe]} = &read_header($file);

    # Set the index range for this file
    for $idim (0..$ndim-1){
        if($ipeD[$idim]==0){
	    $ixminin[$ipe][$idim] = 1;
	}else{
	    # Get processor index for the neighbor in the -idim-th direction
	    my $ipeprev;
            if   ($idim==0){$ipeprev=$ipe-1}
            elsif($idim==1){$ipeprev=$ipe-$npeDin[0]}
            else           {$ipeprev=$ipe-$npeDin[0]*$npeDin[1]}

	    $ixminin[$ipe][$idim] = $ixmaxin[$ipeprev][$idim]+1;
	}
	$ixmaxin[$ipe][$idim] = $ixminin[$ipe][$idim] + $nxin[$ipe][$idim] - 1;
    }

    if($Debug){
	print 
	    "ipeD = @ipeD\n",
	    "nxin     = @{$nxin[$ipe]}\n",
	    "ixminin  = @{$ixminin[$ipe]}\n",
	    "ixmaxin  = @{$ixmaxin[$ipe]}\n";
    }

    &increase_ipeD(@npeDin);
}

# Reduce $lastread to $nsnapshot if necessary
$lastread = $nsnapshot if $lastread > $nsnapshot;

# Calculate global grid size
my @nxall = @{$ixmaxin[$npein-1]};

# Calculate output grid sizes and index ranges
my @nxout;
my @ixminout;
my @ixmaxout;

# Get size of the grid for PE0
if($npeout>1){
    for $idim (0..$ndim-1){
        if($FullGridOut and $npeDout[$idim] > 1){
            # Take ghost cells into account
	    $nxout[0][$idim] = int(($nxall[$idim]-1-2*$GhostCell) / $npeDout[$idim]) + 1+$GhostCell;
	}else{
	    $nxout[0][$idim] = int(($nxall[$idim]-1) / $npeDout[$idim]) + 1;
	}
    }
}else{
    @{$nxout[0]} = @nxall;
}

# Fill up @ipeD with ndim 0-s
for $idim (0..$ndim-1){$ipeD[$idim]=0};

for my $ipe (0..$npeout-1){

    for $idim (0..$ndim-1){
	$ixminout[$ipe][$idim] = $ipeD[$idim] * $nxout[0][$idim] + 1;

	# Subtract ghost cell layers for full grid
	$ixminout[$ipe][$idim] -= $GhostCell*($ipeD[$idim]-1) 
	    if $FullGridOut and $ipeD[$idim]>1;

	$nxout[$ipe][$idim] = ($ipeD[$idim] < $npeDout[$idim]-1) ? 
	    $nxout[0][$idim] : ($nxall[$idim] - $ixminout[$ipe][$idim] + 1);

        # Subtract ghost cell layers for full grid
        $nxout[$ipe][$idim] -= $GhostCell if $FullGridOut
	    and $ipeD[$idim]>0 and $ipeD[$idim]<$npeDout[$idim]-1;

	$ixmaxout[$ipe][$idim] = $ixminout[$ipe][$idim]+$nxout[$ipe][$idim]-1;
    }

    if($Debug){
	print 
	    "ipeD = @ipeD\n",
	    "nxout    = @{$nxout[$ipe]}\n",
	    "ixminout = @{$ixminout[$ipe]}\n",
	    "ixmaxout = @{$ixmaxout[$ipe]}\n";
    }
    &increase_ipeD(@npeDout);
}

# Create output files by copying and merging data from the input files
my $snapshotsizein; # size of one snapshot in current input file
my $snapshotsizeout;# size of one snapshot in current output file
my @nxdata;         # size of the data copied from 1 input to 1 output file
my @shiftin;        # the shift of data in the input file relative to 1,1,1
my @shiftout;       # the shift of data in the output file relative to 1,1,1

my $nxsout;         # number of grid points for the output file
my $nxsin;          # number of grid points for the input file
for my $ipeout (0..$npeout-1){
    my $filepeout = &set_filename($fileout,$ipeout);
    my @headersaved;
    open(OUT,">$filepeout") or die "Could not open output file $filepeout\n";

    # Calculate snapshot size for output file
    $nxsout=1; for $idim (0..$ndim-1){$nxsout *= $nxout[$ipeout][$idim]}
    $snapshotsizeout = &snapshot_size($nxsout);
    print "snapshotsizeout = $snapshotsizeout\n" if $Debug;

    # Find all input files which can contribute to output file
    for my $ipein (0..$npein-1){
	# Check if there is an overlap
	&overlap($ipein,$ipeout) or next;

	# Calculate snapshot size for output file
	$nxsin=1; for $idim (0..$ndim-1){$nxsin *= $nxin[$ipein][$idim]}
	$snapshotsizein = &snapshot_size($nxsin);

	if($Debug){
	    print 
		"ipein $ipein and ipeout $ipeout overlap\n",
		"nxdata   = @nxdata\n",
		"shiftin  = @shiftin\n",
		"shiftout = @shiftout\n",
		"snapshotsizein  = $snapshotsizein\n";
	}

	my $filepein = &set_filename($filein,$ipein);
	open(IN,$filepein);

	print "read range=$firstread:$lastread:$dread\n" if $Debug;
	my $iread;
	my $isave;
	for($iread=$firstread; $iread<=$lastread; $iread+=$dread){
	    $isave++;
	    print "iread=$iread isave=$isave\n" if $Debug;

	    &write_header($ipein,$iread,$ipeout,$isave) unless
		$headersaved[$isave]++;

	    &write_data($ipein,$iread,$ipeout,$isave);
	}
	close(IN);
    }
    close(OUT);
}

##############################################################################
sub get_npeD{

    my $file=shift;
    my @npeD;

    if($file =~ /_np(\d+)/){
	my $npeD = $1;
	while($npeD =~ s/(..)//){push @npeD, $1+0};
    }else{
	@npeD=();
    }
    print "$file: npeD = @npeD\n" if $Debug;

    @npeD;
}

##############################################################################
sub increase_ipeD{
    my @npeD = @_;

    # Increase the dimensional processor index
    $ipeD[0]++;
    for my $idim (0..$ndim-2){
	if($ipeD[$idim] >= $npeD[$idim]){
	    $ipeD[$idim]=0;
	    $ipeD[$idim+1]++;
	}
    }
}

##############################################################################
sub set_filename{
    # Create filename based on base file name and processor number
    # E.g. $file="example_np0203.ini"; ipe=5 --> "example_np0203_005.ini"
    my $file = shift;
    my $ipe  = shift;

    $file =~ s/(_np\d+)(_\d\d\d)?/$1 . '_' . '0' x (3-length($ipe)) . $ipe/e;
    $file;
}

##############################################################################
sub snapshot_size{
    # Calculate the size of a snapshot based on $filetype, $precision, 
    # $headersize and the number of grid points (argument)
    my $nxs = shift;
    if($filetype eq 'ascii'){
	return $headersize + $recordsize*$nxs;
    }else{
	#       header     + data                        + data record markers
	return $headersize + $nxs*($ndim+$nw)*$precision + 8*(1+$nw);
    }
}

##############################################################################
sub overlap{
    # Check if the global index ranges of processors ipein and ipeout overlap
    my $ipein=shift;
    my $ipeout=shift;

    my $idim;

    # return false if the index ranges do not overlap in any dimensions
    for $idim (0..$ndim-1){
	return 0 if
	    $ixminin[$ipein][$idim] > $ixmaxout[$ipeout][$idim] or
	    $ixmaxin[$ipein][$idim] < $ixminout[$ipeout][$idim];
    }

    # Set global variables for data copying
    @nxdata = (1,1,1);
    for $idim (0..$ndim-1){
	my $ixmin = ($ixminin[$ipein][$idim] > $ixminout[$ipeout][$idim])? 
	    $ixminin[$ipein][$idim] : $ixminout[$ipeout][$idim];
	my $ixmax = ($ixmaxin[$ipein][$idim] < $ixmaxout[$ipeout][$idim])? 
	    $ixmaxin[$ipein][$idim] : $ixmaxout[$ipeout][$idim];

	$nxdata[$idim]  = $ixmax - $ixmin + 1;                # size of overlap
	$shiftin[$idim] = $ixmin - $ixminin[$ipein][$idim];   # shift in IN
	$shiftout[$idim]= $ixmin - $ixminout[$ipeout][$idim]; # shift in OUT
    }
    return  1;
}
##############################################################################
sub read_header{
    # Return grid size array from file $file
    # Set global variables $nw, $headersize, $recordsize, $nsnapshot if not set
    my $file=shift;
    my @nx;

    open(FILE,$file) or die "Could not open file $file\n";
    if($filetype eq 'ascii'){
	<FILE>;                           # header line

	my $line = <FILE>;                # it, time, ndim, neqpar, nw

	@nx = split(' ',<FILE>);          # nx

	if(not $nw){                      # if nw is still unknown read it
	    my @numbers = split ' ',$line;
	    my $ndimin = abs($numbers[2]);
	    die "ERROR: ndim in file is $ndimin != $ndim, the ndim for PE-s\n"
		if $ndimin != $ndim;
	    $nw = $numbers[4];
	    $headertail=<FILE>.<FILE>;      # equation parameters
                                            # and the variable names

            $headersize = tell(FILE);       # set the size of header lines
	    $recordsize = ($ndim+$nw)*18+1; # set the size of data records
	    # Calculate number of grid points
	    my $nxs=1; for $idim (0..$ndim-1){$nxs *= $nx[$idim]}
	    $nsnapshot = (-s $file)/($headersize+$nxs*$recordsize);

	    print
		"nw         = $nw\n",
		"headersize = $headersize\n",
		"recordsize = $recordsize\n",
		"nsnapshot  = $nsnapshot\n" if $Debug;

	}
    }else{
	# read 2nd line, get file type and extract ndim
	my $len2; # length of 2nd line in the header
	seek FILE, 500+8, 0;
	read FILE, $len2, 4;
	($len2)=unpack 'L',$len2;
	# double precision binary files have 24 byte long 2nd line
	$precision = ($len2==24)? 8 : 4;
	seek FILE,4+$precision,1;  # skip time step and time
	read FILE,my $line2,12;    # read ndim,neqpar,nw
        seek FILE,8,1;             # skip recordlengths
	read FILE,my $nx,$ndim*4;  # read nxs
	@nx=unpack("l*",$nx);
	if(not $nw){
	    (my $ndimin,my $neqpar,$nw)=unpack "lll",$line2;
	    $ndimin=abs($ndimin);
	    die "ERROR: ndim in file is $ndimin != $ndim, the ndim for PE-s\n"
		if $ndimin != $ndim;
	    # read header tail, ie. equation parameters and variable names
	    read FILE, $headertail, 4 + 8 + $neqpar*$precision + 8 + 500;
	    # Set header and record sizes
	    $headersize = tell(FILE);
	    $recordsize = $precision;
	    # Calculate number of grid points
            my $nxs=1; for $idim (0..$ndim-1){$nxs *= $nx[$idim]}
	    $nsnapshot = (-s $file)/&snapshot_size($nxs);

	    print
		"nw         = $nw\n",
		"headersize = $headersize\n",
		"recordsize = $recordsize\n",
		"nsnapshot  = $nsnapshot\n" if $Debug;

	}
    }
    close(FILE);
    print "$file: nx = @nx\n" if $Debug;
    @nx;
}

##############################################################################
sub write_header{
    # Copy header from IN into OUT but replace nx array with @$nx[$ipeout]

    my $ipein  = shift;
    my $iread  = shift;
    my $ipeout = shift;
    my $isave  = shift;

    # jump to beginning of snapshot in IN and OUT files
    seek IN,  ($iread-1)*$snapshotsizein,  0;
    seek OUT, ($isave-1)*$snapshotsizeout, 0;

    if($filetype eq 'ascii'){
	print OUT scalar <IN>;               # header string
	print OUT scalar <IN>;               # it, time, ndim, neqpar, nw
	printf OUT ("%8d" x $ndim)."\n",@{$nxout[$ipeout]}; # nx (replaced)
	print OUT $headertail;               # eqpar and variable names
    }else{
	read IN,my $headerhead, 500+8+8+$precision+16+4; # header string,
	print OUT $headerhead;                  # and it,time,ndim,neqpar,nw
	print OUT pack "L*",@{$nxout[$ipeout]}; # nx (replaced)
	print OUT $headertail;                  # eqpar and variable names

	# print out record length markers for the x array
	my $xsize = $ndim*$nxsout*$precision;
	print OUT pack "L",$xsize;  # record length at beginning
	seek OUT,$xsize,1;          # skip record
	print OUT pack "L",$xsize;  # record length at end

	# print out record length markers for the nw flow variables
	my $wsize = $nxsout*$precision;
	for (0..$nw-1){
	    print OUT pack "L",$wsize; # record length at beginning
	    seek OUT,$wsize,1;         # skip record
	    print OUT pack "L",$wsize; # record length at end
	}
    }
}

##############################################################################
sub write_data{

    # Copy data of size @nxdata shifted by indexes @shiftin in IN 
    # into OUT shifted by indexes @shiftout

    my $ipein = shift;
    my $iread = shift;
    my $ipeout= shift;
    my $isave = shift;

    my $nxin  = $nxin[$ipein][0];
    my $nxyin = $nxin[$ipein][1] * $nxin;

    my $nxout = $nxout[$ipeout][0];
    my $nxyout= $nxout[$ipeout][1] * $nxout;
    for my $k (0..$nxdata[2]-1){
	for my $j (0..$nxdata[1]-1){

	    my $data;
	    my $offsetin = ($iread-1)*$snapshotsizein + 
		$headersize +
		( ($k+$shiftin[2])*$nxyin +
		  ($j+$shiftin[1])*$nxin +
		  $shiftin[0]
		  ) * $recordsize;

	    my $offsetout = ($isave-1)*$snapshotsizeout +
		$headersize +
		( ($k+$shiftout[2])*$nxyout +
		  ($j+$shiftout[1])*$nxout +
		  $shiftout[0]
		  ) * $recordsize;

	    if($filetype eq 'ascii'){
		seek IN,$offsetin,0;
		read IN,$data,$nxdata[0]*$recordsize;
		seek OUT,$offsetout,0;
		print OUT $data;
	    }else{
		$offsetin  += 4; # length of record marker for x array
		$offsetout += 4; # length of record marker for x array
		for (0..$ndim-1){
		    seek IN,$offsetin,0;
		    read IN,$data,$nxdata[0]*$recordsize;
		    seek OUT,$offsetout,0;
		    print OUT $data;
		    $offsetin  += $precision*$nxsin;    # length of x(...,idim)
		    $offsetout += $precision*$nxsout;
		}
		$offsetin  += 8; # length of record marker for x and first var
		$offsetout += 8; # length of record marker for x and first var
		for (0..$nw-1){
		    seek IN,$offsetin,0;
		    read IN,$data,$nxdata[0]*$recordsize;
		    seek OUT,$offsetout,0;
		    print OUT $data;
		    $offsetin  += $precision*$nxsin +8; # length of var+markers
		    $offsetout += $precision*$nxsout+8;
		}
	    }	    
	}
    }
}
##############################################################################
sub help{

    print "
Convert data files with some MPI distribution into another distribution.

   redistribute.pl [-h] [-D] [-s=SNAPSHOTS] [-f] [-g=GHOSTCELL] filein fileout

-h            print this help message
-D            print debug information
-s=NUMBER     convert snapshot NUMBER (first snapshot is indexed by 1)
-s=FROM:TO    convert snapshots from number FROM to number TO
-s=FROM:TO:STRIDE 
              convert every STRIDE-th snapshot from number FROM to number TO
-f            full grid output with 2 ghost cell layers at the outer boundary
-g=GHOSTCELL  full grid output with GHOSTCELL layers at the outer boundary
filein        name of the input  files without/ignored processor rank '_nnn'
fileout       name of the output files without/ignored processor rank '_nnn'

Examples:

Merge data files on 1 x 2 x 2 PE-s into a single file:

   redistribute.pl fdips_field_np010202.out fdips_field.out

Scatter the first snapshot of a single data file into 3 x 1 PE-s:

   redistribute.pl -s=1 exampleA22.out exampleA22_np0301.out

Rearrange a 2 x 2 distribution into a 1 x 3 distribution with a new name:

   redistribute.pl old_np0202.out new_np0103.out

Rearrange snapshots 2,4,6 of a 2x2x2 distribution into 4x2x1 distribution:

   redistribute.pl -s=2:6:2 old_np020202.out new_np040201.out

Create a file with the default 2 ghost layers at the outer boundary:

   redistribute.pl -f old.out new_np0202.out

Create a file with 4 ghost layers at the outer boundary:

   redistribute.pl -g=4 old.out new_np0202.out

Processor rank in the filenames (if any) is ignored, so these two are the same:

   redistribute.pl old_np020202_000.out new_np040201_003.out
   redistribute.pl old_np020202.out     new_np040201.out

";
    exit;
}
##############################################################################
