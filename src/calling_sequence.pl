#!/usr/bin/perl -s
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

$help = $h;

%filenames = ();

sub find_subroutine {

    local($subroutine,$foundfile,$file,$IsDone,$line,$StartLooking,$IsFound);
    local($deep,$space,$i,$call,$subfile,$sfile,$scall);

    $subroutine = $_[0];
    $deep = $_[1]+0;

    $space = "";
    for ($i=0;$i<$deep;$i++){
	$space = "  $space";
    }
    $scall = $subroutine;
    $scall =~ s/_/\\_/g;
    print "\n\n$space\\item {\\tt $scall} ";

    if ($filenames{$subroutine}) {
	#print "FOUND IN HASH!!\n";
	$foundfile = $filenames{$subroutine};
    } else {

	my @files = <*.f90>;

	foreach $file (@files) {
	    open(INFILE,"<$file");
	    while (<INFILE>){
		if (/^\s*subroutine\s$subroutine\W+/) {
		    if (!$foundfile) {
			$foundfile = $file;
			$filenames{$subroutine} = $foundfile;
		    }
		}
	    }
	    close(INFILE);
	}

    }

    if ($foundfile) {
	$sfile = $foundfile;
	$sfile =~ s/_/\\_/g;
	print "  in file {\\bf $sfile}\n";
	open(INFILE,"<$foundfile");
	$IsDone = 0;
	$StartLooking = 0;
	$IsFound = 0;
	while (!$IsDone) {
	    $line = <INFILE>;
	    if ($line =~ /^\s*subroutine\s$subroutine\W+/) {
		$StartLooking = 1;		
	    }
	    if ($StartLooking) {
		if ($line =~ /call\s(\w*)/) {
		    if (!$IsFound) {
			$IsFound = 1;
			print "$space  \\begin{itemize}\n";
		    }
		    $call = $1;

                    #while (($key, $value) = each(%filenames)){
                    #     print "$key --> $value\n";
                    #}

		    if ($filenames{$call}) {
			#print "FOUND IN HASH!!\n";
			$subfile = $filenames{$call};
		    } else {

			my @files = <*.f90>;
			
			$subfile=0;
			foreach $file (@files) {
			    #print "$file\n";
			    open(INFILE2,"<$file");
			    while (<INFILE2>){
				if (/^\s*subroutine\s$call\W+/) {
				    if (!$subfile) {
					$subfile = $file;
					$filenames{$call} = $subfile;
				    }
				}
			    }
			    close(INFILE2);
			}

		    }

		    $call =~ s/_/\\_/g;
		    print "$space    \\item {\\tt $call}";
		    if ($subfile) {
			$subfile =~ s/_/\\_/g;
			print " in file $subfile";
		    }
		    print "\n";

		}
		$IsDone = 1 if ($line =~ /^\s*end subroutine\s$subroutine\W+/);
		print "$space  \\end{itemize}\n" if ($IsFound && $IsDone);
	    }
	}
	close(INFILE);
    }


}


@files = <*.f90>;
#foreach $file (@files) {
#   print $file . "\n";
#} 

open(MAIN,"<main.f90");

print "{\\tt main program} in {\\bf main.f90}\n\n";
print "\\begin{itemize}\n\n";

$depth = 0;
while (<MAIN>) {

    if (/while/) {
	print "\n\n\\item Loop Start\n";
	print "  \\begin{itemize}\n";
	$depth++;
    }

    if (/call\s(\w*)/) {
	&find_subroutine($1,$depth);
    }

    if (/enddo/) {
	print "\n\n  \\end{itemize} %loop \n";
	print "\n\n \\item Loop End\n";
	$depth--;
    }


}

close(MAIN);

print "\\end{itemize}\n";

exit(1);

