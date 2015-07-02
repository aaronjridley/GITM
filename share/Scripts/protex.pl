#!/usr/bin/perl
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf
#
#  $Id: protex.pl,v 1.4 2013/10/12 04:01:02 kopmanis Exp $
#
#BOP
#
#!ROUTINE: ProTeX v. 2.10 - Translates Prologues to LaTeX
#
#!INTERFACE:
#         protex [-hbACFS] ] [+-nlsxf] [srcfile(s)]
#
#!DESCRIPTION:
#         Perl script to produce a \LaTeX compatible document 
#         from a Fortran and other source code with standard Pro\TeX 
#         prologues. If source files are not specified it
#         reads from stdin; output is always to stdout.
# 
# \noindent        
# {\bf Command Line Switches:} \vspace{0.2cm}
#
# \begin{center}
# \begin{tabular}{|c|l|} \hline \hline
#   -h   & Help mode: list command line options   \\ \hline
#   -b   & Bare mode, meaning no preamble, etc.  \\ \hline
#   +/-n & New Page for each subsection (wastes paper) \\ \hline
#   +/-l & Listing mode, default is prologues only \\ \hline
#   +/-s & Shut-up mode, i.e., ignore any code from BOC to EOC \\ \hline
#   +/-x & No LaTeX mode, i.e., put !DESCRIPTION: in verbatim mode \\ \hline
#   +/-f & No source file info \\ \hline
#   -A   & Ada code \\ \hline
#   -C   & C++ code \\ \hline
#   -F   & F90 code (default) \\ \hline
#   -L   & LaTex file \\ \hline
#   -S   & Shell script \\ \hline \hline
# \end{tabular}
# \end{center}
#
# The options can appear in any order.  The options, -h and -b, affect
# the input from all files listed on command-line input.  Each of the
# remaining options effects only the input from the files listed after
# the option and prior to any overriding option.  The plus sign
# turns off the option.  For example, the command-line input,
# \begin{verbatim}
#      protex -bnS File1 -F File2.f +n File3.f
# \end{verbatim}
# will cause the option, {\tt -n} to affect the input from the files,
# {\tt File} and {\tt File2.f}, but not from {\tt File3.f}.  The
# {\tt -S} option is implemented for {\tt File1} but is overridden by
# the {\tt -F} for files {\tt File2.f} and {\tt File3.f}.
#
#!REVISION HISTORY:
#  20Dec1995  da Silva  First experimental version
#  10Nov1996  da Silva  First internal release (v1.01)
#  28Jun1997  da Silva  Modified so that !DESCRIPTION can appear after
#             !INTERFACE, and !INPUT PARAMETERS etc. changed to italics.
#  02Jul1997  Sawyer    Added shut-up mode
#  20Oct1997  Sawyer    Added support for shell scripts
#  11Mar1998  Sawyer    Added: file name, date in header, C, script support
#  05Aug1998  Sawyer    Fixed LPChang-bug-support-for-files-with-underscores
#  10Oct1998  da Silva  Introduced -f option for removing source file info
#                       from subsection, etc.  Added help (WS).
#  06Dec1999  C. Redder Added LaTeX command "\label{sec:prologues}" just 
#                       after the beginning of the proglogue section.
#  13Dec1999  C. Redder Increased flexbility in command-line
#                       interface.  The options can appear in any
#                       order which will allow the user to implement
#                       options for select files.
#  01Feb1999  C. Redder Added \usepackage commands to preamble of latex
#                       document to include the packages amsmath, epsfig
#                       and hangcaption.
#  10May2000  C. Redder Revised LaTeX command "\label{sec:prologues}"
#                       to "\label{app:ProLogues}"
#  10/10/2002 da Silva  Introduced ARGUMENTS keyword, touch ups.
#  08/09/2003 G. Toth   Using Perl string matching instead of the splitting
#                       into the Fld array, so that comment characters can 
#                       be preceeded with blanks and there is no need to 
#                       repeat the exclamation mark in Fortran 90.
#                       Removed '#end if', '#end for' comments. Replaced 
#                       @author, @affiliation, @title, @date with strings.
#                       Added new key words: 
#                         INPUT ARGUMENTS, 
#                         OUTPUT ARGUMENTS, 
#                         INPUT/OUTPUT ARGUMENTS
#                       Added INPUT/OUTPUT PARAMETERS->ARGUMENTS replacement.
#                       Added \bigskip before CONTENTS (needed when 
#                       CONTENTS (!BOC) immediately follows !DESCRIPTION:
#                       Nicer format for Tex output: newline after 
#                       end{verbatim} etc.
#                       Using single quote instead of double quote when useful.
#  08/27/2003 G. Toth   Added -L option for Latex files, which are simply 
#                       copied. Useful for inserting Latex between source files
#                       Commented out the section Routine/Function...
#  08/07/2004 G. Toth   Corrected the self-description.
#                       Changed the version number to 2.10
#                       Renamed the script from protex to protex.pl
#  08/08/2004 G. Toth   Some clean up of source code quoted from Makefiles.
#                       
#EOP
#----------------------------------------------------------------------------

# Keep this if you don't know what it does...
# -------------------------------------------
  $[ = 1;                 # set array base to 1
  $, = ' ';               # set output field separator
  $\ = "\n";              # set output record separator


### This is just a temporary fix to start modules on a new page
### even if the latex file is put together from individual parts
$not_first=1;

# Set valid options lists
# -----------------------
  $GlobOptions = 'hb';    # Global options (i.e for all files)
  $LangOptions = 'ACLFS'; # Options for setting programming languages
  $SwOptions   = 'flnsx'; # Options that can change for each input 
                          #   file
  $RegOptions  = "$GlobOptions$LangOptions";
                          # Scan for global options until first first
                          #   file is processed.

# Scan for global options
# -----------------------
  $NFiles = 0;
  ARG: foreach $arg (@ARGV) {
      $option = &CheckOpts ( $arg, $RegOptions, $SwOptions ) + 1;
      if ( $option ) {
	  $rc = &GetOpts    ( $arg, $GlobOptions );
	  next ARG;
      } else { 
	  $NFiles++;
      }
  }

# If all input arguments are options, then assume the
# filename, "-", for the standard input
# --------------------------------------------------
  if ( $NFiles == 0 ) { push (@ARGV, "-"); } 

# Implement help option
# ---------------------
  if ( $opt_h ) {
      &print_help;
      exit;
  }

# Optional Prologue Keywords
# --------------------------
  @keys = ( "!INTERFACE:",
            "!USES:",
            "!PUBLIC TYPES:",
            "!PRIVATE TYPES:",
            "!PUBLIC MEMBER FUNCTIONS:",
            "!PRIVATE MEMBER FUNCTIONS:",
            "!PUBLIC DATA MEMBERS:",
            "!PARAMETERS:",
            "!ARGUMENTS:",
            "!DEFINED PARAMETERS:",
            "!INPUT PARAMETERS:",
            "!INPUT/OUTPUT PARAMETERS:",
            "!OUTPUT PARAMETERS:",
            "!INPUT ARGUMENTS:",
            "!OUTPUT ARGUMENTS:",
            "!INPUT/OUTPUT ARGUMENTS:",
            "!RETURN VALUE:",
            "!REQUIREMENTS:",
            "!REVISION HISTORY:",
            "!BUGS:",
            "!SEE ALSO:",
            "!SYSTEM ROUTINES:",
            "!FILES USED:",
            "!REMARKS:",
            "!TO DO:",
            "!CALLING SEQUENCE:",
            "!AUTHOR:",
            "!CALLED FROM:",
            "!LOCAL VARIABLES:" );

# Initialize these for clarity
# ----------------------------
  $intro = 0;             # doing introduction?
  $prologue = 0;          # doing prologue?
  $first = 1;             # first prologue?
  $source = 0;            # source code mode?
  $verb = 0;              # verbatim mode?
  $tpage = 0;             # title page?
  $begdoc = 0;            # has \begin{document} been written?

# Initial LaTeX stuff
# -------------------
  &print_notice();
  &print_preamble();      # \documentclass, text dimensions, etc.
  &print_macros();        # short-hand LaTeX macros

# Main loop -- for each command-line argument
# -------------------------------------------
ARG:
  foreach $arg (@ARGV) {

#    Scan for non-global command-line options
#    ----------------------------------------
     $option = &CheckOpts ( $arg, $RegOptions, $SwOptions, "quiet" ) + 1;
     if ( $option ) {
        &GetOpts  ( $arg, $SwOptions   );
        &SetOpt   ( $arg, $LangOptions );
        next ARG;

}#   end if

#    Determine the type of code, set corresponding search strings
#    ------------------------------------------------------------

     $comment_string = '!'; # if $opt_F  # FORTRAN (default so if is commented out)
     $comment_string = '--'   if $opt_A; # ADA
     $comment_string = '//'   if $opt_C; # C
     $comment_string = '#'    if $opt_S; # Script

     $boi_string = $comment_string.'BOI';
     $eoi_string = $comment_string.'EOI';
     $bop_string = $comment_string.'BOP';
     $eop_string = $comment_string.'EOP';
     $boc_string = $comment_string.'BOC';
     $eoc_string = $comment_string.'EOC';


#    Set file name parameters
#    ------------------------
     $InputFile           = $arg;
     @all_path_components = split( /\//, $InputFile     );
     $FileBaseName        = pop  ( @all_path_components );
     $FileBaseName        =~ s/_/\\_/g;
     if ( $InputFile eq "-" ) {$FileBaseName = "Standard Input";}

#    Set date
#    --------
     $Date                = `date`;

#    Open current file
#    -----------------
     open ( InputFile, "$InputFile" )
          or print STDERR "Unable to open $InputFile: $!";

#    Print page header
#    -----------------
     ### This is an attempt to mark the header with the current
     ### source file. Not perfect yet, need to look at Latex book... 
     #printf "\n\\markboth{Source File: %s}{Source File: %s}\n\n",
     #                               $FileBaseName, $FileBaseName;

LINE:
#    Inner loop --- for processing each line of the input file
#    ---------------------------------------------------------
     while ( <InputFile> ) {

#       LATEX files are simply copied
#       -----------------------------
	if ( $opt_L ){
	     print;
	     next LINE;
	}

        chop;     # strip record separator

#       !PARAMETERS: really mean !ARGUMENTS:
#       ------------------------------------
        s/!(INPUT\/OUTPUT|INPUT|OUTPUT)\s+PARAMETERS:/!$1 ARGUMENTS:/g;

#       Straight quote
#       --------------
        if (/^\s*($comment_string)?\s*!QUOTE:/) {
	    print $';
	    print ' ';
	    next LINE;
	}

#       Handle optional Title Page and Introduction
#       -------------------------------------------
        if (/^\s*$boi_string\b/) {
	    print ' ';
	    $intro = 1;
	    next LINE;
        }

        if (/^\s*($comment_string)?\s*!TITLE:/) {
	    if ( $intro ) {
		$title = $';
		$tpage = 1;
		next LINE;
	    }
	}

        if (/^\s*($comment_string)?\s*!AUTHORS:/) {
	    if ( $intro ) {
		$author = $';
		$tpage = 1;
		next LINE;
	    }
	}

        if (/^\s*($comment_string)?\s*!AFFILIATION:/) {
	    if ( $intro ) {
		$affiliation = $';
		$tpage = 1;
		next LINE;
	    }
	}

        if (/^\s*($comment_string)?\s*!DATE:/) {
	    if ( $intro ) {
		$date = $';
		$tpage = 1;
		next LINE;
	    }
	}

        if (/^\s*($comment_string)?\s*!INTRODUCTION:/) {
	    if ( $intro ) {
		&do_beg();
		print ' ';
		print '%..............................................';
		print "\\section{$'}";
		print ' ';
		next LINE;
	    }
	}


#       End of introduction
#       -------------------
        if (/^\s*$eoi_string\b/) {
           print "\n";
           print '%/////////////////////////////////////////////////////////////',;
           print '\newpage',"\n";
           $intro = 0;
           next LINE;
}#      end if

#       Beginning of prologue
#       ---------------------
        if (/^\s*$bop_string\b/) {
           if ( $source ) { &do_eoc(); }
           print ' ';
           print '%/////////////////////////////////////////////////////////////';
           &do_beg();
           if ($first == 0) {
              ### print "\\newpage";
              print " ";
              print "\\mbox{}\\hrulefill\\ ";
              print " ";}
           else {
# Commented out.
#              unless($opt_b){print "\\section{Routine/Function Prologues} \\label{app:ProLogues}";}
	   }

           $first = 0;
           $prologue = 1;
           $verb = 0;
           $source = 0;
           &set_missing();   # no required keyword yet
           next LINE;
}#      end if

#       A new subroutine/function
#       -------------------------
        if (/^\s*($comment_string)?\s*!ROUTINE:/) { 
	    if ($prologue) {
		$_ = $';
		$name_is = $_;
		s/_/\\_/g;                         # Replace "_" with "\_"
		if ( $opt_n && $not_first ) { 
		    printf "\\newpage\n";
		    ### This is added to mark the header with the current
		    ### source file
		    printf "\n\\markboth{Source File: %s}{Source File: %s}\n\n"
                                    ,$FileBaseName, $FileBaseName;
		}
		unless ($opt_f) {printf "\\subsubsection{%s (Source File: %s)}\n\n", $_, $FileBaseName;}
		else            {printf "\\subsubsection{%s }\n\n", $_;}
		$have_name = 1;
		$not_first = 1;
		next LINE;
	    }
	}

#       A new Module
#       ------------
        if (/^\s*($comment_string)?\s*!MODULE:/) { 
	    if ($prologue) {
		$_ = $';
		$name_is = $_;
		s/_/\\_/g;                         # Replace "_" with "\_"
		if ( $opt_n && $not_first ) {
		    printf "\\newpage\n";
		    ### This is added to mark the header with the current
		    ### source file. 
		    printf "\n\\markboth{Source File: %s}{Source File: %s}\n\n"
                                    ,$FileBaseName, $FileBaseName;
		}
		unless($opt_f) {
#		    printf "\\subsection{Fortran:  Module Interface %s (Source File: %s)}\n\n", $_, $FileBaseName;
		    printf "\\subsection{Module %s (Source File: %s)}\n\n", $_, $FileBaseName;
		}else{
#		    printf "\\subsection{Fortran:  Module Interface %s }\n\n", $_;}
		    printf "\\subsection{Module %s }\n\n", $_;}
		$have_name = 1;
		$have_intf = 1;  # fake it, it does not need one.
		$not_first = 1;
		next LINE;
	    }
	}

#       A new include file
#       ------------------
        if (/^\s*($comment_string)?\s*!INCLUDE:/) { 
	    if ($prologue) {
		$_ = $';
		$name_is = $_;
		s/_/\\_/g;                         # Replace "_" with "\_"
		if ( $opt_n && $not_first ) { printf "\\newpage\n"; }
		unless($opt_f) {
		    printf "\\subsubsection{Include File %s (Source File: %s)}\n\n", $_, $FileBaseName;
		}else{
		    printf "\\subsubsection{Include File %s }\n\n", $_;
		}
		$have_name = 1;
		$have_intf = 1;  # fake it, it does not need one.
		$not_first = 1;
		next LINE;
	    }
	}

#       A new INTERNAL subroutine/function
#       ----------------------------------
        if (/^\s*($comment_string)?\s*!IROUTINE:/) {            # Internal routine
	    if ($prologue) {
		$_ = $';
		$name_is = $_;
		s/_/\\_/g;                        # Replace "_" with "\_"
		printf "\\subsubsection{%s}\n\n", $_;
		$have_name = 1;
		next LINE;
	    }
	}

#       A new CLASS
#       ------------
        if (/^\s*($comment_string)?\s*!CLASS:/) { 
	    if ($prologue) {
		$_ = $';
		$name_is = $_;
		s/_/\\_/g;                         # Replace "_" with "\_"
		if ( $opt_n && $not_first ) { printf "\\newpage\n"; }
		unless($opt_f) {
		    printf "\\subsection{C++:  Class Interface %s (Source File: %s)}\n\n", $_, $FileBaseName;
		}else{
		    printf "\\subsection{C++:  Class Interface %s }\n\n", $_;}
		$have_name = 1;
		$have_intf = 1;  # fake it, it does not need one.
		$not_first = 1;
		next LINE;
	    }
	}

#       A new Method
#       -------------------------
        if (/^\s*($comment_string)?\s*!METHOD:/) { 
	    if ($prologue) {
		$_ = $';
		$name_is = $_;
		s/_/\\_/g;                         # Replace "_" with "\_"
		if ( $opt_n && $not_first ) { printf "\\newpage\n"; }
		unless ($opt_f) {
		    printf "\\subsubsection{%s (Source File: %s)}\n\n", $_, $FileBaseName;
		}else{
		    printf "\\subsubsection{%s }\n\n", $_;
		}
		$have_name = 1;
		$not_first = 1;
		next LINE;
	    }
	}

#       A new function
#       -------------------------
        if (/^\s*($comment_string)?\s*!FUNCTION:/) { 
	    if ($prologue) {
		$_ = $';
		$name_is = $_;
		s/_/\\_/g;                         # Replace "_" with "\_"
		if ( $opt_n && $not_first ) { printf "\\newpage\n"; }
		unless ($opt_f) {
		    printf "\\subsubsection{%s (Source File: %s)}\n\n", $_, $FileBaseName;
		}else{
		    printf "\\subsubsection{%s }\n\n", $_;
		}
		$have_name = 1;
		$not_first = 1;
		next LINE;
	    }
	}

#       Description: what follows will be regular LaTeX (no verbatim)
#       -------------------------------------------------------------
        if (/^\s*($comment_string)?\s*!DESCRIPTION:/) {
           if ($prologue) {
	       if ($verb) {
		   print '\end{verbatim}';
		   print '{\sf DESCRIPTION:}\\\\';
		   $verb = 0; 
	       }else{                          # probably never occurs
	       }
	       if ($opt_x) {
		   printf '\begin{verbatim} ';
		   $verb = 1;
		   $first_verb = 1; 
	       }
	       printf '%s ', $';

	       ### print " ";
	       $have_desc = 1;
	       next LINE;
	   }
       }

#       Handle optional keywords (these will appear as verbatim)
#       --------------------------------------------------------
        if ($prologue) {
	  KEY:foreach $key ( @keys ) {
              if ( /$key/ ) {
		  if ($verb) {
		      print '\end{verbatim}';
		      $verb = 0; 
		  }else{
		      print '';
		      print '\bigskip';
		  }
		  $k = sprintf('%s', $key);
		  $ln = length($k);
		  ###printf "\\subsubsection*{%s}\n", substr($k, 2, $ln - 1);
		  ###printf "{\\Large \\em %s}\n", ucfirst lc substr($k, 2, $ln - 1);
		  $_ = $key;
		  if( /USES/ || /INPUT/ || /OUTPUT/ || /PARAMETERS/ || 
		      /VALUE/  || /ARGUMENTS/ ) {
		      printf "{\\em %s}\n", substr($k, 2, $ln - 1); # italics
		  }else{
		      printf "{\\sf %s}\n", substr($k, 2, $ln - 1); # san serif
		  }

		  printf "\\begin{verbatim} ";
		  $verb = 1;
		  $first_verb = 1;
		  if ( $key eq "!INTERFACE:" )        { $have_intf = 1; }
		  if ( $key eq "!CALLING SEQUENCE:" ) { $have_intf = 1; }
		  if ( $key eq "!REVISION HISTORY:" ) { $have_hist = 1; }
		  next LINE;
	      }
          # end foreach KEY
	  }
	}

#       End of prologue
#       ---------------
        if (/^\s*$eop_string\b/) {
	    if ($verb) {
		print '\end{verbatim}';
		$verb = 0;
	    }
	    $prologue = 0;
#           &check_if_all_there(); # check if all required keyword are there.
	    if ( $opt_l ) {
		s/$eop_string/$boc_string/; # Replace EOP with BOC: code begins
	    }else{
		next LINE;
	    }
}#      end if

        unless ( $opt_s ) {
#
#           Beginning of source code section
#           --------------------------------
	    if (/^\s*$boc_string\b/) {
		print '';
		print '';
		print '%/////////////////////////////////////////////////////////////';
		$first = 0;
		$prologue = 0;
		$source = 1;
		print '';
		print '\bigskip';
		print '{\sf CONTENTS:}';
		print '\begin{verbatim}';
		$verb = 1;
		next LINE;
	    }

#           End of source code
#           ------------------
	    if (/^\s*$eoc_string\b/) {
		&do_eoc();
		$prologue = 0;
		next LINE;
}#         end if
}#      end unless
   
#   Prologue or Introduction, print regular line
#   --------------------------------------------
    if ($prologue or $intro) {
        if ( $verb and /^\s*$comment_string\s*$/ ) {
	    next LINE;                # to eliminate excessive blanks 
	}
	if ( /^\s*$comment_string\s*\\ev\b/ ) {   # special handling
	    $_ = $comment_string . '\end{verbatim}';
	}
	s/^(\s*)$comment_string/$1 /;  # replace comment string with blank

        unless ( $first_verb ) { printf "\n "; }
        printf "%s",$_;
        $first_verb = 0;
        next LINE;
    }

#   Source code: print the full line
#   --------------------------------
    if ($source) {
        # Do not print config lines from Makefiles
        next LINE if /^\s*\@\#\^ CMP/x;
        # replace @echo '...' and @echo "..." with ...
	s/\@echo\s+\'([^\']*)\'/$1/; 
	s/\@echo\s+\"([^\"]*)\"/$1/; 
        # print the line
        print  $_;
        next LINE;
    }

}#   end inner loop for processing each line of the input file
 #   ---------------------------------------------------------

}# end main loop for each command-line argument
 # --------------------------------------------
  print $_;
  if ( $source ) { &do_eoc(); }     
  print '%...............................................................';

  unless ( $opt_b ) {
      print "\\end{document}";
  }


#----------------------------------------------------------------------

  sub CheckOpts
#    Checks options against a given list.  Outputs error message
#    for any invalid option.
#
#    Usage:
#       $rc = &CheckOpts ( options, valid_reg_options,
#                                   valid_sw_options,
#                                   quiet_mode )
#
#       character: options - options to be checked. (e.g. -df+x)  The
#                            list must begin with a positive or
#                            negative sign.  If no sign appears at the
#                            beginning or by itself, then the argument
#                            is not recognized as a list of options.
#       character: valid_reg_options - list of valid regular options.
#                            (i.e. options that are associated only
#                            eith negative sign.)
#       character: valid_sw_options - list of valid switch options.
#                            (i.e. options that can be associated with
#                            either a positive or negative sign.
#       logical:   quiet mode (optional) If true then print no error
#                            messages.
#       integer:   rc      - return code
#                            = -1 if the arguement, options, is
#                               not recognized as a list of options
#                            =  0 if all options are valid.
#                            >  0 for the number of invalid options.
# 
{    local($options,
           $valid_reg_options,
           $valid_sw_options,
           $quiet_mode ) = @_;

     if ( $options eq "+" ||
          $options eq "-" ) {return -1}

     local(@Options) = split( / */, $options );
     if ( $Options[ $[ ] ne "-" &&
          $Options[ $[ ] ne "+" ) {return -1;}

     local($option, $option_sign, $valid_list, $pos);
     local($errs)    = 0;
     foreach $option ( @Options ) {
        if ( $option eq "-" ||
             $option eq "+" ) {$option_sign = $option;}
        else {
           if ( $option_sign eq "-" )
              { $valid_list = $valid_reg_options
                            . $valid_sw_options; }
           else
              { $valid_list = $valid_sw_options; }
           $pos = index ($valid_list,$option); 
           if ( $pos < $[ &&
                $quiet_mode ) {
              $errs++;
              print STDERR "Invalid option: $option_sign$option \n"; 

}#         end if
}#      end if
}#   end foreach
     return $errs;

}#end sub GetOpts

  sub GetOpts
#    Gets options.  If an option is valid,  then opt_[option] is
#    set to 0 or 1 as a side effect if the option is preceeded by
#    a positive or negative sign.
#
#    Usage:
#       $rc = &GetOpts ( options, valid_options )
#
#       character: options - options to be checked. (e.g. -df+x)  The
#                            list must begin with a positive or
#                            negative sign.  If no sign appears at the
#                            beginning or by itself, then the argument
#                            is not recognized as a list of options.
#       character: valid_options - list of valid options (e.g. dfhx)
#       integer:   rc      - return code
#                            = -1 if the arguement, options, is
#                               not recognized as a list of options.
#                            =  0 otherwise
# 
{    local($options,$valid_options) = @_;

     if ( $options eq "+" ||
          $options eq "-" ) {return -1}

     local(@Options)       = split( / */, $options );
     if ( $Options[ $[ ] ne "-" &&
          $Options[ $[ ] ne "+" ) {return -1;}

     local($option, $option_sign);

     foreach $option ( @Options ) {

        if ( $option eq "-" ||
             $option eq "+" ) {
           $option_sign = $option; }

        else {

           if ( index ($valid_options,$option) >= $[ ) {
              if ( $option_sign eq "-" ) {${"opt_$option"} = 1;}
              if ( $option_sign eq "+" ) {${"opt_$option"} = 0;};

}#         end if
}#      end if
}#   end foreach

     return 0;
}#end sub GetOpts

  sub SetOpt
#    Sets option flags.  For the last input option that is in a
#    list, the flag opt_[option] is set to 1 as a side effect.
#    For all other options in the list, opt_[option] is set to 0.
#
#    Usage:
#       $rc = &SetOpt ( options, valid_options )
#
#       character: options - options to be checked. (e.g. -df+x)  The
#                            list must begin with a positive or
#                            negative sign.  If no sign appears at the
#                            beginning or by itself, then the argument
#                            is not recognized as a list of options.
#       character: valid_options - list of valid options (e.g. def )
#       integer:   rc      - return code
#                            = -1 if the arguement, options, is
#                               not recognized as a list of options.
#                            =  0 otherwise
#       Note: For the examples provided for the input arguments,
#             $opt_d = 0, $opt_e = 0, and $opt_f = 1, since the 
#             input option, -f, was the last in the argument,
#             option.
# 
{    local($options,$valid_options) = @_;

     if ( $options eq "+" ||
          $options eq "-" ) {return -1}

     local(@Options)       = split( / */, $options       );
     local(@ValidOptions)  = split( / */, $valid_options );
     if ( $Options[ $[ ] ne "-" &&
          $Options[ $[ ] ne "+" ) {return -1;}

     local($option, $option_sign);

     foreach $option ( @Options ) {
        if ( $option ne "-" &&
             $option ne "+" ) {

           if ( index ($valid_options,$option) >= $[ ) {
              foreach $valid_option (@ValidOptions ) {
                 ${"opt_$valid_option"} = 0;

}#            end foreach
              ${"opt_$option"} = 1;
}#         end if
}#      end if
}#   end foreach

  return 0;
}#end sub SetOpt

sub print_help {

    print "Usage:     protex [-hbACFS] [+-nlsxf] [src_file(s)]";
    print " ";
    print " Options:";
    print "     -h   Help mode: list command line options";
    print "     -b   Bare mode, meaning no preamble, etc."; 
    print "     +-n  New Page for each subsection (wastes paper)";
    print "     +-l  Listing mode, default is prologues only";
    print "     +-s  Shut-up mode, i.e., ignore any code from BOC to EOC";
    print "     +-x  No LaTeX mode, i.e., put !DESCRIPTION: in verbatim mode";
    print "     +-f  No source file info";
    print "     -A   Ada code";
    print "     -C   C++ code";
    print "     -F   F90 code";
    print "     -S   Shell script";
    print " ";
    print "  The options can appear in any order.  The options, -h and -b,";
    print "  affect the input from all files listed on command-line input.";
    print "  Each of the remaining options effects only the input from the";
    print "  files listed after the option and prior to any overriding";
    print "  option.  The plus sign turns off the option."; 
}# end sub print_help

sub print_notice {

    print "%                **** IMPORTANT NOTICE *****" ;
    print "% This LaTeX file has been automatically produced by ProTeX v. 2.1";
    print "% Any changes made to this file will likely be lost next time";
    print "% this file is regenerated from its source.";
    print " ";

}# sub print_notice

sub print_preamble {

  unless ( $opt_b ) {
    print "%------------------------ PREAMBLE --------------------------";
    print "\\documentclass[11pt]{article}";
    print "\\usepackage{amsmath}";
    print "\\usepackage{epsfig}";
    print "\\usepackage{hangcaption}";
    print "\\textheight     9in";
    print "\\topmargin      0pt";
    print "\\headsep        1cm";
    print "\\headheight     0pt";
    print "\\textwidth      6in";
    print "\\oddsidemargin  0in";
    print "\\evensidemargin 0in";
    print "\\marginparpush  0pt";
    print "\\pagestyle{myheadings}";
    print "\\markboth{}{}";
    print "%-------------------------------------------------------------";
}#end unless

    print "\\parskip        0pt";
    print "\\parindent      0pt";
    print "\\baselineskip  11pt";

}# end sub print_preamble

sub print_macros {

    print " ";
    print "%--------------------- SHORT-HAND MACROS ----------------------";
    print "\\def\\bv{\\begin{verbatim}}";
    print "\\def\\ev\{\\end\{verbatim}}";
    print "\\def\\be{\\begin{equation}}";
    print "\\def\\ee{\\end{equation}}";
    print "\\def\\bea{\\begin{eqnarray}}";
    print "\\def\\eea{\\end{eqnarray}}";
    print "\\def\\bi{\\begin{itemize}}";
    print "\\def\\ei{\\end{itemize}}";
    print "\\def\\bn{\\begin{enumerate}}";
    print "\\def\\en{\\end{enumerate}}";
    print "\\def\\bd{\\begin{description}}";
    print "\\def\\ed{\\end{description}}";
    print "\\def\\({\\left (}";
    print "\\def\\){\\right )}";
    print "\\def\\[{\\left [}";
    print "\\def\\]{\\right ]}";
    print "\\def\\<{\\left  \\langle}";
    print "\\def\\>{\\right \\rangle}";
    print "\\def\\cI{{\\cal I}}";
    print "\\def\\diag{\\mathop{\\rm diag}}";
    print "\\def\\tr{\\mathop{\\rm tr}}";
    print "%-------------------------------------------------------------";

}# end sub print_macros

sub do_beg {
    unless ( $opt_b ) {
    if ( $begdoc == 0 ) {
        if ( $tpage ) {
            print "\\title{$title}";
            print "\\author{{\\sc $author}\\\\ {\\em $affiliation}}";
            print "\\date{$date}";
        }
        print "\\begin{document}";
        if ( $tpage ) {
            print "\\maketitle";
        }
        print "\\tableofcontents";
        print "\\newpage";
        $begdoc = 1;
     }
  }
}# end sub do_beg

sub do_eoc {
        print ' ';
        if ($verb) {
            print '\end{verbatim}';
            $verb = 0;
        }
        $source = 0;
}# end sub do_eoc

sub set_missing {

  $have_name = 0;      # have routine name?
  $have_desc = 0;      # have description?
  $have_intf = 0;      # have interface?
  $have_hist = 0;      # have revision history?
  $name_is = "UNKNOWN";

}# end sub set_missing

    
sub check_if_all_there {

$have_name || 
die "ProTeX: invalid prologue, missing !ROUTINE: or !IROUTINE: in <$name_is>";

$have_desc || 
die "ProTeX: invalid prologue, missing !DESCRIPTION: in <$name_is>";

$have_intf || 
die "ProTeX: invalid prologue, missing !INTERFACE: in <$name_is>";

$have_hist || 
 die "ProTeX: invalid prologue, missing !REVISION HISTORY: in <$name_is>";

}# end sub check_if_all_there
