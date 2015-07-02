;  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
;  For more information, see http://csem.engin.umich.edu/tools/swmf
; Purpose:
;
; Reading FITS files into IDL.
;
; This file contains all functions and procedures needed for
; the basic features of the readfits function.
;
; Basic usage:
;
; Array = read_fits('filename.fits' [, header])
;
; where the array is set to the data contained in the fits file, 
; while the optional header parameter returns a string array 
; with the header information.

function gettok,st,char
;+
; NAME:
;   GETTOK
; PURPOSE:
;   Retrieve the first part of the string up to a specified character
; EXPLANATION:
;   GET TOKen - Retrieve first part of string until the character char
;   is encountered.
;
; CALLING SEQUENCE:
;   token = gettok( st, char )
;
; INPUT:
;   char - character separating tokens, scalar string
;
; INPUT-OUTPUT:
;   st - (scalar) string to get token from (on output token is removed)
;
; OUTPUT:
;   token - scalar string value is returned
;
; EXAMPLE:
;   If ST is 'abc=999' then gettok(ST,'=') would return
;   'abc' and ST would be left as '999'
;
; NOTES:
;       A version of GETTOK that accepts vector strings is available for users
;       of IDL V5.3 or later from  http://idlastro.gsfc.nasa.gov/ftp/v53/
; HISTORY
;   version 1  by D. Lindler APR,86
;   Remove leading blanks    W. Landsman (from JKF)    Aug. 1991
;   Converted to IDL V5.0   W. Landsman   September 1997
;-
;----------------------------------------------------------------------
  On_error,2                           ;Return to caller

; if char is a blank treat tabs as blanks

  tab = string(9b)
  while strpos(st,tab) GE 0 do begin    ;Search for tabs
    pos = strpos(st,tab)
    strput,st,' ',pos
  endwhile

  st = strtrim(st,1)              ;Remove leading blanks

; find character in string

  pos = strpos(st,char)
  if pos EQ -1 then begin         ;char not found?
    token = st
    st = ''
    return, token
  endif

; extract token

 token = strmid(st,0,pos)
 len = strlen(st)
 if pos EQ (len-1) then st = '' else st = strmid(st,pos+1,len-pos-1)

;  Return the result.

 return,token
 end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
; NAME:
;     VALID_NUM
; PURPOSE:
;     Check if a string is a valid number representation.
; EXPLANATION:
;     The input string is parsed for characters that may possibly
;     form a valid number.  It is more robust than simply checking
;     for an IDL conversion error because that allows strings such
;     as '22.3qwert' to be returned as the valid number 22.3
;     See also the original NUM_CHK which returns the status in
;     the opposite sense.
;
; CALLING SEQUENCE:
;     IDL> status = valid_num(string  [,value]  [,/integer])
;
; Inputs      : string  -  the string to be tested
;
; Opt. Inputs : None
;
; Outputs     : The function returns 1 for valid, 0 for invalid number
;
; Opt. Outputs: value   - The value the string decodes to.  This will be
;        returned as a double precision number unless /INTEGER
;        is present, in which case a long integer is returned.
;
; Keywords    : Integer   -  if present code checks specfically for an integer.
;
; Calls       : None
;
; Restrictions: None
;
; Category    : Utilities, Numerical
;
; Prev. Hist. : Small changes from NUM_CHK by Andrew Bowen,
;                                             Tessella Support Services, 8/3/93
;
; Written     : CDS version by C D Pike, RAL, 24-May-93
;
; Modified    : Version 1, C D Pike, RAL, 24-May-93
;     Version 2, William Thompson, GSFC, 14 October 1994
;      Added optional output parameter VALUE to allow
;      VALID_NUM to replace STRNUMBER in FITS routines.
;
; Version     : Version 1  24-May-93
;   Converted to IDL V5.0   W. Landsman   September 1997
;-

FUNCTION valid_num, string, value, INTEGER=integer

       ;**** Set defaults for keyword ****
  IF NOT (KEYWORD_SET(integer)) THEN integer=0

       ;**** arrays of legal characters ****
  numbers   = '0123456789'
  signs     = '+-'
  decimal   = '.'
  exponents     = 'ED'

       ;**** trim leading and trailing blanks/compress white ****
       ;**** space and convert any exponents to uppercase.   ****
  numstr = strupcase(strtrim(strcompress(string),2))

       ;**** length of input string ****
  len = strlen(numstr)

  ok = 1

  if integer eq 0 then stage = 1 else stage = 6

  for i = 0, len-1 do begin

    char = strmid(numstr,i,1)

       ;**** the parsing steps 1 to 8 are for floating   ****
       ;**** point, steps 6 to 8, which test for a legal ****
       ;**** exponent, can be used to check for integers ****

;**** The parsing structure is as follows.  Each character in the ****
;**** string is checked against the valid list at the current     ****
;**** stage.  If no match is found an error is reported.  When a  ****
;**** match is found the stage number is updated as indicated     ****
;**** ready for the next character.  The valid end points are     ****
;**** indicated in the diagram.             ****
;
;Stage  1       2     3   4
;
;Valid  sign --> 2   dec-pt    --> 3  digit    --> 5  dec-pt   --> 5
;  "    dec-pt --> 3   digit --> 4      digit   --> 4
;  "    digit  --> 4              exp't  --> 6
;  "                   END
;
;Stage  5       6     7   8
;
;Valid  digit    --> 5  sign --> 7   digit --> 8   digit -->8
;  "    exp't  --> 6    digit  --> 8         END
;  "    END
;

    CASE stage OF

      1 : begin
        if      strpos(signs,char) ge 0    then stage = 2 $
    else if    decimal eq char    then stage = 3 $
    else if    strpos(numbers,char) ge 0     then stage = 4 $
    else    ok = 0
      end

      2 : begin
    if     decimal eq char    then stage = 3 $
    else if    strpos(numbers,char) ge 0     then stage = 4 $
    else    ok = 0
      end

      3 : begin
    if     strpos(numbers,char) ge 0     then stage = 5 $
    else    ok = 0
      end

      4 : begin
    if     decimal eq char    then stage = 5 $
    else if    strpos(numbers,char) ge 0     then stage = 4 $
    else if       strpos(exponents,char) ge 0   then stage = 6 $
    else    ok = 0
      end

      5 : begin
    if     strpos(numbers,char) ge 0     then stage = 5 $
    else if       strpos(exponents,char) ge 0   then stage = 6 $
    else    ok = 0
      end

      6 : begin
        if      strpos(signs,char) ge 0    then stage = 7 $
    else if    strpos(numbers,char) ge 0     then stage = 8 $
    else    ok = 0
      end

      7 : begin
    if     strpos(numbers,char) ge 0     then stage = 8 $
    else    ok = 0
      end

      8 : begin
    if     strpos(numbers,char) ge 0     then stage = 8 $
    else    ok = 0
      end

    ENDCASE

  end

       ;**** check that the string terminated legally ****
       ;**** i.e in stages 4, 5 or 8                  ****
  if (stage ne 4) and (stage ne 5) and (stage ne 8) then ok = 0

       ;**** If requested, then form the value. ****

  if (n_params() eq 2) and ok then begin
    if keyword_set(integer) then value = long(string) else $
       value = double(string)
  endif

       ;**** return error status to the caller ****
  RETURN, ok


END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function SXPAR, hdr, name, abort, COUNT=matches, COMMENT = comments, $
                                  NoContinue = NoContinue, SILENT = silent
;+
; NAME:
;      SXPAR
; PURPOSE:
;      Obtain the value of a parameter in a FITS header
;
; CALLING SEQUENCE:
;      result = SXPAR( Hdr, Name, [ Abort, COUNT=, COMMENT =, /NoCONTINUE  ])   
;
; INPUTS:
;      Hdr =  FITS header array, (e.g. as returned by READFITS) 
;             string array, each element should have a length of 80 characters      
;
;      Name = String name of the parameter to return.   If Name is of the
;             form 'keyword*' then an array is returned containing values of
;             keywordN where N is an integer.  The value of keywordN will be
;             placed in RESULT(N-1).  The data type of RESULT will be the
;             type of the first valid match of keywordN found.
;
; OPTIONAL INPUTS:
;       ABORT - string specifying that SXPAR should do a RETALL
;               if a parameter is not found.  ABORT should contain
;               a string to be printed if the keyword parameter is not found.
;               If not supplied, SXPAR will return quietly with COUNT = 0
;               (and !ERR = -1) if a keyword is not found.
;
; OPTIONAL INPUT KEYWORDS: 
;       /NOCONTINUE = If set, then continuation lines will not be read, even
;                 if present in the header
;       /SILENT - Set this keyword to suppress warning messages about duplicate
;                 keywords in the FITS header.
;
; OPTIONAL OUTPUT KEYWORDS:
;       COUNT - Optional keyword to return a value equal to the number of 
;               parameters found by SXPAR, integer scalar
;
;       COMMENT - Array of comments associated with the returned values
;
; OUTPUTS:
;       Function value = value of parameter in header.
;               If parameter is double precision, floating, long or string,
;               the result is of that type.  Apostrophes are stripped
;               from strings.  If the parameter is logical, 1b is
;               returned for T, and 0b is returned for F.
;               If Name was of form 'keyword*' then a vector of values
;               are returned.
;
; SIDE EFFECTS:
;       !ERR is set to -1 if parameter not found, 0 for a scalar
;       value returned.  If a vector is returned it is set to the
;       number of keyword matches found.    The use of !ERR is deprecated, and
;       instead the COUNT keyword is preferred
;
;       If a keyword (except HISTORY or COMMENT) occurs more than once in a 
;       header, a warning is given, and the *last* occurence is used.
;
; EXAMPLES:
;       Given a FITS header, h, return the values of all the NAXISi values
;       into a vector.    Then place the history records into a string vector.
;
;       IDL> naxisi = sxpar( h ,'NAXIS*')         ; Extract NAXISi value
;       IDL> history = sxpar( h, 'HISTORY' )      ; Extract HISTORY records
;
; PROCEDURE:
;       The first 8 chacters of each element of Hdr are searched for a 
;       match to Name.  The value from the last 20 characters is returned.  
;       An error occurs if there is no parameter with the given name.
;
;       If a numeric value has no decimal point it is returned as type
;       LONG.   If it contains more than 8 numerals, or contains the 
;       characters 'D' or 'E', then it is returned as type DOUBLE.  Otherwise
;       it is returned as type FLOAT.    Very large integer values, outside
;       the range of valid LONG, are returned as DOUBLE.
;
;       If the value is too long for one line, it may be continued on to the
;       the next input card, using the OGIP CONTINUE convention.  For more info,
;       http://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/ofwg_recomm/r13.html
;
;       Complex numbers are recognized as two numbers separated by one or more
;       space characters.
;
;       If a numeric value has no decimal point (or E or D) it is returned as
;       type LONG.  If it contains more than 8 numerals, or contains the
;       character 'D', then it is returned as type DOUBLE.  Otherwise it is
;       returned as type FLOAT.    If an integer is too large to be stored as
;       type LONG, then it is returned as DOUBLE.
;
; NOTES:
;       The functions SXPAR() and FXPAR() are nearly identical, although
;       FXPAR() has slightly more sophisticated parsing.   There is no
;       particular reason for having two nearly identical procedures, but
;       both are too widely used to drop either one.
;
; PROCEDURES CALLED:
;       GETTOK(), VALID_NUM()
; MODIFICATION HISTORY:
;       DMS, May, 1983, STPAR Written.
;       D. Lindler Jan 90 added ABORT input parameter
;       J. Isensee Jul,90 added COUNT keyword
;       W. Thompson, Feb. 1992, added support for FITS complex values.
;       W. Thompson, May 1992, corrected problem with HISTORY/COMMENT/blank
;               keywords, and complex value error correction.
;       W. Landsman, November 1994, fix case where NAME is an empty string 
;       W. Landsman, March 1995,  Added COMMENT keyword, ability to read
;               values longer than 20 character
;       W. Landsman, July 1995, Removed /NOZERO from MAKE_ARRAY call
;       T. Beck May 1998, Return logical as type BYTE
;       W. Landsman May 1998, Make sure integer values are within range of LONG
;       Converted to IDL V5.0, May 1998
;       W. Landsman Feb 1998, Recognize CONTINUE convention 
;       W. Landsman Oct 1999, Recognize numbers such as 1E-10 as floating point
;       W. Landsman Jan 2000, Only accept integer N values when name = keywordN
;       W. Landsman Dec 2001, Optional /SILENT keyword to suppress warnings
;       W. Landsman/D. Finkbeiner  Mar 2002  Make sure extracted vectors 
;             of mixed data type are returned with the highest type.
;-
;----------------------------------------------------------------------
 if N_params() LT 2 then begin
     print,'Syntax -     result =  sxpar( hdr, name, [abort])'
     print,'   Input Keywords:    /NOCONTINUE, /SILENT'
     print,'   Output Keywords:   COUNT=,  COMMENT= '
     return, -1
 endif 

 VALUE = 0
 if N_params() LE 2 then begin
      abort_return = 0
      abort = 'FITS Header'
 end else abort_return = 1
 if abort_return then On_error,1 else On_error,2

;       Check for valid header

  s = size(hdr)         ;Check header for proper attributes.
  if ( s[0] NE 1 ) or ( s[2] NE 7 ) then $
           message,'FITS Header (first parameter) must be a string array'

  nam = strtrim( strupcase(name) )      ;Copy name, make upper case     


;  Determine if NAME is of form 'keyword*'.  If so, then strip off the '*', and
;  set the VECTOR flag.  One must consider the possibility that NAM is an empty
;  string.

   namelength1 = (strlen(nam) - 1 ) > 1         
   if strpos( nam, '*' ) EQ namelength1 then begin    
            nam = strmid( nam, 0, namelength1)  
            vector = 1                  ;Flag for vector output  
            name_length = strlen(nam)   ;Length of name 
            num_length = 8 - name_length        ;Max length of number portion  
            if num_length LE 0 then  $ 
                  message, 'Keyword length must be 8 characters or less'

;  Otherwise, extend NAME with blanks to eight characters.

    endif else begin  
                while strlen(nam) LT 8 do nam = nam + ' ' ;Make 8 chars long
                vector = 0      
    endelse


;  If of the form 'keyword*', then find all instances of 'keyword' followed by
;  a number.  Store the positions of the located keywords in NFOUND, and the
;  value of the number field in NUMBER.

        histnam = (nam eq 'HISTORY ') or (nam eq 'COMMENT ') or (nam eq '') 
        if N_elements(start) EQ 0 then start = -1l
        start = long(start[0])
        if (not vector) and (start GE 0) then begin
            if N_elements(precheck)  EQ 0 then precheck = 5
            if N_elements(postcheck) EQ 0 then postcheck = 20
            nheader = N_elements(hdr)
            mn = (start - precheck)  > 0
            mx = (start + postcheck) < nheader-1
            keyword = strmid(hdr[mn:mx], 0, 8)
        endif else begin
            restart:
            start   = -1l
            keyword = strmid( hdr, 0, 8)
        endelse

        if vector then begin
            nfound = where(strpos(keyword,nam) GE 0, matches)
            if ( matches gt 0 ) then begin
                numst= strmid( hdr[nfound], name_length, num_length)
                number = replicate(-1, matches)
                for i = 0, matches-1 do         $
                    if VALID_NUM( numst[i], num,/INTEGER) then number[i] = num
                igood = where(number GE 0, matches)
                if matches GT 0 then begin
                    nfound = nfound[igood]
                    number = number[igood]
                endif
            endif

;  Otherwise, find all the instances of the requested keyword.  If more than
;  one is found, and NAME is not one of the special cases, then print an error
;  message.

        endif else begin
            nfound = where(keyword EQ nam, matches)
            if (matches EQ 0) and (start GE 0) then goto, RESTART
            if (start GE 0) then nfound = nfound + mn
            if (matches GT 1) and (not histnam) then        $
                if not keyword_set(silent) then $
                message,/informational, 'Warning - keyword ' +   $
                nam + ' located more than once in ' + abort
            if (matches GT 0) then start = nfound[matches-1]
        endelse


; Process string parameter 

 if matches GT 0 then begin
  line = hdr[nfound]
  svalue = strtrim( strmid(line,9,71),2)
  if histnam then $
        value = strtrim(strmid(line,8,71),2) else for i = 0,matches-1 do begin
      if ( strmid(svalue[i],0,1) EQ "'" ) then begin   ;Is it a string?
                  test = strmid( svalue[i],1,strlen( svalue[i] )-1)
                  next_char = 0
                  off = 0
                  value = '' 
          NEXT_APOST:
                  endap = strpos(test, "'", next_char)      ;Ending apostrophe  
                  if endap LT 0 then $ 
                            MESSAGE,'Value of '+name+' invalid in '+abort
                  value = value + strmid( test, next_char, endap-next_char )  

;  Test to see if the next character is also an apostrophe.  If so, then the
;  string isn't completed yet.  Apostrophes in the text string are signalled as
;  two apostrophes in a row.

                 if strmid( test, endap+1, 1) EQ "'" then begin    
                    value = value + "'"
                    next_char = endap+2         
                    goto, NEXT_APOST
                 endif      

; Extract the comment, if any
                
                slash = strpos( test, "/", endap )
                if slash LT 0 then comment = '' else    $
                        comment = strmid( test, slash+1, strlen(test)-slash-1 )

; This is a string that could be continued on the next line.  Check this
; possibility with the following four criteria: *1) Ends with '&'
; (2) Next line is CONTINUE  (3) LONGSTRN keyword is present (recursive call to
; SXPAR) 4. /NOCONTINE is not set

    if not keyword_set(nocontinue) then begin
                off = off + 1
                val = strtrim(value,2)

                if (strlen(val) gt 0) and $
                  (strmid(val, strlen(val)-1, 1) EQ '&') and $
                  (strmid(hdr[nfound[i]+off],0,8) EQ 'CONTINUE') then begin
                   if (size(sxpar(hdr, 'LONGSTRN',/NoCONTINUE)))[1] EQ 7 then begin                    
                  value = strmid(val, 0, strlen(val)-1)
                  test = hdr[nfound[i]+off]
                  test = strmid(test, 8, strlen(test)-8)
                  test = strtrim(test, 2)
                  if strmid(test, 0, 1) NE "'" then message, $
                    'ERROR: Invalidly CONTINUEd string in '+ abort
                  next_char = 1
                  GOTO, NEXT_APOST
                ENDIF
               ENDIF
    ENDIF


; Process non-string value  

          endif else begin

                test = svalue[i]
                slash = strpos( test, "/" )
                if slash GT 0 then begin
                        comment = strmid( test, slash+1, strlen(test)-slash-1 )
                        test = strmid( test, 0, slash )
                end else comment = ''

; Find the first word in TEST.  Is it a logical value ('T' or 'F')

                test2 = test
                value = gettok(test2,' ')
               if ( value EQ 'T' ) then value = 1b else $
               if ( value EQ 'F' ) then value = 0b else begin

;  Test to see if a complex number.  It's  a complex number if the value and
;  the next word, if any, are both valid values.

                if strlen(test2) EQ 0 then goto, NOT_COMPLEX
                value2 = gettok( test2, ' ') 
                if value2 EQ '' then goto, NOT_COMPLEX
                On_ioerror, NOT_COMPLEX
                value2 = float(value2)
                value = complex(value,value2)
                goto, GOT_VALUE

;  Not a complex number.  Decide if it is a floating point, double precision,
;  or integer number.

NOT_COMPLEX:
                On_IOerror, GOT_VALUE
                  if (strpos(value,'.') GE 0) or (strpos(value,'E') GT 0) $
                  or (strpos(value,'D') GE 0) then begin  ;Floating or double?
                      if ( strpos(value,'D') GT 0 ) or $  ;Double?
                         ( strlen(value) GE 8 ) then value = double(value) $
                                                else value = float(value)
                       endif else begin                   ;Long integer
                            lmax = 2.0d^31 - 1.0d
                            lmin = -2.0d31
                            value = double(value)
                            if (value GE lmin) and (value LE lmax) then $
                                value = long(value)
                       endelse

GOT_VALUE:
                On_IOerror, NULL
                endelse
             endelse; if c eq apost

;  Add to vector if required

         if vector then begin
               if ( i EQ 0 ) then begin
                     maxnum = max(number)
                     dtype = size(value,/type)
                     result = make_array( maxnum, TYPE = dtype )
                     comments = strarr( maxnum )
               endif 
               if size(value,/type) GT dtype then begin   ;Do we need to recast?
                    result = result + 0*value
                    dtype = size(value,/type)
               endif
               result[ number[i]-1 ] =  value
               comments[ number[i]-1 ] = comment
          endif else $
                comments = comment
  endfor

  if vector then begin
         !ERR = matches     
         return, result
  endif else !ERR = 0

endif  else  begin    
     if abort_return then message,'Keyword '+nam+' not found in '+abort
     !ERR = -1
endelse     

return, value       

END                 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sxaddpar, Header, Name, Value, Comment, Location, before=before, $
                 savecomment = savecom, after=after , format=format, pdu = pdu
;+
; NAME:
;       SXADDPAR
; PURPOSE:
;       Add or modify a parameter in a FITS header array.
;
; CALLING SEQUENCE:
;       SXADDPAR, Header, Name, Value, [ Comment,  Location, /SaveComment, 
;                               BEFORE =, AFTER = , FORMAT= , /PDU]
;
; INPUTS:
;       Header = String array containing FITS or STSDAS header.    The
;               length of each element must be 80 characters.    If not 
;               defined, then SXADDPAR will create an empty FITS header array.
;
;       Name = Name of parameter. If Name is already in the header the value 
;               and possibly comment fields are modified.  Otherwise a new 
;               record is added to the header.  If name is equal to 'COMMENT'
;               or 'HISTORY' or a blank string then the value will be added to 
;               the record without replacement.  For these cases, the comment 
;               parameter is ignored.
;
;       Value = Value for parameter.  The value expression must be of the 
;               correct type, e.g. integer, floating or string.  String values
;                of 'T' or 'F' are considered logical values.
;
; OPTIONAL INPUT PARAMETERS:
;       Comment = String field.  The '/' is added by this routine.  Added 
;               starting in position 31.    If not supplied, or set equal to 
;               '', or /SAVECOMMENT is set, then the previous comment field is 
;               retained (when found) 
;
;       Location = Keyword string name.  The parameter will be placed before the
;               location of this keyword.    This parameter is identical to
;               the BEFORE keyword and is kept only for consistency with
;               earlier versions of SXADDPAR.
;
; OPTIONAL INPUT KEYWORD PARAMETERS:
;       BEFORE  = Keyword string name.  The parameter will be placed before the
;               location of this keyword.  For example, if BEFORE='HISTORY'
;               then the parameter will be placed before the first history
;               location.  This applies only when adding a new keyword;
;               keywords already in the header are kept in the same position.
;
;       AFTER   = Same as BEFORE, but the parameter will be placed after the
;               location of this keyword.  This keyword takes precedence over
;               BEFORE.
;
;       FORMAT  = Specifies FORTRAN-like format for parameter, e.g. "F7.3".  A
;               scalar string should be used.  For complex numbers the format
;               should be defined so that it can be applied separately to the
;               real and imaginary parts.  If not supplied then the default is
;               'G19.12' for double precision, and 'G14.7' for floating point.
;
;       /PDU    = specifies keyword is to be added to the primary data unit
;               header. If it already exists, it's current value is updated in
;               the current position and it is not moved.
;       /SAVECOMMENT = if set, then any existing comment is retained, i.e. the
;               COMMENT parameter only has effect if the keyword did not 
;               previously exist in the header.
; OUTPUTS:
;       Header = updated FITS header array.
;
; EXAMPLE:
;       Add a keyword 'TELESCOP' with the value 'KPNO-4m' and comment 'Name
;       of Telescope' to an existing FITS header h.
;
;       IDL> sxaddpar, h, 'TELESCOPE','KPNO-4m','Name of Telescope'
; NOTES:
;       The functions SXADDPAR() and FXADDPAR() are nearly identical, with the
;       major difference being that FXADDPAR forces required FITS keywords
;       BITPIX, NAXISi, EXTEND, PCOUNT, GCOUNT to appear in the required order
;       in the header, and FXADDPAR supports the OGIP LongString convention.   
;       There is no particular reason for having two nearly identical 
;       procedures, but both are too widely used to drop either one.
;
;       All HISTORY records are inserted in order at the end of the header.
;
;       All COMMENT records are also inserted in order at the end of the header
;       header, but before the HISTORY records.  The BEFORE and AFTER keywords
;       can override this.
;
;       All records with no keyword (blank) are inserted in order at the end of
;       the header, but before the COMMENT and HISTORY records.  The BEFORE and
;       AFTER keywords can override this.

; RESTRICTIONS:
;       Warning -- Parameters and names are not checked
;               against valid FITS parameter names, values and types.
;
; MODIFICATION HISTORY:
;       DMS, RSI, July, 1983.
;       D. Lindler Oct. 86  Added longer string value capability
;       Converted to NEWIDL  D. Lindler April 90
;       Added Format keyword, J. Isensee, July, 1990
;       Added keywords BEFORE and AFTER. K. Venkatakrishna, May '92
;       Pad string values to at least 8 characters   W. Landsman  April 94
;       Aug 95: added /PDU option and changed routine to update last occurence
;               of an existing keyword (the one SXPAR reads) instead of the
;               first occurence.
;       Comment for string data can start after column 32 W. Landsman June 97
;       Make sure closing quote supplied with string value  W. Landsman  June 98
;       Converted to IDL V5.0    W. Landsman   June 98
;       Increase precision of default formatting of double precision floating
;               point values.   C. Gehman, JPL  September 1998
;       Mar 2000, D. Lindler, Modified to use capital E instead of lower case
;               e for exponential formats.
;       Apr 2000, Make user-supplied format upper-case  W. Landsman 
;       Oct 2001, Treat COMMENT or blank string like HISTORY keyword W. Landsman
;       Jan 2002, Allow BEFORE, AFTER to apply to COMMENT keywords W. Landsman
;       June 2003, Added SAVECOMMENT keyword    W. Landsman
;       Jan 2004, If END is missing, then add it at the end W. Landsman
;       May 2005 Fix SAVECOMMENT error with non-string values W. Landsman
;       Oct 2005 Jan 2004 change made SXADDPAR fail for empty strings W.L. 
;       
;-
 if N_params() LT 3 then begin             ;Need at least 3 parameters
      print,'Syntax - Sxaddpar, Header, Name,  Value, [Comment, Postion'
      print,'                      BEFORE = ,AFTER = , FORMAT =, /SAVECOMMENT]'
      return
 endif

; Define a blank line and the END line

 ENDLINE = 'END' +string(replicate(32b,77))     ;END line
 BLANK = string(replicate(32b,80))             ;BLANK line
;
;  If Location parameter not defined, set it equal to 'END     '
;
 if ( N_params() GT 4 ) then loc = strupcase(location) else $
 if keyword_set( BEFORE) then loc = strupcase(before) else $
 if keyword_set( AFTER)  then loc = strupcase(after) else $
 if keyword_set( PDU) then loc = 'BEGIN EX' else $
                             loc = 'END'

 while strlen(loc) lt 8 do loc = loc + ' '

 if N_params() lt 4 then comment = ''      ;Is comment field specified?

 n = N_elements(header)                  ;# of lines in FITS header
 if (n EQ 0) then begin                  ;header defined?
          header=strarr(10)              ;no, make it.
          header[0]=ENDLINE
          n=10
 endif else begin
          s = size(header)               ;check for string type
              if (s[0] ne 1) or (s[2] ne 7) then $
                  message,'FITS Header (first parameter) must be a string array'
 endelse

;  Make sure Name is 8 characters long

        nn = string(replicate(32b,8))   ;8 char name
        strput,nn,strupcase(name) ;insert name

;  Extract first 8 characters of each line of header, and locate END line

 keywrd = strmid(header,0,8)                 ;Header keywords
 iend = where(keywrd eq 'END     ',nfound)
;
;  If no END, then add it.  Either put it after the last non-null string, or
;  append it to the end.
;
        if nfound EQ 0 then begin
                ii = where(strtrim(header) ne '',nfound)
                ii = max(ii) + 1
                if ii eq n_elements(header) then begin
                        header = [header,endline]
                        n = n+1 
                endif else header[ii] = endline
                keywrd = strmid(header,0,8)
                iend = where(keywrd eq 'END     ',nfound)
        endif
;
        iend = iend[0] > 0                      ;make scalar

;  History, comment and "blank" records are treated differently from the
;  others.  They are simply added to the header array whether there are any
;  already there or not.

 if (nn EQ 'HISTORY ') or (nn EQ 'COMMENT ') or $
    (nn EQ '        ')  then begin             ;add history record?
;
;  If the header array needs to grow, then expand it in increments of 5 lines.
;

     if iend GE (n-1) then begin
                 header = [header,replicate(blank,5)] ;yes, add 5.
                 n = N_elements(header)
      endif

; Format the record

      newline = blank
      strput,newline,nn+string(value),0

;
;  If a history record, then append to the record just before the end.
;
      if nn EQ 'HISTORY ' then begin
             header[iend] = newline             ;add history rec.
             header[iend+1] = endline
;
;  The comment record is placed immediately after the last previous comment
;  record, or immediately before the first history record, unless overridden by
;  either the BEFORE or AFTER keywords.
;
      endif else if nn EQ 'COMMENT ' then begin
            if loc EQ 'END     ' then loc = 'COMMENT '
            iloc = where(keywrd EQ loc, nloc)
            if nloc EQ 0 then iloc = where(keywrd EQ 'HISTORY ', nloc)
            if nloc gt 0 then begin
               i = iloc[nloc-1]
               if keyword_set(after) or (loc EQ 'COMMENT ') then i = i+1 < iend 
               if i gt 0 then header=[header[0:i-1],newline,header[i:n-1]] $
                        else header=[newline,header[0:n-1]]
            endif else begin
                header[iend] = newline
                header[iend+1] = endline
            endelse

;
;  The "blank" record is placed immediately after the last previous "blank"
;  record, or immediately before the first comment or history record, unless
;  overridden by either the BEFORE or AFTER keywords.
;
          ENDIF ELSE BEGIN
            if loc EQ 'END     ' then loc = '       '
            iloc = where(keywrd[0:iend] EQ loc, nloc)
            if nloc gt 0 then begin
               i = iloc[0]
               if keyword_set(after) and loc ne 'HISTORY ' then i = i+1 < iend 
               if i gt 0 then header=[header[0:i-1],newline,header[i:n-1]] $
                        else header=[newline,header[0:n-1]]
            endif else begin
                iloc = where(keywrd EQ 'COMMENT ', nloc)
                if nloc Eq 0 then iloc = where(keywrd EQ 'HISTORY ', nloc)
                if nloc GT 0 then begin
                   i = iloc[0]
                   if i gt 0 then header=[header[0:i-1],newline,header[i:n-1]] $
                        else header=[newline,header[0:n-1]]
                endif else begin
                  header[iend] = newline
                  header[iend+1] = endline
            endelse
            endelse
           endelse
            RETURN
 endif

; Find location to insert keyword.   Save the existing comment if user did
; not supply a new one.   Comment starts after column 32 for numeric data,
; after the slash (but at least after column 20) for string data. 

 ncomment = comment
 ipos  = where(keywrd eq nn,nfound)
 if nfound gt 0 then begin
         i = ipos[nfound-1]
         if comment eq '' or keyword_set(savecom) then begin  ;save comment?
         if strmid(header[i],10,1) NE "'" then $
                 ncomment=strmid(header[i],32,48) else begin
                 slash = strpos(header[i],'/', 20)  
                 if slash NE -1 then $
                        ncomment =  strmid(header[i], slash+1, 80) else $
                        ncomment = string(replicate(32B,80))
                endelse
        endif 
         goto, REPLACE    
 endif

 if loc ne '' then begin
          iloc =  where(keywrd eq loc,nloc)
          if nloc gt 0 then begin
             i = iloc[0]
             if keyword_set(after) and loc ne 'HISTORY ' then i = i+1 < iend 
             if i gt 0 then header=[header[0:i-1],blank,header[i:n-1]] $
                        else header=[blank,header[0:n-1]]
             goto, REPLACE  
          endif
 endif

; At this point keyword and location parameters were not found, so a new
; line is added at the end of the FITS header

        if iend lt (n-1) then begin     ;Not found, add more?
                header[iend+1] = ENDLINE        ;no, already long enough.
                i = iend                ;position to add.
           endif else begin             ;must lengthen.
                header = [header,replicate(blank,5)] ;add an element on the end
                header[n]=ENDLINE               ;save "END"
                i =n-1                  ;add to end
        end

; Now put value into keyword at line i

REPLACE:    
        h=blank                 ;80 blanks
        strput,h,nn+'= '        ;insert name and =.
        apost = "'"             ;quote a quote
        type = size(value)      ;get type of value parameter
        if type[0] ne 0 then $
                message,'Keyword Value (third parameter) must be scalar'

        case type[1] of         ;which type?

7:      begin
          upval = strupcase(value)      ;force upper case.
          if (upval eq 'T') or (upval eq 'F') then begin
                strput,h,upval,29  ;insert logical value.
            end else begin              ;other string?
                if strlen(value) gt 18 then begin       ;long string
                    strput, h, apost + strmid(value,0,68) + apost + $
                        ' /' + ncomment,10
                    header[i] = h
                    return
                endif
                strput, h, apost + value,10       ;insert string val
                strput, h, apost, 11 + (strlen(value)>8)   ;pad string vals
          endelse                                          ;to at least 8 chars
          endcase

5:      BEGIN
        IF (N_ELEMENTS(format) EQ 1) THEN $             ; use format keyword
            v = string(value, FORMAT='('+strupcase(format)+')') $
        ELSE v = STRING(value, FORMAT='(G19.12)')
        s = strlen(v)                                   ; right justify
        strput, h, v, (30-s)>10
        END

 else:  begin
        if (N_elements(format) eq 1) then $            ;use format keyword
            v = string(value, FORMAT='('+strupcase(format)+')' ) else $
            v = strtrim(strupcase(value),2)      
                                      ;convert to string, default format
        s = strlen(v)                 ;right justify
        strput,h,v,(30-s)>10          ;insert
        end
 endcase

 strput,h,' /',30       ;add ' /'
 strput, h, ncomment, 32 ;add comment
 header[i] = h          ;save line

 return
 end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
; NAME:
;       READFITS
; PURPOSE:
;       Read a FITS file into IDL data and header variables. 
; EXPLANATION:
;       READFITS() can also read gzip or Unix compressed FITS files.
;       See http://idlastro.gsfc.nasa.gov/fitsio.html for other ways of
;       reading FITS files with IDL.   
;
; CALLING SEQUENCE:
;       Result = READFITS( Filename/Fileunit,[ Header, heap, /NOSCALE, EXTEN_NO=,
;                     NSLICE=, /SILENT , STARTROW =, NUMROW = , HBUFFER=,
;                     /CHECKSUM, /COMPRESS, /No_Unsigned, NaNVALUE = ]
;
; INPUTS:
;       Filename = Scalar string containing the name of the FITS file  
;                 (including extension) to be read.   If the filename has
;                  a *.gz extension, it will be treated as a gzip compressed
;                  file.   If it has a .Z extension, it will be treated as a
;                  Unix compressed file.
;                                   OR
;       Fileunit - A scalar integer specifying the unit of an already opened
;                  FITS file.  The unit will remain open after exiting 
;                  READFITS().  There are two possible reasons for choosing 
;                  to specify a unit number rather than a file name:
;          (1) For a FITS file with many extensions, one can move to the 
;              desired extensions with FXPOSIT() and then use READFITS().  This
;              is more efficient than repeatedly starting at the beginning of 
;              the file.
;          (2) For reading a FITS file across a Web http: address after opening
;              the unit with the SOCKET procedure (IDL V5.4 or later,
;              Unix and Windows only) 
;
; OUTPUTS:
;       Result = FITS data array constructed from designated record.
;                If the specified file was not found, then Result = -1
;
; OPTIONAL OUTPUT:
;       Header = String array containing the header from the FITS file.
;              If you don't need the header, then the speed may be improved by
;              not supplying this parameter.    Note however, that omitting 
;              the header can imply /NOSCALE, i.e. BSCALE and BZERO values
;              may not be applied.
;       heap = For extensions, the optional heap area following the main
;              data array (e.g. for variable length binary extensions).
;
; OPTIONAL INPUT KEYWORDS:
;       /CHECKSUM - If set, then READFITS() will call FITS_TEST_CHECKSUM to 
;                verify the data integrity if CHECKSUM keywords are present
;                in the FITS header.   Cannot be used with the NSLICE, NUMROW
;                or STARTROW keywords, since verifying the checksum requires 
;               that all the data be read.  See FITS_TEST_CHECKSUM() for more
;               information.
;
;       /COMPRESS - Signal that the file is gzip compressed.  By default, 
;               READFITS will assume that if the file name extension ends in 
;               '.gz' then the file is gzip compressed.   The /COMPRESS keyword
;               is required only if the the gzip compressed file name does not 
;               end in '.gz'
;              
;
;       EXTEN_NO - non-negative scalar integer specifying the FITS extension to
;               read.  For example, specify EXTEN = 1 or /EXTEN to read the 
;               first FITS extension.   
;   
;        HBUFFER - Number of lines in the header, set this to slightly larger
;                than the expected number of lines in the FITS header, to 
;               improve performance when reading very large FITS headers. 
;               Should be a multiple of 36 -- otherwise it will be modified
;               to the next higher multiple of 36.   Default is 180
;
;       /NOSCALE - If present and non-zero, then the ouput data will not be
;                scaled using the optional BSCALE and BZERO keywords in the 
;                FITS header.   Default is to scale.
;
;       /NO_UNSIGNED - By default, if the header indicates an unsigned integer 
;               (BITPIX = 16, BZERO=2^15, BSCALE=1) then FITS_READ will output 
;               an IDL unsigned integer data type (UINT).   But if /NO_UNSIGNED
;               is set, then the data is converted to type LONG.  
;
;       NSLICE - An integer scalar specifying which N-1 dimensional slice of a 
;                N-dimensional array to read.   For example, if the primary 
;                image of a file 'wfpc.fits' contains a 800 x 800 x 4 array, 
;                then 
;
;                 IDL> im = readfits('wfpc.fits',h, nslice=2)
;                           is equivalent to 
;                 IDL> im = readfits('wfpc.fits',h)
;                 IDL> im = im[*,*,2]
;                 but the use of the NSLICE keyword is much more efficient.
;
;       NUMROW -  Scalar non-negative integer specifying the number of rows 
;                 of the image or table extension to read.   Useful when one 
;                 does not want to read the entire image or table.   This
;                 keyword is only for extensions and is ignored for primary
;                 arrays.
;
;       POINT_LUN  -  Position (in bytes) in the FITS file at which to start
;                 reading.   Useful if READFITS is called by another procedure
;                 which needs to directly read a FITS extension.    Should 
;                 always be a multiple of 2880, and not be used with EXTEN_NO
;                 keyword.
;
;       /SILENT - Normally, READFITS will display the size the array at the
;                 terminal.  The SILENT keyword will suppress this
;
;        STARTROW - Non-negative integer scalar specifying the row
;               of the image or extension table at which to begin reading. 
;               Useful when one does not want to read the entire table.  This
;               keyword is ignored when reading a primary data array.
;
;       NaNVALUE - This keyword is included only for backwards compatibility
;                  with routines that require IEEE "not a number" values to be
;                  converted to a regular value.
;
; EXAMPLE:
;       Read a FITS file test.fits into an IDL image array, IM and FITS 
;       header array, H.   Do not scale the data with BSCALE and BZERO.
;
;              IDL> im = READFITS( 'test.fits', h, /NOSCALE)
;
;       If the file contains a FITS extension, it could be read with
;
;              IDL> tab = READFITS( 'test.fits', htab, /EXTEN )
;
;       The function TBGET() can be used for further processing of a binary 
;       table, and FTGET() for an ASCII table.
;       To read only rows 100-149 of the FITS extension,
;
;              IDL> tab = READFITS( 'test.fits', htab, /EXTEN, 
;                                   STARTR=100, NUMR = 50 )
;
;       To read in a file that has been compressed:
;
;              IDL> tab = READFITS('test.fits.gz',h)
;
; ERROR HANDLING:
;       If an error is encountered reading the FITS file, then 
;               (1) the system variable !ERROR_STATE.CODE is set negative 
;                   (via the MESSAGE facility)
;               (2) the error message is displayed (unless /SILENT is set),
;                   and the message is also stored in !!ERROR_STATE.MSG
;               (3) READFITS returns with a value of -1
; RESTRICTIONS:
;       (1) Cannot handle random group FITS
;
; NOTES:
;       (1) If data is stored as integer (BITPIX = 16 or 32), and BSCALE
;       and/or BZERO keywords are present, then the output array is scaled to 
;       floating point (unless /NOSCALE is present) using the values of BSCALE
;       and BZERO.   In the header, the values of BSCALE and BZERO are then 
;       reset to 1. and 0., while the original values are written into the 
;       new keywords O_BSCALE and O_BZERO.     If the BLANK keyword was
;       present, then any input integer values equal to BLANK in the input
;       integer image are unchanged by BSCALE or BZERO
;       
;       (2) The use of the NSLICE keyword is incompatible with the NUMROW
;       or STARTROW keywords.
;
;       (3) READFITS() underwent a substantial rewrite in February 2000 to 
;       take advantage of new features in IDL V5.3
;            1. The /swap_if_little_endian keyword is now used to OPENR rather
;                than calling IEEE_TO_HOST for improved performance
;            2. The /compress keyword is now used with OPENR to allow gzip files
;                to be read on any machine architecture.
;            3. Removed NANvalue keyword, since in V5.3, NaN is recognized on
;                all machine architectures
;            4. Assume unsigned integers are always allowed
;            5. Use STRJOIN to display image size
;            6. Use !ERROR_STATE.MSG rather than !ERR_STRING
;      
;
;       (4) On some Unix shells, one may get a "Broken pipe" message if reading
;        a Unix compressed (.Z) file, and not reading to the end of the file 
;       (i.e. the decompression has not gone to completion).     This is an 
;        informative message only, and should not affect the output of READFITS.   
; PROCEDURES USED:
;       Functions:   SXPAR()
;       Procedures:  SXADDPAR, SXDELPAR
; MINIMUM IDL VERSION:
;       V5.3 (Uses STRJOIN, /COMPRESS keyword to OPENR)
;
; MODIFICATION HISTORY:
;       Original Version written in 1988, W.B. Landsman   Raytheon STX
;       Revision History prior to October 1998 removed          
;       Major rewrite to eliminate recursive calls when reading extensions
;                  W.B. Landsman  Raytheon STX                    October 1998
;       Add /binary modifier needed for Windows  W. Landsman    April 1999
;       Read unsigned datatypes, added /no_unsigned   W. Landsman December 1999
;       Output BZERO = 0 for unsigned data types   W. Landsman   January 2000
;       Update to V5.3 (see notes)  W. Landsman                  February 2000
;       Fixed logic error in use of NSLICE keyword  W. Landsman  March 2000
;       Fixed byte swapping for Unix compress files on little endian machines
;                                    W. Landsman    April 2000
;       Added COMPRESS keyword, catch IO errors W. Landsman September 2000
;       Option to read a unit number rather than file name W.L    October 2001
;       Fix undefined variable problem if unit number supplied W.L. August 2002
;       Don't read entire header unless needed   W. Landsman  Jan. 2003
;       Added HBUFFER keyword    W. Landsman   Feb. 2003
;       Added CHECKSUM keyword   W. Landsman   May 2003
;       Restored NaNVALUE keyword for backwards compatibility,
;               William Thompson, 16-Aug-2004, GSFC
;-
function READFITS, filename, header, heap, CHECKSUM=checksum, $ 
                   COMPRESS = compress, HBUFFER=hbuf, EXTEN_NO = exten_no, $
                   NOSCALE = noscale, NSLICE = nslice, $
                   NO_UNSIGNED = no_unsigned,  NUMROW = numrow, $
                   POINTLUN = pointlun, SILENT = silent, STARTROW = startrow, $
                   NaNvalue = NaNvalue

  On_error,2                    ;Return to user
  On_IOerror, BAD

; Check for filename input

   if N_params() LT 1 then begin                
      print,'Syntax - im = READFITS( filename, [ h, heap, /NOSCALE, /SILENT,'
      print,'                 EXTEN_NO =, STARTROW = , NUMROW=, NSLICE = ,'
      print,'                 HBUFFER = ,/NO_UNSIGNED, /CHECKSUM, /COMPRESS]'
      return, -1
   endif

   unitsupplied = size(filename,/TNAME) NE 'STRING'

; Set default keyword values

   silent = keyword_set( SILENT )
   do_checksum = keyword_set( CHECKSUM )
   if N_elements(exten_no) EQ 0 then exten_no = 0

;  Check if this is a Unix compressed file.   (gzip files are handled 
;  separately using the /compress keyword to OPENR).

    unixZ = 0                
    if unitsupplied then unit = filename else begin
    len = strlen(filename)
    gzip = strmid(filename,len-3,3) EQ '.gz'
    compress = keyword_set(compress) or gzip
    unixZ =  (strmid(filename, len-2, 2) EQ '.Z') and $
             (!VERSION.OS_FAMILY EQ 'unix') 

;  Go to the start of the file.

   openr, unit, filename, ERROR=error,/get_lun,/BLOCK,/binary, $
                COMPRESS = compress, /swap_if_little_endian
   if error NE 0 then begin
        message,/con,' ERROR - Unable to locate file ' + filename
        return, -1
   endif

;  Handle Unix compressed files.   On some Unix machines, users might wish to 
;  force use of /bin/sh in the line spawn, ucmprs+filename, unit=unit,/sh

        if unixZ then begin
                free_lun, unit
                spawn, 'uncompress -c '+filename, unit=unit                 
                gzip = 1b
        endif 
  endelse
  if keyword_set(POINTLUN) then begin
       if gzip then  readu,unit,bytarr(pointlun,/nozero) $
               else point_lun, unit, pointlun
  endif
  doheader = arg_present(header) or do_checksum
  if doheader  then begin
          if N_elements(hbuf) EQ 0 then hbuf = 180 else begin
                  remain = hbuf mod 36
                  if remain GT 0 then hbuf = hbuf + 36-remain
           endelse
  endif else hbuf = 36

  for ext = 0L, exten_no do begin
               
;  Read the next header, and get the number of bytes taken up by the data.

       block = string(replicate(32b,80,36))
       w = [-1]
       if ((ext EQ exten_no) and (doheader)) then header = strarr(hbuf) $
                                             else header = strarr(36)
       headerblock = 0L
       i = 0L      

       while w[0] EQ -1 do begin
          
       if EOF(unit) then begin 
            message,/ CON, $
               'EOF encountered attempting to read extension ' + strtrim(ext,2)
            if not unitsupplied then free_lun,unit
            return,-1
       endif

      readu, unit, block
      headerblock = headerblock + 1
      w = where(strlen(block) NE 80, Nbad)
      if (Nbad GT 0) then begin
           message,'Warning-Invalid characters in header',/INF,NoPrint=Silent
           block[w] = string(replicate(32b, 80))
      endif
      w = where(strmid(block, 0, 8) eq 'END     ', Nend)
      if (headerblock EQ 1) or ((ext EQ exten_no) and (doheader)) then begin
              if Nend GT 0 then  begin
             if headerblock EQ 1 then header = block[0:w[0]]   $
                                 else header = [header[0:i-1],block[0:w[0]]]
             endif else begin
                header[i] = block
                i = i+36
                if i mod hbuf EQ 0 then $
                              header = [header,strarr(hbuf)]
           endelse
          endif
      endwhile

      if (ext EQ 0 ) and (keyword_set(pointlun) EQ 0) then $
             if strmid( header[0], 0, 8)  NE 'SIMPLE  ' then begin
              message,/CON, $
                 'ERROR - Header does not contain required SIMPLE keyword'
                if not unitsupplied then free_lun, unit
                return, -1
      endif

                
; Get parameters that determine size of data region.
                
       bitpix =  sxpar(header,'BITPIX')
       naxis  = sxpar(header,'NAXIS')
       gcount = sxpar(header,'GCOUNT') > 1
       pcount = sxpar(header,'PCOUNT')
                
       if naxis GT 0 then begin 
            dims = sxpar( header,'NAXIS*')           ;Read dimensions
            ndata = dims[0]
            if naxis GT 1 then for i = 2, naxis do ndata = ndata*dims[i-1]
                        
                endif else ndata = 0
                
                nbytes = (abs(bitpix) / 8) * gcount * (pcount + ndata)

;  Move to the next extension header in the file.   Although we could use
;  POINT_LUN with compressed files, a faster way is to simply read into the 
;  file

      if ext LT exten_no then begin
                nrec = long((nbytes + 2879) / 2880)
                if nrec GT 0 then begin     
                if gzip then begin 
                        buf = bytarr(nrec*2880L,/nozero)
                        readu,unit,buf 
                        endif else  begin 
                        point_lun, -unit,curr_pos
                        point_lun, unit,curr_pos + nrec*2880L
                endelse
                endif
       endif
       endfor

 case BITPIX of 
           8:   IDL_type = 1          ; Byte
          16:   IDL_type = 2          ; Integer*2
          32:   IDL_type = 3          ; Integer*4
          64:   IDL_type = 14         ; Integer*8
         -32:   IDL_type = 4          ; Real*4
         -64:   IDL_type = 5          ; Real*8
        else:   begin
                message,/CON, 'ERROR - Illegal value of BITPIX (= ' +  $
                strtrim(bitpix,2) + ') in FITS header'
                if not unitsupplied then free_lun,unit
                return, -1
                end
  endcase     

; Check for dummy extension header

 if Naxis GT 0 then begin 
        Nax = sxpar( header, 'NAXIS*' )   ;Read NAXES
        ndata = nax[0]
        if naxis GT 1 then for i = 2, naxis do ndata = ndata*nax[i-1]

  endif else ndata = 0

  nbytes = (abs(bitpix)/8) * gcount * (pcount + ndata)
 
  if nbytes EQ 0 then begin
        if not SILENT then message, $
                "FITS header has NAXIS or NAXISi = 0,  no data array read",/CON
        if do_checksum then begin
             ;;;result = FITS_TEST_CHECKSUM(header, data, ERRMSG = errmsg)
             if not SILENT then begin
               case result of 
                1: message,/INF,'CHECKSUM keyword in header is verified'
               -1: message,/CON, errmsg
                else: 
                endcase
              endif
        endif
        if not unitsupplied then free_lun, unit
        return,-1
 endif

; Check for FITS extensions, GROUPS

 groups = sxpar( header, 'GROUPS' ) 
 if groups then message,NoPrint=Silent, $
           'WARNING - FITS file contains random GROUPS', /INF

; If an extension, did user specify row to start reading, or number of rows
; to read?

   if not keyword_set(STARTROW) then startrow = 0
   if naxis GE 2 then nrow = nax[1] else nrow = ndata
   if not keyword_set(NUMROW) then numrow = nrow
   if do_checksum then if ((startrow GT 0) or $
      (numrow LT nrow) or (N_elements(nslice) GT 0)) then begin 
      message,/CON, $
      'Warning - CHECKSUM not applied when STARTROW, NUMROW or NSLICE is set'
      do_checksum = 0
   endif 

   if exten_no GT 0 then begin
        xtension = strtrim( sxpar( header, 'XTENSION' , Count = N_ext),2)
        if N_ext EQ 0 then message, /INF, NoPRINT = Silent, $
                'WARNING - Header missing XTENSION keyword'
   endif 

   if (exten_no GT 0) and ((startrow NE 0) or (numrow NE nrow)) then begin
        if startrow GE nax[1] then begin
           message,'ERROR - Specified starting row ' + strtrim(startrow,2) + $
          ' but only ' + strtrim(nax[1],2) + ' rows in extension',/CON
           if not unitsupplied then free_lun,unit
           return,-1
        endif 
        nax[1] = nax[1] - startrow    
        nax[1] = nax[1] < numrow
        sxaddpar, header, 'NAXIS2', nax[1]
        if gzip then begin
                if startrow GT 0 then begin
                        tmp=bytarr(startrow*nax[0],/nozero)
                        readu,unit,tmp
                endif 
        endif else begin 
              point_lun, -unit, pointlun          ;Current position
              point_lun, unit, pointlun + startrow*nax[0]
    endelse
    endif else if (N_elements(NSLICE) EQ 1) then begin
        lastdim = nax[naxis-1]
        if nslice GE lastdim then message,/CON, $
        'ERROR - Value of NSLICE must be less than ' + strtrim(lastdim,2)
        nax = nax[0:naxis-2]
        sxdelpar,header,'NAXIS' + strtrim(naxis,2)
        naxis = naxis-1
        sxaddpar,header,'NAXIS',naxis
        ndata = ndata/lastdim
        nskip = nslice*ndata*abs(bitpix/8) 
        if gzip then  begin 
              if Ndata GT 0 then begin
                  buf = bytarr(nskip,/nozero)
                  readu,unit,buf
               endif   
        endif else begin 
                   point_lun, -unit, currpoint          ;Current position
                   point_lun, unit, currpoint + nskip
        endelse
  endif


  if not (SILENT) then begin   ;Print size of array being read

         if exten_no GT 0 then message, $
                     'Reading FITS extension of type ' + xtension, /INF
         st = 'Now reading ' + strjoin(strtrim(NAX,2),' by ') + ' array'
         if (exten_no GT 0) and (pcount GT 0) then st = st + ' + heap area'
         message,/INF,st   
   endif

; Read Data in a single I/O call.   Only need byteswapping for Unix compress
; files

    data = make_array( DIM = nax, TYPE = IDL_type, /NOZERO)
    readu, unit, data
    ;;; if unixZ then if not is_ieee_big() then ieee_to_host,data
    if (exten_no GT 0) and (pcount GT 0) then begin
        theap = sxpar(header,'THEAP')
        skip = theap - N_elements(data)
        if skip GT 0 then begin 
                temp = bytarr(skip,/nozero)
                readu, unit, skip
        endif
        heap = bytarr(pcount*gcount*abs(bitpix)/8)
        readu, unit, heap
        ;;;if do_checksum then $
        ;;;result = fits_test_checksum(header,[data,heap],ERRMSG=errmsg)
    endif ;;; else if do_checksum then $
          ;;; result = fits_test_checksum(header, data, ERRMSG = errmsg)
    if not unitsupplied then free_lun, unit
    if do_checksum then if not SILENT then begin
        case result of 
        1: message,/INF,'CHECKSUM keyword in header is verified'
       -1: message,/CON, 'CHECKSUM ERROR! ' + errmsg
        else: 
        endcase
    endif

; Scale data unless it is an extension, or /NOSCALE is set
; Use "TEMPORARY" function to speed processing.  

   do_scale = not keyword_set( NOSCALE )
   if (do_scale and (exten_no GT 0)) then do_scale = xtension EQ 'IMAGE' 
   if do_scale then begin

          Nblank = 0
          if bitpix GT 0 then begin
                blank = sxpar( header, 'BLANK', Count = N_blank) 
                if N_blank GT 0 then $ 
                        blankval = where( data EQ blank, Nblank)
          endif

          Bscale = float( sxpar( header, 'BSCALE' , Count = N_bscale))
          Bzero = float( sxpar(header, 'BZERO', Count = N_Bzero ))
 
; Check for unsigned integer (BZERO = 2^15) or unsigned long (BZERO = 2^31)

          if not keyword_set(No_Unsigned) then begin
            no_bscale = (Bscale EQ 1) or (N_bscale EQ 0)
            unsgn_int = (bitpix EQ 16) and (Bzero EQ 32768) and no_bscale
            unsgn_lng = (bitpix EQ 32) and (Bzero EQ 2147483648) and no_bscale
            unsgn = unsgn_int or unsgn_lng
           endif else unsgn = 0

          if unsgn then begin
                 sxaddpar, header, 'BZERO', 0
                 sxaddpar, header, 'O_BZERO', bzero, $
                          'Original Data is unsigned Integer'
                   if unsgn_int then $ 
                        data =  uint(data) - 32768U else $
                   if unsgn_lng then  data = ulong(data) - 2147483648UL 
                
          endif else begin
 
          if N_Bscale GT 0  then $ 
               if ( Bscale NE 1. ) then begin
                   data = temporary(data) * Bscale 
                   sxaddpar, header, 'BSCALE', 1.
                   sxaddpar, header, 'O_BSCALE', Bscale,' Original BSCALE Value'
               endif

         if N_Bzero GT 0  then $
               if (Bzero NE 0) then begin
                     data = temporary( data ) + Bzero
                     sxaddpar, header, 'BZERO', 0.
                     sxaddpar, header, 'O_BZERO', Bzero,' Original BZERO Value'
               endif
        
        endelse

        if (Nblank GT 0) and ((N_bscale GT 0) or (N_Bzero GT 0)) then $
                data[blankval] = blank

        endif

; Return array.  If necessary, first convert NaN values.

        if n_elements(nanvalue) eq 1 then begin
            w = where(finite(data,/nan),count)
            if count gt 0 then data(w) = nanvalue
        endif
        return, data    

; Come here if there was an IO_ERROR
    
 BAD:   print,!ERROR_STATE.MSG
        if (not unitsupplied) and (N_elements(unit) GT 0) then free_lun, unit
        if N_elements(data) GT 0 then return,data else return, -1

 end 
