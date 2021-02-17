;  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
;  For more information, see http://csem.engin.umich.edu/tools/swmf
function ask, what, orig_ans, set_orig = set_orig

  if n_elements(orig_ans) eq 0 then orig_ans = ''

  nAnswers = n_elements(orig_ans)

  if (nAnswers eq 1) then begin

      answer = ''
      read, 'Enter '+what+' ['+orig_ans+'] : ', answer
      if strlen(answer) eq 0 then answer = orig_ans
      if n_elements(set_orig) gt 0 then orig_ans = answer

  endif else begin

      answer = strarr(nAnswers)

      TempAnswer = ''
      for i = 0, nAnswers-1 do begin
          read, 'Enter '+what+' '+tostr(i)+' ['+orig_ans(i)+'] : ', TempAnswer
          if strlen(TempAnswer) eq 0 then TempAnswer = orig_ans(i)
          Answer(i) = TempAnswer
      endfor

  endelse

  return, answer

  end
