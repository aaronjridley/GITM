!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
subroutine merge_str(str1, str2)

  character (len=100) :: str1, str2, temp
  integer :: i, j, k

  i = 1
  do while (iachar(str1(i:i)) /= 32 .and. &
            iachar(str1(i:i)) /= 9  .and. &
            i < 100) 
     i=i+1
  enddo

  j = 1
  do while (iachar(str2(j:j)) /= 32 .and. &
            iachar(str2(j:j)) /= 9  .and. &
            j < 100) 
     j=j+1
  enddo

  temp = str1
  do k = i,100
     temp(k:k) = ' '
  enddo

  if (i+j-1 > 100) j = 100 - i + 1

  temp(i:i+j-1) = str2(1:j)

  str2 = temp

end subroutine merge_str

subroutine strlen(str1, len)

  character (len=100) :: str1
  integer :: len

  len = 1
  do while (iachar(str1(len:len)) /= 32 .and. &
            iachar(str1(len:len)) /= 9  .and. &
            len < 100) 
     len=len+1
  enddo

end subroutine strlen
