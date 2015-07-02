!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModErrors

  integer, parameter :: ecBeforeStartTime_ = 1
  integer, parameter :: ecAfterEndTime_    = 2
  integer, parameter :: ecTimeNotSet_      = 3
  integer, parameter :: ecAllocationError_ = 4
  integer, parameter :: ecPointOutofRange_ = 5
  integer, parameter :: ecFileNotFound_    = 6
  integer, parameter :: ecIMFBzNotSet_     = 7
  integer, parameter :: ecIMFByNotSet_     = 8
  integer, parameter :: ecIMFBxNotSet_     = 9
  integer, parameter :: ecSWVNotSet_       = 10
  integer, parameter :: ecHPINotSet_       = 11
  integer, parameter :: ecKpNotSet_        = 12
  integer, parameter :: ecEFieldModelNotFound_ = 13
  integer, parameter :: ecAuroralModelNotFound_ = 14
  integer, parameter :: ecSWNNotSet_       = 15

  integer, parameter :: nErrorsMax = 1000

  character (len=100), dimension(nErrorsMax) :: cErrorCodes

contains

  subroutine set_error_codes

    cErrorCodes(ecBeforeStartTime_) = "Time requested is before start time"
    cErrorCodes(ecAfterEndTime_)    = "Time requested is after end time"
    cErrorCodes(ecTimeNotSet_)      = "Time not set"
    cErrorCodes(ecAllocationError_) = "Allocation Error"
    cErrorCodes(ecPointOutofRange_) = "Requested point is out of range"
    cErrorCodes(ecFileNotFound_)    = "File not found"
    cErrorCodes(ecIMFBxNotSet_)     = "IMF Bx has not been set"
    cErrorCodes(ecIMFByNotSet_)     = "IMF By has not been set"
    cErrorCodes(ecIMFBzNotSet_)     = "IMF Bz has not been set"
    cErrorCodes(ecSWVNotSet_)       = "Solar Wind V has not been set"
    cErrorCodes(ecHPINotSet_)       = &
         "Hemispheric Power Index has not been set"
    cErrorCodes(ecKpNotSet_)       = "Kp has not been set"
    cErrorCodes(ecEFieldModelNotFound_) = &
         "The Selected Electric Field model is unknown"
    cErrorCodes(ecAuroralModelNotFound_) = &
         "The Selected Auroral model is unknown"
    cErrorCodes(ecSWNNotSet_)       = "Solar Wind N has not been set"

  end subroutine set_error_codes

end module ModErrors
