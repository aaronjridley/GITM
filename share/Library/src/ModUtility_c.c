//  Copyright (C) 2002 Regents of the University of Michigan,
// portions used with permission
// For more information, see http://csem.engin.umich.edu/tools/swmf

#include <sys/stat.h>
#include <errno.h>

/** 

@brief Create a directory with permission "perm".

@param path Path to new directory

@return 0 for success, -1 for failure

*/
int make_dir_c(const char *path, int perm, int *mkdir_errno)
{

  mode_t uperm=perm;

  // Make the directory
  int retval = mkdir(path, uperm);

  *mkdir_errno = errno;

  // It is not an error if the directory already exists
  if(errno == EEXIST){
    retval = 0;
  }

  return retval;

}
