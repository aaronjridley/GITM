#!/bin/bash
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf
#Most of this code is a modification of code from build_visit2_5_1
function initialize_build_libs()
{

# *************************************************************************** #
#                       Section 1, setting up inputs                          #
# --------------------------------------------------------------------------- #
# This section sets up the inputs to the VisIt script.  This is where you can #
# specify which compiler to use, which versions of the third party libraries, #
# etc.  Note that this script is really only known to work with mpicc.          #
# *************************************************************************** #
export MAKE_OPT_FLAGS=-j4
export BUILD_VISIT_ENV=$(env | cut -d'=' -f1 | sort | uniq)

# Can cause problems in some build systems.
unset CDPATH

# Some systems tar command does not support the deep directory hierarchies
# used in Qt, such as AIX. Gnu tar is a good alternative.
### export TAR=/usr/local/bin/tar # Up and Purple
export TAR=tar

export UTILDIR=$(pwd)
export START_DIR=$(pwd)
export OPSYS=${OPSYS:-$(uname -s)}
export PROC=${PROC:-$(uname -p)}
export REL=${REL:-$(uname -r)}
# Determine architecture
if [[ "$OPSYS" == "Darwin" ]]; then
   export ARCH=${ARCH:-"${PROC}-apple-darwin${REL%%.*}"}
#  export VISITARCH=${VISITARCH-${ARCH}}
   export SO_EXT="dylib"
   VER=$(uname -r)
# Check for Panther, because MACOSX_DEPLOYMENT_TARGET will default to 10.1
   if (( ${VER%%.*} < 8 )) ; then
      export MACOSX_DEPLOYMENT_TARGET=10.3
   elif [[ ${VER%%.*} == 8 ]] ; then
      export MACOSX_DEPLOYMENT_TARGET=10.4
   elif [[ ${VER%%.*} == 9 ]] ; then
      export MACOSX_DEPLOYMENT_TARGET=10.5
   elif [[ ${VER%%.*} == 10 ]] ; then
      export MACOSX_DEPLOYMENT_TARGET=10.6
   else
      export MACOSX_DEPLOYMENT_TARGET=10.6
   fi
   export C_COMPILER=${C_COMPILER:-"mpicc"}
   export CXX_COMPILER=${CXX_COMPILER:-"mpic++"}
   export C_OPT_FLAGS=${C_OPT_FLAGS:-"-O2"}
   export CFLAGS=${CFLAGS:-"-fno-common -fexceptions"}
   export CXX_OPT_FLAGS=${CXX_OPT_FLAGS:-"-O2"}
   export CXXFLAGS=${CXXFLAGS:-"-fno-common -fexceptions"}
   export COMPILEF90=$(which mpif90)
   export FCFLAGS=${FCFLAGS:-$CFLAGS}
#    export FCFLAGS=$(make echo_f90_flags)
elif [[ "$OPSYS" == "Linux" ]]; then
   export COMPILEF90=mpif90
   export ARCH=${ARCH:-"linux-$(uname -m)"} # You can change this to say RHEL, SuSE, Fedora.
   export SO_EXT="so"
   if [[ "$(uname -m)" == "x86_64" ]] ; then
      CFLAGS="$CFLAGS -m64 -fPIC"
      FCFLAGS="$FCFLAGS -m64 -fPIC"
      if [[ "$C_COMPILER" == "mpicc" || "$C_COMPILER" == "" ]]; then
          C_OPT_FLAGS="$C_OPT_FLAGS -O2"
      fi
      CXXFLAGS="$CXXFLAGS -m64 -fPIC"
      if [[ "$CXX_COMPILER" == "mpic++" || "$CXX_COMPILER" == "" ]]; then
          CXX_OPT_FLAGS="$CXX_OPT_FLAGS -O2"
      fi
   elif [[ "$(uname -m)" == "ppc64" ]] ; then
      if [[ "$C_COMPILER" == "xlc" ]] ; then
          CFLAGS="$CFLAGS -qpic"
          FCFLAGS="$FCFLAGS -qpic"
          CXXFLAGS="$CXXFLAGS -qpic"
          export CXX_COMPILER=${CXX_COMPILER-"xlC"}
          export MESA_TARGET=${MESA_TARGET-"linux"}
          QT_PLATFORM="linux-xlc" #aix-xlc"
      else
          CFLAGS="$CFLAGS -fPIC"
          FCFLAGS="$FCFLAGS -fPIC"
          if [[ "$C_COMPILER" == "mpicc" || "$C_COMPILER" == "" ]]; then
              C_OPT_FLAGS="$C_OPT_FLAGS -O2"
          fi
          CXXFLAGS="$CXXFLAGS -fPIC"
          if [[ "$CXX_COMPILER" == "mpic++" || "$CXX_COMPILER" == "" ]]; then
              CXX_OPT_FLAGS="$CXX_OPT_FLAGS -O2"
          fi
      fi
   elif [[ "$(uname -m)" == "ia64" ]] ; then
      CFLAGS="$CFLAGS -fPIC"
      FCFLAGS="$FCFLAGS -fPIC"
      if [[ "$C_COMPILER" == "mpicc" || "$C_COMPILER" == "" ]]; then
          C_OPT_FLAGS="$C_OPT_FLAGS -O2"
      fi
      CXXFLAGS="$CXXFLAGS -fPIC"
      if [[ "$CXX_COMPILER" == "mpic++" || "$CXX_COMPILER" == "" ]]; then
          CXX_OPT_FLAGS="$CXX_OPT_FLAGS -O2"
      fi
      QT_PLATFORM="linux-mpic++"
   fi
   export C_COMPILER=${C_COMPILER:-"mpicc"}
   export CXX_COMPILER=${CXX_COMPILER:-"mpic++"}
   export C_OPT_FLAGS=${C_OPT_FLAGS:-"-O2"}
   export CXX_OPT_FLAGS=${CXX_OPT_FLAGS:-"-O2"}
elif [[ "$OPSYS" == "AIX" ]]; then
   export ARCH="aix" # You can change this to say RHEL, SuSE, Fedora, etc.
   export SO_EXT="a"
   export C_COMPILER=${C_COMPILER:-"xlc"}
   export CXX_COMPILER=${CXX_COMPILER:-"xlC"}
   export C_OPT_FLAGS=${C_OPT_FLAGS:-"-O2"}
   export CXX_OPT_FLAGS=${CXX_OPT_FLAGS:-"-O2"}
   export MAKE=${MAKE:-"gmake"}
elif [[ "$OPSYS" == "IRIX64" ]]; then
   export ARCH="irix64" # You can change this to say RHEL, SuSE, Fedora, etc.
   export SO_EXT="so"
   export C_COMPILER=${C_COMPILER:-"mpicc"}
   export COMPILEF90=$(make echo_compilef90)
   export COMPILEF77=$(make echo_compilef77)
   export CXX_COMPILER=${CXX_COMPILER:-"mpic++"}
   export C_OPT_FLAGS=${C_OPT_FLAGS:-"-O2"}
   export CXX_OPT_FLAGS=${CXX_OPT_FLAGS:-"-O2"}
   export MAKE=${MAKE:-"gmake"}
   export MESA_TARGET=${MESA_TARGET:-"irix6-64-dso"}
elif [[ "$OPSYS" == "SunOS" ]]; then
   export ARCH=${ARCH:-"sunos5"}
   export SO_EXT="so"
   export C_COMPILER=${C_COMPILER:-"mpicc"}
   export CXX_COMPILER=${CXX_COMPILER:-"mpic++"}
   export C_OPT_FLAGS=${C_OPT_FLAGS:-"-O2"}
   export CXX_OPT_FLAGS=${CXX_OPT_FLAGS:-"-O2"}
   export MAKE=${MAKE:-"make"}
   export MESA_TARGET=${MESA_TARGET:-"sunos5-mpicc"}
   export QT_PLATFORM="solaris-mpic++"
else
   export ARCH=${ARCH:-"linux-$(uname -m)"} # You can change this to say RHEL, SuSE, Fedora.
   export SO_EXT="so"
   if [[ "$(uname -m)" == "x86_64" ]] ; then
      CFLAGS="$CFLAGS -m64 -fPIC"
      FCFLAGS="$FCFLAGS -m64 -fPIC"
      if [[ "$C_COMPILER" == "mpicc" || "$C_COMPILER" == "" ]]; then
          C_OPT_FLAGS="$C_OPT_FLAGS -O2"
      fi
      CXXFLAGS="$CXXFLAGS -m64 -fPIC"
      if [[ "$CXX_COMPILER" == "mpic++" || "$CXX_COMPILER" == "" ]]; then
          CXX_OPT_FLAGS="$CXX_OPT_FLAGS -O2"
      fi
      QT_PLATFORM="linux-mpic++-64"
   fi
   if [[ "$(uname -m)" == "ia64" ]] ; then
      CFLAGS="$CFLAGS -fPIC"
      FCFLAGS="$FCFLAGS -fPIC"
      if [[ "$C_COMPILER" == "mpicc" || "$C_COMPILER" == "" ]]; then
          C_OPT_FLAGS="$C_OPT_FLAGS -O2"
      fi
      CXXFLAGS="$CXXFLAGS -fPIC"
      if [[ "$CXX_COMPILER" == "mpic++" || "$CXX_COMPILER" == "" ]]; then
          CXX_OPT_FLAGS="$CXX_OPT_FLAGS -O2"
      fi
      QT_PLATFORM="linux-mpic++-64"
   fi
   export C_COMPILER=${C_COMPILER:-"mpicc"}
   export CXX_COMPILER=${CXX_COMPILER:-"mpic++"}
   export C_OPT_FLAGS=${C_OPT_FLAGS:-"-O2"}
   export CXX_OPT_FLAGS=${CXX_OPT_FLAGS:-"-O2"}
fi
export MAKE=${MAKE:-"make"}
export THIRD_PARTY_PATH=${THIRD_PARTY_PATH:-"./visit"}
export GROUP=${GROUP:-"visit"}
#export LOG_FILE=${LOG_FILE:-"${0##*/}_log"}
export SVNREVISION=${SVNREVISION:-"HEAD"}
# Created a temporary value because the user can override most of 
# the components, which for the GUI happens at a later time.
# the tmp value is useful for user feedback.
if [[ $VISITARCH == "" ]] ; then
    export VISITARCHTMP=${ARCH}_${C_COMPILER}
    if [[ "$CXX_COMPILER" == "mpic++" ]] ; then
        VERSION=$(mpic++ -v 2>&1 | grep "mpicc version" | cut -d' ' -f3 | cut -d'.' -f1-2)
        if [[ ${#VERSION} == 3 ]] ; then
            VISITARCHTMP=${VISITARCHTMP}-${VERSION}
        fi
    fi
else
# use environment variable value
    export VISITARCHTMP=$VISITARCH
fi

REDIRECT_ACTIVE="no"
ANY_ERRORS="no"

#initialize VisIt
# bv_visit_initialize

export DO_DEBUG="no"
export ON_DEBUG="off"
parallel="yes"
ON_parallel="on"
export DO_FORTRAN="yes"
export ON_FORTRAN="on"
verify="no"
ON_verify="off"
export DO_STATIC_BUILD="no"
export USE_VISIBILITY_HIDDEN="no"
export VISIT_INSTALL_PREFIX=""
export VISIT_BUILD_MODE="Release"
DOWNLOAD_ONLY="no"

if [[ "$CXX_COMPILER" == "mpic++" ]] ; then
    VERSION=$(mpic++ -v 2>&1 | grep "mpicc version" | cut -d' ' -f3 | cut -d'.' -f1-1)
    if [[ ${VERSION} -ge 4 ]] ; then
        export USE_VISIBILITY_HIDDEN="yes"
    fi
fi
}

# *************************************************************************** #
# Function: uncompress_untar
#                                                                             #
# Purpose: Uncompress and untar the file, checking if GNU tar can be used.    #
#                                                                             #
# Programmer: Thomas R. Treadway                                              #
# Date: Tue May 15 16:48:01 PDT 2007                                          #
#                                                                             #
# *************************************************************************** #

function uncompress_untar
{
    # Check if GNU tar
    if [[ $(echo $1 | egrep "\.gz$" ) != "" ]] ; then
        COMPRESSTYPE="gzip"
    elif [[ $(echo $1 | egrep "\.bz2$" ) != "" ]] ; then
        COMPRESSTYPE="bzip"
    elif [[ $(echo $1 | egrep "\.tgz$" ) != "" ]] ; then
        COMPRESSTYPE="targzip"
    else
        echo "unsupported uncompression method"
        return 1
    fi
    TARVERSION=$($TAR --version >/dev/null 2>&1)
    if [[ $? == 0 ]] ; then
        case $COMPRESSTYPE in
            gzip|targzip) $TAR zxf $1;;
            bzip) $TAR jxf $1;;
        esac
    else
        case $COMPRESSTYPE in
            gzip) 
               gunzip $1
               $TAR xf ${1%.gz}
               ;;
            targzip) 
               gunzip $1
               $TAR xf "${1%.tgz}.tar"
               ;;
            bzip)
               bunzip2 $1
               $TAR xf ${1%.bz2}
               ;;
        esac
    fi
}



# *************************************************************************** #
# Function: prepare_build_dir                                                 #
#                                                                             #
# Purpose: Helper that prepares a build directory from a src file.            #
#                                                                             #
# Returns:                                                                    #
#          -1 on failure                                                      #
#           0 for success without untar                                       #
#           1 for success with untar                                          #
#           2 for failure with checksum                                       #
#                                                                             #
# Programmer: Cyrus Harrison                                                  #
# Date: Thu Nov 13 09:28:26 PST 2008                                          #
#                                                                             #
# *************************************************************************** #
function prepare_build_dir
{
    BUILD_DIR=$1
    SRC_FILE=$2
    
    #optional
    CHECKSUM_TYPE=$3
    CHECKSUM_VALUE=$4

    untarred_src=0
    if [[ -d ${BUILD_DIR} ]] ; then
       echo "Found ${BUILD_DIR} . . ."
       untarred_src=0
    elif [[ -f ${SRC_FILE} ]] ; then
       if [[ $CHECKSUM != "" && $CHECKSUM_TYPE != "" ]]; then
            verify_checksum $CHECKSUM_TYPE $CHECKSUM ${SRC_FILE}
            if [[ $? != 0 ]]; then
                return 2
            fi
       fi
       echo "Unzipping/Untarring ${SRC_FILE} . . ."
       uncompress_untar ${SRC_FILE}
       untarred_src=1
       if [[ $? != 0 ]] ; then
          echo \
"Unable to untar $SRC_FILE  Corrupted file or out of space on device?"
          return -1
       fi
    elif [[ -f ${SRC_FILE%.*} ]] ; then
       echo "Untarring ${SRC_FILE%.*} . . ."
       $TAR xf ${SRC_FILE%.*}
       untarred_src=1
       if [[ $? != 0 ]] ; then
          echo \
"Unable to untar ${SRC_FILE%.*}.  Corrupted file or out of space on device?"
          return -1
       fi
    fi
    
    return $untarred_src
}


 
# *************************************************************************** 
# Function: download_file
#
# Purpose: DONT USE THIS FUNCTION. USE download_file.
# Downloads a file using wget or curl. 
#
# Programmer: Refactored from other sources (Mark C. Miller)
# Creation: February 19, 2009
#
#   Cyrus Harrison, Tue 24 Mar 13:44:31 PST 2009 
#   As an extra guard, check that the downloaded file actually exisits. 
#   (Firewalls can cause strange files to be created.)
# 
#   Cyrus Harrison, Thu Apr  9 19:21:13 PDT 2009
#   Applied patch from Rick Wagner to fix curl downloads on OSX.
#
#   Tom Fogal, Sun Jul 26 17:19:26 MDT 2009
#   Follow redirects.  Don't use a second argument.
#
#   Gunther H. Weber, Fri Oct 23 13:17:34 PDT 2009
#   Specify explicit path to system curl so that we do not use another
#   version without SSL support
#
# *************************************************************************** 

function download_file
{
    if [[ "$OPSYS" == "Darwin" ]]; then
        # MaxOS X comes with curl
        /usr/bin/curl -ksfLO $1
    else
        check_wget
        if [[ $? != 0 ]] ; then
            echo "Need to download $1, but \
                   cannot locate the wget utility to do so."
        fi
        wget $WGET_OPTS -o /dev/null $1
    fi

    if [[ $? == 0 && -e `basename $1` ]] ; then
        echo "Download succeeded: $1"
        return 0
    else    
        echo "Download attempt failed: $1"
        rm -f `basename $1`
        return 1
    fi
}
function bv_szip_info
{
export SZIP_FILE=${SZIP_FILE:-"szip-2.1.tar.gz"}
export SZIP_VERSION=${SZIP_VERSION:-"2.1"}
export SZIP_COMPATIBILITY_VERSION=${SZIP_COMPATIBILITY_VERSION:-"2.0"}
export SZIP_BUILD_DIR=${SZIP_BUILD_DIR:-"szip-2.1"}
export SZIP_URL=${SZIP_URL:-"http://www.hdfgroup.org/ftp/lib-external/szip/${SZIP_VERSION}/src"}
export SZIP_MD5_CHECKSUM="9cc9125a58b905a4148e4e2fda3fabc6"
export SZIP_SHA256_CHECKSUM=""
export SZIPINSTALLDIR=${START_DIR}/szip-${SZIP_VERSION}
export 
cd ${UTILDIR}
if [ -f ${SZIPINSTALLDIR}/${SZIP_FILE} ]; then
    echo 'mv' ${SZIPINSTALLDIR}'/'${SZIP_FILE} ${UTILDIR}
    mv ${SZIPINSTALLDIR}/${SZIP_FILE} ${UTILDIR}
elif [ -f ${START_DIR}/${SZIP_FILE} ]; then
    echo 'mv' ${START_DIR}/${SZIP_FILE} ${UTILDIR}
    mv ${START_DIR}/${SZIP_FILE} ${UTILDIR}
elif [ -f ${UTILDIR}/${SZIP_FILE} ]; then
    echo ${SZIP_FILE} 'in' ${UTILDIR}
else
    echo 'downloading szip'
    download_file ${SZIP_URL}/${SZIP_FILE} ${SZIP_FILE}
fi

if [ -d ${SZIPINSTALLDIR} ]; then
    echo 'rm -rf' ${SZIPINSTALLDIR}
    rm -rf ${SZIPINSTALLDIR}
fi
}

# *************************************************************************** #
#                          Function 8.0, build_szip                           #
# *************************************************************************** #

function build_szip
{
    #
    # Prepare build dir
    #
#     download_file $site/$dfile $dfile
#     mkdir ${SZIPINSTALLDIR}

    prepare_build_dir ${SZIPINSTALLDIR} ${SZIP_FILE}
    untarred_szip=$?
    if [[ $untarred_szip == -1 ]] ; then
       echo "Unable to prepare SZip build directory. Giving Up!"
       return 1
    fi

    #
    echo "Configuring SZIP . . ."
#     mv ${SZIP_BUILD_DIR} ${UTILDIR}
    cd ${SZIPINSTALLDIR} || error "Can't cd to szip build dir."
    echo "Invoking command to configure SZIP"
    cf_szip=""
    if [[ "$DO_STATIC_BUILD" == "yes" ]]; then
        cf_szip="--disable-shared --enable-static"
    fi

    echo ${SZIPINSTALLDIR}
    echo "./configure CXX=\"$CXX_COMPILER\" CC=\"$C_COMPILER\" LIBS=\"-lm\" \
        CFLAGS=\"$CFLAGS $C_OPT_FLAGS\" CXXFLAGS=\"$CXXFLAGS $CXX_OPT_FLAGS\" \
        --prefix=${SZIPINSTALLDIR} ${cf_szip}"

    ./configure CXX="$CXX_COMPILER" CC="$C_COMPILER" LIBS="-lm" \
        CFLAGS="$CFLAGS $C_OPT_FLAGS" CXXFLAGS="$CXXFLAGS $CXX_OPT_FLAGS" \
        --prefix=${SZIPINSTALLDIR} ${cf_szip}

    if [[ $? != 0 ]] ; then
       echo "SZIP configure failed.  Giving up"
       return 1
    fi

    #
    # Build SZIP
    #
    echo "Building SZIP . . . (~1 minutes)"

    $MAKE
    if [[ $? != 0 ]] ; then
       echo "SZIP build failed.  Giving up"
       return 1
    fi
    #
    # Install into the VisIt third party location.
# :noh
#     #
    echo "Installing SZIP . . ."

    $MAKE install
    if [[ $? != 0 ]] ; then
       echo "SZIP install failed.  Giving up"
       return 1
    fi

    if [[ "$DO_STATIC_BUILD" == "no" && "$OPSYS" == "Darwin" ]]; then
        #
        # Make dynamic executable, need to patch up the install path and
        # version information.
        #
        echo "Creating dynamic libraries for SZIP . . ."

## go back to mpicc bacause if "external relocation entries" restFP saveFP
##      /usr/bin/libtool -o libsz.${SO_EXT} -dynamic src/.libs/libsz.a \
##      -lSystem -lz -headerpad_max_install_names \
##      -install_name $SZIPINSTALLDIR/libsz.${SO_EXT} \
##      -compatibility_version $SZIP_COMPATIBILITY_VERSION \
##      -current_version $SZIP_VERSION
        $C_COMPILER -dynamiclib -o libsz.${SO_EXT} src/*.o \
           -Wl,-headerpad_max_install_names \
           -Wl,-twolevel_namespace,-undefined,dynamic_lookup \
           -Wl,-install_name,${SZIPINSTALLDIR}/libsz.${SO_EXT} \
           -Wl,-compatibility_version,$SZIP_COMPATIBILITY_VERSION \
           -Wl,-current_version,$SZIP_VERSION -lSystem
        if [[ $? != 0 ]] ; then
           echo "SZIP dynamic library build failed.  Giving up"
           return 1
        fi

    fi
    
    cd "$START_DIR"
        echo "Done with SZIP"
        mv ${UTILDIR}/${SZIP_FILE} ${SZIPINSTALLDIR}
    return 0
}

function bv_szip_build
{
cd "$START_DIR"
if [[ "$DO_SZIP" == "yes" ]] ; then
    check_if_installed "szip" $SZIP_VERSION
#     if [[ $? == 0 ]] ; then
        echo "Skipping SZIP build.  SZIP is already installed."
#     else
        echo "Building SZIP (~2 minutes)"
        build_szip
        if [[ $? != 0 ]] ; then
            error "Unable to build or install SZIP.  Bailing out."
        fi
#         echo "Done building SZIP"
#     fi
fi
}


function bv_hdf5_info
{
export HDF5_VERSION=${HDF5_VERSION:-"1.8.8"}
export HDF5_FILE=${HDF5_FILE:-"hdf5-${HDF5_VERSION}.tar.gz"}
export HDF5_COMPATIBILITY_VERSION=${HDF5_COMPATIBILITY_VERSION:-"1.8"}
export HDF5_BUILD_DIR=${HDF5_BUILD_DIR:-"hdf5-${HDF5_VERSION}"}
# Note: Versions of HDF5 1.6.5 and earlier DO NOT have last path component
export HDF5_URL=${HDF5_URL:-"http://www.hdfgroup.org/ftp/HDF5/prev-releases/hdf5-${HDF5_VERSION}/src"}
export HDF5_MD5_CHECKSUM="d4ed6892b17e45a59f998b58ade9987d"
export HDF5_SHA256_CHECKSUM=""
export HDF5INSTALLDIR=${START_DIR}/hdf5-${HDF5_VERSION}
# download_file ${HDF5_URL}/${HDF5_FILE} ${HDF5_FILE}
cd ${UTILDIR}
if [ -f ${HDF5INSTALLDIR}/${HDF5_FILE} ]; then
    echo 'mv' ${HDF5INSTALLDIR}'/'${HDF5_FILE} ${UTILDIR}
    mv ${HDF5INSTALLDIR}/${HDF5_FILE} ${UTILDIR}
elif [ -f ${START_DIR}/${HDF5_FILE} ]; then
    echo 'mv' ${START_DIR}/${HDF5_FILE} ${UTILDIR}
    mv ${START_DIR}/${HDF5_FILE} ${UTILDIR}
elif [ -f ${UTILDIR}/${HDF5_FILE} ]; then
    echo ${HDF5_FILE} 'in' ${UTILDIR}
else
    echo 'downloading HDF5'
    download_file ${HDF5_URL}/${HDF5_FILE} ${HDF5_FILE}
fi

if [ -d ${HDF5INSTALLDIR} ]; then
    echo 'rm -rf' ${HDF5INSTALLDIR}
    rm -rf ${HDF5INSTALLDIR}
fi

}

# *************************************************************************** #
#                          Function 8.1, build_hdf5                           #
# *************************************************************************** #

function build_hdf5
{
    #
    # Prepare build dir
    #
    prepare_build_dir $HDF5_BUILD_DIR $HDF5_FILE
    untarred_hdf5=$?
    if [[ $untarred_hdf5 == -1 ]] ; then
       echo "Unable to prepare HDF5 Build Directory. Giving Up"
       return 1
    fi

    #
    echo "Configuring HDF5 . . ."
    cd $HDF5_BUILD_DIR || error "Can't cd to HDF5 build dir."

    # Fix too many continuation lines in fortran/src/H5f90global.f90
    perl -pi -e 's/,\s*\&/\n  INTEGER(HID_T)  \&/ if /H5T_STD_U32LE/ and $.==109' \
	fortran/src/H5f90global.f90

    cf_darwin=""
    if [[ "$OPSYS" == "Darwin" ]]; then
        export DYLD_LIBRARY_PATH=${SZIPINSTALLDIR}/lib:$DYLD_LIBRARY_PATH
    else
        export LD_LIBRARY_PATH=${SZIPINSTALLDIR}/lib:$LD_LIBRARY_PATH
    fi
    if [[ "$DO_STATIC_BUILD" == "yes" ]]; then
            cf_build_type="--disable-shared --enable-static"
        else
            cf_build_type="--enable-shared --disable-static"
    fi
    cf_szip=""
    if test "x${DO_SZIP}" = "xyes"; then
        echo "SZip requested.  Configuring HDF5 with SZip support."
        cf_szip="--with-szlib=${SZIPINSTALLDIR}"
    fi

    # In order to ensure $FORTRANARGS is expanded to build the arguments to
    # configure, we wrap the invokation in 'sh -c "..."' syntax
    echo "Invoking command to configure HDF5"
    # HDF5 is not supported on OSX but it works that is the reason for the --enable-unsupported flag 
    echo "./configure CC="${C_COMPILER} " FC="${COMPILEF90} " --enable-parallel --enable-fortran --enable-unsupported --prefix="${HDF5INSTALLDIR} ${cf_szip}
        ./configure CC=${C_COMPILER} FC=${COMPILEF90} --enable-parallel --enable-fortran --enable-unsupported --prefix=${HDF5INSTALLDIR} ${cf_szip}
    if [[ $? != 0 ]] ; then
       echo "HDF5 configure failed.  Giving up"
       return 1
    fi

    #
    # Build HDF5
    #
    echo "Making HDF5 . . ."
    $MAKE $MAKE_OPT_FLAGS
    if [[ $? != 0 ]] ; then
       echo "HDF5 build failed.  Giving up"
       return 1
    fi
    #
    # Install into the VisIt third party location.
    #
    echo "Installing HDF5 . . ."

    $MAKE install
    if [[ $? != 0 ]] ; then
       echo "HDF5 install failed.  Giving up"
       return 1
    fi

    # comment out arbitrary Fortran compiler flags
    perl -pi -e 's/^H5BLD_FFLAGS/#H5BLD_FFLAGS/' bin/h5pfc

    cd "$START_DIR"

    echo "Done with HDF5"
    return 0
}

function bv_hdf5_build
{
cd "$START_DIR"

if [[ "$DO_HDF5" == "yes" ]] ; then
    check_if_installed "hdf5" $HDF5_VERSION
    if [[ $? == 0 ]] ; then
        echo "Skipping HDF5 build.  HDF5 is already installed."
    else
        echo "Building HDF5 (~15 minutes)"
        build_hdf5
        if [[ $? != 0 ]] ; then
            error "Unable to build or install HDF5.  Bailing out."
        fi
        echo "Done building HDF5"
    fi
fi
}
usage="install_hdf5.sh [-h] -- Script to install hdf5 library in current directory

-h  show this help text
no flags = install hdf5 (and szip) here

1. The proper mpif90 command used for the SWMF should be in the path before running install_hdf5.sh.
2. The install_hdf5.sh script should be run from the directory where you want hdf5 to be installed.
3. Add the hdf5-1.8.8/bin directory to the path (or make links to the executables that will be in the path).
4. Enable HDF5 in the SWMF: Config.pl -hdf5
5. For the NAG compiler only: edit Makefile.conf in the main directory 
   of the installed SWMF/BATSRUS and change the definition of SEARCH to

   SEARCH = -I\${SHAREDIR} \${SEARCH_EXTRA} -I/YOURHDF5PATH/fortran/src 

6. Enjoy"
while getopts ':hs:' option; do
case "$option" in
h) echo "$usage"
exit
;;
?) printf "illegal option: '%s'\n" "$OPTARG" >&2
echo "$usage" >&2
exit 1
;;
esac
done
initialize_build_libs
bv_szip_info
build_szip
bv_hdf5_info
build_hdf5
