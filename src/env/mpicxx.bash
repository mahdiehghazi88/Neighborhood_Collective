#! /cvmfs/soft.computecanada.ca/nix/var/nix/profiles/16.09/bin/bash
#
# Copyright (c) 2001-2014, The Ohio State University. All rights reserved.
#  
# This file is part of the MVAPICH2 software package developed by the team
# members of The Ohio State University's Network-Based Computing Laboratory
# (NBCL), headed by Professor Dhabaleswar K. (DK) Panda.
#
# For detailed copyright and licensing information, please refer to the
# copyright file COPYRIGHT in the top level MVAPICH2 directory.
# 
# (C) 2006 by Argonne National Laboratory.
#     See COPYRIGHT in top-level directory.
#
# mpicxx
# Simple script to compile and/or link MPI programs.
# This script knows the default flags and libraries, and can handle
# alternative C++ compilers and the associated flags and libraries.
# The important terms are:
#    includedir, libdir - Directories containing an *installed* mpich
#    prefix, execprefix - Often used to define includedir and libdir
#    CXX                - C compiler
#    WRAPPER_CXXFLAGS      - Any special flags needed to compile 
#    WRAPPER_LDFLAGS       - Any special flags needed to link
#    WRAPPER_LIBS          - Any special libraries needed in order to link
#    
# We assume that (a) the C++ compiler can both compile and link programs
#
# Handling of command-line options:
#   This is a little tricky because some options may contain blanks.
#
# Special issues with shared libraries - todo
#
# --------------------------------------------------------------------------
# Set the default values of all variables.
#
# Directory locations: Fixed for any MPI implementation
prefix=__PREFIX_TO_BE_FILLED_AT_INSTALL_TIME__
exec_prefix=__EXEC_PREFIX_TO_BE_FILLED_AT_INSTALL_TIME__
sysconfdir=__SYSCONFDIR_TO_BE_FILLED_AT_INSTALL_TIME__
includedir=__INCLUDEDIR_TO_BE_FILLED_AT_INSTALL_TIME__
libdir=__LIBDIR_TO_BE_FILLED_AT_INSTALL_TIME__

# Default settings for compiler, flags, and libraries
CXX="g++"
MVAPICH2_VERSION="2.2"

enable_wrapper_rpath="yes"

# How to pass a linker flag through the compiler.
wl="-Wl,"

# Static library suffix (normally "a").
libext="a"

# Shared library suffix (normally "so").
shlibext="so"

# Format of library name prefix.
libname_spec="lib\$name"

# Library names that the linker finds when passed -lNAME.
library_names_spec="\$libname\$shrext"

# Flag to hardcode $libdir into a binary during linking.
# This must work even if $libdir does not exist.
hardcode_libdir_flag_spec="\${wl}-rpath \${wl}\$libdir -Wl,--enable-new-dtags"

# Whether we need a single -rpath flag with a separated argument.
hardcode_libdir_separator=""

# Set to yes if using DIR/libNAME.so during linking hardcodes DIR into the
# resulting binary.
hardcode_direct="no"

# Set to yes if using the -LDIR flag during linking hardcodes DIR into the
# resulting binary.
hardcode_minus_L="no"


# Internal variables
# Show is set to echo to cause the compilation command to be echoed instead 
# of executed.
Show=
#
# End of initialization of variables
#---------------------------------------------------------------------
# Environment Variables.
# The environment variables MPICH_CXX may be used to override the 
# default choices.
# In addition, if there is a file $sysconfdir/mpicxx-$CXXname.conf, 
# where CXXname is the name of the compiler with all spaces replaced by hyphens
# (e.g., "CC -64" becomes "CC--64", that file is sources, allowing other
# changes to the compilation environment.  See the variables used by the 
# script (defined above)
# Added MPICH_CXX_OLD, MPICH_CXX can be used to prefix CXX
# with external utility, e.g. setenv MPICH_CXX 'eval linkcache $MPICH_CXX_OLD'
if [ -n "$MPICH_CXX" ] ; then
    MPICH_CXX_OLD="$CXX"
    CXX="$MPICH_CXX"
    CXXname=`echo $CXX | sed 's/ /-/g'`
    if [ -s $sysconfdir/mpicxx-$CXXname.conf ] ; then
        . $sysconfdir/mpicxx-$CXXname.conf
    fi
fi
# Allow a profiling option to be selected through an environment variable
if [ -n "$MPICXX_PROFILE" ] ; then
    profConf=$MPICXX_PROFILE
fi
#
# ------------------------------------------------------------------------
# Argument processing.
# This is somewhat awkward because of the handling of arguments within
# the shell.  We want to handle arguments that include spaces without 
# loosing the spacing (an alternative would be to use a more powerful
# scripting language that would allow us to retain the array of values, 
# which the basic (rather than enhanced) Bourne shell does not.  
#
# Look through the arguments for arguments that indicate compile only.
# If these are *not* found, add the library options

linking=yes
allargs=("$@")
argno=0
for arg in "$@" ; do
    # Set addarg to no if this arg should be ignored by the C compiler
    addarg=yes
    case "$arg" in 
 	# ----------------------------------------------------------------
	# Compiler options that affect whether we are linking or no
    -c|-S|-E|-M|-MM)
    # The compiler links by default
    linking=no
    ;;
	# ----------------------------------------------------------------
	# Options that control how we use mpicxx (e.g., -show, 
	# -cxx=* -config=*
    -echo)
    addarg=no
    set -x
    ;;
    -cxx=*)
    CXX=`echo A$arg | sed -e 's/A-cxx=//g'`
    addarg=no
    ;;
    # Backwards compatibility for MPICH1 - scripts
    -CC=*)
    CXX=`echo A$arg | sed -e 's/A-CC=//g'`
    addarg=no
    ;;
    -show)
    addarg=no
    Show=echo
    ;;
    -config=*)
    addarg=no
    CXXname=`echo A$arg | sed -e 's/A-config=//g'`
    if [ -s "$sysconfdir/mpicxx-$CXXname.conf" ] ; then
        . "$sysconfdir/mpicxx-$CXXname.conf"
    else
	echo "Configuration file mpicxx-$CXXname.conf not found"
    fi
    ;;
    -compile-info|-compile_info)
    # -compile_info included for backward compatibility
    Show=echo
    addarg=no
    ;;
    -link-info|-link_info)
    # -link_info included for backward compatibility
    Show=echo
    addarg=no
    ;;
    -v)
    # Pass this argument to the compiler as well.
    echo "mpicxx for MVAPICH2 version $MVAPICH2_VERSION"
    # if there is only 1 argument, it must be -v.
    if [ "$#" -eq "1" ] ; then
        linking=no
    fi
    ;;
    -profile=*)
    # Pass the name of a profiling configuration.  As
    # a special case, lib<name>.so or lib<name>.la may be used
    # if the library is in $libdir
    profConf=`echo A$arg | sed -e 's/A-profile=//g'`
    addarg=no
    # Loading the profConf file is handled below
    ;;
    -nativelinking)
    # Internal option to use native compiler for linking without MPI libraries
    nativelinking=yes
    addarg=no
    ;;
    -itac)
    if [ -n "$VT_ROOT" ] ; then
        if [ -e "$VT_ROOT/bin/itacvars.sh" ] ; then
            VT_LIB_DIR=`bash -c ". $VT_ROOT/bin/itacvars.sh; echo \\$VT_LIB_DIR"`
            VT_ADD_LIBS=`bash -c ". $VT_ROOT/bin/itacvars.sh; echo \\$VT_ADD_LIBS"`
            ITAC_OPTIONS="-L$VT_LIB_DIR -lVT $VT_ADD_LIBS"
        else
            echo "Cannot find itacvars.sh in $VT_ROOT/bin. -itac option ignored."
        fi
    else
        echo "VT_ROOT should be set with the root directory of Intel Trace Analyzer and Collector. -itac option ignored."
    fi
    addarg=no
    ;;
    -help)
    NC=`echo "$CXX" | sed 's%\/% %g' | awk '{print $NF}' -`
    if [ -f "$sysconfdir/mpixxx_opts.conf" ] ; then
        . $sysconfdir/mpixxx_opts.conf
        echo "    -cxx=xxx      - Reset the native compiler to xxx."
    else
        if [ -f "./mpixxx_opts.conf" ] ; then
            . ./mpixxx_opts.conf
            echo "    -cxx=xxx      - Reset the native compiler to xxx."
        fi
    fi
    exit 0
    ;;
    *)
    ;;

    esac
    if [ $addarg = no ] ; then
	unset allargs[$argno]
    fi
    # Some versions of bash do not accept ((argno++))
    argno=`expr $argno + 1`
done

if [ $# -eq 0 ] ; then
    echo "Error: Command line argument is needed!"
    "$0" -help
    exit 1
fi

# -----------------------------------------------------------------------
# Derived variables.  These are assembled from variables set from the
# default, environment, configuration file (if any) and command-line
# options (if any)
cxxlibs=
if [ "mpicxx" != "mpi" ] ; then
    cxxlibs="-lmpicxx"
fi

PROFILE_FOO=
# Handle the case of a profile switch
if [ -n "$profConf" ] ; then
    profConffile=
    if [ -s "$libdir/lib$profConf.a" -o -s "$libdir/lib$profConf.so" ] ; then
	PROFILE_FOO="-l$profConf"
    elif [ -s "$sysconfdir/$profConf.conf" ] ; then
	profConffile="$sysconfdir/$profConf.conf"
    elif [ -s "$profConf.conf" ] ; then
        profConffile="$profConf.conf"
    else
        echo "Profiling configuration file $profConf.conf not found in $sysconfdir"
    fi
    if [ -n "$profConffile" -a -s "$profConffile" ] ; then
	. $profConffile
    fi
fi

# A temporary statement to invoke the compiler
# Place the -L before any args incase there are any mpi libraries in there.
# Eventually, we'll want to move this after any non-MPI implementation 
# libraries
if [ "$linking" = yes ] ; then
    # Attempt to encode rpath info into the executable if the user has not
    # disabled rpath usage and some flavor of rpath makes sense on this
    # platform.
    # TODO configure and config.rpath are computing more sophisticated rpath
    # schemes than this simple one.  Consider updating this logic accordingly.
    if test "X$enable_wrapper_rpath" = "Xyes" ; then
        eval rpath_flags=\"${hardcode_libdir_flag_spec}\"
    else
	rpath_flags=""
    fi

    if [ "$nativelinking" = yes ] ; then
        $Show $CXX $PROFILE_INCPATHS    "${allargs[@]}" -I$includedir
        rc=$?
    else
        $Show $CXX $PROFILE_INCPATHS    "${allargs[@]}" -I$includedir $ITAC_OPTIONS -L$libdir $cxxlibs $PROFILE_PRELIB $PROFILE_FOO $rpath_flags -lmpi  $PROFILE_POSTLIB 
        rc=$?
    fi
else
    $Show $CXX $PROFILE_INCPATHS   "${allargs[@]}" -I$includedir
    rc=$?
fi

exit $rc

