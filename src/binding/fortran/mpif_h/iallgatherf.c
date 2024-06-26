/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*  
 *  (C) 2001 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 *
 * This file is automatically generated by buildiface 
 * DO NOT EDIT
 */
#include "mpi_fortimpl.h"


/* Begin MPI profiling block */
#if defined(USE_WEAK_SYMBOLS) && !defined(USE_ONLY_MPI_NAMES) 
#if defined(HAVE_MULTIPLE_PRAGMA_WEAK)
extern FORT_DLL_SPEC void FORT_CALL MPI_IALLGATHER( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint * );
extern FORT_DLL_SPEC void FORT_CALL mpi_iallgather__( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint * );
extern FORT_DLL_SPEC void FORT_CALL mpi_iallgather( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint * );
extern FORT_DLL_SPEC void FORT_CALL mpi_iallgather_( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint * );

#if defined(F77_NAME_UPPER)
#pragma weak MPI_IALLGATHER = PMPI_IALLGATHER
#pragma weak mpi_iallgather__ = PMPI_IALLGATHER
#pragma weak mpi_iallgather_ = PMPI_IALLGATHER
#pragma weak mpi_iallgather = PMPI_IALLGATHER
#elif defined(F77_NAME_LOWER_2USCORE)
#pragma weak MPI_IALLGATHER = pmpi_iallgather__
#pragma weak mpi_iallgather__ = pmpi_iallgather__
#pragma weak mpi_iallgather_ = pmpi_iallgather__
#pragma weak mpi_iallgather = pmpi_iallgather__
#elif defined(F77_NAME_LOWER_USCORE)
#pragma weak MPI_IALLGATHER = pmpi_iallgather_
#pragma weak mpi_iallgather__ = pmpi_iallgather_
#pragma weak mpi_iallgather_ = pmpi_iallgather_
#pragma weak mpi_iallgather = pmpi_iallgather_
#else
#pragma weak MPI_IALLGATHER = pmpi_iallgather
#pragma weak mpi_iallgather__ = pmpi_iallgather
#pragma weak mpi_iallgather_ = pmpi_iallgather
#pragma weak mpi_iallgather = pmpi_iallgather
#endif



#elif defined(HAVE_PRAGMA_WEAK)

#if defined(F77_NAME_UPPER)
extern FORT_DLL_SPEC void FORT_CALL MPI_IALLGATHER( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint * );

#pragma weak MPI_IALLGATHER = PMPI_IALLGATHER
#elif defined(F77_NAME_LOWER_2USCORE)
extern FORT_DLL_SPEC void FORT_CALL mpi_iallgather__( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint * );

#pragma weak mpi_iallgather__ = pmpi_iallgather__
#elif !defined(F77_NAME_LOWER_USCORE)
extern FORT_DLL_SPEC void FORT_CALL mpi_iallgather( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint * );

#pragma weak mpi_iallgather = pmpi_iallgather
#else
extern FORT_DLL_SPEC void FORT_CALL mpi_iallgather_( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint * );

#pragma weak mpi_iallgather_ = pmpi_iallgather_
#endif

#elif defined(HAVE_PRAGMA_HP_SEC_DEF)
#if defined(F77_NAME_UPPER)
#pragma _HP_SECONDARY_DEF PMPI_IALLGATHER  MPI_IALLGATHER
#elif defined(F77_NAME_LOWER_2USCORE)
#pragma _HP_SECONDARY_DEF pmpi_iallgather__  mpi_iallgather__
#elif !defined(F77_NAME_LOWER_USCORE)
#pragma _HP_SECONDARY_DEF pmpi_iallgather  mpi_iallgather
#else
#pragma _HP_SECONDARY_DEF pmpi_iallgather_  mpi_iallgather_
#endif

#elif defined(HAVE_PRAGMA_CRI_DUP)
#if defined(F77_NAME_UPPER)
#pragma _CRI duplicate MPI_IALLGATHER as PMPI_IALLGATHER
#elif defined(F77_NAME_LOWER_2USCORE)
#pragma _CRI duplicate mpi_iallgather__ as pmpi_iallgather__
#elif !defined(F77_NAME_LOWER_USCORE)
#pragma _CRI duplicate mpi_iallgather as pmpi_iallgather
#else
#pragma _CRI duplicate mpi_iallgather_ as pmpi_iallgather_
#endif

#elif defined(HAVE_WEAK_ATTRIBUTE)
#if defined(F77_NAME_UPPER)
extern FORT_DLL_SPEC void FORT_CALL MPI_IALLGATHER( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("PMPI_IALLGATHER")));
extern FORT_DLL_SPEC void FORT_CALL mpi_iallgather__( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("PMPI_IALLGATHER")));
extern FORT_DLL_SPEC void FORT_CALL mpi_iallgather_( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("PMPI_IALLGATHER")));
extern FORT_DLL_SPEC void FORT_CALL mpi_iallgather( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("PMPI_IALLGATHER")));

#elif defined(F77_NAME_LOWER_2USCORE)
extern FORT_DLL_SPEC void FORT_CALL MPI_IALLGATHER( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_iallgather__")));
extern FORT_DLL_SPEC void FORT_CALL mpi_iallgather__( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_iallgather__")));
extern FORT_DLL_SPEC void FORT_CALL mpi_iallgather_( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_iallgather__")));
extern FORT_DLL_SPEC void FORT_CALL mpi_iallgather( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_iallgather__")));

#elif defined(F77_NAME_LOWER_USCORE)
extern FORT_DLL_SPEC void FORT_CALL MPI_IALLGATHER( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_iallgather_")));
extern FORT_DLL_SPEC void FORT_CALL mpi_iallgather__( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_iallgather_")));
extern FORT_DLL_SPEC void FORT_CALL mpi_iallgather_( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_iallgather_")));
extern FORT_DLL_SPEC void FORT_CALL mpi_iallgather( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_iallgather_")));

#else
extern FORT_DLL_SPEC void FORT_CALL MPI_IALLGATHER( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_iallgather")));
extern FORT_DLL_SPEC void FORT_CALL mpi_iallgather__( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_iallgather")));
extern FORT_DLL_SPEC void FORT_CALL mpi_iallgather_( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_iallgather")));
extern FORT_DLL_SPEC void FORT_CALL mpi_iallgather( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_iallgather")));

#endif
#endif /* HAVE_PRAGMA_WEAK */
#endif /* USE_WEAK_SYMBOLS */
/* End MPI profiling block */


/* These definitions are used only for generating the Fortran wrappers */
#if defined(USE_WEAK_SYMBOLS) && defined(USE_ONLY_MPI_NAMES)
#if defined(HAVE_MULTIPLE_PRAGMA_WEAK)
extern FORT_DLL_SPEC void FORT_CALL MPI_IALLGATHER( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint * );
extern FORT_DLL_SPEC void FORT_CALL mpi_iallgather__( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint * );
extern FORT_DLL_SPEC void FORT_CALL mpi_iallgather( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint * );
extern FORT_DLL_SPEC void FORT_CALL mpi_iallgather_( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint * );

#if defined(F77_NAME_UPPER)
#pragma weak mpi_iallgather__ = MPI_IALLGATHER
#pragma weak mpi_iallgather_ = MPI_IALLGATHER
#pragma weak mpi_iallgather = MPI_IALLGATHER
#elif defined(F77_NAME_LOWER_2USCORE)
#pragma weak MPI_IALLGATHER = mpi_iallgather__
#pragma weak mpi_iallgather_ = mpi_iallgather__
#pragma weak mpi_iallgather = mpi_iallgather__
#elif defined(F77_NAME_LOWER_USCORE)
#pragma weak MPI_IALLGATHER = mpi_iallgather_
#pragma weak mpi_iallgather__ = mpi_iallgather_
#pragma weak mpi_iallgather = mpi_iallgather_
#else
#pragma weak MPI_IALLGATHER = mpi_iallgather
#pragma weak mpi_iallgather__ = mpi_iallgather
#pragma weak mpi_iallgather_ = mpi_iallgather
#endif
#elif defined(HAVE_WEAK_ATTRIBUTE)
#if defined(F77_NAME_UPPER)
extern FORT_DLL_SPEC void FORT_CALL MPI_IALLGATHER( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint * );
extern FORT_DLL_SPEC void FORT_CALL mpi_iallgather__( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("MPI_IALLGATHER")));
extern FORT_DLL_SPEC void FORT_CALL mpi_iallgather_( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("MPI_IALLGATHER")));
extern FORT_DLL_SPEC void FORT_CALL mpi_iallgather( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("MPI_IALLGATHER")));

#elif defined(F77_NAME_LOWER_2USCORE)
extern FORT_DLL_SPEC void FORT_CALL MPI_IALLGATHER( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("mpi_iallgather__")));
extern FORT_DLL_SPEC void FORT_CALL mpi_iallgather__( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint * );
extern FORT_DLL_SPEC void FORT_CALL mpi_iallgather_( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("mpi_iallgather__")));
extern FORT_DLL_SPEC void FORT_CALL mpi_iallgather( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("mpi_iallgather__")));

#elif defined(F77_NAME_LOWER_USCORE)
extern FORT_DLL_SPEC void FORT_CALL MPI_IALLGATHER( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("mpi_iallgather_")));
extern FORT_DLL_SPEC void FORT_CALL mpi_iallgather__( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("mpi_iallgather_")));
extern FORT_DLL_SPEC void FORT_CALL mpi_iallgather_( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint * );
extern FORT_DLL_SPEC void FORT_CALL mpi_iallgather( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("mpi_iallgather_")));

#else
extern FORT_DLL_SPEC void FORT_CALL MPI_IALLGATHER( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("mpi_iallgather")));
extern FORT_DLL_SPEC void FORT_CALL mpi_iallgather__( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("mpi_iallgather")));
extern FORT_DLL_SPEC void FORT_CALL mpi_iallgather_( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("mpi_iallgather")));
extern FORT_DLL_SPEC void FORT_CALL mpi_iallgather( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint * );

#endif
#endif

#endif

/* Map the name to the correct form */
#ifndef MPICH_MPI_FROM_PMPI
#if defined(USE_WEAK_SYMBOLS)
#if defined(HAVE_MULTIPLE_PRAGMA_WEAK)
/* Define the weak versions of the PMPI routine*/
#ifndef F77_NAME_UPPER
extern FORT_DLL_SPEC void FORT_CALL PMPI_IALLGATHER( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint * );
#endif
#ifndef F77_NAME_LOWER_2USCORE
extern FORT_DLL_SPEC void FORT_CALL pmpi_iallgather__( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint * );
#endif
#ifndef F77_NAME_LOWER_USCORE
extern FORT_DLL_SPEC void FORT_CALL pmpi_iallgather_( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint * );
#endif
#ifndef F77_NAME_LOWER
extern FORT_DLL_SPEC void FORT_CALL pmpi_iallgather( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint * );

#endif

#if defined(F77_NAME_UPPER)
#pragma weak pmpi_iallgather__ = PMPI_IALLGATHER
#pragma weak pmpi_iallgather_ = PMPI_IALLGATHER
#pragma weak pmpi_iallgather = PMPI_IALLGATHER
#elif defined(F77_NAME_LOWER_2USCORE)
#pragma weak PMPI_IALLGATHER = pmpi_iallgather__
#pragma weak pmpi_iallgather_ = pmpi_iallgather__
#pragma weak pmpi_iallgather = pmpi_iallgather__
#elif defined(F77_NAME_LOWER_USCORE)
#pragma weak PMPI_IALLGATHER = pmpi_iallgather_
#pragma weak pmpi_iallgather__ = pmpi_iallgather_
#pragma weak pmpi_iallgather = pmpi_iallgather_
#else
#pragma weak PMPI_IALLGATHER = pmpi_iallgather
#pragma weak pmpi_iallgather__ = pmpi_iallgather
#pragma weak pmpi_iallgather_ = pmpi_iallgather
#endif /* Test on name mapping */

#elif defined(HAVE_WEAK_ATTRIBUTE)
#if defined(F77_NAME_UPPER)
extern FORT_DLL_SPEC void FORT_CALL pmpi_iallgather__( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("PMPI_IALLGATHER")));
extern FORT_DLL_SPEC void FORT_CALL pmpi_iallgather_( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("PMPI_IALLGATHER")));
extern FORT_DLL_SPEC void FORT_CALL pmpi_iallgather( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("PMPI_IALLGATHER")));

#elif defined(F77_NAME_LOWER_2USCORE)
extern FORT_DLL_SPEC void FORT_CALL PMPI_IALLGATHER( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_iallgather__")));
extern FORT_DLL_SPEC void FORT_CALL pmpi_iallgather_( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_iallgather__")));
extern FORT_DLL_SPEC void FORT_CALL pmpi_iallgather( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_iallgather__")));

#elif defined(F77_NAME_LOWER_USCORE)
extern FORT_DLL_SPEC void FORT_CALL PMPI_IALLGATHER( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_iallgather_")));
extern FORT_DLL_SPEC void FORT_CALL pmpi_iallgather__( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_iallgather_")));
extern FORT_DLL_SPEC void FORT_CALL pmpi_iallgather( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_iallgather_")));

#else
extern FORT_DLL_SPEC void FORT_CALL PMPI_IALLGATHER( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_iallgather")));
extern FORT_DLL_SPEC void FORT_CALL pmpi_iallgather__( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_iallgather")));
extern FORT_DLL_SPEC void FORT_CALL pmpi_iallgather_( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_iallgather")));

#endif /* Test on name mapping */
#endif /* HAVE_MULTIPLE_PRAGMA_WEAK */
#endif /* USE_WEAK_SYMBOLS */

#ifdef F77_NAME_UPPER
#define mpi_iallgather_ PMPI_IALLGATHER
#elif defined(F77_NAME_LOWER_2USCORE)
#define mpi_iallgather_ pmpi_iallgather__
#elif !defined(F77_NAME_LOWER_USCORE)
#define mpi_iallgather_ pmpi_iallgather
#else
#define mpi_iallgather_ pmpi_iallgather_
#endif /* Test on name mapping */

#ifdef F77_USE_PMPI
/* This defines the routine that we call, which must be the PMPI version
   since we're renaming the Fortran entry as the pmpi version.  The MPI name
   must be undefined first to prevent any conflicts with previous renamings. */
#undef MPI_Iallgather
#define MPI_Iallgather PMPI_Iallgather 
#endif

#else

#ifdef F77_NAME_UPPER
#define mpi_iallgather_ MPI_IALLGATHER
#elif defined(F77_NAME_LOWER_2USCORE)
#define mpi_iallgather_ mpi_iallgather__
#elif !defined(F77_NAME_LOWER_USCORE)
#define mpi_iallgather_ mpi_iallgather
/* Else leave name alone */
#endif


#endif /* MPICH_MPI_FROM_PMPI */

/* Prototypes for the Fortran interfaces */
#include "fproto.h"
FORT_DLL_SPEC void FORT_CALL mpi_iallgather_ ( void*v1, MPI_Fint *v2, MPI_Fint *v3, void*v4, MPI_Fint *v5, MPI_Fint *v6, MPI_Fint *v7, MPI_Fint *v8, MPI_Fint *ierr ){

#ifndef HAVE_MPI_F_INIT_WORKS_WITH_C
    if (MPIR_F_NeedInit){ mpirinitf_(); MPIR_F_NeedInit = 0; }
#endif
    if (v1 == MPIR_F_MPI_IN_PLACE) v1 = MPI_IN_PLACE;
    *ierr = MPI_Iallgather( v1, (int)*v2, (MPI_Datatype)(*v3), v4, (int)*v5, (MPI_Datatype)(*v6), (MPI_Comm)(*v7), (MPI_Request *)(v8) );
}
