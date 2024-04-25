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
extern FORT_DLL_SPEC void FORT_CALL MPI_NEIGHBOR_ALLGATHERV( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint [], MPI_Fint [], MPI_Fint *, MPI_Fint *, MPI_Fint * );
extern FORT_DLL_SPEC void FORT_CALL mpi_neighbor_allgatherv__( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint [], MPI_Fint [], MPI_Fint *, MPI_Fint *, MPI_Fint * );
extern FORT_DLL_SPEC void FORT_CALL mpi_neighbor_allgatherv( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint [], MPI_Fint [], MPI_Fint *, MPI_Fint *, MPI_Fint * );
extern FORT_DLL_SPEC void FORT_CALL mpi_neighbor_allgatherv_( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint [], MPI_Fint [], MPI_Fint *, MPI_Fint *, MPI_Fint * );

#if defined(F77_NAME_UPPER)
#pragma weak MPI_NEIGHBOR_ALLGATHERV = PMPI_NEIGHBOR_ALLGATHERV
#pragma weak mpi_neighbor_allgatherv__ = PMPI_NEIGHBOR_ALLGATHERV
#pragma weak mpi_neighbor_allgatherv_ = PMPI_NEIGHBOR_ALLGATHERV
#pragma weak mpi_neighbor_allgatherv = PMPI_NEIGHBOR_ALLGATHERV
#elif defined(F77_NAME_LOWER_2USCORE)
#pragma weak MPI_NEIGHBOR_ALLGATHERV = pmpi_neighbor_allgatherv__
#pragma weak mpi_neighbor_allgatherv__ = pmpi_neighbor_allgatherv__
#pragma weak mpi_neighbor_allgatherv_ = pmpi_neighbor_allgatherv__
#pragma weak mpi_neighbor_allgatherv = pmpi_neighbor_allgatherv__
#elif defined(F77_NAME_LOWER_USCORE)
#pragma weak MPI_NEIGHBOR_ALLGATHERV = pmpi_neighbor_allgatherv_
#pragma weak mpi_neighbor_allgatherv__ = pmpi_neighbor_allgatherv_
#pragma weak mpi_neighbor_allgatherv_ = pmpi_neighbor_allgatherv_
#pragma weak mpi_neighbor_allgatherv = pmpi_neighbor_allgatherv_
#else
#pragma weak MPI_NEIGHBOR_ALLGATHERV = pmpi_neighbor_allgatherv
#pragma weak mpi_neighbor_allgatherv__ = pmpi_neighbor_allgatherv
#pragma weak mpi_neighbor_allgatherv_ = pmpi_neighbor_allgatherv
#pragma weak mpi_neighbor_allgatherv = pmpi_neighbor_allgatherv
#endif



#elif defined(HAVE_PRAGMA_WEAK)

#if defined(F77_NAME_UPPER)
extern FORT_DLL_SPEC void FORT_CALL MPI_NEIGHBOR_ALLGATHERV( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint [], MPI_Fint [], MPI_Fint *, MPI_Fint *, MPI_Fint * );

#pragma weak MPI_NEIGHBOR_ALLGATHERV = PMPI_NEIGHBOR_ALLGATHERV
#elif defined(F77_NAME_LOWER_2USCORE)
extern FORT_DLL_SPEC void FORT_CALL mpi_neighbor_allgatherv__( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint [], MPI_Fint [], MPI_Fint *, MPI_Fint *, MPI_Fint * );

#pragma weak mpi_neighbor_allgatherv__ = pmpi_neighbor_allgatherv__
#elif !defined(F77_NAME_LOWER_USCORE)
extern FORT_DLL_SPEC void FORT_CALL mpi_neighbor_allgatherv( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint [], MPI_Fint [], MPI_Fint *, MPI_Fint *, MPI_Fint * );

#pragma weak mpi_neighbor_allgatherv = pmpi_neighbor_allgatherv
#else
extern FORT_DLL_SPEC void FORT_CALL mpi_neighbor_allgatherv_( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint [], MPI_Fint [], MPI_Fint *, MPI_Fint *, MPI_Fint * );

#pragma weak mpi_neighbor_allgatherv_ = pmpi_neighbor_allgatherv_
#endif

#elif defined(HAVE_PRAGMA_HP_SEC_DEF)
#if defined(F77_NAME_UPPER)
#pragma _HP_SECONDARY_DEF PMPI_NEIGHBOR_ALLGATHERV  MPI_NEIGHBOR_ALLGATHERV
#elif defined(F77_NAME_LOWER_2USCORE)
#pragma _HP_SECONDARY_DEF pmpi_neighbor_allgatherv__  mpi_neighbor_allgatherv__
#elif !defined(F77_NAME_LOWER_USCORE)
#pragma _HP_SECONDARY_DEF pmpi_neighbor_allgatherv  mpi_neighbor_allgatherv
#else
#pragma _HP_SECONDARY_DEF pmpi_neighbor_allgatherv_  mpi_neighbor_allgatherv_
#endif

#elif defined(HAVE_PRAGMA_CRI_DUP)
#if defined(F77_NAME_UPPER)
#pragma _CRI duplicate MPI_NEIGHBOR_ALLGATHERV as PMPI_NEIGHBOR_ALLGATHERV
#elif defined(F77_NAME_LOWER_2USCORE)
#pragma _CRI duplicate mpi_neighbor_allgatherv__ as pmpi_neighbor_allgatherv__
#elif !defined(F77_NAME_LOWER_USCORE)
#pragma _CRI duplicate mpi_neighbor_allgatherv as pmpi_neighbor_allgatherv
#else
#pragma _CRI duplicate mpi_neighbor_allgatherv_ as pmpi_neighbor_allgatherv_
#endif

#elif defined(HAVE_WEAK_ATTRIBUTE)
#if defined(F77_NAME_UPPER)
extern FORT_DLL_SPEC void FORT_CALL MPI_NEIGHBOR_ALLGATHERV( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint [], MPI_Fint [], MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("PMPI_NEIGHBOR_ALLGATHERV")));
extern FORT_DLL_SPEC void FORT_CALL mpi_neighbor_allgatherv__( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint [], MPI_Fint [], MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("PMPI_NEIGHBOR_ALLGATHERV")));
extern FORT_DLL_SPEC void FORT_CALL mpi_neighbor_allgatherv_( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint [], MPI_Fint [], MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("PMPI_NEIGHBOR_ALLGATHERV")));
extern FORT_DLL_SPEC void FORT_CALL mpi_neighbor_allgatherv( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint [], MPI_Fint [], MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("PMPI_NEIGHBOR_ALLGATHERV")));

#elif defined(F77_NAME_LOWER_2USCORE)
extern FORT_DLL_SPEC void FORT_CALL MPI_NEIGHBOR_ALLGATHERV( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint [], MPI_Fint [], MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_neighbor_allgatherv__")));
extern FORT_DLL_SPEC void FORT_CALL mpi_neighbor_allgatherv__( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint [], MPI_Fint [], MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_neighbor_allgatherv__")));
extern FORT_DLL_SPEC void FORT_CALL mpi_neighbor_allgatherv_( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint [], MPI_Fint [], MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_neighbor_allgatherv__")));
extern FORT_DLL_SPEC void FORT_CALL mpi_neighbor_allgatherv( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint [], MPI_Fint [], MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_neighbor_allgatherv__")));

#elif defined(F77_NAME_LOWER_USCORE)
extern FORT_DLL_SPEC void FORT_CALL MPI_NEIGHBOR_ALLGATHERV( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint [], MPI_Fint [], MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_neighbor_allgatherv_")));
extern FORT_DLL_SPEC void FORT_CALL mpi_neighbor_allgatherv__( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint [], MPI_Fint [], MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_neighbor_allgatherv_")));
extern FORT_DLL_SPEC void FORT_CALL mpi_neighbor_allgatherv_( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint [], MPI_Fint [], MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_neighbor_allgatherv_")));
extern FORT_DLL_SPEC void FORT_CALL mpi_neighbor_allgatherv( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint [], MPI_Fint [], MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_neighbor_allgatherv_")));

#else
extern FORT_DLL_SPEC void FORT_CALL MPI_NEIGHBOR_ALLGATHERV( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint [], MPI_Fint [], MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_neighbor_allgatherv")));
extern FORT_DLL_SPEC void FORT_CALL mpi_neighbor_allgatherv__( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint [], MPI_Fint [], MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_neighbor_allgatherv")));
extern FORT_DLL_SPEC void FORT_CALL mpi_neighbor_allgatherv_( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint [], MPI_Fint [], MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_neighbor_allgatherv")));
extern FORT_DLL_SPEC void FORT_CALL mpi_neighbor_allgatherv( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint [], MPI_Fint [], MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_neighbor_allgatherv")));

#endif
#endif /* HAVE_PRAGMA_WEAK */
#endif /* USE_WEAK_SYMBOLS */
/* End MPI profiling block */


/* These definitions are used only for generating the Fortran wrappers */
#if defined(USE_WEAK_SYMBOLS) && defined(USE_ONLY_MPI_NAMES)
#if defined(HAVE_MULTIPLE_PRAGMA_WEAK)
extern FORT_DLL_SPEC void FORT_CALL MPI_NEIGHBOR_ALLGATHERV( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint [], MPI_Fint [], MPI_Fint *, MPI_Fint *, MPI_Fint * );
extern FORT_DLL_SPEC void FORT_CALL mpi_neighbor_allgatherv__( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint [], MPI_Fint [], MPI_Fint *, MPI_Fint *, MPI_Fint * );
extern FORT_DLL_SPEC void FORT_CALL mpi_neighbor_allgatherv( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint [], MPI_Fint [], MPI_Fint *, MPI_Fint *, MPI_Fint * );
extern FORT_DLL_SPEC void FORT_CALL mpi_neighbor_allgatherv_( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint [], MPI_Fint [], MPI_Fint *, MPI_Fint *, MPI_Fint * );

#if defined(F77_NAME_UPPER)
#pragma weak mpi_neighbor_allgatherv__ = MPI_NEIGHBOR_ALLGATHERV
#pragma weak mpi_neighbor_allgatherv_ = MPI_NEIGHBOR_ALLGATHERV
#pragma weak mpi_neighbor_allgatherv = MPI_NEIGHBOR_ALLGATHERV
#elif defined(F77_NAME_LOWER_2USCORE)
#pragma weak MPI_NEIGHBOR_ALLGATHERV = mpi_neighbor_allgatherv__
#pragma weak mpi_neighbor_allgatherv_ = mpi_neighbor_allgatherv__
#pragma weak mpi_neighbor_allgatherv = mpi_neighbor_allgatherv__
#elif defined(F77_NAME_LOWER_USCORE)
#pragma weak MPI_NEIGHBOR_ALLGATHERV = mpi_neighbor_allgatherv_
#pragma weak mpi_neighbor_allgatherv__ = mpi_neighbor_allgatherv_
#pragma weak mpi_neighbor_allgatherv = mpi_neighbor_allgatherv_
#else
#pragma weak MPI_NEIGHBOR_ALLGATHERV = mpi_neighbor_allgatherv
#pragma weak mpi_neighbor_allgatherv__ = mpi_neighbor_allgatherv
#pragma weak mpi_neighbor_allgatherv_ = mpi_neighbor_allgatherv
#endif
#elif defined(HAVE_WEAK_ATTRIBUTE)
#if defined(F77_NAME_UPPER)
extern FORT_DLL_SPEC void FORT_CALL MPI_NEIGHBOR_ALLGATHERV( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint [], MPI_Fint [], MPI_Fint *, MPI_Fint *, MPI_Fint * );
extern FORT_DLL_SPEC void FORT_CALL mpi_neighbor_allgatherv__( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint [], MPI_Fint [], MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("MPI_NEIGHBOR_ALLGATHERV")));
extern FORT_DLL_SPEC void FORT_CALL mpi_neighbor_allgatherv_( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint [], MPI_Fint [], MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("MPI_NEIGHBOR_ALLGATHERV")));
extern FORT_DLL_SPEC void FORT_CALL mpi_neighbor_allgatherv( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint [], MPI_Fint [], MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("MPI_NEIGHBOR_ALLGATHERV")));

#elif defined(F77_NAME_LOWER_2USCORE)
extern FORT_DLL_SPEC void FORT_CALL MPI_NEIGHBOR_ALLGATHERV( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint [], MPI_Fint [], MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("mpi_neighbor_allgatherv__")));
extern FORT_DLL_SPEC void FORT_CALL mpi_neighbor_allgatherv__( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint [], MPI_Fint [], MPI_Fint *, MPI_Fint *, MPI_Fint * );
extern FORT_DLL_SPEC void FORT_CALL mpi_neighbor_allgatherv_( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint [], MPI_Fint [], MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("mpi_neighbor_allgatherv__")));
extern FORT_DLL_SPEC void FORT_CALL mpi_neighbor_allgatherv( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint [], MPI_Fint [], MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("mpi_neighbor_allgatherv__")));

#elif defined(F77_NAME_LOWER_USCORE)
extern FORT_DLL_SPEC void FORT_CALL MPI_NEIGHBOR_ALLGATHERV( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint [], MPI_Fint [], MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("mpi_neighbor_allgatherv_")));
extern FORT_DLL_SPEC void FORT_CALL mpi_neighbor_allgatherv__( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint [], MPI_Fint [], MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("mpi_neighbor_allgatherv_")));
extern FORT_DLL_SPEC void FORT_CALL mpi_neighbor_allgatherv_( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint [], MPI_Fint [], MPI_Fint *, MPI_Fint *, MPI_Fint * );
extern FORT_DLL_SPEC void FORT_CALL mpi_neighbor_allgatherv( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint [], MPI_Fint [], MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("mpi_neighbor_allgatherv_")));

#else
extern FORT_DLL_SPEC void FORT_CALL MPI_NEIGHBOR_ALLGATHERV( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint [], MPI_Fint [], MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("mpi_neighbor_allgatherv")));
extern FORT_DLL_SPEC void FORT_CALL mpi_neighbor_allgatherv__( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint [], MPI_Fint [], MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("mpi_neighbor_allgatherv")));
extern FORT_DLL_SPEC void FORT_CALL mpi_neighbor_allgatherv_( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint [], MPI_Fint [], MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("mpi_neighbor_allgatherv")));
extern FORT_DLL_SPEC void FORT_CALL mpi_neighbor_allgatherv( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint [], MPI_Fint [], MPI_Fint *, MPI_Fint *, MPI_Fint * );

#endif
#endif

#endif

/* Map the name to the correct form */
#ifndef MPICH_MPI_FROM_PMPI
#if defined(USE_WEAK_SYMBOLS)
#if defined(HAVE_MULTIPLE_PRAGMA_WEAK)
/* Define the weak versions of the PMPI routine*/
#ifndef F77_NAME_UPPER
extern FORT_DLL_SPEC void FORT_CALL PMPI_NEIGHBOR_ALLGATHERV( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint [], MPI_Fint [], MPI_Fint *, MPI_Fint *, MPI_Fint * );
#endif
#ifndef F77_NAME_LOWER_2USCORE
extern FORT_DLL_SPEC void FORT_CALL pmpi_neighbor_allgatherv__( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint [], MPI_Fint [], MPI_Fint *, MPI_Fint *, MPI_Fint * );
#endif
#ifndef F77_NAME_LOWER_USCORE
extern FORT_DLL_SPEC void FORT_CALL pmpi_neighbor_allgatherv_( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint [], MPI_Fint [], MPI_Fint *, MPI_Fint *, MPI_Fint * );
#endif
#ifndef F77_NAME_LOWER
extern FORT_DLL_SPEC void FORT_CALL pmpi_neighbor_allgatherv( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint [], MPI_Fint [], MPI_Fint *, MPI_Fint *, MPI_Fint * );

#endif

#if defined(F77_NAME_UPPER)
#pragma weak pmpi_neighbor_allgatherv__ = PMPI_NEIGHBOR_ALLGATHERV
#pragma weak pmpi_neighbor_allgatherv_ = PMPI_NEIGHBOR_ALLGATHERV
#pragma weak pmpi_neighbor_allgatherv = PMPI_NEIGHBOR_ALLGATHERV
#elif defined(F77_NAME_LOWER_2USCORE)
#pragma weak PMPI_NEIGHBOR_ALLGATHERV = pmpi_neighbor_allgatherv__
#pragma weak pmpi_neighbor_allgatherv_ = pmpi_neighbor_allgatherv__
#pragma weak pmpi_neighbor_allgatherv = pmpi_neighbor_allgatherv__
#elif defined(F77_NAME_LOWER_USCORE)
#pragma weak PMPI_NEIGHBOR_ALLGATHERV = pmpi_neighbor_allgatherv_
#pragma weak pmpi_neighbor_allgatherv__ = pmpi_neighbor_allgatherv_
#pragma weak pmpi_neighbor_allgatherv = pmpi_neighbor_allgatherv_
#else
#pragma weak PMPI_NEIGHBOR_ALLGATHERV = pmpi_neighbor_allgatherv
#pragma weak pmpi_neighbor_allgatherv__ = pmpi_neighbor_allgatherv
#pragma weak pmpi_neighbor_allgatherv_ = pmpi_neighbor_allgatherv
#endif /* Test on name mapping */

#elif defined(HAVE_WEAK_ATTRIBUTE)
#if defined(F77_NAME_UPPER)
extern FORT_DLL_SPEC void FORT_CALL pmpi_neighbor_allgatherv__( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint [], MPI_Fint [], MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("PMPI_NEIGHBOR_ALLGATHERV")));
extern FORT_DLL_SPEC void FORT_CALL pmpi_neighbor_allgatherv_( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint [], MPI_Fint [], MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("PMPI_NEIGHBOR_ALLGATHERV")));
extern FORT_DLL_SPEC void FORT_CALL pmpi_neighbor_allgatherv( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint [], MPI_Fint [], MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("PMPI_NEIGHBOR_ALLGATHERV")));

#elif defined(F77_NAME_LOWER_2USCORE)
extern FORT_DLL_SPEC void FORT_CALL PMPI_NEIGHBOR_ALLGATHERV( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint [], MPI_Fint [], MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_neighbor_allgatherv__")));
extern FORT_DLL_SPEC void FORT_CALL pmpi_neighbor_allgatherv_( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint [], MPI_Fint [], MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_neighbor_allgatherv__")));
extern FORT_DLL_SPEC void FORT_CALL pmpi_neighbor_allgatherv( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint [], MPI_Fint [], MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_neighbor_allgatherv__")));

#elif defined(F77_NAME_LOWER_USCORE)
extern FORT_DLL_SPEC void FORT_CALL PMPI_NEIGHBOR_ALLGATHERV( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint [], MPI_Fint [], MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_neighbor_allgatherv_")));
extern FORT_DLL_SPEC void FORT_CALL pmpi_neighbor_allgatherv__( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint [], MPI_Fint [], MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_neighbor_allgatherv_")));
extern FORT_DLL_SPEC void FORT_CALL pmpi_neighbor_allgatherv( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint [], MPI_Fint [], MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_neighbor_allgatherv_")));

#else
extern FORT_DLL_SPEC void FORT_CALL PMPI_NEIGHBOR_ALLGATHERV( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint [], MPI_Fint [], MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_neighbor_allgatherv")));
extern FORT_DLL_SPEC void FORT_CALL pmpi_neighbor_allgatherv__( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint [], MPI_Fint [], MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_neighbor_allgatherv")));
extern FORT_DLL_SPEC void FORT_CALL pmpi_neighbor_allgatherv_( void*, MPI_Fint *, MPI_Fint *, void*, MPI_Fint [], MPI_Fint [], MPI_Fint *, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_neighbor_allgatherv")));

#endif /* Test on name mapping */
#endif /* HAVE_MULTIPLE_PRAGMA_WEAK */
#endif /* USE_WEAK_SYMBOLS */

#ifdef F77_NAME_UPPER
#define mpi_neighbor_allgatherv_ PMPI_NEIGHBOR_ALLGATHERV
#elif defined(F77_NAME_LOWER_2USCORE)
#define mpi_neighbor_allgatherv_ pmpi_neighbor_allgatherv__
#elif !defined(F77_NAME_LOWER_USCORE)
#define mpi_neighbor_allgatherv_ pmpi_neighbor_allgatherv
#else
#define mpi_neighbor_allgatherv_ pmpi_neighbor_allgatherv_
#endif /* Test on name mapping */

#ifdef F77_USE_PMPI
/* This defines the routine that we call, which must be the PMPI version
   since we're renaming the Fortran entry as the pmpi version.  The MPI name
   must be undefined first to prevent any conflicts with previous renamings. */
#undef MPI_Neighbor_allgatherv
#define MPI_Neighbor_allgatherv PMPI_Neighbor_allgatherv 
#endif

#else

#ifdef F77_NAME_UPPER
#define mpi_neighbor_allgatherv_ MPI_NEIGHBOR_ALLGATHERV
#elif defined(F77_NAME_LOWER_2USCORE)
#define mpi_neighbor_allgatherv_ mpi_neighbor_allgatherv__
#elif !defined(F77_NAME_LOWER_USCORE)
#define mpi_neighbor_allgatherv_ mpi_neighbor_allgatherv
/* Else leave name alone */
#endif


#endif /* MPICH_MPI_FROM_PMPI */

/* Prototypes for the Fortran interfaces */
#include "fproto.h"
FORT_DLL_SPEC void FORT_CALL mpi_neighbor_allgatherv_ ( void*v1, MPI_Fint *v2, MPI_Fint *v3, void*v4, MPI_Fint v5[], MPI_Fint v6[], MPI_Fint *v7, MPI_Fint *v8, MPI_Fint *ierr ){
    *ierr = MPI_Neighbor_allgatherv( v1, (int)*v2, (MPI_Datatype)(*v3), v4, v5, v6, (MPI_Datatype)(*v7), (MPI_Comm)(*v8) );
}
