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
extern FORT_DLL_SPEC void FORT_CALL MPI_ADDRESS( void*, MPI_Fint *, MPI_Fint * );
extern FORT_DLL_SPEC void FORT_CALL mpi_address__( void*, MPI_Fint *, MPI_Fint * );
extern FORT_DLL_SPEC void FORT_CALL mpi_address( void*, MPI_Fint *, MPI_Fint * );
extern FORT_DLL_SPEC void FORT_CALL mpi_address_( void*, MPI_Fint *, MPI_Fint * );

#if defined(F77_NAME_UPPER)
#pragma weak MPI_ADDRESS = PMPI_ADDRESS
#pragma weak mpi_address__ = PMPI_ADDRESS
#pragma weak mpi_address_ = PMPI_ADDRESS
#pragma weak mpi_address = PMPI_ADDRESS
#elif defined(F77_NAME_LOWER_2USCORE)
#pragma weak MPI_ADDRESS = pmpi_address__
#pragma weak mpi_address__ = pmpi_address__
#pragma weak mpi_address_ = pmpi_address__
#pragma weak mpi_address = pmpi_address__
#elif defined(F77_NAME_LOWER_USCORE)
#pragma weak MPI_ADDRESS = pmpi_address_
#pragma weak mpi_address__ = pmpi_address_
#pragma weak mpi_address_ = pmpi_address_
#pragma weak mpi_address = pmpi_address_
#else
#pragma weak MPI_ADDRESS = pmpi_address
#pragma weak mpi_address__ = pmpi_address
#pragma weak mpi_address_ = pmpi_address
#pragma weak mpi_address = pmpi_address
#endif



#elif defined(HAVE_PRAGMA_WEAK)

#if defined(F77_NAME_UPPER)
extern FORT_DLL_SPEC void FORT_CALL MPI_ADDRESS( void*, MPI_Fint *, MPI_Fint * );

#pragma weak MPI_ADDRESS = PMPI_ADDRESS
#elif defined(F77_NAME_LOWER_2USCORE)
extern FORT_DLL_SPEC void FORT_CALL mpi_address__( void*, MPI_Fint *, MPI_Fint * );

#pragma weak mpi_address__ = pmpi_address__
#elif !defined(F77_NAME_LOWER_USCORE)
extern FORT_DLL_SPEC void FORT_CALL mpi_address( void*, MPI_Fint *, MPI_Fint * );

#pragma weak mpi_address = pmpi_address
#else
extern FORT_DLL_SPEC void FORT_CALL mpi_address_( void*, MPI_Fint *, MPI_Fint * );

#pragma weak mpi_address_ = pmpi_address_
#endif

#elif defined(HAVE_PRAGMA_HP_SEC_DEF)
#if defined(F77_NAME_UPPER)
#pragma _HP_SECONDARY_DEF PMPI_ADDRESS  MPI_ADDRESS
#elif defined(F77_NAME_LOWER_2USCORE)
#pragma _HP_SECONDARY_DEF pmpi_address__  mpi_address__
#elif !defined(F77_NAME_LOWER_USCORE)
#pragma _HP_SECONDARY_DEF pmpi_address  mpi_address
#else
#pragma _HP_SECONDARY_DEF pmpi_address_  mpi_address_
#endif

#elif defined(HAVE_PRAGMA_CRI_DUP)
#if defined(F77_NAME_UPPER)
#pragma _CRI duplicate MPI_ADDRESS as PMPI_ADDRESS
#elif defined(F77_NAME_LOWER_2USCORE)
#pragma _CRI duplicate mpi_address__ as pmpi_address__
#elif !defined(F77_NAME_LOWER_USCORE)
#pragma _CRI duplicate mpi_address as pmpi_address
#else
#pragma _CRI duplicate mpi_address_ as pmpi_address_
#endif

#elif defined(HAVE_WEAK_ATTRIBUTE)
#if defined(F77_NAME_UPPER)
extern FORT_DLL_SPEC void FORT_CALL MPI_ADDRESS( void*, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("PMPI_ADDRESS")));
extern FORT_DLL_SPEC void FORT_CALL mpi_address__( void*, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("PMPI_ADDRESS")));
extern FORT_DLL_SPEC void FORT_CALL mpi_address_( void*, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("PMPI_ADDRESS")));
extern FORT_DLL_SPEC void FORT_CALL mpi_address( void*, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("PMPI_ADDRESS")));

#elif defined(F77_NAME_LOWER_2USCORE)
extern FORT_DLL_SPEC void FORT_CALL MPI_ADDRESS( void*, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_address__")));
extern FORT_DLL_SPEC void FORT_CALL mpi_address__( void*, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_address__")));
extern FORT_DLL_SPEC void FORT_CALL mpi_address_( void*, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_address__")));
extern FORT_DLL_SPEC void FORT_CALL mpi_address( void*, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_address__")));

#elif defined(F77_NAME_LOWER_USCORE)
extern FORT_DLL_SPEC void FORT_CALL MPI_ADDRESS( void*, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_address_")));
extern FORT_DLL_SPEC void FORT_CALL mpi_address__( void*, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_address_")));
extern FORT_DLL_SPEC void FORT_CALL mpi_address_( void*, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_address_")));
extern FORT_DLL_SPEC void FORT_CALL mpi_address( void*, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_address_")));

#else
extern FORT_DLL_SPEC void FORT_CALL MPI_ADDRESS( void*, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_address")));
extern FORT_DLL_SPEC void FORT_CALL mpi_address__( void*, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_address")));
extern FORT_DLL_SPEC void FORT_CALL mpi_address_( void*, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_address")));
extern FORT_DLL_SPEC void FORT_CALL mpi_address( void*, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_address")));

#endif
#endif /* HAVE_PRAGMA_WEAK */
#endif /* USE_WEAK_SYMBOLS */
/* End MPI profiling block */


/* These definitions are used only for generating the Fortran wrappers */
#if defined(USE_WEAK_SYMBOLS) && defined(USE_ONLY_MPI_NAMES)
#if defined(HAVE_MULTIPLE_PRAGMA_WEAK)
extern FORT_DLL_SPEC void FORT_CALL MPI_ADDRESS( void*, MPI_Fint *, MPI_Fint * );
extern FORT_DLL_SPEC void FORT_CALL mpi_address__( void*, MPI_Fint *, MPI_Fint * );
extern FORT_DLL_SPEC void FORT_CALL mpi_address( void*, MPI_Fint *, MPI_Fint * );
extern FORT_DLL_SPEC void FORT_CALL mpi_address_( void*, MPI_Fint *, MPI_Fint * );

#if defined(F77_NAME_UPPER)
#pragma weak mpi_address__ = MPI_ADDRESS
#pragma weak mpi_address_ = MPI_ADDRESS
#pragma weak mpi_address = MPI_ADDRESS
#elif defined(F77_NAME_LOWER_2USCORE)
#pragma weak MPI_ADDRESS = mpi_address__
#pragma weak mpi_address_ = mpi_address__
#pragma weak mpi_address = mpi_address__
#elif defined(F77_NAME_LOWER_USCORE)
#pragma weak MPI_ADDRESS = mpi_address_
#pragma weak mpi_address__ = mpi_address_
#pragma weak mpi_address = mpi_address_
#else
#pragma weak MPI_ADDRESS = mpi_address
#pragma weak mpi_address__ = mpi_address
#pragma weak mpi_address_ = mpi_address
#endif
#elif defined(HAVE_WEAK_ATTRIBUTE)
#if defined(F77_NAME_UPPER)
extern FORT_DLL_SPEC void FORT_CALL MPI_ADDRESS( void*, MPI_Fint *, MPI_Fint * );
extern FORT_DLL_SPEC void FORT_CALL mpi_address__( void*, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("MPI_ADDRESS")));
extern FORT_DLL_SPEC void FORT_CALL mpi_address_( void*, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("MPI_ADDRESS")));
extern FORT_DLL_SPEC void FORT_CALL mpi_address( void*, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("MPI_ADDRESS")));

#elif defined(F77_NAME_LOWER_2USCORE)
extern FORT_DLL_SPEC void FORT_CALL MPI_ADDRESS( void*, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("mpi_address__")));
extern FORT_DLL_SPEC void FORT_CALL mpi_address__( void*, MPI_Fint *, MPI_Fint * );
extern FORT_DLL_SPEC void FORT_CALL mpi_address_( void*, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("mpi_address__")));
extern FORT_DLL_SPEC void FORT_CALL mpi_address( void*, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("mpi_address__")));

#elif defined(F77_NAME_LOWER_USCORE)
extern FORT_DLL_SPEC void FORT_CALL MPI_ADDRESS( void*, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("mpi_address_")));
extern FORT_DLL_SPEC void FORT_CALL mpi_address__( void*, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("mpi_address_")));
extern FORT_DLL_SPEC void FORT_CALL mpi_address_( void*, MPI_Fint *, MPI_Fint * );
extern FORT_DLL_SPEC void FORT_CALL mpi_address( void*, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("mpi_address_")));

#else
extern FORT_DLL_SPEC void FORT_CALL MPI_ADDRESS( void*, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("mpi_address")));
extern FORT_DLL_SPEC void FORT_CALL mpi_address__( void*, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("mpi_address")));
extern FORT_DLL_SPEC void FORT_CALL mpi_address_( void*, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("mpi_address")));
extern FORT_DLL_SPEC void FORT_CALL mpi_address( void*, MPI_Fint *, MPI_Fint * );

#endif
#endif

#endif

/* Map the name to the correct form */
#ifndef MPICH_MPI_FROM_PMPI
#if defined(USE_WEAK_SYMBOLS)
#if defined(HAVE_MULTIPLE_PRAGMA_WEAK)
/* Define the weak versions of the PMPI routine*/
#ifndef F77_NAME_UPPER
extern FORT_DLL_SPEC void FORT_CALL PMPI_ADDRESS( void*, MPI_Fint *, MPI_Fint * );
#endif
#ifndef F77_NAME_LOWER_2USCORE
extern FORT_DLL_SPEC void FORT_CALL pmpi_address__( void*, MPI_Fint *, MPI_Fint * );
#endif
#ifndef F77_NAME_LOWER_USCORE
extern FORT_DLL_SPEC void FORT_CALL pmpi_address_( void*, MPI_Fint *, MPI_Fint * );
#endif
#ifndef F77_NAME_LOWER
extern FORT_DLL_SPEC void FORT_CALL pmpi_address( void*, MPI_Fint *, MPI_Fint * );

#endif

#if defined(F77_NAME_UPPER)
#pragma weak pmpi_address__ = PMPI_ADDRESS
#pragma weak pmpi_address_ = PMPI_ADDRESS
#pragma weak pmpi_address = PMPI_ADDRESS
#elif defined(F77_NAME_LOWER_2USCORE)
#pragma weak PMPI_ADDRESS = pmpi_address__
#pragma weak pmpi_address_ = pmpi_address__
#pragma weak pmpi_address = pmpi_address__
#elif defined(F77_NAME_LOWER_USCORE)
#pragma weak PMPI_ADDRESS = pmpi_address_
#pragma weak pmpi_address__ = pmpi_address_
#pragma weak pmpi_address = pmpi_address_
#else
#pragma weak PMPI_ADDRESS = pmpi_address
#pragma weak pmpi_address__ = pmpi_address
#pragma weak pmpi_address_ = pmpi_address
#endif /* Test on name mapping */

#elif defined(HAVE_WEAK_ATTRIBUTE)
#if defined(F77_NAME_UPPER)
extern FORT_DLL_SPEC void FORT_CALL pmpi_address__( void*, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("PMPI_ADDRESS")));
extern FORT_DLL_SPEC void FORT_CALL pmpi_address_( void*, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("PMPI_ADDRESS")));
extern FORT_DLL_SPEC void FORT_CALL pmpi_address( void*, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("PMPI_ADDRESS")));

#elif defined(F77_NAME_LOWER_2USCORE)
extern FORT_DLL_SPEC void FORT_CALL PMPI_ADDRESS( void*, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_address__")));
extern FORT_DLL_SPEC void FORT_CALL pmpi_address_( void*, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_address__")));
extern FORT_DLL_SPEC void FORT_CALL pmpi_address( void*, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_address__")));

#elif defined(F77_NAME_LOWER_USCORE)
extern FORT_DLL_SPEC void FORT_CALL PMPI_ADDRESS( void*, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_address_")));
extern FORT_DLL_SPEC void FORT_CALL pmpi_address__( void*, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_address_")));
extern FORT_DLL_SPEC void FORT_CALL pmpi_address( void*, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_address_")));

#else
extern FORT_DLL_SPEC void FORT_CALL PMPI_ADDRESS( void*, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_address")));
extern FORT_DLL_SPEC void FORT_CALL pmpi_address__( void*, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_address")));
extern FORT_DLL_SPEC void FORT_CALL pmpi_address_( void*, MPI_Fint *, MPI_Fint * ) __attribute__((weak,alias("pmpi_address")));

#endif /* Test on name mapping */
#endif /* HAVE_MULTIPLE_PRAGMA_WEAK */
#endif /* USE_WEAK_SYMBOLS */

#ifdef F77_NAME_UPPER
#define mpi_address_ PMPI_ADDRESS
#elif defined(F77_NAME_LOWER_2USCORE)
#define mpi_address_ pmpi_address__
#elif !defined(F77_NAME_LOWER_USCORE)
#define mpi_address_ pmpi_address
#else
#define mpi_address_ pmpi_address_
#endif /* Test on name mapping */

#ifdef F77_USE_PMPI
/* This defines the routine that we call, which must be the PMPI version
   since we're renaming the Fortran entry as the pmpi version.  The MPI name
   must be undefined first to prevent any conflicts with previous renamings. */
#undef MPI_Address
#define MPI_Address PMPI_Address 
#endif

#else

#ifdef F77_NAME_UPPER
#define mpi_address_ MPI_ADDRESS
#elif defined(F77_NAME_LOWER_2USCORE)
#define mpi_address_ mpi_address__
#elif !defined(F77_NAME_LOWER_USCORE)
#define mpi_address_ mpi_address
/* Else leave name alone */
#endif


#endif /* MPICH_MPI_FROM_PMPI */

/* Prototypes for the Fortran interfaces */
#include "fproto.h"
#include "mpierrs.h"
#include <stdio.h>
#include "mpierror.h"
FORT_DLL_SPEC void FORT_CALL mpi_address_ ( void*v1, MPI_Fint *v2, MPI_Fint *ierr ){
    MPI_Aint a, b;
    *ierr = MPI_Address( v1, &a );

#ifndef HAVE_MPI_F_INIT_WORKS_WITH_C
    if (MPIR_F_NeedInit){ mpirinitf_(); MPIR_F_NeedInit = 0; }
#endif

#ifdef USE_POINTER_FOR_BOTTOM
    b = a;
#else
    b = a - (MPIR_Pint) MPIR_F_MPI_BOTTOM;
#endif
    *v2 = (MPI_Fint)( b );
#ifdef HAVE_AINT_LARGER_THAN_FINT
    /* Check for truncation */
    if ((MPI_Aint)*v2 - b != 0) {
        *ierr = MPIR_Err_create_code( MPI_SUCCESS, MPIR_ERR_RECOVERABLE, 
			  "MPI_Address", __LINE__, MPI_ERR_ARG, "**inttoosmall", 0 );
	(void)MPIR_Err_return_comm( 0, "MPI_Address",  *ierr );
    }
#endif
}
