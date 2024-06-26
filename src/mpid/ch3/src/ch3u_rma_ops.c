/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */
/*
 *  (C) 2001 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */
/* Copyright (c) 2001-2016, The Ohio State University. All rights
 * reserved.
 *
 * This file is part of the MVAPICH2 software package developed by the
 * team members of The Ohio State University's Network-Based Computing
 * Laboratory (NBCL), headed by Professor Dhabaleswar K. (DK) Panda.
 *
 * For detailed copyright and licensing information, please refer to the
 * copyright file COPYRIGHT in the top level MVAPICH2 directory.
 *
 */

#include "mpidrma.h"
#include "coll_shmem.h"

static int enableShortACC=1;

MPIR_T_PVAR_DOUBLE_TIMER_DECL_EXTERN(RMA, rma_rmaqueue_alloc);
MPIR_T_PVAR_DOUBLE_TIMER_DECL_EXTERN(RMA, rma_rmaqueue_set);

#define MPIDI_PASSIVE_TARGET_DONE_TAG  348297
#define MPIDI_PASSIVE_TARGET_RMA_TAG 563924

/* 
 * TODO:
 * Explore use of alternate allocation mechanisms for the RMA queue elements
 * (Because profiling has shown that queue element allocation/deallocation
 * can take a significant amount of time in the RMA operations).
 *    1: Current approach (uses perm memory malloc/free)
 *    2: Preallocate and maintain list (use perm memory malloc, but
 *       free onto window; use first; free on window free)
 *    3: Preallocate and maintain list (use separate memory, but free to
 *       thread/process; free in Finalize handler.  Option to use for
 *       single-threaded to avoid thread overheads)
 * Possible interface
 *    int MPIDI_RMAListAlloc(MPIDI_RMA_Op_t **a,MPID_Win *win)
 *    int MPIDI_RMAListFree(MPIDI_RMA_Op_t *a, MPID_Win *win)
 *    return value is error code (e.g., allocation failure).
 */

#undef FUNCNAME
#define FUNCNAME MPIDI_Win_free
#undef FCNAME
#define FCNAME MPIDI_QUOTE(FUNCNAME)
int MPIDI_Win_free(MPID_Win **win_ptr)
{
    int mpi_errno=MPI_SUCCESS;
    int in_use;
    MPID_Comm *comm_ptr;
    MPIDI_STATE_DECL(MPID_STATE_MPIDI_WIN_FREE);
        
    MPIDI_RMA_FUNC_ENTER(MPID_STATE_MPIDI_WIN_FREE);

    MPIU_ERR_CHKANDJUMP((*win_ptr)->epoch_state != MPIDI_EPOCH_NONE &&
                        /* OSU_MVAPICH */
                        (*win_ptr)->epoch_state != MPIDI_EPOCH_FENCE,
                        /* OSU_MVAPICH */
                        mpi_errno, MPI_ERR_RMA_SYNC, "**rmasync");

    if (!(*win_ptr)->shm_allocated) {
        /* when SHM is allocated, we already waited for operation completion in
         MPIDI_CH3_SHM_Win_free, so we do not need to do it again here. */
        mpi_errno = MPIDI_CH3I_Wait_for_pt_ops_finish(*win_ptr);
        if (mpi_errno) MPIU_ERR_POP(mpi_errno);
    }

#if defined(CHANNEL_MRAIL)
    /*complete any pending outgoing operations*/
    if ((*win_ptr)->fall_back != 1) {
        MPIDI_CH3I_RDMA_finish_rma(*win_ptr);
    }

    /*you can be a passive target and hence should wait for any incoming 
      communication to complete*/
    if ((*win_ptr)->outstanding_rma != 0 || 
        (*win_ptr)->current_lock_type != MPID_LOCK_NONE) {
        MPID_Progress_state progress_state;
        
        MPID_Progress_start(&progress_state);
        while ((*win_ptr)->outstanding_rma != 0 || 
               (*win_ptr)->current_lock_type != MPID_LOCK_NONE) {
            mpi_errno = MPID_Progress_wait(&progress_state);
            /* --BEGIN ERROR HANDLING-- */
            if (mpi_errno != MPI_SUCCESS) {
                MPID_Progress_end(&progress_state);
                MPIU_ERR_SETANDJUMP(mpi_errno,MPI_ERR_OTHER,"**winnoprogress");
            }
            /* --END ERROR HANDLING-- */
        }
        MPID_Progress_end(&progress_state);
    }  

    if ((*win_ptr)->fall_back != 1) {
	MPIDI_CH3I_RDMA_win_free(win_ptr);
    }
    if( (!(*win_ptr)->shm_win_pt2pt) && (*win_ptr)->comm_ptr->dev.ch.shmem_coll_ok == 1) {
        free_2level_comm((*win_ptr)->comm_ptr);
    }
#endif /* defined(CHANNEL_MRAIL) */

#if defined (CHANNEL_PSM)
    MPIU_Free((*win_ptr)->rank_mapping);
#endif /* CHANNEL_PSM */    
    
    comm_ptr = (*win_ptr)->comm_ptr;
    mpi_errno = MPIR_Comm_free_impl(comm_ptr);
    if (mpi_errno) MPIU_ERR_POP(mpi_errno);

    MPIU_Free((*win_ptr)->targets);
    MPIU_Free((*win_ptr)->base_addrs);
    MPIU_Free((*win_ptr)->sizes);
    MPIU_Free((*win_ptr)->disp_units);
    MPIU_Free((*win_ptr)->all_win_handles);
    MPIU_Free((*win_ptr)->pt_rma_puts_accs);

    /* Free the attached buffer for windows created with MPI_Win_allocate() */
    if ((*win_ptr)->create_flavor == MPI_WIN_FLAVOR_ALLOCATE || (*win_ptr)->create_flavor == MPI_WIN_FLAVOR_SHARED) {
        if ((*win_ptr)->shm_allocated == FALSE && (*win_ptr)->size > 0) {
            MPIU_Free((*win_ptr)->base);
        }
    }

    MPIU_Object_release_ref(*win_ptr, &in_use);
    /* MPI windows don't have reference count semantics, so this should always be true */
    MPIU_Assert(!in_use);
    MPIU_Handle_obj_free( &MPID_Win_mem, *win_ptr );

 fn_exit:
    MPIDI_RMA_FUNC_EXIT(MPID_STATE_MPIDI_WIN_FREE);
    return mpi_errno;

 fn_fail:
    goto fn_exit;
}


#undef FUNCNAME
#define FUNCNAME MPIDI_Win_shared_query
#undef FCNAME
#define FCNAME MPIDI_QUOTE(FUNCNAME)
int MPIDI_Win_shared_query(MPID_Win *win_ptr, int target_rank, MPI_Aint *size,
                           int *disp_unit, void *baseptr)
{
    int mpi_errno = MPI_SUCCESS;

    MPIDI_STATE_DECL(MPID_STATE_MPIDI_WIN_SHARED_QUERY);
    MPIDI_RMA_FUNC_ENTER(MPID_STATE_MPIDI_WIN_SHARED_QUERY);

    *(void**) baseptr = win_ptr->base;
    *size             = win_ptr->size;
    *disp_unit        = win_ptr->disp_unit;

 fn_exit:
    MPIDI_RMA_FUNC_EXIT(MPID_STATE_MPIDI_WIN_SHARED_QUERY);
    return mpi_errno;
    /* --BEGIN ERROR HANDLING-- */
 fn_fail:
    goto fn_exit;
    /* --END ERROR HANDLING-- */
}


#undef FUNCNAME
#define FUNCNAME MPIDI_Put
#undef FCNAME
#define FCNAME MPIDI_QUOTE(FUNCNAME)
int MPIDI_Put(const void *origin_addr, int origin_count, MPI_Datatype
            origin_datatype, int target_rank, MPI_Aint target_disp,
            int target_count, MPI_Datatype target_datatype, MPID_Win *win_ptr)
{
    int mpi_errno = MPI_SUCCESS;
    int dt_contig ATTRIBUTE((unused)), rank;
    MPID_Datatype *dtp;
    MPI_Aint dt_true_lb ATTRIBUTE((unused));
    MPIDI_msg_sz_t data_sz;
    MPIDI_VC_t *orig_vc = NULL, *target_vc = NULL;
    MPIDI_STATE_DECL(MPID_STATE_MPIDI_PUT);

#if defined(CHANNEL_MRAIL)
    int transfer_complete = 0;
    MPIDI_msg_sz_t size, target_type_size;
#endif
        
    MPIDI_RMA_FUNC_ENTER(MPID_STATE_MPIDI_PUT);

    if (target_rank == MPI_PROC_NULL) {
        goto fn_exit;
    }

    if (win_ptr->epoch_state == MPIDI_EPOCH_NONE && win_ptr->fence_issued) {
        win_ptr->epoch_state = MPIDI_EPOCH_FENCE;
    }

    MPIU_ERR_CHKANDJUMP(win_ptr->epoch_state == MPIDI_EPOCH_NONE,
                        mpi_errno, MPI_ERR_RMA_SYNC, "**rmasync");

    MPIDI_Datatype_get_info(origin_count, origin_datatype,
			    dt_contig, data_sz, dtp,dt_true_lb); 
    
    if (data_sz == 0) {
	goto fn_exit;
    }

    rank = win_ptr->comm_ptr->rank;
    
    if (win_ptr->shm_allocated == TRUE && target_rank != rank && win_ptr->create_flavor != MPI_WIN_FLAVOR_SHARED) {
        /* check if target is local and shared memory is allocated on window,
           if so, we directly perform this operation on shared memory region. */

        /* FIXME: Here we decide whether to perform SHM operations by checking if origin and target are on
           the same node. However, in ch3:sock, even if origin and target are on the same node, they do
           not within the same SHM region. Here we filter out ch3:sock by checking shm_allocated flag first,
           which is only set to TRUE when SHM region is allocated in nemesis.
           In future we need to figure out a way to check if origin and target are in the same "SHM comm".
        */
        MPIDI_Comm_get_vc(win_ptr->comm_ptr, rank, &orig_vc);
        MPIDI_Comm_get_vc(win_ptr->comm_ptr, target_rank, &target_vc);
    }

    /* If the put is a local operation, do it here */
    if (target_rank == rank || (
#if defined(CHANNEL_MRAIL) || defined(CHANNEL_PSM)
        (win_ptr->use_direct_shm == 1) && 
#endif
        (win_ptr->create_flavor == MPI_WIN_FLAVOR_SHARED ||
        (win_ptr->shm_allocated == TRUE && orig_vc->node_id == target_vc->node_id))))
    {
        mpi_errno = MPIDI_CH3I_Shm_put_op(origin_addr, origin_count, origin_datatype, target_rank,
                                          target_disp, target_count, target_datatype, win_ptr);
        if (mpi_errno) MPIU_ERR_POP(mpi_errno);
    }
    else
    {

#if defined(CHANNEL_MRAIL)
   MPID_Datatype_get_size_macro(target_datatype, target_type_size);
   size = target_count * target_type_size;
   if (MPIR_DATATYPE_IS_PREDEFINED(origin_datatype) 
        && MPIR_DATATYPE_IS_PREDEFINED(target_datatype)
        && win_ptr->fall_back != 1 && win_ptr->enable_fast_path == 1
        && win_ptr->use_rdma_path == 1
        && ((win_ptr->is_active && win_ptr->post_flag[target_rank] == 1)
        || (!win_ptr->is_active && win_ptr->using_lock == 0))
        && size < rdma_large_msg_rail_sharing_threshold)
    {
        transfer_complete = MPIDI_CH3I_RDMA_try_rma_op_fast(MPIDI_RMA_PUT, (void *)origin_addr,
                origin_count, origin_datatype, target_rank, target_disp,
                target_count, target_datatype, NULL, NULL, win_ptr);
    }
    if (transfer_complete) {
        goto fn_exit;
    }
    else 
#endif
    {

        MPIDI_RMA_Ops_list_t *ops_list = MPIDI_CH3I_RMA_Get_ops_list(win_ptr, target_rank);
        MPIDI_RMA_Op_t *new_ptr = NULL;

	/* queue it up */
        MPIR_T_PVAR_TIMER_START(RMA, rma_rmaqueue_alloc);
        mpi_errno = MPIDI_CH3I_RMA_Ops_alloc_tail(ops_list, &new_ptr);
        MPIR_T_PVAR_TIMER_END(RMA, rma_rmaqueue_alloc);
        if (mpi_errno) { MPIU_ERR_POP(mpi_errno); }

	MPIR_T_PVAR_TIMER_START(RMA, rma_rmaqueue_set);
	/* FIXME: For contig and very short operations, use a streamlined op */
	new_ptr->type = MPIDI_RMA_PUT;
        /* Cast away const'ness for the origin address, as the
         * MPIDI_RMA_Op_t structure is used for both PUT and GET like
         * operations */
	new_ptr->origin_addr = (void *) origin_addr;
	new_ptr->origin_count = origin_count;
	new_ptr->origin_datatype = origin_datatype;
	new_ptr->target_rank = target_rank;
	new_ptr->target_disp = target_disp;
	new_ptr->target_count = target_count;
	new_ptr->target_datatype = target_datatype;
	MPIR_T_PVAR_TIMER_END(RMA, rma_rmaqueue_set);

	/* if source or target datatypes are derived, increment their
	   reference counts */
	if (!MPIR_DATATYPE_IS_PREDEFINED(origin_datatype))
	{
	    MPID_Datatype_get_ptr(origin_datatype, dtp);
	    MPID_Datatype_add_ref(dtp);
	}
	if (!MPIR_DATATYPE_IS_PREDEFINED(target_datatype))
	{
	    MPID_Datatype_get_ptr(target_datatype, dtp);
	    MPID_Datatype_add_ref(dtp);
	}
    }
    }

#if defined(CHANNEL_MRAIL) && !defined(_SCHEDULE)
    if (win_ptr->fall_back != 1 && win_ptr->using_lock != 1) {
        MPIDI_CH3I_RDMA_try_rma(win_ptr, target_rank);
    }
#endif /* defined(CHANNEL_MRAIL) && !defined(_SCHEDULE) */

  fn_exit:
    MPIDI_RMA_FUNC_EXIT(MPID_STATE_MPIDI_PUT);    
    return mpi_errno;

    /* --BEGIN ERROR HANDLING-- */
  fn_fail:
    goto fn_exit;
    /* --END ERROR HANDLING-- */
}



#undef FUNCNAME
#define FUNCNAME MPIDI_Get
#undef FCNAME
#define FCNAME MPIDI_QUOTE(FUNCNAME)
int MPIDI_Get(void *origin_addr, int origin_count, MPI_Datatype
            origin_datatype, int target_rank, MPI_Aint target_disp,
            int target_count, MPI_Datatype target_datatype, MPID_Win *win_ptr)
{
    int mpi_errno = MPI_SUCCESS;
    MPIDI_msg_sz_t data_sz;
    int dt_contig ATTRIBUTE((unused)), rank;
    MPI_Aint dt_true_lb ATTRIBUTE((unused));
    MPID_Datatype *dtp;
    MPIDI_VC_t *orig_vc = NULL, *target_vc = NULL;
    MPIDI_STATE_DECL(MPID_STATE_MPIDI_GET);

#if defined(CHANNEL_MRAIL)
    int transfer_complete = 0;
    MPIDI_msg_sz_t size, target_type_size;
#endif
        
    MPIDI_RMA_FUNC_ENTER(MPID_STATE_MPIDI_GET);

    if (target_rank == MPI_PROC_NULL) {
        goto fn_exit;
    }

    if (win_ptr->epoch_state == MPIDI_EPOCH_NONE && win_ptr->fence_issued) {
        win_ptr->epoch_state = MPIDI_EPOCH_FENCE;
    }

    MPIU_ERR_CHKANDJUMP(win_ptr->epoch_state == MPIDI_EPOCH_NONE,
                        mpi_errno, MPI_ERR_RMA_SYNC, "**rmasync");

    MPIDI_Datatype_get_info(origin_count, origin_datatype,
			    dt_contig, data_sz, dtp, dt_true_lb); 

    if (data_sz == 0) {
	goto fn_exit;
    }

    rank = win_ptr->comm_ptr->rank;

    if (win_ptr->shm_allocated == TRUE && target_rank != rank && win_ptr->create_flavor != MPI_WIN_FLAVOR_SHARED) {
        /* check if target is local and shared memory is allocated on window,
           if so, we directly perform this operation on shared memory region. */

        /* FIXME: Here we decide whether to perform SHM operations by checking if origin and target are on
           the same node. However, in ch3:sock, even if origin and target are on the same node, they do
           not within the same SHM region. Here we filter out ch3:sock by checking shm_allocated flag first,
           which is only set to TRUE when SHM region is allocated in nemesis.
           In future we need to figure out a way to check if origin and target are in the same "SHM comm".
        */
        MPIDI_Comm_get_vc(win_ptr->comm_ptr, rank, &orig_vc);
        MPIDI_Comm_get_vc(win_ptr->comm_ptr, target_rank, &target_vc);
    }
    
    /* If the get is a local operation, do it here */
    if (target_rank == rank || (
#if defined(CHANNEL_MRAIL) || defined(CHANNEL_PSM)
        (win_ptr->use_direct_shm == 1) &&
#endif
        (win_ptr->create_flavor == MPI_WIN_FLAVOR_SHARED ||
        (win_ptr->shm_allocated == TRUE && orig_vc->node_id == target_vc->node_id))))
    {
        mpi_errno = MPIDI_CH3I_Shm_get_op(origin_addr, origin_count, origin_datatype, target_rank,
                                          target_disp, target_count, target_datatype, win_ptr);
        if (mpi_errno) MPIU_ERR_POP(mpi_errno);
    }
    else
    {

#if defined(CHANNEL_MRAIL)
   MPID_Datatype_get_size_macro(target_datatype, target_type_size);
   size = target_count * target_type_size;
	if (MPIR_DATATYPE_IS_PREDEFINED(origin_datatype) 
        && MPIR_DATATYPE_IS_PREDEFINED(target_datatype)
        && win_ptr->fall_back != 1 && win_ptr->enable_fast_path == 1
        && win_ptr->use_rdma_path == 1
        && ((win_ptr->is_active && win_ptr->post_flag[target_rank] == 1)
        || (!win_ptr->is_active && win_ptr->using_lock == 0))
        && size < rdma_large_msg_rail_sharing_threshold)
    {
        transfer_complete = MPIDI_CH3I_RDMA_try_rma_op_fast(MPIDI_RMA_GET, origin_addr,
                origin_count, origin_datatype, target_rank, target_disp,
                target_count, target_datatype, NULL, NULL, win_ptr);
    }
    if (transfer_complete) {
        goto fn_exit;
    }
    else 
#endif
    {
        MPIDI_RMA_Ops_list_t *ops_list = MPIDI_CH3I_RMA_Get_ops_list(win_ptr, target_rank);
        MPIDI_RMA_Op_t *new_ptr = NULL;

	/* queue it up */
        MPIR_T_PVAR_TIMER_START(RMA, rma_rmaqueue_alloc);
        mpi_errno = MPIDI_CH3I_RMA_Ops_alloc_tail(ops_list, &new_ptr);
        MPIR_T_PVAR_TIMER_END(RMA, rma_rmaqueue_alloc);
        if (mpi_errno) { MPIU_ERR_POP(mpi_errno); }

	MPIR_T_PVAR_TIMER_START(RMA, rma_rmaqueue_set);
	/* FIXME: For contig and very short operations, use a streamlined op */
	new_ptr->type = MPIDI_RMA_GET;
	new_ptr->origin_addr = origin_addr;
	new_ptr->origin_count = origin_count;
	new_ptr->origin_datatype = origin_datatype;
	new_ptr->target_rank = target_rank;
	new_ptr->target_disp = target_disp;
	new_ptr->target_count = target_count;
	new_ptr->target_datatype = target_datatype;
	MPIR_T_PVAR_TIMER_END(RMA, rma_rmaqueue_set);
	
	/* if source or target datatypes are derived, increment their
	   reference counts */
	if (!MPIR_DATATYPE_IS_PREDEFINED(origin_datatype))
	{
	    MPID_Datatype_get_ptr(origin_datatype, dtp);
	    MPID_Datatype_add_ref(dtp);
	}
	if (!MPIR_DATATYPE_IS_PREDEFINED(target_datatype))
	{
	    MPID_Datatype_get_ptr(target_datatype, dtp);
	    MPID_Datatype_add_ref(dtp);
	}
    }
    }

#if defined(CHANNEL_MRAIL) && !defined(_SCHEDULE)
    if (win_ptr->fall_back != 1 && win_ptr->using_lock != 1) {
        MPIDI_CH3I_RDMA_try_rma(win_ptr, target_rank);
    }
#endif /* defined(CHANNEL_MRAIL) && !defined(_SCHEDULE) */

  fn_exit:
    MPIDI_RMA_FUNC_EXIT(MPID_STATE_MPIDI_GET);
    return mpi_errno;

    /* --BEGIN ERROR HANDLING-- */
  fn_fail:
    goto fn_exit;
    /* --END ERROR HANDLING-- */
}



#undef FUNCNAME
#define FUNCNAME MPIDI_Accumulate
#undef FCNAME
#define FCNAME MPIDI_QUOTE(FUNCNAME)
int MPIDI_Accumulate(const void *origin_addr, int origin_count, MPI_Datatype
                    origin_datatype, int target_rank, MPI_Aint target_disp,
                    int target_count, MPI_Datatype target_datatype, MPI_Op op,
                    MPID_Win *win_ptr)
{
    int mpi_errno=MPI_SUCCESS;
    MPIDI_msg_sz_t data_sz;
    int dt_contig ATTRIBUTE((unused)), rank;
    MPI_Aint dt_true_lb ATTRIBUTE((unused));
    MPID_Datatype *dtp;
    MPIDI_VC_t *orig_vc = NULL, *target_vc = NULL;
    MPIDI_STATE_DECL(MPID_STATE_MPIDI_ACCUMULATE);
    
    MPIDI_RMA_FUNC_ENTER(MPID_STATE_MPIDI_ACCUMULATE);

    if (target_rank == MPI_PROC_NULL) {
        goto fn_exit;
    }

    if (win_ptr->epoch_state == MPIDI_EPOCH_NONE && win_ptr->fence_issued) {
        win_ptr->epoch_state = MPIDI_EPOCH_FENCE;
    }

    MPIU_ERR_CHKANDJUMP(win_ptr->epoch_state == MPIDI_EPOCH_NONE,
                        mpi_errno, MPI_ERR_RMA_SYNC, "**rmasync");

    MPIDI_Datatype_get_info(origin_count, origin_datatype,
			    dt_contig, data_sz, dtp, dt_true_lb);  
    
    if (data_sz == 0) {
	goto fn_exit;
    }

    rank = win_ptr->comm_ptr->rank;
    
    if (win_ptr->shm_allocated == TRUE && target_rank != rank && win_ptr->create_flavor != MPI_WIN_FLAVOR_SHARED) {
        /* check if target is local and shared memory is allocated on window,
           if so, we directly perform this operation on shared memory region. */

        /* FIXME: Here we decide whether to perform SHM operations by checking if origin and target are on
           the same node. However, in ch3:sock, even if origin and target are on the same node, they do
           not within the same SHM region. Here we filter out ch3:sock by checking shm_allocated flag first,
           which is only set to TRUE when SHM region is allocated in nemesis.
           In future we need to figure out a way to check if origin and target are in the same "SHM comm".
        */
        MPIDI_Comm_get_vc(win_ptr->comm_ptr, rank, &orig_vc);
        MPIDI_Comm_get_vc(win_ptr->comm_ptr, target_rank, &target_vc);
    }

    /* Do =! rank first (most likely branch?) */
    if (target_rank == rank || (
#if defined(CHANNEL_MRAIL) || defined(CHANNEL_PSM)
        (win_ptr->use_direct_shm == 1) &&
#endif
        (win_ptr->create_flavor == MPI_WIN_FLAVOR_SHARED ||
        (win_ptr->shm_allocated == TRUE && orig_vc->node_id == target_vc->node_id))))
    {
	mpi_errno = MPIDI_CH3I_Shm_acc_op(origin_addr, origin_count, origin_datatype,
					  target_rank, target_disp, target_count, target_datatype,
					  op, win_ptr);
	if (mpi_errno) MPIU_ERR_POP(mpi_errno);
    }
    else
    {
        MPIDI_RMA_Ops_list_t *ops_list = MPIDI_CH3I_RMA_Get_ops_list(win_ptr, target_rank);
        MPIDI_RMA_Op_t *new_ptr = NULL;

	/* queue it up */
        MPIR_T_PVAR_TIMER_START(RMA, rma_rmaqueue_alloc);
        mpi_errno = MPIDI_CH3I_RMA_Ops_alloc_tail(ops_list, &new_ptr);
        MPIR_T_PVAR_TIMER_END(RMA, rma_rmaqueue_alloc);
        if (mpi_errno) { MPIU_ERR_POP(mpi_errno); }

	/* If predefined and contiguous, use a simplified element */
	if (MPIR_DATATYPE_IS_PREDEFINED(origin_datatype) &&
            MPIR_DATATYPE_IS_PREDEFINED(target_datatype) && enableShortACC) {
	    MPIR_T_PVAR_TIMER_START(RMA, rma_rmaqueue_set);
	    new_ptr->type = MPIDI_RMA_ACC_CONTIG;
	    /* Only the information needed for the contig/predefined acc */
            /* Cast away const'ness for origin_address as
             * MPIDI_RMA_Op_t contain both PUT and GET like ops */
	    new_ptr->origin_addr = (void *) origin_addr;
	    new_ptr->origin_count = origin_count;
	    new_ptr->origin_datatype = origin_datatype;
	    new_ptr->target_rank = target_rank;
	    new_ptr->target_disp = target_disp;
	    new_ptr->target_count = target_count;
	    new_ptr->target_datatype = target_datatype;
	    new_ptr->op = op;
	    MPIR_T_PVAR_TIMER_END(RMA, rma_rmaqueue_set);
	    goto fn_exit;
	}

	MPIR_T_PVAR_TIMER_START(RMA, rma_rmaqueue_set);
	new_ptr->type = MPIDI_RMA_ACCUMULATE;
        /* Cast away const'ness for origin_address as MPIDI_RMA_Op_t
         * contain both PUT and GET like ops */
	new_ptr->origin_addr = (void *) origin_addr;
	new_ptr->origin_count = origin_count;
	new_ptr->origin_datatype = origin_datatype;
	new_ptr->target_rank = target_rank;
	new_ptr->target_disp = target_disp;
	new_ptr->target_count = target_count;
	new_ptr->target_datatype = target_datatype;
	new_ptr->op = op;
	MPIR_T_PVAR_TIMER_END(RMA, rma_rmaqueue_set);
	
	/* if source or target datatypes are derived, increment their
	   reference counts */
	if (!MPIR_DATATYPE_IS_PREDEFINED(origin_datatype))
	{
	    MPID_Datatype_get_ptr(origin_datatype, dtp);
	    MPID_Datatype_add_ref(dtp);
	}
	if (!MPIR_DATATYPE_IS_PREDEFINED(target_datatype))
	{
	    MPID_Datatype_get_ptr(target_datatype, dtp);
	    MPID_Datatype_add_ref(dtp);
	}
    }

 fn_exit:
    MPIDI_RMA_FUNC_EXIT(MPID_STATE_MPIDI_ACCUMULATE);
    return mpi_errno;

    /* --BEGIN ERROR HANDLING-- */
  fn_fail:
    goto fn_exit;
    /* --END ERROR HANDLING-- */
}


#undef FUNCNAME
#define FUNCNAME MPIDI_Alloc_mem
#undef FCNAME
#define FCNAME MPIDI_QUOTE(FUNCNAME)
void *MPIDI_Alloc_mem( size_t size, MPID_Info *info_ptr )
{
    void *ap;
    MPIDI_STATE_DECL(MPID_STATE_MPIDI_ALLOC_MEM);

    MPIDI_FUNC_ENTER(MPID_STATE_MPIDI_ALLOC_MEM);

#if defined (CHANNEL_MRAIL)
    ap = MPIDI_CH3I_Alloc_mem(size, info_ptr);
#else
    ap = MPIU_Malloc(size);
#endif
    
    MPIDI_FUNC_EXIT(MPID_STATE_MPIDI_ALLOC_MEM);
    return ap;
}


#undef FUNCNAME
#define FUNCNAME MPIDI_Free_mem
#undef FCNAME
#define FCNAME MPIDI_QUOTE(FUNCNAME)
int MPIDI_Free_mem( void *ptr )
{
    int mpi_errno = MPI_SUCCESS;
    MPIDI_STATE_DECL(MPID_STATE_MPIDI_FREE_MEM);

    MPIDI_FUNC_ENTER(MPID_STATE_MPIDI_FREE_MEM);

#if defined(CHANNEL_MRAIL)
    MPIDI_CH3I_Free_mem(ptr);
#else
    MPIU_Free(ptr);
#endif
    
    MPIDI_FUNC_EXIT(MPID_STATE_MPIDI_FREE_MEM);
    return mpi_errno;
}
