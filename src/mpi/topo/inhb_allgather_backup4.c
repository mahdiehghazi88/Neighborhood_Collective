/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */
/*
 *  (C) 2012 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include "mpiimpl.h"
#include "topo.h"

/* -- Begin Profiling Symbol Block for routine MPI_Ineighbor_allgather */
#if defined(HAVE_PRAGMA_WEAK)
#pragma weak MPI_Ineighbor_allgather = PMPI_Ineighbor_allgather
#elif defined(HAVE_PRAGMA_HP_SEC_DEF)
#pragma _HP_SECONDARY_DEF PMPI_Ineighbor_allgather  MPI_Ineighbor_allgather
#elif defined(HAVE_PRAGMA_CRI_DUP)
#pragma _CRI duplicate MPI_Ineighbor_allgather as PMPI_Ineighbor_allgather
#elif defined(HAVE_WEAK_ATTRIBUTE)
int MPI_Ineighbor_allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                            void *recvbuf, int recvcount, MPI_Datatype recvtype,
                            MPI_Comm comm, MPI_Request *request)
                            __attribute__((weak,alias("PMPI_Ineighbor_allgather")));
#endif
/* -- End Profiling Symbol Block */

/* Define MPICH_MPI_FROM_PMPI if weak symbols are not supported to build
   the MPI routines */
#ifndef MPICH_MPI_FROM_PMPI
#undef MPI_Ineighbor_allgather
#define MPI_Ineighbor_allgather PMPI_Ineighbor_allgather

int find_incom_tmp_buff_size(int **incom_sched_mat, int num_rows, int size_per_rank)
{
    int i, buff_size;
    buff_size = 0;
    for(i = 0; i < num_rows; i++)
    {
    	if(incom_sched_mat[i][K_vrbl+2]==2)
    		buff_size += K_vrbl * size_per_rank;
    	else if (incom_sched_mat[i][K_vrbl+2]==1)
    		buff_size +=size_per_rank;
       /* if(incom_sched_mat[i][1] != -1) //Initially ON incoming neighbors  //yani be andazeye kolli memory bara t haye mokhtalef memory allocate mokonim
        {
        	if(!incom_sched_mat[i][K_vrbl+2])
            buff_size += K_vrbl * size_per_rank;  //incom_sched_mat[i][2] vaghti be kar mire ke cummulative combining dashte bashim, right?
        	else

        }*/
    }
    return buff_size;
}

int find_in_arr(int *array, int size, int value)
{
    int i, self_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &self_rank);
    for(i = 0; i < size; i++)
    {
    	//if(self_rank==0)
    	//printf("array[%d]=%d\n", i, array[i]);
        if(array[i] == value)
        {
        	//if(self_rank==0)
        		//printf("rank %d is in position %d\n", value, i);
            return i;
        }
    }
   // if(self_rank==0)
    //	printf("rank %d not found\n", value);
    return -1;
}

int find_incom_tmp_buf_offset(int **incom_sched_mat, int nbr_index, int size_per_rank)
{
    int i, offset;
    offset = 0;
    for(i = 0; i < nbr_index; i++)
    {
    	if(incom_sched_mat[i][K_vrbl+2]==2)
    		offset += K_vrbl * size_per_rank;
    	else if (incom_sched_mat[i][K_vrbl+2]==1)
    		offset +=size_per_rank;
       /* if(incom_sched_mat[i][1] !=-1) //Initially ON incoming neighbors
        {
            offset += K_vrbl * size_per_rank;
        }*/
    }
    return offset;
}

#undef FUNCNAME
#define FUNCNAME MPIR_Ineighbor_allgather_SMGM
#undef FCNAME
#define FCNAME MPIU_QUOTE(FUNCNAME)
int MPIR_Ineighbor_allgather_SMGM(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount,
                                MPI_Datatype recvtype, MPID_Comm *comm_ptr, MPID_Sched_t s, int iter)
{

	int mpi_errno = MPI_SUCCESS;
	int tag = -1;
	MPID_Request *reqp = NULL;
	mpi_errno = MPID_Sched_next_tag(comm_ptr, &tag);

	//*request=MPI_REQUEST_NULL;

	//At this point, we already have extracted all information
	//needed for building the schedule. The information is
	//attached to the topo.dist_graph field of the communicator
	//as a SHM_nbr_coll_patt structure. Now, we have to build
	//the schedule out of the information provided by this struct.
	MPIR_Topology *topo_ptr = NULL;
	topo_ptr = MPIR_Topology_get(comm_ptr);
	Common_nbrhood_matrix *cmn_nbh_mat = topo_ptr->topo.dist_graph.shm_nbh_coll_patt->cmn_nbh_mat;
	int **incom_sched_mat = topo_ptr->topo.dist_graph.shm_nbh_coll_patt->incom_sched_mat;

	//printf("MPIR_Ineighbor_allgather_SMGM\n");
	int self_rank = comm_ptr->rank;
	//int *test=(int*) sendbuf;

	int indegree, outdegree;
	int *srcs, *dests;
	MPI_Aint recvtype_extent, recvtype_true_lb, recvtype_true_extent, recvbuf_size;
	MPI_Aint sendtype_extent, sendbuf_size;
	MPI_Aint exchange_recv_count, exchange_send_count, incom_recv_count, out_send_count;
	MPI_Aint exchange_tmp_buf_size, incom_tmp_buf_size;

	MPID_Datatype_get_extent_macro(recvtype, recvtype_extent);
	MPID_Datatype_get_extent_macro(sendtype, sendtype_extent);


	void *exchange_tmp_buf = NULL;
	void *incom_tmp_buf = NULL;

	//exchange_tmp_buf = MPIU_Malloc(K_vrbl * sizeof(void*));

	indegree = topo_ptr->topo.dist_graph.indegree;
	outdegree = topo_ptr->topo.dist_graph.outdegree;
	srcs = topo_ptr->topo.dist_graph.in;
	dests = topo_ptr->topo.dist_graph.out;

	MPID_Datatype_get_extent_macro(recvtype, recvtype_extent);
	MPID_Datatype_get_extent_macro(sendtype, sendtype_extent);
	recvbuf_size = recvcount * recvtype_extent;
	sendbuf_size = sendcount * sendtype_extent;

	MPIR_SCHED_CHKPMEM_DECL(2);
	MPIU_CHKPMEM_DECL(2);
	exchange_tmp_buf_size = K_vrbl*recvbuf_size;
	//int sender_rank=10;

	//int cnt_sent_ind=0, cnt_sent_cmb=0;
	//int cnt_recv_ind=0, cnt_recv_comb=0;


//	if(self_rank==0)
//	printf("exchange_tmp_buf_size=%d, sendcount=%d, sendtype_extent=%d, sizeofint=%d\n", exchange_tmp_buf_size, sendcount, sendtype_extent, sizeof(int));

	incom_tmp_buf_size = find_incom_tmp_buff_size(incom_sched_mat, indegree, recvbuf_size);

	MPIR_SCHED_CHKPMEM_MALLOC(exchange_tmp_buf, void*, exchange_tmp_buf_size, mpi_errno, "exchange_tmp_buf");


	//MPIR_SCHED_CHKPMEM_MALLOC(test, int*, exchange_tmp_buf_size, mpi_errno, "test");

	int i,j;
	//for(i = 0; i < K_vrbl; i++)
		//exchange_tmp_buf[i] = MPIU_Malloc(exchange_tmp_buf_size * sizeof(void));

	MPIR_SCHED_CHKPMEM_MALLOC(incom_tmp_buf, void*, incom_tmp_buf_size, mpi_errno, "incom_tmp_buf");

	//MPID_SCHED_BARRIER(s);
	//printf("rank%d, pid%d: 1.\n", self_rank, getpid());
	  int t=0;
	 // printf("rank%d: cmn_nbh_mat->t=%d\n", self_rank, cmn_nbh_mat->t);
	    for(t = 0; t < cmn_nbh_mat->t; t++) //Building the schedule one time step at a time
	    {
	    	//if(self_rank==0)
	    	//printf("*******rank%d: t=%d\n", self_rank, t);
	    	//printf("rank%d: t=%d\n", self_rank, t);
	    	//Each iteration of this for loop represents
	    	//a single step of the desired final schedule.
	    	//In each iteration, we should find out everything
	    	//that needs to be performed in the given time step (t).
	    	//We have to check to see what should be done at each
	    	//step based on the cmn_nbh_mat as well as the incom_sched_mat.

	    	int k;
	    	int flag=0;

	    	  for(i = 0; i < cmn_nbh_mat->num_rows; i++)
	    	  {
	    		  for(k=0;k<K_vrbl-1;k++)
	    		  {
	    		  if(cmn_nbh_mat->comb_matrix[i][t].opt[k] != IDLE)
	    		  {
	    			  flag=1;
	    			  goto exit_checktloop;
	    		  }
	    		  }
	    	  }

	    	  exit_checktloop:
			  if(flag==1)
			  {
				  for(k=0;k<K_vrbl-1;k++)
				  {
	    			//  if(self_rank==51)
	    				//  printf("rank%d: cmn_nbh_mat->comb_matrix[%d][%d].grp_frnds[%d]=%d\n", self_rank, i, t, k, cmn_nbh_mat->comb_matrix[i][t].grp_frnds[k]);
	    			 // printf("rank%d: i=%d, t=%d\n", self_rank, i, t);
	    			//  for(k=0;k<K_vrbl-1;k++)
	    			  //{
	    			 // if(cmn_nbh_mat->comb_matrix[i][t].opt[k] == RECV && flag_recv==0)
	    			  {
	    				  //flag_recv=1;
	    				  char *rb = ((char *)exchange_tmp_buf)+ k * recvcount * recvtype_extent;

	    				  mpi_errno = MPID_Sched_recv(rb, recvcount, recvtype, cmn_nbh_mat->comb_matrix[i][t].grp_frnds[k], comm_ptr, s);
	    				  if (mpi_errno) MPIU_ERR_POP(mpi_errno);

	    				  int *test_rcv=(int *) rb;
	    				  if(self_rank==0)
	    				  {
	    					  printf("0.rank%d pid%d, iter%d: t=%d, receivebuff from source %d is %d\n", self_rank, getpid(), iter, t, cmn_nbh_mat->comb_matrix[i][t].grp_frnds[k],  recvcount);
	    				  }
	    				//  if(self_rank==0 )
	    					//  cnt_recv_ind++;



	    				//  sender_rank=self_rank;
	    			  }
	    				//  if(self_rank==0)
	    					//  printf("I received from %d at offset %d\n", cmn_nbh_mat->comb_matrix[i][t].grp_frnds[k], k);

	    			//  if(cmn_nbh_mat->comb_matrix[i][t].opt[k] == SEND && flag_send==0)
	    			  {
	    				 // flag_send=1;
	    				  mpi_errno = MPID_Sched_send(sendbuf, sendcount, sendtype, cmn_nbh_mat->comb_matrix[i][t].grp_frnds[k], comm_ptr, s);
	    				  if (mpi_errno) MPIU_ERR_POP(mpi_errno);

	    				  char *copy_from = ((char*)sendbuf);
	    				  char *copy_to = ((char*)exchange_tmp_buf+ (K_vrbl-1) * recvcount * recvtype_extent);
	    				  MPID_Sched_copy(copy_from, sendcount, sendtype, copy_to, sendcount, sendtype, s);

	    				  MPID_SCHED_BARRIER(s);
	    				//  if(cmn_nbh_mat->comb_matrix[i][t].grp_frnds[k]==0)
	    				//	  cnt_sent_ind++;


	    				//  if(cmn_nbh_mat->comb_matrix[i][t].grp_frnds[k]==2)
	    				  if(cmn_nbh_mat->comb_matrix[i][t].grp_frnds[k]==0)
	    					  printf("0.rank%d pid%d iter%d: t=%d, sendbuff to dest %d is %d\n", self_rank, getpid(), iter, t, cmn_nbh_mat->comb_matrix[i][t].grp_frnds[k],  sendcount);

	    			//	  if(cmn_nbh_mat->comb_matrix[i][t].grp_frnds[k]==0)
	    				//  printf("rank %d sent to me\n", self_rank);

	    				//  if(cmn_nbh_mat->comb_matrix[i][t].grp_frnds[k]==0)
	    					//  printf("******rank %d receive: cmn_nbh_mat->comb_matrix[%d][%d].grp_frnds[%d]=%d\n", self_rank, i, t, k, cmn_nbh_mat->comb_matrix[i][t].grp_frnds[k]);

	    			 // }

	    			//  for(k=0;k<K_vrbl-1;k++)
	    			 // {
	    			 // }



	    			  }
	    				  //MPID_SCHED_BARRIER(s);
	    		  }


	    	  MPID_SCHED_BARRIER(s);
	    	  for(i = 0; i < cmn_nbh_mat->num_rows; i++)
	    	  {
	    		  	  	//  printf("rank%d: i=%d, t=%d\n", self_rank, i, t);
	    				  if(cmn_nbh_mat->comb_matrix[i][t].opt[0] == RECV && !cmn_nbh_mat->is_row_offloaded[i])
	    				  {
	    					//  printf("rank%d: 1\n", self_rank);
	    				//  for(j = i; j < cmn_nbh_mat->num_rows; j++)  //why starting from i!?? we are in this loop only once but we want to send to all outgoings that we are responsible for
	    				 // {
	    					 // {
	    						//  printf("rank%d: j=%d, t=%d\n", self_rank, j, t);
	    						 // if(!cmn_nbh_mat->is_row_offloaded[j] && cmn_nbh_mat->comb_matrix[j][t].opt[0] == RECV)
	    						  //{
	    							//  printf("rank%d: 2\n", self_rank);
	    							  mpi_errno = MPID_Sched_send(exchange_tmp_buf, sendcount*K_vrbl,
	    									  sendtype, dests[i], comm_ptr, s);
	    							  if (mpi_errno) MPIU_ERR_POP(mpi_errno);

	    							//  if(dests[i]==0)
	    							 // if(self_rank==1)
	    							//  {
	    							//  printf("1.rank%d pid%d: t=%d, send to dest %d \n", self_rank, getpid(), t, dests[i]);
	    							 // cnt_sent_cmb++;
	    							 // }
	    							    if(dests[i]==0)
	    							  {
	    								  printf("rank%d pid%d iter %d: sent the combined message to rank %d\n", self_rank,getpid(), iter, dests[i]);
	    								 // sender_rank=self_rank;
	    							  }
	    							  //Mark the neighbor OFF as we are done with it
	    							  //cmn_nbh_mat->is_row_offloaded[k] = 1; //DO NOT set this to 1 as it changes the derived pattern for future calls!
	    						 // }
	    					 // }
	    				 // }
	    				  }


	    			  //Schedule a copy from sendbuf to exchange_tmp_buf
	    			//  printf("rank%d: %d done\n", self_rank, t);

	    	  }
	    }

	    	 // MPID_SCHED_BARRIER(s);
	    	//  if(self_rank==0)
	    	/*  {
	    		  for(i=0;i<indegree;i++)
	    		  {
	    			  int ii;
	    			  for(ii=0;ii<K_vrbl+1;ii++)
	    				  printf("rank%d: incom_sched_mat[%d][%d]=%d\n", self_rank, i, ii, incom_sched_mat[i][ii]);
	    		  }

	    	  }*/


	    	  //Checking the incom_sched_mat now
	    	  //Scheduling necessary recv operations
	    	  for(i = 0; i < indegree; i++)
	    	  {
	    		//  if(self_rank==0)
	    			//  printf("rank%d: in[%d][0]=%d\n", self_rank, i, incom_sched_mat[i][0]);
	    		  if(!incom_sched_mat[i][0]) //On incoming neighbors
	    		  {
	    			//  if(self_rank==0)
	    				//  printf("rank%d: in[%d][1]=%d, t=%d\n", self_rank, i, incom_sched_mat[i][1], t);
	    			  if(incom_sched_mat[i][1] == t)
	    			  {

	    				  //Schedule a receive from the corresponding source.
	    				  if(incom_sched_mat[i][K_vrbl+2]==2)
	    				  {
	    				  incom_recv_count = K_vrbl * recvcount;
	    				//  if(self_rank==0)
	    					//  cnt_recv_comb++;
	    				//  if(self_rank==0)
	    				 // printf("rank%d pid%d should receive combined from %d\n", self_rank, getpid(), srcs[i]);
	    				  }
	    				  else if(incom_sched_mat[i][K_vrbl+2]==1)
	    				  {
	    					  incom_recv_count = recvcount;
	    					//  if(self_rank==0)
	    						//  cnt_recv_ind++;
	    					 // printf("rank%d pid%d should receive direct from %d\n", self_rank, getpid(), srcs[i]);
	    				  }

	    				  int recv_offset=find_incom_tmp_buf_offset(incom_sched_mat, i, recvbuf_size);

	    				//  if(self_rank==0)
	    					//  printf("incoming neighbor %d is received in offset %d\n", i, recv_offset);


	    				  char *rb = ((char*)incom_tmp_buf) + recv_offset;
	    				  mpi_errno = MPID_Sched_recv(rb, incom_recv_count, recvtype, srcs[i], comm_ptr, s);
	    				  if (mpi_errno) MPIU_ERR_POP(mpi_errno);



	    				  if(self_rank==0)
	    				  printf("2.rank%d pid%d iter%d: t=%d, receivebuff from source %d is %d\n", self_rank, getpid(), iter, t, srcs[i],  incom_recv_count);

	    				//  if(self_rank==0)
	    					//  printf("rank %d should receive from %d at offset %d\n", self_rank, srcs[i], recv_offset);
	    				  //MPID_SCHED_BARRIER(s); //Can't have this here; could cause deadlock. We have to first receive all messages in a separate tmp buffer.
	    				  //break; //Commented out because it is NOT True that we should receive an incoming message from only one incoming neighbor per each step t.
	    			  }
	    		  }

	    	  }
	    	  MPID_SCHED_BARRIER(s); //Making sure the message has completely been received in tmp bufer before attempting to copy it.*/

	    	  //Scheduling message copies from incom_tmp_buf to final recv_buf
	    	  for(i = 0; i < indegree; i++)
	    	  {
	    		  if(!incom_sched_mat[i][0]) //On incoming neighbors
	    		  {
	    			  if(incom_sched_mat[i][1] == t)
	    			  {
	    				  int nbr_idx, offset1;
	    				  offset1 = find_incom_tmp_buf_offset(incom_sched_mat, i, recvbuf_size);


	    				  if(incom_sched_mat[i][2]>=0) //it means that it has the combined message for other group friends
	    				 {
	    				  for(j = 0; j < K_vrbl; j++)
	    				  {
	    					 // if(self_rank==0)
	    						//  printf("1.\n");


	    					  nbr_idx = find_in_arr(srcs, indegree, incom_sched_mat[i][j+2]);

	    					  //if(self_rank==0)
	    						//  printf("incoming neighbor %d is copied from offset %d\n", i, offset1);
	    					  if(nbr_idx>=0)
	    					  {
	    					  char *copy_from = ((char*)incom_tmp_buf) + offset1 +  j * recvcount * recvtype_extent;
	    					  char *copy_to = ((char*)recvbuf) + nbr_idx * recvcount * recvtype_extent;
	    					  MPID_Sched_copy(copy_from, recvcount, recvtype,
	    							  copy_to, recvcount, recvtype, s);

	    					//  if(self_rank==0)
	    						//  printf("copy from offset %d for src %d to %d\n", offset1 +  j* recvcount * recvtype_extent, srcs[i], nbr_idx * recvcount * recvtype_extent);
	    					  }
	    					  }

	    				  }
	    				  else
	    				  {
	    					  nbr_idx = find_in_arr(srcs, indegree, incom_sched_mat[i][K_vrbl+1]);

	    					  char *copy_from = ((char*)incom_tmp_buf) + offset1 ;
	    					  char *copy_to = ((char*)recvbuf) + nbr_idx * recvcount * recvtype_extent;
	    					  MPID_Sched_copy(copy_from, recvcount, recvtype,
	    							  copy_to, recvcount, recvtype, s);

	    					//  if(self_rank==0)
	    						//  printf("copy from offset %d for src %d to %d\n", offset1 , srcs[i], nbr_idx * recvcount * recvtype_extent);
	    				  }

	    				  //Mark the neighbor OFF as we are done with it
	    				  //incom_sched_mat[i][0] = 1; //DO NOT set this to 1 as it changes the derived pattern for future calls!
	    				  //break; //Commented out because it is NOT True that we should receive an incoming message from only one incoming neighbor per each step t.
	    			  }
	    		  }
	    	  }
	    	  MPID_SCHED_BARRIER(s);
	    }

	  //  printf("rank%d, pid%d: 2.\n", self_rank, getpid());
	  //  printf("rank%d: 2.\n", getpid());
	    /* Now we might have some ON neighbors that we
	         * could not plan to send combined messages to
	         * them. That is, we have to communicate with
	         * them one by one in a simple naive way.
	         */
	        MPID_SCHED_BARRIER(s);//This is needed here because otherwise, a send with a combined message issued above, could be received by a recv operation issued below. This could happen because we have that 'break' after issuing outgoing sends.
	        for(i = 0; i < cmn_nbh_mat->num_rows; i++)
	        {
	            //if(!cmn_nbh_mat->is_row_offloaded[i])
	            if(!cmn_nbh_mat->ignore_row[i]) //out neighbor still active
	            {
	                mpi_errno = MPID_Sched_send(sendbuf, sendcount, sendtype,
	                                            dests[i], comm_ptr, s);
	                if (mpi_errno) MPIU_ERR_POP(mpi_errno);

	               // if(dests[i]==0)
	                //	   cnt_sent_ind++;
	               // if(dests[i]==4)
	                if(dests[i]==0)
	                printf("3.rank%d pid%d iter%d: sends directly to dest %d \n", self_rank, getpid(), iter, dests[i]);
	                //if(dests[i]==0)
	                	//printf("rank%d pid%d send directly to rank %d\n", self_rank, getpid(), dests[i]);
	                //Mark the neighbor OFF as we are done with it
	    			//cmn_nbh_mat->is_row_offloaded[i] = 1; //DO NOT set this to 1 as it changes the derived pattern for future calls!
	            }
	        }
	     //   printf("rank%d: 3.\n", getpid());
	      //  printf("rank%d, pid%d: 3.\n", self_rank, getpid());

	    for(i = 0; i < indegree; i++)
	    {
	    	if(!incom_sched_mat[i][0] && incom_sched_mat[i][1] >= cmn_nbh_mat->t) //On incoming neighbors not covered before
	    	{
	    		//Schedule a receive from the corresponding source.

	    		if(incom_sched_mat[i][K_vrbl+2]==2)
	    		{
	    		incom_recv_count =  K_vrbl * recvcount;
	    	//	if(self_rank==0)
	    		//	cnt_recv_comb++;
	    	//	if(self_rank==1)
	    		// printf("***rank%d pid%d should receive combined from %d\n", self_rank, getpid(), srcs[i]);
	    		}
	    		else if(incom_sched_mat[i][K_vrbl+2]==1)
	    		{
	    			incom_recv_count = recvcount;
	    			//if(self_rank==0)
	    				//cnt_recv_ind++;
	    		//	if(self_rank==0)
	    			// printf("*rank%d pid%d should receive direct from %d\n", self_rank, getpid(), srcs[i]);
	    		}

	    		char *rb = ((char*)incom_tmp_buf) + find_incom_tmp_buf_offset(incom_sched_mat, i, recvbuf_size);
	    		mpi_errno = MPID_Sched_recv(rb, incom_recv_count, recvtype, srcs[i], comm_ptr, s);
	    		if (mpi_errno) MPIU_ERR_POP(mpi_errno);


	    		if(self_rank==0)
	    		printf("4.rank%d pid%d iter%d: receivebuff from source %d is %d\n", self_rank, getpid(), iter, srcs[i],  incom_recv_count);
	    		//if(self_rank==0)
	    		 //printf("*rank %d should receive from %d at offset %d\n", self_rank, srcs[i], find_incom_tmp_buf_offset(incom_sched_mat, i, recvbuf_size));
	    		//MPID_SCHED_BARRIER(s); //Can't have this here; could cause deadlock. We have to first receive all messages in a separate tmp buffer.
	    		//break; //Commented out because it is NOT True that we should receive an incoming message from only one incoming neighbor per each step t.
	    	}
	    }
	    MPID_SCHED_BARRIER(s);

	  //  printf("rank%d, pid%d: 4.\n", self_rank, getpid());
	  //  printf("rank%d: 4.\n", getpid());
	    //Scheduling message copies from incom_tmp_buf to final recv_buf
	    for(i = 0; i < indegree; i++)
	    {
	    	if(!incom_sched_mat[i][0] && incom_sched_mat[i][1] >= cmn_nbh_mat->t) //On incoming neighbors not covered before
	    	{
	    		int nbr_idx, offset1;
	    		offset1 = find_incom_tmp_buf_offset(incom_sched_mat, i, recvbuf_size);

	    		if(incom_sched_mat[i][2]>=0) //it means that it has the combined message for other group friends
	    		{
	    		for(j = 0; j < K_vrbl; j++)
	    		{
	    			//A completely different (potentially recursive) approach is needed for the cumulative case
	    			//if(self_rank==0)
	    				// printf("2.\n");
	    			nbr_idx = find_in_arr(srcs, indegree, incom_sched_mat[i][j + 2]);

	    			if(nbr_idx>=0)
	    			{
	    			char *copy_from = ((char*)incom_tmp_buf) + offset1 + j * recvcount * recvtype_extent;
	    			char *copy_to = ((char*)recvbuf) + nbr_idx * recvcount * recvtype_extent;
	    			MPID_Sched_copy(copy_from, recvcount, recvtype,
	    					copy_to, recvcount, recvtype, s);

	    			//if(self_rank==0)
	    				//printf("*copy from offset %d for src %d to %d\n", offset1 + j * recvcount * recvtype_extent, srcs[i], nbr_idx * recvcount * recvtype_extent);
	    			}
	    		}
	    		}
	    		else
	    		{
	    			nbr_idx = find_in_arr(srcs, indegree, incom_sched_mat[i][K_vrbl+1]);

	    		//	if(self_rank==0)
	    			//	printf("incom_sched_mat[%d][%d]=%d, offset=%d, nbr_idx=%d \n", i,K_vrbl+1, incom_sched_mat[i][K_vrbl+1], offset1, nbr_idx );

	    			char *copy_from = ((char*)incom_tmp_buf) + offset1 ;
	    			char *copy_to = ((char*)recvbuf) + nbr_idx * recvcount * recvtype_extent;
	    			MPID_Sched_copy(copy_from, recvcount, recvtype,
	    					copy_to, recvcount, recvtype, s);


	    			//if(self_rank==0)
	    				//printf("*copy direct from offset %d for src %d to %d\n", offset1 , srcs[i], nbr_idx * recvcount * recvtype_extent);
	    		}

	    		//Mark the neighbor OFF as we are done with it
	    		//incom_sched_mat[i][0] = 1; //DO NOT set this to 1 as it changes the derived pattern for future calls!
	    		//break; //Commented out because it is NOT True that we should receive an incoming message from only one incoming neighbor per each step t.
	    	}
	    }

	    int *data= (int*) recvbuf;
	    if(self_rank==0 && iter==46)
	    {
	    	printf("***rank 0---PID%d received:\n", getpid());
	    	for(i=0;i<8;i++)
	    		printf("%d  ", data[i]);
	    	printf("\n");
	    }


	  //  printf("rank%d, pid%d: 5.\n", self_rank, getpid());
	    MPID_SCHED_BARRIER(s);

	    MPIR_SCHED_CHKPMEM_COMMIT(s);

	 /*   int sum_comb, sum_ind;
	    MPI_Reduce(&cnt_sent_cmb, &sum_comb, 1, MPI_INT, MPI_SUM, 0,
	               MPI_COMM_WORLD);
	    MPI_Reduce(&cnt_sent_ind, &sum_ind, 1, MPI_INT, MPI_SUM, 0,
	   	               MPI_COMM_WORLD);

	    if(self_rank==0)
	    	if(sum_comb!=cnt_recv_comb || sum_ind!=cnt_recv_ind)
	    		printf("ATTENTION!  sent comb=%d, sent ind=%d, recved comb=%d, recved ind=%d \n", sum_comb, sum_ind, cnt_recv_comb, cnt_recv_ind);*/


	 /*   	printf("rank%d: sent comb=%d, sent ind=%d\n", self_rank, sum_comb, sum_ind);

	    if(self_rank==0)
	    	printf("rank%d: recved comb=%d, recved ind=%d \n", self_rank, cnt_recv_comb, cnt_recv_ind);*/

	fn_exit:
	    return mpi_errno;
	fn_fail:
	MPIR_SCHED_CHKPMEM_REAP(s);
	goto fn_exit;
}


/* any non-MPI functions go here, especially non-static ones */

#undef FUNCNAME
#define FUNCNAME MPIR_Ineighbor_allgather_default
#undef FCNAME
#define FCNAME MPIU_QUOTE(FUNCNAME)
int MPIR_Ineighbor_allgather_default(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, MPID_Comm *comm_ptr, MPID_Sched_t s, int iter)
{

	int mpi_errno = MPI_SUCCESS;
	int indegree, outdegree, weighted;
	int k,l;
	int *srcs, *dsts;
	MPI_Aint recvtype_extent;
	//SHM commented out
	//MPIU_CHKLMEM_DECL(2);

	MPID_Datatype_get_extent_macro(recvtype, recvtype_extent);

	/* This is the largest offset we add to recvbuf */
	    MPID_Ensure_Aint_fits_in_pointer(MPI_VOID_PTR_CAST_TO_MPI_AINT recvbuf +
	                                     (comm_ptr->local_size * recvcount * recvtype_extent));

	    //SHM commented out -- Just to be fair to the default approach
	    /*
	    mpi_errno = MPIR_Topo_canon_nhb_count(comm_ptr, &indegree, &outdegree, &weighted);
	    if (mpi_errno) MPIU_ERR_POP(mpi_errno);
	    MPIU_CHKLMEM_MALLOC(srcs, int *, indegree*sizeof(int), mpi_errno, "srcs");
	    MPIU_CHKLMEM_MALLOC(dsts, int *, outdegree*sizeof(int), mpi_errno, "dsts");
	    mpi_errno = MPIR_Topo_canon_nhb(comm_ptr,
	                                    indegree, srcs, MPI_UNWEIGHTED,
	                                    outdegree, dsts, MPI_UNWEIGHTED);
	    if (mpi_errno) MPIU_ERR_POP(mpi_errno);
	    */

	    //SHM
	    MPIR_Topology *topo_ptr = NULL;
	    topo_ptr = MPIR_Topology_get(comm_ptr);
	    indegree = topo_ptr->topo.dist_graph.indegree;
	    outdegree = topo_ptr->topo.dist_graph.outdegree;
	    srcs = topo_ptr->topo.dist_graph.in;
	    dsts = topo_ptr->topo.dist_graph.out;

	    for (k = 0; k < outdegree; ++k) {
	        mpi_errno = MPID_Sched_send(sendbuf, sendcount, sendtype, dsts[k], comm_ptr, s);
	        if (mpi_errno) MPIU_ERR_POP(mpi_errno);
	    }

	    for (l = 0; l < indegree; ++l) {
	        char *rb = ((char *)recvbuf) + l * recvcount * recvtype_extent;
	        mpi_errno = MPID_Sched_recv(rb, recvcount, recvtype, srcs[l], comm_ptr, s);
	        if (mpi_errno) MPIU_ERR_POP(mpi_errno);
	    }

	    MPID_SCHED_BARRIER(s);

	    int i, self_rank;
	    MPI_Comm_rank(MPI_COMM_WORLD, &self_rank);

	    int *data= (int*) recvbuf;
	    if(self_rank==0 && iter==46)
	    {
	    	printf("###rank 0---PID%d received:\n", getpid());
	    	for(i=0;i<8;i++)
	    		printf("%d  ", data[i]);
	    	printf("\n");
	    }

	fn_exit:
	    //SHM commented out
	    //MPIU_CHKLMEM_FREEALL();
	    return mpi_errno;
	fn_fail:
	    goto fn_exit;

  /*  int mpi_errno = MPI_SUCCESS;
    int indegree, outdegree, weighted;
    int k,l;
    int *srcs, *dsts;
    MPI_Aint recvtype_extent;
    MPIU_CHKLMEM_DECL(2);

    MPID_Datatype_get_extent_macro(recvtype, recvtype_extent);

    /* This is the largest offset we add to recvbuf
    MPID_Ensure_Aint_fits_in_pointer(MPI_VOID_PTR_CAST_TO_MPI_AINT recvbuf +
                                     (comm_ptr->local_size * recvcount * recvtype_extent));

    mpi_errno = MPIR_Topo_canon_nhb_count(comm_ptr, &indegree, &outdegree, &weighted);
    if (mpi_errno) MPIU_ERR_POP(mpi_errno);
    MPIU_CHKLMEM_MALLOC(srcs, int *, indegree*sizeof(int), mpi_errno, "srcs");
    MPIU_CHKLMEM_MALLOC(dsts, int *, outdegree*sizeof(int), mpi_errno, "dsts");
    mpi_errno = MPIR_Topo_canon_nhb(comm_ptr,
                                    indegree, srcs, MPI_UNWEIGHTED,
                                    outdegree, dsts, MPI_UNWEIGHTED);
    if (mpi_errno) MPIU_ERR_POP(mpi_errno);

    for (k = 0; k < outdegree; ++k) {
        mpi_errno = MPID_Sched_send(sendbuf, sendcount, sendtype, dsts[k], comm_ptr, s);
        if (mpi_errno) MPIU_ERR_POP(mpi_errno);
    }

    for (l = 0; l < indegree; ++l) {
        char *rb = ((char *)recvbuf) + l * recvcount * recvtype_extent;
        mpi_errno = MPID_Sched_recv(rb, recvcount, recvtype, srcs[l], comm_ptr, s);
        if (mpi_errno) MPIU_ERR_POP(mpi_errno);
    }

    MPID_SCHED_BARRIER(s);

fn_exit:
    MPIU_CHKLMEM_FREEALL();
    return mpi_errno;
fn_fail:
    goto fn_exit;*/
}

#undef FUNCNAME
#define FUNCNAME MPIR_Ineighbor_allgather_impl
#undef FCNAME
#define FCNAME MPIU_QUOTE(FUNCNAME)
int MPIR_Ineighbor_allgather_impl(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, MPID_Comm *comm_ptr, MPI_Request *request, int iter)
{
    int mpi_errno = MPI_SUCCESS;
    int tag = -1;
    MPID_Request *reqp = NULL;
    MPID_Sched_t s = MPID_SCHED_NULL;

    *request = MPI_REQUEST_NULL;

    mpi_errno = MPID_Sched_next_tag(comm_ptr, &tag);
    if (mpi_errno) MPIU_ERR_POP(mpi_errno);
    mpi_errno = MPID_Sched_create(&s);
    if (mpi_errno) MPIU_ERR_POP(mpi_errno);
    MPIU_Assert(comm_ptr->coll_fns != NULL);
    MPIU_Assert(comm_ptr->coll_fns->Ineighbor_allgather != NULL);
  //  mpi_errno = comm_ptr->coll_fns->Ineighbor_allgather(sendbuf, sendcount, sendtype,
    //                                                    recvbuf, recvcount, recvtype,
      //                                                  comm_ptr, s);

    if(nbr_impl == 0)
    {
    mpi_errno = MPIR_Ineighbor_allgather_default(sendbuf, sendcount, sendtype,
    	                                                        recvbuf, recvcount, recvtype,
    	                                                        comm_ptr, s, iter);
    }
    else
    {
    	 mpi_errno = MPIR_Ineighbor_allgather_SMGM(sendbuf, sendcount, sendtype,
    	    	                                                        recvbuf, recvcount, recvtype,
    	    	                                                        comm_ptr, s, iter);
    }



    //SMGM
   /*    if(nbr_impl == 1) //Use SMGM implementation -- This is not the best place for this statement
       {
       	comm_ptr->coll_fns->Ineighbor_allgather = MPIR_Ineighbor_allgather_SMGM;
       }
       else
       {
       	comm_ptr->coll_fns->Ineighbor_allgather = MPIR_Ineighbor_allgather_default;
       }

       mpi_errno = comm_ptr->coll_fns->Ineighbor_allgather(sendbuf, sendcount, sendtype,
                                                           recvbuf, recvcount, recvtype,

                                                         comm_ptr, s);*/
  //  MPID_SCHED_BARRIER(s);

    if (mpi_errno) MPIU_ERR_POP(mpi_errno);

    mpi_errno = MPID_Sched_start(&s, comm_ptr, tag, &reqp);
    if (reqp)
        *request = reqp->handle;
    if (mpi_errno) MPIU_ERR_POP(mpi_errno);



fn_exit:
    return mpi_errno;
fn_fail:
    goto fn_exit;
}

#endif /* MPICH_MPI_FROM_PMPI */

#undef FUNCNAME
#define FUNCNAME MPI_Ineighbor_allgather
#undef FCNAME
#define FCNAME MPIU_QUOTE(FUNCNAME)
/*@
MPI_Ineighbor_allgather - Nonblocking version of MPI_Neighbor_allgather.

Input Parameters:
+ sendbuf - starting address of the send buffer (choice)
. sendcount - number of elements sent to each neighbor (non-negative integer)
. sendtype - data type of send buffer elements (handle)
. recvcount - number of elements received from each neighbor (non-negative integer)
. recvtype - data type of receive buffer elements (handle)
- comm - communicator (handle)

Output Parameters:
+ recvbuf - starting address of the receive buffer (choice)
- request - communication request (handle)

.N ThreadSafe

.N Fortran

.N Errors
@*/
int MPI_Ineighbor_allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request, int iter)
{
    int mpi_errno = MPI_SUCCESS;
    MPID_Comm *comm_ptr = NULL;
    MPID_MPI_STATE_DECL(MPID_STATE_MPI_INEIGHBOR_ALLGATHER);

    MPIU_THREAD_CS_ENTER(ALLFUNC,);
    MPID_MPI_FUNC_ENTER(MPID_STATE_MPI_INEIGHBOR_ALLGATHER);

    /* Validate parameters, especially handles needing to be converted */
#   ifdef HAVE_ERROR_CHECKING
    {
        MPID_BEGIN_ERROR_CHECKS
        {
            MPIR_ERRTEST_DATATYPE(sendtype, "sendtype", mpi_errno);
            MPIR_ERRTEST_DATATYPE(recvtype, "recvtype", mpi_errno);
            MPIR_ERRTEST_COMM(comm, mpi_errno);

            /* TODO more checks may be appropriate */
        }
        MPID_END_ERROR_CHECKS
    }
#   endif /* HAVE_ERROR_CHECKING */

    /* Convert MPI object handles to object pointers */
    MPID_Comm_get_ptr(comm, comm_ptr);

    /* Validate parameters and objects (post conversion) */
#   ifdef HAVE_ERROR_CHECKING
    {
        MPID_BEGIN_ERROR_CHECKS
        {
            if (HANDLE_GET_KIND(sendtype) != HANDLE_KIND_BUILTIN) {
                MPID_Datatype *sendtype_ptr = NULL;
                MPID_Datatype_get_ptr(sendtype, sendtype_ptr);
                MPID_Datatype_valid_ptr(sendtype_ptr, mpi_errno);
                if (mpi_errno != MPI_SUCCESS) goto fn_fail;
                MPID_Datatype_committed_ptr(sendtype_ptr, mpi_errno);
                if (mpi_errno != MPI_SUCCESS) goto fn_fail;
            }

            if (HANDLE_GET_KIND(recvtype) != HANDLE_KIND_BUILTIN) {
                MPID_Datatype *recvtype_ptr = NULL;
                MPID_Datatype_get_ptr(recvtype, recvtype_ptr);
                MPID_Datatype_valid_ptr(recvtype_ptr, mpi_errno);
                if (mpi_errno != MPI_SUCCESS) goto fn_fail;
                MPID_Datatype_committed_ptr(recvtype_ptr, mpi_errno);
                if (mpi_errno != MPI_SUCCESS) goto fn_fail;
            }

            MPID_Comm_valid_ptr( comm_ptr, mpi_errno, FALSE );
            if (mpi_errno != MPI_SUCCESS) goto fn_fail;
            MPIR_ERRTEST_ARGNULL(request, "request", mpi_errno);
            /* TODO more checks may be appropriate (counts, in_place, buffer aliasing, etc) */
        }
        MPID_END_ERROR_CHECKS
    }
#   endif /* HAVE_ERROR_CHECKING */

    /* ... body of routine ...  */

    mpi_errno = MPIR_Ineighbor_allgather_impl(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm_ptr, request, iter);
    if (mpi_errno) MPIU_ERR_POP(mpi_errno);


    /* ... end of body of routine ... */

fn_exit:
    MPID_MPI_FUNC_EXIT(MPID_STATE_MPI_INEIGHBOR_ALLGATHER);
    MPIU_THREAD_CS_EXIT(ALLFUNC,);
    return mpi_errno;

fn_fail:
    /* --BEGIN ERROR HANDLING-- */
#   ifdef HAVE_ERROR_CHECKING
    {
        mpi_errno = MPIR_Err_create_code(
            mpi_errno, MPIR_ERR_RECOVERABLE, FCNAME, __LINE__, MPI_ERR_OTHER,
            "**mpi_ineighbor_allgather", "**mpi_ineighbor_allgather %p %d %D %p %d %D %C %p", sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, request);
    }
#   endif
    mpi_errno = MPIR_Err_return_comm(NULL, FCNAME, mpi_errno);
    goto fn_exit;
    /* --END ERROR HANDLING-- */
}
