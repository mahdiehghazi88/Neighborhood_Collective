/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */
/*
 *
 *  (C) 2009 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include "mpiimpl.h"
#include "topo.h"
#include <math.h>
#include <stddef.h>

#define MAX_COMB_DEGREE 50

/* -- Begin Profiling Symbol Block for routine MPI_Dist_graph_create_adjacent */
#if defined(HAVE_PRAGMA_WEAK)
#pragma weak MPI_Dist_graph_create_adjacent = PMPI_Dist_graph_create_adjacent
#elif defined(HAVE_PRAGMA_HP_SEC_DEF)
#pragma _HP_SECONDARY_DEF PMPI_Dist_graph_create_adjacent  MPI_Dist_graph_create_adjacent
#elif defined(HAVE_PRAGMA_CRI_DUP)
#pragma _CRI duplicate MPI_Dist_graph_create_adjacent as PMPI_Dist_graph_create_adjacent
#elif defined(HAVE_WEAK_ATTRIBUTE)
int MPI_Dist_graph_create_adjacent(MPI_Comm comm_old, int indegree, const int sources[],
                                   const int sourceweights[], int outdegree,
                                   const int destinations[], const int destweights[],
                                   MPI_Info info, int reorder, MPI_Comm *comm_dist_graph) __attribute__((weak,alias("PMPI_Dist_graph_create_adjacent")));
#endif
/* -- End Profiling Symbol Block */

/* Define MPICH_MPI_FROM_PMPI if weak symbols are not supported to build
   the MPI routines */
#ifndef MPICH_MPI_FROM_PMPI
#undef MPI_Dist_graph_create_adjacent
#define MPI_Dist_graph_create_adjacent PMPI_Dist_graph_create_adjacent
/* any utility functions should go here, usually prefixed with PMPI_LOCAL to
 * correctly handle weak symbols and the profiling interface */
#endif

int comp (const void * elem1, const void * elem2) {
    int f = *((int*)elem1);
    int s = *((int*)elem2);
    if (f > s) return  1;
    if (f < s) return -1;
    return 0;
}

#undef FUNCNAME
#define FUNCNAME MPIR_SMGM_Find_Group_Friend
#undef FCNAME
#define FCNAME MPIU_QUOTE(FUNCNAME)
int MPIR_SMGM_Find_Group_Friend(MPID_Comm *comm_ptr, Group_Friendship_Matrix *grp_frnd_mat, Final_Group_Friends *slcted_grp_frnd_ptr, int step) {
    int mpi_errno = MPI_SUCCESS;
    int i, j,p,n, comm_size, self_rank, indegree,outdegree;
    comm_size = comm_ptr->local_size;
    self_rank = comm_ptr->rank;
    MPI_Status status;
    MPI_Status status2;

    int context_offset = (comm_ptr->comm_kind == MPID_INTRACOMM) ?
	                      MPID_CONTEXT_INTRA_COLL : MPID_CONTEXT_INTER_COLL;

    MPID_Request *req_ptr = NULL;
    MPID_Request *req_ptr2 = NULL;

    MPIU_CHKPMEM_DECL(9);

    MPIR_Topology *topo_ptr = NULL;
    topo_ptr = MPIR_Topology_get(comm_ptr);
    if(topo_ptr == NULL) {
        fprintf(stderr, "ERROR: Communicator topology pointer is NULL!\n");
	goto fn_fail;
    }

    Final_Group_Friends slcted_grp_frnd;
    MPIU_CHKPMEM_CALLOC(slcted_grp_frnd.grp_frnds, int*, K_vrbl* sizeof(int), mpi_errno, "slcted_grp_frnd.grp_frnds");
    slcted_grp_frnd.is_grp=0;

    outdegree = topo_ptr->topo.dist_graph.outdegree;

    if(grp_frnd_mat->num_grp_frnds!=0) {
        Transfered_Request sndreq;
	MPIU_CHKPMEM_CALLOC(sndreq.grp, int*, K_vrbl* sizeof(int), mpi_errno, "sndreq.grp");

	F_Group_Friends *avail_grpfrnds;
	MPIU_CHKPMEM_MALLOC(avail_grpfrnds, F_Group_Friends*, sizeof(F_Group_Friends), mpi_errno, "avail_grpfrnds");
	MPIU_CHKPMEM_MALLOC(avail_grpfrnds->grpfriends, int**, K_vrbl * sizeof(int*), mpi_errno, "avail_grpfrnds->grpfriends");

	avail_grpfrnds->num_avail_grpfrnds=  grp_frnd_mat->num_grp_frnds;

	MPIU_CHKPMEM_CALLOC(avail_grpfrnds->num_cmn_brs, int*, avail_grpfrnds->num_avail_grpfrnds * sizeof(int), mpi_errno, "avail_grpfrnds->num_cmn_brs");


	for(i=0;i<K_vrbl;i++) {
	    avail_grpfrnds->grpfriends[i] = MPIU_Malloc(avail_grpfrnds->num_avail_grpfrnds * sizeof(int));
	}

	Chosen_Group_Friends *chosen_grpfrnds;
	MPIU_CHKPMEM_MALLOC(chosen_grpfrnds, Chosen_Group_Friends*, sizeof(Chosen_Group_Friends), mpi_errno, "chosen_grpfrnds");
	MPIU_CHKPMEM_MALLOC(chosen_grpfrnds->grpfriends, int**, K_vrbl * sizeof(int*), mpi_errno, "chosen_grpfrnds->grpfriends");
	MPIU_CHKPMEM_MALLOC(chosen_grpfrnds->source, int*, grp_frnd_mat->num_grp_frnds*K_vrbl * sizeof(int*), mpi_errno, "chosen_grpfrnds->source");


	chosen_grpfrnds->num_chosen_grpfrnds=  0;
	for(i=0;i<K_vrbl;i++) {
	    chosen_grpfrnds->grpfriends[i] = MPIU_Malloc( grp_frnd_mat->num_grp_frnds*K_vrbl* sizeof(int));
	}

	p=0;
	for(j=0;j<grp_frnd_mat->Total_num_grp_frnds;j++) {
	    if(grp_frnd_mat->is_active_grpfrnd[j]) {
	        for(i=0;i<K_vrbl-1;i++) {
                    avail_grpfrnds->grpfriends[i][p]=grp_frnd_mat->grpfriends[i][j];
		}
		avail_grpfrnds->num_cmn_brs[p]=grp_frnd_mat->num_cmn_brs[j];
		p++;
	    }
	}

	if(p!=grp_frnd_mat->num_grp_frnds) printf("rank%d: ERROR: p=%d, num_grp_frnds=%d \n", self_rank, p, grp_frnd_mat->num_grp_frnds);

	for(j=0;j<avail_grpfrnds->num_avail_grpfrnds;j++) avail_grpfrnds->grpfriends[K_vrbl-1][j]=self_rank;


	int max_cmn_nbr_index, max_val=0;
	for(i=0;i<grp_frnd_mat->num_grp_frnds;i++) {
	    if(avail_grpfrnds->num_cmn_brs[i]>max_val) {
	        max_cmn_nbr_index=i;
		max_val=avail_grpfrnds->num_cmn_brs[i];
	    }
        }

	//add the source rank to the group friends
	int order_flag=0;
	for(i=0;i<K_vrbl;i++) {
	    if(self_rank>avail_grpfrnds->grpfriends[i][max_cmn_nbr_index] && order_flag==0) {
	        sndreq.grp[i]=avail_grpfrnds->grpfriends[i][max_cmn_nbr_index];
	    } else {
	        if(order_flag==0) {
		    sndreq.grp[i]=self_rank;
		    order_flag=1;
		} else {
		    sndreq.grp[i]=avail_grpfrnds->grpfriends[i-1][max_cmn_nbr_index];
		}
	    }
	}
	sndreq.msgtype=1;
	for(j=0;j<K_vrbl;j++) {
	    if(sndreq.grp[j]!=self_rank) {
	        mpi_errno= MPID_Send(&sndreq.msgtype, 1, MPI_INT, sndreq.grp[j], 123+step, comm_ptr, context_offset, &req_ptr);
		if (mpi_errno) MPIU_ERR_POP(mpi_errno);

		mpi_errno= MPID_Send(sndreq.grp, K_vrbl, MPI_INT, sndreq.grp[j], 400+step, comm_ptr, context_offset, &req_ptr2);
		if (mpi_errno) MPIU_ERR_POP(mpi_errno);
	    }
	}

	while(1) {
	    Transfered_Request recvreq;
	    MPIU_CHKPMEM_CALLOC(recvreq.grp, int*, K_vrbl* sizeof(int), mpi_errno, "recvreq.grp");

	    MPID_Recv(&recvreq.msgtype, 1, MPI_INT, MPI_ANY_SOURCE, 123+step, comm_ptr, context_offset, &status, &req_ptr);

	    if(recvreq.msgtype==1) {
	        MPID_Recv(recvreq.grp, K_vrbl, MPI_INT, MPI_ANY_SOURCE, 400+step, comm_ptr, context_offset, &status2, &req_ptr2);
		chosen_grpfrnds->source[chosen_grpfrnds->num_chosen_grpfrnds]=status2.MPI_SOURCE;
		for(i=0;i<K_vrbl;i++) {
		    chosen_grpfrnds->grpfriends[i][chosen_grpfrnds->num_chosen_grpfrnds]=recvreq.grp[i];
		}
		    chosen_grpfrnds->num_chosen_grpfrnds++;
	    } else {
                for(j=0;j<grp_frnd_mat->num_grp_frnds;j++) {
		    for(i=0;i<K_vrbl;i++) {
		        if(status.MPI_SOURCE==avail_grpfrnds->grpfriends[i][j]) {
			    avail_grpfrnds->num_cmn_brs[j]=-1;
			    avail_grpfrnds->num_avail_grpfrnds--;
			    i=K_vrbl;
			}
		    }
		}

		int cnt=0;
		for(i=0;i<K_vrbl;i++) {
			if(status.MPI_SOURCE==sndreq.grp[i]) cnt=1;
		}
		if(cnt==1) {//choose a new set of group friends and send them request
		    max_cmn_nbr_index=-1, max_val=0;
		    for(i=0;i<grp_frnd_mat->num_grp_frnds;i++)
		    if(avail_grpfrnds->num_cmn_brs[i]>max_val) {
		        max_cmn_nbr_index=i;
			max_val=avail_grpfrnds->num_cmn_brs[i];
		    }
		    if(max_cmn_nbr_index!=-1) {
		        //add the source rank to the group friends
			int order_flag=0;
			for(i=0;i<K_vrbl;i++) {
			    if(self_rank>avail_grpfrnds->grpfriends[i][max_cmn_nbr_index] && order_flag==0) {
			        sndreq.grp[i]=avail_grpfrnds->grpfriends[i][max_cmn_nbr_index];
			    } else {
				if(order_flag==0) {
				    sndreq.grp[i]=self_rank;
				    order_flag=1;
				} else {
				    sndreq.grp[i]=avail_grpfrnds->grpfriends[i-1][max_cmn_nbr_index];
				}
			    }
			}

			sndreq.msgtype=1;
			for(j=0;j<K_vrbl;j++) {
			    if(sndreq.grp[j]!=self_rank) {
			        MPID_Send(&sndreq.msgtype, 1, MPI_INT, sndreq.grp[j], 123+step, comm_ptr, context_offset, &req_ptr);
				MPID_Send(sndreq.grp, K_vrbl, MPI_INT, sndreq.grp[j], 400+step, comm_ptr, context_offset, &req_ptr2);
			    }
			}
		   } else { // no group friend is found so exit the loop
					goto 	exit_grouping;
		   }
		}
	    }

	    int cnt_R=0;
	    for(i=0;i<K_vrbl;i++) {
	        if(sndreq.grp[i]!=self_rank) {
		    for(j=0;j<chosen_grpfrnds->num_chosen_grpfrnds;j++) {
		        if(chosen_grpfrnds->source[j]==sndreq.grp[i]) {
			    int cnt=0;
			    for(p=0;p<K_vrbl;p++) {
			        if(chosen_grpfrnds->grpfriends[p][j]==sndreq.grp[p]) cnt++;
				}
				if(cnt==K_vrbl) {
                                    j=chosen_grpfrnds->num_chosen_grpfrnds;
				    cnt_R++;
				}
			}
		    }
		}
	    }

	    if(cnt_R==K_vrbl-1) { //we found the group friend, so save it in slcted_grp_frnd and drop the other friends 
	        for(i=0;i<K_vrbl;i++) {
		    slcted_grp_frnd.grp_frnds[i]=sndreq.grp[i];
		    slcted_grp_frnd.is_grp=1;
		}
		int *flag_ranks;
		MPIU_CHKPMEM_CALLOC(flag_ranks, int*, comm_size* sizeof(int), mpi_errno, "flag_ranks");
		for(i=0;i<comm_size;i++) {
		    flag_ranks[i]=0;
		}

		for(i=0;i<K_vrbl;i++) {    //don't drop the ranks that are chosen 
		    flag_ranks[sndreq.grp[i]]=1;
		}

		for(i=0;i<grp_frnd_mat->num_grp_frnds;i++) {
		    if(avail_grpfrnds->num_cmn_brs[i]!=-1) {
		        for(p=0;p<K_vrbl;p++) {
			    if(flag_ranks[avail_grpfrnds->grpfriends[p][i]]==0) {
			        flag_ranks[avail_grpfrnds->grpfriends[p][i]]=1;
				sndreq.msgtype=0;
				MPID_Send(&sndreq.msgtype, 1, MPI_INT, avail_grpfrnds->grpfriends[p][i], 123+step, comm_ptr, context_offset, &req_ptr);
			    }
			}
		    }
		}
		avail_grpfrnds->num_avail_grpfrnds=0;
		goto exit_grouping;
	    }
	}

	exit_grouping:
	*slcted_grp_frnd_ptr=slcted_grp_frnd;
	}

	fn_fail:
		return 0;

}


#undef FUNCNAME
#define FUNCNAME MPIR_Make_GroupFriend_Matrix
#undef FCNAME
#define FCNAME MPIU_QUOTE(FUNCNAME)
int MPIR_Make_GroupFriend_Matrix(MPID_Comm *comm_ptr, Individual_Friendship_Matrix *ind_frnd_mat, Group_Friendship_Matrix **grp_frnd_mat_ptr) {
	int mpi_errno = MPI_SUCCESS;
	int i, j,p,m,n,o,t,u,i1, j1,p1,m1,n1,o1,t1,u1,u2, comm_size, self_rank, indegree,outdegree;
	comm_size = comm_ptr->local_size;
	self_rank = comm_ptr->rank;
	int cnt_cmn_nbrs;
	MPIU_CHKPMEM_DECL(9);

	MPIR_Topology *topo_ptr = NULL;
	topo_ptr = MPIR_Topology_get(comm_ptr);
	if(topo_ptr == NULL) {
	    fprintf(stderr, "ERROR: Communicator topology pointer is NULL!\n");
	    goto fn_fail;
	}

	outdegree = topo_ptr->topo.dist_graph.outdegree;

	Group_Friendship_Matrix *grp_frnd_mat;

	MPIU_CHKPMEM_MALLOC(grp_frnd_mat, Group_Friendship_Matrix*, sizeof(Group_Friendship_Matrix), mpi_errno, "grp_frnd_mat");

	MPIU_CHKPMEM_CALLOC(grp_frnd_mat->grpfrnd_bit_map, int**, outdegree * sizeof(int*), mpi_errno, "grp_frnd_mat->grpfrnd_bit_map");
	MPIU_CHKPMEM_MALLOC(grp_frnd_mat->grpfriends, int**, K_vrbl * sizeof(int*), mpi_errno, "grp_frnd_mat->grpfrnds");
	MPIU_CHKPMEM_MALLOC(grp_frnd_mat->is_active_nbr, int*, outdegree * sizeof(int), mpi_errno, "grp_frnd_mat->is_active_nbr");
	MPIU_CHKPMEM_MALLOC(grp_frnd_mat->num_actv_grp_of_nbr, int*, outdegree * sizeof(int), mpi_errno, "grp_frnd_mat->num_actv_grp_of_nbr");

	for(i=0;i<outdegree;i++) {
	    grp_frnd_mat->is_active_nbr[i]=1;
	    grp_frnd_mat->num_actv_grp_of_nbr[i]=0;
	}

	grp_frnd_mat->num_grp_frnds=0;

	int counter=0;

	if(K_vrbl==2) {
	for(i=0;i<ind_frnd_mat->numfriends;i++) {
		cnt_cmn_nbrs=0;
		for(p=0;p<outdegree;p++) {
		    if(ind_frnd_mat->frndship_bit_map[i][p]) {
		        cnt_cmn_nbrs++;
		    }
		}
		if(cnt_cmn_nbrs>=nbr_frndshp_thr) {
		    grp_frnd_mat->num_grp_frnds++;
		}
	}

	if(grp_frnd_mat->num_grp_frnds>200000) {
	    printf("Not Enough Memory. Increase the Threshold\n");
	    return -1;
	} else if(grp_frnd_mat->num_grp_frnds==0) {
	    printf("num_grp_frnds=0, switch to default\n");
	    return -1;
	}

	//grp_frnd_mat->num_cmn_brs= (int*) malloc(grp_frnd_mat->num_grp_frnds* sizeof(int));
	MPIU_CHKPMEM_CALLOC(grp_frnd_mat->num_cmn_brs, int*,  grp_frnd_mat->num_grp_frnds* sizeof(int), mpi_errno, "grp_frnd_mat->num_cmn_brs");
	MPIU_CHKPMEM_MALLOC(grp_frnd_mat->is_active_grpfrnd, int*, grp_frnd_mat->num_grp_frnds * sizeof(int), mpi_errno, "grp_frnd_mat->is_active_grpfrnd");

	for(i=0;i<grp_frnd_mat->num_grp_frnds;i++) grp_frnd_mat->is_active_grpfrnd[i]=1;

	for(i=0;i<outdegree;i++) {
	    grp_frnd_mat->grpfrnd_bit_map[i] = MPIU_Malloc(grp_frnd_mat->num_grp_frnds * sizeof(int));
	}

	for(i=0;i<K_vrbl;i++) grp_frnd_mat->grpfriends[i] = MPIU_Malloc(grp_frnd_mat->num_grp_frnds * sizeof(int));

	for(i=0;i<outdegree;i++) {
	    for(j=0;j<grp_frnd_mat->num_grp_frnds;j++) grp_frnd_mat->grpfrnd_bit_map[i][j]=0;
	}

	for(j=0;j<grp_frnd_mat->num_grp_frnds;j++) grp_frnd_mat->num_cmn_brs[j]=0;

	int q=0;
	for(i=0;i<ind_frnd_mat->numfriends;i++) {
	    cnt_cmn_nbrs=0;
	    for(p=0;p<outdegree;p++) {
	        if(ind_frnd_mat->frndship_bit_map[i][p]) cnt_cmn_nbrs++;
	    }
	    if(cnt_cmn_nbrs>=nbr_frndshp_thr) {
	        for(p=0;p<outdegree;p++) {
		    grp_frnd_mat->grpfrnd_bit_map[p][q]=ind_frnd_mat->frndship_bit_map[i][p];

		    if(grp_frnd_mat->grpfrnd_bit_map[p][q]) grp_frnd_mat->num_actv_grp_of_nbr[p]++;

		}
		grp_frnd_mat->grpfriends[0][q]=ind_frnd_mat->friends[i];
		q++;
	    }
	}
	} else if(K_vrbl==3) {
	    for(i=0;i<ind_frnd_mat->numfriends;i++) {
	        for(j=i+1;j<ind_frnd_mat->numfriends;j++) {
		    cnt_cmn_nbrs=0;
		    for(p=0;p<outdegree;p++) {
		        if(ind_frnd_mat->frndship_bit_map[i][p]&ind_frnd_mat->frndship_bit_map[j][p]) {
			    cnt_cmn_nbrs++;
			}
		    }

		    if(cnt_cmn_nbrs>=nbr_frndshp_thr) {
		        grp_frnd_mat->num_grp_frnds++;
		    }
		}
	    }

	    if(grp_frnd_mat->num_grp_frnds>200000) {
	        printf("Not Enough Memory. Increase the Threshold\n");
		return -1;
	    } else if(grp_frnd_mat->num_grp_frnds==0) {
		printf("num_grp_frnds=0, switch to default\n");
		return -1;
	    }

	    MPIU_CHKPMEM_CALLOC(grp_frnd_mat->num_cmn_brs, int*,  grp_frnd_mat->num_grp_frnds* sizeof(int), mpi_errno, "grp_frnd_mat->num_cmn_brs");
	    MPIU_CHKPMEM_MALLOC(grp_frnd_mat->is_active_grpfrnd, int*, grp_frnd_mat->num_grp_frnds * sizeof(int), mpi_errno, "grp_frnd_mat->is_active_grpfrnd");

	    for(i=0;i<grp_frnd_mat->num_grp_frnds;i++) { 
	        grp_frnd_mat->is_active_grpfrnd[i]=1;
	    }

	    for(i=0;i<outdegree;i++) {
	        grp_frnd_mat->grpfrnd_bit_map[i] = MPIU_Calloc(grp_frnd_mat->num_grp_frnds, sizeof(int));
	    }

	    for(i=0;i<K_vrbl;i++) {
	    grp_frnd_mat->grpfriends[i] = MPIU_Malloc(grp_frnd_mat->num_grp_frnds * sizeof(int));
	    }

	    for(j=0;j<grp_frnd_mat->num_grp_frnds;j++) {
	        for(i=0;i<outdegree;i++) {
		    grp_frnd_mat->grpfrnd_bit_map[i][j]=0;
		}
	    }
	    for(j=0;j<grp_frnd_mat->num_grp_frnds;j++) {
	        grp_frnd_mat->num_cmn_brs[j]=0;
	    }

	    int q=0;
	    for(i=0;i<ind_frnd_mat->numfriends;i++) {
                for(j=i+1;j<ind_frnd_mat->numfriends;j++) {
		    cnt_cmn_nbrs=0;
		    for(p=0;p<outdegree;p++) {
		        if(ind_frnd_mat->frndship_bit_map[i][p]&ind_frnd_mat->frndship_bit_map[j][p]) {
						cnt_cmn_nbrs++;
			}
		    }
		    if(cnt_cmn_nbrs>=nbr_frndshp_thr) {
		        for(p=0;p<outdegree;p++) {
			    grp_frnd_mat->grpfrnd_bit_map[p][q]=ind_frnd_mat->frndship_bit_map[i][p]&ind_frnd_mat->frndship_bit_map[j][p];
			    if(grp_frnd_mat->grpfrnd_bit_map[p][q]) {
			        grp_frnd_mat->num_actv_grp_of_nbr[p]++;
			    }
			}
			grp_frnd_mat->grpfriends[0][q]=ind_frnd_mat->friends[i];
			grp_frnd_mat->grpfriends[1][q]=ind_frnd_mat->friends[j];
			q++;
		    }
		}
	    }
	} else if(K_vrbl==4) {
	    for(i=0;i<ind_frnd_mat->numfriends;i++) {
	        for(j=i+1;j<ind_frnd_mat->numfriends;j++) {
		    for(m=j+1;m<ind_frnd_mat->numfriends;m++) {
		        cnt_cmn_nbrs=0;
			for(p=0;p<outdegree;p++) {
			    if(ind_frnd_mat->frndship_bit_map[i][p]&ind_frnd_mat->frndship_bit_map[j][p]&ind_frnd_mat->frndship_bit_map[m][p]) {
			        cnt_cmn_nbrs++;
			    }
			}
			if(cnt_cmn_nbrs>=nbr_frndshp_thr) {
			    grp_frnd_mat->num_grp_frnds++;
			}
		    }
		}
	    }

	    if(grp_frnd_mat->num_grp_frnds>200000) {
	        printf("Not Enough Memory. Increase the Threshold\n");
		return -1;
	    } else if(grp_frnd_mat->num_grp_frnds==0) {
		printf("num_grp_frnds=0, switch to default\n");
		return -1;
	    }

	    MPIU_CHKPMEM_CALLOC(grp_frnd_mat->num_cmn_brs, int*,  grp_frnd_mat->num_grp_frnds* sizeof(int), mpi_errno, "grp_frnd_mat->num_cmn_brs");
	    MPIU_CHKPMEM_MALLOC(grp_frnd_mat->is_active_grpfrnd, int*, grp_frnd_mat->num_grp_frnds * sizeof(int), mpi_errno, "grp_frnd_mat->is_active_grpfrnd");

	    for(i=0;i<grp_frnd_mat->num_grp_frnds;i++) {
	        grp_frnd_mat->is_active_grpfrnd[i]=1;
	    }

	    for(i=0;i<outdegree;i++){
	        grp_frnd_mat->grpfrnd_bit_map[i] = MPIU_Calloc(grp_frnd_mat->num_grp_frnds, sizeof(int));
	    }

	    for(i=0;i<K_vrbl;i++) {
	        grp_frnd_mat->grpfriends[i] = MPIU_Malloc(grp_frnd_mat->num_grp_frnds * sizeof(int));
	    }

	    for(j=0;j<grp_frnd_mat->num_grp_frnds;j++) {
	        for(i=0;i<outdegree;i++) {
		    grp_frnd_mat->grpfrnd_bit_map[i][j]=0;
		}
	    }
	    for(j=0;j<grp_frnd_mat->num_grp_frnds;j++) {
	        grp_frnd_mat->num_cmn_brs[j]=0;
	    }

	    int q=0;
	    for(i=0;i<ind_frnd_mat->numfriends;i++) {
		for(j=i+1;j<ind_frnd_mat->numfriends;j++) {
		    for(m=j+1;m<ind_frnd_mat->numfriends;m++) {
		        cnt_cmn_nbrs=0;
			for(p=0;p<outdegree;p++) {
			    if(ind_frnd_mat->frndship_bit_map[i][p]&ind_frnd_mat->frndship_bit_map[j][p]&ind_frnd_mat->frndship_bit_map[m][p]) {
			        cnt_cmn_nbrs++;
			    }
			}
			if(cnt_cmn_nbrs>=nbr_frndshp_thr) {
			    for(p=0;p<outdegree;p++) {
			        grp_frnd_mat->grpfrnd_bit_map[p][q]=ind_frnd_mat->frndship_bit_map[i][p]&ind_frnd_mat->frndship_bit_map[j][p]&ind_frnd_mat->frndship_bit_map[m][p];
				if(grp_frnd_mat->grpfrnd_bit_map[p][q]) {
				    grp_frnd_mat->num_actv_grp_of_nbr[p]++;
				}
			    }
			    grp_frnd_mat->grpfriends[0][q]=ind_frnd_mat->friends[i];
			    grp_frnd_mat->grpfriends[1][q]=ind_frnd_mat->friends[j];
			    grp_frnd_mat->grpfriends[2][q]=ind_frnd_mat->friends[m];
			    q++;
			}
		    }
		}
	    }
	} 

	for(j=0;j<grp_frnd_mat->num_grp_frnds;j++) {
	    for(i=0;i<outdegree;i++)
	    grp_frnd_mat->num_cmn_brs[j]=grp_frnd_mat->num_cmn_brs[j]+grp_frnd_mat->grpfrnd_bit_map[i][j];
	}

	grp_frnd_mat->Total_num_grp_frnds= grp_frnd_mat->num_grp_frnds;
	*grp_frnd_mat_ptr=grp_frnd_mat;

	fn_fail:
	return 0;
}


//SMGM
#undef FUNCNAME
#define FUNCNAME MPIR_Make_Friendship_Matrix
#undef FCNAME
#define FCNAME MPIU_QUOTE(FUNCNAME)
int MPIR_Make_Friendship_Matrix(MPID_Comm *comm_ptr, Common_nbrhood_matrix *cmn_nbh_mat, Individual_Friendship_Matrix **ind_frnd_mat_ptr) {
	int mpi_errno = MPI_SUCCESS;
	int i, j,comm_size, self_rank, indegree,outdegree;
	comm_size = comm_ptr->local_size;
	self_rank = comm_ptr->rank;

	MPIR_Topology *topo_ptr = NULL;
	topo_ptr = MPIR_Topology_get(comm_ptr);
	if(topo_ptr == NULL) {
	    fprintf(stderr, "ERROR: Communicator topology pointer is NULL!\n");
	    goto fn_fail;
	}
	indegree = topo_ptr->topo.dist_graph.indegree;
	outdegree = topo_ptr->topo.dist_graph.outdegree;

	int num_frnds=0;

	int *glob_frndshp_arr;
	glob_frndshp_arr = MPIU_Malloc(comm_size * sizeof(int));
	for(i = 0; i < comm_size; i++) {
	    glob_frndshp_arr[i] = 0;
	}

	int mygrp_lwr_rng, mygrp_upr_rng;

	if(topo_aware==1) {
	    //determine the processes on the same node
	    MPI_Comm shmcomm;
	    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0,
		                    MPI_INFO_NULL, &shmcomm);
	    int shmrank, shmsz;
	    MPI_Comm_rank(shmcomm, &shmrank);
	    MPI_Comm_size(shmcomm, &shmsz);

	    int is_rank0, nodes;
	    is_rank0 = (shmrank == 0) ? 1 : 0;
	    MPI_Allreduce(&is_rank0, &nodes, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	    int proc_per_node=comm_size/nodes;
	    int num_large_grp=comm_size%nodes;
	    int num_small_grp=nodes-num_large_grp;

	    int grpnum;
	    if(self_rank< (num_small_grp*proc_per_node)) {
	        grpnum= self_rank/proc_per_node;
	        mygrp_lwr_rng=(grpnum*proc_per_node);
	        mygrp_upr_rng=(grpnum*proc_per_node)+proc_per_node-1;
	    } else {
	        int tmp= (self_rank - (proc_per_node*num_small_grp))/(proc_per_node+1);
	        grpnum=num_small_grp+tmp;
	        mygrp_lwr_rng= (num_small_grp*proc_per_node)+ (grpnum-num_small_grp)* (proc_per_node+1);
	        mygrp_upr_rng= mygrp_lwr_rng+proc_per_node;
	    }
	} else {
	    mygrp_lwr_rng=0;
	    mygrp_upr_rng=comm_size-1;
	}

	for(i = 0; i < cmn_nbh_mat->num_rows; i++) {
	    for(j = 0; j < cmn_nbh_mat->row_sizes[i]; j++) {
	        if(cmn_nbh_mat->matrix[i][j] != self_rank && cmn_nbh_mat->matrix[i][j]>=mygrp_lwr_rng && cmn_nbh_mat->matrix[i][j]<=mygrp_upr_rng) {
		    glob_frndshp_arr[cmn_nbh_mat->matrix[i][j]]++;
		    if(glob_frndshp_arr[cmn_nbh_mat->matrix[i][j]] == nbr_frndshp_thr) {
		        num_frnds++;
		    }
		}
	    }
	}

	MPIU_CHKPMEM_DECL(9);

	Individual_Friendship_Matrix *ind_frnd_mat;
	MPIU_CHKPMEM_MALLOC(ind_frnd_mat, Individual_Friendship_Matrix*, sizeof(Individual_Friendship_Matrix), mpi_errno, "ind_frnd_mat");
	MPIU_CHKPMEM_CALLOC(ind_frnd_mat->friends, int*,  num_frnds* sizeof(int), mpi_errno, "ind_frnd_mat->friends");
	MPIU_CHKPMEM_CALLOC(ind_frnd_mat->nbrs, int*,  outdegree* sizeof(int), mpi_errno, "ind_frnd_mat->nbrs");
	MPIU_CHKPMEM_MALLOC(ind_frnd_mat->frndship_bit_map, int**, num_frnds * sizeof(int*), mpi_errno, "ind_frnd_mat->frndship_bit_map");

	ind_frnd_mat->numfriends=num_frnds;
	for(i=0;i<num_frnds;i++)
	ind_frnd_mat->frndship_bit_map[i] = MPIU_Malloc(outdegree * sizeof(int));

	j=0;
	for(i=0;i<comm_size;i++) {
	    if(glob_frndshp_arr[i]>=nbr_frndshp_thr && i!=self_rank && i>=mygrp_lwr_rng && i<=mygrp_upr_rng) {
	        ind_frnd_mat->friends[j]=i;
		j++;
	    }
	}

	ind_frnd_mat->nbrs = topo_ptr->topo.dist_graph.out;

	int p;
	for(p=0;p<num_frnds;p++) {
	    for(i=0;i<outdegree; i++) {
	        ind_frnd_mat->frndship_bit_map[p][i]=0;
	    }
	}

	for(p=0;p<num_frnds;p++) {
	    for(i = 0; i < outdegree; i++) {
	        for(j = 0; j < cmn_nbh_mat->row_sizes[i]; j++) {
		    if(cmn_nbh_mat->matrix[i][j]==ind_frnd_mat->friends[p]) {
		        ind_frnd_mat->frndship_bit_map[p][i]=1;
 		    }
		}
	    }
	}

	*ind_frnd_mat_ptr= ind_frnd_mat;


	fn_fail:
	return 0;
}

#undef FUNCNAME
#define FUNCNAME MPIR_Get_inNbrs_of_outNbrs
#undef FCNAME
#define FCNAME MPIU_QUOTE(FUNCNAME)
int MPIR_Get_inNbrs_of_outNbrs(MPID_Comm *comm_ptr, Common_nbrhood_matrix **cmn_nbh_mat_ptr) {
    //This function returns the matrix that gives the
    //incoming neighbors of each of my outgoing neighbors.
    //Number of rows will be equal to my outdegree, and
    //the number of elements in row i will be equal to
    //the indegree of my ith outgoing neighbor.

    int mpi_errno = MPI_SUCCESS;
    int indegree, outdegree, comm_size, i, j, out_idx, in_idx, all_reqs_idx, reqs_max_size, context_offset;
    comm_size = comm_ptr->local_size;
    MPIR_Topology *topo_ptr = NULL;
    topo_ptr = MPIR_Topology_get(comm_ptr);
    if(topo_ptr == NULL) {
        fprintf(stderr, "ERROR: Communicator topology pointer is NULL!\n");
        goto fn_fail;
    }
    indegree = topo_ptr->topo.dist_graph.indegree;
    outdegree = topo_ptr->topo.dist_graph.outdegree;

    MPIU_CHKPMEM_DECL(9);

    /* Getting the indegree of my outgoing neighbors.
     * To this end, everybody should send its indegree
     * to each of its incoming neighbors, and also receive
     * the indegree from each of its outgoing neighbors.
    */
    context_offset = (comm_ptr->comm_kind == MPID_INTRACOMM) ?
                      MPID_CONTEXT_INTRA_COLL : MPID_CONTEXT_INTER_COLL;
    reqs_max_size = indegree + outdegree;
    all_reqs_idx = 0;
    MPI_Request *all_reqs = MPIU_Malloc(reqs_max_size * sizeof(MPI_Request));
    MPID_Request *req_ptr = NULL;

    //Recv buffer
    int *outnbrs_indegree;
    MPIU_CHKPMEM_MALLOC(outnbrs_indegree, int*, outdegree * sizeof(int), mpi_errno, "outnbrs_indegree");

    //Send the indegree to each incoming neighbor
    for(in_idx = 0; in_idx < indegree; in_idx++) //for each of my incoming neighbors {
        mpi_errno = MPID_Isend(&indegree, 1, MPI_INT, topo_ptr->topo.dist_graph.in[in_idx],
		                       1000, comm_ptr, context_offset, &req_ptr);
	if (mpi_errno) MPIU_ERR_POP(mpi_errno);
	all_reqs[all_reqs_idx++] = req_ptr->handle;
    }
    for(out_idx = 0; out_idx < outdegree; out_idx++) //for each of my outgoing neighbors {
	mpi_errno = MPID_Irecv(&outnbrs_indegree[out_idx], 1, MPI_INT,
		                       topo_ptr->topo.dist_graph.out[out_idx],
		                       1000, comm_ptr, context_offset, &req_ptr);
	if (mpi_errno) MPIU_ERR_POP(mpi_errno);
	all_reqs[all_reqs_idx++] = req_ptr->handle;
    }

    mpi_errno = MPIR_Waitall_impl(all_reqs_idx, all_reqs, MPI_STATUS_IGNORE);
    if (mpi_errno) MPIU_ERR_POP(mpi_errno);
    all_reqs_idx = 0; //set index back to zero for future use

#ifdef SHM_DEBUG
	print_vec(comm_ptr->rank, outdegree, outnbrs_indegree, "My outnbrs_indegree is");
#endif
	/** Done with getting the indegree of my outgoing neighbors **/

	/** Going to get the incoming neighbors of my outgoing neighbors **/
	/* I should send the list of my incoming neighbors to each of my
	 * incoming neighbors. Plus, I should receive a list (of incoming
	 * neighbors) from each of my outgoing neighbors.
	 */

	//Build recv buffers first which would be the result matrix
	Common_nbrhood_matrix *cmn_nbh_mat;
	MPIU_CHKPMEM_MALLOC(cmn_nbh_mat, Common_nbrhood_matrix*, sizeof(Common_nbrhood_matrix), mpi_errno, "cmn_nbh_mat");
	MPIU_CHKPMEM_CALLOC(cmn_nbh_mat->ignore_row, int*, outdegree * sizeof(int), mpi_errno, "cmn_nbh_mat->ignore_row");
	MPIU_CHKPMEM_CALLOC(cmn_nbh_mat->is_row_offloaded, int*, outdegree * sizeof(int), mpi_errno, "cmn_nbh_mat->is_row_offloaded");
	MPIU_CHKPMEM_MALLOC(cmn_nbh_mat->my_innbrs_bitmap, int*, indegree * sizeof(int), mpi_errno, "cmn_nbh_mat->my_innbrs_bitmap");
	MPIU_CHKPMEM_MALLOC(cmn_nbh_mat->outnbrs_innbrs_bitmap, int**, outdegree * sizeof(int*), mpi_errno, "cmn_nbh_mat->outnbrs_innbrs_bitmap");
	MPIU_CHKPMEM_MALLOC(cmn_nbh_mat->matrix, int**, outdegree * sizeof(int*), mpi_errno, "cmn_nbh_mat->matrix");
	MPIU_CHKPMEM_CALLOC(cmn_nbh_mat->comb_matrix, Comb_element**, outdegree* sizeof(Comb_element*), mpi_errno, "cmn_nbh_mat->comb_matrix");
	//MPIU_CHKPMEM_CALLOC(cmn_nbh_mat->num_entries_in_row, int*, outdegree* sizeof(int), mpi_errno, "cmn_nbh_mat->num_entries_in_row");

	cmn_nbh_mat->num_rows = outdegree;
        cmn_nbh_mat->indegree = indegree;
        cmn_nbh_mat->num_elements = 0; //figure out later
        cmn_nbh_mat->t = 0; // Keeping track of pairing steps order
        cmn_nbh_mat->row_sizes = outnbrs_indegree;

        for(i=0;i<outdegree;i++) {
    	cmn_nbh_mat->comb_matrix[i] = MPIU_Malloc(MAX_COMB_DEGREE * sizeof(Comb_element));
        }

	for(in_idx = 0; in_idx < indegree; in_idx++) //for each of my incoming neighbors {
	    mpi_errno = MPID_Isend(topo_ptr->topo.dist_graph.in, indegree, MPI_INT,
	                           topo_ptr->topo.dist_graph.in[in_idx],
	                           2000, comm_ptr, context_offset, &req_ptr);
            if (mpi_errno) MPIU_ERR_POP(mpi_errno);
            all_reqs[all_reqs_idx++] = req_ptr->handle;
		cmn_nbh_mat->my_innbrs_bitmap[in_idx] = 1; //Set all incoming neighbors to active
	}
	for(out_idx = 0; out_idx < outdegree; out_idx++) //for each of my outgoing neighbors {
	    cmn_nbh_mat->ignore_row[out_idx]=0;
	    cmn_nbh_mat->is_row_offloaded[out_idx]=0;
	    cmn_nbh_mat->outnbrs_innbrs_bitmap[out_idx] = MPIU_Malloc(outnbrs_indegree[out_idx] * sizeof(int));
	    for(j = 0; j < outnbrs_indegree[out_idx]; j++) {
	        cmn_nbh_mat->outnbrs_innbrs_bitmap[out_idx][j] = 1; //Set all incoming neighbors of all outgoing neighbors to active
	    }

	    cmn_nbh_mat->matrix[out_idx] = MPIU_Malloc(outnbrs_indegree[out_idx] * sizeof(int));
	    mpi_errno = MPID_Irecv(cmn_nbh_mat->matrix[out_idx], outnbrs_indegree[out_idx],
	                           MPI_INT, topo_ptr->topo.dist_graph.out[out_idx],
                                   2000, comm_ptr, context_offset, &req_ptr);
            if (mpi_errno) MPIU_ERR_POP(mpi_errno);
            all_reqs[all_reqs_idx++] = req_ptr->handle;
		cmn_nbh_mat->num_elements += cmn_nbh_mat->row_sizes[out_idx];
	}

	mpi_errno = MPIR_Waitall_impl(all_reqs_idx, all_reqs, MPI_STATUS_IGNORE);
	if (mpi_errno) MPIU_ERR_POP(mpi_errno);
	all_reqs_idx = 0; //set index back to zero for future use

	*cmn_nbh_mat_ptr = cmn_nbh_mat;

    fn_exit:
        MPIU_Free(all_reqs);

        return 0;
    fn_fail:
        MPIU_CHKPMEM_REAP();
        goto fn_exit;
}

int find_cmn_nbrs(MPID_Comm *comm_ptr, Group_Friendship_Matrix* grp_frnd_mat, Final_Group_Friends slcted_grp_frnd, int outdegree, int *dests, Common_neighbor** cmn_nbrs_ptr) {
	int i,j,k;
	int self_rank = comm_ptr->rank;


	int slct_grp_indx;

	for(i=0;i<grp_frnd_mat->Total_num_grp_frnds;i++) {
	    if(grp_frnd_mat->is_active_grpfrnd[i]) {
		k=0;
		for(j=0;j<K_vrbl;j++) {
		    if(slcted_grp_frnd.grp_frnds[j]!=self_rank) {
		        if(slcted_grp_frnd.grp_frnds[j]==grp_frnd_mat->grpfriends[k][i]) {
			    k++;
			} else {
			    j=K_vrbl; //else move to the next group friend in grp_frnd_mat
			}
		    }
		}
		if(k==K_vrbl-1) {
		    slct_grp_indx=i;
		    goto slct_grp_found;
		}
	    }
	}
	slct_grp_found:

	k=0;
	Common_neighbor *cmn_nbrs = MPIU_Malloc(grp_frnd_mat->num_cmn_brs[slct_grp_indx] * sizeof(Common_neighbor));
	for(i=0;i<outdegree;i++) {
	    if(grp_frnd_mat->grpfrnd_bit_map[i][slct_grp_indx]==1) {
	        cmn_nbrs[k].index=i;
		cmn_nbrs[k].rank= dests[i];
		k++;
	    }
	}

	*cmn_nbrs_ptr=cmn_nbrs;

	return slct_grp_indx;
}



int add_frnd_to_comb_matrix(Common_nbrhood_matrix *cmn_nbh_mat, int row_idx, int frnd, Operation opt, int grp_frnd_num) {
	 int t = cmn_nbh_mat->t;

	 if(t >= MAX_COMB_DEGREE) {
	     fprintf(stderr, "ERROR: No more space to add a new scheduling step to the comb_matrix! t = %d\n", t);
	     return -1;
	 }

	 cmn_nbh_mat->comb_matrix[row_idx][t].grp_frnds[grp_frnd_num] = frnd;
	 cmn_nbh_mat->comb_matrix[row_idx][t].opt[grp_frnd_num] = opt;

	 return 0;
}

#undef FUNCNAME
#define FUNCNAME MPIR_Update_grp_frnd_mat
#undef FCNAME
#define FCNAME MPIU_QUOTE(FUNCNAME)
int MPIR_Update_grp_frnd_mat(Common_nbrhood_matrix* cmn_nbh_mat, Group_Friendship_Matrix *grp_frnd_mat, MPID_Comm *comm_ptr) {
	int mpi_errno = MPI_SUCCESS;
	int indegree, outdegree, comm_size, i, j,k, out_idx, in_idx, all_reqs_idx, reqs_max_size, context_offset;
	int num_removed_nbrs_in_update=0;
	comm_size = comm_ptr->local_size;
	int self_rank = comm_ptr->rank;
	MPIR_Topology *topo_ptr = NULL;
	topo_ptr = MPIR_Topology_get(comm_ptr);
	if(topo_ptr == NULL) {
	    fprintf(stderr, "ERROR: Communicator topology pointer is NULL!\n");
	    return -1;
	}
	indegree = topo_ptr->topo.dist_graph.indegree;
	outdegree = topo_ptr->topo.dist_graph.outdegree;

	context_offset = (comm_ptr->comm_kind == MPID_INTRACOMM) ?
			MPID_CONTEXT_INTRA_COLL : MPID_CONTEXT_INTER_COLL;
	reqs_max_size = indegree + outdegree;
	all_reqs_idx = 0;
	MPI_Request *all_reqs = MPIU_Malloc(reqs_max_size * sizeof(MPI_Request));
	MPID_Request *req_ptr = NULL;

	MPIU_CHKPMEM_DECL(9);


	//Sending my own innbrs bitmap to each not-ignored incoming neighbor
	for(in_idx = 0; in_idx < indegree; in_idx++) //for each of my incoming neighbors {
		if(!cmn_nbh_mat->my_innbrs_bitmap[in_idx]) continue;
		mpi_errno = MPID_Isend(cmn_nbh_mat->my_innbrs_bitmap, indegree, MPI_INT,
				topo_ptr->topo.dist_graph.in[in_idx], 1100,
				comm_ptr, context_offset, &req_ptr);

		if (mpi_errno) MPIU_ERR_POP(mpi_errno);
		all_reqs[all_reqs_idx++] = req_ptr->handle;
	}

	int *recv_flag;
	MPIU_CHKPMEM_MALLOC(recv_flag, int*, outdegree*sizeof(int), mpi_errno, "recv_flag");
	for(i=0;i<outdegree;i++)
		recv_flag[i]=0;
	//Receiving the innbrs bitmap of each of my not-ignored outgoing neighbors
	for(out_idx = 0; out_idx < outdegree && grp_frnd_mat->num_grp_frnds != 0; out_idx++) { //for each of my outgoing neighbors
	    if(!grp_frnd_mat->is_active_nbr[out_idx]) continue;
	    mpi_errno = MPID_Irecv(cmn_nbh_mat->outnbrs_innbrs_bitmap[out_idx],
				cmn_nbh_mat->row_sizes[out_idx], MPI_INT,
				topo_ptr->topo.dist_graph.out[out_idx],
				1100, comm_ptr, context_offset, &req_ptr);
	    recv_flag[out_idx]=1;
	    if (mpi_errno) MPIU_ERR_POP(mpi_errno);
	    all_reqs[all_reqs_idx++] = req_ptr->handle;
	}

	mpi_errno = MPIR_Waitall_impl(all_reqs_idx, all_reqs, MPI_STATUS_IGNORE);
	if (mpi_errno) MPIU_ERR_POP(mpi_errno);
	all_reqs_idx = 0; //set index back to zero for future use
	/** Done with getting the innbrs bitmaps of my outgoing neighbors **/

	//Update the grp_frnd_bit_map based on the updated outnbrs_innbrs_bitmap matrix
	int num_nonactive_grps=0;
	i=0;
	for(out_idx = 0; out_idx < outdegree ; out_idx++) {
	    num_nonactive_grps=0;
	    for(in_idx = 0; in_idx < cmn_nbh_mat->row_sizes[out_idx]; in_idx++) {
	        if(cmn_nbh_mat->outnbrs_innbrs_bitmap[out_idx][in_idx]==0) {
		    if(grp_frnd_mat->is_active_nbr[out_idx]) {
		        for(i=0;i<grp_frnd_mat->Total_num_grp_frnds;i++) {
			    if(grp_frnd_mat->is_active_grpfrnd[i]) {
			        if(grp_frnd_mat->grpfrnd_bit_map[out_idx][i]) {
				    for(k=0;k<K_vrbl-1;k++) {
				        if(grp_frnd_mat->grpfriends[k][i]==cmn_nbh_mat->matrix[out_idx][in_idx]) {
					    grp_frnd_mat->grpfrnd_bit_map[out_idx][i]=0;
						grp_frnd_mat->num_actv_grp_of_nbr[out_idx]--;
						if(grp_frnd_mat->num_actv_grp_of_nbr[out_idx]==0) {
						    grp_frnd_mat->is_active_nbr[out_idx]=2;
						    num_removed_nbrs_in_update++;
						}

						grp_frnd_mat->num_cmn_brs[i]--;

						if(grp_frnd_mat->num_cmn_brs[i]<nbr_frndshp_thr && grp_frnd_mat->is_active_grpfrnd[i]) { 
						    //we add the second condition to make sure this grp_frnd has not been removed previously
						    grp_frnd_mat->is_active_grpfrnd[i]=0;
						    grp_frnd_mat->num_grp_frnds--;
						    num_nonactive_grps++;
						}
					    }
					}
				     } else {
				         num_nonactive_grps++;
				     }
				} else {
				    num_nonactive_grps++;
				}
			    }
			}
		    }
		}
	    }


	fn_fail:
	return num_removed_nbrs_in_update;
}

#undef FUNCNAME
#define FUNCNAME MPIR_Build_SHM_nbh_coll_patt
#undef FCNAME
#define FCNAME MPIU_QUOTE(FUNCNAME)
int MPIR_Build_SHM_nbh_coll_patt(MPID_Comm *comm_ptr) {
    int mpi_errno = MPI_SUCCESS;
    int i, j,k,p, all_reqs_idx, reqs_max_size, context_offset;
    int self_rank = comm_ptr->rank;
    int comm_size = comm_ptr->local_size;
    int steps=0;

    int ON = 1;
    int OFF = 0;

    //Getting the topology pointer
    MPIR_Topology *topo_ptr = NULL;
    topo_ptr = MPIR_Topology_get(comm_ptr);
    if(topo_ptr == NULL) {
        fprintf(stderr, "ERROR: Communicator topology pointer is NULL!\n");
        goto fn_fail;
    }
    int indegree = topo_ptr->topo.dist_graph.indegree;
    int outdegree = topo_ptr->topo.dist_graph.outdegree;
    int *srcs = topo_ptr->topo.dist_graph.in;
    int *dests = topo_ptr->topo.dist_graph.out;

    reqs_max_size = indegree + outdegree;
    all_reqs_idx = 0;
    MPI_Request *all_reqs = MPIU_Malloc(reqs_max_size * sizeof(MPI_Request));
    MPID_Request *req_ptr = NULL;

    MPIU_CHKPMEM_DECL(9);

    context_offset = (comm_ptr->comm_kind == MPID_INTRACOMM) ?
                      MPID_CONTEXT_INTRA_COLL : MPID_CONTEXT_INTER_COLL;

    //Extract common neighborhood matrix
    Common_nbrhood_matrix *cmn_nbh_mat;
    mpi_errno = MPIR_Get_inNbrs_of_outNbrs(comm_ptr, &cmn_nbh_mat);

    for(i=0;i<outdegree;i++) {
        for(j = 0; j < MAX_COMB_DEGREE; j++) {
            MPIU_CHKPMEM_CALLOC(cmn_nbh_mat->comb_matrix[i][j].grp_frnds, int*, K_vrbl* sizeof(int), mpi_errno, "cmn_nbh_mat->comb_matrix[i][j].grp_frnds");
            MPIU_CHKPMEM_CALLOC(cmn_nbh_mat->comb_matrix[i][j].opt, Operation*, K_vrbl* sizeof(Operation), mpi_errno, "cmn_nbh_mat->comb_matrix[i][j].opt");
        }
    }

    for(i=0;i<outdegree;i++) {
        for(j = 0; j < MAX_COMB_DEGREE; j++) {
     	    for(k=0;k<K_vrbl;k++) {
     	        cmn_nbh_mat->comb_matrix[i][j].opt[k] = IDLE;
     	    	cmn_nbh_mat->comb_matrix[i][j].grp_frnds[k] = -1;
     	    }
     	}
    }


    if(mpi_errno) MPIU_ERR_POP(mpi_errno);

    Individual_Friendship_Matrix *ind_frnd_mat;
    MPIR_Make_Friendship_Matrix(comm_ptr, cmn_nbh_mat, &ind_frnd_mat);


    Group_Friendship_Matrix *grp_frnd_mat;
    int hyperedge= MPIR_Make_GroupFriend_Matrix(comm_ptr,ind_frnd_mat, &grp_frnd_mat);
    if(hyperedge==-1)
    	return -1;

    Final_Group_Friends slcted_grp_frnd;

    do {

    MPIR_SMGM_Find_Group_Friend(comm_ptr, grp_frnd_mat, &slcted_grp_frnd, steps);

    if(slcted_grp_frnd.is_grp==1 && grp_frnd_mat->num_grp_frnds!=0) {
    	Common_neighbor *cmn_nbrs;
    	int slct_grp_indx= find_cmn_nbrs(comm_ptr, grp_frnd_mat, slcted_grp_frnd, outdegree, dests, &cmn_nbrs);

    	int num_cmn_nbrs= grp_frnd_mat->num_cmn_brs[slct_grp_indx];

	//mask the group friends so as NOT to pair with it again
	// Do this by in-activing the corresponding column in grp_frnd_mat
    	if(grp_frnd_mat->is_active_grpfrnd[slct_grp_indx]==0) {
    	    printf("rank%d: ERROR: a group friend is selected twice\n", self_rank);
    	} else {
    	    grp_frnd_mat->is_active_grpfrnd[slct_grp_indx]=0;
        }
    	grp_frnd_mat->num_grp_frnds--;

    	//find the index of self_rank in the group friend
    	int my_indx;
    	for(i=0;i<K_vrbl;i++) {
    	    if(slcted_grp_frnd.grp_frnds[i]==self_rank) {
    	        my_indx=i;
	    }
	}

    	int aux_indx=0;
    	int *gap;

    	gap=MPIU_Malloc(K_vrbl * sizeof(int));

    	for(i=0;i<K_vrbl;i++) {
    	    gap[i]=num_cmn_nbrs/K_vrbl;
	}

    	int reminder= num_cmn_nbrs%K_vrbl;

    	for(i=0;i<reminder;i++) {
    	    gap[i]=gap[i]+1;
	}


    	int i_grpfrnd=0;
    	for(i=0;i<K_vrbl;i++) {
    		if(my_indx==i) { //keep the common neighbors
    		    int start_on_indx= aux_indx;
    		    int end_on_indx= aux_indx+gap[i];
    		    aux_indx=end_on_indx;

    		    for(j=start_on_indx;j<end_on_indx;j++) {
    		        grp_frnd_mat->is_active_nbr[cmn_nbrs[j].index]=0;
    			cmn_nbh_mat->ignore_row[cmn_nbrs[j].index]=1;

    			//update grp_frnd_mat considering the removed common neighbor. Update the group friends that have the same common neighbor
    			for(p=0;p<grp_frnd_mat->Total_num_grp_frnds;p++) {
    			    if(grp_frnd_mat->is_active_grpfrnd[p]) {
    			        if(grp_frnd_mat->grpfrnd_bit_map[cmn_nbrs[j].index][p]==1) {
    				    grp_frnd_mat->grpfrnd_bit_map[cmn_nbrs[j].index][p]=0;
    				    grp_frnd_mat->num_actv_grp_of_nbr[cmn_nbrs[j].index]--;

    				    grp_frnd_mat->num_cmn_brs[p]--;
    				    if(grp_frnd_mat->num_cmn_brs[p]<nbr_frndshp_thr) {
    				        grp_frnd_mat->is_active_grpfrnd[p]=0;
    					grp_frnd_mat->num_grp_frnds--;
    				    }
    			        }
    			    }
    			}

    			mpi_errno = MPID_Isend(&OFF, 1, MPI_INT, dests[cmn_nbrs[j].index], 1000, comm_ptr, context_offset, &req_ptr);
    			if (mpi_errno) MPIU_ERR_POP(mpi_errno);
    			all_reqs[all_reqs_idx++] = req_ptr->handle;

    			int k_grpfrnd=0;
    			for(k=0;k<K_vrbl;k++) {
    			    if(slcted_grp_frnd.grp_frnds[k]!=self_rank) {
    			        add_frnd_to_comb_matrix(cmn_nbh_mat, cmn_nbrs[j].index, slcted_grp_frnd.grp_frnds[k], RECV, k_grpfrnd);
    				k_grpfrnd++;
    			    }
    			}
    		    }
    		} else {  //offload common neighbors to other friends
    		    int start_off_indx= aux_indx;
    		    int end_off_indx= aux_indx+gap[i];
    		    aux_indx=end_off_indx;

    		    for(j=start_off_indx;j<end_off_indx;j++) {
    		        grp_frnd_mat->is_active_nbr[cmn_nbrs[j].index]=0;
    			cmn_nbh_mat->is_row_offloaded[cmn_nbrs[j].index]=1;
    			cmn_nbh_mat->ignore_row[cmn_nbrs[j].index]=1;

    			//update grp_frnd_mat considering the removed common neighbor. Update the group friends that have the same common neighbor
    			for(p=0;p<grp_frnd_mat->Total_num_grp_frnds;p++) {
    			    if(grp_frnd_mat->is_active_grpfrnd[p]) {
    			        if(grp_frnd_mat->grpfrnd_bit_map[cmn_nbrs[j].index][p]==1) {
    			            grp_frnd_mat->grpfrnd_bit_map[cmn_nbrs[j].index][p]=0;
    				    grp_frnd_mat->num_actv_grp_of_nbr[cmn_nbrs[j].index]--;

    				    grp_frnd_mat->num_cmn_brs[p]--;
    				    if(grp_frnd_mat->num_cmn_brs[p]<nbr_frndshp_thr) {
    				        grp_frnd_mat->is_active_grpfrnd[p]=0;
    				        grp_frnd_mat->num_grp_frnds--;
    				    }
    			        }
    			    }
    			}

    			mpi_errno = MPID_Isend(&OFF, 1, MPI_INT, dests[cmn_nbrs[j].index], 1000, comm_ptr, context_offset, &req_ptr);
    			if (mpi_errno) MPIU_ERR_POP(mpi_errno);
    			all_reqs[all_reqs_idx++] = req_ptr->handle;

    			//add_frnd_to_comb_matrix(cmn_nbh_mat, cmn_nbrs[j].index, slcted_grp_frnd.grp_frnds[i], SEND, i_grpfrnd);
    			int k_grpfrnd=0;
    			for(k=0;k<K_vrbl;k++) {
    			    if(slcted_grp_frnd.grp_frnds[k]!=self_rank) {
    			        if(k==i) {
    				    add_frnd_to_comb_matrix(cmn_nbh_mat, cmn_nbrs[j].index, slcted_grp_frnd.grp_frnds[k], SEND, k_grpfrnd);
				} else {
    				    add_frnd_to_comb_matrix(cmn_nbh_mat, cmn_nbrs[j].index, slcted_grp_frnd.grp_frnds[k], HALT, k_grpfrnd);
    				    k_grpfrnd++;
				}
    			    }
    			}
    		    }
    		    i_grpfrnd++;
    		}
    	    }
        }


    for(i=0;i<outdegree;i++) {
        if(grp_frnd_mat->is_active_nbr[i]==1) {
    	    if(grp_frnd_mat->num_grp_frnds==0) {
    	        /* We send OFF because this rank will be quitting
    		 * the main while loop with num_frnd == 0, and so
    		 * others should not be expecting to receive any
    		 * more notifications from it.
    		 */
    		 mpi_errno = MPID_Isend(&OFF, 1, MPI_INT, dests[i], 1000,
    					comm_ptr, context_offset, &req_ptr);
    		 grp_frnd_mat->is_active_nbr[i]=0;
    		 if (mpi_errno) MPIU_ERR_POP(mpi_errno);
    		 all_reqs[all_reqs_idx++] = req_ptr->handle;
    	    } else {
    	         mpi_errno = MPID_Isend(&ON, 1, MPI_INT, dests[i], 1000,
    					comm_ptr, context_offset, &req_ptr);
    		 if (mpi_errno) MPIU_ERR_POP(mpi_errno);
    		 all_reqs[all_reqs_idx++] = req_ptr->handle;
    	    }
    	}
    }

    int *recv_flag;
    MPIU_CHKPMEM_MALLOC(recv_flag, int*, indegree*sizeof(int), mpi_errno, "recv_flag");
    for(i=0;i<indegree;i++) {
    	recv_flag[i]=0;
    }

    for(i = 0; i < indegree; i++) { //for each of my still-active incoming neighbors
    	if(cmn_nbh_mat->my_innbrs_bitmap[i] == 1) {
    	    mpi_errno = MPID_Irecv(&(cmn_nbh_mat->my_innbrs_bitmap[i]), 1, MPI_INT,
    				srcs[i], 1000, comm_ptr, context_offset, &req_ptr);
    	    recv_flag[i]=1;
    	    if (mpi_errno) MPIU_ERR_POP(mpi_errno);
    	    all_reqs[all_reqs_idx++] = req_ptr->handle;
    	}
    }
    mpi_errno = MPIR_Waitall_impl(all_reqs_idx, all_reqs, MPI_STATUS_IGNORE);
    if (mpi_errno) MPIU_ERR_POP(mpi_errno);
    all_reqs_idx = 0; //set index back to zero for future use


    for(i=0;i<indegree;i++) {
        recv_flag[i]=0; //for future use
    }

    int num_removed_nbrs_in_update=0, *rcv_num_removed_nbrs_in_update;
    MPIU_CHKPMEM_CALLOC(rcv_num_removed_nbrs_in_update, int*, indegree* sizeof(int), mpi_errno, "rcv_num_removed_nbrs_in_update");

    for(i=0;i<indegree;i++) rcv_num_removed_nbrs_in_update[i]=0;

    for(i=0;i<indegree;i++) rcv_num_removed_nbrs_in_update[i]=0;

    num_removed_nbrs_in_update= MPIR_Update_grp_frnd_mat(cmn_nbh_mat, grp_frnd_mat, comm_ptr);

    //perform Alltoall to see if a new neighbor is removed from any of the processes

    for(i=0;i<outdegree;i++) {
        if(grp_frnd_mat->is_active_nbr[i]) {
    	    mpi_errno = MPID_Isend(&num_removed_nbrs_in_update, 1, MPI_INT, dests[i], 88,
    				comm_ptr, context_offset, &req_ptr);
    	    if (mpi_errno) MPIU_ERR_POP(mpi_errno);
    	    all_reqs[all_reqs_idx++] = req_ptr->handle;
    	}
    }

    for(i = 0; i < indegree; i++) //for each of my still-active incoming neighbors {
    	if(cmn_nbh_mat->my_innbrs_bitmap[i] == 1) {
    	    mpi_errno = MPID_Irecv(&rcv_num_removed_nbrs_in_update[i], 1, MPI_INT,
    				srcs[i], 88, comm_ptr, context_offset, &req_ptr);
    	    recv_flag[i]=1;
    	    if (mpi_errno) MPIU_ERR_POP(mpi_errno);
    	    all_reqs[all_reqs_idx++] = req_ptr->handle;
    	}
    }
    mpi_errno = MPIR_Waitall_impl(all_reqs_idx, all_reqs, MPI_STATUS_IGNORE);
    if (mpi_errno) MPIU_ERR_POP(mpi_errno);
    all_reqs_idx = 0;


    for(i=0;i<indegree;i++) recv_flag[i]=0; //for future use


    for(i=0;i<outdegree;i++) {
    	if(grp_frnd_mat->is_active_nbr[i]==2) {
    	    mpi_errno = MPID_Isend(&OFF, 1, MPI_INT, dests[i], 44,
    				comm_ptr, context_offset, &req_ptr);

    	    grp_frnd_mat->is_active_nbr[i]=0;
    	    if (mpi_errno) MPIU_ERR_POP(mpi_errno);
    	    all_reqs[all_reqs_idx++] = req_ptr->handle;
    	} else if(grp_frnd_mat->is_active_nbr[i]==1 && num_removed_nbrs_in_update!=0) {
    	    mpi_errno = MPID_Isend(&ON, 1, MPI_INT, dests[i], 44,
    				comm_ptr, context_offset, &req_ptr);

    	    if (mpi_errno) MPIU_ERR_POP(mpi_errno);
    	    all_reqs[all_reqs_idx++] = req_ptr->handle;
    	}
     }

    for(i=0;i<indegree;i++) {
        if(rcv_num_removed_nbrs_in_update[i]) {
    		mpi_errno = MPID_Irecv(&(cmn_nbh_mat->my_innbrs_bitmap[i]), 1, MPI_INT,
    				srcs[i], 44, comm_ptr, context_offset, &req_ptr);
    		recv_flag[i]=1;
    		if (mpi_errno) MPIU_ERR_POP(mpi_errno);
    		all_reqs[all_reqs_idx++] = req_ptr->handle;
    	}
    }

    mpi_errno = MPIR_Waitall_impl(all_reqs_idx, all_reqs, MPI_STATUS_IGNORE);
    if (mpi_errno) MPIU_ERR_POP(mpi_errno);
    all_reqs_idx = 0;

    if(grp_frnd_mat->num_grp_frnds==0) {
    	for(i=0;i<outdegree;i++) {
    	    if(grp_frnd_mat->is_active_nbr[i]) {
    	        /* We send OFF because this rank will be quitting
    	       	 * the main while loop with num_frnd == 0, and so
    	       	 * others should not be expecting to receive any
    	       	 * more notifications from it.
    	       	 */
    	       	 mpi_errno = MPID_Isend(&OFF, 1, MPI_INT, dests[i], 1000,
    	       			comm_ptr, context_offset, &req_ptr);

    	       	 grp_frnd_mat->is_active_nbr[i]=0;
    	       	 if (mpi_errno) MPIU_ERR_POP(mpi_errno);
    	       	 all_reqs[all_reqs_idx++] = req_ptr->handle;
    	    }
    	}
    }

    steps++;
    cmn_nbh_mat->t++;

    } while (grp_frnd_mat->num_grp_frnds>0);

    /* Once a rank gets out of the while loop above due
       * to num_frnds == 0, it should still issue the recv
       * operations corresponding to its still-active incoming
       * neighbors because those neighbors will be sending
       * notifications to this rank until they get out of the
       * loop too.
       */
      int have_atleast_one_active_in_nbr;
      do {
          int *recv_flag;
    	  MPIU_CHKPMEM_MALLOC(recv_flag, int*, indegree*sizeof(int), mpi_errno, "recv_flag");
    	  for(i=0;i<indegree;i++) {
    	      recv_flag[i]=0;
          }

          have_atleast_one_active_in_nbr = 0;
          for(i = 0; i < indegree; i++) {
              //receive from the ith incoming neighbor
              if(cmn_nbh_mat->my_innbrs_bitmap[i]) {
                  mpi_errno = MPID_Irecv(&(cmn_nbh_mat->my_innbrs_bitmap[i]), 1, MPI_INT,
                                         srcs[i], 1000, comm_ptr, context_offset, &req_ptr);
                  recv_flag[i]=1;
                  if (mpi_errno) MPIU_ERR_POP(mpi_errno);
                  all_reqs[all_reqs_idx++] = req_ptr->handle;
                  have_atleast_one_active_in_nbr = 1;
              }
          }
          mpi_errno = MPIR_Waitall_impl(all_reqs_idx, all_reqs, MPI_STATUS_IGNORE);
          if (mpi_errno) MPIU_ERR_POP(mpi_errno);
          all_reqs_idx = 0; //set index back to zero for future use

          if(have_atleast_one_active_in_nbr) {
          int num_removed_nbrs_in_update=0, *rcv_num_removed_nbrs_in_update;
          MPIU_CHKPMEM_CALLOC(rcv_num_removed_nbrs_in_update, int*, indegree* sizeof(int), mpi_errno, "rcv_num_removed_nbrs_in_update");

          for(i=0;i<indegree;i++) rcv_num_removed_nbrs_in_update[i]=0;

          num_removed_nbrs_in_update= MPIR_Update_grp_frnd_mat(cmn_nbh_mat, grp_frnd_mat, comm_ptr);

          //perform Alltoall to see if a new neighbor is removed from any of the processes

          for(i=0;i<outdegree;i++) {
              if(grp_frnd_mat->is_active_nbr[i]) {
                  mpi_errno = MPID_Isend(&num_removed_nbrs_in_update, 1, MPI_INT, dests[i], 88,
        				  comm_ptr, context_offset, &req_ptr);
        	  if (mpi_errno) MPIU_ERR_POP(mpi_errno);
        	  all_reqs[all_reqs_idx++] = req_ptr->handle;
              }
          }

          for(i = 0; i < indegree; i++) { //for each of my still-active incoming neighbors
              if(cmn_nbh_mat->my_innbrs_bitmap[i] == 1) {
                  mpi_errno = MPID_Irecv(&rcv_num_removed_nbrs_in_update[i], 1, MPI_INT,
        				  srcs[i], 88, comm_ptr, context_offset, &req_ptr);
        	  recv_flag[i]=1;
        	  if (mpi_errno) MPIU_ERR_POP(mpi_errno);
        	  all_reqs[all_reqs_idx++] = req_ptr->handle;
              }
          }
          mpi_errno = MPIR_Waitall_impl(all_reqs_idx, all_reqs, MPI_STATUS_IGNORE);
          if (mpi_errno) MPIU_ERR_POP(mpi_errno);
          all_reqs_idx = 0;

          for(i=0;i<indegree;i++) recv_flag[i]=0; //for future use

          for(i=0;i<outdegree;i++) {
              if(grp_frnd_mat->is_active_nbr[i]==2) {
                  mpi_errno = MPID_Isend(&OFF, 1, MPI_INT, dests[i], 44,
        			comm_ptr, context_offset, &req_ptr);

        	  grp_frnd_mat->is_active_nbr[i]=0;
              	  if (mpi_errno) MPIU_ERR_POP(mpi_errno);
              	  all_reqs[all_reqs_idx++] = req_ptr->handle;
              } else if(grp_frnd_mat->is_active_nbr[i]==1 && num_removed_nbrs_in_update!=0) {
        	  mpi_errno = MPID_Isend(&ON, 1, MPI_INT, dests[i], 1000,
        			comm_ptr, context_offset, &req_ptr);

        	  if (mpi_errno) MPIU_ERR_POP(mpi_errno);
              	  all_reqs[all_reqs_idx++] = req_ptr->handle;
              }
          }

          for(i=0;i<indegree;i++) {
              if(rcv_num_removed_nbrs_in_update[i]) {
                  mpi_errno = MPID_Irecv(&(cmn_nbh_mat->my_innbrs_bitmap[i]), 1, MPI_INT,
        				  srcs[i], 44, comm_ptr, context_offset, &req_ptr);
        	  recv_flag[i]=1;
        	  if (mpi_errno) MPIU_ERR_POP(mpi_errno);
        	  all_reqs[all_reqs_idx++] = req_ptr->handle;
              }
          }

          mpi_errno = MPIR_Waitall_impl(all_reqs_idx, all_reqs, MPI_STATUS_IGNORE);
          if (mpi_errno) MPIU_ERR_POP(mpi_errno);
          all_reqs_idx = 0;
          }
      } while(have_atleast_one_active_in_nbr);

      int **sched_msg;
      sched_msg = NULL;
      sched_msg = MPIU_Malloc(outdegree * sizeof(int*));

      for(i = 0; i < outdegree; i++) {
    	  sched_msg[i] = MPIU_Malloc((K_vrbl+3) * sizeof(int));
      }

      for(i=0;i<outdegree;i++) {
          for(j=0;j<K_vrbl+3;j++) {
             sched_msg[i][j]=-1;
      }

      for(i = 0; i < outdegree; i++) {
          sched_msg[i][0] = cmn_nbh_mat->is_row_offloaded[i];

    	  if(!cmn_nbh_mat->is_row_offloaded[i]) {
    	      sched_msg[i][1] = cmn_nbh_mat->t;
    	      sched_msg[i][K_vrbl+2] =1; //it determines how much memory should be allocated for incom_buffer in scheduling. in other words, it determined if the received message is combined
    	  } else {  //Again we do not consider the last element for the offloaded rows
    	      /* The value of 't' does not not mean anything for offloaded rows,
    	       * but we set it to -1 to recognize offloaded neighbors later on
    	       * while building the schedule as we modify the ON/OFF field there.
    	       * I know, terrible design! but I need to get results quickly.
    	       */
    	       sched_msg[i][1] = -1;
    	       sched_msg[i][K_vrbl+2] =0;
    	  }

    	  int indx=2;
    	  for(j = 0; j < cmn_nbh_mat->t; j++) {
    	      for(k=0;k<K_vrbl;k++) {
    	          if(cmn_nbh_mat->comb_matrix[i][j].opt[k] == RECV) {
    		      sched_msg[i][K_vrbl+2] =2;
    		      sched_msg[i][1] = j;
    		      sched_msg[i][indx] = cmn_nbh_mat->comb_matrix[i][j].grp_frnds[k];
    		      indx++;
    		  }
    	      }
    	  }
      }

      //Communicating the extracted scheduling information
      int **sched_recv_buffs;
      sched_recv_buffs = MPIU_Malloc(indegree * sizeof(int*));
      for(i = 0; i < indegree; i++) {
    	  sched_recv_buffs[i] = MPIU_Malloc((K_vrbl+3) * sizeof(int));
      }

      int ii;
      for(i=0;i<indegree;i++)
    	  for(ii=0;ii<K_vrbl+3;ii++)
    		  sched_recv_buffs[i][ii]=-1;

      for(i = 0; i < outdegree; i++) {
    	  mpi_errno = MPID_Isend(sched_msg[i], (K_vrbl+3), MPI_INT,
    			  dests[i], 2000, comm_ptr, context_offset, &req_ptr);

    	  if (mpi_errno) MPIU_ERR_POP(mpi_errno);
    	  all_reqs[all_reqs_idx++] = req_ptr->handle;
    	  if (mpi_errno) MPIU_ERR_POP(mpi_errno);
      }
      for(i = 0; i < indegree; i++) {
    	  mpi_errno = MPID_Irecv(sched_recv_buffs[i], (K_vrbl+3), MPI_INT,
    			  srcs[i], 2000, comm_ptr, context_offset, &req_ptr);

    	  if (mpi_errno) MPIU_ERR_POP(mpi_errno);
    	  all_reqs[all_reqs_idx++] = req_ptr->handle;
      }
      mpi_errno = MPIR_Waitall_impl(all_reqs_idx, all_reqs, MPI_STATUS_IGNORE);
      if (mpi_errno) MPIU_ERR_POP(mpi_errno);
      all_reqs_idx = 0; //set index back to zero for future use

      for(i = 0; i < indegree; i++)
    	  sched_recv_buffs[i][K_vrbl+1]=srcs[i];

      //Attaching the received sched_recv_buff and cmn_nbh_mat to the topology of the communicator
      MPIU_CHKPMEM_MALLOC(topo_ptr->topo.dist_graph.shm_nbh_coll_patt, SHM_nbh_coll_patt*, sizeof(SHM_nbh_coll_patt), mpi_errno, "topo_ptr->topo.dist_graph.shm_nbh_coll_patt");
      topo_ptr->topo.dist_graph.shm_nbh_coll_patt->cmn_nbh_mat = cmn_nbh_mat;
      topo_ptr->topo.dist_graph.shm_nbh_coll_patt->incom_sched_mat = sched_recv_buffs;

    fn_fail:
    return 0;
}

#undef FUNCNAME
#define FUNCNAME MPI_Dist_graph_create_adjacent
#undef FCNAME
#define FCNAME MPIU_QUOTE(FUNCNAME)
/*@
MPI_Dist_graph_create_adjacent - returns a handle to a new communicator to
which the distributed graph topology information is attached.

Input Parameters:
+ comm_old - input communicator (handle)
. indegree - size of sources and sourceweights arrays (non-negative integer)
. sources - ranks of processes for which the calling process is a
            destination (array of non-negative integers)
. sourceweights - weights of the edges into the calling
                  process (array of non-negative integers or MPI_UNWEIGHTED)
. outdegree - size of destinations and destweights arrays (non-negative integer)
. destinations - ranks of processes for which the calling process is a
                 source (array of non-negative integers)
. destweights - weights of the edges out of the calling process
                (array of non-negative integers or MPI_UNWEIGHTED)
. info - hints on optimization and interpretation of weights (handle)
- reorder - the ranks may be reordered (true) or not (false) (logical)

Output Parameters:
. comm_dist_graph - communicator with distributed graph topology (handle)

.N ThreadSafe

.N Fortran

.N Errors
.N MPI_SUCCESS
.N MPI_ERR_ARG
.N MPI_ERR_OTHER
@*/
int MPI_Dist_graph_create_adjacent(MPI_Comm comm_old,
                                   int indegree,  const int sources[],
                                   const int sourceweights[],
                                   int outdegree,  const int destinations[],
                                   const int destweights[],
                                   MPI_Info info, int reorder, MPI_Comm *comm_dist_graph) {
    int       mpi_errno = MPI_SUCCESS;
    MPID_Comm *comm_ptr = NULL;
    MPID_Comm *comm_dist_graph_ptr = NULL;
    MPIR_Topology *topo_ptr = NULL;
    MPIR_Dist_graph_topology *dist_graph_ptr = NULL;
    MPIU_CHKPMEM_DECL(5);
    MPID_MPI_STATE_DECL(MPID_STATE_MPI_DIST_GRAPH_CREATE_ADJACENT);

    MPIR_ERRTEST_INITIALIZED_ORDIE();

    MPIU_THREAD_CS_ENTER(ALLFUNC,);
    MPID_MPI_FUNC_ENTER(MPID_STATE_MPI_DIST_GRAPH_CREATE_ADJACENT);

    /* Validate parameters, especially handles needing to be converted */
#   ifdef HAVE_ERROR_CHECKING
    {
        MPID_BEGIN_ERROR_CHECKS;
        {
            MPIR_ERRTEST_COMM(comm_old, mpi_errno);
            MPIR_ERRTEST_INFO_OR_NULL(info, mpi_errno);
        }
        MPID_END_ERROR_CHECKS;
    }
#   endif

    /* Convert MPI object handles to object pointers */
    MPID_Comm_get_ptr(comm_old, comm_ptr);

    /* Validate parameters and objects (post conversion) */
#   ifdef HAVE_ERROR_CHECKING
    {
        MPID_BEGIN_ERROR_CHECKS;
        {
            /* Validate comm_ptr */
            MPID_Comm_valid_ptr( comm_ptr, mpi_errno, FALSE );
            if (mpi_errno != MPI_SUCCESS) goto fn_fail;
            /* If comm_ptr is not valid, it will be reset to null */
            if (comm_ptr) {
                MPIR_ERRTEST_COMM_INTRA(comm_ptr, mpi_errno);
            }

            MPIR_ERRTEST_ARGNEG(indegree, "indegree", mpi_errno);
            MPIR_ERRTEST_ARGNEG(outdegree, "outdegree", mpi_errno);

            if (indegree > 0) {
                MPIR_ERRTEST_ARGNULL(sources, "sources", mpi_errno);
                if (sourceweights == MPI_UNWEIGHTED && destweights != MPI_UNWEIGHTED) {
                    MPIU_ERR_SET(mpi_errno, MPI_ERR_TOPOLOGY, "**unweightedboth");
                    goto fn_fail;
                }
                /* TODO check ranges for array elements too (**argarrayneg / **rankarray)*/
            }
            if (outdegree > 0) {
                MPIR_ERRTEST_ARGNULL(destinations, "destinations", mpi_errno);
                if (destweights == MPI_UNWEIGHTED && sourceweights != MPI_UNWEIGHTED) {
                    MPIU_ERR_SET(mpi_errno, MPI_ERR_TOPOLOGY, "**unweightedboth");
                    goto fn_fail;
                }
            }
            MPIR_ERRTEST_ARGNULL(comm_dist_graph, "comm_dist_graph", mpi_errno);
        }
        MPID_END_ERROR_CHECKS;
    }
#   endif /* HAVE_ERROR_CHECKING */

    /* ... body of routine ...  */

    /* Implementation based on Torsten Hoefler's reference implementation
     * attached to MPI-2.2 ticket #33. */
    *comm_dist_graph = MPI_COMM_NULL;

    /* following the spirit of the old topo interface, attributes do not
     * propagate to the new communicator (see MPI-2.1 pp. 243 line 11) */
    mpi_errno = MPIR_Comm_copy(comm_ptr, comm_ptr->local_size, &comm_dist_graph_ptr);
    if (mpi_errno) MPIU_ERR_POP(mpi_errno);

    /* Create the topology structure */
    MPIU_CHKPMEM_MALLOC(topo_ptr, MPIR_Topology *, sizeof(MPIR_Topology), mpi_errno, "topo_ptr");
    topo_ptr->kind = MPI_DIST_GRAPH;
    dist_graph_ptr = &topo_ptr->topo.dist_graph;
    dist_graph_ptr->indegree = indegree;
    dist_graph_ptr->in = NULL;
    dist_graph_ptr->in_weights = NULL;
    dist_graph_ptr->outdegree = outdegree;
    dist_graph_ptr->out = NULL;
    dist_graph_ptr->out_weights = NULL;
    dist_graph_ptr->is_weighted = (sourceweights != MPI_UNWEIGHTED);

    MPIU_CHKPMEM_MALLOC(dist_graph_ptr->in, int *, indegree*sizeof(int), mpi_errno, "dist_graph_ptr->in");
    MPIU_CHKPMEM_MALLOC(dist_graph_ptr->out, int *, outdegree*sizeof(int), mpi_errno, "dist_graph_ptr->out");
    MPIU_Memcpy(dist_graph_ptr->in, sources, indegree*sizeof(int));
    MPIU_Memcpy(dist_graph_ptr->out, destinations, outdegree*sizeof(int));

    if (dist_graph_ptr->is_weighted) {
        MPIU_CHKPMEM_MALLOC(dist_graph_ptr->in_weights, int *, indegree*sizeof(int), mpi_errno, "dist_graph_ptr->in_weights");
        MPIU_CHKPMEM_MALLOC(dist_graph_ptr->out_weights, int *, outdegree*sizeof(int), mpi_errno, "dist_graph_ptr->out_weights");
        MPIU_Memcpy(dist_graph_ptr->in_weights, sourceweights, indegree*sizeof(int));
        MPIU_Memcpy(dist_graph_ptr->out_weights, destweights, outdegree*sizeof(int));
    }

    mpi_errno = MPIR_Topology_put(comm_dist_graph_ptr, topo_ptr);
    if (mpi_errno) MPIU_ERR_POP(mpi_errno);

    MPIU_OBJ_PUBLISH_HANDLE(*comm_dist_graph, comm_dist_graph_ptr->handle);
    MPIU_CHKPMEM_COMMIT();

    double start_time, end_time;
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    start_time = MPI_Wtime();

    if(nbr_impl!=0)
    {

    int SHM_out= MPIR_Build_SHM_nbh_coll_patt(comm_dist_graph_ptr);
    if(SHM_out==-1)
    	nbr_impl=0;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    end_time = MPI_Wtime();

    if(my_rank == 0)
    	printf("Total overhead = %lf  (s)\n", end_time - start_time);

    /* ... end of body of routine ... */
  fn_exit:
    MPID_MPI_FUNC_EXIT(MPID_STATE_MPI_DIST_GRAPH_CREATE_ADJACENT);
    MPIU_THREAD_CS_EXIT(ALLFUNC,);
    return mpi_errno;

    /* --BEGIN ERROR HANDLING-- */
  fn_fail:
    MPIU_CHKPMEM_REAP();
#ifdef HAVE_ERROR_CHECKING
    mpi_errno = MPIR_Err_create_code(
        mpi_errno, MPIR_ERR_RECOVERABLE, FCNAME, __LINE__, MPI_ERR_OTHER,
        "**mpi_dist_graph_create_adjacent",
        "**mpi_dist_graph_create_adjacent %C %d %p %p %d %p %p %I %d %p",
        comm_old, indegree, sources, sourceweights,
        outdegree, destinations, destweights,
        info, reorder, comm_dist_graph);
#endif
    mpi_errno = MPIR_Err_return_comm(comm_ptr, FCNAME, mpi_errno);
    goto fn_exit;
    /* --END ERROR HANDLING-- */
}

