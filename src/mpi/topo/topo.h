/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */
/*
 *
 *  (C) 2001 by Argonne National Laboratory.

 *      See COPYRIGHT in top-level directory.
 */

#define MAX_COMB_DEGREE 50

typedef enum Operation {
    SEND,
    RECV,
    IDLE,
	HALT
}Operation;

typedef struct Comb_element {
    int *grp_frnds;
    Operation *opt;
   // int *num_entries_in_row;
}Comb_element;

typedef struct Common_nbrhood_matrix {
    int **matrix;
    int **outnbrs_innbrs_bitmap; //Not currently an actual bitmap
    int *row_sizes; //indegree of each outgoing neighbor
    int *ignore_row; //This determines whether the neighbor is considered in the process of finding friends.
    int *is_row_offloaded; //This determines whether a neighbor should be communicated with at the end. I.e., whether it has been offloaded or not.
    int *my_innbrs_bitmap; //Not currently an actual bitmap
    int indegree;
    int num_rows; //outdegree
    int num_elements; //sum of all row_sizes
    int t;
    Comb_element **comb_matrix; //For each outgoing neighbor, represents the list of ranks with whom I should combine my message before sending it out.
   // int *comb_matrix_num_entries_in_row;
}Common_nbrhood_matrix;


//SMGM
typedef struct Scheduler_In_Info{
	int incom_flag;
	int time_stp;
	int *combnd_prcs;
}Scheduler_In_Info;


typedef struct Individual_Friendship_Matrix {
    int **frndship_bit_map;
    int *friends;
    int *nbrs;
    int numfriends;
}Individual_Friendship_Matrix;

typedef struct Group_Friendship_Matrix {
    int **grpfrnd_bit_map;
    int **grpfriends;
    int *num_actv_grp_of_nbr;
   // int *grpnbrs;
    int *is_active_grpfrnd; //determines is the corresponding group friends are active. Initially all group friends are active
    int *is_active_nbr; //determines is the corresponding neighbor is active. Initially all neighbors are active. if it is set to 2, it means that in is deactivated in update function but it should send OFF to its incoming neighbors before being set to zero.
    int *num_cmn_brs; //number of common neighbors for each group friends
    long num_grp_frnds; //number of active group friends
    long Total_num_grp_frnds; //total number of group friends
}Group_Friendship_Matrix;

typedef struct F_Group_Friends { //this is F in the algorithm
    int **grpfriends;
    int *num_cmn_brs; //number of common neighbors for each group friends
    long num_avail_grpfrnds;
}F_Group_Friends;

typedef struct Chosen_Group_Friends { //this is R in the algorithm
	int   *source;
    int **grpfriends;
    int num_chosen_grpfrnds;
}Chosen_Group_Friends;

typedef struct Transfered_Request_s {
	int  msgtype;  // it specifies whether the request is drop(=0) or grant(=1)
    int *grp;
}Transfered_Request;

typedef struct Final_Group_Friends {
	int  is_grp;
    int *grp_frnds;
}Final_Group_Friends;

typedef struct SHM_nbh_coll_patt {
    Common_nbrhood_matrix *cmn_nbh_mat;
    int **incom_sched_mat;

}SHM_nbh_coll_patt;

//SMGM

typedef struct Common_neighbor{
    int rank;
    int index;
}Common_neighbor;

typedef struct MPIR_Graph_topology {
  int nnodes;
  int nedges;
  int *index;
  int *edges;
} MPIR_Graph_topology;

typedef struct MPIR_Cart_topology {
  int nnodes;     /* Product of dims[*], gives the size of the topology */
  int ndims;
  int *dims;
  int *periodic;
  int *position;
} MPIR_Cart_topology;

typedef struct MPIR_Dist_graph_topology {
    int indegree;
    int *in;
    int *in_weights;
    int outdegree;
    int *out;
    int *out_weights;
    int is_weighted;
    //SHM added
    SHM_nbh_coll_patt *shm_nbh_coll_patt;
    MPID_Sched_t shm_nbh_coll_sched;
    void *sched_mem_to_free[2];
} MPIR_Dist_graph_topology;

typedef struct MPIR_Topology { 
  MPIR_Topo_type kind;
  union topo { 
    MPIR_Graph_topology graph;
    MPIR_Cart_topology  cart;
    MPIR_Dist_graph_topology dist_graph;
  } topo;
} MPIR_Topology;

MPIR_Topology *MPIR_Topology_get( MPID_Comm * );
int MPIR_Topology_put( MPID_Comm *, MPIR_Topology * );
int MPIR_Cart_create( MPID_Comm *, int, const int [], 
		      const int [], int, MPI_Comm * );
int MPIR_Graph_create( MPID_Comm *, int, 
		       const int[], const int[], int, 
		       MPI_Comm *);
int MPIR_Dims_create( int, int, int * );
int MPIR_Graph_map( const MPID_Comm *, int, const int[], const int[], int* );
int MPIR_Cart_map( const MPID_Comm *, int, const int[],  const int[], int* );

/* Returns the canonicalized count of neighbors for the given topology as though
 * MPI_Dist_graph_neighbors_count were called with a distributed graph topology,
 * even if the given topology is actually Cartesian or Graph.  Useful for
 * implementing neighborhood collective operations. */
int MPIR_Topo_canon_nhb_count(MPID_Comm *comm_ptr, int *indegree, int *outdegree, int *weighted);

/* Returns the canonicalized list of neighbors for a given topology, separated
 * into inbound and outbound edges.  Equivalent to MPI_Dist_graph_neighbors but
 * works for any topology type by canonicalizing according to the rules in
 * Section 7.6 of the MPI-3.0 standard. */
int MPIR_Topo_canon_nhb(MPID_Comm *comm_ptr,
                        int indegree, int sources[], int inweights[],
                        int outdegree, int dests[], int outweights[]);

#define MAX_CART_DIM 16
