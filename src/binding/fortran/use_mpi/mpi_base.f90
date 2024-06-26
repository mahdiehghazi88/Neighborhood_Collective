!     -*- Mode: Fortran; -*-
!  (C) 2008 by Argonne National Laboratory.
!      See COPYRIGHT in top-level directory.
       MODULE MPI_BASE
       IMPLICIT NONE
!      This module was created by the script buildiface
       INTERFACE
      SUBROUTINE MPI_TYPE_CREATE_DARRAY(size,rank,ndims,array_of_gsizes,&
                  array_of_distribs,array_of_dargs,array_of_psizes,order,&
                  oldtype,newtype,ierror)
           INTEGER size
           INTEGER rank
           INTEGER ndims
           INTEGER array_of_gsizes(*)
           INTEGER array_of_distribs(*)
           INTEGER array_of_dargs(*)
           INTEGER array_of_psizes(*)
           INTEGER order
           INTEGER oldtype
           INTEGER newtype
           INTEGER ierror
      END SUBROUTINE MPI_TYPE_CREATE_DARRAY

      SUBROUTINE MPI_TYPE_EXTENT(datatype,extent,ierror)
           INTEGER datatype
           INTEGER extent
           INTEGER ierror
      END SUBROUTINE MPI_TYPE_EXTENT

      SUBROUTINE MPI_TYPE_GET_NAME(datatype,type_name,resultlen,ierror)
           INTEGER datatype
           CHARACTER (LEN=*) type_name
           INTEGER resultlen
           INTEGER ierror
      END SUBROUTINE MPI_TYPE_GET_NAME

      SUBROUTINE MPI_WIN_LOCK(lock_type,rank,assert,win,ierror)
           INTEGER lock_type
           INTEGER rank
           INTEGER assert
           INTEGER win
           INTEGER ierror
      END SUBROUTINE MPI_WIN_LOCK

      SUBROUTINE MPI_COMM_SPLIT(comm,color,key,newcomm,ierror)
           INTEGER comm
           INTEGER color
           INTEGER key
           INTEGER newcomm
           INTEGER ierror
      END SUBROUTINE MPI_COMM_SPLIT

      SUBROUTINE MPI_WIN_COMPLETE(win,ierror)
           INTEGER win
           INTEGER ierror
      END SUBROUTINE MPI_WIN_COMPLETE

      SUBROUTINE MPI_GROUP_SIZE(group,size,ierror)
           INTEGER group
           INTEGER size
           INTEGER ierror
      END SUBROUTINE MPI_GROUP_SIZE

      SUBROUTINE MPI_FILE_GET_ERRHANDLER(file,errhandler,ierror)
           INTEGER file
           INTEGER errhandler
           INTEGER ierror
      END SUBROUTINE MPI_FILE_GET_ERRHANDLER

      SUBROUTINE MPI_REQUEST_FREE(request,ierror)
           INTEGER request
           INTEGER ierror
      END SUBROUTINE MPI_REQUEST_FREE

      SUBROUTINE MPI_TYPE_CREATE_HINDEXED_BLOCK(count,blocklength,&
                  array_of_displacements,oldtype,newtype,ierror)
           USE MPI_CONSTANTS,ONLY:MPI_ADDRESS_KIND
           INTEGER count
           INTEGER blocklength
           INTEGER(KIND=MPI_ADDRESS_KIND) array_of_displacements(*)
           INTEGER oldtype
           INTEGER newtype
           INTEGER ierror
      END SUBROUTINE MPI_TYPE_CREATE_HINDEXED_BLOCK

      SUBROUTINE MPI_BARRIER(comm,ierror)
           INTEGER comm
           INTEGER ierror
      END SUBROUTINE MPI_BARRIER

      SUBROUTINE MPI_TYPE_COMMIT(datatype,ierror)
           INTEGER datatype
           INTEGER ierror
      END SUBROUTINE MPI_TYPE_COMMIT

      SUBROUTINE MPI_GROUP_RANGE_EXCL(group,n,ranges,newgroup,ierror)
           INTEGER group
           INTEGER n
           INTEGER ranges(3,*)
           INTEGER newgroup
           INTEGER ierror
      END SUBROUTINE MPI_GROUP_RANGE_EXCL

      SUBROUTINE MPI_REQUEST_GET_STATUS(request,flag,status,ierror)
           USE MPI_CONSTANTS,ONLY:MPI_STATUS_SIZE
           INTEGER request
           LOGICAL flag
           INTEGER status(MPI_STATUS_SIZE)
           INTEGER ierror
      END SUBROUTINE MPI_REQUEST_GET_STATUS

      SUBROUTINE MPI_COMM_SPAWN_MULTIPLE(count,array_of_commands,array_of_argv,&
                  array_of_maxprocs,array_of_info,root,comm,intercomm,&
                  array_of_errcodes,ierror)
           INTEGER count
           CHARACTER (LEN=*) array_of_commands(*)
           CHARACTER (LEN=*) array_of_argv(count,*)
           INTEGER array_of_maxprocs(*)
           INTEGER array_of_info(*)
           INTEGER root
           INTEGER comm
           INTEGER intercomm
           INTEGER array_of_errcodes(*)
           INTEGER ierror
      END SUBROUTINE MPI_COMM_SPAWN_MULTIPLE

      SUBROUTINE MPI_TYPE_GET_EXTENT(datatype,lb,extent,ierror)
           USE MPI_CONSTANTS,ONLY:MPI_ADDRESS_KIND
           INTEGER datatype
           INTEGER(KIND=MPI_ADDRESS_KIND) lb
           INTEGER(KIND=MPI_ADDRESS_KIND) extent
           INTEGER ierror
      END SUBROUTINE MPI_TYPE_GET_EXTENT

      SUBROUTINE MPI_INFO_GET_VALUELEN(info,key,valuelen,flag,ierror)
           INTEGER info
           CHARACTER (LEN=*) key
           INTEGER valuelen
           LOGICAL flag
           INTEGER ierror
      END SUBROUTINE MPI_INFO_GET_VALUELEN

      SUBROUTINE MPI_OP_CREATE(user_fn,commute,op,ierror)
           EXTERNAL user_fn
           LOGICAL commute
           INTEGER op
           INTEGER ierror
      END SUBROUTINE MPI_OP_CREATE

      SUBROUTINE MPI_TYPE_CREATE_STRUCT(count,array_of_blocklengths,&
                  array_of_displacements,array_of_types,newtype,ierror)
           USE MPI_CONSTANTS,ONLY:MPI_ADDRESS_KIND
           INTEGER count
           INTEGER array_of_blocklengths(*)
           INTEGER(KIND=MPI_ADDRESS_KIND) array_of_displacements(*)
           INTEGER array_of_types(*)
           INTEGER newtype
           INTEGER ierror
      END SUBROUTINE MPI_TYPE_CREATE_STRUCT

      SUBROUTINE MPI_WIN_GET_GROUP(win,group,ierror)
           INTEGER win
           INTEGER group
           INTEGER ierror
      END SUBROUTINE MPI_WIN_GET_GROUP

      SUBROUTINE MPI_GROUP_COMPARE(group1,group2,result,ierror)
           INTEGER group1
           INTEGER group2
           INTEGER result
           INTEGER ierror
      END SUBROUTINE MPI_GROUP_COMPARE

      SUBROUTINE MPI_CART_SHIFT(comm,direction,disp,rank_source,rank_dest,&
                  ierror)
           INTEGER comm
           INTEGER direction
           INTEGER disp
           INTEGER rank_source
           INTEGER rank_dest
           INTEGER ierror
      END SUBROUTINE MPI_CART_SHIFT

      SUBROUTINE MPI_COMM_GROUP(comm,group,ierror)
           INTEGER comm
           INTEGER group
           INTEGER ierror
      END SUBROUTINE MPI_COMM_GROUP

      SUBROUTINE MPI_WIN_CALL_ERRHANDLER(win,errorcode,ierror)
           INTEGER win
           INTEGER errorcode
           INTEGER ierror
      END SUBROUTINE MPI_WIN_CALL_ERRHANDLER

      SUBROUTINE MPI_LOOKUP_NAME(service_name,info,port_name,ierror)
           CHARACTER (LEN=*) service_name
           INTEGER info
           CHARACTER (LEN=*) port_name
           INTEGER ierror
      END SUBROUTINE MPI_LOOKUP_NAME

      SUBROUTINE MPI_INFO_FREE(info,ierror)
           INTEGER info
           INTEGER ierror
      END SUBROUTINE MPI_INFO_FREE

      SUBROUTINE MPI_GRAPH_GET(comm,maxindex,maxedges,indx,edges,ierror)
           INTEGER comm
           INTEGER maxindex
           INTEGER maxedges
           INTEGER indx(*)
           INTEGER edges(*)
           INTEGER ierror
      END SUBROUTINE MPI_GRAPH_GET

      SUBROUTINE MPI_STATUS_SET_ELEMENTS(status,datatype,count,ierror)
           USE MPI_CONSTANTS,ONLY:MPI_STATUS_SIZE
           INTEGER status(MPI_STATUS_SIZE)
           INTEGER datatype
           INTEGER count
           INTEGER ierror
      END SUBROUTINE MPI_STATUS_SET_ELEMENTS

      SUBROUTINE MPI_COMM_SET_INFO(comm,info,ierror)
           INTEGER comm
           INTEGER info
           INTEGER ierror
      END SUBROUTINE MPI_COMM_SET_INFO

      SUBROUTINE MPI_WIN_FREE(win,ierror)
           INTEGER win
           INTEGER ierror
      END SUBROUTINE MPI_WIN_FREE

      SUBROUTINE MPI_PACK_EXTERNAL_SIZE(datarep,incount,datatype,size,ierror)
           USE MPI_CONSTANTS,ONLY:MPI_ADDRESS_KIND
           CHARACTER (LEN=*) datarep
           INTEGER incount
           INTEGER datatype
           INTEGER(KIND=MPI_ADDRESS_KIND) size
           INTEGER ierror
      END SUBROUTINE MPI_PACK_EXTERNAL_SIZE

      SUBROUTINE MPI_OPEN_PORT(info,port_name,ierror)
           INTEGER info
           CHARACTER (LEN=*) port_name
           INTEGER ierror
      END SUBROUTINE MPI_OPEN_PORT

      SUBROUTINE MPI_COMM_SPLIT_TYPE(comm,split_type,key,info,newcomm,ierror)
           INTEGER comm
           INTEGER split_type
           INTEGER key
           INTEGER info
           INTEGER newcomm
           INTEGER ierror
      END SUBROUTINE MPI_COMM_SPLIT_TYPE

      SUBROUTINE MPI_FILE_CREATE_ERRHANDLER(file_errhandler_fn,errhandler,&
                  ierror)
           INTERFACE 
       SUBROUTINE file_errhandler_fn(vv0,vv1)
       INTEGER vv0,vv1
       END SUBROUTINE
       END INTERFACE
           INTEGER errhandler
           INTEGER ierror
      END SUBROUTINE MPI_FILE_CREATE_ERRHANDLER

      SUBROUTINE MPI_COMM_CREATE_GROUP(comm,group,tag,newcomm,ierror)
           INTEGER comm
           INTEGER group
           INTEGER tag
           INTEGER newcomm
           INTEGER ierror
      END SUBROUTINE MPI_COMM_CREATE_GROUP

      SUBROUTINE MPI_ATTR_DELETE(comm,keyval,ierror)
           INTEGER comm
           INTEGER keyval
           INTEGER ierror
      END SUBROUTINE MPI_ATTR_DELETE

      SUBROUTINE MPI_WIN_LOCK_ALL(assert,win,ierror)
           INTEGER assert
           INTEGER win
           INTEGER ierror
      END SUBROUTINE MPI_WIN_LOCK_ALL

      SUBROUTINE MPI_IBARRIER(comm,request,ierror)
           INTEGER comm
           INTEGER request
           INTEGER ierror
      END SUBROUTINE MPI_IBARRIER

      SUBROUTINE MPI_TYPE_GET_CONTENTS(datatype,max_integers,max_addresses,&
                  max_datatypes,array_of_integers,array_of_addresses,&
                  array_of_datatypes,ierror)
           USE MPI_CONSTANTS,ONLY:MPI_ADDRESS_KIND
           INTEGER datatype
           INTEGER max_integers
           INTEGER max_addresses
           INTEGER max_datatypes
           INTEGER array_of_integers(*)
           INTEGER(KIND=MPI_ADDRESS_KIND) array_of_addresses(*)
           INTEGER array_of_datatypes(*)
           INTEGER ierror
      END SUBROUTINE MPI_TYPE_GET_CONTENTS

      SUBROUTINE MPI_WIN_FLUSH_ALL(win,ierror)
           INTEGER win
           INTEGER ierror
      END SUBROUTINE MPI_WIN_FLUSH_ALL

      SUBROUTINE MPI_START(request,ierror)
           INTEGER request
           INTEGER ierror
      END SUBROUTINE MPI_START

      SUBROUTINE MPI_FINALIZE(ierror)
           INTEGER ierror
      END SUBROUTINE MPI_FINALIZE

      SUBROUTINE MPI_DIST_GRAPH_NEIGHBORS_COUNT(comm,indegree,outdegree,&
                  weighted,ierror)
           INTEGER comm
           INTEGER indegree
           INTEGER outdegree
           LOGICAL weighted
           INTEGER ierror
      END SUBROUTINE MPI_DIST_GRAPH_NEIGHBORS_COUNT

      SUBROUTINE MPI_COMM_GET_PARENT(parent,ierror)
           INTEGER parent
           INTEGER ierror
      END SUBROUTINE MPI_COMM_GET_PARENT

      SUBROUTINE MPI_FINALIZED(flag,ierror)
           LOGICAL flag
           INTEGER ierror
      END SUBROUTINE MPI_FINALIZED

      SUBROUTINE MPI_INTERCOMM_MERGE(intercomm,high,newintracomm,ierror)
           INTEGER intercomm
           LOGICAL high
           INTEGER newintracomm
           INTEGER ierror
      END SUBROUTINE MPI_INTERCOMM_MERGE

      SUBROUTINE MPI_INFO_GET_NTHKEY(info,n,key,ierror)
           INTEGER info
           INTEGER n
           CHARACTER (LEN=*) key
           INTEGER ierror
      END SUBROUTINE MPI_INFO_GET_NTHKEY

      SUBROUTINE MPI_TYPE_MATCH_SIZE(typeclass,size,datatype,ierror)
           INTEGER typeclass
           INTEGER size
           INTEGER datatype
           INTEGER ierror
      END SUBROUTINE MPI_TYPE_MATCH_SIZE

      SUBROUTINE MPI_INITIALIZED(flag,ierror)
           LOGICAL flag
           INTEGER ierror
      END SUBROUTINE MPI_INITIALIZED

      SUBROUTINE MPI_TYPE_CONTIGUOUS(count,oldtype,newtype,ierror)
           INTEGER count
           INTEGER oldtype
           INTEGER newtype
           INTEGER ierror
      END SUBROUTINE MPI_TYPE_CONTIGUOUS

      SUBROUTINE MPI_TYPE_UB(datatype,displacement,ierror)
           INTEGER datatype
           INTEGER displacement
           INTEGER ierror
      END SUBROUTINE MPI_TYPE_UB

      SUBROUTINE MPI_INFO_DUP(info,newinfo,ierror)
           INTEGER info
           INTEGER newinfo
           INTEGER ierror
      END SUBROUTINE MPI_INFO_DUP

      SUBROUTINE MPI_WIN_DELETE_ATTR(win,win_keyval,ierror)
           INTEGER win
           INTEGER win_keyval
           INTEGER ierror
      END SUBROUTINE MPI_WIN_DELETE_ATTR

      SUBROUTINE MPI_INFO_GET_NKEYS(info,nkeys,ierror)
           INTEGER info
           INTEGER nkeys
           INTEGER ierror
      END SUBROUTINE MPI_INFO_GET_NKEYS

      SUBROUTINE MPI_WIN_CREATE_DYNAMIC(info,comm,win,ierror)
           INTEGER info
           INTEGER comm
           INTEGER win
           INTEGER ierror
      END SUBROUTINE MPI_WIN_CREATE_DYNAMIC

      SUBROUTINE MPI_GROUP_EXCL(group,n,ranks,newgroup,ierror)
           INTEGER group
           INTEGER n
           INTEGER ranks(*)
           INTEGER newgroup
           INTEGER ierror
      END SUBROUTINE MPI_GROUP_EXCL

      SUBROUTINE MPI_WIN_GET_INFO(win,info_used,ierror)
           INTEGER win
           INTEGER info_used
           INTEGER ierror
      END SUBROUTINE MPI_WIN_GET_INFO

      SUBROUTINE MPI_COMM_DELETE_ATTR(comm,comm_keyval,ierror)
           INTEGER comm
           INTEGER comm_keyval
           INTEGER ierror
      END SUBROUTINE MPI_COMM_DELETE_ATTR

      SUBROUTINE MPI_GET_COUNT(status,datatype,count,ierror)
           USE MPI_CONSTANTS,ONLY:MPI_STATUS_SIZE
           INTEGER status(MPI_STATUS_SIZE)
           INTEGER datatype
           INTEGER count
           INTEGER ierror
      END SUBROUTINE MPI_GET_COUNT

      SUBROUTINE MPI_ADD_ERROR_CLASS(errorclass,ierror)
           INTEGER errorclass
           INTEGER ierror
      END SUBROUTINE MPI_ADD_ERROR_CLASS

      SUBROUTINE MPI_COMM_SET_NAME(comm,comm_name,ierror)
           INTEGER comm
           CHARACTER (LEN=*) comm_name
           INTEGER ierror
      END SUBROUTINE MPI_COMM_SET_NAME

      SUBROUTINE MPI_COMM_DISCONNECT(comm,ierror)
           INTEGER comm
           INTEGER ierror
      END SUBROUTINE MPI_COMM_DISCONNECT

      SUBROUTINE MPI_ADD_ERROR_CODE(errorclass,errorcode,ierror)
           INTEGER errorclass
           INTEGER errorcode
           INTEGER ierror
      END SUBROUTINE MPI_ADD_ERROR_CODE

      SUBROUTINE MPI_COMM_GET_ERRHANDLER(comm,errhandler,ierror)
           INTEGER comm
           INTEGER errhandler
           INTEGER ierror
      END SUBROUTINE MPI_COMM_GET_ERRHANDLER

      SUBROUTINE MPI_COMM_CREATE(comm,group,newcomm,ierror)
           INTEGER comm
           INTEGER group
           INTEGER newcomm
           INTEGER ierror
      END SUBROUTINE MPI_COMM_CREATE

      SUBROUTINE MPI_COMM_REMOTE_SIZE(comm,size,ierror)
           INTEGER comm
           INTEGER size
           INTEGER ierror
      END SUBROUTINE MPI_COMM_REMOTE_SIZE

      SUBROUTINE MPI_PROBE(source,tag,comm,status,ierror)
           USE MPI_CONSTANTS,ONLY:MPI_STATUS_SIZE
           INTEGER source
           INTEGER tag
           INTEGER comm
           INTEGER status(MPI_STATUS_SIZE)
           INTEGER ierror
      END SUBROUTINE MPI_PROBE

      SUBROUTINE MPI_WIN_WAIT(win,ierror)
           INTEGER win
           INTEGER ierror
      END SUBROUTINE MPI_WIN_WAIT

      SUBROUTINE MPI_TYPE_CREATE_SUBARRAY(ndims,array_of_sizes,&
                  array_of_subsizes,array_of_starts,order,oldtype,newtype,&
                  ierror)
           INTEGER ndims
           INTEGER array_of_sizes(*)
           INTEGER array_of_subsizes(*)
           INTEGER array_of_starts(*)
           INTEGER order
           INTEGER oldtype
           INTEGER newtype
           INTEGER ierror
      END SUBROUTINE MPI_TYPE_CREATE_SUBARRAY

      SUBROUTINE MPI_TYPE_SIZE_X(datatype,size,ierror)
           USE MPI_CONSTANTS,ONLY:MPI_COUNT_KIND
           INTEGER datatype
           INTEGER(KIND=MPI_COUNT_KIND) size
           INTEGER ierror
      END SUBROUTINE MPI_TYPE_SIZE_X

      SUBROUTINE MPI_TYPE_FREE(datatype,ierror)
           INTEGER datatype
           INTEGER ierror
      END SUBROUTINE MPI_TYPE_FREE

      SUBROUTINE MPI_GRAPHDIMS_GET(comm,nnodes,nedges,ierror)
           INTEGER comm
           INTEGER nnodes
           INTEGER nedges
           INTEGER ierror
      END SUBROUTINE MPI_GRAPHDIMS_GET

      SUBROUTINE MPI_FILE_CALL_ERRHANDLER(fh,errorcode,ierror)
           INTEGER fh
           INTEGER errorcode
           INTEGER ierror
      END SUBROUTINE MPI_FILE_CALL_ERRHANDLER

      SUBROUTINE MPI_TYPE_DELETE_ATTR(datatype,type_keyval,ierror)
           INTEGER datatype
           INTEGER type_keyval
           INTEGER ierror
      END SUBROUTINE MPI_TYPE_DELETE_ATTR

      SUBROUTINE MPI_TYPE_CREATE_HINDEXED(count,array_of_blocklengths,&
                  array_of_displacements,oldtype,newtype,ierror)
           USE MPI_CONSTANTS,ONLY:MPI_ADDRESS_KIND
           INTEGER count
           INTEGER array_of_blocklengths(*)
           INTEGER(KIND=MPI_ADDRESS_KIND) array_of_displacements(*)
           INTEGER oldtype
           INTEGER newtype
           INTEGER ierror
      END SUBROUTINE MPI_TYPE_CREATE_HINDEXED

      SUBROUTINE MPI_STATUS_SET_ELEMENTS_X(status,datatype,count,ierror)
           USE MPI_CONSTANTS,ONLY:MPI_COUNT_KIND, MPI_STATUS_SIZE
           INTEGER status(MPI_STATUS_SIZE)
           INTEGER datatype
           INTEGER(KIND=MPI_COUNT_KIND) count
           INTEGER ierror
      END SUBROUTINE MPI_STATUS_SET_ELEMENTS_X

      SUBROUTINE MPI_TYPE_INDEXED(count,array_of_blocklengths,&
                  array_of_displacements,oldtype,newtype,ierror)
           INTEGER count
           INTEGER array_of_blocklengths(*)
           INTEGER array_of_displacements(*)
           INTEGER oldtype
           INTEGER newtype
           INTEGER ierror
      END SUBROUTINE MPI_TYPE_INDEXED

      SUBROUTINE MPI_GRAPH_NEIGHBORS_COUNT(comm,rank,nneighbors,ierror)
           INTEGER comm
           INTEGER rank
           INTEGER nneighbors
           INTEGER ierror
      END SUBROUTINE MPI_GRAPH_NEIGHBORS_COUNT

      SUBROUTINE MPI_KEYVAL_FREE(keyval,ierror)
           INTEGER keyval
           INTEGER ierror
      END SUBROUTINE MPI_KEYVAL_FREE

      SUBROUTINE MPI_GROUP_DIFFERENCE(group1,group2,newgroup,ierror)
           INTEGER group1
           INTEGER group2
           INTEGER newgroup
           INTEGER ierror
      END SUBROUTINE MPI_GROUP_DIFFERENCE

      SUBROUTINE MPI_COMM_DUP(comm,newcomm,ierror)
           INTEGER comm
           INTEGER newcomm
           INTEGER ierror
      END SUBROUTINE MPI_COMM_DUP

      SUBROUTINE MPI_ERROR_CLASS(errorcode,errorclass,ierror)
           INTEGER errorcode
           INTEGER errorclass
           INTEGER ierror
      END SUBROUTINE MPI_ERROR_CLASS

      SUBROUTINE MPI_COMM_FREE_KEYVAL(comm_keyval,ierror)
           INTEGER comm_keyval
           INTEGER ierror
      END SUBROUTINE MPI_COMM_FREE_KEYVAL

      SUBROUTINE MPI_GET_LIBRARY_VERSION(version,resultlen,ierror)
           CHARACTER (LEN=*) version
           INTEGER resultlen
           INTEGER ierror
      END SUBROUTINE MPI_GET_LIBRARY_VERSION

      SUBROUTINE MPI_WIN_FLUSH_LOCAL(rank,win,ierror)
           INTEGER rank
           INTEGER win
           INTEGER ierror
      END SUBROUTINE MPI_WIN_FLUSH_LOCAL

      SUBROUTINE MPI_GROUP_INTERSECTION(group1,group2,newgroup,ierror)
           INTEGER group1
           INTEGER group2
           INTEGER newgroup
           INTEGER ierror
      END SUBROUTINE MPI_GROUP_INTERSECTION

      SUBROUTINE MPI_CARTDIM_GET(comm,ndims,ierror)
           INTEGER comm
           INTEGER ndims
           INTEGER ierror
      END SUBROUTINE MPI_CARTDIM_GET

      SUBROUTINE MPI_WIN_GET_ERRHANDLER(win,errhandler,ierror)
           INTEGER win
           INTEGER errhandler
           INTEGER ierror
      END SUBROUTINE MPI_WIN_GET_ERRHANDLER

      SUBROUTINE MPI_CANCEL(request,ierror)
           INTEGER request
           INTEGER ierror
      END SUBROUTINE MPI_CANCEL

      SUBROUTINE MPI_WIN_POST(group,assert,win,ierror)
           INTEGER group
           INTEGER assert
           INTEGER win
           INTEGER ierror
      END SUBROUTINE MPI_WIN_POST

      SUBROUTINE MPI_TEST_CANCELLED(status,flag,ierror)
           USE MPI_CONSTANTS,ONLY:MPI_STATUS_SIZE
           INTEGER status(MPI_STATUS_SIZE)
           LOGICAL flag
           INTEGER ierror
      END SUBROUTINE MPI_TEST_CANCELLED

      SUBROUTINE MPI_ADD_ERROR_STRING(errorcode,string,ierror)
           INTEGER errorcode
           CHARACTER (LEN=*) string
           INTEGER ierror
      END SUBROUTINE MPI_ADD_ERROR_STRING

      SUBROUTINE MPI_GET_ELEMENTS(status,datatype,count,ierror)
           USE MPI_CONSTANTS,ONLY:MPI_STATUS_SIZE
           INTEGER status(MPI_STATUS_SIZE)
           INTEGER datatype
           INTEGER count
           INTEGER ierror
      END SUBROUTINE MPI_GET_ELEMENTS

      SUBROUTINE MPI_PACK_SIZE(incount,datatype,comm,size,ierror)
           INTEGER incount
           INTEGER datatype
           INTEGER comm
           INTEGER size
           INTEGER ierror
      END SUBROUTINE MPI_PACK_SIZE

      SUBROUTINE MPI_ERRHANDLER_GET(comm,errhandler,ierror)
           INTEGER comm
           INTEGER errhandler
           INTEGER ierror
      END SUBROUTINE MPI_ERRHANDLER_GET

      SUBROUTINE MPI_WIN_SET_INFO(win,info,ierror)
           INTEGER win
           INTEGER info
           INTEGER ierror
      END SUBROUTINE MPI_WIN_SET_INFO

      SUBROUTINE MPI_OP_COMMUTATIVE(op,commute,ierror)
           INTEGER op
           LOGICAL commute
           INTEGER ierror
      END SUBROUTINE MPI_OP_COMMUTATIVE

      SUBROUTINE MPI_TYPE_LB(datatype,displacement,ierror)
           INTEGER datatype
           INTEGER displacement
           INTEGER ierror
      END SUBROUTINE MPI_TYPE_LB

      SUBROUTINE MPI_GROUP_RANGE_INCL(group,n,ranges,newgroup,ierror)
           INTEGER group
           INTEGER n
           INTEGER ranges(3,*)
           INTEGER newgroup
           INTEGER ierror
      END SUBROUTINE MPI_GROUP_RANGE_INCL

      SUBROUTINE MPI_TYPE_GET_TRUE_EXTENT(datatype,true_lb,true_extent,ierror)
           USE MPI_CONSTANTS,ONLY:MPI_ADDRESS_KIND
           INTEGER datatype
           INTEGER(KIND=MPI_ADDRESS_KIND) true_lb
           INTEGER(KIND=MPI_ADDRESS_KIND) true_extent
           INTEGER ierror
      END SUBROUTINE MPI_TYPE_GET_TRUE_EXTENT

      SUBROUTINE MPI_IS_THREAD_MAIN(flag,ierror)
           LOGICAL flag
           INTEGER ierror
      END SUBROUTINE MPI_IS_THREAD_MAIN

      SUBROUTINE MPI_WIN_FREE_KEYVAL(win_keyval,ierror)
           INTEGER win_keyval
           INTEGER ierror
      END SUBROUTINE MPI_WIN_FREE_KEYVAL

      SUBROUTINE MPI_QUERY_THREAD(provided,ierror)
           INTEGER provided
           INTEGER ierror
      END SUBROUTINE MPI_QUERY_THREAD

      SUBROUTINE MPI_ERRHANDLER_CREATE(function,errhandler,ierror)
           INTERFACE 
       SUBROUTINE function(vv0,vv1)
       INTEGER vv0,vv1
       END SUBROUTINE
       END INTERFACE
           INTEGER errhandler
           INTEGER ierror
      END SUBROUTINE MPI_ERRHANDLER_CREATE

      SUBROUTINE MPI_TYPE_GET_TRUE_EXTENT_X(datatype,lb,extent,ierror)
           USE MPI_CONSTANTS,ONLY:MPI_COUNT_KIND
           INTEGER datatype
           INTEGER(KIND=MPI_COUNT_KIND) lb
           INTEGER(KIND=MPI_COUNT_KIND) extent
           INTEGER ierror
      END SUBROUTINE MPI_TYPE_GET_TRUE_EXTENT_X

      SUBROUTINE MPI_IMPROBE(source,tag,comm,flag,message,status,ierror)
           USE MPI_CONSTANTS,ONLY:MPI_STATUS_SIZE
           INTEGER source
           INTEGER tag
           INTEGER comm
           LOGICAL flag
           INTEGER message
           INTEGER status(MPI_STATUS_SIZE)
           INTEGER ierror
      END SUBROUTINE MPI_IMPROBE

      SUBROUTINE MPI_COMM_REMOTE_GROUP(comm,group,ierror)
           INTEGER comm
           INTEGER group
           INTEGER ierror
      END SUBROUTINE MPI_COMM_REMOTE_GROUP

      SUBROUTINE MPI_COMM_COMPARE(comm1,comm2,result,ierror)
           INTEGER comm1
           INTEGER comm2
           INTEGER result
           INTEGER ierror
      END SUBROUTINE MPI_COMM_COMPARE

      SUBROUTINE MPI_INFO_GET(info,key,valuelen,value,flag,ierror)
           INTEGER info
           CHARACTER (LEN=*) key
           INTEGER valuelen
           CHARACTER (LEN=*) value
           LOGICAL flag
           INTEGER ierror
      END SUBROUTINE MPI_INFO_GET

      SUBROUTINE MPI_TYPE_VECTOR(count,blocklength,stride,oldtype,newtype,&
                  ierror)
           INTEGER count
           INTEGER blocklength
           INTEGER stride
           INTEGER oldtype
           INTEGER newtype
           INTEGER ierror
      END SUBROUTINE MPI_TYPE_VECTOR

      SUBROUTINE MPI_COMM_DUP_WITH_INFO(comm,info,newcomm,ierror)
           INTEGER comm
           INTEGER info
           INTEGER newcomm
           INTEGER ierror
      END SUBROUTINE MPI_COMM_DUP_WITH_INFO

      SUBROUTINE MPI_WIN_SET_ERRHANDLER(win,errhandler,ierror)
           INTEGER win
           INTEGER errhandler
           INTEGER ierror
      END SUBROUTINE MPI_WIN_SET_ERRHANDLER

      SUBROUTINE MPI_COMM_SPAWN(command,argv,maxprocs,info,root,comm,intercomm,&
                  array_of_errcodes,ierror)
           CHARACTER (LEN=*) command
           CHARACTER (LEN=*) argv(*)
           INTEGER maxprocs
           INTEGER info
           INTEGER root
           INTEGER comm
           INTEGER intercomm
           INTEGER array_of_errcodes(*)
           INTEGER ierror
      END SUBROUTINE MPI_COMM_SPAWN

      SUBROUTINE MPI_GROUP_FREE(group,ierror)
           INTEGER group
           INTEGER ierror
      END SUBROUTINE MPI_GROUP_FREE

      SUBROUTINE MPI_COMM_SET_ERRHANDLER(comm,errhandler,ierror)
           INTEGER comm
           INTEGER errhandler
           INTEGER ierror
      END SUBROUTINE MPI_COMM_SET_ERRHANDLER

      SUBROUTINE MPI_WIN_TEST(win,flag,ierror)
           INTEGER win
           LOGICAL flag
           INTEGER ierror
      END SUBROUTINE MPI_WIN_TEST

      SUBROUTINE MPI_WIN_FLUSH_LOCAL_ALL(win,ierror)
           INTEGER win
           INTEGER ierror
      END SUBROUTINE MPI_WIN_FLUSH_LOCAL_ALL

      SUBROUTINE MPI_GRAPH_MAP(comm,nnodes,indx,edges,newrank,ierror)
           INTEGER comm
           INTEGER nnodes
           INTEGER indx(*)
           INTEGER edges(*)
           INTEGER newrank
           INTEGER ierror
      END SUBROUTINE MPI_GRAPH_MAP

      SUBROUTINE MPI_PUBLISH_NAME(service_name,info,port_name,ierror)
           CHARACTER (LEN=*) service_name
           INTEGER info
           CHARACTER (LEN=*) port_name
           INTEGER ierror
      END SUBROUTINE MPI_PUBLISH_NAME

      SUBROUTINE MPI_TYPE_GET_EXTENT_X(datatype,lb,extent,ierror)
           USE MPI_CONSTANTS,ONLY:MPI_COUNT_KIND
           INTEGER datatype
           INTEGER(KIND=MPI_COUNT_KIND) lb
           INTEGER(KIND=MPI_COUNT_KIND) extent
           INTEGER ierror
      END SUBROUTINE MPI_TYPE_GET_EXTENT_X

      SUBROUTINE MPI_TYPE_CREATE_F90_REAL(precision,range,newtype,ierror)
           INTEGER precision
           INTEGER range
           INTEGER newtype
           INTEGER ierror
      END SUBROUTINE MPI_TYPE_CREATE_F90_REAL

      SUBROUTINE MPI_GROUP_UNION(group1,group2,newgroup,ierror)
           INTEGER group1
           INTEGER group2
           INTEGER newgroup
           INTEGER ierror
      END SUBROUTINE MPI_GROUP_UNION

      SUBROUTINE MPI_COMM_ACCEPT(port_name,info,root,comm,newcomm,ierror)
           CHARACTER (LEN=*) port_name
           INTEGER info
           INTEGER root
           INTEGER comm
           INTEGER newcomm
           INTEGER ierror
      END SUBROUTINE MPI_COMM_ACCEPT

      SUBROUTINE MPI_WIN_GET_NAME(win,win_name,resultlen,ierror)
           INTEGER win
           CHARACTER (LEN=*) win_name
           INTEGER resultlen
           INTEGER ierror
      END SUBROUTINE MPI_WIN_GET_NAME

      SUBROUTINE MPI_INFO_CREATE(info,ierror)
           INTEGER info
           INTEGER ierror
      END SUBROUTINE MPI_INFO_CREATE

      SUBROUTINE MPI_TYPE_CREATE_F90_INTEGER(range,newtype,ierror)
           INTEGER range
           INTEGER newtype
           INTEGER ierror
      END SUBROUTINE MPI_TYPE_CREATE_F90_INTEGER

      SUBROUTINE MPI_TYPE_SET_NAME(datatype,type_name,ierror)
           INTEGER datatype
           CHARACTER (LEN=*) type_name
           INTEGER ierror
      END SUBROUTINE MPI_TYPE_SET_NAME

      SUBROUTINE MPI_GROUP_INCL(group,n,ranks,newgroup,ierror)
           INTEGER group
           INTEGER n
           INTEGER ranks(*)
           INTEGER newgroup
           INTEGER ierror
      END SUBROUTINE MPI_GROUP_INCL

      SUBROUTINE MPI_COMM_CONNECT(port_name,info,root,comm,newcomm,ierror)
           CHARACTER (LEN=*) port_name
           INTEGER info
           INTEGER root
           INTEGER comm
           INTEGER newcomm
           INTEGER ierror
      END SUBROUTINE MPI_COMM_CONNECT

      SUBROUTINE MPI_COMM_CREATE_ERRHANDLER(comm_errhandler_fn,errhandler,&
                  ierror)
           INTERFACE 
       SUBROUTINE comm_errhandler_fn(vv0,vv1)
       INTEGER vv0,vv1
       END SUBROUTINE
       END INTERFACE
           INTEGER errhandler
           INTEGER ierror
      END SUBROUTINE MPI_COMM_CREATE_ERRHANDLER

      SUBROUTINE MPI_ERROR_STRING(errorcode,string,resultlen,ierror)
           INTEGER errorcode
           CHARACTER (LEN=*) string
           INTEGER resultlen
           INTEGER ierror
      END SUBROUTINE MPI_ERROR_STRING

      SUBROUTINE MPI_TYPE_STRUCT(count,array_of_blocklengths,&
                  array_of_displacements,array_of_types,newtype,ierror)
           INTEGER count
           INTEGER array_of_blocklengths(*)
           INTEGER array_of_displacements(*)
           INTEGER array_of_types(*)
           INTEGER newtype
           INTEGER ierror
      END SUBROUTINE MPI_TYPE_STRUCT

      SUBROUTINE MPI_TYPE_CREATE_INDEXED_BLOCK(count,blocklength,&
                  array_of_displacements,oldtype,newtype,ierror)
           INTEGER count
           INTEGER blocklength
           INTEGER array_of_displacements(*)
           INTEGER oldtype
           INTEGER newtype
           INTEGER ierror
      END SUBROUTINE MPI_TYPE_CREATE_INDEXED_BLOCK

      SUBROUTINE MPI_TYPE_CREATE_HVECTOR(count,blocklength,stride,oldtype,&
                  newtype,ierror)
           USE MPI_CONSTANTS,ONLY:MPI_ADDRESS_KIND
           INTEGER count
           INTEGER blocklength
           INTEGER(KIND=MPI_ADDRESS_KIND) stride
           INTEGER oldtype
           INTEGER newtype
           INTEGER ierror
      END SUBROUTINE MPI_TYPE_CREATE_HVECTOR

      SUBROUTINE MPI_ABORT(comm,errorcode,ierror)
           INTEGER comm
           INTEGER errorcode
           INTEGER ierror
      END SUBROUTINE MPI_ABORT

      SUBROUTINE MPI_TYPE_FREE_KEYVAL(type_keyval,ierror)
           INTEGER type_keyval
           INTEGER ierror
      END SUBROUTINE MPI_TYPE_FREE_KEYVAL

      SUBROUTINE MPI_INTERCOMM_CREATE(local_comm,local_leader,peer_comm,&
                  remote_leader,tag,newintercomm,ierror)
           INTEGER local_comm
           INTEGER local_leader
           INTEGER peer_comm
           INTEGER remote_leader
           INTEGER tag
           INTEGER newintercomm
           INTEGER ierror
      END SUBROUTINE MPI_INTERCOMM_CREATE

      SUBROUTINE MPI_COMM_RANK(comm,rank,ierror)
           INTEGER comm
           INTEGER rank
           INTEGER ierror
      END SUBROUTINE MPI_COMM_RANK

      SUBROUTINE MPI_COMM_IDUP(comm,newcomm,request,ierror)
           INTEGER comm
           INTEGER newcomm
           INTEGER request
           INTEGER ierror
      END SUBROUTINE MPI_COMM_IDUP

      SUBROUTINE MPI_STATUS_SET_CANCELLED(status,flag,ierror)
           USE MPI_CONSTANTS,ONLY:MPI_STATUS_SIZE
           INTEGER status(MPI_STATUS_SIZE)
           LOGICAL flag
           INTEGER ierror
      END SUBROUTINE MPI_STATUS_SET_CANCELLED

      SUBROUTINE MPI_FILE_SET_ERRHANDLER(file,errhandler,ierror)
           INTEGER file
           INTEGER errhandler
           INTEGER ierror
      END SUBROUTINE MPI_FILE_SET_ERRHANDLER

      SUBROUTINE MPI_UNPUBLISH_NAME(service_name,info,port_name,ierror)
           CHARACTER (LEN=*) service_name
           INTEGER info
           CHARACTER (LEN=*) port_name
           INTEGER ierror
      END SUBROUTINE MPI_UNPUBLISH_NAME

      SUBROUTINE MPI_INFO_DELETE(info,key,ierror)
           INTEGER info
           CHARACTER (LEN=*) key
           INTEGER ierror
      END SUBROUTINE MPI_INFO_DELETE

      SUBROUTINE MPI_TYPE_CREATE_RESIZED(oldtype,lb,extent,newtype,ierror)
           USE MPI_CONSTANTS,ONLY:MPI_ADDRESS_KIND
           INTEGER oldtype
           INTEGER(KIND=MPI_ADDRESS_KIND) lb
           INTEGER(KIND=MPI_ADDRESS_KIND) extent
           INTEGER newtype
           INTEGER ierror
      END SUBROUTINE MPI_TYPE_CREATE_RESIZED

      SUBROUTINE MPI_ERRHANDLER_SET(comm,errhandler,ierror)
           INTEGER comm
           INTEGER errhandler
           INTEGER ierror
      END SUBROUTINE MPI_ERRHANDLER_SET

      SUBROUTINE MPI_TYPE_DUP(oldtype,newtype,ierror)
           INTEGER oldtype
           INTEGER newtype
           INTEGER ierror
      END SUBROUTINE MPI_TYPE_DUP

      SUBROUTINE MPI_WIN_FLUSH(rank,win,ierror)
           INTEGER rank
           INTEGER win
           INTEGER ierror
      END SUBROUTINE MPI_WIN_FLUSH

      SUBROUTINE MPI_WAIT(request,status,ierror)
           USE MPI_CONSTANTS,ONLY:MPI_STATUS_SIZE
           INTEGER request
           INTEGER status(MPI_STATUS_SIZE)
           INTEGER ierror
      END SUBROUTINE MPI_WAIT

      SUBROUTINE MPI_INFO_SET(info,key,value,ierror)
           INTEGER info
           CHARACTER (LEN=*) key
           CHARACTER (LEN=*) value
           INTEGER ierror
      END SUBROUTINE MPI_INFO_SET

      SUBROUTINE MPI_TEST(request,flag,status,ierror)
           USE MPI_CONSTANTS,ONLY:MPI_STATUS_SIZE
           INTEGER request
           LOGICAL flag
           INTEGER status(MPI_STATUS_SIZE)
           INTEGER ierror
      END SUBROUTINE MPI_TEST

      SUBROUTINE MPI_COMM_GET_NAME(comm,comm_name,resultlen,ierror)
           INTEGER comm
           CHARACTER (LEN=*) comm_name
           INTEGER resultlen
           INTEGER ierror
      END SUBROUTINE MPI_COMM_GET_NAME

      SUBROUTINE MPI_COMM_FREE(comm,ierror)
           INTEGER comm
           INTEGER ierror
      END SUBROUTINE MPI_COMM_FREE

      SUBROUTINE MPI_COMM_GET_INFO(comm,info,ierror)
           INTEGER comm
           INTEGER info
           INTEGER ierror
      END SUBROUTINE MPI_COMM_GET_INFO

      SUBROUTINE MPI_IPROBE(source,tag,comm,flag,status,ierror)
           USE MPI_CONSTANTS,ONLY:MPI_STATUS_SIZE
           INTEGER source
           INTEGER tag
           INTEGER comm
           LOGICAL flag
           INTEGER status(MPI_STATUS_SIZE)
           INTEGER ierror
      END SUBROUTINE MPI_IPROBE

      SUBROUTINE MPI_GET_PROCESSOR_NAME(name,resultlen,ierror)
           CHARACTER (LEN=*) name
           INTEGER resultlen
           INTEGER ierror
      END SUBROUTINE MPI_GET_PROCESSOR_NAME

      SUBROUTINE MPI_TOPO_TEST(comm,status,ierror)
           INTEGER comm
           INTEGER status
           INTEGER ierror
      END SUBROUTINE MPI_TOPO_TEST

      SUBROUTINE MPI_OP_FREE(op,ierror)
           INTEGER op
           INTEGER ierror
      END SUBROUTINE MPI_OP_FREE

      SUBROUTINE MPI_COMM_SIZE(comm,size,ierror)
           INTEGER comm
           INTEGER size
           INTEGER ierror
      END SUBROUTINE MPI_COMM_SIZE

      SUBROUTINE MPI_WIN_UNLOCK(rank,win,ierror)
           INTEGER rank
           INTEGER win
           INTEGER ierror
      END SUBROUTINE MPI_WIN_UNLOCK

      SUBROUTINE MPI_ERRHANDLER_FREE(errhandler,ierror)
           INTEGER errhandler
           INTEGER ierror
      END SUBROUTINE MPI_ERRHANDLER_FREE

      SUBROUTINE MPI_WIN_UNLOCK_ALL(win,ierror)
           INTEGER win
           INTEGER ierror
      END SUBROUTINE MPI_WIN_UNLOCK_ALL

      SUBROUTINE MPI_TYPE_HINDEXED(count,array_of_blocklengths,&
                  array_of_displacements,oldtype,newtype,ierror)
           INTEGER count
           INTEGER array_of_blocklengths(*)
           INTEGER array_of_displacements(*)
           INTEGER oldtype
           INTEGER newtype
           INTEGER ierror
      END SUBROUTINE MPI_TYPE_HINDEXED

      SUBROUTINE MPI_MPROBE(source,tag,comm,message,status,ierror)
           USE MPI_CONSTANTS,ONLY:MPI_STATUS_SIZE
           INTEGER source
           INTEGER tag
           INTEGER comm
           INTEGER message
           INTEGER status(MPI_STATUS_SIZE)
           INTEGER ierror
      END SUBROUTINE MPI_MPROBE

      SUBROUTINE MPI_WIN_SET_NAME(win,win_name,ierror)
           INTEGER win
           CHARACTER (LEN=*) win_name
           INTEGER ierror
      END SUBROUTINE MPI_WIN_SET_NAME

      SUBROUTINE MPI_TYPE_SIZE(datatype,size,ierror)
           INTEGER datatype
           INTEGER size
           INTEGER ierror
      END SUBROUTINE MPI_TYPE_SIZE

      SUBROUTINE MPI_WIN_START(group,assert,win,ierror)
           INTEGER group
           INTEGER assert
           INTEGER win
           INTEGER ierror
      END SUBROUTINE MPI_WIN_START

      SUBROUTINE MPI_WIN_CREATE_ERRHANDLER(win_errhandler_fn,errhandler,ierror)
           INTERFACE 
       SUBROUTINE win_errhandler_fn(vv0,vv1)
       INTEGER vv0,vv1
       END SUBROUTINE
       END INTERFACE
           INTEGER errhandler
           INTEGER ierror
      END SUBROUTINE MPI_WIN_CREATE_ERRHANDLER

      SUBROUTINE MPI_WIN_FENCE(assert,win,ierror)
           INTEGER assert
           INTEGER win
           INTEGER ierror
      END SUBROUTINE MPI_WIN_FENCE

      SUBROUTINE MPI_TYPE_GET_ENVELOPE(datatype,num_integers,num_addresses,&
                  num_datatypes,combiner,ierror)
           INTEGER datatype
           INTEGER num_integers
           INTEGER num_addresses
           INTEGER num_datatypes
           INTEGER combiner
           INTEGER ierror
      END SUBROUTINE MPI_TYPE_GET_ENVELOPE

      SUBROUTINE MPI_GREQUEST_COMPLETE(request,ierror)
           INTEGER request
           INTEGER ierror
      END SUBROUTINE MPI_GREQUEST_COMPLETE

      SUBROUTINE MPI_GET_VERSION(version,subversion,ierror)
           INTEGER version
           INTEGER subversion
           INTEGER ierror
      END SUBROUTINE MPI_GET_VERSION

      SUBROUTINE MPI_TYPE_HVECTOR(count,blocklength,stride,oldtype,newtype,&
                  ierror)
           INTEGER count
           INTEGER blocklength
           INTEGER stride
           INTEGER oldtype
           INTEGER newtype
           INTEGER ierror
      END SUBROUTINE MPI_TYPE_HVECTOR

      SUBROUTINE MPI_WIN_SYNC(win,ierror)
           INTEGER win
           INTEGER ierror
      END SUBROUTINE MPI_WIN_SYNC

      SUBROUTINE MPI_GET_ELEMENTS_X(status,datatype,count,ierror)
           USE MPI_CONSTANTS,ONLY:MPI_COUNT_KIND, MPI_STATUS_SIZE
           INTEGER status(MPI_STATUS_SIZE)
           INTEGER datatype
           INTEGER(KIND=MPI_COUNT_KIND) count
           INTEGER ierror
      END SUBROUTINE MPI_GET_ELEMENTS_X

      SUBROUTINE MPI_COMM_CALL_ERRHANDLER(comm,errorcode,ierror)
           INTEGER comm
           INTEGER errorcode
           INTEGER ierror
      END SUBROUTINE MPI_COMM_CALL_ERRHANDLER

      SUBROUTINE MPI_COMM_TEST_INTER(comm,flag,ierror)
           INTEGER comm
           LOGICAL flag
           INTEGER ierror
      END SUBROUTINE MPI_COMM_TEST_INTER

      SUBROUTINE MPI_COMM_JOIN(fd,intercomm,ierror)
           INTEGER fd
           INTEGER intercomm
           INTEGER ierror
      END SUBROUTINE MPI_COMM_JOIN

      SUBROUTINE MPI_CLOSE_PORT(port_name,ierror)
           CHARACTER (LEN=*) port_name
           INTEGER ierror
      END SUBROUTINE MPI_CLOSE_PORT

      SUBROUTINE MPI_TYPE_CREATE_F90_COMPLEX(precision,range,newtype,ierror)
           INTEGER precision
           INTEGER range
           INTEGER newtype
           INTEGER ierror
      END SUBROUTINE MPI_TYPE_CREATE_F90_COMPLEX

      SUBROUTINE MPI_GROUP_RANK(group,rank,ierror)
           INTEGER group
           INTEGER rank
           INTEGER ierror
      END SUBROUTINE MPI_GROUP_RANK


        SUBROUTINE MPI_INIT(ierror)
        INTEGER ierror
        END SUBROUTINE MPI_INIT

        SUBROUTINE MPI_INIT_THREAD(v0,v1,ierror)
        INTEGER v0, v1, ierror
        END SUBROUTINE MPI_INIT_THREAD

        FUNCTION MPI_WTIME()
            REAL*8 MPI_WTIME
        END FUNCTION MPI_WTIME
!
        FUNCTION MPI_WTICK()
            REAL*8 MPI_WTICK
        END FUNCTION MPI_WTICK

! style:PMPIuse:PMPI_WTIME:3 sig:0
        FUNCTION PMPI_WTIME()
            REAL*8 PMPI_WTIME
        END FUNCTION PMPI_WTIME
!
! style:PMPIuse:PMPI_WTICK:3 sig:0
        FUNCTION PMPI_WTICK()
            REAL*8 PMPI_WTICK
        END FUNCTION PMPI_WTICK

        SUBROUTINE MPI_NULL_DELETE_FN(COMM, KEYVAL, ATTRIBUTE_VAL,&
          EXTRA_STATE, IERROR)
            USE MPI_CONSTANTS,ONLY: MPI_ADDRESS_KIND
            INTEGER COMM, KEYVAL, IERROR
            INTEGER(KIND=MPI_ADDRESS_KIND) ATTRIBUTE_VAL, EXTRA_STATE
        END SUBROUTINE MPI_NULL_DELETE_FN

        SUBROUTINE MPI_DUP_FN(OLDCOMM, KEYVAL, EXTRA_STATE,&
          ATTRIBUTE_VAL_IN, ATTRIBUTE_VAL_OUT, FLAG, IERROR)
            USE MPI_CONSTANTS,ONLY: MPI_ADDRESS_KIND
            INTEGER OLDCOMM, KEYVAL, IERROR
            INTEGER(KIND=MPI_ADDRESS_KIND) EXTRA_STATE, ATTRIBUTE_VAL_IN, ATTRIBUTE_VAL_OUT
            LOGICAL FLAG
        END SUBROUTINE MPI_DUP_FN

        SUBROUTINE MPI_NULL_COPY_FN(OLDCOMM, KEYVAL, EXTRA_STATE,&
          ATTRIBUTE_VAL_IN, ATTRIBUTE_VAL_OUT, FLAG, IERROR)
            USE MPI_CONSTANTS,ONLY: MPI_ADDRESS_KIND
            INTEGER OLDCOMM, KEYVAL, IERROR
            INTEGER(KIND=MPI_ADDRESS_KIND) EXTRA_STATE, ATTRIBUTE_VAL_IN, ATTRIBUTE_VAL_OUT
            LOGICAL FLAG
        END SUBROUTINE MPI_NULL_COPY_FN

        SUBROUTINE MPI_COMM_NULL_DELETE_FN(COMM, COMM_KEYVAL, ATTRIBUTE_VAL,&
          EXTRA_STATE, IERROR)
            USE MPI_CONSTANTS,ONLY: MPI_ADDRESS_KIND
            INTEGER COMM, COMM_KEYVAL, IERROR
            INTEGER(KIND=MPI_ADDRESS_KIND) ATTRIBUTE_VAL, EXTRA_STATE
        END SUBROUTINE MPI_COMM_NULL_DELETE_FN

        SUBROUTINE MPI_COMM_DUP_FN(OLDCOMM, COMM_KEYVAL, EXTRA_STATE,&
          ATTRIBUTE_VAL_IN, ATTRIBUTE_VAL_OUT, FLAG, IERROR)
            USE MPI_CONSTANTS,ONLY: MPI_ADDRESS_KIND
            INTEGER OLDCOMM, COMM_KEYVAL, IERROR
            INTEGER(KIND=MPI_ADDRESS_KIND) EXTRA_STATE, ATTRIBUTE_VAL_IN, ATTRIBUTE_VAL_OUT
            LOGICAL FLAG
        END SUBROUTINE MPI_COMM_DUP_FN

        SUBROUTINE MPI_COMM_NULL_COPY_FN(OLDCOMM, COMM_KEYVAL, EXTRA_STATE,&
          ATTRIBUTE_VAL_IN, ATTRIBUTE_VAL_OUT, FLAG, IERROR)
            USE MPI_CONSTANTS,ONLY: MPI_ADDRESS_KIND
            INTEGER OLDCOMM, COMM_KEYVAL, IERROR
            INTEGER(KIND=MPI_ADDRESS_KIND) EXTRA_STATE, ATTRIBUTE_VAL_IN, ATTRIBUTE_VAL_OUT
            LOGICAL FLAG
        END SUBROUTINE MPI_COMM_NULL_COPY_FN

        SUBROUTINE MPI_TYPE_NULL_DELETE_FN(DATATYPE, TYPE_KEYVAL, ATTRIBUTE_VAL,&
          EXTRA_STATE, IERROR)
            USE MPI_CONSTANTS,ONLY: MPI_ADDRESS_KIND
            INTEGER DATATYPE, TYPE_KEYVAL, IERROR
            INTEGER(KIND=MPI_ADDRESS_KIND) ATTRIBUTE_VAL, EXTRA_STATE
        END SUBROUTINE MPI_TYPE_NULL_DELETE_FN

        SUBROUTINE MPI_TYPE_DUP_FN(OLDTYPE, TYPE_KEYVAL, EXTRA_STATE,&
          ATTRIBUTE_VAL_IN, ATTRIBUTE_VAL_OUT, FLAG, IERROR)
            USE MPI_CONSTANTS,ONLY: MPI_ADDRESS_KIND
            INTEGER OLDTYPE, TYPE_KEYVAL, IERROR
            INTEGER(KIND=MPI_ADDRESS_KIND) EXTRA_STATE, ATTRIBUTE_VAL_IN, ATTRIBUTE_VAL_OUT
            LOGICAL FLAG
        END SUBROUTINE MPI_TYPE_DUP_FN

        SUBROUTINE MPI_TYPE_NULL_COPY_FN(OLDTYPE, TYPE_KEYVAL, EXTRA_STATE,&
          ATTRIBUTE_VAL_IN, ATTRIBUTE_VAL_OUT, FLAG, IERROR)
            USE MPI_CONSTANTS,ONLY: MPI_ADDRESS_KIND
            INTEGER OLDTYPE, TYPE_KEYVAL, IERROR
            INTEGER(KIND=MPI_ADDRESS_KIND) EXTRA_STATE, ATTRIBUTE_VAL_IN, ATTRIBUTE_VAL_OUT
            LOGICAL FLAG
        END SUBROUTINE MPI_TYPE_NULL_COPY_FN

        SUBROUTINE MPI_WIN_NULL_DELETE_FN(WIN, WIN_KEYVAL, ATTRIBUTE_VAL,&
          EXTRA_STATE, IERROR)
            USE MPI_CONSTANTS,ONLY: MPI_ADDRESS_KIND
            INTEGER WIN, WIN_KEYVAL, IERROR
            INTEGER(KIND=MPI_ADDRESS_KIND) ATTRIBUTE_VAL, EXTRA_STATE
        END SUBROUTINE MPI_WIN_NULL_DELETE_FN

        SUBROUTINE MPI_WIN_DUP_FN(OLDWIN, WIN_KEYVAL, EXTRA_STATE,&
          ATTRIBUTE_VAL_IN, ATTRIBUTE_VAL_OUT, FLAG, IERROR)
            USE MPI_CONSTANTS,ONLY: MPI_ADDRESS_KIND
            INTEGER OLDWIN, WIN_KEYVAL, IERROR
            INTEGER(KIND=MPI_ADDRESS_KIND) EXTRA_STATE, ATTRIBUTE_VAL_IN, ATTRIBUTE_VAL_OUT
            LOGICAL FLAG
        END SUBROUTINE MPI_WIN_DUP_FN

        SUBROUTINE MPI_WIN_NULL_COPY_FN(OLDWIN, WIN_KEYVAL, EXTRA_STATE,&
          ATTRIBUTE_VAL_IN, ATTRIBUTE_VAL_OUT, FLAG, IERROR)
            USE MPI_CONSTANTS,ONLY: MPI_ADDRESS_KIND
            INTEGER OLDWIN, WIN_KEYVAL, IERROR
            INTEGER(KIND=MPI_ADDRESS_KIND) EXTRA_STATE, ATTRIBUTE_VAL_IN, ATTRIBUTE_VAL_OUT
            LOGICAL FLAG
        END SUBROUTINE MPI_WIN_NULL_COPY_FN
       END INTERFACE
       END MODULE MPI_BASE
