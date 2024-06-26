/* Copyright (c) 2001-2016, The Ohio State University. All rights
 * reserved.
 *
 * This file is part of the MVAPICH2 software package developed by the
 * team members of The Ohio State University's Network-Based Computing
 * Laboratory (NBCL), headed by Professor Dhabaleswar K. (DK) Panda.
 *
 * For detailed copyright and licensing information, please refer to the
 * copyright file COPYRIGHT in the top level MVAPICH2 directory.
 */

#define GEN2__INTEL_XEON_E5_2680_24__MLX_CX_FDR__24PPN {\
	{		\
	48,		\
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},		\
	15,		\
	{		\
	{1, &MPIR_Alltoall_RD_MV2},		\
	{2, &MPIR_Alltoall_RD_MV2},		\
	{4, &MPIR_Alltoall_bruck_MV2},		\
	{8, &MPIR_Alltoall_bruck_MV2},		\
	{16, &MPIR_Alltoall_bruck_MV2},		\
	{32, &MPIR_Alltoall_bruck_MV2},		\
	{64, &MPIR_Alltoall_bruck_MV2},		\
	{128, &MPIR_Alltoall_bruck_MV2},		\
	{256, &MPIR_Alltoall_bruck_MV2},		\
	{512, &MPIR_Alltoall_bruck_MV2},		\
	{1024, &MPIR_Alltoall_bruck_MV2},		\
	{2048, &MPIR_Alltoall_bruck_MV2},		\
	{4096, &MPIR_Alltoall_bruck_MV2},		\
	{8192, &MPIR_Alltoall_inplace_MV2},		\
	{16384, &MPIR_Alltoall_inplace_MV2},		\
	{32768, &MPIR_Alltoall_inplace_MV2}		\
	}		\
	},		\
	{		\
	96,		\
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0},		\
	15,		\
	{		\
	{1, &MPIR_Alltoall_RD_MV2},		\
	{2, &MPIR_Alltoall_bruck_MV2},		\
	{4, &MPIR_Alltoall_bruck_MV2},		\
	{8, &MPIR_Alltoall_bruck_MV2},		\
	{16, &MPIR_Alltoall_bruck_MV2},		\
	{32, &MPIR_Alltoall_bruck_MV2},		\
	{64, &MPIR_Alltoall_bruck_MV2},		\
	{128, &MPIR_Alltoall_bruck_MV2},		\
	{256, &MPIR_Alltoall_bruck_MV2},		\
	{512, &MPIR_Alltoall_bruck_MV2},		\
	{1024, &MPIR_Alltoall_bruck_MV2},		\
	{2048, &MPIR_Alltoall_bruck_MV2},		\
	{4096, &MPIR_Alltoall_bruck_MV2},		\
	{8192, &MPIR_Alltoall_Scatter_dest_MV2},		\
	{16384, &MPIR_Alltoall_pairwise_MV2},		\
	{32768, &MPIR_Alltoall_pairwise_MV2}		\
	}		\
	},		\
	{		\
	192,		\
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},		\
	15,		\
	{		\
	{1, &MPIR_Alltoall_bruck_MV2},		\
	{2, &MPIR_Alltoall_bruck_MV2},		\
	{4, &MPIR_Alltoall_bruck_MV2},		\
	{8, &MPIR_Alltoall_bruck_MV2},		\
	{16, &MPIR_Alltoall_bruck_MV2},		\
	{32, &MPIR_Alltoall_bruck_MV2},		\
	{64, &MPIR_Alltoall_bruck_MV2},		\
	{128, &MPIR_Alltoall_bruck_MV2},		\
	{256, &MPIR_Alltoall_bruck_MV2},		\
	{512, &MPIR_Alltoall_bruck_MV2},		\
	{1024, &MPIR_Alltoall_bruck_MV2},		\
	{2048, &MPIR_Alltoall_bruck_MV2},		\
	{4096, &MPIR_Alltoall_pairwise_MV2},		\
	{8192, &MPIR_Alltoall_pairwise_MV2},		\
	{16384, &MPIR_Alltoall_pairwise_MV2},		\
	{32768, &MPIR_Alltoall_pairwise_MV2}		\
	}		\
	},		\
	{		\
	384,		\
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},		\
	15,		\
	{		\
	{1, &MPIR_Alltoall_bruck_MV2},		\
	{2, &MPIR_Alltoall_bruck_MV2},		\
	{4, &MPIR_Alltoall_bruck_MV2},		\
	{8, &MPIR_Alltoall_bruck_MV2},		\
	{16, &MPIR_Alltoall_bruck_MV2},		\
	{32, &MPIR_Alltoall_bruck_MV2},		\
	{64, &MPIR_Alltoall_bruck_MV2},		\
	{128, &MPIR_Alltoall_bruck_MV2},		\
	{256, &MPIR_Alltoall_bruck_MV2},		\
	{512, &MPIR_Alltoall_bruck_MV2},		\
	{1024, &MPIR_Alltoall_bruck_MV2},		\
	{2048, &MPIR_Alltoall_bruck_MV2},		\
	{4096, &MPIR_Alltoall_pairwise_MV2},		\
	{8192, &MPIR_Alltoall_pairwise_MV2},		\
	{16384, &MPIR_Alltoall_pairwise_MV2},		\
	{32768, &MPIR_Alltoall_pairwise_MV2}		\
	}		\
	},		\
	{		\
	768,		\
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0},		\
	15,		\
	{		\
	{1, &MPIR_Alltoall_bruck_MV2},		\
	{2, &MPIR_Alltoall_bruck_MV2},		\
	{4, &MPIR_Alltoall_bruck_MV2},		\
	{8, &MPIR_Alltoall_bruck_MV2},		\
	{16, &MPIR_Alltoall_bruck_MV2},		\
	{32, &MPIR_Alltoall_bruck_MV2},		\
	{64, &MPIR_Alltoall_bruck_MV2},		\
	{128, &MPIR_Alltoall_bruck_MV2},		\
	{256, &MPIR_Alltoall_bruck_MV2},		\
	{512, &MPIR_Alltoall_bruck_MV2},		\
	{1024, &MPIR_Alltoall_bruck_MV2},		\
	{2048, &MPIR_Alltoall_pairwise_MV2},		\
	{4096, &MPIR_Alltoall_pairwise_MV2},		\
	{8192, &MPIR_Alltoall_Scatter_dest_MV2},		\
	{16384, &MPIR_Alltoall_Scatter_dest_MV2},		\
	{32768, &MPIR_Alltoall_pairwise_MV2}		\
	}		\
	}		\
}
