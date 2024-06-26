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
#include "mpichconf.h"
#include "mpidi_ch3_impl.h"
#include <mpimem.h>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <netdb.h>
#include <sys/mman.h>
#include <errno.h>
#include <string.h>
#include <assert.h>
#include "upmi.h"
#include "mpiutil.h"
#include "hwloc_bind.h"
#if defined(HAVE_LIBIBVERBS)
#include <hwloc/openfabrics-verbs.h>
#endif
#if defined(CHANNEL_MRAIL)
#include "smp_smpi.h"
#include "rdma_impl.h"
#endif /*defined(CHANNEL_MRAIL)*/
#include "mv2_arch_hca_detect.h"
#include "debug_utils.h"

/* CPU Mapping related definitions */

#define CONFIG_FILE "/proc/cpuinfo"
#define MAX_LINE_LENGTH 512
#define MAX_NAME_LENGTH 64

int mv2_my_cpu_id = -1;
int mv2_my_sock_id = -1;
int mv2_my_async_cpu_id = -1;
int *local_core_ids = NULL;
int mv2_user_defined_mapping = FALSE;


unsigned int mv2_enable_affinity = 1;
unsigned int mv2_enable_leastload = 0;
unsigned int mv2_hca_aware_process_mapping = 1;

typedef enum {
    CPU_FAMILY_NONE = 0,
    CPU_FAMILY_INTEL,
    CPU_FAMILY_AMD,
} cpu_type_t;

int CLOVERTOWN_MODEL = 15;
int HARPERTOWN_MODEL = 23;
int NEHALEM_MODEL = 26;

int ip = 0;
unsigned long *core_mapping = NULL;
int *obj_tree = NULL;

policy_type_t mv2_binding_policy;
level_type_t mv2_binding_level;
hwloc_topology_t topology;

static int INTEL_XEON_DUAL_MAPPING[] = { 0, 1, 0, 1 };

/* ((0,1),(4,5))((2,3),(6,7)) */
static int INTEL_CLOVERTOWN_MAPPING[] = { 0, 0, 1, 1, 0, 0, 1, 1 };

/* legacy ((0,2),(4,6))((1,3),(5,7)) */
static int INTEL_HARPERTOWN_LEG_MAPPING[] = { 0, 1, 0, 1, 0, 1, 0, 1 };

/* common ((0,1),(2,3))((4,5),(6,7)) */
static int INTEL_HARPERTOWN_COM_MAPPING[] = { 0, 0, 0, 0, 1, 1, 1, 1 };

/* legacy (0,2,4,6)(1,3,5,7) with hyperthreading */
static int INTEL_NEHALEM_LEG_MAPPING[] =
    { 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1 };

/* common (0,1,2,3)(4,5,6,7) with hyperthreading */
static int INTEL_NEHALEM_COM_MAPPING[] =
    { 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1 };

static int AMD_OPTERON_DUAL_MAPPING[] = { 0, 0, 1, 1 };
static int AMD_BARCELONA_MAPPING[] = { 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3 };

extern int use_hwloc_cpu_binding;

char *s_cpu_mapping = NULL;
static char *custom_cpu_mapping = NULL;
int s_cpu_mapping_line_max = _POSIX2_LINE_MAX;
static int custom_cpu_mapping_line_max = _POSIX2_LINE_MAX;
char *cpu_mapping = NULL;

int ib_socket_bind = 0;

#if defined(CHANNEL_MRAIL)
int get_ib_socket(struct ibv_device * ibdev)
{
    hwloc_cpuset_t set = hwloc_bitmap_alloc();
    hwloc_obj_t osdev = NULL;
    char string[256];
    int socket_id = 0;

    if (NULL == set) {
        goto free_and_return;
    }

    if (hwloc_ibv_get_device_cpuset(topology, ibdev, set)) {
        goto free_and_return;
    }

    osdev = hwloc_get_obj_inside_cpuset_by_type(topology, set,
            HWLOC_OBJ_SOCKET, 0);

    if (NULL == osdev) {
        goto free_and_return;
    }

    /*
     * The hwloc object "string" will have the form "Socket#n" so we are
     * looking at the 8th char to detect which socket is.
     */
    hwloc_obj_type_snprintf(string, sizeof(string), osdev, 1);
    return osdev->os_index;

free_and_return:
    hwloc_bitmap_free(set);
    return socket_id;
}
#endif /* defined(CHANNEL_MRAIL) */

static int first_num_from_str(char **str)
{
    int val = atoi(*str);
    while (isdigit(**str)) {
        (*str)++;
    }
    return val;
}

static inline int compare_float(const float a, const float b)
{
    const float precision = 0.00001;
    if ((a - precision) < b && (a + precision) > b) {
        return 1;
    } else {
        return 0;
    }
}

static int pid_filter(const struct dirent *dir_obj)
{
    int i;
    int length = strlen(dir_obj->d_name);

    for (i = 0; i < length; i++) {
        if (!isdigit(dir_obj->d_name[i])) {
            return 0;
        }
    }
    return 1;
}

static void find_parent(hwloc_obj_t obj, hwloc_obj_type_t type, hwloc_obj_t * parent)
{
    if ((type == HWLOC_OBJ_CORE) || (type == HWLOC_OBJ_SOCKET)
        || (type == HWLOC_OBJ_NODE)) {
        if (obj->parent->type == type) {
            *parent = obj->parent;
            return;
        } else {
            find_parent(obj->parent, type, parent);
        }
    } else {
        return;
    }
}

static void find_leastload_node(obj_attribute_type * tree, hwloc_obj_t original,
                                hwloc_obj_t * result)
{
    int i, j, k, per, ix, depth_nodes, num_nodes, depth_sockets, num_sockets;
    hwloc_obj_t obj, tmp;

    depth_nodes = hwloc_get_type_depth(topology, HWLOC_OBJ_NODE);
    num_nodes = hwloc_get_nbobjs_by_depth(topology, depth_nodes);

    /* One socket includes multi numanodes. */
    if ((original->type == HWLOC_OBJ_SOCKET)) {
        depth_sockets = hwloc_get_type_depth(topology, HWLOC_OBJ_SOCKET);
        num_sockets = hwloc_get_nbobjs_by_depth(topology, depth_sockets);
        per = num_nodes / num_sockets;
        ix = (original->logical_index) * per;
        if (per == 1) {
            *result = tree[depth_nodes * num_nodes + ix].obj;
        } else {
            i = depth_nodes * num_nodes + ix;
            for (k = 0; k < (per - 1); k++) {
                j = i + k + 1;
                i = (tree[i].load > tree[j].load) ? j : i;
            }
            *result = tree[i].obj;
        }
    } else if (original->type == HWLOC_OBJ_MACHINE) {
        tmp = NULL;
        for (k = 0; k < num_nodes; k++) {
            obj = hwloc_get_obj_by_depth(topology, depth_nodes, k);
            if (tmp == NULL) {
                tmp = obj;
            } else {
                i = depth_nodes * num_nodes + tmp->logical_index;
                j = depth_nodes * num_nodes + obj->logical_index;
                if (tree[i].load > tree[j].load)
                    tmp = obj;
            }
        }
        *result = tmp;
    } else {
        *result = NULL;
    }
    return;
}

static void find_leastload_socket(obj_attribute_type * tree, hwloc_obj_t original,
                                  hwloc_obj_t * result)
{
    int i, j, k, per, ix, depth_sockets, num_sockets, depth_nodes, num_nodes;
    hwloc_obj_t obj, tmp;

    depth_sockets = hwloc_get_type_depth(topology, HWLOC_OBJ_SOCKET);
    num_sockets = hwloc_get_nbobjs_by_depth(topology, depth_sockets);

    /* One numanode includes multi sockets. */
    if ((original->type == HWLOC_OBJ_NODE)) {
        depth_nodes = hwloc_get_type_depth(topology, HWLOC_OBJ_NODE);
        num_nodes = hwloc_get_nbobjs_by_depth(topology, depth_nodes);
        per = num_sockets / num_nodes;
        ix = (original->logical_index) * per;
        if (per == 1) {
            *result = tree[depth_sockets * num_sockets + ix].obj;
        } else {
            i = depth_sockets * num_sockets + ix;
            for (k = 0; k < (per - 1); k++) {
                j = i + k + 1;
                i = (tree[i].load > tree[j].load) ? j : i;
            }
            *result = tree[i].obj;
        }
    } else if (original->type == HWLOC_OBJ_MACHINE) {
        tmp = NULL;
        for (k = 0; k < num_sockets; k++) {
            obj = hwloc_get_obj_by_depth(topology, depth_sockets, k);
            if (tmp == NULL) {
                tmp = obj;
            } else {
                i = depth_sockets * num_sockets + tmp->logical_index;
                j = depth_sockets * num_sockets + obj->logical_index;
                if (tree[i].load > tree[j].load)
                    tmp = obj;
            }
        }
        *result = tmp;
    } else {
        *result = NULL;
    }
    return;
}

static void find_leastload_core(obj_attribute_type * tree, hwloc_obj_t original,
                                hwloc_obj_t * result)
{
    int i, j, k, per, ix;
    int depth_cores, num_cores, depth_sockets, num_sockets, depth_nodes, num_nodes;

    depth_cores = hwloc_get_type_depth(topology, HWLOC_OBJ_CORE);
    num_cores = hwloc_get_nbobjs_by_depth(topology, depth_cores);

    /* Core may have Socket or Numanode as direct parent. */
    if ((original->type == HWLOC_OBJ_NODE)) {
        depth_nodes = hwloc_get_type_depth(topology, HWLOC_OBJ_NODE);
        num_nodes = hwloc_get_nbobjs_by_depth(topology, depth_nodes);
        per = num_cores / num_nodes;
        ix = (original->logical_index) * per;
        if (per == 1) {
            *result = tree[depth_cores * num_cores + ix].obj;
        } else {
            i = depth_cores * num_cores + ix;
            for (k = 0; k < (per - 1); k++) {
                j = i + k + 1;
                i = (tree[i].load > tree[j].load) ? j : i;
            }
            *result = tree[i].obj;
        }
    } else if (original->type == HWLOC_OBJ_SOCKET) {
        depth_sockets = hwloc_get_type_depth(topology, HWLOC_OBJ_SOCKET);
        num_sockets = hwloc_get_nbobjs_by_depth(topology, depth_sockets);
        per = num_cores / num_sockets;
        ix = (original->logical_index) * per;
        if (per == 1) {
            *result = tree[depth_cores * num_cores + ix].obj;
        } else {
            i = depth_cores * num_cores + ix;
            for (k = 0; k < (per - 1); k++) {
                j = i + k + 1;
                i = (tree[i].load > tree[j].load) ? j : i;
            }
            *result = tree[i].obj;
        }
    } else {
        *result = NULL;
    }
    return;
}

static void find_leastload_pu(obj_attribute_type * tree, hwloc_obj_t original,
                              hwloc_obj_t * result)
{
    int i, j, k, per, ix, depth_pus, num_pus, depth_cores, num_cores;

    depth_pus = hwloc_get_type_depth(topology, HWLOC_OBJ_PU);
    num_pus = hwloc_get_nbobjs_by_depth(topology, depth_pus);

    /* Assume: pu only has core as direct parent. */
    if ((original->type == HWLOC_OBJ_CORE)) {
        depth_cores = hwloc_get_type_depth(topology, HWLOC_OBJ_CORE);
        num_cores = hwloc_get_nbobjs_by_depth(topology, depth_cores);
        per = num_pus / num_cores;
        ix = (original->logical_index) * per;
        if (per == 1) {
            *result = tree[depth_pus * num_pus + ix].obj;
        } else {
            i = depth_pus * num_pus + ix;
            for (k = 0; k < (per - 1); k++) {
                j = i + k + 1;
                i = (tree[i].load > tree[j].load) ? j : i;
            }
            *result = tree[i].obj;
        }
    } else {
        *result = NULL;
    }
    return;
}


static void update_obj_attribute(obj_attribute_type * tree, int ix, hwloc_obj_t obj,
                                 int cpuset, float load)
{
    tree[ix].obj = obj;
    if (!(cpuset < 0)) {
        CPU_SET(cpuset, &(tree[ix].cpuset));
    }
    tree[ix].load += load;
}

static void insert_load(obj_attribute_type * tree, hwloc_obj_t pu, int cpuset, float load)
{
    int k, depth_pus, num_pus = 0;
    int depth_cores, depth_sockets, depth_nodes, num_cores = 0, num_sockets =
        0, num_nodes = 0;
    hwloc_obj_t parent;

    depth_pus = hwloc_get_type_or_below_depth(topology, HWLOC_OBJ_PU);
    num_pus = hwloc_get_nbobjs_by_depth(topology, depth_pus);

    depth_nodes = hwloc_get_type_depth(topology, HWLOC_OBJ_NODE);
    if (depth_nodes != HWLOC_TYPE_DEPTH_UNKNOWN) {
        num_nodes = hwloc_get_nbobjs_by_depth(topology, depth_nodes);
    }
    depth_sockets = hwloc_get_type_depth(topology, HWLOC_OBJ_SOCKET);
    if (depth_sockets != HWLOC_TYPE_DEPTH_UNKNOWN) {
        num_sockets = hwloc_get_nbobjs_by_depth(topology, depth_sockets);
    }
    depth_cores = hwloc_get_type_depth(topology, HWLOC_OBJ_CORE);
    if (depth_cores != HWLOC_TYPE_DEPTH_UNKNOWN) {
        num_cores = hwloc_get_nbobjs_by_depth(topology, depth_cores);
    }

    /* Add obj, cpuset and load for HWLOC_OBJ_PU */
    k = depth_pus * num_pus + pu->logical_index;
    update_obj_attribute(tree, k, pu, cpuset, load);
    /* Add cpuset and load for HWLOC_OBJ_CORE */
    if (depth_cores != HWLOC_TYPE_DEPTH_UNKNOWN) {
        find_parent(pu, HWLOC_OBJ_CORE, &parent);
        k = depth_cores * num_cores + parent->logical_index;
        update_obj_attribute(tree, k, parent, cpuset, load);
    }
    /* Add cpuset and load for HWLOC_OBJ_SOCKET */
    if (depth_sockets != HWLOC_TYPE_DEPTH_UNKNOWN) {
        find_parent(pu, HWLOC_OBJ_SOCKET, &parent);
        k = depth_sockets * num_sockets + parent->logical_index;
        update_obj_attribute(tree, k, parent, cpuset, load);
    }
    /* Add cpuset and load for HWLOC_OBJ_NODE */
    if (depth_nodes != HWLOC_TYPE_DEPTH_UNKNOWN) {
        find_parent(pu, HWLOC_OBJ_NODE, &parent);
        k = depth_nodes * num_nodes + parent->logical_index;
        update_obj_attribute(tree, k, parent, cpuset, load);
    }
    return;
}

static void cac_load(obj_attribute_type * tree, cpu_set_t cpuset)
{
    int i, j, depth_pus, num_pus;
    float proc_load;
    int num_processes = 0;
    hwloc_obj_t obj;

    depth_pus = hwloc_get_type_or_below_depth(topology, HWLOC_OBJ_PU);
    num_pus = hwloc_get_nbobjs_by_depth(topology, depth_pus);

    for (i = 0; i < num_pus; i++) {
        if (CPU_ISSET(i, &cpuset)) {
            num_processes++;
        }
    }

    /* Process is running on num_processes cores; for each core, the load is proc_load. */
    proc_load = 1 / num_processes;

    /*
     * num_objs is HWLOC_OBJ_PU number, and system CPU number;
     * also HWLOC_OBJ_CORE number when HT disabled or without HT.
     */

    for (i = 0; i < num_pus; i++) {
        if (CPU_ISSET(i, &cpuset)) {
            for (j = 0; j < num_pus; j++) {
                obj = hwloc_get_obj_by_depth(topology, depth_pus, j);
                if (obj->os_index == i) {
                    insert_load(tree, obj, i, proc_load);
                }
            }
        }
    }
    return;
}

static void insert_core_mapping(int ix, hwloc_obj_t pu, obj_attribute_type * tree)
{
    core_mapping[ix] = pu->os_index;
    /* This process will be binding to one pu/core.
     * The load for this pu/core is 1; and not update cpuset.
     */
    insert_load(tree, pu, -1, 1);
    return;
}

void map_scatter_load(obj_attribute_type * tree)
{
    int k;
    int depth_cores, depth_sockets, depth_nodes, num_cores = 0;
    hwloc_obj_t root, node, sockets, core_parent, core, result;

    root = hwloc_get_root_obj(topology);

    depth_nodes = hwloc_get_type_depth(topology, HWLOC_OBJ_NODE);

    depth_sockets = hwloc_get_type_depth(topology, HWLOC_OBJ_SOCKET);

    depth_cores = hwloc_get_type_depth(topology, HWLOC_OBJ_CORE);
    if (depth_cores != HWLOC_TYPE_DEPTH_UNKNOWN) {
        num_cores = hwloc_get_nbobjs_by_depth(topology, depth_cores);
    }

    k = 0;
    /*Assume: there is always existing SOCKET, but not always existing NUMANODE(like Clovertown). */
    while (k < num_cores) {
        if (depth_nodes == HWLOC_TYPE_DEPTH_UNKNOWN) {
            find_leastload_socket(tree, root, &result);
        } else {
            if ((depth_nodes) < (depth_sockets)) {
                find_leastload_node(tree, root, &result);
                node = result;
                find_leastload_socket(tree, node, &result);
            } else {
                find_leastload_socket(tree, root, &result);
                sockets = result;
                find_leastload_node(tree, sockets, &result);
            }
        }
        core_parent = result;
        find_leastload_core(tree, core_parent, &result);
        core = result;
        find_leastload_pu(tree, core, &result);
        insert_core_mapping(k, result, tree);
        k++;
    }
}

void map_bunch_load(obj_attribute_type * tree)
{
    int i, j, k, per = 0;
    int per_socket_node, depth_pus, num_pus = 0;
    float current_socketornode_load = 0, current_core_load = 0;
    int depth_cores, depth_sockets, depth_nodes, num_cores = 0, num_sockets =
        0, num_nodes = 0;
    hwloc_obj_t root, node, sockets, core_parent, core, pu, result;

    root = hwloc_get_root_obj(topology);

    depth_nodes = hwloc_get_type_depth(topology, HWLOC_OBJ_NODE);
    if (depth_nodes != HWLOC_TYPE_DEPTH_UNKNOWN) {
        num_nodes = hwloc_get_nbobjs_by_depth(topology, depth_nodes);
    }

    depth_sockets = hwloc_get_type_depth(topology, HWLOC_OBJ_SOCKET);
    if (depth_sockets != HWLOC_TYPE_DEPTH_UNKNOWN) {
        num_sockets = hwloc_get_nbobjs_by_depth(topology, depth_sockets);
    }

    depth_cores = hwloc_get_type_depth(topology, HWLOC_OBJ_CORE);
    if (depth_cores != HWLOC_TYPE_DEPTH_UNKNOWN) {
        num_cores = hwloc_get_nbobjs_by_depth(topology, depth_cores);
    }

    depth_pus = hwloc_get_type_depth(topology, HWLOC_OBJ_PU);
    if (depth_pus != HWLOC_TYPE_DEPTH_UNKNOWN) {
        num_pus = hwloc_get_nbobjs_by_depth(topology, depth_pus);
    }

    k = 0;
    /*Assume: there is always existing SOCKET, but not always existing NUMANODE(like Clovertown). */
    while (k < num_cores) {
        if (depth_nodes == HWLOC_TYPE_DEPTH_UNKNOWN) {
            find_leastload_socket(tree, root, &result);
            core_parent = result;
            per = num_cores / num_sockets;
            for (i = 0; (i < per) && (k < num_cores); i++) {
                find_leastload_core(tree, core_parent, &result);
                core = result;
                find_leastload_pu(tree, core, &result);
                pu = result;
                if (i == 0) {
                    current_core_load =
                        tree[depth_pus * num_pus + pu->logical_index].load;
                    insert_core_mapping(k, pu, tree);
                    k++;
                } else {
                    if (compare_float
                        (tree[depth_pus * num_pus + pu->logical_index].load,
                         current_core_load)) {
                        insert_core_mapping(k, pu, tree);
                        k++;
                    }
                }
            }
        } else {
            if ((depth_nodes) < (depth_sockets)) {
                find_leastload_node(tree, root, &result);
                node = result;
                per_socket_node = num_sockets / num_nodes;
                for (j = 0; (j < per_socket_node) && (k < num_cores); j++) {
                    find_leastload_socket(tree, node, &result);
                    sockets = result;
                    if (j == 0) {
                        current_socketornode_load =
                            tree[depth_sockets * num_sockets +
                                 sockets->logical_index].load;
                        per = num_cores / num_sockets;
                        for (i = 0; (i < per) && (k < num_cores); i++) {
                            find_leastload_core(tree, sockets, &result);
                            core = result;
                            find_leastload_pu(tree, core, &result);
                            pu = result;
                            if (i == 0) {
                                current_core_load =
                                    tree[depth_pus * num_pus + pu->logical_index].load;
                                insert_core_mapping(k, pu, tree);
                                k++;
                            } else {
                                if (compare_float
                                    (tree[depth_pus * num_pus + pu->logical_index].load,
                                     current_core_load)) {
                                    insert_core_mapping(k, pu, tree);
                                    k++;
                                }
                            }
                        }
                    } else {
                        if (compare_float
                            (tree
                             [depth_sockets * num_sockets + sockets->logical_index].load,
                             current_socketornode_load)) {
                            for (i = 0; (i < per) && (k < num_cores); i++) {
                                find_leastload_core(tree, sockets, &result);
                                core = result;
                                find_leastload_pu(tree, core, &result);
                                pu = result;
                                if (i == 0) {
                                    current_core_load =
                                        tree[depth_pus * num_pus +
                                             pu->logical_index].load;
                                    insert_core_mapping(k, pu, tree);
                                    k++;
                                } else {
                                    if (compare_float
                                        (tree
                                         [depth_pus * num_pus + pu->logical_index].load,
                                         current_core_load)) {
                                        insert_core_mapping(k, pu, tree);
                                        k++;
                                    }
                                }
                            }

                        }
                    }
                }
            } else {    // depth_nodes > depth_sockets
                find_leastload_socket(tree, root, &result);
                sockets = result;
                per_socket_node = num_nodes / num_sockets;
                for (j = 0; (j < per_socket_node) && (k < num_cores); j++) {
                    find_leastload_node(tree, sockets, &result);
                    node = result;
                    if (j == 0) {
                        current_socketornode_load =
                            tree[depth_nodes * num_nodes + node->logical_index].load;
                        per = num_cores / num_sockets;
                        for (i = 0; (i < per) && (k < num_cores); i++) {
                            find_leastload_core(tree, node, &result);
                            core = result;
                            find_leastload_pu(tree, core, &result);
                            pu = result;
                            if (i == 0) {
                                current_core_load =
                                    tree[depth_pus * num_pus + pu->logical_index].load;
                                insert_core_mapping(k, pu, tree);
                                k++;
                            } else {
                                if (compare_float
                                    (tree[depth_pus * num_pus + pu->logical_index].load,
                                     current_core_load)) {
                                    insert_core_mapping(k, pu, tree);
                                    k++;
                                }
                            }
                        }
                    } else {
                        if (compare_float
                            (tree[depth_nodes * num_nodes + node->logical_index].load,
                             current_socketornode_load)) {
                            for (i = 0; (i < per) && (k < num_cores); i++) {
                                find_leastload_core(tree, node, &result);
                                core = result;
                                find_leastload_pu(tree, core, &result);
                                pu = result;
                                if (i == 0) {
                                    current_core_load =
                                        tree[depth_pus * num_pus +
                                             pu->logical_index].load;
                                    insert_core_mapping(k, pu, tree);
                                    k++;
                                } else {
                                    if (compare_float
                                        (tree
                                         [depth_pus * num_pus + pu->logical_index].load,
                                         current_core_load)) {
                                        insert_core_mapping(k, pu, tree);
                                        k++;
                                    }
                                }
                            }
                        }
                    }
                }
            }   /* depth_nodes > depth_sockets */
        }
    }   /* while */
}

/*
 * Compare two hwloc_obj_t of type HWLOC_OBJ_PU according to sibling_rank, used with qsort
 */
static int cmpproc_smt(const void *a, const void *b)
{
    hwloc_obj_t pa = *(hwloc_obj_t *) a;
    hwloc_obj_t pb = *(hwloc_obj_t *) b;
    return (pa->sibling_rank ==
            pb->sibling_rank) ? pa->os_index - pb->os_index : pa->sibling_rank -
        pb->sibling_rank;
}

static int cmpdepth_smt(const void *a, const void *b)
{
    ancestor_type pa = *(ancestor_type *) a;
    ancestor_type pb = *(ancestor_type *) b;
    if ((pa.ancestor)->depth > (pb.ancestor)->depth) {
        return -1;
    } else if ((pa.ancestor)->depth < (pb.ancestor)->depth) {
        return 1;
    } else {
        return 0;
    }
}

static int cmparity_smt(const void *a, const void *b)
{
    ancestor_type pa = *(ancestor_type *) a;
    ancestor_type pb = *(ancestor_type *) b;
    if ((pa.ancestor)->arity > (pb.ancestor)->arity) {
        return -1;
    } else if ((pa.ancestor)->arity < (pb.ancestor)->arity) {
        return 1;
    } else {
        return 0;
    }
}

static void get_first_obj_bunch(hwloc_obj_t * result)
{
    hwloc_obj_t *objs;
    ancestor_type *array;
    int i, j, k, num_objs, num_ancestors;

    if ((num_objs = hwloc_get_nbobjs_by_type(topology, HWLOC_OBJ_PU)) <= 0) {
        return;
    }

    if ((objs = (hwloc_obj_t *) MPIU_Malloc(num_objs * sizeof(hwloc_obj_t))) == NULL) {
        return;
    }

    for (i = 0; i < num_objs; i++) {
        objs[i] = hwloc_get_obj_by_type(topology, HWLOC_OBJ_PU, i);
    }

    num_ancestors = num_objs * (num_objs - 1) / 2;

    if ((array =
         (ancestor_type *) MPIU_Malloc(num_ancestors * sizeof(ancestor_type))) == NULL) {
        return;
    }

    k = 0;
    for (i = 0; i < (num_objs - 1); i++) {
        for (j = i + 1; j < num_objs; j++) {
            array[k].obja = objs[i];
            array[k].objb = objs[j];
            array[k].ancestor = hwloc_get_common_ancestor_obj(topology, objs[i], objs[j]);
            k++;
        }
    }

    qsort(array, num_ancestors, sizeof(ancestor_type), cmpdepth_smt);

    for (i = 0; i < (num_ancestors - 1); i++) {
        if ((array[i + 1].ancestor)->depth < (array[i].ancestor)->depth) {
            break;
        }
    }

    qsort(array, (i + 1), sizeof(ancestor_type), cmparity_smt);

    *result = array[0].obja;

    MPIU_Free(objs);
    MPIU_Free(array);
    return;
}

static void get_first_socket_bunch(hwloc_obj_t * result, hwloc_obj_type_t binding_level)
{
    hwloc_obj_t *objs;
    ancestor_type *array;
    int i, j, k, num_objs, num_ancestors;

    if ((num_objs = hwloc_get_nbobjs_by_type(topology, binding_level)) <= 0) {
        return;
    }

    if ((objs = (hwloc_obj_t *) MPIU_Malloc(num_objs * sizeof(hwloc_obj_t))) == NULL) {
        return;
    }

    for (i = 0; i < num_objs; i++) {
        objs[i] = hwloc_get_obj_by_type(topology, binding_level, i);
    }

    num_ancestors = num_objs * (num_objs - 1) / 2;

    if ((array =
         (ancestor_type *) MPIU_Malloc(num_ancestors * sizeof(ancestor_type))) == NULL) {
        return;
    }

    k = 0;
    for (i = 0; i < (num_objs - 1); i++) {
        for (j = i + 1; j < num_objs; j++) {
            array[k].obja = objs[i];
            array[k].objb = objs[j];
            array[k].ancestor = hwloc_get_common_ancestor_obj(topology, objs[i], objs[j]);
            k++;
        }
    }

    qsort(array, num_ancestors, sizeof(ancestor_type), cmpdepth_smt);

    for (i = 0; i < (num_ancestors - 1); i++) {
        if ((array[i + 1].ancestor)->depth < (array[i].ancestor)->depth) {
            break;
        }
    }

    if (i < num_ancestors - 1)
        qsort(array, (i + 1), sizeof(ancestor_type), cmparity_smt);

    *result = array[0].obja;

    MPIU_Free(objs);
    MPIU_Free(array);
    return;
}

/*
 * Yields "scatter" affinity scenario in core_mapping.
 */
void map_scatter_core(int num_cpus)
{
    hwloc_obj_t *objs, obj, a;
    unsigned *pdist, maxd;
    int i, j, ix, jp, d, s;

    /* Init and load HWLOC_OBJ_PU objects */
    if ((objs = (hwloc_obj_t *) MPIU_Malloc(num_cpus * sizeof(hwloc_obj_t *))) == NULL)
        return;

    obj = NULL;
    i = 0;
    while ((obj = hwloc_get_next_obj_by_type(topology, HWLOC_OBJ_PU, obj)) != NULL)
        objs[i++] = obj;
    if (i != num_cpus) {
        MPIU_Free(objs);
        return;
    }

    /* Sort HWLOC_OBJ_PU objects according to sibling_rank */
    qsort(objs, num_cpus, sizeof(hwloc_obj_t *), cmpproc_smt);

    /* Init cumulative distances */
    if ((pdist = (unsigned *) MPIU_Malloc(num_cpus * sizeof(unsigned))) == NULL) {
        MPIU_Free(objs);
        return;
    }

    /* Loop over objects, ix is index in objs where sorted objects start */
    ix = num_cpus;
    s = -1;
    while (ix > 0) {
        /* If new group of SMT processors starts, zero distances */
        if (s != objs[0]->sibling_rank) {
            s = objs[0]->sibling_rank;
            for (j = 0; j < ix; j++)
                pdist[j] = 0;
        }
        /*
         * Determine object that has max. distance to all already stored objects.
         * Consider only groups of SMT processors with same sibling_rank.
         */
        maxd = 0;
        jp = 0;
        for (j = 0; j < ix; j++) {
            if ((j) && (objs[j - 1]->sibling_rank != objs[j]->sibling_rank))
                break;
            if (pdist[j] > maxd) {
                maxd = pdist[j];
                jp = j;
            }
        }

        /* Rotate found object to the end of the list, map out found object from distances */
        obj = objs[jp];
        for (j = jp; j < num_cpus - 1; j++) {
            objs[j] = objs[j + 1];
            pdist[j] = pdist[j + 1];
        }
        objs[j] = obj;
        ix--;

        /*
         * Update cumulative distances of all remaining objects with new stored one.
         * If two HWLOC_OBJ_PU objects don't share a common ancestor, the topology is broken.
         * Our scheme cannot be used in this case.
         */
        for (j = 0; j < ix; j++) {
            if ((a = hwloc_get_common_ancestor_obj(topology, obj, objs[j])) == NULL) {
                MPIU_Free(pdist);
                MPIU_Free(objs);
                return;
            }
            d = objs[j]->depth + obj->depth - 2 * a->depth;
            pdist[j] += d * d;
        }
    }

    /* Collect os_indexes into core_mapping */
    for (i = 0; i < num_cpus; i++) {
        core_mapping[i] = objs[i]->os_index;
    }

    MPIU_Free(pdist);
    MPIU_Free(objs);
    return;
}

void map_scatter_socket(int num_sockets, hwloc_obj_type_t binding_level)
{
    hwloc_obj_t *objs, obj, a;
    unsigned *pdist, maxd;
    int i, j, ix, jp, d, s, num_cores;

    /* Init and load HWLOC_OBJ_SOCKET or HWLOC_OBJ_NODE objects */
    if ((objs = (hwloc_obj_t *) MPIU_Malloc(num_sockets * sizeof(hwloc_obj_t *))) == NULL)
        return;

    if ((num_cores = hwloc_get_nbobjs_by_type(topology, HWLOC_OBJ_CORE)) <= 0) {
        return;
    }

    obj = NULL;
    i = 0;
    while ((obj = hwloc_get_next_obj_by_type(topology, binding_level, obj)) != NULL)
        objs[i++] = obj;
    if (i != num_sockets) {
        MPIU_Free(objs);
        return;
    }

    /* Sort HWLOC_OBJ_SOCKET or HWLOC_OBJ_NODE objects according to sibling_rank */
    qsort(objs, num_sockets, sizeof(hwloc_obj_t *), cmpproc_smt);

    /* Init cumulative distances */
    if ((pdist = (unsigned *) MPIU_Malloc(num_sockets * sizeof(unsigned))) == NULL) {
        MPIU_Free(objs);
        return;
    }

    /* Loop over objects, ix is index in objs where sorted objects start */
    ix = num_sockets;
    s = -1;
    while (ix > 0) {
        /* If new group of SMT processors starts, zero distances */
        if (s != objs[0]->sibling_rank) {
            s = objs[0]->sibling_rank;
            for (j = 0; j < ix; j++)
                pdist[j] = 0;
        }
        /*
         * Determine object that has max. distance to all already stored objects.
         * Consider only groups of SMT processors with same sibling_rank.
         */
        maxd = 0;
        jp = 0;
        for (j = 0; j < ix; j++) {
            if ((j) && (objs[j - 1]->sibling_rank != objs[j]->sibling_rank))
                break;
            if (pdist[j] > maxd) {
                maxd = pdist[j];
                jp = j;
            }
        }

        /* Rotate found object to the end of the list, map out found object from distances */
        obj = objs[jp];
        for (j = jp; j < num_sockets - 1; j++) {
            objs[j] = objs[j + 1];
            pdist[j] = pdist[j + 1];
        }
        objs[j] = obj;
        ix--;

        /*
         * Update cumulative distances of all remaining objects with new stored one.
         * If two HWLOC_OBJ_SOCKET or HWLOC_OBJ_NODE objects don't share a common ancestor, the topology is broken.
         * Our scheme cannot be used in this case.
         */
        for (j = 0; j < ix; j++) {
            if ((a = hwloc_get_common_ancestor_obj(topology, obj, objs[j])) == NULL) {
                MPIU_Free(pdist);
                MPIU_Free(objs);
                return;
            }
            d = objs[j]->depth + obj->depth - 2 * a->depth;
            pdist[j] += d * d;
        }
    }

    /* Collect os_indexes into core_mapping */
    for (i = 0, j = 0; i < num_cores; i++, j++) {
        if (j == num_sockets) {
            j = 0;
        }
        core_mapping[i] = hwloc_bitmap_to_ulong((hwloc_const_bitmap_t) (objs[j]->cpuset));
    }

    MPIU_Free(pdist);
    MPIU_Free(objs);
    return;
}

 /*
  * Yields "bunch" affinity scenario in core_mapping.
  */
void map_bunch_core(int num_cpus)
{
    hwloc_obj_t *objs, obj, a;
    unsigned *pdist, mind;
    int i, j, ix, jp, d, s, num_cores, num_pus;

    /* Init and load HWLOC_OBJ_PU objects */
    if ((objs = (hwloc_obj_t *) MPIU_Malloc(num_cpus * sizeof(hwloc_obj_t *))) == NULL)
        return;

    obj = NULL;
    i = 0;

    if ((num_cores = hwloc_get_nbobjs_by_type(topology, HWLOC_OBJ_CORE)) <= 0) {
        MPIU_Free(objs);
        return;
    }

    if ((num_pus = hwloc_get_nbobjs_by_type(topology, HWLOC_OBJ_PU)) <= 0) {
        MPIU_Free(objs);
        return;
    }

    /* SMT Disabled */
    if (num_cores == num_pus) {

        get_first_obj_bunch(&obj);

        if (obj == NULL) {
            MPIU_Free(objs);
            return;
        }

        objs[i] = obj;
        i++;

        while ((obj = hwloc_get_next_obj_by_type(topology, HWLOC_OBJ_PU, obj)) != NULL) {
            objs[i] = obj;
            i++;
        }

        obj = NULL;
        while (i != num_cpus) {
            obj = hwloc_get_next_obj_by_type(topology, HWLOC_OBJ_PU, obj);
            objs[i++] = obj;
        }

        if (i != num_cpus) {
            MPIU_Free(objs);
            return;
        }

    } else {    /* SMT Enabled */

        while ((obj = hwloc_get_next_obj_by_type(topology, HWLOC_OBJ_PU, obj)) != NULL)
            objs[i++] = obj;

        if (i != num_cpus) {
            MPIU_Free(objs);
            return;
        }

        /* Sort HWLOC_OBJ_PU objects according to sibling_rank */
        qsort(objs, num_cpus, sizeof(hwloc_obj_t *), cmpproc_smt);
    }

    /* Init cumulative distances */
    if ((pdist = (unsigned *) MPIU_Malloc(num_cpus * sizeof(unsigned))) == NULL) {
        MPIU_Free(objs);
        return;
    }

    /* Loop over objects, ix is index in objs where sorted objects start */
    ix = num_cpus;
    s = -1;
    while (ix > 0) {
        /* If new group of SMT processors starts, zero distances */
        if (s != objs[0]->sibling_rank) {
            s = objs[0]->sibling_rank;
            for (j = 0; j < ix; j++)
                pdist[j] = UINT_MAX;
        }
        /*
         * Determine object that has min. distance to all already stored objects.
         * Consider only groups of SMT processors with same sibling_rank.
         */
        mind = UINT_MAX;
        jp = 0;
        for (j = 0; j < ix; j++) {
            if ((j) && (objs[j - 1]->sibling_rank != objs[j]->sibling_rank))
                break;
            if (pdist[j] < mind) {
                mind = pdist[j];
                jp = j;
            }
        }

        /* Rotate found object to the end of the list, map out found object from distances */
        obj = objs[jp];
        for (j = jp; j < num_cpus - 1; j++) {
            objs[j] = objs[j + 1];
            pdist[j] = pdist[j + 1];
        }
        objs[j] = obj;
        ix--;

        /*
         * Update cumulative distances of all remaining objects with new stored one.
         * If two HWLOC_OBJ_PU objects don't share a common ancestor, the topology is broken.
         * Our scheme cannot be used in this case.
         */
        for (j = 0; j < ix; j++) {
            if ((a = hwloc_get_common_ancestor_obj(topology, obj, objs[j])) == NULL) {
                MPIU_Free(pdist);
                MPIU_Free(objs);
                return;
            }
            d = objs[j]->depth + obj->depth - 2 * a->depth;
            pdist[j] += d * d;
        }
    }

    /* Collect os_indexes into core_mapping */
    for (i = 0; i < num_cpus; i++) {
        core_mapping[i] = objs[i]->os_index;
    }

    MPIU_Free(pdist);
    MPIU_Free(objs);
    return;
}

int check_num_child(hwloc_obj_t obj)
{
    int i = 0, k, num_cores;

    if ((num_cores = hwloc_get_nbobjs_by_type(topology, HWLOC_OBJ_CORE)) <= 0) {
        return 0;
    }

    for (k = 0; k < num_cores; k++) {
        if (hwloc_bitmap_isset((hwloc_const_bitmap_t) (obj->cpuset), k)) {
            i++;
        }
    }

    return i;
}

void map_bunch_socket(int num_sockets, hwloc_obj_type_t binding_level)
{
    hwloc_obj_t *objs, obj, a;
    unsigned *pdist, mind;
    int i, j, ix, jp, d, s, num_cores, num_pus;

    /* Init and load HWLOC_OBJ_PU objects */
    if ((objs = (hwloc_obj_t *) MPIU_Malloc(num_sockets * sizeof(hwloc_obj_t *))) == NULL)
        return;

    obj = NULL;
    i = 0;

    if ((num_cores = hwloc_get_nbobjs_by_type(topology, HWLOC_OBJ_CORE)) <= 0) {
        MPIU_Free(objs);
        return;
    }

    if ((num_pus = hwloc_get_nbobjs_by_type(topology, HWLOC_OBJ_PU)) <= 0) {
        MPIU_Free(objs);
        return;
    }

    /* SMT Disabled */
    if (num_cores == num_pus) {

        get_first_socket_bunch(&obj, binding_level);

        if (obj == NULL) {
            MPIU_Free(objs);
            return;
        }

        objs[i] = obj;
        i++;

        while ((obj = hwloc_get_next_obj_by_type(topology, binding_level, obj)) != NULL) {
            objs[i] = obj;
            i++;
        }

        obj = NULL;
        while (i != num_sockets) {
            obj = hwloc_get_next_obj_by_type(topology, binding_level, obj);
            objs[i++] = obj;
        }

        if (i != num_sockets) {
            MPIU_Free(objs);
            return;
        }

    } else {    /* SMT Enabled */

        while ((obj = hwloc_get_next_obj_by_type(topology, binding_level, obj)) != NULL)
            objs[i++] = obj;

        if (i != num_sockets) {
            MPIU_Free(objs);
            return;
        }

        /* Sort HWLOC_OBJ_SOCKET or HWLOC_OBJ_NODE objects according to sibling_rank */
        qsort(objs, num_sockets, sizeof(hwloc_obj_t *), cmpproc_smt);

    }

    /* Init cumulative distances */
    if ((pdist = (unsigned *) MPIU_Malloc(num_sockets * sizeof(unsigned))) == NULL) {
        MPIU_Free(objs);
        return;
    }

    /* Loop over objects, ix is index in objs where sorted objects start */
    ix = num_sockets;
    s = -1;
    while (ix > 0) {
        /* If new group of SMT processors starts, zero distances */
        if (s != objs[0]->sibling_rank) {
            s = objs[0]->sibling_rank;
            for (j = 0; j < ix; j++)
                pdist[j] = UINT_MAX;
        }
        /*
         * Determine object that has min. distance to all already stored objects.
         * Consider only groups of SMT processors with same sibling_rank.
         */
        mind = UINT_MAX;
        jp = 0;
        for (j = 0; j < ix; j++) {
            if ((j) && (objs[j - 1]->sibling_rank != objs[j]->sibling_rank))
                break;
            if (pdist[j] < mind) {
                mind = pdist[j];
                jp = j;
            }
        }

        /* Rotate found object to the end of the list, map out found object from distances */
        obj = objs[jp];
        for (j = jp; j < num_sockets - 1; j++) {
            objs[j] = objs[j + 1];
            pdist[j] = pdist[j + 1];
        }
        objs[j] = obj;
        ix--;

        /*
         * Update cumulative distances of all remaining objects with new stored one.
         * If two HWLOC_OBJ_SOCKET or HWLOC_OBJ_NODE objects don't share a common ancestor, the topology is broken.
         * Our scheme cannot be used in this case.
         */
        for (j = 0; j < ix; j++) {
            if ((a = hwloc_get_common_ancestor_obj(topology, obj, objs[j])) == NULL) {
                MPIU_Free(pdist);
                MPIU_Free(objs);
                return;
            }
            d = objs[j]->depth + obj->depth - 2 * a->depth;
            pdist[j] += d * d;
        }
    }

    /* Collect os_indexes into core_mapping */
    int num_child_in_socket[num_sockets];

    for (i = 0; i < num_sockets; i++) {
        num_child_in_socket[i] = check_num_child(objs[i]);
    }

    for (i = 1; i < num_sockets; i++)
        num_child_in_socket[i] += num_child_in_socket[i - 1];

    for (i = 0, j = 0; i < num_cores; i++) {
        if (i == num_child_in_socket[j]) {
            j++;
        }
        core_mapping[i] = hwloc_bitmap_to_ulong((hwloc_const_bitmap_t) (objs[j]->cpuset));
    }

    MPIU_Free(pdist);
    MPIU_Free(objs);
    return;
}

static int num_digits(unsigned long numcpus)
{
    int n_digits = 0;
    while (numcpus > 0) {
        n_digits++;
        numcpus /= 10;
    }
    return n_digits;
}

int get_cpu_mapping_hwloc(long N_CPUs_online, hwloc_topology_t tp)
{
    unsigned topodepth = -1, depth = -1;
    int num_processes = 0, rc = 0, i;
    int num_sockets = 0;
    int num_numanodes = 0;
    int num_cpus = 0;
    char *s;
    struct dirent **namelist;
    pid_t pid;
    obj_attribute_type *tree = NULL;
    char *value;

    /* Determine topology depth */
    topodepth = hwloc_topology_get_depth(tp);
    if (topodepth == HWLOC_TYPE_DEPTH_UNKNOWN) {
        fprintf(stderr, "Warning: %s: Failed to determine topology depth.\n", __func__);
        return (topodepth);
    }

    /* Count number of (logical) processors */
    depth = hwloc_get_type_depth(tp, HWLOC_OBJ_PU);

    if (depth == HWLOC_TYPE_DEPTH_UNKNOWN) {
        fprintf(stderr, "Warning: %s: Failed to determine number of processors.\n",
                __func__);
        return (depth);
    }
    if ((num_cpus = hwloc_get_nbobjs_by_type(tp, HWLOC_OBJ_PU)) <= 0) {
        fprintf(stderr, "Warning: %s: Failed to determine number of processors.\n",
                __func__);
        return -1;
    }

    /* Count number of sockets */
    depth = hwloc_get_type_depth(tp, HWLOC_OBJ_SOCKET);
    if (depth == HWLOC_TYPE_DEPTH_UNKNOWN) {
        fprintf(stderr, "Warning: %s: Failed to determine number of sockets.\n",
                __func__);
        return (depth);
    } else {
        num_sockets = hwloc_get_nbobjs_by_depth(tp, depth);
    }

    /* Count number of numanodes */
    depth = hwloc_get_type_depth(tp, HWLOC_OBJ_NODE);
    if (depth == HWLOC_TYPE_DEPTH_UNKNOWN) {
        num_numanodes = -1;
    } else {
        num_numanodes = hwloc_get_nbobjs_by_depth(tp, depth);
    }

    if (s_cpu_mapping == NULL) {
        /* We need to do allocate memory for the custom_cpu_mapping array
         * and determine the current load on the different cpu's only
         * when the user has not specified a mapping string. If the user
         * has provided a mapping string, it overrides everything.
         */
        /*TODO: might need a better representation as number of cores per node increases */
        unsigned long long_max = ULONG_MAX;
        int n_digits = num_digits(long_max);
        custom_cpu_mapping =
            MPIU_Malloc(sizeof(char) * num_cpus * (n_digits + 1) + 1);
        if (custom_cpu_mapping == NULL) {
            goto error_free;
        }
        MPIU_Memset(custom_cpu_mapping, 0,
                    sizeof(char) * num_cpus * (n_digits + 1) + 1);
        core_mapping = (unsigned long *) MPIU_Malloc(num_cpus * sizeof(unsigned long));
        if (core_mapping == NULL) {
            goto error_free;
        }
        for (i = 0; i < num_cpus; i++) {
            core_mapping[i] = -1;
        }

        tree = MPIU_Malloc(num_cpus * topodepth * sizeof(obj_attribute_type));
        if (tree == NULL) {
            goto error_free;
        }
        for (i = 0; i < num_cpus * topodepth; i++) {
            tree[i].obj = NULL;
            tree[i].load = 0;
            CPU_ZERO(&(tree[i].cpuset));
        }

        if (!(obj_tree = (int *) MPIU_Malloc(num_cpus * topodepth * sizeof(*obj_tree)))) {
            goto error_free;
        }
        for (i = 0; i < num_cpus * topodepth; i++) {
            obj_tree[i] = -1;
        }

        ip = 0;

        /* MV2_ENABLE_LEASTLOAD: map_bunch/scatter or map_bunch/scatter_load */
        if ((value = getenv("MV2_ENABLE_LEASTLOAD")) != NULL) {
            mv2_enable_leastload = atoi(value);
            if (mv2_enable_leastload != 1) {
                mv2_enable_leastload = 0;
            }
        }

        /* MV2_ENABLE_LEASTLOAD=1, map_bunch_load or map_scatter_load is used */
        if (mv2_enable_leastload == 1) {
            /*
             * Get all processes' pid and cpuset.
             * Get numanode, socket, and core current load according to processes running on it.
             */
            num_processes = scandir("/proc", &namelist, pid_filter, alphasort);
            if (num_processes < 0) {
                fprintf(stderr, "Warning: %s: Failed to scandir /proc.\n", __func__);
                return -1;
            } else {
                int status;
                cpu_set_t pid_cpuset;
                CPU_ZERO(&pid_cpuset);

                /* Get cpuset for each running process. */
                for (i = 0; i < num_processes; i++) {
                    pid = atol(namelist[i]->d_name);
                    status = sched_getaffinity(pid, sizeof(pid_cpuset), &pid_cpuset);
                    /* Process completed. */
                    if (status < 0) {
                        continue;
                    }
                    cac_load(tree, pid_cpuset);
                }
                while (num_processes--) {
                    MPIU_Free(namelist[num_processes]);
                }
                MPIU_Free(namelist);
            }

            if (mv2_binding_policy == POLICY_SCATTER) {
                map_scatter_load(tree);
            } else if (mv2_binding_policy == POLICY_BUNCH) {
                map_bunch_load(tree);
            } else {
                goto error_free;
            }
        } else {
            /* MV2_ENABLE_LEASTLOAD != 1 or MV2_ENABLE_LEASTLOAD == NULL, map_bunch or map_scatter is used */
            if (mv2_binding_policy == POLICY_SCATTER) {
                /* Scatter */
                hwloc_obj_type_t binding_level = HWLOC_OBJ_SOCKET;
                if (mv2_binding_level == LEVEL_SOCKET) {
                    map_scatter_socket(num_sockets, binding_level);
                } else if (mv2_binding_level == LEVEL_NUMANODE) {
                    if (num_numanodes == -1) {
                        /* There is not numanode, fallback to socket */
                        map_scatter_socket(num_sockets, binding_level);
                    } else {
                        binding_level = HWLOC_OBJ_NODE;
                        map_scatter_socket(num_numanodes, binding_level);
                    }
                } else {
                    map_scatter_core(num_cpus);
                }

            } else if (mv2_binding_policy == POLICY_BUNCH) {
                /* Bunch */
                hwloc_obj_type_t binding_level = HWLOC_OBJ_SOCKET;
                if (mv2_binding_level == LEVEL_SOCKET) {
                    map_bunch_socket(num_sockets, binding_level);
                } else if (mv2_binding_level == LEVEL_NUMANODE) {
                    if (num_numanodes == -1) {
                        /* There is not numanode, fallback to socket */
                        map_bunch_socket(num_sockets, binding_level);
                    } else {
                        binding_level = HWLOC_OBJ_NODE;
                        map_bunch_socket(num_numanodes, binding_level);
                    }
                } else {
                    map_bunch_core(num_cpus);
                }
            } else {
                goto error_free;
            }
        }

        /* Assemble custom_cpu_mapping string */
        s = custom_cpu_mapping;
        for (i = 0; i < num_cpus; i++) {
            sprintf(s, "%lu:", core_mapping[i]);
            s = custom_cpu_mapping + strlen(custom_cpu_mapping);
        }
        i = strlen(custom_cpu_mapping);
        if (i) {
            custom_cpu_mapping[i - 1] = '\0';
        }
    }

    /* Done */
    rc = MPI_SUCCESS;

  error_free:
    if (core_mapping != NULL) {
        MPIU_Free(core_mapping);
    }
    if (tree != NULL) {
        MPIU_Free(tree);
    }
    if (obj_tree) {
        MPIU_Free(obj_tree);
    }

    MPIU_DBG_MSG_FMT(OTHER, TYPICAL,
                     (MPIU_DBG_FDEST,
                      "num_cpus=%d, num_sockets=%d, custom_cpu_mapping=\"%s\"", num_cpus,
                      num_sockets, custom_cpu_mapping));

    return rc;
}


int get_cpu_mapping(long N_CPUs_online)
{
    char line[MAX_LINE_LENGTH];
    char input[MAX_NAME_LENGTH];
    char bogus1[MAX_NAME_LENGTH];
    char bogus2[MAX_NAME_LENGTH];
    char bogus3[MAX_NAME_LENGTH];
    int physical_id;            //return value
    int mapping[N_CPUs_online];
    int core_index = 0;
    cpu_type_t cpu_type = 0;
    int model;
    int vendor_set = 0, model_set = 0, num_cpus = 0;

    FILE *fp = fopen(CONFIG_FILE, "r");
    if (fp == NULL) {
        printf("can not open cpuinfo file \n");
        return 0;
    }

    MPIU_Memset(mapping, 0, sizeof(mapping));
    custom_cpu_mapping = (char *) MPIU_Malloc(sizeof(char) * N_CPUs_online * 2);
    if (custom_cpu_mapping == NULL) {
        return 0;
    }
    MPIU_Memset(custom_cpu_mapping, 0, sizeof(char) * N_CPUs_online * 2);

    while (!feof(fp)) {
        MPIU_Memset(line, 0, MAX_LINE_LENGTH);
        fgets(line, MAX_LINE_LENGTH, fp);

        MPIU_Memset(input, 0, MAX_NAME_LENGTH);
        sscanf(line, "%s", input);

        if (!vendor_set) {
            if (strcmp(input, "vendor_id") == 0) {
                MPIU_Memset(input, 0, MAX_NAME_LENGTH);
                sscanf(line, "%s%s%s", bogus1, bogus2, input);

                if (strcmp(input, "AuthenticAMD") == 0) {
                    cpu_type = CPU_FAMILY_AMD;
                } else {
                    cpu_type = CPU_FAMILY_INTEL;
                }
                vendor_set = 1;
            }
        }

        if (!model_set) {
            if (strcmp(input, "model") == 0) {
                sscanf(line, "%s%s%d", bogus1, bogus2, &model);
                model_set = 1;
            }
        }

        if (strcmp(input, "physical") == 0) {
            sscanf(line, "%s%s%s%d", bogus1, bogus2, bogus3, &physical_id);
            mapping[core_index++] = physical_id;
        }
    }

    num_cpus = core_index;
    if (num_cpus == 4) {
        if ((memcmp(INTEL_XEON_DUAL_MAPPING, mapping, sizeof(int) * num_cpus) == 0)
            && (cpu_type == CPU_FAMILY_INTEL)) {
            strcpy(custom_cpu_mapping, "0:2:1:3");
        } else
            if ((memcmp(AMD_OPTERON_DUAL_MAPPING, mapping, sizeof(int) * num_cpus) == 0)
                && (cpu_type == CPU_FAMILY_AMD)) {
            strcpy(custom_cpu_mapping, "0:1:2:3");
        }
    } else if (num_cpus == 8) {
        if (cpu_type == CPU_FAMILY_INTEL) {
            if (model == CLOVERTOWN_MODEL) {
                if (memcmp(INTEL_CLOVERTOWN_MAPPING, mapping, sizeof(int) * num_cpus) ==
                    0) {
                    strcpy(custom_cpu_mapping, "0:1:4:5:2:3:6:7");
                }
            } else if (model == HARPERTOWN_MODEL) {
                if (memcmp(INTEL_HARPERTOWN_LEG_MAPPING, mapping, sizeof(int) * num_cpus)
                    == 0) {
                    strcpy(custom_cpu_mapping, "0:1:4:5:2:3:6:7");
                } else
                    if (memcmp
                        (INTEL_HARPERTOWN_COM_MAPPING, mapping,
                         sizeof(int) * num_cpus) == 0) {
                    strcpy(custom_cpu_mapping, "0:4:2:6:1:5:3:7");
                }
            } else if (model == NEHALEM_MODEL) {
                if (memcmp(INTEL_NEHALEM_LEG_MAPPING, mapping, sizeof(int) * num_cpus) ==
                    0) {
                    strcpy(custom_cpu_mapping, "0:2:4:6:1:3:5:7");
                } else
                    if (memcmp(INTEL_NEHALEM_COM_MAPPING, mapping, sizeof(int) * num_cpus)
                        == 0) {
                    strcpy(custom_cpu_mapping, "0:4:1:5:2:6:3:7");
                }
            }
        }
    } else if (num_cpus == 16) {
        if (cpu_type == CPU_FAMILY_INTEL) {
            if (model == NEHALEM_MODEL) {
                if (memcmp(INTEL_NEHALEM_LEG_MAPPING, mapping, sizeof(int) * num_cpus) ==
                    0) {
                    strcpy(custom_cpu_mapping, "0:2:4:6:1:3:5:7:8:10:12:14:9:11:13:15");
                } else
                    if (memcmp(INTEL_NEHALEM_COM_MAPPING, mapping, sizeof(int) * num_cpus)
                        == 0) {
                    strcpy(custom_cpu_mapping, "0:4:1:5:2:6:3:7:8:12:9:13:10:14:11:15");
                }
            }
        } else if (cpu_type == CPU_FAMILY_AMD) {
            if (memcmp(AMD_BARCELONA_MAPPING, mapping, sizeof(int) * num_cpus) == 0) {
                strcpy(custom_cpu_mapping, "0:1:2:3:4:5:6:7:8:9:10:11:12:13:14:15");
            }
        }
    }
    fclose(fp);

    return MPI_SUCCESS;
}

#if defined(CHANNEL_MRAIL)
int get_socket_id (int ib_socket, int cpu_socket, int num_sockets,
        tab_socket_t * tab_socket)
{
    extern int rdma_local_id, rdma_num_hcas;

    int rdma_num_proc_per_hca;
    int offset_id;
    int j;
    int socket_id = ib_socket;
    int delta = cpu_socket / tab_socket[ib_socket].num_hca;

    rdma_num_proc_per_hca = rdma_num_local_procs / rdma_num_hcas;

    if (rdma_num_local_procs % rdma_num_hcas) {
        rdma_num_proc_per_hca++;
    }

    offset_id = rdma_local_id % rdma_num_proc_per_hca;

    if (offset_id < delta) {
        return ib_socket;
    }

    for (j = 0; j < num_sockets - 1; j++) {
        socket_id = tab_socket[ib_socket].closest[j];

        if (tab_socket[socket_id].num_hca == 0) {
            offset_id -= delta;

            if (offset_id < delta) {
                return socket_id;
            }
        }
    }

    /*
     * Couldn't find a free socket, spread remaining processes
     */
    return rdma_local_id % num_sockets;
}

#undef FUNCNAME
#define FUNCNAME mv2_get_cpu_core_closest_to_hca
#undef FCNAME
#define FCNAME MPIDI_QUOTE(FUNCNAME)
int mv2_get_cpu_core_closest_to_hca(int my_local_id, int total_num_cores,
                                    int num_sockets, int depth_sockets)
{
    int i = 0, k = 0;
    int ib_hca_selected = 0;
    int selected_socket = 0;
    int cores_per_socket = 0;
    tab_socket_t *tab_socket = NULL;
    int linelen = strlen(custom_cpu_mapping);

    if (linelen < custom_cpu_mapping_line_max) {
        custom_cpu_mapping_line_max = linelen;
    }

    cores_per_socket = total_num_cores / num_sockets;

    /*
     * Make ib_hca_selected global or make this section a function
     */
    if (FIXED_MAPPING == rdma_rail_sharing_policy) {
        ib_hca_selected = rdma_process_binding_rail_offset /
                            rdma_num_rails_per_hca;
    } else {
        ib_hca_selected = 0;
    }

    tab_socket = (tab_socket_t*)MPIU_Malloc(num_sockets * sizeof(tab_socket_t));
    if (NULL == tab_socket) {
        fprintf(stderr, "could not allocate the socket table\n");
        return -1;
    }

    for (i = 0; i < num_sockets; i++) {
        tab_socket[i].num_hca = 0;

        for(k = 0; k < num_sockets; k++) {
            tab_socket[i].closest[k] = -1;
        }
    }

    for (i = 0; i < rdma_num_hcas; i++) {
        struct ibv_device * ibdev = mv2_MPIDI_CH3I_RDMA_Process.ib_dev[i];
        int socket_id = get_ib_socket(ibdev);
        /*
         * Make this information available globally
         */
        if (i == ib_hca_selected) {
            ib_socket_bind = socket_id;
        }
        tab_socket[socket_id].num_hca++;
    }

    hwloc_obj_t obj_src;
    hwloc_obj_t objs[num_sockets];
    char string[20];

    for (i = 0; i < num_sockets; i++) {
        obj_src = hwloc_get_obj_by_type(topology, HWLOC_OBJ_SOCKET,i);
        hwloc_get_closest_objs(topology, obj_src, (hwloc_obj_t *)&objs,
                                num_sockets - 1);

        for (k = 0; k < num_sockets - 1; k++) {
            hwloc_obj_type_snprintf(string, sizeof(string),
                                objs[k], 1);
            tab_socket[i].closest[k] = objs[k]->os_index;
        }
    }

    selected_socket = get_socket_id(ib_socket_bind, cores_per_socket,
                                    num_sockets, tab_socket);
    MPIU_Free(tab_socket);

    return selected_socket;
}
#endif /* defined(CHANNEL_MRAIL) */

#undef FUNCNAME
#define FUNCNAME mv2_get_assigned_cpu_core
#undef FCNAME
#define FCNAME MPIDI_QUOTE(FUNCNAME)
int mv2_get_assigned_cpu_core(int my_local_id, char *cpu_mapping, int max_cpu_map_len, char *tp_str)
{
    int i = 0, j = 0;
    char *cp = NULL;
    char *tp = cpu_mapping;
    long N_CPUs_online = sysconf(_SC_NPROCESSORS_ONLN);

    while (*tp != '\0') {
        i = 0;
        cp = tp;

        while (*cp != '\0' && *cp != ':' && i < max_cpu_map_len) {
            ++cp;
            ++i;
        }

        if (j == my_local_id) {
            strncpy(tp_str, tp, i);
            if (atoi(tp) < 0 ||
                ((mv2_binding_level == LEVEL_CORE) && (atoi(tp) >= N_CPUs_online))) {
                fprintf(stderr,
                        "Warning! : Core id %d does not exist on this architecture! \n",
                        atoi(tp));
                fprintf(stderr, "CPU Affinity is undefined \n");
                mv2_enable_affinity = 0;
                return -1;
            }
            tp_str[i] = '\0';
            return 0;
        }

        if (*cp == '\0') {
            break;
        }

        tp = cp;
        ++tp;
        ++j;
    }

    return -1;
}

#if defined(CHANNEL_MRAIL)
#undef FUNCNAME
#define FUNCNAME smpi_set_progress_thread_affinity
#undef FCNAME
#define FCNAME MPIDI_QUOTE(FUNCNAME)
int smpi_set_progress_thread_affinity()
{
    int mpi_errno = MPI_SUCCESS;
    hwloc_cpuset_t cpuset;

    /* Alloc cpuset */
    cpuset = hwloc_bitmap_alloc();
    /* Set cpuset to mv2_my_async_cpu_id */
    hwloc_bitmap_set(cpuset, mv2_my_async_cpu_id);
    /* Attachement progress thread to mv2_my_async_cpu_id */
    hwloc_set_thread_cpubind(topology, pthread_self(), cpuset, 0);
    /* Free cpuset */
    hwloc_bitmap_free(cpuset);

    return mpi_errno;
}

#undef FUNCNAME
#define FUNCNAME smpi_identify_allgather_local_core_ids
#undef FCNAME
#define FCNAME MPIDI_QUOTE(FUNCNAME)
int smpi_identify_allgather_local_core_ids(MPIDI_PG_t * pg)
{
    int mpi_errno = MPI_SUCCESS;
    int p = 0;
    MPIDI_VC_t *vc = NULL;
    MPI_Request *request = NULL;
    MPI_Status *status= NULL;
    int errflag = 0;
    MPI_Comm comm;
    MPID_Comm *comm_ptr=NULL;

    MPID_Comm_get_ptr(MPI_COMM_WORLD, comm_ptr );
    comm = comm_ptr->handle;

    /* Allocate memory */
    local_core_ids = MPIU_Malloc(g_smpi.num_local_nodes * sizeof(int));
    if (local_core_ids== NULL) {
        ibv_error_abort(GEN_EXIT_ERR, "Failed to allocate memory for local_core_ids\n");
    }
    request = MPIU_Malloc(g_smpi.num_local_nodes * 2 * sizeof(MPI_Request));
    if (request == NULL) {
        ibv_error_abort(GEN_EXIT_ERR, "Failed to allocate memory for requests\n");
    }
    status = MPIU_Malloc(g_smpi.num_local_nodes * 2 * sizeof(MPI_Status));
    if (request == NULL) {
        ibv_error_abort(GEN_EXIT_ERR, "Failed to allocate memory for statuses\n");
    }
    /* Perform intra-node allgather */
    for (p = 0; p < g_smpi.num_local_nodes; ++p) {
        MPIDI_PG_Get_vc(pg, g_smpi.l2g_rank[p], &vc);
        if (vc->smp.local_nodes >= 0) {
            mpi_errno = MPIC_Irecv((void*)&local_core_ids[vc->smp.local_nodes],
                                    1, MPI_INT, vc->pg_rank, MPIR_ALLGATHER_TAG,
                                    comm, &request[g_smpi.num_local_nodes+p]);
            if (mpi_errno) {
                MPIU_ERR_POP(mpi_errno);
            }
            mpi_errno = MPIC_Isend((void*)&mv2_my_cpu_id, 1, MPI_INT, vc->pg_rank,
                                    MPIR_ALLGATHER_TAG, comm, &request[p], &errflag);
            if (mpi_errno) {
                MPIU_ERR_POP(mpi_errno);
            }
        }
    }
    /* Wait for intra-node allgather to finish */
    mpi_errno = PMPI_Waitall(g_smpi.num_local_nodes*2, request, status);
    if (mpi_errno) {
        MPIU_ERR_POP(mpi_errno);
    }

fn_exit:
    if (request) {
        MPIU_Free(request);
    }
    if (status) {
        MPIU_Free(status);
    }
    return mpi_errno;
fn_fail:
    goto fn_exit;
}

#undef FUNCNAME
#define FUNCNAME smpi_identify_free_cores
#undef FCNAME
#define FCNAME MPIDI_QUOTE(FUNCNAME)
int smpi_identify_free_cores(hwloc_cpuset_t *sock_cpuset, hwloc_cpuset_t *free_sock_cpuset)
{
    int i = 0;
    int mpi_errno = MPI_SUCCESS;
    int num_sockets = -1;
    int depth_sockets = -1;
    hwloc_obj_t socket = NULL;
    hwloc_cpuset_t my_cpuset = NULL;
    char cpu_str[128];

    /* Alloc cpuset */
    my_cpuset = hwloc_bitmap_alloc();
    *sock_cpuset = hwloc_bitmap_alloc();
    /* Clear CPU set */
    hwloc_bitmap_zero(my_cpuset);
    hwloc_bitmap_zero(*sock_cpuset);
    /* Set cpuset to mv2_my_cpu_id */
    hwloc_bitmap_set(my_cpuset, mv2_my_cpu_id);

    depth_sockets   = hwloc_get_type_depth(topology, HWLOC_OBJ_SOCKET);
    num_sockets     = hwloc_get_nbobjs_by_depth(topology, depth_sockets);

    for (i = 0; i < num_sockets; ++i) {
        socket = hwloc_get_obj_by_depth(topology, depth_sockets, i);
        /* Find the list of CPUs we're allowed to use in the socket */
        hwloc_bitmap_and(*sock_cpuset, socket->online_cpuset, socket->allowed_cpuset);
        /* Find the socket the core I'm bound to resides on */
        if (hwloc_bitmap_intersects(my_cpuset, *sock_cpuset)) {
            /* Create a copy to identify list of free coress */
            *free_sock_cpuset = hwloc_bitmap_dup(*sock_cpuset);
            /* Store my sock ID */
            mv2_my_sock_id = i;
            break;
        }
    }
    if (i == num_sockets) {
        mpi_errno = MPI_ERR_OTHER;
        MPIU_ERR_POP(mpi_errno);
    } else {
        /* Remove cores used by processes from list of available cores */
        for (i = 0; i < g_smpi.num_local_nodes; ++i) {
            hwloc_bitmap_clr(*free_sock_cpuset, local_core_ids[i]);
        }
        hwloc_bitmap_snprintf(cpu_str, 128, *free_sock_cpuset);
        PRINT_DEBUG(DEBUG_INIT_verbose, "Free sock_cpuset = %s\n", cpu_str);
    }

    if (my_cpuset) {
        hwloc_bitmap_free(my_cpuset);
    }
fn_fail:
    return mpi_errno;
}

#undef FUNCNAME
#define FUNCNAME smpi_identify_core_for_async_thread
#undef FCNAME
#define FCNAME MPIDI_QUOTE(FUNCNAME)
int smpi_identify_core_for_async_thread(MPIDI_PG_t * pg)
{
    int i = 0;
    int mpi_errno = MPI_SUCCESS;
    hwloc_cpuset_t sock_cpuset = NULL;
    hwloc_cpuset_t free_sock_cpuset = NULL;

    /* Gather IDs of cores local processes are bound to */
    mpi_errno = smpi_identify_allgather_local_core_ids(pg);
    if (mpi_errno) {
        MPIU_ERR_POP(mpi_errno);
    }
    /* Identify my socket and cores available in my socket */
    mpi_errno = smpi_identify_free_cores(&sock_cpuset, &free_sock_cpuset);
    if (mpi_errno) {
        MPIU_ERR_POP(mpi_errno);
    }
    /* Identify core to be used for async thread */
    if (!hwloc_bitmap_iszero(free_sock_cpuset)) {
        for (i = 0; i < g_smpi.num_local_nodes; ++i) {
            /* If local process 'i' is on a core on my socket */
            if (hwloc_bitmap_isset(sock_cpuset, local_core_ids[i])) {
                mv2_my_async_cpu_id = hwloc_bitmap_next(free_sock_cpuset, mv2_my_async_cpu_id);
                if (i == g_smpi.my_local_id) {
                    break;
                }
            }
        }
        /* Ensure async thread gets bound to a core */
        while (mv2_my_async_cpu_id < 0) {
            mv2_my_async_cpu_id = hwloc_bitmap_next(free_sock_cpuset, mv2_my_async_cpu_id);
        }
    }
    PRINT_DEBUG(DEBUG_INIT_verbose>0, "[local_rank: %d]: sock_id = %d, cpu_id = %d, async_cpu_id = %d\n",
                    g_smpi.my_local_id, mv2_my_sock_id, mv2_my_cpu_id, mv2_my_async_cpu_id);

fn_exit:
    /* Free temporary memory */
    if (local_core_ids) {
        MPIU_Free(local_core_ids);
    }
    /* Free cpuset */
    if (sock_cpuset) {
        hwloc_bitmap_free(sock_cpuset);
    }
    if (free_sock_cpuset) {
        hwloc_bitmap_free(free_sock_cpuset);
    }
    return mpi_errno;

fn_fail:
    goto fn_exit;
}
#endif /*defined(CHANNEL_MRAIL)*/

#undef FUNCNAME
#define FUNCNAME smpi_setaffinity
#undef FCNAME
#define FCNAME MPIDI_QUOTE(FUNCNAME)
int smpi_setaffinity(int my_local_id)
{
    int selected_socket = 0;
    int mpi_errno = MPI_SUCCESS;

    hwloc_cpuset_t cpuset;
    MPIDI_STATE_DECL(MPID_STATE_SMPI_SETAFFINITY);
    MPIDI_FUNC_ENTER(MPID_STATE_SMPI_SETAFFINITY);

#if !defined(CHANNEL_MRAIL)
    mv2_hca_aware_process_mapping = 0;
#endif
    mpi_errno = hwloc_topology_init(&topology);
    hwloc_topology_set_flags(topology, HWLOC_TOPOLOGY_FLAG_IO_DEVICES);
    if (mpi_errno != 0) {
        mv2_enable_affinity = 0;
    }

    if (mv2_enable_affinity > 0) {
        long N_CPUs_online = sysconf(_SC_NPROCESSORS_ONLN);

        if (N_CPUs_online < 1) {
            MPIU_ERR_SETFATALANDJUMP2(mpi_errno,
                                      MPI_ERR_OTHER,
                                      "**fail", "%s: %s", "sysconf",
                                      strerror(errno));
        }

        hwloc_topology_load(topology);
        cpuset = hwloc_bitmap_alloc();

        /* Call the cpu_mapping function to find out about how the
         * processors are numbered on the different sockets.
         * The hardware information gathered from this function
         * is required to determine the best set of intra-node thresholds.
         * However, since the user has specified a mapping pattern,
         * we are not going to use any of our proposed binding patterns
         */
        mpi_errno = get_cpu_mapping_hwloc(N_CPUs_online, topology);
        if (mpi_errno != MPI_SUCCESS) {
            /* In case, we get an error from the hwloc mapping function */
            mpi_errno = get_cpu_mapping(N_CPUs_online);
        }

        if (s_cpu_mapping) {
            /* If the user has specified how to map the processes, use it */
            char tp_str[s_cpu_mapping_line_max + 1];

            mpi_errno = mv2_get_assigned_cpu_core(my_local_id, s_cpu_mapping,
                                                    s_cpu_mapping_line_max, tp_str);
            if (mpi_errno != 0) {
                fprintf(stderr, "Error parsing CPU mapping string\n");
                mv2_enable_affinity = 0;
                MPIU_Free(s_cpu_mapping);
                s_cpu_mapping = NULL;
                goto fn_fail;
            }

            // parsing of the string
            char *token = tp_str;
            int cpunum = 0;
            while (*token != '\0') {
                if (isdigit(*token)) {
                    cpunum = first_num_from_str(&token);
                    if (cpunum >= N_CPUs_online) {
                        fprintf(stderr,
                                "Warning! : Core id %d does not exist on this architecture! \n",
                                cpunum);
                        fprintf(stderr, "CPU Affinity is undefined \n");
                        mv2_enable_affinity = 0;
                        MPIU_Free(s_cpu_mapping);
                        goto fn_fail;
                    }
                    hwloc_bitmap_set(cpuset, cpunum);
                    mv2_my_cpu_id = cpunum;
                    PRINT_DEBUG(DEBUG_SHM_verbose>0, "Set mv2_my_cpu_id = %d\n", mv2_my_cpu_id);
                } else if (*token == ',') {
                    token++;
                } else if (*token == '-') {
                    token++;
                    if (!isdigit(*token)) {
                        fprintf(stderr,
                                "Warning! : Core id %c does not exist on this architecture! \n",
                                *token);
                        fprintf(stderr, "CPU Affinity is undefined \n");
                        mv2_enable_affinity = 0;
                        MPIU_Free(s_cpu_mapping);
                        goto fn_fail;
                    } else {
                        int cpuend = first_num_from_str(&token);
                        if (cpuend >= N_CPUs_online || cpuend < cpunum) {
                            fprintf(stderr,
                                    "Warning! : Core id %d does not exist on this architecture! \n",
                                    cpuend);
                            fprintf(stderr, "CPU Affinity is undefined \n");
                            mv2_enable_affinity = 0;
                            MPIU_Free(s_cpu_mapping);
                            goto fn_fail;
                        }
                        int cpuval;
                        for (cpuval = cpunum + 1; cpuval <= cpuend; cpuval++)
                            hwloc_bitmap_set(cpuset, cpuval);
                    }
                } else if (*token != '\0') {
                    fprintf(stderr,
                            "Warning! Error parsing the given CPU mask! \n");
                    fprintf(stderr, "CPU Affinity is undefined \n");
                    mv2_enable_affinity = 0;
                    MPIU_Free(s_cpu_mapping);
                    goto fn_fail;
                }
            }
            // then attachement
            hwloc_set_cpubind(topology, cpuset, 0);

            MPIU_Free(s_cpu_mapping);
            s_cpu_mapping = NULL;
        } else {
            /* The user has not specified how to map the processes,
             * use the data available in /proc/cpuinfo file to decide
             * on the best cpu mapping pattern
             */
            if (mpi_errno != MPI_SUCCESS || custom_cpu_mapping == NULL) {
                /* For some reason, we were not able to retrieve the cpu mapping
                 * information. We are falling back on the linear mapping.
                 * This may not deliver the best performace
                 */
                hwloc_bitmap_only(cpuset, my_local_id % N_CPUs_online);
                mv2_my_cpu_id = (my_local_id % N_CPUs_online);
                PRINT_DEBUG(DEBUG_SHM_verbose>0, "Set mv2_my_cpu_id = %d\n", mv2_my_cpu_id);
                hwloc_set_cpubind(topology, cpuset, 0);
            } else {
                /*
                 * We have all the information that we need. We will bind the
                 * processes to the cpu's now
                 */
                char tp_str[custom_cpu_mapping_line_max + 1];

                mpi_errno = mv2_get_assigned_cpu_core(my_local_id, custom_cpu_mapping,
                        custom_cpu_mapping_line_max, tp_str);
                if (mpi_errno != 0) {
                    fprintf(stderr, "Error parsing CPU mapping string\n");
                    mv2_enable_affinity = 0;
                    goto fn_fail;
                }

                int cores_per_socket = 0;
#if defined(CHANNEL_MRAIL)
                if (!SMP_ONLY && !mv2_user_defined_mapping) {
                    char *value = NULL;
                    if ((value = getenv("MV2_HCA_AWARE_PROCESS_MAPPING")) != NULL) {
                        mv2_hca_aware_process_mapping = !!atoi(value);
                    }
                    if (likely(mv2_hca_aware_process_mapping)) {
                        int num_cpus = hwloc_get_nbobjs_by_type(topology, HWLOC_OBJ_PU);
                        int depth_sockets = hwloc_get_type_depth(topology, HWLOC_OBJ_SOCKET);
                        int num_sockets = hwloc_get_nbobjs_by_depth(topology, depth_sockets);

                        selected_socket = mv2_get_cpu_core_closest_to_hca(my_local_id, num_cpus,
                                num_sockets, depth_sockets);
                        if (selected_socket < 0) {
                            fprintf(stderr, "Error getting closest socket\n");
                            mv2_enable_affinity = 0;
                            goto fn_fail;
                        }
                        cores_per_socket = num_cpus/num_sockets;
                    }
                }
#endif /* defined(CHANNEL_MRAIL) */

                if (mv2_binding_level == LEVEL_CORE) {
                    if (
#if defined(CHANNEL_MRAIL)
                        SMP_ONLY ||
#endif
                        mv2_user_defined_mapping || !mv2_hca_aware_process_mapping
                       )
                    {
                        hwloc_bitmap_only(cpuset, atol(tp_str));
                        mv2_my_cpu_id = atol(tp_str);
                        PRINT_DEBUG(DEBUG_SHM_verbose>0, "Set mv2_my_cpu_id = %d\n", mv2_my_cpu_id);
                    } else {
                        hwloc_bitmap_only(cpuset,
                                (atol(tp_str) % cores_per_socket)
                                + (selected_socket * cores_per_socket));
                        mv2_my_cpu_id = ((atol(tp_str) % cores_per_socket)
                                        + (selected_socket * cores_per_socket));
                        PRINT_DEBUG(DEBUG_SHM_verbose>0, "Set mv2_my_cpu_id = %d\n", mv2_my_cpu_id);
                    }
                } else {
                    if (
#if defined(CHANNEL_MRAIL)
                        SMP_ONLY ||
#endif
                        mv2_user_defined_mapping || !mv2_hca_aware_process_mapping
                        ) {
                        hwloc_bitmap_from_ulong(cpuset, atol(tp_str));
                    } else {
                        hwloc_bitmap_from_ulong(cpuset,
                                (atol(tp_str) % cores_per_socket)
                                + (selected_socket * cores_per_socket));
                    }
                }
                hwloc_set_cpubind(topology, cpuset, 0);
            }

            MPIU_Free(custom_cpu_mapping);
        }
        /* Free cpuset */
        hwloc_bitmap_free(cpuset);
    }

  fn_exit:
    MPIDI_FUNC_EXIT(MPID_STATE_SMPI_SETAFFINITY);
    return mpi_errno;

  fn_fail:
    goto fn_exit;
}

void mv2_show_cpu_affinity(MPIDI_PG_t * pg)
{
    int i, j, num_cpus, my_rank, pg_size;
    int mpi_errno = MPI_SUCCESS, errflag = 0;
    char buf[512];
    cpu_set_t *allproc_cpu_set;
    MPID_Comm *comm_world = NULL;
    MPIDI_VC_t *vc = NULL;

    comm_world = MPIR_Process.comm_world;
    pg_size = comm_world->local_size;
    my_rank = comm_world->rank;

    allproc_cpu_set = (cpu_set_t *) MPIU_Malloc(sizeof(cpu_set_t) * pg_size);
    CPU_ZERO(&allproc_cpu_set[my_rank]);
    sched_getaffinity(0, sizeof(cpu_set_t), &allproc_cpu_set[my_rank]);

    mpi_errno = MPIR_Allgather_impl(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, allproc_cpu_set,
                                    sizeof(cpu_set_t), MPI_BYTE, comm_world, &errflag);
    if (mpi_errno != MPI_SUCCESS) {
        fprintf(stderr, "MPIR_Allgather_impl returned error");
        return;
    }
    if (my_rank == 0) {
        num_cpus = sysconf(_SC_NPROCESSORS_CONF);
        fprintf(stderr, "-------------CPU AFFINITY-------------\n");
        for (i = 0; i < pg_size; i++) {
            MPIDI_Comm_get_vc(comm_world, i, &vc);
            if (vc->smp.local_rank != -1) {
                MPIU_Memset(buf, 0, sizeof(buf));
                for (j = 0; j < num_cpus; j++) {
                    if (CPU_ISSET(j, &allproc_cpu_set[vc->pg_rank])) {
                        sprintf((char *) (buf + strlen(buf)), "%3d", j);
                    }
                }
                fprintf(stderr, "RANK:%d  CPU_SET: %s\n", i, buf);
            }
        }
        fprintf(stderr, "-------------------------------------\n");
    }
    MPIU_Free(allproc_cpu_set);
}

#undef FUNCNAME
#define FUNCNAME MPIDI_CH3I_set_affinity
#undef FCNAME
#define FCNAME MPIDI_QUOTE(FUNCNAME)
int MPIDI_CH3I_set_affinity(MPIDI_PG_t * pg, int pg_rank)
{
    char *value;
    int mpi_errno = MPI_SUCCESS;
    int my_local_id;
    MPIDI_VC_t *vc;

    MPIDI_STATE_DECL(MPID_STATE_MPIDI_CH3I_SET_AFFINITY);
    MPIDI_FUNC_ENTER(MPID_STATE_MPIDI_CH3I_SET_AFFINITY);

    if ((value = getenv("MV2_ENABLE_AFFINITY")) != NULL) {
        mv2_enable_affinity = atoi(value);
    }

    if (mv2_enable_affinity && (value = getenv("MV2_CPU_MAPPING")) != NULL) {
        /* Affinity is on and the user has supplied a cpu mapping string */
        int linelen = strlen(value);
        if (linelen < s_cpu_mapping_line_max) {
            s_cpu_mapping_line_max = linelen;
        }
        s_cpu_mapping =
            (char *) MPIU_Malloc(sizeof(char) * (s_cpu_mapping_line_max + 1));
        strncpy(s_cpu_mapping, value, s_cpu_mapping_line_max);
        s_cpu_mapping[s_cpu_mapping_line_max] = '\0';
        mv2_user_defined_mapping = TRUE;
    }

    if (mv2_enable_affinity && (value = getenv("MV2_CPU_MAPPING")) == NULL) {
        /* Affinity is on and the user has not specified a mapping string */
        if ((value = getenv("MV2_CPU_BINDING_POLICY")) != NULL) {
            /* User has specified a binding policy */
            if (!strcmp(value, "bunch") || !strcmp(value, "BUNCH")) {
                mv2_binding_policy = POLICY_BUNCH;
            } else if (!strcmp(value, "scatter") || !strcmp(value, "SCATTER")) {
                mv2_binding_policy = POLICY_SCATTER;
            } else {
                MPIU_ERR_SETFATALANDJUMP1(mpi_errno, MPI_ERR_OTHER,
                                          "**fail", "**fail %s",
                                          "CPU_BINDING_PRIMITIVE: Policy should be bunch or scatter.");
            }
            mv2_user_defined_mapping = TRUE;
        } else {
            /* User has not specified a binding policy.
             * We are going to do "bunch" binding, by default  */
            mv2_binding_policy = POLICY_BUNCH;
        }
    }

    if (mv2_enable_affinity && (value = getenv("MV2_CPU_MAPPING")) == NULL) {
        /* Affinity is on and the user has not specified a mapping string */
        if ((value = getenv("MV2_CPU_BINDING_LEVEL")) != NULL) {
            /* User has specified a binding level */
            if (!strcmp(value, "core") || !strcmp(value, "CORE")) {
                mv2_binding_level = LEVEL_CORE;
            } else if (!strcmp(value, "socket") || !strcmp(value, "SOCKET")) {
                mv2_binding_level = LEVEL_SOCKET;
            } else if (!strcmp(value, "numanode") || !strcmp(value, "NUMANODE")) {
                mv2_binding_level = LEVEL_NUMANODE;
            } else {
                MPIU_ERR_SETFATALANDJUMP1(mpi_errno, MPI_ERR_OTHER,
                                          "**fail", "**fail %s",
                                          "CPU_BINDING_PRIMITIVE: Level should be core, socket, or numanode.");
            }
            if (MV2_ARCH_INTEL_XEON_PHI_7250 == mv2_get_arch_type() &&
                    mv2_binding_level != LEVEL_CORE) {
                if (MPIDI_Process.my_pg_rank == 0) {
                    fprintf(stderr, "CPU_BINDING_PRIMITIVE: Only core level binding supported for this architecture.\n");
                }
                mpi_errno = MPI_ERR_OTHER;
                goto fn_fail;
            }
            mv2_user_defined_mapping = TRUE;
        } else {
            /* User has not specified a binding level.
             * We are going to do "core" binding, by default  */
            mv2_binding_level = LEVEL_CORE;
        }
    }

    /* Get my VC */
    MPIDI_PG_Get_vc(pg, pg_rank, &vc);
    my_local_id = vc->smp.local_rank;

    if (mv2_enable_affinity) {
        mpi_errno = smpi_setaffinity(my_local_id);
        if (mpi_errno != MPI_SUCCESS) {
            MPIU_ERR_POP(mpi_errno);
        }
    }
  fn_exit:
    MPIDI_FUNC_EXIT(MPID_STATE_MPIDI_CH3I_SET_AFFINITY);
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}
