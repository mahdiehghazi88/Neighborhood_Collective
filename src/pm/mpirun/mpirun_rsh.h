#ifndef _MPIRUN_RSH_H
#define _MPIRUN_RSH_H 1
/*
 * Copyright (C) 1999-2001 The Regents of the University of California
 * (through E.O. Lawrence Berkeley National Laboratory), subject to
 * approval by the U.S. Department of Energy.
 *
 * Use of this software is under license. The license agreement is included
 * in the file MVICH_LICENSE.TXT.
 *
 * Developed at Berkeley Lab as part of MVICH.
 *
 * Authors: Bill Saphir      <wcsaphir@lbl.gov>
 *          Michael Welcome  <mlwelcome@lbl.gov>
 */

/* Copyright (c) 2001-2014, The Ohio State University. All rights
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

#define BASE_ENV_LEN    17
#define COMMAND_LEN     20000

#define RSH_CMD         "/usr/bin/rsh"
#define SSH_CMD         "/cvmfs/soft.computecanada.ca/custom/bin/ssh"
#define ENV_CMD         "/cvmfs/soft.computecanada.ca/nix/var/nix/profiles/16.09/bin/env"
#define DBG_CMD         "/cvmfs/soft.computecanada.ca/nix/var/nix/profiles/16.09/bin/gdb"
#define XTERM_CMD       "/cvmfs/soft.computecanada.ca/nix/var/nix/profiles/16.09/bin/xterm"
#define SHELL_CMD       "/cvmfs/soft.computecanada.ca/nix/var/nix/profiles/16.09/bin/bash"
#define TOTALVIEW_CMD   "/usr/totalview/bin/totalview"

#define SSH_ARG         "-q"
#define SHELL_ARG       "-c"

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <stdio.h>
#include <stdarg.h>
#include <unistd.h>
#ifdef MAC_OSX
#include <sys/wait.h>
#else
#include <wait.h>
#endif
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <fcntl.h>
#include <getopt.h>
#include <netdb.h>
#include <assert.h>
#include <libgen.h>
#include "mpirun_util.h"
#include <mpichconf.h>

#define PRINT_MVAPICH2_VERSION() printf("Version: mvapich2-" MVAPICH2_VERSION "\n")

#define PMGR_VERSION PMGR_COLLECTIVE

/*This list is used for dpm to take care of the mpirun_rsh started.*/
typedef struct list_pid_mpirun {
    pid_t pid;
    struct list_pid_mpirun *next;
} list_pid_mpirun_t;

extern int NSPAWNS;

#define RUNNING(i) ((plist[i].state == P_STARTED ||                 \
            plist[i].state == P_CONNECTED ||                        \
            plist[i].state == P_RUNNING) ? 1 : 0)

/* other information: a.out and rank are implicit. */

#define SEPARATOR ':'

#ifndef PARAM_GLOBAL
#define PARAM_GLOBAL "/etc/mvapich.conf"
#endif

#endif

int handle_spawn_req(int readsock);
void mpispawn_checkin(int);

/* vi:set sw=4 sts=4 tw=80: */
