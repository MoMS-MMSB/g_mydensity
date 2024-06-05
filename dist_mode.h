#ifndef _dist_mode_h
#define _dist_mode_h

#include <math.h>

#include <gromacs/statutil.h>
#include <gromacs/macros.h>
#include <gromacs/smalloc.h>
#include <gromacs/typedefs.h>
#include <gromacs/gmx_fatal.h>
#include <gromacs/pbc.h>
#include <gromacs/xvgr.h>
#include <gromacs/vec.h>
#include <gromacs/physics.h>

#include "distances.h"

#define PI (3.141592653589793)

typedef struct DistMode {
    real **data;    
    int  length;
    FILE *out_dist;
    real width;
    int axis[2];
    real box_width;
    int nframes;
    int ngroups;
    char dens;
    atom_id *ref_index;
    int ref_size;
    real max_dist;
    real height;
    gmx_bool b3D;
    gmx_bool bCOM;
    real ref_mass;
    rvec *com;
} DistMode; 

DistMode *build_dist(int length, int normal_axis, int ngroups, char dens,
        const char *dist_fn, output_env_t oenv,
        const char *index_fn, t_topology *top, const char **legend,
        gmx_bool b3D, gmx_bool bCOM);

void clean_dist(DistMode *dist_store);

void dist_start_frame(DistMode *dist_store, matrix box, rvec *x,
        t_topology *top);

void dist_end_frame(DistMode *dist_store, int adt);

void dist_store(DistMode *dist, int group, int atom, rvec *x, t_pbc *pbc,
        real mass);

void dist_end(DistMode *dist_store);

#endif
