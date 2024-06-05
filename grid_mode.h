#ifndef _grid_mode_h
#define _grid_mode_h

#include <math.h>

#include <gromacs/macros.h>
#include <gromacs/smalloc.h>
#include <gromacs/typedefs.h>
#include <gromacs/gmx_fatal.h>
#include <gromacs/pbc.h>
#include <gromacs/physics.h>
#include <gromacs/futil.h>

#include "matrix.h"

/** Store the height field of each leaflet and the membrane thickness as grids
 *
 * The sampling for each grid is also stored for averaging purposes and to
 * filter low sampling cells.
 *
 * Grids and sampling are in arrays describing, in order, the first and second
 * lealet and the thickness.
 *
 * The shape of the grids is also stored to avoid looking out of boundaries.
 */
typedef struct GridHeight {
    real ***grids;    
    int  shape[2];
    FILE *out_grid;
    real width[2];
    int axis[3];
    real box_width[2];
    int nframes;
    int ngroups;
    real invvol;
    char dens;
} GridHeight;

GridHeight *build_grids(int shape[2], int normal_axis, int ngroups,
        const char *grid_fn, char dens);

void clean_grids(GridHeight *grid_store);

void grid_start_frame(GridHeight *grid_store, matrix box);

void grid_store(GridHeight *grid, int group, rvec atom, t_pbc *pbc, real mass);

void grid_end(GridHeight *grid_store);

#endif /*  _grid_mode_h */
