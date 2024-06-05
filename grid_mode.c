#include "grid_mode.h"

/** Contruct an instance of GridHeight
 *
 * All the dimensions described in the "shape" array have to be greater than 0.
 */
GridHeight *build_grids(int shape[2], int normal_axis, int ngroups,
        const char *grid_fn, char dens) {
    GridHeight *grid_store;
    int grid, i;

    /* Check dimensions */
    if (shape[0] <= 0 || shape[1] <= 0) {
        fprintf(stderr,
                "I can not build a grid with this dimensions: (%d, %d)\n",
                shape[0], shape[1]);
        exit(1);
    }

    /* Allocate empty structure */
    snew(grid_store, 1);
    
    grid_store->nframes = 0;
    grid_store->ngroups = ngroups;
    grid_store->invvol = 0;
    grid_store->dens = dens;
    for (i=0; i<2; ++i) {
        /* Store the shape */
        grid_store->shape[i] = shape[i];
        /* Set the slice width to 0, just in case */
        grid_store->width[i] = 0;
        /* Set box_width sommation to 0 */
        grid_store->box_width[i] = 0.0;
    }
    /* Define the axis */
    grid_store->axis[0] = normal_axis;
    switch (normal_axis) {
        case 0:
            grid_store->axis[1] = 1; grid_store->axis[2] = 2;
            break;
        case 1:
            grid_store->axis[1] = 0; grid_store->axis[2] = 2;
            break;
        case 2:
            grid_store->axis[1] = 0; grid_store->axis[2] = 1;
            break;
        default:
            gmx_fatal(FARGS,"Invalid axes. Terminating. \n");
    }

    /* Allocate the grids */
    snew(grid_store->grids, ngroups);
    for (grid = 0; grid < ngroups; ++grid) {
        (grid_store->grids)[grid] = realMatrix(shape[0], shape[1], 0.0);
    }

    /* Open the files */
    grid_store->out_grid = ffopen(grid_fn, "w");
    if (grid_store->out_grid == NULL) {
        fprintf(stderr, "Error oppenning %s for grid mode\n", grid_fn);
        exit(1);
    }

    return grid_store;
}

/** Clean an instance of GridHeight
 */
void clean_grids(GridHeight *grid_store) {
    int grid = 0;
    if (grid_store) {
        for (grid = 0; grid < grid_store->ngroups; ++grid) {
            deleteRealMat(grid_store->grids[grid], grid_store->shape[0]);
        }
        sfree(grid_store->grids);
        ffclose(grid_store->out_grid);
        sfree(grid_store);
    }
}

void grid_start_frame(GridHeight *grid_store, matrix box) {
    int i = 0;
    int axis = 0;
    if (grid_store) {
        grid_store->nframes += 1;
        for (i=0; i<2; ++i) {
            axis = grid_store->axis[i+1];
            grid_store->width[i] = box[axis][axis]/grid_store->shape[i];
            grid_store->box_width[i] += box[axis][axis];
        }
        grid_store->invvol = (grid_store->shape[0] * grid_store->shape[1])/
            (box[XX][XX] * box[YY][YY] * box[ZZ][ZZ]);
    }
}

void grid_store(GridHeight *grid, int group, rvec atom, t_pbc *pbc, real mass) {
    int slice[2] = {0, 0};
    int i = 0;
    if (grid) {
        put_atom_in_box((real (*)[3])pbc->box,atom);
        for (i =0; i<2; ++i) {
            slice[i] = atom[grid->axis[i+1]]/grid->width[i];
        }
        grid->grids[group][slice[0]][slice[1]] += mass*grid->invvol;
    }
}

void grid_end(GridHeight *grid_store) {
    char labels[] = "XYZ";
    int i, j, group;
    if (grid_store) {
        for (group = 0; group < grid_store->ngroups; ++group) {
            for (i=0; i < grid_store->shape[0]; ++i) {
                for (j=0; j < grid_store->shape[1]; ++j) {
                    grid_store->grids[group][i][j] /= grid_store->nframes;
                    if (grid_store->dens == 'm') {
                        grid_store->grids[group][i][j] *= AMU/(NANO*NANO*NANO);
                    }
                }
            }
        }
        /* Write the output */
        fprintf(grid_store->out_grid, "@xwidth %7.3f\n",
                grid_store->box_width[0]/grid_store->nframes);
        fprintf(grid_store->out_grid, "@ywidth %7.3f\n",
                grid_store->box_width[1]/grid_store->nframes);
        fprintf(grid_store->out_grid, "@xlabel %c (nm)\n",
                labels[grid_store->axis[1]]);
        fprintf(grid_store->out_grid, "@ylabel %c (nm)\n",
                labels[grid_store->axis[2]]);
        fprintf(grid_store->out_grid,
                "@legend Partial mass density (kg/m^3)\n");
        for (group = 0; group < grid_store->ngroups; ++group) {
            for (i=0; i < grid_store->shape[0]; ++i) {
                for (j=0; j < grid_store->shape[1]; ++j) {
                    if (j > 0) {
                        fprintf(grid_store->out_grid, "\t");
                    }
                    fprintf(grid_store->out_grid, "%7.3f",
                            grid_store->grids[group][i][j]);
                }
                fprintf(grid_store->out_grid, "\n");
            }
            fprintf(grid_store->out_grid, "&&\n");
        }
    }
}


