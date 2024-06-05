#include "dist_mode.h"
/** Calculate the minimum function between an atom and a set of atoms
 *
 * Parameters :
 *  - idx_fix : the index of the reference atom in the table of coordinates (x)
 *  - nmobile : the number of atoms in the set of atoms
 *  - x       : the table of atom coordinates
 *  - index   : the indices of the atom of the set in x
 *  - pbc     : the periodic box
 *  - axis    : the axis to exclude from distance calculation
 * Returns :
 *  - the minimum distance between the reference group and the atom of interest
 */
real _min_distance(int idx_fix, int nmobile, rvec *x, atom_id *index,
        t_pbc *pbc, int axis) {
    real rd;
    rvec reference = {0,0,0}, interest = {0,0,0};
    rvec dx;
    int i, j;
    real rd_min=GMX_REAL_MAX;
    for (i=0; i<DIM; ++i) {
        if (i != axis) {
            interest[i] = x[idx_fix][i];
        }
    }
    for (j=0;(j<nmobile);j++) {
        for (i=0; i<DIM; ++i) {
            if (i != axis) {
                reference[i] = x[index[j]][i];
            }
        }
        pbc_dx(pbc,interest,reference,dx);
        rd=norm(dx);
        if (rd < rd_min) {
            rd_min=rd;
        }
    }
    return rd_min;
}

DistMode *build_dist(int length, int normal_axis, int ngroups, char dens,
        const char *dist_fn, output_env_t oenv,
        const char *index_fn, t_topology *top, const char **legend,
        gmx_bool b3D, gmx_bool bCOM) {
    DistMode *dist_store;
    int prof, i;
    atom_id **index;
    int *isize;
    char **grpnames;

    /* Check dimensions */
    if (length <= 0) {
        fprintf(stderr,
                "I can not build a distance profile with this length: %d\n",
                length);
        exit(1);
    }

    /* Allocate empty structure */
    snew(dist_store, 1);
    
    dist_store->nframes = 0;
    /* Store the shape */
    dist_store->length = length;
    /* Set the slice width to 0, just in case */
    dist_store->width = 0;
    /* Set box_width sommation to 0 */
    dist_store->box_width = 0.0;
    /* Define the axis */
    dist_store->axis[0] = normal_axis;
    if (b3D)
        dist_store->axis[1] = -1;
    else
        dist_store->axis[1] = normal_axis;
    dist_store->ngroups = ngroups;
    dist_store->max_dist = 0;
    dist_store->height = 0;
    dist_store->dens = dens;
    dist_store->b3D = b3D;

    /* Allocate the profiles */
    snew(dist_store->data, ngroups);
    for (prof = 0; prof < ngroups; ++prof) {
        snew(dist_store->data[prof], length);
        for (i=0; i<length; ++i) {
            dist_store->data[prof][i] = 0;
        }
    }

    /* Get the reference group index */
    snew(index, 1);
    snew(isize, 1);
    snew(grpnames, 1);
    printf("Select reference group for distance calcultation:\n");
    get_index(&(top->atoms), index_fn, 1, isize,index,grpnames);
    dist_store->ref_index = index[0];
    dist_store->ref_size = isize[0];

    dist_store->com = NULL;
    dist_store->ref_mass = 0;
    dist_store->bCOM = bCOM;
    if (bCOM) {
        dist_store->ref_mass = get_mass(index[0], isize[0], top);
    }

    /* Open the output file */
    dist_store->out_dist = xvgropen(dist_fn,"Density",
            "Distance from Protein (nm)","Density (kg/m^3)",oenv);
    xvgr_legend(dist_store->out_dist, dist_store->ngroups, legend, oenv);
    return dist_store;
}

void clean_dist(DistMode *dist_store) {
    int prof = 0;
    if (dist_store) {
        for (prof = 0; prof < dist_store->ngroups; ++prof) {
            sfree(dist_store->data[prof]);
        }
        sfree(dist_store->data);
        sfree(dist_store->ref_index);
        fclose(dist_store->out_dist);
        sfree(dist_store);
    }
}

void dist_start_frame(DistMode *dist_store, matrix box, rvec *x,
        t_topology *top) {
    int i = 0;
    real max_dist = INT_MAX;
    if (dist_store) {
        /* Find what the maximum distance is */
        for (i=0; i<DIM; ++i) {
            if ((dist_store->b3D || i != dist_store->axis[0]) 
                    && box[i][i]/2 < max_dist) {
                max_dist = box[i][i]/2;
            }
        }
        dist_store->nframes += 1;
        dist_store->width = max_dist/dist_store->length;
        dist_store->box_width += max_dist;
        dist_store->max_dist = max_dist;
        dist_store->height = box[dist_store->axis[0]][dist_store->axis[0]];
        if (dist_store->bCOM) {
            dist_store->com = center_of_mass(dist_store->ref_index, 
                    dist_store->ref_size, x, top, dist_store->ref_mass);
            make_2D(*dist_store->com, dist_store->axis[1], *dist_store->com);
        }
    }
}

void dist_store(DistMode *dist, int group, int atom, rvec *x, t_pbc *pbc,
        real mass) {
    int slice = 0;
    int i = 0;
    rvec pointA;
    real distance = 0;
    real vslice = 0;
    real r1 = 0, r2 = 0;
    real height = 0;
    if (dist) {
        height = dist->height;
        if (dist->bCOM) {
            make_2D(x[atom], dist->axis[1], pointA);
            distance = get_distance(pointA, *dist->com, pbc);
        }
        else {
            distance = min_dist(x[atom], dist->ref_index, dist->ref_size, x,
                    pbc, dist->axis[1]);
        }
        slice = distance/dist->width;
        if (slice < dist->length) {
            r1 = dist->max_dist * ((float)slice / (float)dist->length);
            r2 = dist->max_dist * ((float)(slice + 1) / (float)dist->length);
            if (dist->b3D)
                vslice = (4.0/3.0) * PI  * (r2*r2*r2 - r1*r1*r1);
            else
                vslice = height * PI * (r2*r2 - r1*r1);
            /*printf("1/vslice: %f, r1: %f, r2: %f, height: %f, vol: %f\n", 1/vslice, r1, r2, height, vslice);*/
            /*printf("mass: %f\n", mass);*/
            /*printf("slice: %d, vslice: %f\n", slice, vslice);*/
            dist->data[group][slice] += mass/vslice;
            /*dist->data[group][slice] += 1;*/
        }
    }
}

void dist_end(DistMode *dist_store) {
    if (dist_store) {
        int i, group;
        real bin_size = 0;
        dist_store->box_width /= dist_store->nframes;
        bin_size = dist_store->box_width/dist_store->length;
        /* Write the output */
        for (i=0; i<dist_store->length; ++i) {
            fprintf(dist_store->out_dist, "%7.3f", i*bin_size);
            for (group=0; group < dist_store->ngroups; ++group) {
                dist_store->data[group][i] /= dist_store->nframes;
                if (dist_store->dens == 'm') {
                    dist_store->data[group][i] *= AMU/(NANO*NANO*NANO);
                }
                fprintf(dist_store->out_dist, "\t%7.3f",
                        dist_store->data[group][i]);
            }
            fprintf(dist_store->out_dist, "\n");
        }
    }
}
