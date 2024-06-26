/*
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Green Red Orange Magenta Azure Cyan Skyblue
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <math.h>
#include <ctype.h>

#include "string.h"
#include <gromacs/string2.h>

#include <gromacs/sysstuff.h>
#include <gromacs/typedefs.h>
#include <gromacs/smalloc.h>
#include <gromacs/macros.h>
#include <gromacs/gstat.h>
#include <gromacs/vec.h>
#include <gromacs/xvgr.h>
#include <gromacs/pbc.h>
#include <gromacs/copyrite.h>
#include <gromacs/futil.h>
#include <gromacs/statutil.h>
#include <gromacs/index.h>
#include <gromacs/tpxio.h>
#include <gromacs/physics.h>
#include <gromacs/gmx_ana.h>

#include "grid_mode.h"
#include "dist_mode.h"

typedef struct {
  char *atomname;
  int nr_el;
} t_electron;

/****************************************************************************/
/* This program calculates the partial density across the box.              */
/* Peter Tieleman, Mei 1995                                                 */
/****************************************************************************/

/* used for sorting the list */
int compare(void *a, void *b)
{
  t_electron *tmp1,*tmp2;
  tmp1 = (t_electron *)a; tmp2 = (t_electron *)b;

  return strcmp(tmp1->atomname,tmp2->atomname);
}

int get_electrons(t_electron **eltab, const char *fn)
{
  char buffer[256];  /* to read in a line   */
  char tempname[80]; /* buffer to hold name */
  int tempnr; 

  FILE *in;
  int nr;            /* number of atomstypes to read */
  int i;

  if ( !(in = ffopen(fn,"r")))
    gmx_fatal(FARGS,"Couldn't open %s. Exiting.\n",fn);

  if(NULL==fgets(buffer, 255, in))
  {
      gmx_fatal(FARGS,"Error reading from file %s",fn);
  }
 
  if (sscanf(buffer, "%d", &nr) != 1)
    gmx_fatal(FARGS,"Invalid number of atomtypes in datafile\n");

  snew(*eltab,nr);

  for (i=0;i<nr;i++) {
    if (fgets(buffer, 255, in) == NULL)
      gmx_fatal(FARGS,"reading datafile. Check your datafile.\n");
    if (sscanf(buffer, "%s = %d", tempname, &tempnr) != 2)
      gmx_fatal(FARGS,"Invalid line in datafile at line %d\n",i+1);
    (*eltab)[i].nr_el = tempnr;
    (*eltab)[i].atomname = strdup(tempname);
  }
  ffclose(in);
  
  /* sort the list */
  fprintf(stderr,"Sorting list..\n");
  qsort ((void*)*eltab, nr, sizeof(t_electron), 
	 (int(*)(const void*, const void*))compare);

  return nr;
}

void center_coords(t_atoms *atoms,matrix box,rvec x0[],int axis)
{
  int  i,m;
  real tmass,mm;
  rvec com,shift,box_center;
  
  tmass = 0;
  clear_rvec(com);
  for(i=0; (i<atoms->nr); i++) {
    mm     = atoms->atom[i].m;
    tmass += mm;
    for(m=0; (m<DIM); m++) 
      com[m] += mm*x0[i][m];
  }
  for(m=0; (m<DIM); m++) 
    com[m] /= tmass;
  calc_box_center(ecenterDEF,box,box_center);
  rvec_sub(box_center,com,shift);
  shift[axis] -= box_center[axis];
  
  for(i=0; (i<atoms->nr); i++) 
    rvec_dec(x0[i],shift);
}

void calc_electron_density(const char *fn, atom_id **index, int gnx[], 
			   real ***slDensity, int *nslices, t_topology *top,
			   int ePBC,
			   int axis, int nr_grps, real *slWidth, 
			   t_electron eltab[], int nr,gmx_bool bCenter,
                           const output_env_t oenv)
{
  rvec *x0;              /* coordinates without pbc */
  matrix box;            /* box (3x3) */
  double invvol;
  int natoms;            /* nr. atoms in trj */
  t_trxstatus *status;  
  int i,n,               /* loop indices */
      nr_frames = 0,     /* number of frames */
      slice;             /* current slice */
  t_electron *found;     /* found by bsearch */
  t_electron sought;     /* thingie thought by bsearch */
  gmx_rmpbc_t  gpbc=NULL;
  t_pbc *pbc = NULL;
 
  real t, 
        z;

  if (axis < 0 || axis >= DIM) {
    gmx_fatal(FARGS,"Invalid axes. Terminating\n");
  }

  if ((natoms = read_first_x(oenv,&status,fn,&t,&x0,box)) == 0)
    gmx_fatal(FARGS,"Could not read coordinates from statusfile\n");
  
  if (! *nslices)
    *nslices = (int)(box[axis][axis] * 10); /* default value */
  fprintf(stderr,"\nDividing the box in %d slices\n",*nslices);

  snew(*slDensity, nr_grps);
  for (i = 0; i < nr_grps; i++)
    snew((*slDensity)[i], *nslices);
 
  if (ePBC != epbcNONE)
      snew(pbc,1);
  else
      pbc = NULL;
  gpbc = gmx_rmpbc_init(&top->idef,ePBC,top->atoms.nr,box);
  /*********** Start processing trajectory ***********/
  do {
      if (pbc) {
          set_pbc(pbc,ePBC,box);
          /* make molecules whole again */
          gmx_rmpbc(gpbc,natoms,box,x0);
      }

    if (bCenter)
      center_coords(&top->atoms,box,x0,axis);
    
    *slWidth = box[axis][axis]/(*nslices);
    invvol = *nslices/(box[XX][XX]*box[YY][YY]*box[ZZ][ZZ]);

    for (n = 0; n < nr_grps; n++) {      
      for (i = 0; i < gnx[n]; i++) {   /* loop over all atoms in index file */
	  z = x0[index[n][i]][axis];
	  while (z < 0) 
	    z += box[axis][axis];
	  while (z > box[axis][axis])
	    z -= box[axis][axis];
      
	  /* determine which slice atom is in */
	  slice = (z / (*slWidth)); 
	  sought.nr_el = 0;
	  sought.atomname = strdup(*(top->atoms.atomname[index[n][i]]));

	  /* now find the number of electrons. This is not efficient. */
	  found = (t_electron *)
	    bsearch((const void *)&sought,
		    (const void *)eltab, nr, sizeof(t_electron), 
		    (int(*)(const void*, const void*))compare);

	  if (found == NULL)
	    fprintf(stderr,"Couldn't find %s. Add it to the .dat file\n",
		    *(top->atoms.atomname[index[n][i]]));
	  else  
	    (*slDensity)[n][slice] += (found->nr_el - 
				       top->atoms.atom[index[n][i]].q)*invvol;
	  free(sought.atomname);
	}
    }
      nr_frames++;
  } while (read_next_x(oenv,status,&t,natoms,x0,box));
  gmx_rmpbc_done(gpbc);

  /*********** done with status file **********/
  close_trj(status);
  
/* slDensity now contains the total number of electrons per slice, summed 
   over all frames. Now divide by nr_frames and volume of slice 
*/

  fprintf(stderr,"\nRead %d frames from trajectory. Counting electrons\n",
	  nr_frames);

  for (n =0; n < nr_grps; n++) {
    for (i = 0; i < *nslices; i++)
      (*slDensity)[n][i] /= nr_frames;
  }

  sfree(x0);  /* free memory used by coordinate array */
}

void calc_density(const char *fn, atom_id **index, int gnx[], 
		  real ***slDensity, int *nslices, t_topology *top, int ePBC,
		  int axis, int nr_grps, real *slWidth, gmx_bool bCenter,
                  const output_env_t oenv, GridHeight *grid, DistMode *dist)
{
  rvec *x0;              /* coordinates without pbc */
  matrix box;            /* box (3x3) */
  double invvol;
  int natoms;            /* nr. atoms in trj */
  t_trxstatus *status;  
  int  **slCount,         /* nr. of atoms in one slice for a group */
      i,j,n,               /* loop indices */
      teller = 0,      
      ax1=0, ax2=0,
      nr_frames = 0,     /* number of frames */
      slice;             /* current slice */
  real t, 
        z;
  char *buf;             /* for tmp. keeping atomname */
  gmx_rmpbc_t  gpbc=NULL;

  t_pbc *pbc;

  if (axis < 0 || axis >= DIM) {
    gmx_fatal(FARGS,"Invalid axes. Terminating\n");
  }

  if ((natoms = read_first_x(oenv,&status,fn,&t,&x0,box)) == 0)
    gmx_fatal(FARGS,"Could not read coordinates from statusfile\n");
  
  if (! *nslices) {
    *nslices = (int)(box[axis][axis] * 10); /* default value */
    fprintf(stderr,"\nDividing the box in %d slices\n",*nslices);
  }
  
  snew(*slDensity, nr_grps);
  for (i = 0; i < nr_grps; i++)
    snew((*slDensity)[i], *nslices);

  if (ePBC != epbcNONE)
      snew(pbc,1);
  else
      pbc = NULL;

  gpbc = gmx_rmpbc_init(&top->idef,ePBC,top->atoms.nr,box);
  /*********** Start processing trajectory ***********/
  do {
      if (pbc) {
          set_pbc(pbc,ePBC,box);
          /* make molecules whole again */
          gmx_rmpbc(gpbc,natoms,box,x0);
      }

    if (bCenter)
      center_coords(&top->atoms,box,x0,axis);
   
    grid_start_frame(grid, box);
    dist_start_frame(dist, box, x0, top);

    *slWidth = box[axis][axis]/(*nslices);
    invvol = *nslices/(box[XX][XX]*box[YY][YY]*box[ZZ][ZZ]);
    teller++;
    
    for (n = 0; n < nr_grps; n++) {      
        for (i = 0; i < gnx[n]; i++) {   /* loop over all atoms in index file */
            grid_store(grid, n, x0[index[n][i]], pbc,
                    top->atoms.atom[index[n][i]].m);
            dist_store(dist, n, index[n][i], x0, pbc,
                    top->atoms.atom[index[n][i]].m);
            z = x0[index[n][i]][axis];
            while (z < 0) 
                z += box[axis][axis];
            while (z > box[axis][axis])
                z -= box[axis][axis];

            /* determine which slice atom is in */
            slice = (int)(z / (*slWidth)); 
            (*slDensity)[n][slice] += top->atoms.atom[index[n][i]].m*invvol;
            /*printf("invvol: %f, vol: %f\n", invvol, 1/invvol);*/
            /*printf("mass prof: %f\n", top->atoms.atom[index[n][i]].m);*/
        }
    }
    nr_frames++;
  } while (read_next_x(oenv,status,&t,natoms,x0,box));
  gmx_rmpbc_done(gpbc);

  grid_end(grid);
  dist_end(dist);

  /*********** done with status file **********/
  close_trj(status);
  
  /* slDensity now contains the total mass per slice, summed over all
     frames. Now divide by nr_frames and volume of slice 
     */
  
  fprintf(stderr,"\nRead %d frames from trajectory. Calculating density\n",
	  nr_frames);

  for (n =0; n < nr_grps; n++) {
    for (i = 0; i < *nslices; i++) {
      (*slDensity)[n][i] /= nr_frames;
    }
  }

  sfree(x0);  /* free memory used by coordinate array */
}

void plot_density(real *slDensity[], const char *afile, int nslices,
		  int nr_grps, char *grpname[], real slWidth, 
		  const char **dens_opt,
		  gmx_bool bSymmetrize, const output_env_t oenv)
{
  FILE  *den;
  const char *ylabel=NULL;
  int   slice, n;
  real  ddd;

  switch (dens_opt[0][0]) {
  case 'm': ylabel = "Density (kg m\\S-3\\N)"; break;
  case 'n': ylabel = "Number density (nm\\S-3\\N)"; break;
  case 'c': ylabel = "Charge density (e nm\\S-3\\N)"; break;
  case 'e': ylabel = "Electron density (e nm\\S-3\\N)"; break;
  }
  
  den = xvgropen(afile, "Partial densities", "Box (nm)", ylabel,oenv);

  xvgr_legend(den,nr_grps,(const char**)grpname,oenv);

  for (slice = 0; (slice < nslices); slice++) { 
    fprintf(den,"%12g  ", slice * slWidth);
    for (n = 0; (n < nr_grps); n++) {
      if (bSymmetrize)
	ddd = (slDensity[n][slice]+slDensity[n][nslices-slice-1])*0.5;
      else
	ddd = slDensity[n][slice];
      if (dens_opt[0][0] == 'm')
	fprintf(den,"   %12g", ddd*AMU/(NANO*NANO*NANO));
      else
	fprintf(den,"   %12g", ddd);
    }
    fprintf(den,"\n");
  }

  ffclose(den);
}
 
int gmx_mydensity(int argc,char *argv[])
{
  const char *desc[] = {
    "Compute partial densities across the box, using an index file.[PAR]",
    "For the total density of NPT simulations, use [TT]g_energy[tt] instead.",
    "[PAR]",
    "Densities are in kg/m^3, and number densities or electron densities can also be",
    "calculated. For electron densities, a file describing the number of",
    "electrons for each type of atom should be provided using [TT]-ei[tt].",
    "It should look like:[BR]",
    "   [TT]2[tt][BR]",
    "   [TT]atomname = nrelectrons[tt][BR]",
    "   [TT]atomname = nrelectrons[tt][BR]",
    "The first line contains the number of lines to read from the file.",
    "There should be one line for each unique atom name in your system.",
    "The number of electrons for each atom is modified by its atomic",
    "partial charge.",
    "[PAR]",
    "WARNING: This is a modified version of g_density. It allows to calculate partial density landscapes on a grid (using the [TT]-og[tt] option) and partial density profile as a function of the distance from a group (using the [TT]-od[tt] option). In the latter case, distances are calculated in the plane normal to the axis given with the [TT]-d[tt] option. To get the distances in 3D, use the [TT]-3d[tt] option."
  };

  output_env_t oenv;
  static const char *dens_opt[] = 
    { NULL, "mass", "number", "charge", "electron", NULL };
  static int  axis = 2;          /* normal to memb. default z  */
  static const char *axtitle="Z"; 
  static int  nslices = 50;      /* nr of slices defined       */
  static int  nslices2 = -1;      /* nr of slices defined       */
  static int  ngrps   = 1;       /* nr. of groups              */
  static gmx_bool bSymmetrize=FALSE;
  static gmx_bool bCenter=FALSE;
  static gmx_bool b3D=TRUE;
  static gmx_bool bCOM=FALSE;
  t_pargs pa[] = {
    { "-d", FALSE, etSTR, {&axtitle}, 
      "Take the normal on the membrane in direction X, Y or Z." },
    { "-sl",  FALSE, etINT, {&nslices},
      "Divide the box in #nr slices." },
    { "-sl2",  FALSE, etINT, {&nslices2},
      "Divide the box second dimension in #nr slices." },
    { "-dens",    FALSE, etENUM, {dens_opt},
      "Density"},
    { "-ng",       FALSE, etINT, {&ngrps},
      "Number of groups to compute densities of" },
    { "-symm",    FALSE, etBOOL, {&bSymmetrize},
      "Symmetrize the density along the axis, with respect to the center. Useful for bilayers." },
    { "-center",  FALSE, etBOOL, {&bCenter},
      "Shift the center of mass along the axis to zero. This means if your axis is Z and your box is bX, bY, bZ, the center of mass will be at bX/2, bY/2, 0."},
    { "-3d",  FALSE, etBOOL, {&b3D}, "Calculate distance in 3D instead of 2D."},
    /*
     * { "-com",  FALSE, etBOOL, {&bCOM},
     *   "Use distance to the center of mass instead of minimum distance."}
     */
  };

  const char *bugs[] = {
    "When calculating electron densities, atomnames are used instead of types. This is bad.",
  };
  
  real **density;      /* density per slice          */
  real slWidth;          /* width of one slice         */
  char **grpname;        /* groupnames                 */
  int  nr_electrons;     /* nr. electrons              */
  int  *ngx;             /* sizes of groups            */
  t_electron *el_tab;    /* tabel with nr. of electrons*/
  t_topology *top;       /* topology 		       */ 
  int  ePBC;
  atom_id   **index;     /* indices for all groups     */
  int  i;

  GridHeight *grid_store = NULL;
  DistMode *dist_store = NULL;

  t_filenm  fnm[] = {    /* files for g_density 	  */
    { efTRX, "-f", NULL,  ffREAD },  
    { efNDX, NULL, NULL,  ffOPTRD }, 
    { efTPX, NULL, NULL,  ffREAD },    	    
    { efDAT, "-ei", "electrons", ffOPTRD }, /* file with nr. of electrons */
    { efXVG,"-o","density",ffWRITE }, 	    
    { efDAT,"-og","density_grid",ffOPTWR }, 	    
    { efDAT,"-od","density_dist",ffOPTWR }, 	    
  };
  
#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);

  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME | PCA_BE_NICE,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,asize(bugs),bugs,
                    &oenv);

  if (bSymmetrize && !bCenter) {
    fprintf(stderr,"Can not symmetrize without centering. Turning on -center\n");
    bCenter = TRUE;
  }
  /* Calculate axis */
  axis = toupper(axtitle[0]) - 'X';
  
  top = read_top(ftp2fn(efTPX,NFILE,fnm),&ePBC);     /* read topology file */
  if (dens_opt[0][0] == 'n') {
    for(i=0; (i<top->atoms.nr); i++)
      top->atoms.atom[i].m = 1;  
  } else if (dens_opt[0][0] == 'c') {
    for(i=0; (i<top->atoms.nr); i++)
      top->atoms.atom[i].m = top->atoms.atom[i].q;  
  }

  snew(grpname,ngrps);
  snew(index,ngrps);
  snew(ngx,ngrps);
 
  get_index(&top->atoms,ftp2fn_null(efNDX,NFILE,fnm),ngrps,ngx,index,grpname); 

  if (dens_opt[0][0] == 'e') {
    nr_electrons =  get_electrons(&el_tab,ftp2fn(efDAT,NFILE,fnm));
    fprintf(stderr,"Read %d atomtypes from datafile\n", nr_electrons);

    calc_electron_density(ftp2fn(efTRX,NFILE,fnm),index, ngx, &density, 
			  &nslices, top, ePBC, axis, ngrps, &slWidth, el_tab, 
			  nr_electrons,bCenter,oenv);
  } else {
    if (opt2fn_null("-og",NFILE,fnm)) {
        if (nslices2 <= 0) {
            nslices2 = nslices;
        }
        grid_store = build_grids((int[2]){nslices, nslices2}, axis, ngrps,
                opt2fn("-og",NFILE,fnm), dens_opt[0][0]);
    }
    if (opt2bSet("-od", NFILE, fnm)) {
        dist_store = build_dist(nslices, axis, ngrps, dens_opt[0][0],
                opt2fn("-od",NFILE,fnm), oenv, ftp2fn(efNDX,NFILE,fnm), top,
                (const char **)grpname, b3D, bCOM);
    }
    calc_density(ftp2fn(efTRX,NFILE,fnm),index, ngx, &density, &nslices, top, 
		 ePBC, axis, ngrps, &slWidth, bCenter,oenv, grid_store, dist_store); 
	clean_grids(grid_store);
	clean_dist(dist_store);
  }
  
  plot_density(density, opt2fn("-o",NFILE,fnm),
	       nslices, ngrps, grpname, slWidth, dens_opt,
	       bSymmetrize,oenv);
  
  do_view(oenv,opt2fn("-o",NFILE,fnm), "-nxy");       /* view xvgr file */
  thanx(stderr);
  return 0;
}
    
int
main(int argc, char *argv[])
{
  gmx_mydensity(argc,argv);
  return 0;
}


