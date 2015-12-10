/* Conversion to arbitrary, specified supercells, with grid interpolation */

/* Copyright (c) 2007 MJ Rutter 
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License version 2
 * as published by the Free Software Foundation.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA  02110-1301, USA.
 */ 

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include "c2xsf.h"
#include "c2xsf_extern.h"

void to235(int *i);
void grid_interp(struct grid *grid,int old_fft[3],
                    double old_basis[3][3],double old_recip[3][3]);

#ifdef QSORT
/* Function to sort atoms by type, third co-ord, then second, then first */
int atom_sort(const void *a, const void *b){ /* This ghastly declaration */
  struct atom *a1 = (struct atom *) a;       /* stops compilers moaning */
  struct atom *a2 = (struct atom *) b;       /* about qsort's prototype */

  if (a1->atno < a2->atno) return(1);
  if (a1->atno > a2->atno) return(-1);
  if (a1->frac[2] < a2->frac[2]-1e-6) return(-1);
  if (a1->frac[2] > a2->frac[2]+1e-6) return(1);
  if (a1->frac[1] < a2->frac[1]-1e-6) return(-1);
  if (a1->frac[1] > a2->frac[1]+1e-6) return(1);
  if (a1->frac[0] < a2->frac[0]-1e-6) return(-1);
  if (a1->frac[0] > a2->frac[0]+1e-6) return(1);
  return(0); /* All co-ords equal! */
}
#endif

void super(double new_basis[3][3], int rhs){
  int i,j,k,l,na,at,old_fft[3];
  double old_basis[3][3],old_recip[3][3],old_vol,dtmp;
  double new_in_old[3][3],abc[6];
  struct atom *old_atoms;
  struct grid *gptr;
  int old_natoms;
  int scan_min[3],scan_max[3];
  double corner, fscan_min[3],fscan_max[3],disp[3],new_abs[3],new_frac[3];
  double fft_res;

  /* Save old cell */

  for(i=0;i<3;i++){
    for(j=0;j<3;j++){
      old_basis[i][j]=basis[i][j];
      old_recip[i][j]=recip[i][j];
    }
  }

  old_atoms=atoms;
  old_natoms=natoms;
  old_vol=cell_vol;

  /* Make new cell */

  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      basis[i][j]=new_basis[i][j];

  real2rec(); /* This will also update cell_vol */

  if (cell_vol<0){
    if (rhs){ /* Wrong hand set. Reverse two... */
      for(i=0;i<3;i++){
        dtmp=basis[1][i];
        basis[1][i]=basis[2][i];
        basis[2][i]=dtmp;
      }
      real2rec();
    }
  } else cell_vol=fabs(cell_vol);

  if (debug) fprintf(stderr,"New cell volume %f (%g times old)\n",
               cell_vol,cell_vol/old_vol);
  if (debug>1){
      fprintf(stderr,"New basis set\n");
      for(i=0;i<=2;i++)
        fprintf(stderr,"%f %f %f\n",basis[i][0],basis[i][1],basis[i][2]);
  }

  if (debug>2){
    fprintf(stderr,"New reciprocal basis set\n");
      for(i=0;i<=2;i++)
        fprintf(stderr,"%f %f %f\n",recip[i][0],recip[i][1],recip[i][2]);
  }

  /* Worry about atoms */

  natoms=old_natoms*cell_vol/old_vol+0.5;
  if (fabs(natoms-old_natoms*cell_vol/old_vol)>1e-6){
    fprintf(stderr,"Impossible cell transformation leads to %f atoms\n",
            old_natoms*cell_vol/old_vol);
    exit(1);
  }
  
  if(!(atoms=malloc(natoms*sizeof(struct atom))))
    error_exit("Malloc error in super");

  /* Find extent of new corners in old cell */

  for(i=0;i<3;i++) fscan_min[i]=fscan_max[i]=0;

  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      new_in_old[i][j]=basis[i][0]*old_recip[j][0]+
                       basis[i][1]*old_recip[j][1]+
                       basis[i][2]*old_recip[j][2];

  for(i=0;i<2;i++)
    for(j=0;j<2;j++)
      for(k=0;k<2;k++)
        for(l=0;l<3;l++){
          corner=i*new_in_old[0][l]+j*new_in_old[1][l]+k*new_in_old[2][l];
          if (corner<fscan_min[l]) fscan_min[l]=corner;
          if (corner>fscan_max[l]) fscan_max[l]=corner;
        }

  for(i=0;i<3;i++){
    scan_max[i]=fscan_max[i]+0.01;
    scan_min[i]=fscan_min[i]-0.99;
    if (debug>2) fprintf(stderr,"i=%d scan_min=%d scan_max=%d\n",i,
                         scan_min[i],scan_max[i]);
  }

  /* Now loop over required array of cells */

  na=0;
  for(i=scan_min[0];i<=scan_max[0];i++){
    for(j=scan_min[1];j<=scan_max[1];j++){
      for(k=scan_min[2];k<=scan_max[2];k++){
        for(l=0;l<3;l++)
          disp[l]=i*old_basis[0][l]+j*old_basis[1][l]+k*old_basis[2][l];
        for(at=0;at<old_natoms;at++){
          for(l=0;l<3;l++) new_abs[l]=old_atoms[at].abs[l]+disp[l];
  /* reduce this atom to new basis set */
          for(l=0;l<3;l++){
            new_frac[l]=new_abs[0]*recip[l][0]+
                        new_abs[1]*recip[l][1]+
                        new_abs[2]*recip[l][2];
            new_frac[l]=fmod(new_frac[l]+15,1.0);
          }

  /* add to our list if we have not yet seen it */
          for(l=0;l<na;l++)
            if(((atoms[l].frac[0]-new_frac[0])*(atoms[l].frac[0]-new_frac[0])+
                (atoms[l].frac[1]-new_frac[1])*(atoms[l].frac[1]-new_frac[1])+
                (atoms[l].frac[2]-new_frac[2])*(atoms[l].frac[2]-new_frac[2]))
                <1e-8) break;
          if (l==na){
            if (na>=natoms) error_exit("Too many new atoms found in super.c");
            for(l=0;l<3;l++){
              atoms[na].frac[l]=new_frac[l];
              atoms[na].abs[l]=new_frac[0]*basis[0][l]+
                               new_frac[1]*basis[1][l]+
                               new_frac[2]*basis[2][l];
            }           
            atoms[na].atno=old_atoms[at].atno;
            na++;
          }
        }
      }
    }
  }

#ifdef QSORT
  qsort(atoms,(size_t)na,sizeof(struct atom),atom_sort);
#endif

  if (debug>1) fprintf(stderr,"New cell: na=%d, natoms=%d\n",na,natoms);

  /* Worry about grids */

  gptr=&grid1;

  while((gptr)&&(gptr->data)){
    fft_res=0;
    for(i=0;i<3;i++){
      dtmp=gptr->size[i]/sqrt(old_basis[i][0]*old_basis[i][0]+
                              old_basis[i][1]*old_basis[i][1]+
                              old_basis[i][2]*old_basis[i][2]);
      if (dtmp>fft_res) fft_res=dtmp;
    }
    cart2abc(basis,abc,0);
    for(i=0;i<3;i++){
      old_fft[i]=gptr->size[i];
      gptr->size[i]=abc[i]*fft_res;
      to235(gptr->size+i);
    }
    if (debug) fprintf(stderr,"New FFT grid is %d %d %d\n",
                       gptr->size[0],gptr->size[1],gptr->size[2]);
    grid_interp(gptr,old_fft,old_basis,old_recip);
    gptr=gptr->next;
  }
}

void to235(int *i){
/* Force i to next highest int whose only factors are 2, 3 and 5 */
  int tmp;

  if ((*i)<1) *i=2;
  tmp=*i;
  while (tmp%2==0) tmp/=2;
  while (tmp%3==0) tmp/=3;
  while (tmp%5==0) tmp/=5;

  if (tmp==1) return;

  (*i)++;
  to235(i);
}

void grid_interp(struct grid *grid,int old_fft[3],
                    double old_basis[3][3],double old_recip[3][3]){
  int i,j,k,l,ii,jj,kk,pt[3][2];
  double pabs[3],pfrac[3],ifrac[3][2],*gnew;
  double di,dj,dk,x,sum;

  sum=0;

  if (!(gnew=malloc(sizeof(double)*grid->size[0]*grid->size[1]*grid->size[2])))
     error_exit("Malloc error in grid_interp");

  for(i=0;i<grid->size[0];i++){
    di=(double)i/grid->size[0];
    for(j=0;j<grid->size[1];j++){
      dj=(double)j/grid->size[1];
      for(k=0;k<grid->size[2];k++){
        dk=(double)k/grid->size[2];
        for(l=0;l<3;l++)
          pabs[l]=di*basis[0][l]+dj*basis[1][l]+dk*basis[2][l];
        /* convert to old basis */
        for(l=0;l<3;l++){
          pfrac[l]=pabs[0]*old_recip[l][0]+pabs[1]*old_recip[l][1]+
                   pabs[2]*old_recip[l][2];
        /* reduce to unit cell */
          pfrac[l]=fmod(pfrac[l]+15,1.0);
        /* convert to old grid */
          pfrac[l]=pfrac[l]*old_fft[l];
        /* and to int and fractional parts */
          pt[l][0]=pfrac[l];
          ifrac[l][1]=pfrac[l]-pt[l][0];
        /* and the other side of the grid "cube" */
          pt[l][1]=(pt[l][0]+1)%old_fft[l];
          ifrac[l][0]=1-ifrac[l][1];
        }
        /* Trilinear interpolation */
#define OX(x,y,z) grid->data[((x)*old_fft[1]+(y))*old_fft[2]+(z)]
        x=0;
        for(ii=0;ii<2;ii++)
          for(jj=0;jj<2;jj++)
            for(kk=0;kk<2;kk++)
              x+=ifrac[0][ii]*ifrac[1][jj]*ifrac[2][kk]*
                   OX(pt[0][ii],pt[1][jj],pt[2][kk]);
        gnew[((i*grid->size[1])+j)*grid->size[2]+k]=x;
        sum+=x;
      }
    }
  }

  free(grid->data);
  grid->data=gnew;

  if (debug>1) fprintf(stderr,"On new grid sum=%g int=%g\n",sum,
                     sum*cell_vol/(grid->size[0]*grid->size[1]*grid->size[2]));
}
