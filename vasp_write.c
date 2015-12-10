/* Write a VASP-style output file, one density only
 */

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

#include "c2xsf.h"
#include "c2xsf_extern.h"

void vasp_write(FILE* outfile, struct grid *g){
  int i,j,k;
  double *dptr1;
  int nspec,*natomsp,*spatno;

  fprintf(outfile,"Kyptonite\n1.0000\n");  /* No idea */

  for(i=0;i<3;i++) fprintf(outfile," %f %f %f\n",
                                basis[i][0],
                                basis[i][1],
                                basis[i][2]);

  /* Now we need to know the number of species.
     It must be fewer than the number of atoms...
     This is horribly inefficient, but I intend to
     use it for ethene only... */

  nspec=0;
  natomsp=malloc(natoms*sizeof(int));
  if (!natomsp) error_exit("Malloc error in vasp_write");
  spatno=malloc(natoms*sizeof(int));
  if (!spatno) error_exit("Malloc error in vasp_write");


  for(i=0;i<natoms;i++){
    for(j=0;j<nspec;j++) if (atoms[i].atno==spatno[j]) break;
    if (j==nspec){  /* new species */
      spatno[j]=atoms[i].atno;
      natomsp[j]=1;
      nspec++;
    }else{          /* existing species */
      natomsp[j]++;
    }
  }

  /* Must not have trailing space on this line... */
  for(i=0;i<nspec;i++) fprintf(outfile," %d",natomsp[i]);
  fprintf(outfile,"\n");

  fprintf(outfile,"Direct\n");  /* No idea */

  /* Write out atoms, sorted by species */
  for(i=0;i<nspec;i++)
    for(j=0;j<natoms;j++)
      if (atoms[j].atno==spatno[i])
        fprintf(outfile," %f %f %f\n",atoms[j].frac[0],
                atoms[j].frac[1],atoms[j].frac[2]);

  /* And now write density */

  if (g->data==0) return;

  fprintf(outfile,"\n%d %d %d\n",g->size[0],g->size[1],g->size[2]);

  for(k=0;k<g->size[2];k++){
    for(j=0;j<g->size[1];j++){
      dptr1=g->data+k+j*g->size[2];
      for(i=0;i<g->size[0];i++)
        fprintf(outfile,"%f\n",*(dptr1+i*g->size[1]*g->size[2]));
    }
  }
}
