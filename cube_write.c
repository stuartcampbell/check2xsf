/* Write a Gaussian cube file, one density only
 *       Units must be Bohr, not Angstoms
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

void cube_write(FILE* outfile, struct grid *g){
  int i,j,k;
  double x,y,z,*dptr1,*dptr2;

  if (!g->data){
    fprintf(stderr,"Cannot write .cube file when no 3D data requested.\n");
    exit(1);
  }

  fprintf(outfile,"DENSITY:\n\n");

  fprintf(outfile,"%d 0.0 0.0 0.0\n",natoms);

  for(i=0;i<3;i++) fprintf(outfile,"%d %f %f %f\n",g->size[i],
                                basis[i][0]/g->size[i]/BOHR,
                                basis[i][1]/g->size[i]/BOHR,
                                basis[i][2]/g->size[i]/BOHR);

/* Need to write coords in Cartesian basis */
  for(i=0;i<natoms;i++){
    x=atoms[i].abs[0]/BOHR;
    y=atoms[i].abs[1]/BOHR;
    z=atoms[i].abs[2]/BOHR;
    fprintf(outfile,"%d 0.0 %f %f %f\n",atoms[i].atno,x,y,z);
  }

  dptr2=g->data;
  for(k=0;k<g->size[0];k++){
    for(j=0;j<g->size[1];j++){
      dptr1=dptr2+((k*g->size[1])+j)*g->size[2];
      for(i=0;i<g->size[2];i++)
        fprintf(outfile,"%f\n",*(dptr1+i));
    }
  }
}
