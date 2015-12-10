/* Write an xsf file (XCrySDen). */

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

/* The format is:

CRYSTAL
PRIMVEC
[lattice_vector_1_x] [y] [z]
[lattice_vector_2_x] [y] [z]
[lattice_vector_3_x] [y] [z]
CONVVEC
[lattice_vector_1_x] [y] [z]
[lattice_vector_2_x] [y] [z]
[lattice_vector_3_x] [y] [z]
PRIMCOORD
[number_of_atoms] 1
[atom_symbol] [pos_x] [pos_y] [pos_z]
....
BEGIN_BLOCK_DATAGRID_3D
Densities
BEGIN_DATAGRID_3D_chden
[grid_points_vec1+1] [grid_points_vec2+1] [grid_points_vec3+1]
[offset_x] [offset_y] [offset_z] 
[grid_vector_1_x] [y] [z]
[grid_vector_2_x] [y] [z]
[grid_vector_3_x] [y] [z]
[data x=0 y=0 z=0]
[data x=1 y=0 z=0]
...
[data x=grid_points_vec1 y=0 z=0]
[data x=0 y=0 z=0]
[data x=0 y=1 z=0]
...
END_DATAGRID_3D
BEGIN_DATAGRID_3D_spin (optional, as
END_DATAGRID_3D         for chden above)
END_BLOCK_DATAGRID_3D

For added confusion, and to convert from Fortran to C ordering, this
code outputs the grid vectors with x and z exchanged, and thus exchanges
x and z when writing the data...
*/

#include<stdio.h>
#include<stdlib.h>

#include "c2xsf.h"
#include "c2xsf_extern.h"

void xsf_write(FILE* outfile){
  int i,j,k;
  double x,y,z,*dptr1,*dptr2;
  struct grid *gptr;

  if (molecule==0){
    fprintf(outfile,"CRYSTAL\nPRIMVEC\n");
    for(i=0;i<=2;i++) fprintf(outfile,"%f %f %f\n",basis[i][0],
                              basis[i][1],basis[i][2]);
    fprintf(outfile,"CONVVEC\n");
    for(i=0;i<=2;i++) fprintf(outfile,"%f %f %f\n",basis[i][0],
                              basis[i][1],basis[i][2]);
    fprintf(outfile,"PRIMCOORD\n");
    fprintf(outfile,"%d 1\n",natoms);
  }else{
    fprintf(outfile,"MOLECULE\nATOMS\n");
  }

/* Need to write coords in Cartesian basis */
  for(i=0;i<natoms;i++){
    x=atoms[i].abs[0];
    y=atoms[i].abs[1];
    z=atoms[i].abs[2];
    if (forces)
      fprintf(outfile,"%d %f %f %f %f %f %f\n",atoms[i].atno,x,y,z,
              atoms[i].force[0],atoms[i].force[1],atoms[i].force[2]);
    else
      fprintf(outfile,"%d %f %f %f\n",atoms[i].atno,x,y,z);
  }

  gptr=&grid1;
  if((gptr)&&(gptr->data)){
    fprintf(outfile,"BEGIN_BLOCK_DATAGRID_3D\n"
                    "Densities\n");
    while((gptr)&&(gptr->data)){
      if (debug>1) fprintf(stderr,"Writing %s\n",gptr->name);
      fprintf(outfile,"BEGIN_DATAGRID_3D_%s\n",gptr->name);
      fprintf(outfile,"%d %d %d\n",gptr->size[2]+1,
                      gptr->size[1]+1,gptr->size[0]+1);
      fprintf(outfile,"0.0 0.0 0.0\n");
      for(i=2;i>=0;i--) fprintf(outfile,"%f %f %f\n",basis[i][0],
                                  basis[i][1],basis[i][2]);

      dptr2=gptr->data;
      for(k=0;k<=gptr->size[0];k++){
        for(j=0;j<=gptr->size[1];j++){
          dptr1=dptr2+
           (((k%gptr->size[0])*gptr->size[1])+(j%gptr->size[1]))*gptr->size[2];
          for(i=0;i<=gptr->size[2];i++){
            fprintf(outfile,"%f\n",*(dptr1+i%gptr->size[2]));
          }
        }
      }
      fprintf(outfile,"END_DATAGRID_3D\n");
      gptr=gptr->next;
    }

    fprintf(outfile,"END_BLOCK_DATAGRID_3D\n");
  }
}
