

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

#include <stdio.h>

#include "c2xsf.h"
#include "c2xsf_extern.h"

void print_globals(int level){
  int i,j;

  if (level==0) return;

  if (level>=1){
    fprintf(stderr,"Cell volume %f\n",cell_vol);
    fprintf(stderr,"natoms      %d\n",natoms);
    fprintf(stderr,"Basis set\n");
    for(i=0;i<=2;i++)
      fprintf(stderr,"%f %f %f\n",basis[i][0],basis[i][1],basis[i][2]);
    fprintf(stderr,"First FFT grid     %d %d %d\n",
                      grid1.size[0],grid1.size[1],grid1.size[2]);
    fprintf(stderr,"spins=%d\n",spins);
    if (level>=2){
      fprintf(stderr,"verbosity=%d\n",debug);
      fprintf(stderr,"Grid data %s\n",grid1.data?"exists":"absent");
      fprintf(stderr,"Ionic positions, fractional\n");
      for(i=0;i<natoms;i++)
        fprintf(stderr,"%d %f %f %f\n",atoms[i].atno,atoms[i].frac[0],
                       atoms[i].frac[1],atoms[i].frac[2]);
      if (level>=3){
        if(forces){
          fprintf(stderr,"Ionic forces\n");
          for(i=0;i<natoms;i++)
            fprintf(stderr,"%d %f %f %f\n",atoms[i].atno,
                  atoms[i].force[0],atoms[i].force[1],atoms[i].force[2]);
        }
        if (level>=4){
          fprintf(stderr,"Ionic positions, absolute\n");
          for(i=0;i<natoms;i++)
            fprintf(stderr,"%d %f %f %f\n",atoms[i].atno,atoms[i].abs[0],
                       atoms[i].abs[1],atoms[i].abs[2]);
          fprintf(stderr,"Reciprocal basis set\n");
          for(i=0;i<=2;i++)
            fprintf(stderr,"%f %f %f\n",recip[i][0],recip[i][1],recip[i][2]);
        }
      }
    }
  }
}
