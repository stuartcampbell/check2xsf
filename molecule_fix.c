/* If passed a non-null pointer, assume is pointer to three integers
 * Shift everything by these numbers of FFT grid cells
 *
 * If pointer is null, attempt to work out what shift would be required
 * in x,y and z independently.
 *
 * Assumes all 3D data are on same FFT grid
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

void molecule_fix(int* m_abc){
  int i,j,k,ii,jj,kk,off1,off2;
  int fft[3],shift[3];
  double min,max,*grid3;
  struct grid *gptr;

  if (grid1.data){
    fft[0]=grid1.size[0];
    fft[1]=grid1.size[1];
    fft[2]=grid1.size[2];
  }
  else fft[0]=fft[1]=fft[2]=1000;


  if(!m_abc){ /* We are expected to determine the shift automatically */

    for(k=0;k<3;k++){
      min=1;
      max=-1;
      for(i=0;i<natoms;i++){
        if (atoms[i].frac[k]>max)
          max=atoms[i].frac[k];
        if (atoms[i].frac[k]<min)
          min=atoms[i].frac[k];
      }
      if (((max-min)<=1)&&((min<0)||(max>1))){  /* We should shift */
        if (!m_abc){
          m_abc=malloc(3*sizeof(int));
          if (!m_abc) error_exit("Malloc error for three ints!");
          for(kk=0;kk<k;kk++) m_abc[kk]=0;
        }
        if ((max<=0.5)&&(min>=-0.5))
          m_abc[k]=fft[k]/2;
        else
          m_abc[k]=(0.5-0.5*(max+min))*fft[k]+0.99;
      }
    }
  }

  if(!m_abc){
    if (debug>1) fprintf(stderr,"Shift requested, but none found\n");
    return;
  }

  if (debug>1)
      fprintf(stderr,"Will translate grid by (%d,%d,%d) gridpoints\n"
                     " which is (%f,%f,%f)\n",
              m_abc[0],m_abc[1],m_abc[2],(double)m_abc[0]/fft[0],
              (double)m_abc[1]/fft[1],(double)m_abc[2]/fft[2]);

  /* Ions are easy */

  for(i=0;i<natoms;i++)
    for(k=0;k<3;k++)
      atoms[i].frac[k]+=(double)m_abc[k]/fft[k];

  addabs();

    /* Grids are harder */

  gptr=&grid1;
  while((gptr)&&(gptr->data)){

    grid3=malloc(sizeof(double)*gptr->size[0]*gptr->size[1]*gptr->size[2]);
    if (!grid3) error_exit("Memory allocation error in grid shift");
    for(k=0;k<3;k++){
      shift[k]=(gptr->size[k]-(m_abc[k]*gptr->size[k])/fft[k])%gptr->size[k];
      while(shift[k]<0) shift[k]+=gptr->size[k];
    }

    for(k=0;k<gptr->size[0];k++){
      kk=(k+shift[0])%gptr->size[0];
      for(j=0;j<gptr->size[1];j++){
        jj=(j+shift[1])%gptr->size[1];
        off1=(k*gptr->size[1]+j)*gptr->size[2];
        off2=(kk*gptr->size[1]+jj)*gptr->size[2];
        for(i=0;i<gptr->size[2];i++){
          ii=(i+shift[2])%gptr->size[2];
          *(grid3+off1+i)=*(gptr->data+off2+ii);
        }
      }
    }
    free(gptr->data);
    gptr->data=grid3;
    gptr=gptr->next;
  }
}
