/* Function is passed two vectors. It applies the rotation
 * which moves first to second to the basis set, returning
 * that in new_cell. For obfuscation, the input vectors are
 * also passed in new_cell.
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

#include <stdio.h>
#include <math.h>

#include "c2xsf.h"
#include "c2xsf_extern.h"

void rotation(double new_cell[3][3]){
  double v1[3],v2[3],e[3],angle,tmp;
  double rot_mat[3][3];
  int i,j,k;
  
  /* Extract two vectors in absolute co-ords*/
  for(i=0;i<3;i++){
    v1[i]=new_cell[0][i];
    v2[i]=new_cell[1][i];
  }

  /* euler is v1xv2 */

  e[0]= v1[1]*v2[2]-v1[2]*v2[1];
  e[1]=-v1[0]*v2[2]+v1[2]*v2[0];
  e[2]= v1[0]*v2[1]-v1[1]*v2[0];

  /* Which needs normalising */

  tmp=sqrt(e[0]*e[0]+e[1]*e[1]+e[2]*e[2]);

  if (tmp<1e-20) error_exit("Impossibly small cross product in rotation");

  e[0]/=tmp;
  e[1]/=tmp;
  e[2]/=tmp;

  /* angle is arccos(v1.v2 / mod(v1).mod(v2)) */
  angle=acos((v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])/
            sqrt((v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2]) *
                 (v2[0]*v2[0]+v2[1]*v2[1]+v2[2]*v2[2])));

  if (debug>1) fprintf(stderr,"Rotating by %g degrees about (%g,%g,%g)\n",
                        angle*180/M_PI,e[0],e[1],e[2]);

  /* Wikipedia says that the rotation matrix is:
   *
   *  I cos(theta) + (1-cos(theta))ee^t - E sin(theta)
   *
   *                 ( 0  -e3   e2 )
   * where       E = ( e3  0   -e1 )
   *                 ( -e2 e1   0  )
   *
   * with e=(e1,e2,e3) being the Euler vector (axis of rotation)
   */

  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      rot_mat[i][j]=0;

  for(i=0;i<3;i++) rot_mat[i][i]=cos(angle);

  tmp=1-cos(angle);

  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      rot_mat[i][j]+=tmp*e[i]*e[j];

  tmp=sin(angle);

  rot_mat[0][1]-=tmp*e[2];
  rot_mat[0][2]+=tmp*e[1];
  rot_mat[1][0]+=tmp*e[2];
  rot_mat[1][2]-=tmp*e[0];
  rot_mat[2][0]-=tmp*e[1];
  rot_mat[2][1]+=tmp*e[0];

  if (debug>2){
    fprintf(stderr,"Rotation matrix\n");
    for(i=0;i<=2;i++)
      fprintf(stderr,"%f %f %f\n",
              rot_mat[i][0],rot_mat[i][1],rot_mat[i][2]);
  }

  /* Now apply rotation to basis */

  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      new_cell[i][j]=0;

  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      for(k=0;k<3;k++)
        new_cell[i][j]+=rot_mat[j][k]*basis[i][k];

  if (debug>1){
    fprintf(stderr,"New basis set\n");
    for(i=0;i<=2;i++)
      fprintf(stderr,"%f %f %f\n",
              new_cell[i][0],new_cell[i][1],new_cell[i][2]);
  }

  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      basis[i][j]=new_cell[i][j];

  /* Everything remains as was in relative co-ordinates, but the
   * absolute co-ords of the atoms need recalculating, as does the
   * reciprocal basis set
   */

  addabs();
  real2rec();

}
