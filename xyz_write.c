/* Write an xyz format file */

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

/* Each frame (one here) has the format of
 * NATOMS
 * Comment line (blank here)
 * Atomic symbol x y z (repeated natoms times)
 */

#include<stdio.h>
#include<stdlib.h>

#include "c2xsf.h"
#include "c2xsf_extern.h"

void xyz_write(FILE* outfile, struct grid *g){
  int i;

  fprintf(outfile,"%d\n\n",natoms);

  for(i=0;i<natoms;i++)
    fprintf(outfile,"%3s %8f %8f %8f\n",
                     atno2sym(atoms[i].atno),atoms[i].abs[0],
                     atoms[i].abs[1],atoms[i].abs[2]);
}
