/* Write a PDB file */


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

#include "c2xsf.h"
#include "c2xsf_extern.h"

void cell_write(FILE* outfile){
  int i;

  fprintf(outfile,"%%block LATTICE_CART\nang\n");

  for(i=0;i<3;i++)
    fprintf(outfile,"%f %f %f\n",
                  basis[i][0],basis[i][1],basis[i][2]);
  fprintf(outfile,"%%endblock LATTICE_CART\n\n");

  fprintf(outfile,"%%block POSITIONS_FRAC\n");
  for(i=0;i<natoms;i++)
    fprintf(outfile,"%3s %14.9f %14.9f %14.9f\n",
                     atno2sym(atoms[i].atno),atoms[i].frac[0],
                     atoms[i].frac[1],atoms[i].frac[2]);
  fprintf(outfile,"%%endblock POSITIONS_FRAC\n");
}


void cell_write_abc(FILE* outfile){
  int i;
  double abc[6];

  fprintf(outfile,"%%block LATTICE_ABC\nang\n");
  cart2abc(basis,abc,1);
  fprintf(outfile,"%f %f %f\n%f %f %f\n",
                  abc[0],abc[1],abc[2],abc[3],abc[4],abc[5]);
  fprintf(outfile,"%%endblock LATTICE_ABC\n\n");

  fprintf(outfile,"%%block POSITIONS_FRAC\n");
  for(i=0;i<natoms;i++)
    fprintf(outfile,"%3s %14.9f %14.9f %14.9f\n",
                     atno2sym(atoms[i].atno),atoms[i].frac[0],
                     atoms[i].frac[1],atoms[i].frac[2]);
  fprintf(outfile,"%%endblock POSITIONS_FRAC\n");
}

void cell_write_abs(FILE* outfile){
  int i;

  fprintf(outfile,"%%block LATTICE_CART\nang\n");

  for(i=0;i<3;i++)
    fprintf(outfile,"%f %f %f\n",
                  basis[i][0],basis[i][1],basis[i][2]);
  fprintf(outfile,"%%endblock LATTICE_CART\n\n");

  fprintf(outfile,"%%block POSITIONS_ABS\n");
  for(i=0;i<natoms;i++)
    fprintf(outfile,"%3s %14.9f %14.9f %14.9f\n",
                     atno2sym(atoms[i].atno),atoms[i].abs[0],
                     atoms[i].abs[1],atoms[i].abs[2]);
  fprintf(outfile,"%%endblock POSITIONS_ABS\n");
}


void cell_write_abc_abs(FILE* outfile){
  int i;
  double abc[6];

  fprintf(outfile,"%%block LATTICE_ABC\nang\n");
  cart2abc(basis,abc,1);
  fprintf(outfile,"%f %f %f\n%f %f %f\n",
                  abc[0],abc[1],abc[2],abc[3],abc[4],abc[5]);
  fprintf(outfile,"%%endblock LATTICE_ABC\n\n");

  fprintf(outfile,"%%block POSITIONS_ABS\n");
  for(i=0;i<natoms;i++)
    fprintf(outfile,"%3s %14.9f %14.9f %14.9f\n",
                     atno2sym(atoms[i].atno),atoms[i].abs[0],
                     atoms[i].abs[1],atoms[i].abs[2]);
  fprintf(outfile,"%%endblock POSITIONS_ABS\n");
}

