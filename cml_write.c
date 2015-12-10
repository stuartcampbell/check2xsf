/* Someone will want CML output, even if I cannot imagine why.
 * This might amuse him.
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
#include<stdlib.h> /* calloc */

#include "c2xsf.h"
#include "c2xsf_extern.h"

#define MAX_ELS 103

void cml_write(FILE* outfile){
  int i,*n_in_el;
  double abc[6];

  n_in_el=calloc(MAX_ELS+1,sizeof(int));
  if (!n_in_el) error_exit("Calloc error in cml_write");

  fprintf(outfile,"<?xml version=\"1.0\"?>\n"
         "<molecule xmlns=\"http://www.xml-cml.org/schema\">\n");

  cart2abc(basis,abc,1);

  fprintf(outfile,
      "<crystal xmlns:cmldict=\"http://www.xml-cml.org/dict/cmlDict\">\n");
  fprintf(outfile,"  <scalar title=\"a\">%f</scalar>\n",abc[0]);
  fprintf(outfile,"  <scalar title=\"b\">%f</scalar>\n",abc[1]);
  fprintf(outfile,"  <scalar title=\"c\">%f</scalar>\n",abc[2]);
  fprintf(outfile,"  <scalar title=\"alpha\">%f</scalar>\n",abc[3]);
  fprintf(outfile,"  <scalar title=\"beta\">%f</scalar>\n",abc[4]);
  fprintf(outfile,"  <scalar title=\"gamma\">%f</scalar>\n",abc[5]);
  fprintf(outfile,"</crystal>\n");

  fprintf(outfile,"<atomArray>\n");
  for(i=0;i<natoms;i++)
      fprintf(outfile,"  <atom id=\"%s%d\" xFract=\"%f\" yFract=\"%f\""
              " zFract=\"%f\" elementType=\"%s\"/>\n",atno2sym(atoms[i].atno),
              ++n_in_el[min(atoms[i].atno,MAX_ELS)],
              atoms[i].frac[0],atoms[i].frac[1],atoms[i].frac[2],
              atno2sym(atoms[i].atno));
  fprintf(outfile,"</atomArray>\n");
  fprintf(outfile,"</molecule>\n");

  free(n_in_el);
}
