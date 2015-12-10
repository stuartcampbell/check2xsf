/* Convert CASTEP .check file to XCrySDen .xsf file */


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
#include <stdlib.h>
#include <unistd.h>  /* isatty */
#include <math.h> /* fabs */
#include <string.h>

#include "c2xsf.h"

void help(void);

void xsf_write(FILE* outfile);
void cube_write(FILE* outfile, struct grid *g);
void xplor_write(FILE* outfile, struct grid *g);
void xplor_fudge(struct grid *g);
void pdb_write(FILE* outfile);
void cell_write(FILE* outfile);
void cell_write_abc(FILE* outfile);
void cell_write_abs(FILE* outfile);
void cell_write_abc_abs(FILE* outfile);
void dx_write(FILE* outfile, struct grid *g);
void vasp_write(FILE* outfile, struct grid *g);
void xyz_write(FILE* outfile, struct grid *g);
void cml_write(FILE* outfile);
void fdf_write(FILE* fdf, char* filename, struct grid *g);

void cell_read(FILE* infile);
void check_read(FILE* infile);
void pdb_read(FILE* infile);

void molecule_fix(int* m_abc);
void super(double new_basis[3][3], int rhs);
void rotation(double new_basis[3][3]);

/* Global variables for system description */

double (*basis)[3],recip[3][3]; /* Basis sets */
double cell_vol,dtmp;
int natoms,spins;
int debug,molecule,format,forces,flags;
char *band_range="-";
char *kpt_range="1";
char *spin_range="-";
struct atom *atoms;
struct grid grid1, *grid_next;

int main(int argc, char **argv)
{
  int i,j,k,l,opt=1,expand,rotate,half_shift=0;
  char *optp,*file1=NULL,*file2=NULL;
  FILE *infile,*outfile;
  double abc[6],new_cell[3][3],new_cell_rel[3][3];
  int *m_abc;
  struct grid *gptr;

/* Initialise all pointers to NULL, etc. */

  atoms=NULL;
  basis=NULL;
  m_abc=NULL;
  natoms=forces=0;
  for(i=0;i<3;i++){
    grid1.size[i]=0;
    for(j=0;j<3;j++) recip[i][j]=0.0;
  }
  spins=1;
  grid_next=&grid1;
  grid1.next=NULL;
  grid1.data=NULL;
  flags=0;

  opt=1;
  optp=argv[opt];
  debug=molecule=format=expand=rotate=0;

  while (opt<argc){
    switch(*optp){
    case 0:
      opt++;
      optp=argv[opt];
      break;
    case '-':
      if(*(optp+1)=='-'){ /* Gnu-style "--" option */
        if (!strcmp(optp,"--xsf")) format=XSF;
        else if (!strcmp(optp,"--cube")) format=CUBE;
        else if (!strcmp(optp,"--xplor")) format=XPLOR;
        else if (!strcmp(optp,"--pdb")) format=PDB;
        else if (!strcmp(optp,"--pdbn")) {format=PDB; flags|=PDBN;}
        else if (!strcmp(optp,"--cell")) format=CELL;
        else if (!strcmp(optp,"--cell_abc")) format=CELL_ABC;
        else if (!strcmp(optp,"--cell_abs")) format=CELL_ABS;
        else if (!strcmp(optp,"--cell_abc_abs")) format=CELL_ABC_ABS;
        else if (!strcmp(optp,"--dx")) format=DX;
        else if (!strcmp(optp,"--vasp")) format=VASP;
        else if (!strcmp(optp,"--xyz")) format=XYZ;
        else if (!strcmp(optp,"--cml")) format=CML;
        else if (!strcmp(optp,"--fdf")) format=FDF;
        else if (!strcmp(optp,"--null")) format=CNULL;
        else if (!strcmp(optp,"--help")) help();
       else {
          fprintf(stderr,"Invalid option %s.\n%s -h for usage.\n",
                   optp,argv[0]);
          exit(1);
        }
        opt++;
        optp=argv[opt]; 
        break;
      }
      while(*(++optp)){
        switch(*optp){
        case 'v':
          debug++;
          break;
        case 'h':
          help();
          break;
        case 'V':
          printf("check2xsf version " C2XSF_VER "\n");
          exit(0);
          break;
        case 'a':
          rotate=1;
          break;
        case 'A':
          flags|=ACCUMULATE;
          break;
        case 'b':
          flags|=BANDS;
          if(*(optp+1)=='='){
            band_range=optp+2;
            while(*((++optp)+1));
          }
          break;
        case 'B':
          flags|=BANDDEN;
          if(*(optp+1)=='='){
            band_range=optp+2;
            while(*((++optp)+1));
          }
          break;
        case 'k':
          if(*(optp+1)=='='){
            kpt_range=optp+2;
            while(*((++optp)+1));
          }
          else error_exit("malformed option -k");
          break;
        case 'c':
          flags|=CHDEN;
          break;
        case 'H':
          half_shift=1;
          break;
        case 'i':
          flags|=BANDIMAG;
          break;
        case 'm':
          if(*(optp+1)=='='){
            m_abc=malloc(3*sizeof(int));
            if (!m_abc) error_exit("malloc error for three ints!");
            if (sscanf(optp+2,"%d,%d,%d",m_abc,m_abc+1,m_abc+2)!=3)
              error_exit("malformed option -m=");
            while(*((++optp)+1));
          }
          molecule=1;
          break;
        case 'p':
          flags|=BANDPHASE;
          break;
        case 'r':
          flags|=BANDREAL;
          break;
        case 'R':
          flags|=RAW;
          break;
        case 's':
          flags|=SPINDEN;
          break;
        case 'S':
          if(*(optp+1)=='='){
            spin_range=optp+2;
            while(*((++optp)+1));
          }
          else error_exit("malformed option -S\n");
          break;
        case 't':
          if(*(optp+1)=='='){
            if (sscanf(optp+2,"(%lf,%lf,%lf)(%lf,%lf,%lf)",
               &new_cell_rel[0][0],&new_cell_rel[0][1],&new_cell_rel[0][2],
               &new_cell_rel[1][0],&new_cell_rel[1][1],&new_cell_rel[1][2])!=6)
               error_exit("malformed option -t=");
             while(*((++optp)+1));
             expand=5;
             break;
          }
          else error_exit("malformed option -t=");
        case 'T':
          if(*(optp+1)=='='){
            if (sscanf(optp+2,"(%lf,%lf,%lf)(%lf,%lf,%lf)",
                       &new_cell[0][0],&new_cell[0][1],&new_cell[0][2],
                       &new_cell[1][0],&new_cell[1][1],&new_cell[1][2])!=6)
               error_exit("malformed option -T=");
             while(*((++optp)+1));
             expand=4;
             break;
          }
          else error_exit("malformed option -T=");
        case 'x':
          if(*(optp+1)=='='){
            if (sscanf(optp+2,"(%lf,%lf,%lf)(%lf,%lf,%lf)(%lf,%lf,%lf)",
               &new_cell_rel[0][0],&new_cell_rel[0][1],&new_cell_rel[0][2],
               &new_cell_rel[1][0],&new_cell_rel[1][1],&new_cell_rel[1][2],
               &new_cell_rel[2][0],&new_cell_rel[2][1],&new_cell_rel[2][2])!=9)
               error_exit("malformed option -x=");
            while(*((++optp)+1));
            expand=2;
            break;
          }
          expand=1;
          break;
        case 'X':
          if(*(optp+1)=='='){
            if (sscanf(optp+2,"(%lf,%lf,%lf)(%lf,%lf,%lf)(%lf,%lf,%lf)",
                       &new_cell[0][0],&new_cell[0][1],&new_cell[0][2],
                       &new_cell[1][0],&new_cell[1][1],&new_cell[1][2],
                       &new_cell[2][0],&new_cell[2][1],&new_cell[2][2])!=9)
               error_exit("malformed option -X=");
             while(*((++optp)+1));
             expand=3;
             break;
          }
          else error_exit("malformed option -X=");
        default:
          fprintf(stderr,"Invalid option %c.\n%s -h for usage.\n",
                   *optp,argv[0]);
           exit(1);
         }
       }
       opt++;
       optp=argv[opt];
       break;
     default:
       if (!file1) file1=argv[opt];
       else if (!file2) file2=argv[opt];
       else{
         fprintf(stderr,"Unexpected argument %s\n%s -h for usage.\n",
                 argv[opt],argv[0]);
         exit(1);
       }
       opt++;
       optp=argv[opt];
       break;
     }
  }

  if (!file1) error_exit("no input file specified.");

  if ((!file2)&&(isatty(fileno(stdout)))&&(format!=CNULL))
    error_exit("refusing to output to a terminal");

  infile=fopen(file1,"rb");
  if(!infile){
    fprintf(stderr,"Error, unable to open %s for reading.\n",file1);
    exit(1);
  }

  if (file2){
    outfile=fopen(file2,"wb");
    if(!outfile){
      fprintf(stderr,"Error, unable to open %s for writing.\n",file2);
      exit(1);
    }
  }
  else
    outfile=stdout;

  i=strlen(file1);
  if((i>4)&&(!strcmp(file1+i-4,".pdb")))
    pdb_read(infile);
  else{
    i=fgetc(infile);
    rewind(infile);
    if ((i==0)||(i==30)||(i==10))  /* three possible first bytes of a */
      check_read(infile);          /* .check or .castep_bin file */
    else
      cell_read(infile);
  }

  if (cell_vol<0) cell_vol=fabs(cell_vol);

  /* Check that we have most of what we need */

  if (!atoms) error_exit("no atoms found!");

  if (debug>=1){
    if (debug>1){
      fprintf(stderr,"Basis set\n");
      for(i=0;i<=2;i++)
        fprintf(stderr,"%f %f %f\n",basis[i][0],basis[i][1],basis[i][2]);
    }
    if (debug>2){
      fprintf(stderr,"Reciprocal basis set\n");
      for(i=0;i<=2;i++)
        fprintf(stderr,"%f %f %f\n",recip[i][0],recip[i][1],recip[i][2]);
    }
    fprintf(stderr,"Cell volume %f\n",cell_vol);
    fprintf(stderr,"natoms      %d\n",natoms);
    if (debug>2){
      fprintf(stderr,"Ionic positions, fractional\n");
      for(i=0;i<natoms;i++)
        fprintf(stderr,"%d %f %f %f\n",atoms[i].atno,atoms[i].frac[0],
                       atoms[i].frac[1],atoms[i].frac[2]);
      fprintf(stderr,"Ionic positions, absolute\n");
      for(i=0;i<natoms;i++)
        fprintf(stderr,"%d %f %f %f\n",atoms[i].atno,atoms[i].abs[0],
                       atoms[i].abs[1],atoms[i].abs[2]);
      if(forces){
        fprintf(stderr,"Ionic forces\n");
        for(i=0;i<natoms;i++)
        fprintf(stderr,"%d %f %f %f\n",atoms[i].atno,
                  atoms[i].force[0],atoms[i].force[1],atoms[i].force[2]);
      }
    }
    fprintf(stderr,"First FFT grid     %d %d %d\n",
                      grid1.size[0],grid1.size[1],grid1.size[2]);
    fprintf(stderr,"spins=%d\n",spins);
    gptr=&grid1;
    while((gptr)&&(gptr->data)){
      fprintf(stderr,"Found 3D data for %s\n",gptr->name);
      if (debug>1){
        double min,max,sum;
        min=1e20;
        max=-1e20;
        sum=0;
        for(i=0;i<gptr->size[0]*gptr->size[1]*gptr->size[2];i++){
          sum=sum+gptr->data[i];
          if (gptr->data[i]>max) max=gptr->data[i];
          if (gptr->data[i]<min) min=gptr->data[i];
        }
        fprintf(stderr,"  min=%g  max=%g  sum=%g  int=%g\n",min,max,
                sum,sum*cell_vol/(gptr->size[0]*gptr->size[1]*gptr->size[2]));
      }
      gptr=gptr->next;
    }
  }

  if (molecule) molecule_fix(m_abc);

  if (half_shift) xplor_fudge(&grid1);

  if (expand){
    switch(expand){
    case 5: /* First vector given in relative terms */
      for(i=0;i<3;i++){
        new_cell[0][i]=new_cell_rel[0][0]*basis[0][i]+
                       new_cell_rel[0][1]*basis[1][i]+
                       new_cell_rel[0][2]*basis[2][i];
        new_cell[1][i]=new_cell_rel[1][i];
      }
      rotation(new_cell);
      break;
    case 4:
      rotation(new_cell);
      break;
    case 3:  /* New basis given explicitly in absolute terms */
      cart2abc(new_cell,abc,0);
      super(new_cell,0);
      break;
    case 2:  /* New basis given explicitly in relative terms */
      for(i=0;i<3;i++)
        for(j=0;j<3;j++){
          new_cell[i][j]=0;
          for(k=0;k<3;k++)
            new_cell[i][j]+=new_cell_rel[i][k]*basis[k][j];
        }
      cart2abc(new_cell,abc,0);
      super(new_cell,0);
      break;
    case 1:  /* We are expected to guess the appropriate new basis */
      cart2abc(basis,abc,0);
      if (debug>1) fprintf(stderr,"LatticeABC=%f %f %f %f %f %f\n",
                       abc[0],abc[1],abc[2],abc[3],abc[4],abc[5]);

/* We know how to convert a=b=c, alpha=beta=gamma=60, to cubic */

      if ((fabs(abc[0]-abc[1])<1e-4*abc[0])&&
          (fabs(abc[0]-abc[2])<1e-4*abc[0])&&
          (fabs(abc[3]-60)<1e-4)&&
          (fabs(abc[4]-60)<1e-4)&&
          (fabs(abc[5]-60)<1e-4)){
        for(i=0;i<3;i++){
          j=(i+1)%3;
          k=(i+2)%3;
          for(l=0;l<3;l++)
            new_cell[i][l]=basis[i][l]+basis[j][l]-basis[k][l];
        }

        cart2abc(new_cell,abc,0);
        super(new_cell,1);
      }

/* We also know how to convert a=b=c, alpha=beta=gamma!=60 to hexagonal */
      else if ((fabs(abc[0]-abc[1])<1e-4*abc[0])&&
               (fabs(abc[0]-abc[2])<1e-4*abc[0])&&
               (fabs(abc[4]-abc[3])<1e-4*abc[3])&&
               (fabs(abc[5]-abc[3])<1e-4*abc[3])){
        for(i=0;i<3;i++)
          for(j=0;j<3;j++) new_cell[i][j]=0;

        for(i=0;i<3;i++)
          for(j=0;j<3;j++) new_cell[0][j]+=basis[i][j];

        for(i=0;i<3;i++){
          new_cell[1][i]=basis[0][i]-basis[1][i];
          new_cell[2][i]=basis[0][i]-basis[2][i];
        }

        cart2abc(new_cell,abc,0);

/* If this worked, alpha will be 60 */

        if (fabs(abc[3]-60)<1e-4) super(new_cell,1);
        else fprintf(stderr,"Expansion (-x) specified, but failed\n");
      }

    /* Finally, we know something else about alpha=beta=gamma=60 */
      else if ((fabs(abc[3]-60)<1e-4)&&
          (fabs(abc[4]-60)<1e-4)&&
          (fabs(abc[5]-60)<1e-4)){
        double amax,amin;
        amin=min(abc[0],min(abc[1],abc[2]));
        if ((fabs(floor(abc[0]/amin+0.5)-abc[0]/amin)<1e-4)&&
            (fabs(floor(abc[1]/amin+0.5)-abc[1]/amin)<1e-4)&&
            (fabs(floor(abc[2]/amin+0.5)-abc[2]/amin)<1e-4)){ /* sizes are in
                                                           * integer ratio */
          amax=max(abc[0],max(abc[1],abc[2]));
          if ((fabs(floor(amax/abc[0]+0.5)-amax/abc[0])<1e-4)&&
              (fabs(floor(amax/abc[1]+0.5)-amax/abc[1])<1e-4)&&
              (fabs(floor(amax/abc[2]+0.5)-amax/abc[2])<1e-4)){
            for(i=0;i<3;i++){
              j=(i+1)%3;
              k=(i+2)%3;
              for(l=0;l<3;l++)
                new_cell[i][l]=(amax/abc[i])*basis[i][l]+
                               (amax/abc[j])*basis[j][l]-
                               (amax/abc[k])*basis[k][l];
            }

            cart2abc(new_cell,abc,0);
            super(new_cell,1);
          }
          else fprintf(stderr,"Expansion (-x) specified, but too complicated\n");
        }
        else fprintf(stderr,"Expansion (-x) specified, but not found.\n");
      }
      else fprintf(stderr,"Expansion (-x) specified, but not found\n");
    } /* select(expand) ... */
  }  /* if (expand) ... */

  if (rotate) cart2abc(basis,abc,1);

  switch(format){
    case XSF:
      xsf_write(outfile);
      break;
    case CUBE:
      cube_write(outfile,&grid1);
      break;
    case XPLOR:
      xplor_write(outfile,&grid1);
      break;
    case PDB:
      pdb_write(outfile);
      break;
    case CELL:
      cell_write(outfile);
      break;
    case CELL_ABC:
      cell_write_abc(outfile);
      break;
    case CELL_ABS:
      cell_write_abs(outfile);
      break;
    case CELL_ABC_ABS:
      cell_write_abc_abs(outfile);
      break;
    case DX:
      dx_write(outfile,&grid1);
      break;
    case VASP:
      vasp_write(outfile,&grid1);
      break;
    case XYZ:
      xyz_write(outfile,&grid1);
      break;
    case CML:
      cml_write(outfile);
      break;
    case FDF:
      fdf_write(outfile,file2,&grid1);
      break;
    case CNULL:
      break;

    default:
      fprintf(stderr,"This cannot happen. Sorry\n");
  }

  return(0);
}

void help(void){
  printf("Usage: check2xsf [-aAbBckmRsSvx] [--FORMAT] infile [outfile]\n\n"
         "-a           rotate as though outputing in abc format\n"
         "-A           accumulate (sum) requested bands\n"
         "-b[=range]   include bands (as psi)\n"
         "-B[=range]   include bands (as densities)\n"
         "-c           include charge density\n"
         "-H           shift atoms by half a grid cell\n"
         "-k=range     include given k-points (default 1) for bands\n"
         "-m[=a,b,c]   assume input is molecule, not crystal, and move by\n"
         "               given nos of grid cells, or move automatically\n"
         "-R           don't rescale densities\n"
         "-s           include spin densities\n"
         "-S=range     include given spins (0 or 1) for bands\n"
         "-t=(x1,y1,z1)(x2,y2,z2)\n"
         "             rotate coords so 1st vector becomes 2nd\n"
         "-T=(x1,y1,z1)(x2,y2,z2)\n"
         "             ditto, but first vec expressed in absolute coords\n"
         "-v           be verbose (may be repeated)\n"
         "-x           expand rhombohedral cell to cubic/hexagonal "
         "automatically\n");
  printf("-x=(x1,y1,z1)(x2,y2,z2)(x3,y3,z3)\n"
         "             re-express in new basis given in terms of old\n"
         "-X=(x1,y1,z1)(x2,y2,z2)(x3,y3,z3)\n"
         "             re-express in new basis given in absolute terms\n\n");
  printf("FORMAT is one of: xsf       XCrySDen (default)\n"
         "                  cell      CASTEP .cell, cartesian and fractional\n"
         "                  cell_abc                abc and fractional\n"
         "                  cell_abs                cartesian and absolute\n"
         "                  cell_abc_abs            abc and absolute\n"
         "                  cml       Chemical Markup Language\n"
         "                  cube      Gaussian cube\n"
         "                  dx        OpenDX\n"
         "                  fdf       Flexible Data Format (Siesta) (beta)\n"
         "                  null      Discard output\n"
         "                  pdb       PDB\n"
         "                  pdbn      PDB with atoms numbered\n"
         "                  vasp      VASP output\n"
         "                  xplor     Xplor\n"
         "                  xyz       XYZ\n\n");
  printf("range specifies band numbers as \"a,b-c,d\"\n"
         "-b and -B are mutually exclusive. Only one of the x and t"
         " options may be given\n\n"
         "Input files ending .pdb are assumed to be in pdb format, otherwise\n"
         "automatic detection of .cell or .check input. Compatible with\n"
         ".check files from CASTEP 3.0 to 4.1 (and perhaps beyond).\n\n"
         "Version " C2XSF_VER ", (c) MJ Rutter 2007,"
         " licenced under the GPL v2\n\n");
  exit(0);
}

void error_exit(char *msg){
  fprintf(stderr,"Aborting: %s\n",msg);
  exit(1);
}


