/* Read some useful data from a CASTEP .check file
 *
 * Cope with either endianness
 *
 * Make various assumptions about the format...
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
#include<stdlib.h> /* malloc */
#include<string.h> /* memcpy */
#include<errno.h>
#include<math.h>   /*sqrt*/

#include "c2xsf.h"
#include "c2xsf_extern.h"

#define H_eV 27.21138342902473  /* Right or wrong, it's Castep's value
                                 * physics.nist.gov gives 27.2113845(23)
                                 */

/* HEAD_LEN must be >=30 */
#define HEAD_LEN 30

/* scan for token x, read data of length len_x to targ, allocating
 * if required. Reverse endian if len_x==4
 */
#define SCAN(x,len_x,targ) \
  if (match(head,x)){ \
    fread(&tmp,4,1,infile); \
    if (endian) reverse4(&tmp); \
    if (tmp!=len_x) error_exit("Unexpected data in " #x); \
    if ((!targ)&&!(targ=malloc(len_x))) error_exit("Malloc error in " #x); \
    fread(targ,len_x,1,infile); \
    if (endian && (len_x==4)) reverse4(targ); \
    fseek(infile,4,SEEK_CUR); \
    if (debug>2) fprintf(stderr,"Found %s, %d bytes\n", #x, len_x); \
  }

int match(char *s1, char *s2);
int inrange(int x, char *range);
static void reverse4(void *data);
static void reverse4n(int *data,int n);
static void reverse8n(double *data,int n);
static void wave_read(FILE *infile, int kpt, int spin);

void fft3d(double *c, int *ngptar, int dir);


/* variables shared with wavefunction reader */

static int endian;
static int *nkpts,nbands;
static double *cr_cell_vol;

int check_read(FILE* infile){
  int tmp;
  char head[HEAD_LEN+1];
  int i,j,k,start,density,na,fft[3],ion_sym_len,ver_maj,ver_min;
  double *dptr1,*dptr2,*column,(*castep_basis)[3],conv;
  int *nsp,*nsp_max,*nionsp;
  char *ion_sym;
  double *ion_pos,*ion_force;
  struct grid *gptr,*gptr2;
  int section;

  ion_sym_len=8;
  spins=1;
  endian=0;
  castep_basis=NULL;
  ion_pos=ion_force=NULL;
  nsp=nsp_max=nionsp=NULL;
  ion_sym=NULL;
  gptr=gptr2=NULL;
  section=1;

  if(debug>2) fprintf(stderr,"check_read called with flags=%d\n",flags);

  /* The first record is a string of length 30. Being Fortran, the first
   * item will therefore be an integer, 4 bytes, of value 30. If we
   * have an endian problem, we will see this as an integer of value
   * 30*(1<<24)
   *
   * Unless we have the new "castep_bin" file format, in which case
   * the first record has a length of 10
   */

  fread(&tmp,4,1,infile);

  if ((tmp!=30)&&(tmp!=10)){ /* Oh dear, wrong magic */
    if ((tmp==(30*(1<<24)))||(tmp==(10*(1<<24)))){
      endian=1;
    }else{
      fprintf(stderr,"Error: unable to read input as a CASTEP .check file.\n");
      exit(1);
    }
  }

  rewind(infile);
  start=1;
  density=0;

/* FORTRAN stores records as  int32     length
 *                           <length>   data
 *                            int32     length  (i.e. repeated of first)
 *
 * consider check file as five sections: start to 1st END_CELL_GLOBAL
 *                                       rest to the start of wavefuncts
 *                                       wavefunctions
 *                                       charge (& spin) density
 *                                       forces
 *
 * note that the cell appears twice, and we want the first copy
 * note that the wavefunctions will be absent in a .castep_bin file
 */

  while(fread(&tmp,4,1,infile)){
    if (endian) reverse4(&tmp);
    i=(tmp<HEAD_LEN)? tmp : HEAD_LEN;
    fread(head,i,1,infile);
    head[i]=0;
    if (debug>3) fprintf(stderr,"%d: %.4s\n",tmp,head);
    if (i<tmp) fseek(infile,tmp-i,SEEK_CUR);
    fseek(infile,4,SEEK_CUR); /* capricious */

    /* Need to find version of file so that we know how long the ionic
     * species symbol records will be. Version is stored in ASCII
     * as major.minor immediately after the BEGIN_PARAMETERS_DUMP
     * marker
     */

    switch (section){
    case(1):
      if(match(head,"BEGIN_PARAMETERS_DUMP")){
        fread(&tmp,4,1,infile);
        if (endian) reverse4(&tmp);
        i=(tmp<HEAD_LEN)? tmp : HEAD_LEN;
        fread(head,i,1,infile);
        fseek(infile,4,SEEK_CUR);
        head[i]=0;
        j=0;
        ver_maj=0;
        /* Swallow leading spaces */
        while(head[j]==' ') j++;
        while((head[j])&&(head[j]!='.')&&(j<i))
          ver_maj=10*ver_maj+head[j++]-'0';
        if (head[j]=='.') j++;
        ver_min=0;
        while((head[j])&&(head[j]!=' ')&&(j<i))
          ver_min=10*ver_min+head[j++]-'0';
        if(ver_maj<4) ion_sym_len=3;
        if (debug>1) fprintf(stderr,"Castep version %d.%d\n",ver_maj,ver_min);
      }

      SCAN("CELL%REAL_LATTICE",72,castep_basis);
      SCAN("CELL%NUM_SPECIES",4,nsp);
      SCAN("CELL%MAX_IONS_IN_SPECIES",4,nsp_max);
      SCAN("CELL%NUM_IONS_IN_SPECIES",(4*(*nsp)),nionsp);
      SCAN("CELL%IONIC_POSITIONS",(24*(*nsp)*(*nsp_max)),ion_pos);
      SCAN("CELL%SPECIES_SYMBOL",(ion_sym_len*(*nsp)),ion_sym);
      SCAN("CELL%VOLUME",8,cr_cell_vol);
      SCAN("NKPTS",4,nkpts);
      if (match(head,"BEGIN_ELECTRONIC")){
        fseek(infile,3*12+4,SEEK_CUR);
        fread(&spins,4,1,infile);
        if (endian) reverse4(&spins);
        if (spins!=2) spins=1;
        fseek(infile,8,SEEK_CUR);
        fread(&nbands,4,1,infile);
        if (endian) reverse4(&nbands);
        fseek(infile,4,SEEK_CUR);
        if (debug>2) fprintf(stderr,"Spins=%d  nbands=%d\n",spins,nbands);
      }
      if (match(head,"END_CELL_GLOBAL")) {
        /* We have an urgent need for the volume so that densities can
         * be scaled as soon as they are read
         */
        if (endian) reverse8n(cr_cell_vol,1);
        (*cr_cell_vol)*=BOHR*BOHR*BOHR; /* CASTEP uses Bohrs,
                                           we use Angstroms */
        section=2;
      }
      break;
    case 2:      /* matches select(section) */
      if (match(head,"END_CELL_GLOBAL")){
        if (debug>1){
          double energy;
          int wave_ok,den_ok;
          fseek(infile,4,SEEK_CUR);
          fread(&wave_ok,4,1,infile);
          fseek(infile,8,SEEK_CUR);
          fread(&den_ok,4,1,infile);
          fseek(infile,8,SEEK_CUR);
          fread(&energy,8,1,infile);
          if (endian) reverse8n(&energy,1);
          fprintf(stderr,"Total Energy %.6f eV\n",energy*H_eV);
          fseek(infile,8,SEEK_CUR);
          fread(&energy,8,1,infile);
          if (endian) reverse8n(&energy,1);
          fprintf(stderr,"Fermi Energy %.6f eV\n",energy*H_eV);
          fseek(infile,4,SEEK_CUR);
        }
        section=3;
      }
      break;
    case 3: case 4: /* matches select(section) */
      if (match(head,"wave")){
        if (flags&BANDS) wave_read(infile, 1, 0);
        section=4;
      }
      if (density==0){
        if (tmp==12){ /* we might have an FFT grid */
          memcpy(fft,head,12);
          fread(&tmp,4,1,infile);
          if (endian) {
            reverse4(&tmp);
            reverse4(fft);
            reverse4(fft+1);
            reverse4(fft+2);
          }
          if (debug>2) fprintf(stderr,"Potential FFT grid     %d %d %d\n",
                 fft[0],fft[1],fft[2]);
          if (tmp==(spins*16*fft[2]+8)){ /* it might well be an FFT grid
               It will be stored as complex, we will read a complex column,
               then store as real */
            if (debug>2)
              fprintf(stderr,"We think we have found a charge density grid\n");
            if (flags&CHDEN){
              gptr=grid_next;
              grid_next=malloc(sizeof(struct grid));
              if (!grid_next) {
                fprintf(stderr,"Malloc error\n");
                exit(1);
              }
              grid_next->next=NULL;
              grid_next->data=NULL;
              gptr->next=grid_next;
              gptr->name="Density";
              for(i=0;i<3;i++) gptr->size[i]=fft[i];
              gptr->data=malloc(8*fft[0]*fft[1]*fft[2]);
              if (!gptr->data){
                fprintf(stderr,"Error allocating density grid\n");
                exit(1);
              }
            }
            if ((spins==2)&&(flags&SPINDEN)){
              gptr2=grid_next;
              grid_next=malloc(sizeof(struct grid));
              if (!grid_next) {
                fprintf(stderr,"Malloc error\n");
                exit(1);
              }
              grid_next->next=NULL;
              grid_next->data=NULL;
              gptr2->next=grid_next;
              gptr2->name="Spin";
              for(i=0;i<3;i++) gptr2->size[i]=fft[i];
              gptr2->data=malloc(8*fft[0]*fft[1]*fft[2]);
              if (!gptr2->data){
                fprintf(stderr,"Error allocating density grid\n");
                exit(1);
              }
            }
            column=malloc(16*fft[2]);
            while(tmp==(spins*16*fft[2]+8)){
              fread(&i,4,1,infile);
              if (endian) reverse4(&i);
              fread(&j,4,1,infile);
              if (endian) reverse4(&j);
              fread(column,16*fft[2],1,infile);
              if (gptr){
                dptr1=gptr->data+((i-1)*fft[1]+(j-1))*fft[2];
                dptr2=column;
                for(k=0;k<fft[2];k++){*dptr1++=*dptr2;dptr2+=2;}
              }
              if (spins==2){
                if (gptr2){
                  fread(column,16*fft[2],1,infile);
                  dptr1=gptr2->data+((i-1)*fft[1]+(j-1))*fft[2];
                  dptr2=column;
                  for(k=0;k<fft[2];k++){*dptr1++=*dptr2;dptr2+=2;}
                }
                else fseek(infile,16*fft[2],SEEK_CUR);
              }
              fseek(infile,4,SEEK_CUR);
              fread(&tmp,4,1,infile);
              if (endian) reverse4(&tmp);
            }
            density=1;
            free(column);
            fseek(infile,-4,SEEK_CUR);
            /* Correct endianness */
            if (endian){
              if (gptr) reverse8n(gptr->data,fft[0]*fft[1]*fft[2]);
              if (gptr2) reverse8n(gptr2->data,fft[0]*fft[1]*fft[2]);
            }
            /* Scale */
            if (((flags&RAW)==0)&&((gptr)||(gptr2))){
              conv=1/(*cr_cell_vol);
              if (debug>2) fprintf(stderr,"Rescaling densities by %f\n",conv);
              if (gptr)
                for(i=0;i<fft[0]*fft[1]*fft[2];i++) gptr->data[i]*=conv;
              if (gptr2)
                for(i=0;i<fft[0]*fft[1]*fft[2];i++) gptr2->data[i]*=conv;
            }
            section=5;
          }else{
            fseek(infile,4+tmp,SEEK_CUR);
          }
        }
      }
      break;
    case 5:  
        SCAN("FORCES",24*(*nsp)*(*nsp_max),ion_force);
        break;
    } /* end select(section) */
  } /* end while(fread(&tmp,4...) */

/* Now sort out the endianness of everything we've read */

  if (endian){
    reverse8n(castep_basis,9);
/* remember single 4 byte objects will have been converted anyway */
    if (*nsp>1) reverse4n(nionsp,*nsp);
    reverse8n(ion_pos,3*(*nsp)*(*nsp_max));
    if (ion_force) reverse8n(ion_force,3*(*nsp)*(*nsp_max));
  }

/* change basis from Fortran ordering to C ordering */

  if(!(basis=malloc(72))) error_exit("Error in malloc for basis");
  for(i=0;i<=2;i++)
    for(j=0;j<=2;j++) basis[i][j]=castep_basis[j][i];

  /* Convert basis to Angstoms */
  for(i=0;i<3;i++)
    for(j=0;j<3;j++) basis[i][j]=basis[i][j]*BOHR;

/* Add reciprocal basis and volume */

  real2rec();

  if ((fabs(cell_vol-(*cr_cell_vol))/cell_vol)>0.005*cell_vol)
    fprintf(stderr,"Warning: calc cell vol %f, CASTEP cell vol %f\n",
            cell_vol,*cr_cell_vol);

/* Pack ions */

  natoms=0;
  for(i=0;i<*nsp;i++) natoms+=nionsp[i];
  if (debug>1) fprintf(stderr,"%d atoms found\n",natoms);
  if (debug>2)
    for(i=0;i<*nsp;i++)
      fprintf(stderr,"Species %d atoms %d\n",i,nionsp[i]);

  if (!(atoms=malloc(natoms*sizeof(struct atom))))
    error_exit("Malloc error in check_read");

  na=0;
  for(i=0;i<*nsp;i++){
    for(j=0;j<nionsp[i];j++){
      atoms[na].atno=atsym2no(ion_sym+ion_sym_len*i);
      if((!atoms[na].atno)&&(debug))
        fprintf(stderr,"Warning: atom symbol %s converted to 0\n",
                ion_sym+ion_sym_len*i);
      atoms[na].frac[0]=ion_pos[3*(i*(*nsp_max)+j)];  
      atoms[na].frac[1]=ion_pos[3*(i*(*nsp_max)+j)+1];  
      atoms[na].frac[2]=ion_pos[3*(i*(*nsp_max)+j)+2];  
      if (ion_force){
        atoms[na].force[0]=ion_force[3*(i*(*nsp_max)+j)];  
        atoms[na].force[1]=ion_force[3*(i*(*nsp_max)+j)+1];  
        atoms[na].force[2]=ion_force[3*(i*(*nsp_max)+j)+2];  
      }
      else for(k=0;k<3;k++)atoms[na].force[k]=0;
      na++;
    }
  }
  if (ion_force) forces=1;

  addabs();

  if(match(head,"END")) return(1);
  return(0);

}

int match(char *s1, char *s2){
/* Returns one if s2 is the initial substring of s1 */
/* Why not use "#define match(s1,s2) (!strncmp(s1,s2,strlen(s2)))"? */
  int i;

  i=0;
  while(s1[i] && s2[i] && s1[i]==s2[i]) i++;

  if(!s2[i]) return(1);
  else return(0);
}

static void reverse4(void *data){
/* reverse endian a single 4 byte int */
   int out;
   char *p1,*p2;

   p1=(char*)data;
   p2=(char*)&out;

   p2=p2+3;

   *(p2--)=*(p1++);
   *(p2--)=*(p1++);
   *(p2--)=*(p1++);
   *(p2--)=*(p1++);

   *((int*)data)=out;
}

static void reverse4n(int *data, int n){
/* reverse endian n words of 4 byte data */
   int i,out;
   char *p1,*p2;

   for(i=0;i<n;i++){
     p1=(char*)(data+i);
     p2=(char*)&out;
  
     p2=p2+3;
  
     *(p2--)=*(p1++);
     *(p2--)=*(p1++);
     *(p2--)=*(p1++);
     *(p2--)=*(p1++);
   
     *(data+i)=out;
   }
}

static void reverse8n(double *data, int n){
/* reverse endian n words of 8 byte data */
   int i;
   double out;
   char *p1,*p2;

   for(i=0;i<n;i++){
     p1=(char*)(data+i);
     p2=(char*)&out;
  
     p2=p2+7;
  
     *(p2--)=*(p1++);
     *(p2--)=*(p1++);
     *(p2--)=*(p1++);
     *(p2--)=*(p1++);
     *(p2--)=*(p1++);
     *(p2--)=*(p1++);
     *(p2--)=*(p1++);
     *(p2--)=*(p1++);
   
     *(data+i)=out;
   }
}

static void wave_read(FILE *infile, int kpt, int spin){
  int tmp,ns,b,i,k,nplwv,fft[3],ffft[3],fft_pts,dummy[4],offset;
  int *pwgrid=NULL;
  double *dptr1,*dptr2,*dptr3,sum,phase_r,phase_i,phase_r2,phase_i2,dtmp;
  double max,min,conv;
  double kpoint[3],phi;
  int ii,jj,kk,ind;
  struct grid *g;

  dptr2=NULL;
  g=grid_next;

/* Assume we have just read the key "wave"
 *
 * we are therefore faced with:
 * ngx,ngy,ngz
 * max_plane_waves,nbands,nkpts,nspins
 * do over spins
 *   do over kpoints
 *     kpoint[3],num_plane_waves
 *     do i=1,3
 *       plane_wave_to_grid_point_component_i_mapping_array[num_plane_waves]
 *     do over bands
 *       wavefunction(complex)[num_plane_waves]
 * do over kpoints
 *   kpoint[3]
 *   do over spins
 *     band_occupancies[nbands]
 *     band_evalues[nbands]
 *
 * I am not convinced that the two "do over kpoints" loops will always see
 * the kpoints in the same order. As they explicitly record the kpoint
 * co-ordinates, there is no reason for confusion...
 */

  fread(&tmp,4,1,infile);
  if (endian) reverse4(&tmp);

  if (tmp!=12) error_exit("Error parsing wavefunction");

  fread(fft,12,1,infile);
  fseek(infile,4,SEEK_CUR);
  if (endian) {
     reverse4(fft);
     reverse4(fft+1);
     reverse4(fft+2);
  }
  fft_pts=fft[0]*fft[1]*fft[2];
  if (debug) fprintf(stderr,"Wavefunction grid     %d %d %d\n",
                 fft[0],fft[1],fft[2]);


  fread(&tmp,4,1,infile);
  if (endian) reverse4(&tmp);

  if (tmp!=16) error_exit("Error parsing wavefunction");
  fread(dummy,16,1,infile);
  if (endian) reverse4n(dummy,4);
  if (debug>2) fprintf(stderr,"Read: %d %d %d %d\n",dummy[0],
                        dummy[1],dummy[2],dummy[3]);
  fseek(infile,4,SEEK_CUR); 

  for(ns=0;ns<spins;ns++){
    for(k=1;k<=*nkpts;k++){
/* record of kpoint[3],nplwv */
    fseek(infile,4,SEEK_CUR);
    fread(kpoint,3*8,1,infile);
    fread(&nplwv,4,1,infile);
    if (endian) reverse8n(kpoint,3);
    if (endian) reverse4(&nplwv);
    if(debug>2) fprintf(stderr,"kpoint no %d (%f,%f,%f) nplwv=%d\n",
                 k,kpoint[0],kpoint[1],kpoint[2],nplwv);
    fseek(infile,4,SEEK_CUR);
/* Read component to FFT grid mapping */
      fread(&tmp,4,1,infile);
      if (endian) reverse4(&tmp);
      if(4*nplwv!=tmp) error_exit("Error parsing band");
      if(debug>2) fprintf(stderr,"Spin=%d kpt=%d nplwv=%d\n",ns,k,nplwv);
      if(inrange(k,kpt_range)&&inrange(ns,spin_range)){
        if (!(pwgrid=malloc(12*nplwv))) error_exit("Malloc error for pwgrid");
        fread(pwgrid,tmp,1,infile);
        fseek(infile,8,SEEK_CUR);
        fread(pwgrid+nplwv,tmp,1,infile);
        fseek(infile,8,SEEK_CUR);
        fread(pwgrid+2*nplwv,tmp,1,infile);
        fseek(infile,4,SEEK_CUR);
        if(endian) reverse4n(pwgrid,3*nplwv);
        if (debug>2) fprintf(stderr,"Read pwgrid\n");
/* CASTEP stores these reciprocal space coeffs as -n/2 to n/2, whereas
   we want 0 to n/2, -n/2 to -1, so .. */

        for(i=0;i<nplwv;i++){
          if(pwgrid[i]<0) pwgrid[i]+=fft[0];
          if(pwgrid[nplwv+i]<0) pwgrid[nplwv+i]+=fft[1];
          if(pwgrid[2*nplwv+i]<0) pwgrid[2*nplwv+i]+=fft[2];
        }
      }
      else fseek(infile,tmp*3+5*4,SEEK_CUR);
      for(b=1;b<=nbands;b++){
        fread(&tmp,4,1,infile);
        if (debug>2) fprintf(stderr,"Start band %d\n",b);
        if (endian) reverse4(&tmp);
        if (tmp!=16*nplwv) error_exit("Error parsing wavefunction band");
        if (inrange(k,kpt_range)&&inrange(ns,spin_range)&&
            inrange(b,band_range)){
          if (debug>2) fprintf(stderr,"starting band read\n");
          if (!(dptr1=malloc(16*nplwv)))
            error_exit("Malloc error in band read");
          fread(dptr1,16*nplwv,1,infile);
          fseek(infile,4,SEEK_CUR);
          if(endian) reverse8n(dptr1,2*nplwv);
          if (!(dptr2=malloc(16*fft_pts)))
              error_exit("Malloc error in band read");
          for(i=0;i<2*fft_pts;i++) dptr2[i]=0.0;
          sum=0;
          for(i=0;i<nplwv;i++){
/*!!*/
            offset=pwgrid[i+2*nplwv]+fft[2]*(pwgrid[i+nplwv]+fft[1]*pwgrid[i]);
            if ((offset<0)||(offset>fft_pts)){
              fprintf(stderr,"Impossible offset in wave_read off=%d i=%d\n",
                      offset,i);
              exit(1);
            }
            sum+=dptr1[2*i]*dptr1[2*i]+dptr1[2*i+1]*dptr1[2*i+1];
            dptr2[2*offset]=dptr1[2*i];
            dptr2[2*offset+1]=dptr1[2*i+1];
           }
          free(dptr1);
          if (debug>2) fprintf(stderr,"Before FFT g=0 component is %g+%gi\n",
                    dptr2[0],dptr2[1]);
          if (debug>2) fprintf(stderr,"And normalisation is %g\n",sum);
/* A FORTRAN data order ... */
          ffft[0]=fft[2];
          ffft[1]=fft[1];
          ffft[2]=fft[0];
          fft3d(dptr2,ffft,-1);

          if (!(dptr3=malloc(8*fft_pts)))
             error_exit("Malloc error in band read");

          if (((flags&BANDDEN)==BANDS)&&((kpoint[0]!=0)||
               (kpoint[1]!=0)||(kpoint[2]!=0))){ /* want psi,
                                                    but not at gamma! */
            if (debug)
              fprintf(stderr,"unwinding psi for non-gamma k-point...\n");
            for(ii=0;ii<fft[0];ii++){
              for(jj=0;jj<fft[1];jj++){
                for(kk=0;kk<fft[2];kk++){
                  phi=2*M_PI*((ii*kpoint[0])/fft[0]+
                              (jj*kpoint[1])/fft[1]+
                              (kk*kpoint[2])/fft[2]);
                  phase_r=cos(phi);
                  phase_i=sin(phi);
                  ind=2*(kk+fft[2]*(jj+ii*fft[1]));
                  dtmp=dptr2[ind];
                  dptr2[ind]=phase_r*dptr2[ind]-phase_i*dptr2[ind+1];
                  dptr2[ind+1]=phase_r*dptr2[ind+1]+phase_i*dtmp;
                }
              }
            }
          }
          phase_r=phase_i=phase_r2=phase_i2=0;
          for(i=0;i<fft_pts;i++){
            if (dptr2[2*i]>0){
              phase_r+=dptr2[2*i];
              phase_i-=dptr2[2*i+1];
            }else{
              phase_r2-=dptr2[2*i];
              phase_i2+=dptr2[2*i+1];
            }
          } 
          phase_r+=phase_r2;
          phase_i+=phase_i2;
          dtmp=sqrt(phase_r*phase_r+phase_i*phase_i);
          phase_r/=dtmp;
          phase_i/=dtmp;
          ii=0;
          max=-1e300;
          min=1e300;
          for (i=0;i<fft_pts;i++){
            if (flags&BANDPHASE){
              dptr3[i]=atan2(dptr2[2*i+1],dptr2[2*i]);
            }
            else if (flags&BANDREAL){
              dptr3[i]=dptr2[2*i];
            }
            else if (flags&BANDIMAG){
              dptr3[i]=dptr2[2*i+1];
            }
            else
            if ((flags&BANDDEN)==BANDDEN)
              dptr3[i]=dptr2[2*i]*dptr2[2*i]+dptr2[2*i+1]*dptr2[2*i+1];
            else{
              dptr3[i]=dptr2[2*i]*phase_r-dptr2[2*i+1]*phase_i;
              dtmp=dptr2[2*i]*phase_i+dptr2[2*i+1]*phase_r;
              if((fabs(dtmp)>.05))ii++;
            }
            sum+=dptr3[i];
            if(dptr3[i]<min) min=dptr3[i];
            if(dptr3[i]>max) max=dptr3[i];
          }
          if (debug>2) fprintf(stderr,"Min=%g Max=%g Sum=%g\n",min,max,sum);
          if (ii>0) fprintf(stderr,"Warning: %d components with imaginary"
                                   " part >0.05\n",ii);
          free(dptr2);
          /* Do we need to rescale? */
          if (((flags&RAW)==0)&&((flags&BANDPHASE)==0)){ /* Yes */
            if ((flags&BANDDEN)==BANDDEN) conv=1/(*cr_cell_vol);
            else conv=1/sqrt(*cr_cell_vol);
            if (debug>2) fprintf(stderr,"Scaling wavefun by %f\n",conv);
            for(i=0;i<fft_pts;i++) dptr3[i]*=conv;
          }
          if (!(flags&ACCUMULATE)){
            g->data=dptr3;
            for(i=0;i<3;i++) g->size[i]=fft[i];
            g->name=malloc(40);
            if (spins==2)
              sprintf(g->name,"band_s%d_k%d_b%d",spin,k,b);
            else
              sprintf(g->name,"band_k%d_b%d",k,b);
            g->next=malloc(sizeof(struct grid));
            grid_next=g->next;
            g=g->next;
            g->data=NULL;
            g->next=NULL;
          }else{  /* Are accumulating */
            if (!g->data){  /* This is the first set */
              g->data=dptr3;
              for(i=0;i<3;i++) g->size[i]=fft[i];
              g->name=malloc(40);
              sprintf(g->name,"bands"); /* Don't move to a new grid */
            }else{
              for(i=0;i<fft_pts;i++) g->data[i]+=dptr3[i];
              free(dptr3);
            }
          }
        }
        else fseek(infile,16*nplwv+4,SEEK_CUR);
      }
    }
  }
  if (pwgrid) free(pwgrid);
  if (debug>1){ /* Try reporting evals and occupancies */
    if (!(dptr1=malloc(8*nbands))) error_exit("Malloc error for occs\n");
    if (!(dptr2=malloc(8*nbands))) error_exit("Malloc error for evals\n");
    fprintf(stderr,"kpoint, band, spin, occupancy, evalue (eV):\n");
    for (k=0;k<*nkpts;k++){
      fread(&tmp,4,1,infile);
      if (endian) reverse4(&tmp);
      if (tmp!=24) error_exit("Error parsing end of wavefunction");
      fread(kpoint,3*8,1,infile);
      if (endian) reverse8n(kpoint,3);
      fseek(infile,4,SEEK_CUR); 
      for(ns=0;ns<spins;ns++){
        fread(&tmp,4,1,infile);
        if (endian) reverse4(&tmp);
        if (tmp!=8*nbands) error_exit("Error parsing end of wavefunction");
        fread(dptr1,8*nbands,1,infile);
        if (endian) reverse8n(dptr1,nbands);
        fseek(infile,4,SEEK_CUR);
        fread(&tmp,4,1,infile);
        if (endian) reverse4(&tmp);
        if (tmp!=8*nbands) error_exit("Error parsing end of wavefunction");
        fread(dptr2,8*nbands,1,infile);
        if (endian) reverse8n(dptr2,nbands);
        fseek(infile,4,SEEK_CUR);
        for(b=0;b<nbands;b++){
          fprintf(stderr,"( %8f %8f %8f ) %d %d  %10f   %14f\n",
                  kpoint[0],kpoint[1],kpoint[2],b+1,ns,dptr1[b],dptr2[b]*H_eV);
        }
      }
    }
  }

}

int inrange(int x, char *range){
/* determine whether x is within a range given by comma-separated
 * whitespace-free list of ints and hyphen-separated ranges
 */
  char *cptr;
  int test1,test2;

  cptr=range;
  if (!*cptr) return(0);
  while(*cptr){
    test1=0;
    while ((*cptr>='0')&&(*cptr<='9'))
      test1=10*test1+(*(cptr++)-'0');
    if (*cptr=='-'){ /* a range */
      cptr++;
      if((*cptr>='0')&&(*cptr<='9')){
        test2=0;
        while ((*cptr>='0')&&(*cptr<='9'))
          test2=10*test2+(*(cptr++)-'0');
      } else test2=1<<30;
    }else test2=test1;

    if ((*cptr!=',')&&(*cptr!=0)) error_exit("Parse error in inrange");

    if ((x>=test1)&&(x<=test2)) return(1);

    if (*cptr==',') cptr++;
  }
  return(0);
}
