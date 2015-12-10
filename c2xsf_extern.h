extern double (*basis)[3],recip[3][3]; /* Basis sets */
extern double cell_vol;
extern int natoms,spins;
extern int debug,molecule,flags,forces;
extern char *band_range,*kpt_range,*spin_range;
extern struct atom *atoms;
extern struct grid grid1, *grid_next;
