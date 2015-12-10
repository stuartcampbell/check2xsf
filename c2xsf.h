#define C2XSF_VER "1.0"

/* Global variables for system description */

#define BOHR 0.529177249

/* constants for the file formats which can be written */
#define XSF 0
#define CUBE 1
#define XPLOR 2
#define PDB 3
#define CELL 4
#define CELL_ABC 5
#define CELL_ABS 6
#define CELL_ABC_ABS 7
#define DX 8
#define VASP 9
#define XYZ 10
#define CML 11
#define FDF 12
#define CNULL 13

/* flags for reading and output */
#define CHDEN 1
#define SPINDEN 2
#define BANDS 4
#define BANDDEN 12
#define BANDPHASE 16
#define BANDREAL 32
#define BANDIMAG 64
#define ACCUMULATE 128
#define RAW 256
#define PDBN 512

/* And for laziness when reading .cell or .pdb file */
#define MAX_ATOMS 2000


struct atom
   {unsigned int atno; double abs[3]; double frac[3]; double force[3];};
struct grid {char *name; int size[3]; double *data; struct grid *next;};

void error_exit(char* msg);
void real2rec(void);
void addfrac(void);
void addabs(void);
void abc2cart(double *abc, double basis[3][3]);
void cart2abc(double basis[3][3], double *abc, int fix);
void print_globals(int level);

unsigned int atsym2no(char* sym);
char* atno2sym(unsigned no);

#ifndef min
#define min(a,b) ((a)<(b)?(a):(b))
#endif

#ifndef max
#define max(a,b) ((a)>(b)?(a):(b))
#endif

