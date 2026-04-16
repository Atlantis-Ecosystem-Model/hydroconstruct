/*********************************************************************

 File:           mainexchange.h

 Created:        Mon Oct 20 16:57:26 EST 2003

 Author:         Beth Fulton
 CSIRO Division of Marine Research

 Purpose:        Include file for exchange calculator to accompany box models

 Revisions:      28/11/2004 updated to allow multiple input files and
 temperature and salinity profiles

 12-03-2009 Bec Gorton
 Added additional parameters required for loading Cams data.
 Additional arrays:
 int ***col_faceid;
 double ***col_faceProportion;
 double **depth_layer;
 also added the layerDepthIfname field to Stuff to store
 the name of the depth layer information file.

 15-04-2009 Bec Gorton
 Added the max_nconn variable to the stuff structure.

 02-06-2009 Bec Gorton
 Added support for the new ECCAL hydro file format. Defined the new input file type eccal_flat_type and
 added a new function definition get_Eccal_flat_hydro() to deal with this format.

 14-09-2009 Bec Gorton
 Added the vertexchange name so that we can read in the vertical exchanges from the
 ECCAL input file.


 *********************************************************************/
#include <sjwlib.h>

/* Type of floating point values for tracers, epi variables */
#define FPTYPE  double

/* Useful macros */
#if !defined(min)
#define min(x,y)  ( (x) < (y) ? (x) : (y) )
#endif

#if !defined(min)
#define max(x,y)  ( (x) > (y) ? (x) : (y) )
#endif

/* Miscellaneous defines */
#define STSLEN    (512)     /* Max string lengths */

#if !defined(MAXDOUBLE)
#define MAXDOUBLE (1.7976931348623158e+308)
#endif

#if !defined(MINDOUBLE)
#define MINDOUBLE (2.2250738585072014e-308)
#endif

#if !defined(MAXFLOAT)
#define MAXFLOAT (3.402823466e+38)
#endif

#if !defined(MINFLOAT)
#define MINFLOAT (1.175494351e-38)
#endif

#if !defined(MAXINT)
#define MAXINT (2147483647)
#endif

#define tiny 0.0000000000000001
#define buffer_rounding 0.0001

#define temp_id 0
#define salt_id 1
#define ph_id 2
#define swr_id 3
#define exchange_id 4
#define padding_max_id 5

#define num_TSP_values 4

#define time_id 0 // like seconds since start
#define TofY_id 1 // time of year
#define done_id 2 // check to see if this slot has been done.
#define altface_id 0
#define bgmface_id 1
#define altbox_id 2
#define bgmbox_id 3
#define nalt_data_id 4

#define NOAA_type_id 0		/* Original Al Herman data type */
#define flat_type_id 1
#define cdf_type_id 2
#define col_face_type 3
#define ROMS_cdf_type 4		/* AL Herman cdf ROMS related files */
#define col_data_type 5		/* Files from Mexican oceanographers for Cam */
#define allinROMS_cdf_type 6
#define eccal_flat_type 7
#define newROMS_cdf_type 8
#define new_col_data_type 9	/* New format from CAM  - hydro file includes depth info.*/
#define guam_data_type 10	/* New format from CAM  - hydro file includes depth info.*/


#define lbox_id 0
#define rbox_id 1
#define fceid_id 2
#define px1_id 3
#define py1_id 4
#define px2_id 5
#define py2_id 6
#define length_id 7

#define num_faceinfo 8

#define metres_in_degree 111120

#define box_id 0
#define circum_id 1

#define to_north 0
#define to_south 1
#define same_lat 2

#define to_west 0
#define to_east 1
#define same_long 2

/*********************************************************************
 The face structure. The face geometry is defined by the two points
 p1 and p2. The positive direction for fluxes etc is defined
 to be to the right when looking from p1 to p2. I'm not sure whether
 this is compatible with John Hunters definition or not.
 *********************************************************************/
typedef struct {
	/* These values don't change during a run, and are read
	 * and written by readFaceGeom() and writeFaceGeom()
	 */
	int n; /* face number (starts at 0) */
	dpoint p1; /* first point coordinates */
	dpoint p2; /* second point coordinates */
	dpoint pll1; /* first point coordinates in latitude and longitude */
	dpoint pll2; /* second point coordinates in latitude and longtitude */
	double len; /* length of face */
	double cos; /* cos angle from face normal to x-axis */
	double sin; /* sin angle from face normal to x-axis */
	int ibl; /* index of box to left */
	int ibr; /* index of box to right */
} Face;

/*********************************************************************
 The hydrodynamic model data input structure
 *********************************************************************/
typedef struct {
	/* Input file management */
	char inputDir[500]; /* path to the input files - only used for col_dat_type */
	char **fname; /* List of file names */
	int nfiles; /* Number of input files */
	double missing_data; /* missing data value (flows < -(this value) will be ignored */
	int fid; /* netCDF id for current input file */
	int t_did; /* t dimension id */
	int b_did; /* b dimension id */
	int f_did; /* f dimension id */
	int z_did; /* z dimension id */
	int t_vid; /* time variable id */
	int srcb_vid; /* source box variable id */
	int destb_vid; /* destination box variable id */
	int fce_vid; /* face_id variable id */
	int fcex1_vid; /* pt1 x coordinate variable id */
	int fcey1_vid; /* pt1 y coordinate variable id */
	int fcex2_vid; /* pt2 x coordinate variable id */
	int fcey2_vid; /* pt2 y coordinate variable id */
	int e_vid; /* exchange variable id */

	char t_units[STSLEN]; /* Hydrodynamic time units */
	char e_units[STSLEN]; /* Hydrodynamic exchange units */
	double dt; /* Hydrodynamic time step */
	long nstep; /* Number of time steps in file */

} HydroData;

/*********************************************************************
 The temperature and salinity data input structure
 *********************************************************************/
typedef struct {
	/* Input file management */
	char inputDir[500]; /* path to the input files - only used for col_dat_type */
	char **tsfname; /* List of file names */
	int ntsfiles; /* Number of input files */
	int temp_missing_data; /* missing data value (temperatures < -(this value) will be ignored */
	int salt_missing_data; /* missing data value (salinities < -(this value) will be ignored */
	int ph_missing_data; /* Data of this values will be ignored */
	int tsfid; /* netCDF id for current input file */
	int t_did; /* t dimension id */
	int b_did; /* b dimension id */
	int z_did; /* z dimension id */
	int t_vid; /* time variable id */
	int b_vid; /* box variable id */
	int salt_vid; /* salinity variable id */
	int temp_vid; /* temperature box variable id */
	int ve_vid; /* vertical exchange variable id */
	int swr_vid; /* swr box variable id */

	char t_units[STSLEN]; /* Salinity and temperature hydrodynamic time units */
	char ve_units[STSLEN]; /* Hydrodynamic vertical exchange units */
	double tsdt; /* Salinity and temperature hydrodynamic time step */
	long tsnstep; /* Number of time steps in salinity and temperature raw data file */

} TSData;


typedef struct{

	int originalLeaveBox;
	int originalDestBox;
	int targetLeaveBox;
	int targetDestBox;
} SubsBoxInfo;
/*********************************************************************
 The Box Model structure
 *********************************************************************/
typedef struct {
	/* Command line arguments/parameter file names */
	char params[7*STSLEN];

	/* I/O file info for this model */
	char geomIfname[STSLEN]; /* geometry input file name */
	char llgeomIfname[STSLEN]; /* lat-long faces geometry input file name */
	char runprmIfname[STSLEN]; /* run parameter file */
	char salinityOfname[STSLEN]; /* exchange output file name */
	char temperatureOfname[STSLEN]; /* exchange output file name */
	char transportOfname[STSLEN]; /* exchange output file name */
	char pHOfname[STSLEN]; /* exchange output file name */
	char swrOfname[STSLEN]; /* exchange output file name */

	int doPH;	/* True if we are dealing with pH as well as temp and salt */
	int doSWR;

	char lookupIfname[STSLEN]; /* lookup table filename for when using Al Herman type data */
	char layerDepthIfname[STSLEN]; /* Lookup file for the layer depth info */
	char tsTEMPinitfname[STSLEN]; /* lookup table of temperature iniital conditions for use with
	 Al Herman type temperature flux data */
	char tsSALTinitfname[STSLEN]; /* lookup table of salinity iniital conditions for use with
	 Al Herman type salinity flux data */
	int numoutfile; /* Number of output files */

	/* Geometry data for this model */
	int nbox; /* number of boxes */
	int nface; /* number of open faces */
	int nbgmface; /* number of faces in bgm file (for mapping Jeff Dunn to bgm faces) */
	int nfacegeo; /* number of faces in geometry file */
	int wcnz; /* number of layers in water column */
	int ndest; /* maximum number of destinations per box */
	int fndest; /* maximum number of destinations actually visited per box */
	Face *faces; /* array of open faces */
	double *boxarea; /* area of each box */
	double *boxtotdepth; /* total depth of each box */
	double *nominal_dz; /* default size of each water column layer */
	double *boxscale; /* scalar of flows for each box */
	double **boxvol; /* volume per layer per box */
	double **boxdz; /* depth per layer per box */
	double *boxinsidex; /* x coordinate of "inside" point for each box */
	double *boxinsidey; /* y coordinate of "inside" point for each box */
	double *boxht; /* height (north-south on y coordinate) in metres of each box */
	double *boxwd; /* width (west-east on x coordinate) in metres of each box */
	double *boxnscoefft; /* coefficient of scaling for area_correction in north-south decomposition */
	double *boxwecoefft; /* coefficient of scaling for area_correction in east-west decomposition */
	int **ibox; /* id numbers of neighbopuring boxes for each box */
	int **iface; /* id numbers of neighbopuring boxes for each box */
	int **boxns; /* flag indicating whether connecting box to north or south of calling box */
	int **boxwe; /* flag indicating whether connecting box to east or west of calling box */
	int *boxnconn; /* number of boxes adjacent to each box */
	int *boxnz; /* number of layers in each box */
	int *boundaries; /* boundary state flags for each box (whether a boundary box,
	 and if so what type, reflective, absorptive etc) */

	int *faceVertIndex; /* Used to store the index of each box face in the vert list */
	/* Estuary parameters */
	int num_estuaries; /* Number of estuaries to calculate */
	int *estnum; /* ID numbers of faces marking estuary mouthes */
	double *estinflux; /* Fluxes into the estuary mouth */
	double *estoutflux; /* Fluxes out of the estuary mouth */

	/* Run parameters */
	char t_units[STSLEN]; /* time units */

	int rewind; /* flag to indicate whether to rewind file between line reads or not
	 0 = no, 1 = yes (fast if no, safer if yes) */
	int generic; /* 1 = generic calculations , 0 = hardwired/readin */
	int tstart; /* start time (days) */
	int tstop; /* stop time (days) */
	int num_entries;
	double dt; /* time step */
	double datadt; /* time step of data - for col_data_type */
	int nt; /* total number of steps */
	int uniflow; /* flag indicating whether unidorectional flows used */
	int sdiff; /* flag indicating whether or not slow diffusion allowed */
	int vdiff; /* flag indicating whether or not assumed vertical diffusion allowed */
	int bdiff; /* flag indicating whether or not there is (minimal) back
	 flow diffusion allowed */
	int recycle_flow; /* flag toindicate flows recycled (alternated) through time */
	int area_correct_flow; /* flag to indicate whether to area correct flow
	 to try and avoid hyperdiffusion issue */
	int area_correct_vflow; /* flag to indicate whether to area correct vertical flows
	 to try and avoid hyperdiffusion issue (less of a problem vertically
	 as those depth bands are in metres so usually completely mixed */
	int tsflagflux; /* flag to indicate whether absolute values or fluxes read in for
	 temperature and salinity */
	int input_style; /* 0 = Al Herman style, 1 = Jeff Dunn style input files */
	int unit_type; /* unit type, 0 = sverdrups, 1 = m3s-1 */
	int altfaces; /* number of alternative faces in lookup table for Al Herman -> bgm ids */
	int notbgmface; /* number indicating face is in Al Herman data but not bgm file (i.e. boundary face) */
	int n_inline; /* number of "datapoints" per line (i.e. flows per face per line and repeat through time
	 or all flows per face for the entire time, each entry in separate column */
	int ndata; /* total number of data time steps in complete set of
	 Al Herman style data file */
	int max_nconn; /* The maximum number of connections between boxes */
	int nhdsteps; /* minimum number of exchange entries to consider */
	int nvhdsteps; /* minimum number of vertical exchange entries to consider */
	int ntssteps; /* minimum number of temperature and salinity
	 entries to consider */
	int Exmin_id; /* ID of minimum date readin */
	int Exmax_id; /* ID of maximum date readin */
	int minExyear; /* earliest year of exchange data */
	int monthEx; /* month of earliest exchange data entry */
	int dayEx; /* day of earliest exchange data entry */
	int maxExyear; /* last year of exchange data */
	int mmonthEx; /* month of last exchange data entry */
	int mdayEx; /* day of last exchange data entry */
	int TSmin_id; /* ID of minimum date readin */
	int TSmax_id; /* ID of maximum date readin */
	int minTSyear; /* earliest year of temperature and salinity data */
	int monthTS; /* month of earliest temperature and salinity data entry */
	int dayTS; /* day of earliest temperature and salinity data entry */
	int maxTSyear; /* last year of temperature and salinity data */
	int mmonthTS; /* month of last temperature and salinity data entry */
	int mdayTS; /* day of last temperature and salinity data entry */

	int REFyear; /* year used as reference year - exchanges etc are referenced
	 from 01-01 00:00:00 (MM-DD HH:MM:SS) of this year */
	int start_day; /* Julian day in year of first line of data in Al Herman style data file */
	int timecorrect_on; /* Whether need to correct dates due to min date read */
	int reset_time; /* Whether to reset min time to match t=0 rather than min value in input files */
	double time_reset; /* Adjustment needed to reset time */

	int verbose_exchange; /* Flag to indicate want detailed exchange calc notes */
	int verbose_pad; /* Flag to indicate want detailed time padding notes */
	int verbose_box;
	int pad_time;

	/* Raw data input structure */
	HydroData hd; /* Horizontal fluxes */
	HydroData vhd; /* Vertical fluxes */
	TSData ts; /* Temperature and salinity data */

	/* Diagnostic test info */
	int nitcount; /* Number of tests to make on temperature */

	/* Names of the variables - this is used to deal with the many different versions of ROMS data we have */
	char 		timeVariableName[STSLEN];
	char 		transportVariableName[STSLEN];
	char 		vertExchangeVariableName[STSLEN];
	char 		tempVariableName[STSLEN];
	char 		saltVariableName[STSLEN];
	char 		boxVariableName[STSLEN];
	char		secIDVariableName[STSLEN];
	char		zbinVariableName[STSLEN];


	/* Information about using data from one box for another box where these is no flow */
	int 	numSubBoxes;
	SubsBoxInfo 	*subBoxInfo;
} Stuff;

/*********************************************************************
 Global variables
 *********************************************************************/
int *boxmap;
int *boxndest;
int *checkedthis;
int *checkedthisbox;
int *lookuplayer;
int *facemap;
int *faceid;
int *next_test;
int *test_ok;
int **boxdest;
int *boxLookup; /** Array to store the mappings between the input hydro file and the bgm file */
int ***col_faceid; /* This array provides linkage between the bgm face indices and the one in the flux files from Cam */
double ***col_faceProportion; /* This array stores the proportion of the flux that will be going to the face in the col_faceid array
 This allows for the fact that some of the oceanography faces span two bgm faces - this allows us to
 proportionally split the flux */
double **depth_layer; /** The depth of each layer - this is read in from the layerDepthIfname input file when
 the input type == col_data_type.
 These depth values are used to distribute the horizontal flux values loaded through the water column
 */

int **paddeddates;
int **lookuptable;
int **leftrightcheck;
int **topbotcheck;
int ***layerdest;
int ***tempfaceinfo;
double *lookuplength;
double *lookupdir;
double **exchangedate;
double **tempsaltdate;
double **lookupdz;
double ***boxstats;
double ***faceinfo;
double ***TempSaltInit;
double ****finalexchanges;
double ****rawtempsalt;
double ****rawexchanges;
double ****rawvertexchanges;
double ****finaltempsalt;
double ****tottempsaltflux;

extern double reference_time;
extern int verbose;
/*********************************************************************
 Protytypes
 *********************************************************************/
double GetTime(char *d);

void CalcTempSalt(Stuff *st);
void Calc_Wc_Vols(Stuff *st);
void close_hydro(Stuff *st);
void close_tsfile(Stuff *st);

void ExchangeStuff(Stuff *st);
void freeStuff(Stuff *st);
void freeGeom(Stuff *st);

void get_alt_data(Stuff *st);
void get_alt_temp_salt(Stuff *st);
void get_alt_vert_data(Stuff *st);
void get_cdf_hydro(Stuff *st);
void get_cdf_temp_and_salt(Stuff *st);
void get_flat_hydro(Stuff *st);
void get_Eccal_flat_hydro(Stuff *st);
void get_flat_temp_and_salt(Stuff *st);
void get_hydro(Stuff *st);
void get_netcdf_dimensions(Stuff *st);
void load_bgm_mappings(Stuff *str);
void get_ROMS_cdf_hydro(Stuff *st);
void get_ROMS_cdf_temp_and_salt(Stuff *st);
void get_ROMS_netcdf_dimensions(Stuff *st);
void get_2in1ROMS_cdf_hydro(Stuff *st);
void get_NEW2in1ROMS_cdf_hydro(Stuff *st);
void get_2in1ROMS_netcdf_dimensions(Stuff *st);
void get_temp_and_salt(Stuff *st);
void Get_Date(int startday, int startmonth, int startyear, int data_dt,
		int *day, int *month, int *year);
void Get_Info(Stuff *st, double valface, double valbox, double vallayer,
		int *bgm_face, int *bgm_box, int *bgm_layer, double *dirscale,
		double *len);

void Load_Face_Lookup_Table(Stuff *st);
void Load_Lookup_Table(Stuff *st);

void map_bgm_netcdf_faces(Stuff *st, int i);
void map_ROMS_netcdf_faces(Stuff *st);

void open_hydro(Stuff *st, char *name, int filenum);
void open_2in1ROMS_hydro(Stuff *st, char *name, int filenum);
void open_2in1ROMS_ts(Stuff *st, char *name, int filenum);
void open_ROMS_hydro(Stuff *st, char *name, int filenum);
void open_ROMStsfile(Stuff *st, char *name, int filenum);
void open_tsfile(Stuff *st, char *name, int filenum);

void PadTimeSeries(Stuff *st);

void readRunParams(char *name, Stuff *st);
void readBoxModel(char *name, long dump, Stuff *st);
void readBoxModelGeom(char *name, Stuff *st);
void readBoxGeomll(char *name, Stuff *st);
void readBoxGeom(char *name, Stuff *st);
void readFaceGeom(FILE *fp, Face *f);

void setup_NorthSouth(Stuff *st);

void write_exchange_out(Stuff *st);
void write_temp_salt_out(Stuff *st);
void writeExchangeOut(FILE *fid, Stuff *st, int numid);
void WriteNewOutput(Stuff *st);


void checkTransport(Stuff *st, int indexed_countend, int indexed_start, int numid);
void writeTimeSteps(Stuff *st, FILE *fid, int countStart, int countend);

void get_guam_hydro(Stuff *st);
void get_guam_temp_and_salt(Stuff *st);
