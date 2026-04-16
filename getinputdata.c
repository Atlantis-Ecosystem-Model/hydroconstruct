/**********************************************************************

 getinputdata.c

 by:	Beth Fulton
 date:	28/11/2004

 comments: This file is the code that reads in run parameters and
 calls appropriate form of readin for the raw hydrodynamics
 flows, and temperature and salinity profiles

 assumptions:

 revisions:
 22/6/2005 Split out the NOAA data version to getaltinputdata.c

 2/3/2006 Split out CSIRO flat and netcdf data versions to
 getflatinputdata.c and getcdfinputdata.c

 11/12/2006 added the code to read in excel estimated flows
 (under col_face_type input_type). Also streamlined the setting
 of the nhdsteps and ndata values.

 18-03-2009 Bec Gorton

 Added code to allocate the depth_layer array.
 When the input_style == col_daa_type the st->tsflagflux = 0; is set to false.


 15-04-2009 Bec Gorton
 Added the max_nconn variable to the stuff structure.

 02-06-2009 Bec Gorton
 Added support for the new format of the ECCAL input files.
 Added the getEccalFlatFormatData.c file and added code in readRunParams to deal with this new format.

 14-09-2009 Bec Gorton
 Changed the HydroConstruct code to store the st->boxnz as the actual number of
 layers in a box not the max index of the layers.
 *************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <math.h>
#include <string.h>
#include <sjwlib.h>
#include "mainexchange.h"
//#include "windows.h"

#ifdef WIN32
#include "windows.h"
#else
#include <fcntl.h>
#include <sys/types.h>   /* for pid_t pid;                                  */
#include <signal.h>      /* for kill(pid, SIGTERM); function                */
#include <dirent.h>      /* for DIR  *d = opendir("/proc"); declaration     */
#include <sys/stat.h>    /* for struct stat sb;                             */

typedef struct dirent Dirent;
#endif

#define MAX_PATH 1024
#define MAXFILES 10000
void dirwalk(char *dir, char **dirNames, int *index, char *filePart);

int parseInputDir(char *dirName, char **fileNames, char *filePart);

/***************************************************************************************************
 Routine to get the exchange data out of the raw data file
 ****************************************************************************************************/
void readRunParams(char *name, Stuff *st) {
	FILE *fp;
	FILE *fp2;
	char key[STSLEN];
	char *boxnzi = "numlayers";
	char *boundaryi = "boundaries";
	char *boxscalei = "box_scaling";
	char *nom_dzi = "default_layer_dz";
	char *nriverii = "river_influxes";
	char *nriveroi = "river_outfluxes";
	char *riveridi = "river_ids";
	int i, nline, k = 0; // make max_nconn global.
	double *readinnz;
	double *boundary;
	double *boxscalar;
	double *riverinflux;
	double *riveroutflux;
	double *riverfcenum;

	st->max_nconn = 0;

	/* Open the parameter file */
	if ((fp = fopen(name, "r")) == NULL)
		quit("hydro_exchange_init: Can't open %s\n", st->runprmIfname);

	/* Read verbose setting */
	readkeyprm_i(fp, "verbose", &verbose);
	readkeyprm_i(fp, "verbose_exchange", &st->verbose_exchange);
	readkeyprm_i(fp, "verbose_pad", &st->verbose_pad);
	readkeyprm_i(fp, "verbose_box", &st->verbose_box);
	readkeyprm_i(fp, "rewind", &st->rewind);
	readkeyprm_i(fp, "generic", &st->generic);

	readkeyprm_i(fp, "pad_time", &st->pad_time);

	if (verbose)
		printf("Reading exchange calculation parameters\n");

	/* Read geometry filenames */
	readkeyprm_s(fp, "geofile", st->geomIfname);
	readBoxGeom(st->geomIfname, st);
	readkeyprm_s(fp, "llgeofile", st->llgeomIfname);
	readBoxGeomll(st->llgeomIfname, st);
	setup_NorthSouth(st);

	/* Whether need to recycle flows */
	readkeyprm_i(fp, "recycle_flow", &st->recycle_flow);

	/* Whether need to recycle flows */
	readkeyprm_i(fp, "unidirectional_flow", &st->uniflow);
	if (st->uniflow)
		warn("Unidirectional flows chosen, this is fine if inputs are net flows\n");

	/* Slow diffusion enabled (replaces 0.0000 flows with 0.0001 flows) */
	readkeyprm_i(fp, "slow_diffusion", &st->sdiff);

	/* Assumed vertical diffusion enabled */
	readkeyprm_i(fp, "vert_diffusion", &st->vdiff);

	/* Assumed (minimal) back diffusion enabled */
	readkeyprm_i(fp, "back_diffusion", &st->bdiff);

	/* Read reference year */
	readkeyprm_i(fp, "reference_year", &st->REFyear);

	/* Read time file starts */
	readkeyprm_i(fp, "tstart", &st->tstart);

	/* Read period of time to be covered */
	readkeyprm_i(fp, "tstop", &st->tstop);
	st->tstop++;

	/* Read ntimestep of output file */
	readkeyprm_d(fp, "dt", &st->dt);

	if(st->dt > 86400){
		st->num_entries = (int)(st->tstop / (st->dt/86400));
	}else{
		st->num_entries = (int)(st->tstop * (86400 / st->dt)) - 1;
	}


	/* Multiple by 86400 / dt as may be tidal entries not daily
	 and need expanded arrays to match */
	st->nt = (int) ((st->tstop - st->tstart) * (86400 / st->dt));


	if(st->nt < 0){
		printf("st->tstop = %d\n", st->tstop);
		printf("st->tstart = %d\n", st->tstart);
		printf("st->dt = %e\n", st->dt);

		abort();
	}

	/* Check for reset time flag */
	readkeyprm_i(fp, "reset_time", &st->reset_time);

	/* Read number of water column layers */
	readkeyprm_i(fp, "wcnz", &st->wcnz);
	printf("st->wcnz = %d\n", st->wcnz);

	/* Read in number of water column layers per box */
	st->nominal_dz = (double *) d_alloc1d(st->wcnz);

	i = st->wcnz;
	readkeyprm_darray(fp, nom_dzi, &st->nominal_dz, &i);

	printf("st->nbox = %d\n", st->nbox);
	readinnz = (double *) d_alloc1d(st->nbox);
	boundary = (double *) d_alloc1d(st->nbox);
	boxscalar = (double *) d_alloc1d(st->nbox);
	st->boxnz = (int *) i_alloc1d(st->nbox);
	st->boundaries = (int *) i_alloc1d(st->nbox);
	st->boxscale = (double *) d_alloc1d(st->nbox);
	TempSaltInit = (double ***) d_alloc3d(2, st->wcnz, st->nbox);
	lookupdz = (double **) d_alloc2d(st->wcnz, st->nbox);

	st->max_nconn = 0;
	for (i = 0; i < st->nbox; i++) {
		if (st->boxnconn[i] > st->max_nconn)
			st->max_nconn = st->boxnconn[i];
	}
	col_faceid = (int ***) i_alloc3d(2, st->max_nconn, st->nbox);
	for (i = 0; i < st->nbox; i++) {
		for (k = 0; k < st->max_nconn; k++) {
			col_faceid[i][k][0] = -1;
			col_faceid[i][k][1] = -1;
		}
	}

	col_faceProportion = (double ***) d_alloc3d(2, st->max_nconn, st->nbox);
	for (i = 0; i < st->nbox; i++) {
		for (k = 0; k < st->max_nconn; k++) {
			col_faceProportion[i][k][0] = 1.0;
			col_faceProportion[i][k][1] = 1.0;
		}
	}

	depth_layer = (double **) d_alloc2d(st->wcnz, st->nbox);

	i = st->nbox;
	readkeyprm_darray(fp, boxnzi, &readinnz, &i);
	readkeyprm_darray(fp, boundaryi, &boundary, &i);
	readkeyprm_darray(fp, boxscalei, &boxscalar, &i);

	if (verbose)
		printf("Flow scaling: ");

	readkeyprm_i(fp, "area_correct_flow", &st->area_correct_flow);
	readkeyprm_i(fp, "area_correct_vflow", &st->area_correct_vflow);

	for (i = 0; i < st->nbox; i++) {
		// Initalise lookupdz
		for (k = 0; k < st->wcnz; k++)
			lookupdz[i][k] = 0;

		// No longer subtract 1 here.
		st->boxnz[i] = (int) (readinnz[i]);
		if (st->boxnz[i] < 0)
			st->boxnz[i] = 0;
		st->boundaries[i] = (int) (boundary[i]);
		st->boxscale[i] = boxscalar[i];

		if (verbose)
			printf(" %f", st->boxscale[i]);
	}

	if (verbose)
		printf("\n");

	/* Read number of destination cells */
	readkeyprm_i(fp, "ndest", &st->ndest);

	/* Read flow mising data value */
	readkeyprm_d(fp, "missing_data", &st->hd.missing_data);
	st->hd.missing_data = -st->hd.missing_data;

	/* Read number of destination cells */
	readkeyprm_i(fp, "n_inline", &st->n_inline);

	/* Read the list of files */
	readkeyprm_i(fp, "input_type", &st->input_style);

	if (!st->input_style && st->uniflow)
		warn("Do you really mean unidirectional flows for Al Herman-like gross flow input data?\n");
	else if (st->input_style == 1 && !st->uniflow)
		warn("Read-in of Jeff Dunn data type not set up for gross flow data as yet, recode\n");

	/**
	 * If this is a mexico data type then we need to parse the input dir name and
	 * build the fname array from the directory contents.
	 */
	if (st->input_style == col_data_type || st->input_style == new_col_data_type) {
		readkeyprm_s(fp, "HydroInputDir", st->hd.inputDir);

		/* Allocate memory for file names if necessary */
		if (!st->hd.fname)
			st->hd.fname = c_alloc2d(STSLEN, MAXFILES);

		/* Now parse that directory */
		st->hd.nfiles = parseInputDir(st->hd.inputDir, st->hd.fname, "flux.dat");
	} else {
		/* Read the number of files */
		readkeyprm_i(fp, "nhdfiles", &st->hd.nfiles);
		if (st->hd.nfiles < 1)
			quit("Must be at least 1 exchange input file (%d in %s) to convert\n", st->hd.nfiles, st->runprmIfname);

		/* Allocate memory for file names if necessary */
		if (!st->hd.fname)
			st->hd.fname = c_alloc2d(STSLEN, st->hd.nfiles);

		if (verbose)
			printf("Reading flow data filenames\n");
		for (i = 0; i < st->hd.nfiles; i++) {
			sprintf(key, "trans%d.name", i);
			readkeyprm_s(fp, key, st->hd.fname[i]);
		}
	}

	st->timecorrect_on = 0;
	st->minExyear = MAXINT;
	st->monthEx = 0;
	st->dayEx = 0;
	st->Exmin_id = 0;
	st->maxExyear = 0;
	st->mmonthEx = 0;
	st->mdayEx = 0;
	st->Exmax_id = 0;

	/* Units of scale identifier */
	readkeyprm_i(fp, "unit_type", &st->unit_type);

	printf("st->unit_type = %d\n", st->unit_type);

	/* Number of output files to create for the hydrodynamics */
	readkeyprm_i(fp, "numoutfile", &st->numoutfile);

	/* Estuarine flows */
	readkeyprm_i(fp, "num_estuaries", &st->num_estuaries);
	if (st->num_estuaries > 0)
		i = st->num_estuaries;
	else
		i = 1;
	riverfcenum = (double *) d_alloc1d(i);
	riverinflux = (double *) d_alloc1d(i);
	riveroutflux = (double *) d_alloc1d(i);
	st->estinflux = (double *) d_alloc1d(i);
	st->estoutflux = (double *) d_alloc1d(i);
	st->estnum = (int *) i_alloc1d(i);

	readkeyprm_darray(fp, riveridi, &riverfcenum, &i);
	readkeyprm_darray(fp, nriverii, &riverinflux, &i);
	readkeyprm_darray(fp, nriveroi, &riveroutflux, &i);

	for (i = 0; i < st->num_estuaries; i++) {
		st->estnum[i] = (int) (riverfcenum[i]);
		st->estinflux[i] = riverinflux[i];
		st->estoutflux[i] = riveroutflux[i];
	}

	/* Read the number of files */
	readkeyprm_i(fp, "nvhdfiles", &st->vhd.nfiles);
	if (st->hd.nfiles < 1)
		warn("Must be at least 1 vertical exchange input file (%d in %s) to convert\n", st->vhd.nfiles, st->runprmIfname);

	/* Allocate memory for file names if necessary */
	if (!st->vhd.fname)
		st->vhd.fname = c_alloc2d(STSLEN, st->vhd.nfiles);

	if (verbose)
		printf("Reading vertical flow data filenames\n");
	for (i = 0; i < st->vhd.nfiles; i++) {
		sprintf(key, "vtrans%d.name", i);
		readkeyprm_s(fp, key, st->vhd.fname[i]);
	}

	/* Read flow mising data value */
	readkeyprm_i(fp, "temp_missing_data", &st->ts.temp_missing_data);
	st->ts.temp_missing_data = -st->ts.temp_missing_data;
	readkeyprm_i(fp, "salt_missing_data", &st->ts.salt_missing_data);
	st->ts.salt_missing_data = -st->ts.salt_missing_data;
	readkeyprm_i(fp, "ph_missing_data", &st->ts.ph_missing_data);
	st->ts.ph_missing_data = -st->ts.ph_missing_data;

	if(st->input_style == NOAA_type_id){
		/* Get temperature and salinity initial conditions file */
		sprintf(key, "tsTEMPinitfname");
		readkeyprm_s(fp, key, st->tsTEMPinitfname);
		sprintf(key, "tsSALTinitfname");
		readkeyprm_s(fp, key, st->tsSALTinitfname);
	}

	/* Read the number of temperature and salinity files */

	/**
	 * If this is a mexico data type then we need to parse the input dir name and
	 * build the fname array from the directory contents.
	 */
	if (st->input_style == col_data_type || st->input_style == new_col_data_type) {

		/* Allocate memory for file names if necessary */
		if (!st->ts.tsfname)
			st->ts.tsfname = c_alloc2d(STSLEN, MAXFILES);

		readkeyprm_s(fp, "TempInputDir", st->ts.inputDir);

		/* Allocate memory for file names if necessary */
		if (!st->hd.fname)
			st->hd.fname = c_alloc2d(STSLEN, MAXFILES);

		/* Now parse that directory */
		if (verbose)
			printf("Reading temperature and salinity data filenames\n");
		st->ts.ntsfiles = parseInputDir(st->ts.inputDir, st->ts.tsfname, "avg.dat");

		st->tsflagflux = 0;
	} else {

		readkeyprm_i(fp, "tsflagflux", &st->tsflagflux);
		readkeyprm_i(fp, "ntsfiles", &st->ts.ntsfiles);
		if (st->ts.ntsfiles < 1)
			warn("No temperature or salinity data files entered for conversion in %s\n", st->runprmIfname);

		/* Allocate memory for file names if necessary */
		if (!st->ts.tsfname)
			st->ts.tsfname = c_alloc2d(STSLEN, st->ts.ntsfiles);

		/* Read the list of files */
		if (verbose)
			printf("Reading temperature and salinity data filenames\n");
		for (i = 0; i < st->ts.ntsfiles; i++) {
			sprintf(key, "tempsalt%d.name", i);
			readkeyprm_s(fp, key, st->ts.tsfname[i]);
		}
	}

	st->minTSyear = MAXINT;
	st->monthTS = 0;
	st->dayTS = 0;
	st->Exmin_id = 0;
	st->maxTSyear = 0;
	st->mmonthTS = 0;
	st->mdayTS = 0;
	st->TSmax_id = 0;
		/* Read in maximum line number for Al Herman style input files */
	switch (st->input_style) {
	case NOAA_type_id:
		/* If using Al Herman style data */
		readkeyprm_i(fp, "ndata", &st->ndata);
		readkeyprm_i(fp, "altfaces", &st->altfaces);
		readkeyprm_i(fp, "notbgmface", &st->notbgmface);
		readkeyprm_i(fp, "start_day", &st->start_day);
		readkeyprm_s(fp, "lookup_table", st->lookupIfname);



		/* Get final number of data points */
		st->ndata *= st->hd.nfiles;

		/* Iterate so can have value matching max number */
		st->altfaces++;
		break;
	case flat_type_id:
		/* Standard flat file hydro input data */
		if (st->wcnz == 1) {
			/* CARS seasonal data - one average year entered so ndata = 4 */
			st->nhdsteps = 4;
		} else if (st->wcnz == 4) {
			/* Original Jeff Dunn flat format with one entry per watercolumn
			 level - currently only coded for 4 watercolumn layers. Number
			 of data entries = number of files as assumes one timestep per file
			 */
			st->nhdsteps = st->hd.nfiles;
		} else {
			st->nhdsteps = 1; // Will fail check below anyway
		}

		st->ndata = st->nhdsteps;
		st->altfaces = 0;
		st->notbgmface = -99999;
		st->start_day = 0;

		strcpy(st->lookupIfname, "");
		//sprintf(st->lookupIfname, "");

		/* Check for code restrictions to do with number of water column layers when not using Al Herman type data entry */
		if (st->wcnz != 1 && st->wcnz != 4)
			quit("Code currently assumes 1 OR 4 water column layers, recode otherwise\n");

		break;
	case cdf_type_id:
		/* Standard netcdf file hydro input data */
		get_netcdf_dimensions(st);
		st->notbgmface = -99999;
		st->start_day = 0;
		//readkeyprm_s(fp, "lookup_table", st->lookupIfname);

		break;
	case col_face_type:
		/* Excel created flat flow file - originally for the clarence river */
		st->nhdsteps = 0;
		for (i = 0; i < st->hd.nfiles; i++) {
			/* Open the exchange file */
			if ((fp2 = fopen(st->hd.fname[i], "r")) == NULL)
				quit("hydro_exchange_dataread: Can't open %s\n", st->hd.fname[i]);
			/* rewind file */
			fseek(fp2, 0L, SEEK_SET);
			sprintf(key, "nline");
			skipToKeyEnd(fp2, key);
			if (fscanf(fp2, "%d", &nline) != 1)
				quit("the number of lines must be marked at the head of the file as 'nline x' where x is the number of lines of data to be read in\n");
			st->nhdsteps += (nline + 1);
			/* Close the exchange file */
			fclose(fp2);
		}
		st->ndata = st->nhdsteps;
		st->altfaces = 0;
		st->notbgmface = -99999;
		st->start_day = 0;
		strcpy(st->lookupIfname, "");
		break;
	case ROMS_cdf_type:
		/* AL Herman cdf ROMS related files */
		get_ROMS_netcdf_dimensions(st);
		st->start_day = 0;
		strcpy(st->lookupIfname, "");
		break;
	case allinROMS_cdf_type:
		strcpy(st->timeVariableName, "ocean_time");
		strcpy(st->transportVariableName, "trans");
		strcpy(st->tempVariableName, "tempmean");
		strcpy(st->saltVariableName, "saltmean");
		strcpy(st->boxVariableName, "boxid");
		strcpy(st->zbinVariableName, "zbin");
		strcpy(st->secIDVariableName, "secid");

		/* AL Herman cdf ROMS related files */
		get_2in1ROMS_netcdf_dimensions(st);

		st->start_day = 0;
		strcpy(st->lookupIfname, "");
		break;
	case newROMS_cdf_type:
		strcpy(st->timeVariableName, "TAX");
		strcpy(st->transportVariableName, "TRANS");
		strcpy(st->tempVariableName, "TEMPMEAN");
		strcpy(st->saltVariableName, "SALTMEAN");
		strcpy(st->boxVariableName, "BOXID");
		strcpy(st->zbinVariableName, "ZBIN");
		strcpy(st->secIDVariableName, "SECID");
		strcpy(st->vertExchangeVariableName, "W");


		/* AL Herman cdf ROMS related files */
		get_2in1ROMS_netcdf_dimensions(st);

		st->start_day = 0;
		strcpy(st->lookupIfname, "");
		break;

	case col_data_type:
	case new_col_data_type:

		/* Files from Mexican oceanographers for Cam */
		readkeyprm_i(fp, "ndata", &st->ndata);
		readkeyprm_i(fp, "notbgmface", &st->notbgmface);
		readkeyprm_i(fp, "start_day", &st->start_day);
		readkeyprm_s(fp, "lookup_table", st->lookupIfname);
		readkeyprm_s(fp, "layer_depth_table", st->layerDepthIfname);

		readkeyprm_d(fp, "data_dt", &st->datadt);
		if (st->datadt > 28)
			quit(
					"to date code assumes this data will be entered at less than monthly intervals if this is not the case (e.g. here have dt = %d) Get_Date() needs to be recoded\n",
					st->datadt);

		printf("st->hd.nfiles = %d\n", st->hd.nfiles);
		printf("st->datadt = %f\n", st->datadt);
		st->nhdsteps = (int) (st->hd.nfiles * 31 * 86400 / st->dt);
		printf("st->nhdsteps = %d\n", st->nhdsteps);

		break;
	case eccal_flat_type:
		/* This might need a rethink */
		st->nhdsteps = (int) (st->hd.nfiles * 31 * 86400 / st->dt);
		st->ndata = st->nhdsteps;
		printf("st->nhdsteps = %d\n", st->nhdsteps);
		break;
	case guam_data_type:
		st->nhdsteps = st->nt;
		st->ndata = st->nt;

		readkeyprm_s(fp, "lookup_table", st->lookupIfname);

		printf("st->nhdstep = %d\n", st->nhdsteps);
		break;
	default:
		quit("No code for this input_style %d option as yet, try 0 (NOAA type), 1 (flat CSIRO format) or 2 (netcdf CSIRO format)\n", st->input_style);
	}

	/* Close the parameter file */
	fclose(fp);

	/* Free local vectors */
	d_free1d(readinnz);
	d_free1d(boundary);
	d_free1d(boxscalar);
	d_free1d(riverinflux);
	d_free1d(riveroutflux);
	d_free1d(riverfcenum);

	return;
}

/***************************************************************************************************
 Routine to get the exchange data out of the raw data file
 ****************************************************************************************************/
void get_hydro(Stuff *st) {
	char key[STSLEN];

	sprintf(key, "%d-01-01 00:00:00", st->REFyear);
	reference_time = GetTime(key);

	warn("Dates re-referenced with %d as year 0\n", st->REFyear);

	if (verbose)
		printf("Read in raw data\n");

	/* Call apropriate hydro and temperature and salinity data read-in */
	switch (st->input_style) {
	case NOAA_type_id:
		/* If using Al Herman style data */
		get_alt_data(st);
		get_alt_vert_data(st);
		get_alt_temp_salt(st);
		break;
	case flat_type_id:
		/* Standard flat file hydro input data */
		get_flat_hydro(st);
		get_flat_temp_and_salt(st);
		break;
	case cdf_type_id:
		/* Load the bgm mappings */
		load_bgm_mappings(st);
		/* Standard netcdf file hydro input data */
		get_cdf_hydro(st);
		get_cdf_temp_and_salt(st);
		break;
	case col_face_type:
		/* Flat form of CSIRO data originally developed for flows estimated in excel */
		get_flat_hydro(st);
		break;
	case ROMS_cdf_type:
		/* AL Herman cdf ROMS related files - reading in cdf files */
		get_ROMS_cdf_hydro(st);
		get_ROMS_cdf_temp_and_salt(st);
		break;
	case allinROMS_cdf_type:
		/* AL Herman cdf ROMS related files which contain all info in one file - reading in cdf files */
		get_2in1ROMS_cdf_hydro(st);
		break;
	case newROMS_cdf_type:
		get_NEW2in1ROMS_cdf_hydro(st);

		/* This is a fudge */
		st->numSubBoxes = 5;
		st->subBoxInfo = (SubsBoxInfo *) malloc(sizeof(SubsBoxInfo) * st->numSubBoxes);

		st->subBoxInfo[0].originalDestBox = 25;
		st->subBoxInfo[0].originalLeaveBox = 24;
		st->subBoxInfo[0].targetDestBox = 18;
		st->subBoxInfo[0].targetLeaveBox = 17;

		st->subBoxInfo[1].originalDestBox = 54;
		st->subBoxInfo[1].originalLeaveBox = 53;
		st->subBoxInfo[1].targetDestBox = 61;
		st->subBoxInfo[1].targetLeaveBox = 60;

		st->subBoxInfo[2].originalDestBox = 72;
		st->subBoxInfo[2].originalLeaveBox = 73;
		st->subBoxInfo[2].targetDestBox = 65;
		st->subBoxInfo[2].targetLeaveBox = 66;

		st->subBoxInfo[3].originalDestBox = 73;
		st->subBoxInfo[3].originalLeaveBox = 74;
		st->subBoxInfo[3].targetDestBox = 66;
		st->subBoxInfo[3].targetLeaveBox = 67;

		st->subBoxInfo[4].originalDestBox = 74;
		st->subBoxInfo[4].originalLeaveBox = 75;
		st->subBoxInfo[4].targetDestBox = 67;
		st->subBoxInfo[4].targetLeaveBox = 68;

		break;
	case col_data_type:
	case new_col_data_type:
		/* Files from Mexican oceanographers for Cam - reading in flat files */
		get_flat_hydro(st);
		get_flat_temp_and_salt(st);
		break;
	case guam_data_type:
			Load_Face_Lookup_Table(st);
			/* Files from Mexican oceanographers for Cam - reading in flat files */
			get_guam_hydro(st);
			get_guam_temp_and_salt(st);
			break;
	case eccal_flat_type:
		get_Eccal_flat_hydro(st);
		break;
	default:
		quit("No code for this input_style option as yet, try 0 (NOAA type), 1 (flat CSIRO format) or 2 (netcdf CSIRO format)\n");
	}

	return;

}

/***************************************************************************************************
 Routine to convert time in ISO date/time format (YYYY-MM-DD HH:MM:SS) to seconds
 ****************************************************************************************************/

double GetTime(char *d) {
	double j;
	double t;

	/* Convert date to Julian days */
	j = DateToJulian(d);
	/* Calculate time in seconds */
	t = 86400.0 * j;

	return (t);
}

/***************************************************************************************************
 * \brief Routine to find the input files in the given input dir.
 *
 *	This can be used for the hydro and temp_salt input files.
 *	The actually implementation of dirwalk that is called depends on the OS.
 *
 *
 ****************************************************************************************************/
int parseInputDir(char *dirName, char **fileNames, char *filePart) {
	int fileCount = 0;

	dirwalk(dirName, fileNames, &fileCount, filePart);

	return fileCount;
}

#ifdef WIN32

/**
 * \brief The Windows version of the dirwalk function.
 *
 *	This will recursively traverse the input directory.
 *	When a file is found its name is stored in the dirName array and the index value is
 *	incremented.
 *	When a directory is found dirwalk is called to traverse the files in that dir.
 *
 *
 */
void dirwalk(char *_current, char **dirNames, int *index, char *filePart)
{
	char DirName[MAX_PATH];
	HANDLE Hnd;
	WIN32_FIND_DATA WFD;
	static int _dlevel = 0;

	/*  Set the new current directory */
	SetCurrentDirectory( _current );

	/* Starts the search */
	Hnd = FindFirstFile( "*.*", &WFD );

	//  loop to get all inside the current directory
	while ( FindNextFile( Hnd, &WFD ) )
	{

		if (strcmp(WFD.cFileName, "..") && strcmp(WFD.cFileName, ".") )
		{
			//Get the current directory
			GetCurrentDirectory( MAX_PATH, DirName );

			if(WFD.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)
			{
				//       Put a "\" if necessary
				if ( strncmp( &DirName[strlen(DirName)-1], "\\", 1 ) )
				(void) strcat( DirName, "\\" );

				//       Create a new path
				(void) strcat( DirName, WFD.cFileName );

				//       Add 1 to level counter
				_dlevel++;

				//       Make a new call to itself
				dirwalk( DirName, dirNames, index, filePart);

				//       Go back one level
				SetCurrentDirectory( ".." );
			} else {
				if(strstr(WFD.cFileName, filePart) != NULL){
					//print the file name
					sprintf(dirNames[*index], "%s\\%s", DirName, WFD.cFileName );
					printf("fileName = %s\n", dirNames[*index]);
					(*index)++;
				}
			}
		}
	}
	// End the search to this call
	(void) FindClose( Hnd );
	return;
}
#else
/* dirwalk: apply fcn to all files in dir */
void dirwalk(char *dir, char **dirNames, int *index, char *filePart) {
	char name[MAX_PATH];
	Dirent *dp;
	DIR *dfd;

	struct stat stbuf;

	if ((dfd = opendir(dir)) == NULL) {
		quit("dirwalk: can 't open %s\n", dir);
		return;
	}

	while ((dp = readdir(dfd)) != NULL) {
		if (strcmp(dp -> d_name, ".") == 0 || strcmp(dp -> d_name, "..") == 0)
			continue; /* skip self and parent directory */

		if (strlen(dir) + strlen(dp -> d_name) + 2 > sizeof(name))
			fprintf(stderr, "dirwalk: name %s/%s too long\n", dir, dp->d_name);
		else {
			sprintf(name, "%s/%s", dir, dp -> d_name);

			if (stat(name, &stbuf) == -1) {
				fprintf(stderr, "fsize: can 't access %s\n", name);
				return;
			}

			if ((stbuf.st_mode & S_IFMT) == S_IFDIR)
				dirwalk(name, dirNames, index, filePart);
			else {

				if(strstr(name, filePart) != NULL){
					//print the file name
					sprintf(dirNames[*index], "%s", name);
					printf("fileName = %s\n", dirNames[*index]);
					(*index)++;
				}
			}

		} /* if-then-else */
	} /* while loop   */

	closedir(dfd);
}

#endif
