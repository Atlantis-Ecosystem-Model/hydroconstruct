/**********************************************************************

 get2in1ROMScdfinputdata.c

 by:	Beth Fulton
 date:	2/4/2009

 comments: This file is the code that reads in netcdf file version of
 raw hydrodynamics flows

 assumptions: assumes raw input is a series of data entries in netcdf file


 Changes:

 15-04-2009 Bec Gorton
 Changed the valboxs array to be a FPTYPE instead of an int array (netcdf lib was spitting when run with valgrind).
 Also fixed a number of other small array indexing bugs.

 27-04-2009 Bec Gorton
 Got rid of some incorrect code that was setting the exchange count array vlaues to all 1 - so no exchange values
 were being read in.

 28-04-2009 Bec Gorton
 Changed the code that sets the rawexchange values to use += instead of = as the fluxes loaded from these input
 files are not net values.

 02-06-20009 Bec Gorton
 Changed the error message displayed when the input file can't be read to be more descriptive.

 26-06-2009 Bec Gorton
 Added support for the new input format - its very similar the only changes are:
 - the variables have different names
 - the transport, salt and temp values are all floats so i need to be careful about the
 alloc functions.
 - The transport variable dimensions are different to the old format.
 - old format: double trans(ocean_time, boxid, zbin, secid) ;
 - new format: float TRANS(TAX, ZBIN, BOXID, SECID) ;

 - The temp variable is also different:
 - old format double tempmean(ocean_time, boxid, zbin)
 - new format float TEMPMEAN(TAX, ZBIN, BOXID) ;


 14-09-2009 Bec Gorton
 Added support to read in the vertical flux values.
 Also change the code to correctly match the values read in from the netcdf files
 with the faces in the bgm file.

 24-09-2009 Bec Gorton
 Added code to correct for the direction of the flow in the input file.
 Flows should be generated from left to right across a face - the input file
 just has flows out of a box. So the new scale value will multiply the exchange values
 by -1 if the direction is wrong.
 *************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <math.h>
#include <string.h>
#include <netcdf.h>
#include <sjwlib.h>
#include "mainexchange.h"

/***************************************************************************************************
 Routine to get the exchange data out of the raw data file - Al Herman ROMS data.
 Involves reading data from netcdf file to buffer
 ****************************************************************************************************/
void get_2in1ROMS_cdf_hydro(Stuff *st) {
	int i, k, lastk, entry_i, nz, ij, final_fce, b;
	double time_of_year, totnum, this_tofy, mintime;
	FPTYPE ****valexchange;
	FPTYPE ***valtemp;
	FPTYPE ***valsalt;
	FPTYPE *valtime;
	FPTYPE *valboxs;
	long start[4];
	long count[4];
	long startts[3];
	long countts[3];
	long startt[1];
	long countt[1];
	long startbx[1];
	long countbx[1];
	long nt, fce, lvl, bx = 0;

	ncopts = NC_VERBOSE | NC_FATAL;

	if (verbose)
		printf("Loading raw data from netcdf files\n");

	/* Give code warning */
	if (st->hd.nfiles > 1)
		warn("%d raw hydro files listed, code may not cope properly with more than one - in particular the face mapping may be wrong\n", st->hd.nfiles);

	/* Loop through the exchange files */
	lastk = 0;
	for (i = 0; i < st->hd.nfiles; i++) {
		if (verbose)
			printf("Reading raw exchange file %s\n", st->hd.fname[i]);

		/* Open the exchange file */
		open_2in1ROMS_hydro(st, st->hd.fname[i], i);

		if (verbose > 1)
			printf("Checking for dimensions\n");

		/* Check for dimensions in file */

		ncdiminq(st->hd.fid, ncdimid(st->hd.fid, st->timeVariableName), NULL, &nt);
		ncdiminq(st->hd.fid, ncdimid(st->hd.fid, st->boxVariableName), NULL, &bx);
		ncdiminq(st->hd.fid, ncdimid(st->hd.fid, st->zbinVariableName), NULL, &lvl);
		ncdiminq(st->hd.fid, ncdimid(st->hd.fid, st->secIDVariableName), NULL, &fce);

		if (fce > st->max_nconn) {
			warn("Less faces in netcdf file %s (nface = %d) than in bgm file %s (nface = %d)\n", st->hd.fname[i], fce, st->geomIfname, st->nface);
		}

		/* Allocate temporary storage for one tracer */
		switch (sizeof(FPTYPE)) {
		case sizeof(float):
			valexchange = (FPTYPE ****) f_alloc4d(fce, lvl, bx, nt);
			valtemp = (FPTYPE ***) f_alloc3d(lvl, bx, nt);
			valsalt = (FPTYPE ***) f_alloc3d(lvl, bx, nt);
			valtime = (FPTYPE *) f_alloc1d(nt);
			valboxs = (FPTYPE *) f_alloc1d(bx);
			break;
		case sizeof(double):
			valexchange = (FPTYPE ****) d_alloc4d(fce, lvl, bx, nt);
			valtemp = (FPTYPE ***) d_alloc3d(lvl, bx, nt);
			valsalt = (FPTYPE ***) d_alloc3d(lvl, bx, nt);
			valtime = (FPTYPE *) d_alloc1d(nt);
			valboxs = (FPTYPE *) d_alloc1d(bx);
			break;
		default:
			quit("get_cdf_hydro: Unknown size for FPTYPE\n");
		}

		/* Set indices for reading tracers */
		if (verbose > 1)
			printf("Set start and count\n");

		start[0] = 0;
		start[1] = 0;
		start[2] = 0;
		start[3] = 0;

		count[0] = nt;
		count[1] = bx;
		count[2] = lvl;
		count[3] = fce;

		if (verbose > 1)
			printf("Set startts and countts\n");

		startts[0] = 0;
		startts[1] = 0;
		startts[2] = 0;

		countts[0] = nt;
		countts[1] = bx;
		countts[2] = lvl;

		startbx[0] = 0;
		countbx[0] = bx;

		if (verbose > 1)
			printf("Set startbx and countbx\n");

		startt[0] = 0;
		countt[0] = nt;

		if (verbose > 1)
			printf("Set startt and countt\n");

		/* Read data - exchanges, dates and box-ids associated with the faces */

		if (verbose > 1)
			printf("Get valtime (id: %d)\n", st->hd.t_vid);
		ncvarget(st->hd.fid, st->hd.t_vid, startt, countt, valtime);

		if (verbose > 1)
			printf("Get valbox (id : %d)\n", st->hd.srcb_vid);
		ncvarget(st->hd.fid, st->hd.srcb_vid, startbx, countbx, valboxs);

		if (verbose > 1)
			printf("Get valtemp (id: %d)\n", st->ts.temp_vid);
		ncvarget(st->hd.fid, st->ts.temp_vid, startts, countts, valtemp[0][0]);

		if (verbose > 1)
			printf("Get valsalt (id : %d)\n", st->ts.salt_vid);
		ncvarget(st->hd.fid, st->ts.salt_vid, startts, countts, valsalt[0][0]);

		if (verbose > 1)
			printf("Get valexchange (id: %d)\n", st->hd.e_vid);
		ncvarget(st->hd.fid, st->hd.e_vid, start, count, valexchange[0][0][0]);

		/** Transfer to rawdata arrays **/
		/* Transfer boxes and faces */
		for (b = 0; b < bx; b++) {
			for (k = 0; k < fce; k++) {
				tempfaceinfo[(int) valboxs[b]][k][fceid_id] = k;
			}
		}

		/* Transfer timestamps - assumes this is in seconds */
		if (st->reset_time) {
			mintime = MAXDOUBLE;
			for (k = 0; k < nt; k++) {
				if (valtime[k] < mintime)
					mintime = valtime[k];
			}
			st->time_reset = mintime * st->dt;
		}

		for (k = 0; k < nt; k++) {
			exchangedate[k + lastk][time_id] = valtime[k] * st->dt - st->time_reset;

			/* Multiple by 86400 / dt as may be tidal entries so must expand
			 list to fit these extra entries in */
			time_of_year = valtime[k] * (86400 / st->dt);
			totnum = 365.0 * (86400 / st->dt);

			if (time_of_year / totnum < 1.0)
				this_tofy = floor(time_of_year / (86400 / st->dt) + 0.5);
			else
				this_tofy = floor((time_of_year / totnum - floor(time_of_year / totnum)) * totnum + 0.5);

			exchangedate[k + lastk][TofY_id] = this_tofy;

			if (st->verbose_pad)
				printf("k: %d, time: %e, tofy: %e valtime: %e\n", k, exchangedate[k + lastk][time_id], exchangedate[k + lastk][TofY_id], valtime[k]);

			/* Save time and date for temp-salt data too */
			tempsaltdate[k + lastk][TofY_id] = exchangedate[k + lastk][TofY_id];
			tempsaltdate[k + lastk][time_id] = exchangedate[k + lastk][time_id];
		}


		/* Transfer exchanges */
		entry_i = 0;
		printf("st->hd.missing_data = %e\n", st->hd.missing_data);

		for (k = 0; k < (int) nt; k++) {
			for (b = 0; b < bx; b++) {
				for (nz = 0; nz < lvl; nz++) {
					for (ij = 0; ij < st->boxnconn[b]; ij++) {
						final_fce = tempfaceinfo[b][ij][rbox_id]; // Actually this is a face ID (overloading rbox_id as just a placeholder index

						if (final_fce < 0 || final_fce > st->nface)
							quit("Trying to assign an exchange to a non-existent face - secid: %d, box: %d, final_fce: %d\n", ij, b, final_fce);

						if (valexchange[k][b][nz][ij] <= (st->hd.missing_data + 1)) {
							/* DO nothing as data missing */
							if (!st->sdiff)
								rawexchanges[k + lastk][nz][final_fce][entry_i] += 0;
							else
								rawexchanges[k + lastk][nz][final_fce][entry_i] += 0.00000001;
						} else {
							rawexchanges[k + lastk][nz][final_fce][entry_i] += valexchange[k][b][nz][ij];
						}
					}

					if (valtemp[k][b][nz] <= ((float) (st->hd.missing_data + 1))) {
						/* DO nothing as data missing */
						if (verbose > 1)
							printf("missing temp data\n");
					} else
						rawtempsalt[k + lastk][nz][b][temp_id] = valtemp[k][b][nz];

					if (valsalt[k][b][nz] <= ((float) (st->hd.missing_data + 1))) {
						/* DO nothing as data missing */
						if (verbose > 1)
							printf("missing salt data\n");
					} else
						rawtempsalt[k + lastk][nz][b][salt_id] = valsalt[k][b][nz];

					if (verbose > 1)
						printf("Storing rawtemp[%d][%d][%d] = %e\n", k + lastk, nz, b, rawtempsalt[k + lastk][nz][b][temp_id]);

				}
			}
		}

		/* Set lastk */
		lastk += nt;

		/* Free temporary storage */
		switch (sizeof(FPTYPE)) {
		case sizeof(float):
			f_free4d((float ****) valexchange);
			f_free3d((float ***) valtemp);
			f_free3d((float ***) valsalt);
			f_free1d((float *) valtime);
			break;
		case sizeof(double):
			d_free4d((double ****) valexchange);
			d_free3d((double ***) valtemp);
			d_free3d((double ***) valsalt);
			d_free1d((double *) valtime);
			break;
		default:
			quit("get_cdf_hydro: Unknown size for FPTYPE\n");
		}
		i_free1d((int *) valboxs);

		/* Close hydro file */
		close_hydro(st);
	}

	/* Map netcdf box assignment to bgm assignment */
	map_ROMS_netcdf_faces(st);

	if (verbose)
		printf("Finished reading raw exchange files\n");

	return;
}

/***************************************************************************************************
 Routine to get the exchange data out of the raw data file - Al Herman ROMS data.
 Involves reading data from netcdf file to buffer
 ****************************************************************************************************/
void get_NEW2in1ROMS_cdf_hydro(Stuff *st) {
	int i, k, lastk, entry_i, nz, ij, final_fce, b;
	double time_of_year, totnum, this_tofy, mintime, dtmult, scale = 0.0;
	float ****valexchange;
	float ***vertexchange;
	float ***valtemp;
	float ***valsalt;
	FPTYPE *valtime;
	FPTYPE *valboxs;
	long start[4];
	long count[4];
	long startts[3];
	long countts[3];
	long startt[1];
	long countt[1];
	long startbx[1];
	long countbx[1];
	long nt, fce, lvl, bx = 0;
	int netcdfVertID;

	ncopts = NC_VERBOSE | NC_FATAL;

	if (verbose)
		printf("Loading raw data from netcdf files\n");

	/* Give code warning */
	if (st->hd.nfiles > 1)
		warn("%d raw hydro files listed, code may not cope properly with more than one - in particular the face mapping may be wrong\n", st->hd.nfiles);

	/* Loop through the exchange files */
	lastk = 0;
	for (i = 0; i < st->hd.nfiles; i++) {
		if (verbose)
			printf("Reading raw exchange file %s\n", st->hd.fname[i]);

		/* Open the exchange file */
		open_2in1ROMS_hydro(st, st->hd.fname[i], i);

		if (verbose > 1)
			printf("Checking for dimensions\n");

		/* Check for dimensions in file */
		ncdiminq(st->hd.fid, ncdimid(st->hd.fid, st->timeVariableName), NULL, &nt);
		ncdiminq(st->hd.fid, ncdimid(st->hd.fid, st->boxVariableName), NULL, &bx);
		ncdiminq(st->hd.fid, ncdimid(st->hd.fid, st->zbinVariableName), NULL, &lvl);
		ncdiminq(st->hd.fid, ncdimid(st->hd.fid, st->secIDVariableName), NULL, &fce);

		if (fce > st->max_nconn) {
			warn("Less faces in netcdf file %s (nface = %d) than in bgm file %s (nface = %d)\n", st->hd.fname[i], fce, st->geomIfname, st->nface);
		}

		valexchange = (float ****) f_alloc4d(fce, bx, lvl, nt);
		vertexchange = (float ***) f_alloc3d(bx, lvl, nt);
		valtemp = (float ***) f_alloc3d(bx, lvl, nt);
		valsalt = (float ***) f_alloc3d(bx, lvl, nt);

		/* Allocate temporary storage for one tracer */
		switch (sizeof(FPTYPE)) {
		case sizeof(float):
			valtime = (FPTYPE *) f_alloc1d(nt);
			valboxs = (FPTYPE *) f_alloc1d(bx);
			break;
		case sizeof(double):
			valtime = (FPTYPE *) d_alloc1d(nt);
			valboxs = (FPTYPE *) d_alloc1d(bx);
			break;
		default:
			quit("get_cdf_hydro: Unknown size for FPTYPE\n");
		}

		/* Set indices for reading tracers */
		if (verbose > 1)
			printf("Set start and count\n");

		start[0] = 0;
		start[1] = 0;
		start[2] = 0;
		start[3] = 0;

		count[0] = nt;
		count[1] = lvl;
		count[2] = bx;
		count[3] = fce;

		if (verbose > 1)
			printf("Set startts and countts\n");

		startts[0] = 0;
		startts[1] = 0;
		startts[2] = 0;

		countts[0] = nt;
		countts[1] = lvl;
		countts[2] = bx;

		startbx[0] = 0;
		countbx[0] = bx;

		if (verbose > 1)
			printf("Set startbx and countbx\n");

		startt[0] = 0;
		countt[0] = nt;

		if (verbose > 1)
			printf("Set startt and countt\n");

		/* Read data - exchanges, dates and box-ids associated with the faces */

		if (verbose > 1)
			printf("Get valtime (id: %d)\n", st->hd.t_vid);
		ncvarget(st->hd.fid, st->hd.t_vid, startt, countt, valtime);

		if (verbose > 1)
			printf("Get valbox (id : %d)\n", st->hd.srcb_vid);
		ncvarget(st->hd.fid, st->hd.srcb_vid, startbx, countbx, valboxs);

		if (verbose > 1)
			printf("Get valtemp (id: %d)\n", st->ts.temp_vid);
		ncvarget(st->hd.fid, st->ts.temp_vid, startts, countts, valtemp[0][0]);

		if (verbose > 1)
			printf("Get valsalt (id : %d)\n", st->ts.salt_vid);
		ncvarget(st->hd.fid, st->ts.salt_vid, startts, countts, valsalt[0][0]);

		if (verbose > 1)
			printf("Get vertexchange (id : %d)\n", st->ts.ve_vid);
		ncvarget(st->hd.fid, st->ts.ve_vid, startts, countts, vertexchange[0][0]);


		if (verbose > 1)
			printf("Get valexchange (id: %d)\n", st->hd.e_vid);
		ncvarget(st->hd.fid, st->hd.e_vid, start, count, valexchange[0][0][0]);

		/** Transfer to rawdata arrays **/
		/* Transfer boxes and faces */
		for (b = 0; b < bx; b++) {
			for (k = 0; k < fce; k++) {
				tempfaceinfo[(int) valboxs[b]][k][fceid_id] = k;
			}
		}

		/* Transfer timestamps - assumes this is in seconds */
		dtmult = 3600; // As Al Herman data in hours not seconds
		if (st->reset_time) {
			mintime = MAXDOUBLE;
			for (k = 0; k < nt; k++) {
				if (valtime[k] < mintime)
					mintime = valtime[k];
			}
			st->time_reset = mintime * dtmult;
		}

		printf("st->time_reset = %e\n", st->time_reset);
		for (k = 0; k < nt; k++) {
			exchangedate[k + lastk][time_id] = valtime[k] * dtmult - st->time_reset;

			/* Multiple by 86400 / dt as may be tidal entries so must expand
			 list to fit these extra entries in */
			//			time_of_year = valtime[k] *  (86400 / st->dt);
			//			totnum = 365.0 * (86400 / st->dt);
			//
			//			if( time_of_year/totnum < 1.0 )
			//				this_tofy = floor(time_of_year/(86400 / st->dt) + 0.5);
			//			else
			//				this_tofy = floor((time_of_year/totnum - floor(time_of_year/totnum))*totnum + 0.5);

			/* Add 15 days to allow for the data started on the 15th day of this year - should perhaps add this as a variable to the input file.*/
			time_of_year = (valtime[k] + 15* 24 ) * (dtmult / st->dt);

			totnum = 365.0 * (86400 / st->dt);

			if (time_of_year / totnum < 1.0)
				this_tofy = floor(time_of_year / (86400 / st->dt) + 0.5);
			else
				this_tofy = floor((time_of_year / totnum - floor(time_of_year / totnum)) * totnum + 0.5);

			exchangedate[k + lastk][TofY_id] = this_tofy;

			if (st->verbose_pad)
				printf("k: %d, time: %e, tofy: %e valtime: %e\n", k, exchangedate[k + lastk][time_id], exchangedate[k + lastk][TofY_id], valtime[k]);

			/* Save time and date for temp-salt data too */
			tempsaltdate[k + lastk][TofY_id] = exchangedate[k + lastk][TofY_id];
			tempsaltdate[k + lastk][time_id] = exchangedate[k + lastk][time_id];
		}


		/* Transfer exchanges */
		entry_i = 0;
		printf("st->hd.missing_data = %e\n", st->hd.missing_data);

		for (k = 0; k < (int) nt; k++) {
			for (b = 0; b < bx; b++) {
				if (verbose > 1)
					printf("Box %d\n", b);
				for (nz = 0; nz < lvl; nz++) {

					/* Need to read in the actual internal boundary faces */
					for (ij = 0; ij < st->boxnconn[b]; ij++) {

						final_fce = tempfaceinfo[b][ij][rbox_id]; // Actually this is a face ID (overloading rbox_id as just a placeholder index

						/* Get the index of this face in the input data  - there is a value in the input file for all
						 * faces - not just those that are inportant like in the bgm file
						 */
						netcdfVertID = st->faceVertIndex[final_fce];

						if (netcdfVertID == -1)
							quit("Could not match the netcdf connection to the bgm connection. Box = %d, nconn = %d\n", b, ij);

						if (final_fce < 0 || final_fce > st->nface)
							quit("Trying to assign an exchange to a non-existent face - secid: %d, box: %d, final_fce: %d\n", ij, b, final_fce);

						if (verbose > 1)
							printf("nconn = %d, face = %d, netcdfVertID = %d\n", ij, final_fce, netcdfVertID);

						/* Need to check that the direction of the flow is correct.
						 * The data read in is positive into the box but we want positive values from
						 * left to right accross the face.
						 * So need to work out if 'into the box' is the same as left to right.
						 */
						if(st->faces[final_fce].ibl == b)
							scale = 1.0;
						else
							scale = -1.0;

						if (valexchange[k][nz][b][netcdfVertID] <= ((float) (st->hd.missing_data + 1))) {
							/* DO nothing as data missing */
							if (!st->sdiff)
								rawexchanges[k + lastk][nz][final_fce][entry_i] += 0;
							else
								rawexchanges[k + lastk][nz][final_fce][entry_i] += 0.00000001;
						} else {
							rawexchanges[k + lastk][nz][final_fce][entry_i] += valexchange[k][nz][b][netcdfVertID] * scale;
						}

						if (verbose > 1)
							printf("Storing rawexchanges[%d][%d][%d][%d] = %e\n", k + lastk, nz, final_fce, entry_i,
									rawexchanges[k + lastk][nz][final_fce][entry_i]);

					}

					if (vertexchange[k][nz][b] <= ((float) (st->hd.missing_data + 1))) {
						/* DO nothing as data missing */
						if (verbose > 1)
							printf("missing vertexchange data\n");
					} else
						rawvertexchanges[k + lastk][nz][b][entry_i] = vertexchange[k][nz][b];

					if (valtemp[k][nz][b] <= ((float) (st->hd.missing_data + 1))) {
						/* DO nothing as data missing */
						if (verbose > 1)
							printf("missing temp data\n");
					} else
						rawtempsalt[k + lastk][nz][b][temp_id] = valtemp[k][nz][b];

					if (valsalt[k][nz][b] <= ((float) (st->hd.missing_data + 1))) {
						/* DO nothing as data missing */
						if (verbose > 1)
							printf("missing salt data\n");
					} else {
						rawtempsalt[k + lastk][nz][b][salt_id] = valsalt[k][nz][b];
					}

					if (verbose > 1) {
						printf("Storing rawtemp[%d][%d][%d] = %e\n", k + lastk, nz, b, rawtempsalt[k + lastk][nz][b][temp_id]);
						printf("Storing rawsalt[%d][%d][%d] = %e\n", k + lastk, nz, b, rawtempsalt[k + lastk][nz][b][salt_id]);
					}

				}
			}
		}

		/* Set lastk */
		lastk += nt;

		/* Free temporary storage */
		f_free4d((float ****) valexchange);
		f_free3d((float ***) vertexchange);
		f_free3d((float ***) valtemp);
		f_free3d((float ***) valsalt);

		switch (sizeof(FPTYPE)) {
		case sizeof(float):
			f_free1d((float *) valtime);
			f_free1d((float *) valboxs);
			break;
		case sizeof(double):
			d_free1d((double *) valtime);
			d_free1d((double *) valboxs);
			break;
		default:
			quit("get_cdf_hydro: Unknown size for FPTYPE\n");
		}

		/* Close hydro file */
		close_hydro(st);
	}

	/* Map netcdf box assignment to bgm assignment */
	map_ROMS_netcdf_faces(st);

	if (verbose)
		printf("Finished reading raw exchange files\n");

	return;
}

/****************************************************************************************************
 Routine to open netcdf file containing raw hydrodynamic fluxes
 *****************************************************************************************************/
void open_2in1ROMS_hydro(Stuff *st, char *name, int filenum) {
	int ndims;
	int nvars;
	int natts;
	int recdim;
	long n;
	nc_type daty;
	int dims[MAX_NC_DIMS];
	char hdu[STSLEN];
	char stu[STSLEN];

	/* Set netCDF library error handling */
	ncopts = NC_VERBOSE;

	if (verbose > 1)
		printf("open_2in1ROMS_hydro: opening the file\n");

	/* Open the file */
	if ((st->hd.fid = ncopen(name, NC_NOWRITE)) < 0)
		quit("open_2in1ROMS_hydro: Can't open hydrodynamic model input data file %s\n", name);

	/* Inquire about this file */
	ncopts = NC_VERBOSE | NC_FATAL;
	ncinquire(st->hd.fid, &ndims, &nvars, &natts, &recdim);
	if (ndims < 3)
		quit("open_hydro: not enough dimensions in %s\n", name);

	if (verbose > 1)
		printf("open_hydro: check dimensions\n");

	/* Check dimensions are as expected */
	printf("check for time for 2-in-1 file\n");

	if ((st->hd.t_did = ncdimid(st->hd.fid, st->timeVariableName)) == -1)
		quit("open_hydro: no t dimension in %s\n", name);
	if (st->hd.t_did != recdim)
		warn("open_hydro: (recdim: %d, hd.t_did: %d) t dimension not unlimited in %s\n", recdim, st->hd.t_did, name);

	printf("check for boxes\n");

	if ((st->hd.b_did = ncdimid(st->hd.fid, st->boxVariableName)) == -1) {
		quit("open_hydro: no boxes dimension in %s\n", name);
	}

	printf("check for level\n");

	if ((st->hd.z_did = ncdimid(st->hd.fid, st->zbinVariableName)) == -1)
		quit("open_hydro: no levels dimension in %s\n", name);

	printf("check for faces\n");

	if ((st->hd.f_did = ncdimid(st->hd.fid, st->secIDVariableName)) == -1) {
		quit("open_hydro: no faces dimension in %s\n", name);
	}

	/* Get dimension sizes. Note check against geometry now happens in get_cdf_hydro() */

	if (verbose > 1)
		printf("open_hydro: get dimensions.\n");

	ncdiminq(st->hd.fid, st->hd.b_did, NULL, &n);
	ncdiminq(st->hd.fid, st->hd.z_did, NULL, &n);

	/* Check that time units and steps match this model */
	if (verbose > 1)
		printf("open_hydro: get units\n");

	if (verbose > 1) {
		printf("time variable name = %s\n", st->timeVariableName);

	}
	// st->hd.t_vid = ncvarid(st->hd.fid,st->timeVariableName);
	st->hd.t_vid = ncvarid(st->hd.fid, "TAX");

	printf("st->hd.t_vid  = %d\n", st->hd.t_vid);
	if (st->hd.t_vid < 0)
		quit("open_hydro: no t variable in %s\n", name);

	memset(st->hd.t_units, 0, STSLEN);
	ncattget(st->hd.fid, st->hd.t_vid, "units", st->hd.t_units);
	sscanf(st->hd.t_units, "%s", hdu);
	if (filenum == 0)
		sprintf(st->t_units, "%s", st->hd.t_units);
	sscanf(st->t_units, "%s", stu);
	if (strcmp(hdu, stu) != 0) {
		quit("open_hydro: Time units (%s) in %s don't match model time units (%s) loaded from %s\n", hdu, name, stu, st->hd.fname[0]);
	}
	st->hd.dt = 604800.0;

	/* Find out how many time steps are in the file */
	ncdiminq(st->hd.fid, st->hd.t_did, NULL, &st->hd.nstep);

	if (verbose > 1) {
		printf("Box variable name = %s\n", st->boxVariableName);

	}

	/* Get other variable ids */
	st->hd.srcb_vid = ncvarid(st->hd.fid, st->boxVariableName);
	printf("st->hd.srcb_vid = %d\n", st->hd.srcb_vid);
	if (st->hd.srcb_vid < 0)
		quit("open_hydro: no %s variable in %s\n", st->boxVariableName, name);

	st->hd.e_vid = ncvarid(st->hd.fid, st->transportVariableName);
	if (st->hd.e_vid < 0)
		quit("open_hydro: no transport variable in %s\n", name);

	/* Check variable types and dimensions for the transport data */
	if (verbose > 1)
		printf("open_hydro: get exchange data ids\n");

	ncvarinq(st->hd.fid, st->hd.e_vid, NULL, &daty, &ndims, dims, &natts);
	if (nctypelen(daty) != sizeof(FPTYPE))
		warn("open_hydro: Type of exchange variable doesn't match model (ok if float in file and double FPTYPE)\n");

	if (st->input_style == newROMS_cdf_type) {
		if (ndims != 4 || dims[0] != st->hd.t_did || dims[1] != st->hd.z_did || dims[2] != st->hd.b_did || dims[3] != st->hd.f_did)
			quit("open_hydro: transport variable has incorrect dimensions\n");

	} else {
		if (ndims != 4 || dims[0] != st->hd.t_did || dims[1] != st->hd.b_did || dims[2] != st->hd.z_did || dims[3] != st->hd.f_did)
			quit("open_hydro: transport variable has incorrect dimensions\n");

	}

	ncattget(st->hd.fid, st->hd.e_vid, "units", st->hd.e_units);
	if (verbose > 1)
		printf("open_hydro: get exchange units\n");

	sscanf(st->hd.e_units, "%s", stu);
	if (strstr("Sverdrup", stu) != NULL)
		st->unit_type = 0;
	else
		st->unit_type = 1;

	/** Now do temperature and salinity information **/
	/* Inquire about this file */
	printf("open_tsfile: inquire about units\n");

	/* Set netCDF library error handling */
	ncopts = NC_VERBOSE;

	/* Open the file */
	if ((st->ts.tsfid = ncopen(name, NC_NOWRITE)) < 0)
		quit("open_2in1ROMS_hydro: Can't open temperature and salinity model input data file %s\n", name);

	/* Inquire about this file */
	ncopts = NC_VERBOSE | NC_FATAL;
	ncinquire(st->ts.tsfid, &ndims, &nvars, &natts, &recdim);
	if (ndims < 3)
		quit("open_tsfile: not enough dimensions in %s\n", name);

	/* Check dimensions are as expected */
	if ((st->ts.t_did = ncdimid(st->ts.tsfid, st->timeVariableName)) == -1)
		quit("open_tsfile: no t dimension in %s\n", name);
	if (st->ts.t_did != recdim)
		warn("open_tsfile: (recdim: %d, ts.t_did: %d) t dimension not unlimited in %s\n", recdim, st->ts.t_did, name);
	if ((st->ts.b_did = ncdimid(st->ts.tsfid, st->boxVariableName)) == -1)
		quit("open_tsfile: no boxes dimension in %s\n", name);
	if ((st->ts.z_did = ncdimid(st->ts.tsfid, st->zbinVariableName)) == -1)
		quit("open_tsfile: no levels dimension in %s\n", name);

	/* Get dimension sizes. Note check against geometry now happens in get_cdf_hydro() */
	ncdiminq(st->ts.tsfid, st->ts.b_did, NULL, &n);
	ncdiminq(st->ts.tsfid, st->ts.z_did, NULL, &n);

	/* Check that time units and steps match this model */
	st->ts.t_vid = ncvarid(st->ts.tsfid, st->timeVariableName);
	if (st->ts.t_vid < 0)
		quit("open_tsfile: no t variable in %s\n", name);
	memset(st->ts.t_units, 0, STSLEN);
	ncattget(st->ts.tsfid, st->ts.t_vid, "units", st->ts.t_units);
	sscanf(st->ts.t_units, "%s", hdu);
	if (!filenum)
		sprintf(st->t_units, "%s", st->ts.t_units);
	sscanf(st->t_units, "%s", stu);
	if (strcmp(hdu, stu) != 0) {
		quit("open_tsfile: Time units (%s) in %s don't match model time units (%s) loaded from %s\n", hdu, name, stu, st->ts.tsfname[0]);
	}

	/* Check for timestep size */
	st->ts.tsdt = 604800.0;

	/* Find out how many time steps are in the file */
	ncdiminq(st->ts.tsfid, st->ts.t_did, NULL, &st->ts.tsnstep);

	/* Get other variable ids */
	st->ts.b_vid = ncvarid(st->ts.tsfid, st->boxVariableName);
	if (st->ts.b_vid < 0)
		quit("open_tsfile: no boxes variable in %s\n", name);

	st->ts.temp_vid = ncvarid(st->ts.tsfid, st->tempVariableName);
	if (st->ts.temp_vid < 0)
		quit("open_tsfile: no temperature variable in %s\n", name);

	st->ts.salt_vid = ncvarid(st->ts.tsfid, st->saltVariableName);
	if (st->ts.salt_vid < 0)
		quit("open_tsfile: no salinity variable in %s\n", name);

	/* The vertical exchange */
	st->ts.ve_vid = ncvarid(st->ts.tsfid, st->vertExchangeVariableName);
	if( st->ts.ve_vid < 0 )
	warn("open_tsfile: no verticalflux variable in %s\n",name);


	/* Check variable types and dimensions for the vertical transport data */
	ncvarinq(st->ts.tsfid,st->ts.ve_vid,NULL,&daty,&ndims,dims,&natts);
	if( nctypelen(daty) != sizeof(FPTYPE) )
		warn("open_tsfile: Type of verticalflux variable doesn't match model (OK if float in file and double FPTYPE)\n");

	if (st->input_style == newROMS_cdf_type) {
		if (ndims != 3 || dims[0] != st->ts.t_did || dims[1] != st->ts.z_did || dims[2] != st->ts.b_did)
			quit("open_hydro: transport variable has incorrect dimensions\n");

	} else {
		if (ndims != 3 || dims[0] != st->ts.t_did || dims[1] != st->ts.b_did || dims[2] != st->ts.z_did)
			quit("open_hydro: verticalflux variable has incorrect dimensions\n");

	}

	/* Can't check the units as the file has no unit value */
//
//	ncattget(st->ts.tsfid,st->ts.ve_vid,"units",st->ts.ve_units);
//	sscanf(st->ts.ve_units,"%s",stu);
//	if(strstr("Sverdrup",stu) != NULL ){
//		if(st->unit_type)
//			quit("Flux data mismatch - horizontal fluxes in m3 and vertical in Sverdrup\n");
//	} else {
//		if(!st->unit_type)
//			quit("Flux data mismatch - horizontal fluxes in Sverdrup and vertical in m3\n");
//	}

	/* Reset netCDF error handling */
	ncopts = NC_VERBOSE | NC_FATAL;

	return;
}

/***************************************************************************************************
 Routine to calculate number of data points, faces and levels to be loaded from netcdf file
 ****************************************************************************************************/
void get_2in1ROMS_netcdf_dimensions(Stuff *st) {
	int i, totndata, maxfaces, minlevel, maxlevel, maxboxes;
	long nt, fce, lvl, bx;

	/* Initialise counters */
	totndata = 0;
	maxfaces = 0;
	minlevel = MAXINT;
	maxlevel = 0;
	maxboxes = 0;

	for (i = 0; i < st->hd.nfiles; i++) {
		if (verbose)
			printf("Reading raw exchange file to get dimensions %s\n", st->hd.fname[i]);

		/* Open the exchange file */
		open_2in1ROMS_hydro(st, st->hd.fname[i], i);

		/* Check for dimensions in file */

		printf("dimensions: getting %s\n", st->timeVariableName);
		ncdiminq(st->hd.fid, ncdimid(st->hd.fid, st->timeVariableName), NULL, &nt);

		printf("dimensions: getting %s\n", st->boxVariableName);
		ncdiminq(st->hd.fid, ncdimid(st->hd.fid, st->boxVariableName), NULL, &bx);

		printf("dimensions: getting %s\n", st->zbinVariableName);
		ncdiminq(st->hd.fid, ncdimid(st->hd.fid, st->zbinVariableName), NULL, &lvl);

		printf("dimensions: getting %s\n", st->secIDVariableName);
		ncdiminq(st->hd.fid, ncdimid(st->hd.fid, st->secIDVariableName), NULL, &fce);

		/* Check against counters */
		totndata += nt;
		if (fce > maxfaces)
			maxfaces = fce;
		if (bx > maxboxes)
			maxboxes = bx;
		if (lvl < minlevel)
			minlevel = lvl;
		if (lvl > maxlevel)
			maxlevel = lvl;

		/* Close hydro file */
		close_hydro(st);

	}

	/* Check aganist stored values */
	st->ndata = totndata;
	st->nhdsteps = st->ndata;

	printf("ndata: %d, nface: %d vs maxface: %d, nbox: %d vs maxbox: %d\n", st->ndata, st->nface, maxfaces, st->nbox, maxboxes);

	if (st->nface < maxfaces) {
		warn("One of the raw hydro files contains more faces than the bgm file - reset nface\n");
		st->nface = maxfaces;
	}
	if (st->nbox < maxboxes) {
		warn("One of the raw hydro files contains more boxes than the bgm file - reset nbox\n");
		st->nbox = maxboxes;
	}
	if (st->wcnz > minlevel) {
		warn("One of the raw hydro files contains less water column layers (%d) than the parameter file specifies (%d)\n", minlevel, st->wcnz);
	}
	if (st->wcnz < maxlevel) {
		warn("One of the raw hydro files contains less water column layers (%d) than the parameter file specifies (%d)\n", maxlevel, st->wcnz);
	}

	return;
}
