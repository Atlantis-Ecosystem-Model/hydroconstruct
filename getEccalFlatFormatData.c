/**********************************************************************

 getflatinputdata.c

 by:	Beth Fulton
 date:	2/3/2006

 comments: This file is the code that reads in flat file version of
 raw hydrodynamics flows, and temperature and salinity profiles

 assumptions: assumes raw input is a series of text files with entry values
 per depth per face in sverdups for flow, deg C for temperature
 and PSU for salinity

 revisions   11/12/2006 added the code to read in excel estimated flows
 (under col_face_type input_type). Also added the Get_Date()
 routine.

 18-03-2009 Bec Gorton
 Added code to read in data for cam - this is when the input_stype == col_data_type.
 Added code to read in the dept_layer information so that the horizontal
 flux values can be distributed in the water column.

 *************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <math.h>
#include <string.h>
#include <sjwlib.h>
#include "mainexchange.h"


/***************************************************************************************************
 Routine to get the exchange data out of the raw data file - Jeff Dunn style flat input (not netcdf files)
 ****************************************************************************************************/
void get_Eccal_flat_hydro(Stuff *st) {
	char key[STSLEN];
	FILE *fp;
	int i, year, month, day, hour, minute, second, ti, fce_id, dbox, layer;
	double current_time, time_of_year, calc_time,
		dp1, dp2, dp3, dp4;
	int entry_i = 0;
	int buflen = 500;
	char buf[500];
	double lastTimeUsed = 0.0;
	double vertFlux, horizontalFlux, temp, salt;

	st->minExyear = MAXINT;
	ti = 0;

	/* Initialise face mapping - assumed to use bgm set as data not coming from netcdf file format */
	for (i = 0; i < st->nface; i++)
		facemap[i] = i;

	/* Loop through the exchange files */
	for (i = 0; i < st->hd.nfiles; i++) {
		if (verbose)
			printf("Reading raw exchange file %s\n", st->hd.fname[i]);

		/* Open the exchange file */
		if ((fp = fopen(st->hd.fname[i], "r")) == NULL)
			quit("hydro_exchange_dataread: Can't open '%s'\n", st->hd.fname[i]);
		/*
		 File format is:
		 t = time in yrday (1 = 01-jan-1958)
		 face = face index # (as in your figure)
		 polygon = polygon index #
		 layer = layer index # (layer 1 is the topmost layer)
		 trans = transport across the face in m^3/s (positive values are out of
		 the box)
		 w = transport across the bottom of the box in m^3/s (positive values are
		 upward velocity)
		 tempmean = temperature of the box in deg C
		 saltmean = salinity of the box in psu

		 0.7000000E+01 1. 1. 1. -0.1000000E+35 -0.1000000E+35 -0.1000000E+35 -0.1000000E+35

		 */

		/* Read in each line */
		while (fgets(buf, buflen, fp) != NULL) {

			printf("buf = %s\n", buf);
			if (sscanf(buf, "%lf %lf %lf %lf %lf %lf %lf %lf", &dp1, &dp2, &dp3, &dp4, &horizontalFlux, &vertFlux, &temp, &salt) != 8)
				quit("not enough entries for '%s' - with %lf %lf %lf %lf %lf %lf %lf %lf\n",  buf, dp1, dp2, dp3, dp4, horizontalFlux, vertFlux, temp, salt);

//			if (sscanf(buf, "%lf %lf %lf %lf %lf %lf %lf %lf", &time, &faceid, &boxid, &binid, &vertTransport, &horizonTransport,&tempmean, &saltmean ) != 8)
//				quit("not enough entries for '%s' - with %lf %lf %lf %lf %lf %lf %lf %lf\n", buf, time, faceid, boxid, binid. vertTransport, horizonTransport, tempmean, saltmean);
//
//
			current_time = dp1;
			/* Need to work out what timestep this is */
			if(current_time > lastTimeUsed){
				lastTimeUsed = current_time;
				ti++;
			}
			exchangedate[ti][time_id] = current_time;

			/* Calculate the time of year value for possible padding */

			// Get the year.
			//void    todat(double j, int *y, int *mo, int *d, int *h, int *mi, int *s);
			todat(dp1, &year, &month, &day, &hour, &minute, &second);
			/* Get TofY in days */
			sprintf(key, "%d-01-01 00:00:00", year);
			calc_time = GetTime(key);

			time_of_year = floor(((current_time - (calc_time - reference_time)) / st->dt) + 0.5);
			exchangedate[ti][TofY_id] = time_of_year;

			/* Save time and date for temp-salt data too */
			tempsaltdate[ti][TofY_id] = exchangedate[ti][TofY_id];
			tempsaltdate[ti][time_id] = exchangedate[ti][time_id];

			/* Store the raw exhanges, temp and salt values */
			dbox = (int) (dp3); // BoxID from column 1
			fce_id = (int) (dp2); // FaceID from column 2 and then use lookup table to do the translation
			layer = (int) (dp4);


			//rawexchanges = (double ****)d_alloc4d(2,st->nface,st->wcnz,st->nhdsteps); // two entries per face as can flow two ways
			if(horizontalFlux <= (st->hd.missing_data + 1)) {
				/* DO nothing as data missing */
				if (!st->sdiff)
					rawexchanges[ti][layer][fce_id][entry_i] += 0;
				else
					rawexchanges[ti][layer][fce_id][entry_i] += 0.00000001;
			} else{
				rawexchanges[ti][layer][fce_id][entry_i] += horizontalFlux;
			}

			if(temp < (st->hd.missing_data + 1)) {
				/* DO nothing as data missing */
				if(verbose > 1)
					printf("missing temp data\n");
			} else
				rawtempsalt[ti][layer][dbox][temp_id] = temp;

			if(salt < (st->hd.missing_data + 1)) {
				/* DO nothing as data missing */
				if(verbose > 1)
					printf("missing salt data\n");
			} else
				rawtempsalt[ti][layer][dbox][salt_id] = salt;

			if(verbose > 1){
				printf("Storing rawtemp[%d][%d][%d] = %e\n", ti, dbox, layer, rawtempsalt[ti][layer][dbox][temp_id]);
				printf("Storing rawsalt[%d][%d][%d] = %e\n", ti, dbox, layer, rawtempsalt[ti][layer][dbox][salt_id]);
			}

		}

		/* Close the exchange file */
		fclose(fp);
	}
}


