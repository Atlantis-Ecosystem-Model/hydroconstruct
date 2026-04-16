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

void Load_Layer_Table(Stuff *st);
void getDatesFromFileName(char *fileName, int *year, int *month, int *day);
void getDatesFromPathName(char *fileName, int *year, int *month, int *day);
int getFaceIndex(Stuff *st, int baseboxID, int adjBox);
int getMonthIndex(char *monthStr);
void trim(char *s);
char *months[12] = { "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec" };

/* This file was originally for Cam in Mexico and assumes
	 the hydro file has the format:

	 Polygon             Face               Time step      Water flux rate

	 0.                  1.                   1                 -732463.
	 0.                  1.                   2                 -536770.

	 */

static int read_Hydro_col_data(Stuff *st, FILE *fp, int *data_dt, int *old_nline, int year, int month, int day){

	char key[STSLEN];
	int fce, dp, nline, start_day, start_month,
			start_year, ti, fce_id, dbox, newdatid;
	double current_time, time_of_year, calc_time, dp1, dp2, dp3, dp4;
	int entry_i = 0;
	int buflen = 500;
	char buf[500];
	double prop, totalFlux, sum;

	*data_dt = (int)(st->datadt);
	/*if(st->wcnz > 1)
	 quit("The input_type 5 assumes only one watercolumn layer for now - recode if more are needed\n");
	 */
	/* Read in data */
	nline = 0;
	//printf("got ot here\n");
	/* Read in the first line and ignore */

	if(fgets(buf, buflen, fp) == NULL){
		quit("error in read_Hydro_col_data(). Problem loaded data form file\n");
	}

	/* Read in each line */
	while (fgets(buf, buflen, fp) != NULL) {

		if (strstr(buf, "Polygon") != NULL || strstr(buf, "No.") != NULL || strstr(buf, "transect") != NULL || strlen(buf) == 0 || strcmp(
				buf, "\n") == 0)
			continue;


		//printf("buf = %s\n", buf);
		if (sscanf(buf, "%lf %lf %lf %lf", &dp1, &dp2, &dp3, &dp4) != 4)
			quit("not enough entries for '%s' - with %lf %lf %lf %lf\n", buf, dp1, dp2, dp3, dp4);

		newdatid = (int) (dp3);
		ti = *old_nline + newdatid - 1; // New time step of entry based on base date + value from column 3 (minus 1 as values start at 1 not 0



		if (exchangedate[ti][done_id] < 1) {
			/* Get current day, month and year - Assumes one timestep (12 hours) worth of entry per line */
			start_day = day;
			start_month = month;
			start_year = year;
			Get_Date(start_day, start_month, start_year, *data_dt, &day, &month, &year);

			/* Get date into ISO date/time format (YYYY-MM-DD HH:MM:SS).*/
			if (month < 10) {
				if (day < 10) {
					sprintf(key, "%d-0%d-0%d 00:00:00", year, month, day);
				} else {
					sprintf(key, "%d-0%d-%d 00:00:00", year, month, day);
				}
			} else {
				if (day < 10) {
					sprintf(key, "%d-%d-0%d 00:00:00", year, month, day);
				} else {
					sprintf(key, "%d-%d-%d 00:00:00", year, month, day);
				}
			}


			/* Julian seconds */
			calc_time = GetTime(key);
			current_time = calc_time - reference_time;
			exchangedate[ti][time_id] = current_time;

			/* Get TofY in days */
			sprintf(key, "%d-01-01 00:00:00", year);
			calc_time = GetTime(key);
			/* Use st->dt not 86400 so enough entries to fit tidal cycle within day, i.e.
			 expand list to fit (i.e. have 2 entries per day etc) */
			time_of_year = floor(((current_time - (calc_time - reference_time)) / st->dt) + 0.5);
			exchangedate[ti][TofY_id] = time_of_year;

			if (verbose) {
				printf("current_time: %e, calc_time: %e, time_of_year: %e (ref_time: %e) - ",  current_time, calc_time,
						time_of_year, reference_time);
				printf("which matches %e (TofY = %e)\n", exchangedate[ti][time_id], exchangedate[ti][TofY_id]);
			}
			exchangedate[ti][done_id] = 1;
			nline++;
		}

		dbox = (int) (dp1); // BoxID from column 1
		fce = (int) (dp2); // FaceID from column 2 and then use lookup table to do the translation
		fce_id = col_faceid[dbox][fce][0];
		prop = col_faceProportion[dbox][fce][0];
		sum= 0;

		if (fce_id > -1){// Don't store flows from boundary box out of the system
			if(verbose > 3){
				printf("st->boxnz[%d] = %d, st->boxtotdepth[dbox] = %e, total flux = %e \n",
						dbox,
						st->boxnz[dbox],
						st->boxtotdepth[dbox],
						dp4 * prop);
			}
			if(st->boxnz[dbox] > 0){
				totalFlux = dp4 * prop;

				/*
				 * We need to allocate the horizontal flux across this face through the water column
				 * The flux at a layer = totalFlux * layer_depth/total_depth;
				 */
				for(dp = 0; dp < st->boxnz[dbox]; dp++){

					rawexchanges[ti][dp][fce_id][entry_i] = totalFlux * depth_layer[dbox][dp]/ st->boxtotdepth[dbox];
					sum += totalFlux * depth_layer[dbox][dp]/ st->boxtotdepth[dbox];
					if(verbose > 3){
						printf("st->boxdz[%d][%d] = %e, rawexchanges[%d][%d][%d][%d] = %.20e\n",
								dbox, dp, depth_layer[dbox][dp],
								ti, dp, fce_id, entry_i, rawexchanges[ti][dp][fce_id][entry_i]);
					}
				}
				if(fabs(sum - totalFlux) > tiny){
					quit("get_flat_hydro: ti = %d, box = %d diff = %e allocation of horizontal flux to layers in water column did not meet the summation check. Change verbose to > 3 to see intermediate output\n", ti, dbox, (sum - totalFlux));
				}

			}
			else{
				/* don't do anything if this box has a zero depth */
				rawexchanges[ti][0][fce_id][entry_i] = 0.0;
				if(verbose > 3)
					printf("rawexchanges[%d][0][%d][%d]= %e\n", ti, fce_id, entry_i, rawexchanges[ti][0][fce_id][entry_i]);
			}
		} // Data from column 4

		/* Now do the second face if its not -1 */
		if (col_faceid[dbox][fce][1] > 0) {
			fce_id = col_faceid[dbox][fce][1];
			prop = col_faceProportion[dbox][fce][1];
			sum = 0.0;
			if (fce_id > -1){ // Don't store flows from boundary box out of the system
				if(verbose > 3){
					printf("st->boxnz[%d] = %d, st->boxtotdepth[dbox] = %e, total flux = %e \n",
							dbox,
							st->boxnz[dbox],
							st->boxtotdepth[dbox],
							dp4 * prop);
				}
				if(st->boxnz[dbox] > 0){
					totalFlux = dp4 * prop;
					/*
					 * We need to allocate the horizontal flux across this face through the water column
					 * The flux at a layer = totalFlux * layer_depth/total_depth;
					 */

					for(dp = 0; dp < st->boxnz[dbox]; dp++){

						rawexchanges[ti][dp][fce_id][entry_i] = totalFlux * depth_layer[dbox][dp]/ st->boxtotdepth[dbox];
						sum += rawexchanges[ti][dp][fce_id][entry_i];
						if(verbose > 3){
							printf("st->boxdz[%d][%d] = %e, rawexchanges[ti][dp][fce_id][entry_i] = %e\n",
									dbox, dp, depth_layer[dbox][dp],
									rawexchanges[ti][dp][fce_id][entry_i]);
						}
					}
					if(fabs(sum - totalFlux) > tiny)
						quit("get_flat_hydro: ti = %d, box = %d diff = %e allocation of horizontal flux to layers in water column did not meet the summation check. Change verbose to > 3 to see intermediate output\n", ti, dbox, (sum - totalFlux));
				}else{
					/* don't do anything if this box has a zero depth */
					rawexchanges[ti][0][fce_id][entry_i] = 0.0;
				}
			} // Data from column 4
		}
	}

	*old_nline += nline;

	return 1;
}

/* This file was originally for Cam in Mexico and assumes
	 the hydro file has the format:

	        Polygon           Face      Time Step          Depth           Flux
         	 number         number         (12)hr          Layer         [m3/s]
             0.              1.              1.              1.       5633.179
             0.              1.              1.              2.       9210.793
             0.              1.              1.              3.      40082.371
             0.              1.              1.              4.     120441.215
             0.              1.              1.              5.       9175.180
             0.              1.              1.              6.          0.000
	 */



static int read_Hydro_New_col_data(Stuff *st, FILE *fp, int *data_dt, int old_nline, int year, int month, int day){

	char key[STSLEN];
	int fce, dp, nline, ij = 0, start_day, start_month,
			start_year, ti, fce_id, dbox, newdatid;
	double current_time, time_of_year, calc_time, dp1, dp2, dp3, dp4, dp5;
	int entry_i = 0;
	int buflen = 500;
	char buf[500];

	*data_dt = (int)(st->datadt);
	/*if(st->wcnz > 1)
	 quit("The input_type 5 assumes only one watercolumn layer for now - recode if more are needed\n");
	 */
	/* Read in data */
	nline = 0;
	//printf("got ot here\n");
	/* Read in the first line and ignore */
	if(fgets(buf, buflen, fp) == NULL){
		quit("error in read_Hydro_New_col_data(). Problem loaded data form file\n");
	}

	/* Read in each line */
	while (fgets(buf, buflen, fp) != NULL) {

		if (strstr(buf, "Polygon") != NULL || strstr(buf, "No.") != NULL || strstr(buf, "transect") != NULL || strstr(buf, "number") != NULL|| strlen(buf) == 0 || strcmp(
				buf, "\n") == 0 || strstr(buf, "Box_ID") != NULL)
			continue;

		if (sscanf(buf, "%lf %lf %lf %lf %lf", &dp1, &dp2, &dp3, &dp4, &dp5) != 5)
			quit("not enough entries for '%s' - with %lf %lf %lf %lf %lf\n", buf, dp1, dp2, dp3, dp4, dp5);

		newdatid = (int) (dp3);
		ti = old_nline + newdatid - 1; // New time step of entry based on base date + value from column 3 (minus 1 as values start at 1 not 0)

		if (exchangedate[ti][done_id] < 1) {
			/* Get current day, month and year - Assumes one timestep (12 hours) worth of entry per line */
			start_day = day;
			start_month = month;
			start_year = year;

			Get_Date(start_day, start_month, start_year, *data_dt, &day, &month, &year);

			/* Get date into ISO date/time format (YYYY-MM-DD HH:MM:SS).*/
			if (month < 10) {
				if (day < 10) {
					sprintf(key, "%d-0%d-0%d 00:00:00", year, month, day);
				} else {
					sprintf(key, "%d-0%d-%d 00:00:00", year, month, day);
				}
			} else {
				if (day < 10) {
					sprintf(key, "%d-%d-0%d 00:00:00", year, month, day);
				} else {
					sprintf(key, "%d-%d-%d 00:00:00", year, month, day);
				}
			}


			/* Julian seconds */
			calc_time = GetTime(key);
			current_time = calc_time - reference_time;
			exchangedate[ti][time_id] = current_time;

			/* Get TofY in days */
			sprintf(key, "%d-01-01 00:00:00", year);
			calc_time = GetTime(key);
			/* Use st->dt not 86400 so enough entries to fit tidal cycle within day, i.e.
			 expand list to fit (i.e. have 2 entries per day etc) */
			time_of_year = floor(((current_time - (calc_time - reference_time)) / st->dt) + 0.5);
			exchangedate[ti][TofY_id] = time_of_year;

			if (verbose) {
				printf("i: %d, current_time: %e, calc_time: %e, time_of_year: %e (ref_time: %e) - ", ij, current_time, calc_time,
						time_of_year, reference_time);
				printf("which matches %e (TofY = %e)\n", exchangedate[ti][time_id], exchangedate[ti][TofY_id]);
			}
			exchangedate[ti][done_id] = 1;
			nline++;
		}

		dbox = (int) (dp1); // BoxID from column 1
		fce = (int) (dp2) - 1; // FaceID from column 2 and then use lookup table to do the translation
		dp = (int)dp4 - 1;

		fce_id = col_faceid[dbox][fce][0];

		if(fce_id == 54 && ti == 0){
			printf("st->boxnz[%d] = %d, st->boxtotdepth[dbox] = %e, total flux = %e, entry_i = %d \n",
									dbox,
									st->boxnz[dbox],
									st->boxtotdepth[dbox],
									dp5, entry_i);
		}

		if (fce_id > -1){// Don't store flows from boundary box out of the system
			if(verbose > 3){
				printf("st->boxnz[%d] = %d, st->boxtotdepth[dbox] = %e, total flux = %e \n",
						dbox,
						st->boxnz[dbox],
						st->boxtotdepth[dbox],
						dp5			);
			}

			rawexchanges[ti][dp][fce_id][entry_i] += dp5;
			if(verbose > 3)
				printf("rawexchanges[%d][0][%d][%d]= %e\n", ti, fce_id, entry_i, rawexchanges[ti][0][fce_id][entry_i]);

			if(fce_id == 54 && ti == 0){
				printf("rawexchanges[%d][%d][%d][%d] = %e\n", ti, dp, fce_id, entry_i, rawexchanges[ti][dp][fce_id][entry_i]);
			}
		} // Data from column 4

	}
	old_nline += nline;

	return 1;
}


/* This assumes file format is

 date year month day
 dt t
 nline y
 nface_id z
 face_id i j k...... z
 line1 val11 val21 val31 val41...... valn1
 |
 linei val1i val2i val3i val4i...... valni
 |
 |
 liney val1z val2z val3z val4z...... valn

 where val1 is the entry for face1, val2 is the entry for face 2 etc.

 */
static int read_Hydro_col_face(Stuff *st, FILE *fp, int *data_dt, int *old_nline, int year, int month, int day){

	char key[STSLEN];
	int fce, dp, nline, ij, start_day, start_month,
			start_year, ti, fce_id, nface_id;
	double current_time, time_of_year, calc_time, val1;
	int entry_i = 0;


	if (st->wcnz > 1)
		quit("The input_type 3 assumes only one watercolumn layer for now - recode if more are needed\n");

	/* rewind file and find time step between data entries (rather than assume daily step) */
	fseek(fp, 0L, SEEK_SET);

	sprintf(key, "dt");
	skipToKeyEnd(fp, key);
	if (fscanf(fp, "%d", data_dt) != 1)
		quit(
				"the timestep between data (i.e. if daily or weekly etc) must be marked in days at the head of the file as 'dt x' where x is the number of days between data points (if daily = 1, if weekly = 7 etc)\n");
	if (*data_dt > 28)
		quit(
				"to date code assumes this data will be entered at less than monthly intervals if this is not the case (e.g. here have dt = %d) Get_Date() needs to be recoded\n",
				*data_dt);

	/* rewind file and find number of lines of data to read in */
	fseek(fp, 0L, SEEK_SET);

	sprintf(key, "nline");
	skipToKeyEnd(fp, key);
	if (fscanf(fp, "%d", &nline) != 1)
		quit("the number of lines must be marked at the head of the file as 'nline x' where x is the number of lines of data to be read in\n");

	/* rewind file and find number of columns of data to read in */
	fseek(fp, 0L, SEEK_SET);

	sprintf(key, "nface_id");
	skipToKeyEnd(fp, key);
	if (fscanf(fp, "%d", &nface_id) != 1)
		quit(
				"the number of faces to readin must be marked at the head of the file as 'nface_id x' where x is the number of columns of data to be read in\n");

	/* Initialise fceid array */
	faceid = (int *) i_alloc1d(nface_id);

	/* rewind file */
	fseek(fp, 0L, SEEK_SET);

	/* Load face ids */
	sprintf(key, "face_id");
	skipToKeyEnd(fp, key);
	for (fce = 0; fce < nface_id; fce++) {
		faceid[fce] = -1;
		if (fscanf(fp, "%d", &dp) != 1) {
			if (!fce)
				quit("Failed to fine any face id data\n");
			else
				quit("Only found values for %d columns - please fill out the rest of the face ids (nface_id = %d)\n", fce, nface_id);
		} else {
			/* Found face_id so store it */
			faceid[fce] = dp;
		}
	}

	for (ij = 1; ij <= nline; ij++) {
		ti = *old_nline + ij;

		/* Get current day, month and year - Assumes one days worth of entry per line */
		start_day = day;
		start_month = month;
		start_year = year;
		Get_Date(start_day, start_month, start_year, *data_dt, &day, &month, &year);

		/* Get date into ISO date/time format (YYYY-MM-DD HH:MM:SS).*/
		if (month < 10) {
			if (day < 10) {
				sprintf(key, "%d-0%d-0%d 00:00:00", year, month, day);
			} else {
				sprintf(key, "%d-0%d-%d 00:00:00", year, month, day);
			}
		} else {
			if (day < 10) {
				sprintf(key, "%d-%d-0%d 00:00:00", year, month, day);
			} else {
				sprintf(key, "%d-%d-%d 00:00:00", year, month, day);
			}
		}

		/* Julian seconds */
		calc_time = GetTime(key);
		current_time = calc_time - reference_time;
		exchangedate[ti][time_id] = current_time;

		/* Get TofY in days */
		sprintf(key, "%d-01-01 00:00:00", year);
		calc_time = GetTime(key);
		/* Use st->dt not 86400 so enough entries to fit tidal cycle within day, i.e.
		 expand list to fit (i.e. have 2 entries per day etc) */
		time_of_year = floor(((current_time - (calc_time - reference_time)) / st->dt) + 0.5);
		exchangedate[ti][TofY_id] = time_of_year;

		if (verbose) {
			printf("i: %d, current_time: %e, calc_time: %e, time_of_year: %e (ref_time: %e) - ", ij, current_time, calc_time, time_of_year,
					reference_time);
			printf("which matches %e (TofY = %e)\n", exchangedate[ti][time_id], exchangedate[ti][TofY_id]);
		}

		/* Now read face flow data */
		sprintf(key, "line%d", ij);
		skipToKeyEnd(fp, key);
		dp = 0;
		for (fce = 0; fce < nface_id; fce++) {
			if (fscanf(fp, "%lf", &val1) != 1) {
				if (!fce)
					quit("Failed to fine any data on line %d\n", ij);
				else
					quit("Only found values for faces up to %d - please fill out the rest of the data columns (nface_id = %d)\n", fce,
							nface_id);
			} else {
				/* Found data so store it */
				if (st->sdiff) {
					if (!val1)
						val1 = 0.00000001;
				}
				fce_id = faceid[fce];
				rawexchanges[ti][dp][fce_id][entry_i] = val1;
			}
		}
	}
	*old_nline += nline;

	return 1;

}

/* For flat CSIRO entry assumes one data set (and date) per file, the
	 alternative format below assumes one per line so must split the code here
	 */
static int read_Hydro_csiro(Stuff *st, FILE *fp, int *data_dt, int *old_nline, char *key, int year, int month, int day, int i){


	int fce, faced, dp;
	double current_time, time_of_year, calc_time, dp1, dp2, dp3, dp4, val1, val2, val3, val4;
	int entry_i = 0;

	*old_nline = 0;
	*data_dt = 1;

	/* Get date into ISO date/time format (YYYY-MM-DD HH:MM:SS).*/
	if (month < 10) {
		if (day < 10) {
			sprintf(key, "%d-0%d-0%d 00:00:00", year, month, day);
		} else {
			sprintf(key, "%d-0%d-%d 00:00:00", year, month, day);
		}
	} else {
		if (day < 10) {
			sprintf(key, "%d-%d-0%d 00:00:00", year, month, day);
		} else {
			sprintf(key, "%d-%d-%d 00:00:00", year, month, day);
		}
	}

	/* Julian seconds */
	calc_time = GetTime(key);
	current_time = calc_time - reference_time;
	exchangedate[i][time_id] = current_time;

	if (verbose)
		printf("Readin date %s ", key);

	/* Get TofY in days */
	sprintf(key, "%d-01-01 00:00:00", year);
	calc_time = GetTime(key);
	/* Use st->dt not 86400 so enough entries to fit tidal cycle within day, i.e.
	 expand list to fit (i.e. have 2 entries per day etc) */
	time_of_year = floor(((current_time - (calc_time - reference_time)) / st->dt) + 0.5);
	exchangedate[i][TofY_id] = time_of_year;

	if (verbose) {
		printf("i: %d, current_time: %e, calc_time: %e, time_of_year: %e (ref_time: %e) - ", i, current_time, calc_time, time_of_year,
				reference_time);
		printf("which matches %e (TofY = %e)\n", exchangedate[i][time_id], exchangedate[i][TofY_id]);
	}

	/* rewind file */
	fseek(fp, 0L, SEEK_SET);

	/* Read raw exchanges - assuming entries are one face per line */
	for (fce = 0; fce < st->nface; fce++) {
		faced = fce;
		dp1 = 0;
		dp2 = 0;
		dp3 = 0;
		dp4 = 0;

		sprintf(key, "face%d", faced);
		skipToKeyEnd(fp, key);

		if (st->n_inline > 1) {
			/* Read raw exchanges - assumes all entries in one file, with each column = a different time period */
			for (fce = 0; fce < st->nface; fce++) {
				for (dp = 0; dp < st->wcnz; dp++) {
					sprintf(key, "f%dd%d", fce, dp);
					skipToKeyEnd(fp, key);
					if (fscanf(fp, "%lf %lf %lf %lf", &val1, &val2, &val3, &val4) != 4)
						quit("not enough entries for this face (%d) and water layer (%d) combination\n", fce, dp);

					rawexchanges[0][dp][fce][entry_i] = val1;
					rawexchanges[1][dp][fce][entry_i] = val2;
					rawexchanges[2][dp][fce][entry_i] = val3;
					rawexchanges[3][dp][fce][entry_i] = val4;

				}
			}

			if (i > 0)
				quit("With mutiple columns of flow data per line, code currently assumes a single file. Recode if desire otherwise\n");
		} else {
			/* Read raw exchanges - if single data time period per line */
			switch (st->wcnz) {
			case 1:
				if (fscanf(fp, "%lf", &dp1) != 1)
					quit("not enough entries for %s - with %lf\n", key, dp1);

				if (st->sdiff) {
					if (!dp1)
						dp1 = 0.00000001;
				}

				if (st->verbose_exchange) {
					printf("read f%d: ", faced);
					printf("%lf\n", dp1);
				}

				dp = 0;
				if (dp1 < st->hd.missing_data) {
					if (!st->sdiff)
						rawexchanges[i][dp][fce][entry_i] = 0;
					else
						rawexchanges[i][dp][fce][entry_i] = 0.00000001;
				} else
					rawexchanges[i][dp][fce][entry_i] = dp1;
				break;
			case 4:
				if (fscanf(fp, "%lf %lf %lf %lf", &dp1, &dp2, &dp3, &dp4) != 4)
					quit("not enough entries for %s - with %lf %lf %lf %lf\n", key, dp1, dp2, dp3, dp4);

				if (st->sdiff) {
					if (!dp1)
						dp1 = 0.00000001;
					if (!dp2)
						dp2 = 0.00000001;
					if (!dp3)
						dp3 = 0.00000001;
					if (!dp4)
						dp4 = 0.00000001;
				}

				if (st->verbose_exchange) {
					printf("read f%d: ", faced);
					printf("%lf %lf %lf %lf\n", dp1, dp2, dp3, dp4);
				}

				dp = 0;
				if (dp1 < st->hd.missing_data) {
					if (!st->sdiff)
						rawexchanges[i][dp][fce][entry_i] = 0;
					else
						rawexchanges[i][dp][fce][entry_i] = 0.00000001;
				} else
					rawexchanges[i][dp][fce][entry_i] = dp1;

				dp = 1;
				if (dp2 < st->hd.missing_data) {
					if (!st->sdiff)
						rawexchanges[i][dp][fce][entry_i] = 0;
					else
						rawexchanges[i][dp][fce][entry_i] = 0.00000001;
				} else
					rawexchanges[i][dp][fce][entry_i] = dp2;

				dp = 2;
				if (dp3 < st->hd.missing_data) {
					if (!st->sdiff)
						rawexchanges[i][dp][fce][entry_i] = 0;
					else
						rawexchanges[i][dp][fce][entry_i] = 0.00000001;
				} else
					rawexchanges[i][dp][fce][entry_i] = dp3;

				dp = 3;
				if (dp4 < st->hd.missing_data) {
					if (!st->sdiff)
						rawexchanges[i][dp][fce][entry_i] = 0;
					else
						rawexchanges[i][dp][fce][entry_i] = 0.00000001;
				} else
					rawexchanges[i][dp][fce][entry_i] = dp4;
				break;
			default:
				quit("Currently only set up for 1 or 4 water column layers recode for anything else\n");
				break;
			}
		}
	}

	return 1;

}

/***************************************************************************************************
 Routine to get the exchange data out of the raw data file - Jeff Dunn style flat input (not netcdf files)
 ****************************************************************************************************/
void get_flat_hydro(Stuff *st) {
	char key[STSLEN];
	FILE *fp;
	int i, year, month, day, old_nline;
	int data_dt = 1;
	double current_time, calc_time;
	int timeIndex;

	st->minExyear = MAXINT;

	if (st->input_style == col_data_type){
		Load_Face_Lookup_Table(st);
		Load_Layer_Table(st);
	}

	if (st->input_style == new_col_data_type){
		Load_Face_Lookup_Table(st);
	}

	/* Initialise face mapping - assumed to use bgm set as data not coming from netcdf file format */
	for (i = 0; i < st->nface; i++)
		facemap[i] = i;

	/* Loop through the exchange files */
	old_nline = 0;
	for (i = 0; i < st->hd.nfiles; i++) {
		if (verbose)
			printf("Reading raw exchange file %s\n", st->hd.fname[i]);

		/* Open the exchange file */
		if ((fp = fopen(st->hd.fname[i], "r")) == NULL)
			quit("hydro_exchange_dataread: Can't open '%s'\n", st->hd.fname[i]);

		if (st->input_style == col_data_type ) {
			getDatesFromFileName(st->hd.fname[i], &year, &month, &day);
		}else if (st->input_style == new_col_data_type){
			getDatesFromPathName(st->hd.fname[i], &year, &month, &day);
		} else {

			/* rewind file */
			fseek(fp, 0L, SEEK_SET);

			/* Read time */
			sprintf(key, "date");
			skipToKeyEnd(fp, key);
			if (fscanf(fp, "%d %d %d", &year, &month, &day) != 3)
				quit("the date of the flow must be give as 'date year month day' (e.g. year 2003 5 31)\n");
		}
		if (year < st->minExyear) {
			st->minExyear = year;
			st->monthEx = month;
			st->dayEx = day;
			st->Exmax_id = i;
			if (st->minExyear < st->REFyear)
				st->timecorrect_on = 1;
			else
				st->timecorrect_on = 2;
		}

		if (year > st->maxExyear) {
			st->maxExyear = year;
			st->mmonthEx = month;
			st->mdayEx = day;
			st->Exmax_id = i;

			if (st->verbose_pad)
				printf("TSMax_id set to %d - yr: %d, mth: %d, day: %d\n ", st->Exmax_id, year, month, day);

		} else if (year == st->maxExyear) {
			if (month > st->mmonthEx) {
				st->maxExyear = year;
				st->mmonthEx = month;
				st->mdayEx = day;
				st->Exmax_id = i;

				if (st->verbose_pad)
					printf("TSMax_id set to %d - yr: %d, mth: %d, day: %d\n ", st->Exmax_id, year, month, day);
			} else if ((month == st->mmonthEx) && (day > st->mdayEx)) {
				st->maxExyear = year;
				st->mmonthEx = month;
				st->mdayEx = day;
				st->Exmax_id = i;

				if (st->verbose_pad)
					printf("TSMax_id set to %d - yr: %d, mth: %d, day: %d\n ", st->Exmax_id, year, month, day);
			}
		}

		switch (st->input_style) {

		case col_face_type:

			/* This assumes file format is

			 date year month day
			 dt t
			 nline y
			 nface_id z
			 face_id i j k...... z
			 line1 val11 val21 val31 val41...... valn1
			 |
			 linei val1i val2i val3i val4i...... valni
			 |
			 |
			 liney val1z val2z val3z val4z...... valn

			 where val1 is the entry for face1, val2 is the entry for face 2 etc.

			 */
			read_Hydro_col_face(st, fp, &data_dt, &old_nline, year, month, day);
			break;

		case col_data_type:
			/* This file was originally for Cam in Mexico and assumes
			 the hydro file has the format:

			 Polygon             Face               Time step      Water flux rate

			 0.                  1.                   1                 -732463.
			 0.                  1.                   2                 -536770.

			 */

			read_Hydro_col_data(st, fp, &data_dt, &old_nline, year, month, day);

			break;

		case new_col_data_type:


			/*
			 *
				Polygon           Face      Time Step          Depth           Flux
         	 	 number         number         (12)hr          Layer         [m3/s]
             	 	 0.              1.              1.              1.       5633.179
             	 	 0.              1.              1.              2.       9210.793
             	 	 0.              1.              1.              3.      40082.371
             	 	 0.              1.              1.              4.     120441.215
             	 	 0.              1.              1.              5.       9175.180

             */

			/* Get date into ISO date/time format (YYYY-MM-DD HH:MM:SS).*/
			if (month < 10) {
				if (day < 10) {
					sprintf(key, "%d-0%d-0%d 00:00:00", year, month, day);
				} else {
					sprintf(key, "%d-0%d-%d 00:00:00", year, month, day);
				}
			} else {
				if (day < 10) {
					sprintf(key, "%d-%d-0%d 00:00:00", year, month, day);
				} else {
					sprintf(key, "%d-%d-%d 00:00:00", year, month, day);
				}
			}


			/* Julian seconds */
			calc_time = GetTime(key);
			current_time = calc_time - reference_time;


			timeIndex = (int)(current_time/st->dt);

			read_Hydro_New_col_data(st, fp, &data_dt, timeIndex, year, month, day);
			break;
		default:

			/* For flat CSIRO entry assumes one data set (and date) per file, the
			 alternative format below assumes one per line so must split the code here
			 */

			read_Hydro_csiro(st, fp, &data_dt, &old_nline, key, year, month, day, i);
			break;
		}
		/* Close the exchange file */
		fclose(fp);
	}

	if (verbose)
		printf("Finished reading raw exchange files\n");
	return;
}

/*******************************************************************************
 Routine to determine month of the year for input data
 */
void Get_Date(int startday, int startmonth, int startyear, int data_dt, int *day, int *month, int *year) {
	int newday;
	int newmonth = startmonth;
	int newyear = startyear;

	newday = startday + data_dt;

	/* Check if startyear is a leap year */
	if ((startyear % 4) == 0) {
		/** startyear = Leap year **/
		/* Determine newday, newmonth and newyear */
		switch (startmonth) {
		case 1:
		case 3:
		case 5:
		case 7:
		case 8:
		case 10:
			/* Month end check on 31 days in the month */
			if (newday > 31) {
				newday -= 31;
				newmonth = startmonth + 1;
			}
			break;
		case 2:
			/* Month end check on 29 days in the month */
			if (newday > 29) {
				newday -= 29;
				newmonth = startmonth + 1;
			}
			break;
		case 4:
		case 6:
		case 9:
		case 11:
			/* Month end check on 30 days in the month */
			if (newday > 30) {
				newday -= 30;
				newmonth = startmonth + 1;
			}
			break;
		case 12:
			/* Year end check */
			if (newday > 31) {
				newday -= 31;
				newmonth = 1;
				newyear = startyear + 1;
			}
			break;
		default:
			quit("If got to a month == 0 or month > 12 then something is a miss! mpnth = %d\n", month);
			break;
		}
	} else {
		/** startyear != leap year */
		/* Determine newday, newmonth and newyear */
		switch (startmonth) {
		case 1:
		case 3:
		case 5:
		case 7:
		case 8:
		case 10:
			/* Month end check on 31 days in the month */
			if (newday > 31) {
				newday -= 31;
				newmonth = startmonth + 1;
			}
			break;
		case 2:
			/* Month end check on 28 days in the month */
			if (newday > 28) {
				newday -= 29;
				newmonth = startmonth + 1;
			}
			break;
		case 4:
		case 6:
		case 9:
		case 11:
			/* Month end check on 30 days in the month */
			if (newday > 30) {
				newday -= 30;
				newmonth = startmonth + 1;
			}
			break;
		case 12:
			/* Year end check */
			if (newday > 31) {
				newday -= 31;
				newmonth = 1;
				newyear = startyear + 1;
			}
			break;
		default:
			quit("If got to a month == 0 or month > 12 then something is a miss! month = %d\n", month);
			break;
		}
	}

	*day = newday;
	*month = newmonth;
	*year = newyear;

	return;
}


/***************************************************************************************************
 Routine to get the temperature and salinity data out of the raw data files - - Jeff Dunn style flat
 input (not netcdf files)
 ****************************************************************************************************/
void get_flat_temp_and_salt(Stuff *st) {
	char key[STSLEN];
	FILE *fp;
	int i, b, boxed, dp, year, month, day, dpA, dpB, dpC, dpD, old_nline, nline, newdatid, ti, k;
	double dpt1, dpt2, dpt3, dpt4, dps1, dps2, dps3, dps4, current_time, time_of_year, calc_time,
			dp1, dp2, dp3, dp4, ph;
	int entry_i = 0;
	int buflen = 500;
	char buf[500];
	int maxTime = -1;
	double depth;
	int lineIndex, timeIndex;

	dpD = dpC = dpA = dpB = 0;

	entry_i = 0;
	old_nline = 0;

	if (verbose)
		printf("Read in raw temperature and salinity profiles \n");

	printf("st->ts.ntsfiles = %d\n",  st->ts.ntsfiles);

	/* Loop through the exchange files */
	for (i = 0; i < st->ts.ntsfiles; i++) {
	//	if (verbose)
			printf("Reading raw temperature_salinity file %s\n", st->ts.tsfname[i]);

		/* Open the exchange file */
		if ((fp = fopen(st->ts.tsfname[i], "r")) == NULL)
			quit("hydro_exchange_dataread: Can't open %s\n", st->ts.tsfname[i]);

		if (st->input_style == col_data_type){
			getDatesFromFileName(st->ts.tsfname[i], &year, &month, &day);
		}else if (st->input_style == new_col_data_type){
			getDatesFromPathName(st->ts.tsfname[i], &year, &month, &day);

		} else {
			/* rewind file */
			fseek(fp, 0L, SEEK_SET);

			/* Read time */
			sprintf(key, "date");
			skipToKeyEnd(fp, key);
			if (fscanf(fp, "%d %d %d", &year, &month, &day) != 3)
				quit("the date of the temperature or salinity profile must be give as 'date year month day' (e.g. year 2003 5 31)\n");
		}

		if (year < st->minTSyear) {
			st->minTSyear = year;
			st->monthTS = month;
			st->dayTS = day;
			st->TSmin_id = i;
			if (st->minTSyear < st->REFyear)
				st->timecorrect_on = 1;
			else
				st->timecorrect_on = 2;
		}

		if ((year > st->maxTSyear)) {
			st->maxTSyear = year;
			st->mmonthTS = month;
			st->mdayTS = day;
			st->TSmax_id = i;

			if (st->verbose_pad)
				printf("TSMax_id set to %d - yr: %d, mth: %d, day: %d\n ", st->TSmax_id, year, month, day);

		} else if (year == st->maxTSyear) {
			if (month > st->mmonthTS) {
				st->maxTSyear = year;
				st->mmonthTS = month;
				st->mdayTS = day;
				st->TSmax_id = i;

				if (st->verbose_pad)
					printf("TSMax_id set to %d - yr: %d, mth: %d, day: %d\n ", st->TSmax_id, year, month, day);
			} else if ((month == st->mmonthTS) && (day > st->mdayTS)) {
				st->maxTSyear = year;
				st->mmonthTS = month;
				st->mdayTS = day;
				st->TSmax_id = i;

				if (st->verbose_pad)
					printf("TSMax_id set to %d - yr: %d, mth: %d, day: %d\n ", st->TSmax_id, year, month, day);
			}
		}

		/* Get date into ISO date/time format (YYYY-MM-DD HH:MM:SS).*/
		if (month < 10) {
			if (day < 10) {
				sprintf(key, "%d-0%d-0%d 00:00:00", year, month, day);
			} else {
				sprintf(key, "%d-0%d-%d 00:00:00", year, month, day);
			}
		} else {
			if (day < 10) {
				sprintf(key, "%d-%d-0%d 00:00:00", year, month, day);
			} else {
				sprintf(key, "%d-%d-%d 00:00:00", year, month, day);
			}
		}

		/* Julian seconds */
		calc_time = GetTime(key);
		current_time = calc_time - reference_time;
		tempsaltdate[i][time_id] = (int) (current_time);

		timeIndex = (int)(current_time/st->dt);

		printf("tempsaltdate[%d][time_id] = %e\n", i, tempsaltdate[i][time_id]);

		if (verbose)
			printf("Readin date %s ", key);

		/* Get TofY in days */
		sprintf(key, "%d-01-01 00:00:00", year);
		calc_time = GetTime(key);
		/* Use st->dt not 86400 so enough entries to fit tidal cycle within day, i.e.
		 expand list to fit (i.e. have 2 entries per day etc) */
		time_of_year = floor(((current_time - (calc_time - reference_time)) / st->dt) + 0.5);
		tempsaltdate[i][TofY_id] = (int) (time_of_year);

		if (verbose) {
			printf("i: %d, tscurrent_time: %e, tscalc_time: %e, tstime_of_year: %e (ref_time: %e) - ", i, current_time, calc_time, time_of_year,
					reference_time);
			printf("which matches %e (TofY = %e)\n", tempsaltdate[i][time_id], tempsaltdate[i][TofY_id]);
		}

		switch (st->input_style) {
		case col_data_type:	/* Intentional follow thur */
		case new_col_data_type :
		case guam_data_type:
			/* This file was originally for Cam in Mexico and assumes
			 the hydro file has the format

			 date year month day
			 dt t
			 nline y
			 LineNum Time_step  Polygon  Depth   Vertical_velocity[10^-5*m/s] Average_Temp[Celsius]  Average_Salinity[PartsPerThousand]
			 linei	   i		j			x			vali					valti					valsi

			 BEC TO FIX - CAN NOT TAG HEAD OF EACH FILE AS TOO MANY SO HOW GET TAG OUT OF FILE NAME? HOW KNOW FILES TO LOAD?
			 */

			ph = -1;
			nline = 0;
			maxTime = -1;
			/* Read in the first line and ignore */
			if(fgets(buf, buflen, fp) == NULL){
				quit("error in get_flat_temp_and_salt(). Problem loaded data form file %s\n", st->ts.tsfname[i]);
			}
			lineIndex = 0;
			/* Read in each line */
			while (fgets(buf, buflen, fp) != NULL) {

				trim(buf);
				//printf("buf = %s\n", buf);

				if (strstr(buf, "Polygon") != NULL || strstr(buf, "transect") != NULL || strlen(buf) == 0 || strcmp(buf, "\n") == 0 || strstr(buf,
						"Layer") != NULL || strstr(buf, "PartsPerThousand") != NULL || strstr(buf, "PartsPer1000") != NULL)

					continue;

				if(st->doPH){
					if (sscanf(buf, "%lf %lf %lf %lf %lf %lf %lf", &dp1, &dp2, &dp3, &dp4, &dpt1, &dps1, &ph) != 7){
						fprintf(stderr, "Line %d - %s in file %s is not the correct format\n", lineIndex, buf, st->ts.tsfname[i]);
						quit("\n\n get_flat_temp_and_salt: not enough entries for %s - with %lf %lf %lf %lf %lf %lf %lf in file %s\n", buf, dp1, dp2, dp3, dp4, dpt1,
							dps1, ph, st->ts.tsfname[i]);
					}

				}else{
					if (sscanf(buf, "%lf %lf %lf %lf %lf %lf", &dp1, &dp2, &dp3, &dp4, &dpt1, &dps1) != 6){
						fprintf(stderr, "Line %d - %s in file %s is not the correct format\n", lineIndex, buf, st->ts.tsfname[i]);
						quit("\n\n get_flat_temp_and_salt: not enough entries for %s - with %lf %lf %lf %lf %lf %lf in file %s\n", buf, dp1, dp2, dp3, dp4, dpt1,
								dps1, st->ts.tsfname[i]);
					}
				}

				lineIndex++;

				newdatid = (int) (dp1);
				ti = timeIndex + newdatid - 1; // New time step of entry based no base date + value from column 1 (minus 1 as values start at 1 not 0)

				if(ti  >= st->nvhdsteps){
					printf("line %s\n", buf);
					quit("The number of timesteps has been calculated as %d but in the file %s we have the %dth data point\n", st->nvhdsteps, st->ts.tsfname[i], ti);
				}
				b = (int) (dp2); // BoxID from column 2


				/* Find layer from Depth in column 3 */
				dp = st->boxnz[b];

				depth = 0;
				for (k = 0; k < st->wcnz; k++) {
					depth  = depth +  st->nominal_dz[k];
					//printf("st->boxdz[b][%d] = %f, depth= %e\n", k, st->nominal_dz[k], depth);
					if (dp3 <= depth){
						dpA = dp - 1 - k; // Tip layers upside down (so layer 0 next to sediment) in calcexchange
						dpA = k; // Tip layers upside down (so layer 0 next to sediment) in calcexchange
						break;
					}
				}

				if(dpA > (dp - 1)){
				//	abort();
					continue;
				}


				/* Enter raw temperature and salinity values */
				if (dpt1 < st->ts.temp_missing_data)
					rawtempsalt[ti][dpA][b][temp_id] = 0;
				else
					rawtempsalt[ti][dpA][b][temp_id] = dpt1;

				if (dps1 < st->ts.salt_missing_data)
					rawtempsalt[ti][dpA][b][salt_id] = 0;
				else
					rawtempsalt[ti][dpA][b][salt_id] = dps1;

				if(st->doPH){
					if (ph < st->ts.ph_missing_data)
						rawtempsalt[ti][dpA][b][ph_id] = 0;
					else
						rawtempsalt[ti][dpA][b][ph_id] = ph;
				}

				/* Enter raw vertical exchange values */
				if (dps1 < st->hd.missing_data){
					rawvertexchanges[ti][dpA][b][entry_i] = 0;
				}else{
					rawvertexchanges[ti][dpA][b][entry_i] = dp4 * 10e-5;
				}

				if (newdatid > maxTime) {
					nline++;
					maxTime = newdatid;
				}
			}
			old_nline += nline;
			break;
		default:
			/* Read raw exchanges */
			for (b = 0; b < st->nbox; b++) {
				boxed = b + 1;
				sprintf(key, "box%d", boxed);
				skipToKeyEnd(fp, key);
				if (fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf", &dpt1, &dpt2, &dpt3, &dpt4, &dps1, &dps2, &dps3, &dps4) != 8)
					quit("not enough entries for this box (%d) - there must be 4 depths\n", b);

				dp = st->boxnz[b];

				switch (dp) {
				case 0:
					dpA = 0;
					dpB = 1;
					dpC = 2;
					dpD = 3;
					break;
				case 1:
					dpA = 1;
					dpB = 0;
					dpC = 2;
					dpD = 3;
					break;
				case 2:
					dpA = 2;
					dpB = 1;
					dpC = 0;
					dpD = 3;
					break;
				case 3:
					dpA = 3;
					dpB = 2;
					dpC = 1;
					dpD = 0;
					break;
				}

				if (dpt1 < st->ts.temp_missing_data)
					rawtempsalt[i][dpA][b][temp_id] = 0;
				else
					rawtempsalt[i][dpA][b][temp_id] = dpt1;
				if (dps1 < st->ts.salt_missing_data)
					rawtempsalt[i][dpA][b][salt_id] = 0;
				else
					rawtempsalt[i][dpA][b][salt_id] = dps1;

				dp = 1;
				if (dpt2 < st->ts.temp_missing_data)
					rawtempsalt[i][dpB][b][temp_id] = 0;
				else
					rawtempsalt[i][dpB][b][temp_id] = dpt2;
				if (dps2 < st->ts.salt_missing_data)
					rawtempsalt[i][dpB][b][salt_id] = 0;
				else
					rawtempsalt[i][dpB][b][salt_id] = dps2;

				dp = 2;
				if (dpt3 < st->ts.temp_missing_data)
					rawtempsalt[i][dpC][b][temp_id] = 0;
				else
					rawtempsalt[i][dpC][b][temp_id] = dpt3;
				if (dps3 < st->ts.salt_missing_data)
					rawtempsalt[i][dpC][b][salt_id] = 0;
				else
					rawtempsalt[i][dpC][b][salt_id] = dps3;

				dp = 3;
				if (dpt4 < st->ts.temp_missing_data)
					rawtempsalt[i][dpD][b][temp_id] = 0;
				else
					rawtempsalt[i][dpD][b][temp_id] = dpt4;
				if (dps4 < st->ts.salt_missing_data)
					rawtempsalt[i][dpD][b][salt_id] = 0;
				else
					rawtempsalt[i][dpD][b][salt_id] = dps4;

			}
			break;
		}

		/* Close the exchange file */
		fclose(fp);
		//quit("");
	}

	if (verbose)
		printf("Finished reading raw temperature and salinity files\n");

	return;
}

/****************************************************************************************************
 Routine to check look-up table for conversion of faces in Mexican flat files to bgm face ids
 *****************************************************************************************************/
void Load_Face_Lookup_Table(Stuff *st) {
	int buflen = 500;
	char buf[500];
	FILE *fp;
	int baseboxID, fce;
	double proportion;
	int box1, box2;
	/* Read and load lookup table */
	if (verbose)
		printf("Reading lookup table file %s\n", st->lookupIfname);

	/* Open the lookup table file */
	if ((fp = fopen(st->lookupIfname, "r")) == NULL)
		quit("Load_Face_Lookup_Table - Lookup table file: Can't open %s\n", st->lookupIfname);

	/* rewind file */
	fseek(fp, 0L, SEEK_SET);

	/* Read in data */
	while (fgets(buf, buflen, fp) != NULL) {
		//printf("buf = %s\n", buf);
		if (strstr(buf, "Polygon") != NULL)
			continue;

		/* Check to see if the adjBox value contains a '&' - is so then we need to find both faces */
		/* Need to get the face from the bgm value */

		if (sscanf(buf, "%d,%d,%d&%d,%lf", &baseboxID, &fce, &box1, &box2, &proportion) == 5) {
			/* get each of the indices */
			col_faceid[baseboxID][fce][0] = getFaceIndex(st, baseboxID, box1);
			col_faceProportion[baseboxID][fce][0] = proportion;

			col_faceid[baseboxID][fce][1] = getFaceIndex(st, baseboxID, box2);
			col_faceProportion[baseboxID][fce][1] = 1.0 - proportion;
		} else if (sscanf(buf, "%d,%d,%d,%lf", &baseboxID, &fce, &box1, &proportion) == 4) {
			col_faceid[baseboxID][fce][0] = getFaceIndex(st, baseboxID, box1);
			col_faceProportion[baseboxID][fce][0] = 1.0;

			col_faceid[baseboxID][fce][1] = -100;
			col_faceProportion[baseboxID][fce][1] = 1.0;


			//printf("col_faceid[%d][%d][0] = %d\n",baseboxID, fce, col_faceid[baseboxID][fce][0]);

		} else
			quit("not enough entries for %s - with %d %d %s %f\n", buf, baseboxID, fce, box1, proportion);
	}

	/* Close the exchange file */
	fclose(fp);

	return;
}

/****************************************************************************************************
 Routine to load the depth_layer input file data into the depth_layer array.
 *****************************************************************************************************/
void Load_Layer_Table(Stuff *st) {
	int buflen = 500;
	char buf[500];
	FILE *fp;
	int dp;
	char	seps[] = ",";
	char *varStr;
	int box;
	double depth;

	/* Read and load lookup table */
	if (verbose)
		printf("Reading layer depth file %s\n", st->layerDepthIfname);

	/* Open the lookup table file */
	if ((fp = fopen(st->layerDepthIfname, "r")) == NULL)
		quit("Lookup layer depth file: Can't open %s\n", st->layerDepthIfname);

	/* rewind file */
	fseek(fp, 0L, SEEK_SET);

	/* Read in data */
	while (fgets(buf, buflen, fp) != NULL) {
		//printf("buf = %s\n", buf);

		varStr = strtok ( buf, seps );
		box = atoi(varStr);

		for(dp= 0; dp < st->wcnz; dp++){
			varStr = strtok ( NULL, seps );
			if(varStr == NULL || strstr(varStr, "-") != NULL){
				depth_layer[box][dp] = 0;
			}else{
				depth = atof(varStr);
				depth_layer[box][dp] = depth;

			}
			if(verbose > 3)
				printf("depth_layer[%d][%d] = %e\n", box, dp, depth_layer[box][dp]);

		}

	}
	/* Close the exchange file */
	fclose(fp);
//
	return;
}


/**
 *	\brief Will parse th egiven file name and pull out the month and year, the day value is set to 1.
 *
 *	This is used when input_style == col_data_type.
 *	File name should be like:
 *		1987_Apr_flux.dat or
 *		1986_Dec_Ave.dat
 *
 *
 */
void getDatesFromFileName(char *fileName, int *year, int *month, int *day) {
	char *strPtr;
	char monthStr[50];

	printf("fileName- %s\n", fileName);
#ifdef WIN32
	/* Get the last instance of '\' */
	strPtr = strrchr(fileName, '\\');
	/* Now we should just have the file name */
	sscanf(strPtr, "\\%d_", year);
#else
	/* Get the last instance of '\' */
	strPtr = strrchr(fileName, '/');
	/* Now we should just have the file name */
	sscanf(strPtr, "/%d_", year);
#endif

	strPtr = strstr(strPtr, "_");
	strPtr++;
	strcpy(monthStr, strPtr);
	*(strstr(monthStr, "_")) = '\0';
	*month = getMonthIndex(monthStr);
	*day = 1;
}



/**
 *	\brief Will parse th egiven file name and pull out the month and year, the day value is set to 1.
 *
 *	This is used when input_style == new_col_data_type.
 *	File name should be like:
 *		/home/bec/BackedUp/Code/Matlab/trunk/Private/GOM/gom_all_out_2010_ias/2010-01/fileName
 *
 *	So we want to pull the year and month from the path instead of the filename.
 *
 */
void getDatesFromPathName(char *fileName, int *year, int *month, int *day) {
	char *strPtr;

	printf("getDatesFromPathName - fileName- %s\n", fileName);
#ifdef WIN32


	strPtr = strrchr(fileName, '\\');
	printf("strPtr = %s\n", strPtr);
	//strPtr++;
	*strPtr = '\0';

	printf("fileName = %s\n", fileName);
	/* Now find the next last one */
	strPtr = strrchr(fileName, '\\');

	printf("strPtr = %s\n",strPtr);

	/* Now we should just have the file name */
	sscanf(strPtr, "\\%d-", year);

	printf("year = %d\n", *year);

#else
	/* Get the last instance of '\' */
	strPtr = strrchr(fileName, '/');
	printf("strPtr = %s\n", strPtr);
	//strPtr++;
	*strPtr = '\0';

	printf("fileName = %s\n", fileName);
	/* Now find the next last one */
	strPtr = strrchr(fileName, '/');

	/* Now we should just have the file name */
	sscanf(strPtr, "/%d_", year);
#endif

	strPtr = strstr(strPtr, "-");
	strPtr++;

	printf("strPtr = %s\n", strPtr);

	sscanf(strPtr, "%d", month);
	month++;

	//strcpy(monthStr, strPtr);
	//*(strstr(monthStr, "_")) = '\0';
	//*month = getMonthIndex(monthStr);
	*day = 1;
}



int getFaceIndex(Stuff *st, int baseboxID, int adjBox) {
	int index;

	//printf("adjBox = %d\n", adjBox);

	for (index = 0; index < st->boxnconn[baseboxID]; index++) {
		//printf("st->ibox[baseboxID][index] = %d\n", st->ibox[baseboxID][index]);
		if (st->ibox[baseboxID][index] == adjBox) {
			return st->iface[baseboxID][index];
		}
	}
	return -100;
}
int getMonthIndex(char *monthStr) {
	int index;

	printf("monthStr = %s\n", monthStr);
	for (index = 0; index < 12; index++)
		if (strcmp(monthStr, months[index]) == 0)
			return index + 1;

	return index;
}

/**
 * Trim leading and trailing white space.
 */

void trim(char *s)
{
	// Trim spaces and tabs from beginning:
	int i=0,j;
	while((s[i]==' ')||(s[i]=='\t') ||(s[i]=='\n')) {
		i++;
	}
	if(i>0) {
		for(j=0;j<(int)strlen(s);j++) {
			s[j]=s[j+i];
		}
		s[j]='\0';
	}

	// Trim spaces and tabs from end:
	i=strlen(s)-1;

	while((s[i]==' ')||(s[i]=='\t') || (s[i]=='\n')|| (s[i]=='\r')) {
		i--;
	}
	if(i<((int)strlen(s)-1)) {
		s[i + 1]='\0';
	}
}

