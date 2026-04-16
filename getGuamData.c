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
void get_guam_hydro(Stuff *st) {
	FILE *fp;
	int i, old_nline;
	int fce, dp, nline, ti, fce_id, dbox, newdatid;
	double dp1, dp2, dp3, dp4, dp5, dummy;
	int entry_i = 0;
	int buflen = 500;
	char buf[500];

	st->minExyear = MAXINT;

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


		/* Read in data */
		nline = 0;
		//printf("got ot here\n");
		/* Read in the first line and ignore */
		if(fgets(buf, buflen, fp) == NULL){
			quit("error in read_Hydro_New_col_data(). Problem loaded data form file\n");
		}

		/* Read in each line */
		while (fgets(buf, buflen, fp) != NULL) {
			//printf("buf = %s\n", buf);
			if (strstr(buf, "Polygon") != NULL || strstr(buf, "No.") != NULL || strstr(buf, "transect") != NULL || strstr(buf, "number") != NULL|| strlen(buf) == 0 || strcmp(
					buf, "\n") == 0 || strstr(buf, "Box_ID") != NULL)
				continue;

			if (sscanf(buf, "%lf,%lf,%lf,%lf,%lf,%lf", &dummy, &dp1, &dp2, &dp3, &dp4, &dp5) != 6)
				quit("not enough entries for '%s' - with %lf %lf %lf %lf %lf %lf\n", buf, dummy, dp1, dp2, dp3, dp4, dp5);

			newdatid = (int) (dp1);
			ti = old_nline + newdatid - 1; // New time step of entry based on base date + value from column 3 (minus 1 as values start at 1 not 0)

			exchangedate[ti][time_id] = ti;
			exchangedate[ti][done_id] = 1;
			exchangedate[ti][TofY_id] = ti;


			nline++;
			dbox = (int) (dp3); // BoxID from column 1
			fce = (int) (dp4) - 1; // FaceID from column 2 and then use lookup table to do the translation
			dp = (int)dp2 - 1;

			fce_id = col_faceid[dbox][fce][0];

			if (fce_id > -1){// Don't store flows from boundary box out of the system
				if(verbose > 3){
					printf("st->boxnz[%d] = %d, st->boxtotdepth[dbox] = %e, total flux = %e \n",
							dbox,
							st->boxnz[dbox],
							st->boxtotdepth[dbox],
							dp5			);
				}

				rawexchanges[ti][dp][fce_id][entry_i] = dp5;
				if(verbose > 3)
					printf("box %d, link %d, rawexchanges[%d][0][%d][%d]= %e\n", dbox, fce, ti, fce_id, entry_i, rawexchanges[ti][0][fce_id][entry_i]);

			} // Data from column 4

		}

		old_nline += nline;

		/* Close the exchange file */
		fclose(fp);
	}

	if (verbose)
		printf("Finished reading raw exchange files\n");
	return;
}


 /***************************************************************************************************
  Routine to get the temperature and salinity data out of the raw data files - - Jeff Dunn style flat
  input (not netcdf files)
  ****************************************************************************************************/
 void get_guam_temp_and_salt(Stuff *st) {
 	FILE *fp;
 	int i, b, old_nline, nline, newdatid, ti, dpA, entry_i;
 	double dpt1, dps1, dp1, dp2, dp3, dp4, dummy;
 	int buflen = 500;
 	char buf[500];
 	int maxTime = -1;
 	//double depth;
 	int lineIndex;

 	dpA = 0;

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


 		tempsaltdate[i][time_id] = i;
 		tempsaltdate[i][TofY_id] = i;


		/* This file was originally for Cam in Mexico and assumes
		 the hydro file has the format

		 date year month day
		 dt t
		 nline y
		 LineNum Time_step  Polygon  Depth   Vertical_velocity[10^-5*m/s] Average_Temp[Celsius]  Average_Salinity[PartsPerThousand]
		 linei	   i		j			x			vali					valti					valsi

		 BEC TO FIX - CAN NOT TAG HEAD OF EACH FILE AS TOO MANY SO HOW GET TAG OUT OF FILE NAME? HOW KNOW FILES TO LOAD?
		 */

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
			if (strstr(buf, "depthLayer") != NULL || strstr(buf, "transect") != NULL || strlen(buf) == 0 || strcmp(buf, "\n") == 0 || strstr(buf,
					"Layer") != NULL || strstr(buf, "PartsPerThousand") != NULL || strstr(buf, "PartsPer1000") != NULL)

				continue;

			if (sscanf(buf, "%lf,%lf,%lf,%lf,%lf,%lf,%lf", &dummy, &dp1, &dp2, &dp3, &dpt1, &dps1, &dp4) != 7){
				fprintf(stderr, "Line %d - %s in file %s is not the correct format\n", lineIndex, buf, st->ts.tsfname[i]);
				quit("\n\n get_guam_temp_and_salt: not enough entries for %s - with %lf %lf %lf %lf %lf %lf in file %s\n", buf, dp1, dp2, dp3, dp4, dpt1,
						dps1, st->ts.tsfname[i]);
			}
			lineIndex++;

			newdatid = (int) (dp1);
			ti = old_nline + newdatid - 1; // New time step of entry based no base date + value from column 1 (minus 1 as values start at 1 not 0)

			if(ti  >= st->nvhdsteps){
				printf("line %s\n", buf);
				quit("The number of timesteps has been calculated as %d but in the file %s we have the %dth data point\n", st->nvhdsteps, st->ts.tsfname[i], ti);
			}
			b = (int) (dp3); // BoxID from column 2


			/* Find layer from Depth in column 3 */
			//dp = st->boxnz[b];
			dpA = (int)dp2 - 1;

//			depth = 0;
//			for (k = 0; k < st->wcnz; k++) {
//				depth  = depth +  st->nominal_dz[k];
//				//printf("st->boxdz[b][%d] = %f, depth= %e\n", k, st->nominal_dz[k], depth);
//				if (dp3 <= depth){
//					dpA = dp - 1 - k; // Tip layers upside down (so layer 0 next to sediment) in calcexchange
//					dpA = k; // Tip layers upside down (so layer 0 next to sediment) in calcexchange
//					break;
//				}
//			}
//
//			if(dpA > (dp - 1)){
//			//	abort();
//				continue;
//			}

			/* Enter raw temperature and salinity values */
			if (dpt1 < st->ts.temp_missing_data)
				rawtempsalt[ti][dpA][b][temp_id] = 0;
			else
				rawtempsalt[ti][dpA][b][temp_id] = dpt1;




			if (dps1 < st->ts.salt_missing_data)
				rawtempsalt[ti][dpA][b][salt_id] = 0;
			else
				rawtempsalt[ti][dpA][b][salt_id] = dps1;

			if(b == 16){
				rawtempsalt[ti][1][b][temp_id] = dpt1;
				rawtempsalt[ti][1][b][salt_id] = dps1;
			}


			/* Enter raw vertical exchange values */
			if (dps1 < st->hd.missing_data){
				rawvertexchanges[ti][dpA][b][entry_i] = 0;
			}else{
				rawvertexchanges[ti][dpA][b][entry_i] = dp4 ;//* 10e5;
			}

			//rawvertexchanges[ti][dpA][b][entry_i] = 0;

			if (newdatid > maxTime) {
				nline++;
				maxTime = newdatid;
			}
		}
		old_nline += nline;

 		/* Close the exchange file */
 		fclose(fp);
 	}

 	if (verbose)
 		printf("Finished reading raw temperature and salinity files\n");

 	return;
 }

