/**********************************************************************

 getinputdata.c

 by:	Beth Fulton
 date:	2/3/2006

 comments: This file is the code that reads in netcdf file version of
 raw hydrodynamics flows, and temperature and salinity profiles

 assumptions: assumes raw input is a series of data entries in netcdf file

 Changes:

 12-05-2009 Bec Gorton

 The new netcdf files form Jeff contain the bgm face id's so there is no need to
 use the face geometry to match the input file faces with the bgm faces.
 So all of the facepx facepy code has been removed.

 02-06-2009 Bec Gorton
 Got rid of some commented out code thats not used anymore.
 03-06-2009 Bec Gorton
 Added the load_bgm_mappings function to load the mappings between the hydro box numbers
 and the bgm box numbers - this functionality was added to allow users to change the box
 numbers after the hydro files have been generated.

 07-06-2009 Bec Gorton
 Changed the code that reads in the input files to always allocate floats for the exchange array.
 The values in the input file are always floats - the valtime array is still allocated based on the current
 definition of FPTYPE.
 *************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <math.h>
#include <string.h>
#include <netcdf.h>
#include <sjwlib.h>
#include "mainexchange.h"
#include "unixUtils.h"

/***************************************************************************************************
 Routine to get the exchange data out of the raw data file - Jeff Dunn style netcdf input (not netcdf files)
 Involves reading data from netcdf file to buffer
 ****************************************************************************************************/
void get_cdf_hydro(Stuff *st) {
	int i, k, lastk, entry_i, nz, ij;
	double time_of_year, totnum, this_tofy, mintime;
	float ***valexchange;
	FPTYPE *valtime;
	int *valboxs;
	int *valboxd;
	int *valfaceid;
	long start[3];
	long count[3];
	long startt[1];
	long countt[1];
	long startfce[1];
	long countfce[1];
	nc_type datatype;
	long nt, fce, lvl;

	ncopts = NC_VERBOSE | NC_FATAL;

	if (verbose)
		printf("Loading raw data from netcdf files\n");

	/* Give code warning */
	if (st->hd.nfiles > 1)
		warn(
				"%d raw hydro files listed, code may not cope properly with more than one - in particular the face mapping may be wrong\n",
				st->hd.nfiles);

	/* Loop through the exchange files */
	lastk = 0;
	for (i = 0; i < st->hd.nfiles; i++) {
		if (verbose)
			printf("Reading raw exchange file %s\n", st->hd.fname[i]);

		/* Open the exchange file */
		open_hydro(st, st->hd.fname[i], i);

		if (verbose > 1)
			printf("Checking for dimensions\n");

		/* Check for dimensions in file */
		ncdiminq(st->hd.fid, ncdimid(st->hd.fid, "time"), NULL, &nt);
		ncdiminq(st->hd.fid, ncdimid(st->hd.fid, "faces"), NULL, &fce);
		ncdiminq(st->hd.fid, ncdimid(st->hd.fid, "level"), NULL, &lvl);

		if (fce > st->nbgmface) {
			warn(
					"More faces in netcdf file %s (nface = %d) than in bgm file %s (nface = %d) - reconcile and retry\n",
					st->hd.fname[i], fce, st->geomIfname, st->nface);
		} else if (fce < st->nbgmface) {
			warn(
					"Less faces in netcdf file %s (nface = %d) than in bgm file %s (nface = %d)\n",
					st->hd.fname[i], fce, st->geomIfname, st->nface);
		}

		/* Allocate temporary storage for one tracer */
		nc_inq_vartype(st->hd.fid, st->hd.fcex1_vid, &datatype);

		if (datatype != NC_FLOAT)
			quit("The data type in the netcdf must be NC_FLOAT\n");

		valexchange = (float ***) f_alloc3d(lvl, fce, nt);

		switch (sizeof(FPTYPE)) {
		case sizeof(float):
			valtime = (FPTYPE *) f_alloc1d(nt);
			break;
		case sizeof(double):
			valtime = (FPTYPE *) d_alloc1d(nt);
			break;
		default:
			quit("get_cdf_hydro: Unknown size for FPTYPE\n");
		}
		valboxs = (int *) i_alloc1d(fce);
		valboxd = (int *) i_alloc1d(fce);
		valfaceid = (int *) i_alloc1d(fce);

		/* Set indices for reading tracers */
		start[0] = 0;
		start[1] = 0;
		start[2] = 0;
		count[0] = nt;
		count[1] = fce;
		count[2] = lvl;

		startfce[0] = 0;
		countfce[0] = fce;

		startt[0] = 0;
		countt[0] = nt;

		/* Read data - exchanges, dates and box-ids associated with the faces */
		ncvarget(st->hd.fid, st->hd.t_vid, startt, countt, valtime);

		ncvarget(st->hd.fid, st->hd.fce_vid, startfce, countfce, valfaceid);
		ncvarget(st->hd.fid, st->hd.srcb_vid, startfce, countfce, valboxs);
		ncvarget(st->hd.fid, st->hd.destb_vid, startfce, countfce, valboxd);
		ncvarget(st->hd.fid, st->hd.e_vid, start, count, valexchange[0][0]);

		printf("st->dt = %e\n", st->dt);

		/* Transfer to rawdata arrays */
		/* Transfer timestamps - while says seconds in input file actually in days since reference date */
		if (st->reset_time) {
			mintime = MAXDOUBLE;
			for (k = 0; k < nt; k++) {
				if (valtime[k] < mintime)
					mintime = valtime[k];
			}
			//st->time_reset = mintime * st->dt;
			st->time_reset = mintime * 86400;
		}

		printf("st->time_reset = %e\n", st->time_reset);
		for (k = 0; k < nt; k++) {

			if (st->verbose_pad)
				printf(
						"k = %d,  valtime[k] = %e, st->dt = %e, st->time_reset = %e\n",
						k, valtime[k], st->dt, st->time_reset);
			//exchangedate[k + lastk][time_id] = valtime[k] * st->dt - st->time_reset;
			exchangedate[k + lastk][time_id] = valtime[k] * 86400
					- st->time_reset;

			/* Multiple by 86400 / dt as may be tidal entries so must expand
			 list to fit these extra entries in */
			if(st->dt < 86400){
				time_of_year = valtime[k] * (86400 / st->dt);
				totnum = 365.0 * (86400 / st->dt);

				if (time_of_year / totnum < 1.0)
					this_tofy = floor(time_of_year / (86400 / st->dt) + 0.5);
				else
					this_tofy = floor((time_of_year / totnum - floor(time_of_year
							/ totnum)) * totnum + 0.5);

			}else{
				time_of_year = valtime[k];// * (st->dt / 86400.0);
				totnum = 365.0 * (st->dt/86400.0);


				if ((time_of_year / totnum) < 1.0)
					this_tofy = floor(time_of_year / (st->dt/86400.0) + 0.5);
				else
					this_tofy = floor((time_of_year / totnum - floor(time_of_year
							/ totnum)) * totnum + 0.5);

			}
			exchangedate[k + lastk][TofY_id] = this_tofy;

			if (st->verbose_pad)
				printf("k: %d, time: %e, tofy: %e valtime: %e\n", k,
						exchangedate[k + lastk][time_id], exchangedate[k
								+ lastk][TofY_id], valtime[k]);
		}

		/* Transfer boxes */
		for (k = 0; k < fce; k++) {
			/* Swap the box values to the bgm box values using the preloaded boxLookup array.
			 * These values have been loaded in the load_bgm_mappings function.
			 * The boxLookup array simply stores the bgm box id in the hydro box value index.
			 * */
			faceinfo[i][k][lbox_id] = boxLookup[valboxs[k]];
			faceinfo[i][k][rbox_id] = boxLookup[valboxd[k]];


			faceinfo[i][k][fceid_id] = valfaceid[k];

			facemap[k] = valfaceid[k];
			//facemap[k] = k;
		}

		/* Transfer exchanges */
		entry_i = 0;
		for (k = 0; k < nt; k++) {
			for (nz = 0; nz < lvl; nz++) {
				for (ij = 0; ij < fce; ij++) {


					if(valexchange[k][ij][nz] > 1e15){
						printf("horizontal exchange value is too large. valexchange[%d][%d][%d] = %e\n", k, ij, nz, valexchange[k][ij][nz]);
						quit("");
					}

					if (isnan(valexchange[k][ij][nz]) || !_finite(valexchange[k][ij][nz]))
						rawexchanges[k + lastk][nz][ij][entry_i] = 0.0;
					else
						rawexchanges[k + lastk][nz][ij][entry_i]
								= valexchange[k][ij][nz];

					if(ij == 113){
						printf("rawexchanges[k + lastk][nz][ij][entry_i] = %e\n",rawexchanges[k + lastk][nz][ij][entry_i]);
					}
				}
			}
		}

		/* Set lastk */
		lastk += nt;

		/* Free temporary storage */
		switch (sizeof(FPTYPE)) {
		case sizeof(float):
			f_free1d((float *) valtime);
			break;
		case sizeof(double):
			d_free1d((double *) valtime);
			break;
		default:
			quit("get_cdf_hydro: Unknown size for FPTYPE\n");
		}

		f_free3d((float ***) valexchange);
		i_free1d((int *) valboxs);
		i_free1d((int *) valboxd);
		i_free1d((int *) valfaceid);

		/* Close hydro file */
		close_hydro(st);

	}

	if (verbose)
		printf("Finished reading raw exchange files\n");

	return;
}

/****************************************************************************************************
 Routine to map box and face ids from netcdf file to matching ones in bgm file
 *****************************************************************************************************/
void map_bgm_netcdf_faces(Stuff *st, int i) {
	int k, fce, actfce, bgmlbox, bgmrbox, cdflbox,
			cdfrbox, p1forp1, p1forp2, p2forp2, p2forp1,
			positive_equiv;
	//, matchscrew = 0;
	//, foundmatch = 0;

	/* Do face assignment */
	for (fce = 0; fce < st->nface; fce++) {
		facemap[fce] = -1;
		//foundmatch = 0;
		//matchscrew = 0;
		positive_equiv = -1;
		for (k = 0; k < st->nbgmface; k++) {
			/* Initialise markers */
			p1forp1 = 0;
			p1forp2 = 0;
			p2forp1 = 0;
			p2forp2 = 0;

			if (verbose > 1) {
				printf("fce (%e, %e)(%e, %e) vs bgmfce (%e, %e)(%e, %e)\n",
						faceinfo[i][fce][px1_id], faceinfo[i][fce][py1_id],
						faceinfo[i][fce][px2_id], faceinfo[i][fce][py2_id],
						st->faces[k].pll1.x, st->faces[k].pll1.y,
						st->faces[k].pll2.x, st->faces[k].pll2.y);
			}

			/* As ordering not always identical first check against one pair then the other
			 then decide if ordering needs to be reversed or not */
			/* Check p1 vs pll1 */
			if (fabs(faceinfo[i][fce][px1_id] - st->faces[k].pll1.x)
					< buffer_rounding) {
				if (fabs(faceinfo[i][fce][py1_id] - st->faces[k].pll1.y)
						< buffer_rounding) {
					p1forp1 = 1;
					//if(fce == 13)
					printf("p1p1 with %d\n", k);
				}
			}

			/* Check p1 vs pll2 */
			if (fabs(faceinfo[i][fce][px1_id] - st->faces[k].pll2.x)
					< buffer_rounding) {
				if (fabs(faceinfo[i][fce][py1_id] - st->faces[k].pll2.y)
						< buffer_rounding) {
					p1forp2 = 1;
					//if(fce == 13)
					printf("p1p2 with %d\n", k);
				}
			}

			/* Check p2 vs pll1 */
			if (fabs(faceinfo[i][fce][px2_id] - st->faces[k].pll1.x)
					< buffer_rounding) {
				if (fabs(faceinfo[i][fce][py2_id] - st->faces[k].pll1.y)
						< buffer_rounding) {
					p2forp1 = 1;
					//if(fce == 13)
					printf("p2p1 with %d\n", k);
				}
			}

			/* Check p2 vs pll2 */
			if (fabs(faceinfo[i][fce][px2_id] - st->faces[k].pll2.x)
					< buffer_rounding) {
				if (fabs(faceinfo[i][fce][py2_id] - st->faces[k].pll2.y)
						< buffer_rounding) {
					p2forp2 = 1;
					//if(fce == 13)
					printf("p2p2 with %d\n", k);
				}
			}

			if (p1forp1 && p2forp2) {
				/* All ordinates match so assign positive equivalency */
				if (facemap[fce] > -1) {
					printf(
							"Already assigned face %d to bgmface: %d now reassigning to %d - is this right?\n",
							fce, facemap[fce], k);
					//matchscrew = 1;
				}
				facemap[fce] = k;
				positive_equiv = 1;
				//foundmatch = 1;

				printf("Assigning face %d to bgmface: %d\n", fce, k);

			} else if (p1forp2 && p2forp1) {
				/* All ordinates match but points back to front from one file to the other (p1 in one
				 file is p2 in the other) then assign negative equivalency */
				if (facemap[fce] > -1) {
					printf(
							"Already assigned face %d to bgmface: %d now reassigning to %d - is this right?\n",
							fce, facemap[fce], k);
					//matchscrew = 1;
				}
				facemap[fce] = k;
				positive_equiv = 0;
				//foundmatch = 1;
				printf("Assigning face %d to bgmface: %d\n", fce, k);
			}

			/**
			 if(foundmatch or matchscrew){
			 if(verbose > 1){
			 printf("for file %d, checking fce %d p1(%e,%e) and p2(%e,%e) vs bgmfce %d pll1(%e,%e) and pll2(%e,%e)\n",
			 i, fce, faceinfo[i][fce][px1_id], faceinfo[i][fce][py1_id], faceinfo[i][fce][px2_id], faceinfo[i][fce][py2_id], k, st->faces[k].pll1.x, st->faces[k].pll1.y, st->faces[k].pll2.x, st->faces[k].pll2.y);
			 printf("do p1-p1 xcheck: %e, ycheck: %e (vs %e)\n", fabs(faceinfo[i][fce][px1_id] - st->faces[k].pll1.x), fabs(faceinfo[i][fce][py1_id] - st->faces[k].pll1.y), buffer_rounding);
			 printf("do p1-p2 xcheck: %e, ycheck: %e (vs %e)\n", fabs(faceinfo[i][fce][px1_id] - st->faces[k].pll2.x), fabs(faceinfo[i][fce][py1_id] - st->faces[k].pll2.y), buffer_rounding);
			 printf("do p2-p1 xcheck: %e, ycheck: %e (vs %e)\n", fabs(faceinfo[i][fce][px2_id] - st->faces[k].pll1.x), fabs(faceinfo[i][fce][py2_id] - st->faces[k].pll1.y), buffer_rounding);
			 printf("do p2-p2 xcheck: %e, ycheck: %e (vs %e)\n", fabs(faceinfo[i][fce][px2_id] - st->faces[k].pll2.x), fabs(faceinfo[i][fce][py2_id] - st->faces[k].pll2.y), buffer_rounding);
			 }
			 }
			 **/
		}
		if (facemap[fce] < 0) {
			warn("Face %d does not have a cdf <-> bgm mapping\n", fce);

			if (verbose > 1) {
				for (k = 0; k < st->nbgmface; k++) {
					printf(
							"for file %d, checking fce %d p1(%e,%e) and p2(%e,%e) vs bgmfce %d pll1(%e,%e) and pll2(%e,%e)\n",
							i, fce, faceinfo[i][fce][px1_id],
							faceinfo[i][fce][py1_id], faceinfo[i][fce][px2_id],
							faceinfo[i][fce][py2_id], k, st->faces[k].pll1.x,
							st->faces[k].pll1.y, st->faces[k].pll2.x,
							st->faces[k].pll2.y);
					printf(
							"do p1-p1 xcheck: %e, ycheck: %e (vs %e)\n",
							fabs(faceinfo[i][fce][px1_id] - st->faces[k].pll1.x),
							fabs(faceinfo[i][fce][py1_id] - st->faces[k].pll1.y),
							buffer_rounding);
					printf(
							"do p1-p2 xcheck: %e, ycheck: %e (vs %e)\n",
							fabs(faceinfo[i][fce][px1_id] - st->faces[k].pll2.x),
							fabs(faceinfo[i][fce][py1_id] - st->faces[k].pll2.y),
							buffer_rounding);
					printf(
							"do p2-p1 xcheck: %e, ycheck: %e (vs %e)\n",
							fabs(faceinfo[i][fce][px2_id] - st->faces[k].pll1.x),
							fabs(faceinfo[i][fce][py2_id] - st->faces[k].pll1.y),
							buffer_rounding);
					printf(
							"do p2-p2 xcheck: %e, ycheck: %e (vs %e)\n",
							fabs(faceinfo[i][fce][px2_id] - st->faces[k].pll2.x),
							fabs(faceinfo[i][fce][py2_id] - st->faces[k].pll2.y),
							buffer_rounding);
				}
			}
			//quit("");
		}

		if (verbose > 1) {
			printf("Map fce %d <-> bgmfce: %d\n", fce, facemap[fce]);
		}

		/* Do box assignment */
		actfce = facemap[fce];
		if (actfce < 0)
			/* If not a real face nothing to do so skip ahead */
			continue;

		/* Match faces */
		bgmlbox = st->faces[actfce].ibl;
		bgmrbox = st->faces[actfce].ibr;

		printf("bgmlbox = %d\n", bgmlbox);
		printf("bgmrbox = %d\n", bgmrbox);

		cdflbox = (int) (faceinfo[i][fce][lbox_id]);
		cdfrbox = (int) (faceinfo[i][fce][rbox_id]);

		printf("cdflbox = %d\n", cdflbox);
		printf("cdfrbox = %d\n", cdfrbox);

		if (verbose > 1) {
			if (!positive_equiv) {
				printf(
						"+ve equiv: %d fce%d (l-%d, r-%d) vs bgmfce%d (l-%d, r-%d)\n",
						positive_equiv, fce, cdflbox, cdfrbox, actfce, bgmlbox,
						bgmrbox);
			} else {
				printf(
						"-ve equiv: %d fce%d (l-%d, r-%d) vs bgmfce%d (l-%d, r-%d)\n",
						positive_equiv, fce, cdflbox, cdfrbox, actfce, bgmlbox,
						bgmrbox);

			}
			//printf("--> fp(%e, %e | %e, %e) vs bp(%e, %e | %e, %e)\n",
			//	faceinfo[i][fce][px1_id], faceinfo[i][fce][py1_id], faceinfo[i][fce][px2_id], faceinfo[i][fce][py2_id], st->faces[actfce].pll1.x, st->faces[actfce].pll1.y, st->faces[actfce].pll2.x, st->faces[actfce].pll2.y);

		}

		/* Left side box */
		if (positive_equiv && boxmap[cdflbox] < 0) {
			/* Positive equivalency and previously unassigned */
			printf("Assigning lbox %d to lbgmbox: %d\n", cdflbox, bgmlbox);
			boxmap[cdflbox] = bgmlbox;
		} else if (positive_equiv && boxmap[cdflbox] != bgmlbox) {
			/* Positive equivalence, but previously assigned */
			printf(
					"Already assigned lbox %d to lbgmbox: %d now reassigning to %d - is this right?\n",
					cdflbox, boxmap[cdflbox], bgmlbox);
			printf(
					"+ve_equiv: %d, fce%d p(%e, %e | %e, %e) vs bgmfce%d p (%e, %e | %e, %e)\n",
					positive_equiv, fce, faceinfo[i][fce][px1_id],
					faceinfo[i][fce][py1_id], faceinfo[i][fce][px2_id],
					faceinfo[i][fce][py2_id], actfce, st->faces[actfce].pll1.x,
					st->faces[actfce].pll1.y, st->faces[actfce].pll2.x,
					st->faces[actfce].pll2.y);
			boxmap[cdflbox] = bgmlbox;
		} else if (!positive_equiv && boxmap[cdfrbox] < 0) {
			/* Positive equivalency and previously unassigned */
			printf("Assigning rbox %d to lbgmbox: %d\n", cdfrbox, bgmlbox);
			boxmap[cdfrbox] = bgmlbox;
		} else if (!positive_equiv && boxmap[cdfrbox] != bgmlbox) {
			/* Positive equivalence, but previously assigned */
			printf(
					"Already assigned rbox %d to lbgmbox: %d now reassigning to %d - is this right?\n",
					cdfrbox, boxmap[cdfrbox], bgmlbox);
			printf(
					"+ve_equiv: %d, fce%d p(%e, %e | %e, %e) vs bgmfce%d p (%e, %e | %e, %e)\n",
					positive_equiv, fce, faceinfo[i][fce][px1_id],
					faceinfo[i][fce][py1_id], faceinfo[i][fce][px2_id],
					faceinfo[i][fce][py2_id], actfce, st->faces[actfce].pll1.x,
					st->faces[actfce].pll1.y, st->faces[actfce].pll2.x,
					st->faces[actfce].pll2.y);
			boxmap[cdfrbox] = bgmlbox;
		}

		/* Right side box */
		if (positive_equiv && boxmap[cdfrbox] < 0) {
			/* Positive equivalency and previously unassigned */
			//printf("Assigning rbox %d to rbgmbox: %d\n", cdfrbox, bgmrbox);
			boxmap[cdfrbox] = bgmrbox;
		} else if (positive_equiv && boxmap[cdfrbox] != bgmrbox) {
			/* Positive equivalence, but previously assigned */
			//printf("Already assigned rbox %d to rbgmbox: %d now reassigning to %d - is this right?\n", cdfrbox, boxmap[cdfrbox], bgmrbox);
			//printf("+ve_equiv: %d, fce%d p(%e, %e | %e, %e) vs bgmfce%d p (%e, %e | %e, %e)\n",
			//	positive_equiv, fce, faceinfo[i][fce][px1_id], faceinfo[i][fce][py1_id], faceinfo[i][fce][px2_id], faceinfo[i][fce][py2_id], actfce, st->faces[actfce].pll1.x, st->faces[actfce].pll1.y, st->faces[actfce].pll2.x, st->faces[actfce].pll2.y);
			boxmap[cdfrbox] = bgmrbox;
		} else if (!positive_equiv && boxmap[cdflbox] < 0) {
			/* Positive equivalency and previously unassigned */
			//printf("Assigning lbox %d to rbgmbox: %d\n", cdflbox, bgmrbox);
			boxmap[cdflbox] = bgmrbox;
		} else if (!positive_equiv && boxmap[cdflbox] != bgmrbox) {
			/* Positive equivalence, but previously assigned */
			//printf("Already assigned lbox %d to rbgmbox: %d now reassigning to %d - is this right?\n", cdflbox, boxmap[cdflbox], bgmrbox);
			//printf("+ve_equiv: %d, fce%d p(%e, %e | %e, %e) vs bgmfce%d p (%e, %e | %e, %e)\n",
			//	positive_equiv, fce, faceinfo[i][fce][px1_id], faceinfo[i][fce][py1_id], faceinfo[i][fce][px2_id], faceinfo[i][fce][py2_id], actfce, st->faces[actfce].pll1.x, st->faces[actfce].pll1.y, st->faces[actfce].pll2.x, st->faces[actfce].pll2.y);
			boxmap[cdflbox] = bgmrbox;
		}
		printf("boxmap[%d] = %d\n", cdflbox, boxmap[cdflbox]);
		printf("boxmap[%d] = %d\n", cdfrbox, boxmap[cdfrbox]);

		//quit("");
	}

	/* Check box assignments */
	for (k = 0; k < st->nbox; k++) {
		if (boxmap[k] < 0)
			warn(
					"Box %d does not have a cdf <-> bgm mapping, boxmap[k] = %d\n",
					k, boxmap[k]);
	}

	return;
}

/****************************************************************************************************
 Routine to open netcdf file containing raw hydrodynamic fluxes
 *****************************************************************************************************/
void open_hydro(Stuff *st, char *name, int filenum) {
	int ndims;
	int nvars;
	int natts;
	int recdim;
	int len;
	long n;
	nc_type daty;
	int dims[MAX_NC_DIMS];
	char hdu[STSLEN];
	char stu[STSLEN];

	/* Set netCDF library error handling */
	ncopts = NC_VERBOSE;

	if (verbose > 1)
		printf("open_hydro: opening the file\n");

	/* Open the file */
	if ((st->hd.fid = ncopen(name, NC_NOWRITE)) < 0)
		quit("open_hydro: Can't open hydrodynamic model input data file %s\n",
				name);

	/* Inquire about this file */
	ncopts = NC_VERBOSE | NC_FATAL;
	ncinquire(st->hd.fid, &ndims, &nvars, &natts, &recdim);
	if (ndims < 3)
		quit("open_hydro: not enough dimensions in %s\n", name);

	if (verbose > 1)
		printf("open_hydro: check dimensions\n");

	/* Check dimensions are as expected */
	printf("check for time in standard CSIRO cdf\n");

	if ((st->hd.t_did = ncdimid(st->hd.fid, "time")) == -1)
		quit("open_hydro: no t dimension in %s\n", name);
	if (st->hd.t_did != recdim)
		quit("open_hydro: t dimension not unlimited in %s\n", name);

	printf("check for faces\n");

	if ((st->hd.b_did = ncdimid(st->hd.fid, "faces")) == -1) {
		quit("open_hydro: no faces dimension in %s\n", name);
	}

	printf("check for level\n");

	if ((st->hd.z_did = ncdimid(st->hd.fid, "level")) == -1)
		quit("open_hydro: no levels dimension in %s\n", name);

	/* Get dimension sizes. Note check against geometry now happens in get_cdf_hydro() */

	if (verbose > 1)
		printf("open_hydro: get dimnensionn\n");

	ncdiminq(st->hd.fid, st->hd.b_did, NULL, &n);
	ncdiminq(st->hd.fid, st->hd.z_did, NULL, &n);

	/* Check that time units and steps match this model */

	if (verbose > 1)
		printf("open_hydro: get units\n");

	st->hd.t_vid = ncvarid(st->hd.fid, "time");
	if (st->hd.t_vid < 0)
		quit("open_hydro: no t variable in %s\n", name);
	memset(st->hd.t_units, 0, STSLEN);
	ncattget(st->hd.fid, st->hd.t_vid, "units", st->hd.t_units);
	sscanf(st->hd.t_units, "%s", hdu);
	if (!filenum)
		sprintf(st->t_units, "%s", st->hd.t_units);
	sscanf(st->t_units, "%s", stu);
	if (strcmp(hdu, stu) != 0) {
		quit(
				"open_hydro: Time units (%s) in %s don't match model time units (%s) loaded from %s\n",
				hdu, name, stu, st->hd.fname[0]);
	}

	/* Check for timestep size */
	if (verbose > 1)
		printf("open_hydro: get timestep\n");

	ncattinq(st->hd.fid, st->hd.t_vid, "dt", &daty, &len);
	if (nctypelen(daty) != sizeof(st->hd.dt))
		quit("open_hydro: dt attribute wrong type\n");
	ncattget(st->hd.fid, st->hd.t_vid, "dt", &st->hd.dt);
	if (st->hd.dt <= 0.0)
		quit("open_hydro: dt attribute must have positive value\n");

	/* Find out how many time steps are in the file */
	ncdiminq(st->hd.fid, st->hd.t_did, NULL, &st->hd.nstep);

	/* Get other variable ids */
	if (verbose > 1)
		printf("open_hydro: get variable ids\n");

	st->hd.srcb_vid = ncvarid(st->hd.fid, "source_boxid");
	if (st->hd.srcb_vid < 0)
		quit("open_hydro: no source_boxid variable in %s\n", name);
	st->hd.destb_vid = ncvarid(st->hd.fid, "dest_boxid");
	if (st->hd.destb_vid < 0)
		quit("open_hydro: no dest_boxid variable in %s\n", name);

	st->hd.fce_vid = ncvarid(st->hd.fid, "faces");
	if (st->hd.fce_vid < 0)
		quit("open_hydro: no faces variable in %s\n", name);

	st->hd.fcex1_vid = ncvarid(st->hd.fid, "pt1_x");
	if (st->hd.fcex1_vid < 0)
		quit("open_hydro: no pt1_x variable in %s\n", name);
	st->hd.fcey1_vid = ncvarid(st->hd.fid, "pt1_y");
	if (st->hd.fcey1_vid < 0)
		quit("open_hydro: no pt1_y variable in %s\n", name);
	st->hd.fcex2_vid = ncvarid(st->hd.fid, "pt2_x");
	if (st->hd.fcex2_vid < 0)
		quit("open_hydro: no pt2_x variable in %s\n", name);
	st->hd.fcey2_vid = ncvarid(st->hd.fid, "pt2_y");
	if (st->hd.fcey2_vid < 0)
		quit("open_hydro: no pt2_y variable in %s\n", name);

	st->hd.e_vid = ncvarid(st->hd.fid, "transport");
	if (st->hd.e_vid < 0)
		quit("open_hydro: no transport variable in %s\n", name);

	/* Check variable types and dimensions for the transport data */
	if (verbose > 1)
		printf("open_hydro: get exchange data ids\n");

	ncvarinq(st->hd.fid, st->hd.e_vid, NULL, &daty, &ndims, dims, &natts);
	if (nctypelen(daty) != sizeof(FPTYPE))
		warn(
				"open_hydro: Type of exchange variable doesn't match model (ok if float in file and double FPTYPE)\n");
	if (ndims != 3 || dims[0] != st->hd.t_did || dims[1] != st->hd.b_did
			|| dims[2] != st->hd.z_did)
		quit("open_hydro: transport variable has incorrect dimensions\n");
	ncattget(st->hd.fid, st->hd.e_vid, "units", st->hd.e_units);

	if (verbose > 1)
		printf("open_hydro: get exchange units\n");

	sscanf(st->hd.e_units, "%s", stu);
	if (strstr("Sverdrup", stu) != NULL)
		st->unit_type = 0;
	else
		st->unit_type = 1;

	/* Reset netCDF error handling */
	ncopts = NC_VERBOSE | NC_FATAL;

	return;
}

/****************************************************************************************************
 Routine to close netcdf file containing raw hydrodynamic fluxes
 *****************************************************************************************************/
void close_hydro(Stuff *st) {
	/* Close file */
	if (st->hd.fid >= 0)
		ncclose(st->hd.fid);
	st->hd.fid = -1;

	return;
}

/***************************************************************************************************
 Routine to calculate number of data points, faces and levels to be loaded from netcdf file
 ****************************************************************************************************/
void get_netcdf_dimensions(Stuff *st) {
	int i, totndata, maxfaces, minlevel, maxlevel, totnvdata, maxboxes;
	long nt, fce, lvl, nvt, vlvl, nvb;

	/* Initialise counters */
	totndata = 0;
	maxfaces = 0;
	minlevel = MAXINT;
	maxlevel = 0;

	for (i = 0; i < st->hd.nfiles; i++) {
		if (verbose)
			printf("Reading raw exchange file to get dimensions %s\n",
					st->hd.fname[i]);

		/* Open the exchange file */
		open_hydro(st, st->hd.fname[i], i);

		/* Check for dimensions in file */
		ncdiminq(st->hd.fid, ncdimid(st->hd.fid, "time"), NULL, &nt);
		ncdiminq(st->hd.fid, ncdimid(st->hd.fid, "faces"), NULL, &fce);
		ncdiminq(st->hd.fid, ncdimid(st->hd.fid, "level"), NULL, &lvl);

		/* Check against counters */
		totndata += nt;
		if (fce > maxfaces)
			maxfaces = fce;
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

	printf("ndata: %d, nface: %d vs maxface: %d\n", st->ndata, st->nface,
			maxfaces);

	if (st->nface < maxfaces) {
		warn(
				"One of the raw hydro files contains more faces than the bgm file - reset nface\n");
		st->nface = maxfaces;
	}
	if (st->wcnz > minlevel) {
		warn(
				"One of the raw hydro files contains less water column layers (%d) than the parameter file specifies (%d)\n",
				minlevel, st->wcnz);
	}
	if (st->wcnz < maxlevel) {
		warn(
				"One of the raw hydro files contains less water column layers (%d) than the parameter file specifies (%d)\n",
				maxlevel, st->wcnz);
	}

	/* Get vertical exchange dimensions */
	totnvdata = 0;
	minlevel = MAXINT;
	maxlevel = 0;
	maxboxes = 0;
	for (i = 0; i < st->ts.ntsfiles; i++) {
		if (verbose)
			printf("Reading raw exchange file to get dimensions %s\n",
					st->ts.tsfname[i]);

		/* Open the exchange file */
		open_tsfile(st, st->ts.tsfname[i], i);

		/* Check for dimensions in file */
		ncdiminq(st->ts.tsfid, ncdimid(st->ts.tsfid, "time"), NULL, &nvt);
		ncdiminq(st->ts.tsfid, ncdimid(st->ts.tsfid, "level"), NULL, &vlvl);
		ncdiminq(st->ts.tsfid, ncdimid(st->ts.tsfid, "boxes"), NULL, &nvb);

		/* Check against counters */
		totnvdata += nvt;
		if (nvb > maxboxes)
			maxboxes = nvb;
		if (vlvl < minlevel)
			minlevel = vlvl;
		if (vlvl > maxlevel)
			maxlevel = vlvl;

		/* Close hydro file */
		close_tsfile(st);
	}

	/* Check aganist stored values */
	st->nvhdsteps = totnvdata;

	//st->nvhdsteps = 0;
	printf("set nhd: %d, nvhd: %d\n", st->nhdsteps, st->nvhdsteps);

	if (st->nbox < maxboxes) {
		warn(
				"One of the raw vert hydro files contains more boxes than the bgm file - reset nbox\n");
		st->nbox = maxboxes;
	}
	if (st->wcnz > minlevel) {
		warn(
				"One of the raw vert hydro files contains less water column layers (%d) than the parameter file specifies (%d)\n",
				minlevel, st->wcnz);
	}
	if (st->wcnz < maxlevel) {
		warn(
				"One of the raw vert hydro files contains less water column layers (%d) than the parameter file specifies (%d)\n",
				maxlevel, st->wcnz);
	}

	return;
}

/**
 *
 * \brief Read in the mappings between the input hydro file boxes and the bgm file.
 *
 * This allows the user to swap box numbers.
 *
 *
 *
 */
void load_bgm_mappings(Stuff *st) {

	FILE *fp;
	int hydroBoxNumber, bgmBoxNumber;
	int buflen = 500;
	char buf[500];

	/* Read and load lookup table */
	if (verbose)
		printf("Reading lookup table file %s\n", st->lookupIfname);

	/* Initialise look-up box array */
	boxLookup = (int *) i_alloc1d(st->nbox);

	/* Open the lookup table file */
	if ((fp = fopen(st->lookupIfname, "r")) == NULL){
		warn("load_bgm_mappings - Lookup table file: Can't open %s. Values loaded from the input file will be used.\n",
				st->lookupIfname);
		for(hydroBoxNumber = 0; hydroBoxNumber < st->nbox; hydroBoxNumber++)
			boxLookup[hydroBoxNumber] = hydroBoxNumber;

		return;
	}
	/* read in each line */
	while (fgets(buf, buflen, fp) != NULL) {

		//printf("buf = %s\n",buf);
		/* Ignore the first line */
		if (strstr(buf, "BGMBox") != NULL) {
			continue;
		}

		/* Read data line */
		if (sscanf(buf, "%d,%d", &hydroBoxNumber, &bgmBoxNumber) != 2)
			quit("not enough entries for %s - with %d, %d\n", buf,
					hydroBoxNumber, bgmBoxNumber);

		if (verbose > 3)
			printf("readin %s, hydroBoxNumber-%d, bgmBoxNumber-%d\n", buf,
					hydroBoxNumber, bgmBoxNumber);

		if (hydroBoxNumber >= st->nbox)
			quit(
					"There are only %d boxes in the model - your box_conversion input has a hydro box number of %d?\n",
					st->nbox, hydroBoxNumber);

		if (bgmBoxNumber >= st->nbox)
			quit(
					"There are only %d boxes in the model - your box_conversion input has a bgm box number of %d?\n",
					st->nbox, bgmBoxNumber);

		/* All is ok - store the numbers */
		boxLookup[hydroBoxNumber] = bgmBoxNumber;
	}

	/* Close the conversion file */
	fclose(fp);
}
