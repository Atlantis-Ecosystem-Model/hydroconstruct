/**********************************************************************

    getROMScdfinputdata.c

	by:	Beth Fulton
	date:	2/3/2009

	comments: This file is the code that reads in netcdf file version of
			raw hydrodynamics flows

    assumptions: assumes raw input is a series of data entries in netcdf file


    Changes:


		27-04-2009 Bec Gorton
		Changed the map_ROMS_netcdf_faces function to set the facemap values. These were just -1 so the calcexchange code
		was not running properly as it was never finding a destination face.

		02-06-2009 Bec Gorton
		Changed the error message that is printed when the input file is not found to
		reflect the name of the calling function.
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
void get_ROMS_cdf_hydro(Stuff *st)
{
	int i, k, lastk, entry_i, nz, ij, final_fce, b;
	double time_of_year, totnum, this_tofy, mintime, dtmult;
    FPTYPE ****valexchange;
    FPTYPE *valtime;
    int *valboxs;
    int *valfaceid;
    long start[4];
    long count[4];
	long startt[1];
	long countt[1];
	long startfce[1];
	long countfce[1];
	long startbx[1];
	long countbx[1];

	long nt, fce, lvl, bx;

    ncopts = NC_VERBOSE | NC_FATAL;

	if(verbose)
		printf("Loading raw data from netcdf files\n");

	/* Give code warning */
	if(st->hd.nfiles > 1)
		warn("%d raw hydro files listed, code may not cope properly with more than one - in particular the face mapping may be wrong\n",st->hd.nfiles);

	/* Loop through the exchange files */
	lastk = 0;
	for(i=0; i<st->hd.nfiles; i++) {
		if(verbose)
			printf("Reading raw exchange file %s\n",st->hd.fname[i]);

		/* Open the exchange file */
	    open_ROMS_hydro(st,st->hd.fname[i],i);

		if(verbose > 1)
			printf("Checking for dimensions\n");

	    /* Check for dimensions in file */
		ncdiminq(st->hd.fid,ncdimid(st->hd.fid,"T12"),NULL,&nt);
		ncdiminq(st->hd.fid,ncdimid(st->hd.fid,"SECID"),NULL,&fce);
		ncdiminq(st->hd.fid,ncdimid(st->hd.fid,"BOXID"),NULL,&bx);
		ncdiminq(st->hd.fid,ncdimid(st->hd.fid,"ZBIN"),NULL,&lvl);

		if(fce > st->nbgmface){
			warn("More faces in netcdf file %s (nface = %d) than in bgm file %s (nface = %d) - reconcile and retry\n",
				st->hd.fname[i], fce, st->geomIfname, st->nface);
		} else if (fce < st->nbgmface){
			warn("Less faces in netcdf file %s (nface = %d) than in bgm file %s (nface = %d)\n",
				st->hd.fname[i], fce, st->geomIfname, st->nface);
		}

	    /* Allocate temporary storage for one tracer */
		switch( sizeof(FPTYPE) ) {
		case sizeof(float):
			valexchange = (FPTYPE ****)f_alloc4d(fce,lvl,bx,nt);
			valtime = (FPTYPE *)f_alloc1d(nt);
			break;
		case sizeof(double):
			valexchange = (FPTYPE ****)d_alloc4d(fce,lvl,bx,nt);
			valtime = (FPTYPE *)d_alloc1d(nt);
			break;
		default:
			quit("get_cdf_hydro: Unknown size for FPTYPE\n");
		}
		valboxs = (int *)i_alloc1d(bx);
		valfaceid = (int *)i_alloc1d(fce);

		/* Set indices for reading tracers */
		start[0] = 0;
		start[1] = 0;
		start[2] = 0;
		start[3] = 0;
		count[0] = nt;
		count[1] = bx;
		count[2] = lvl;
		count[3] = fce;

		startbx[0] = 0;
		countbx[0] = bx;

		startfce[0] = 0;
		countfce[0] = fce;

		startt[0] = 0;
		countt[0] = nt;

		/* Read data - exchanges, dates and box-ids associated with the faces */
		ncvarget(st->hd.fid,st->hd.t_vid,startt,countt,valtime);
		ncvarget(st->hd.fid,st->hd.srcb_vid,startbx,countbx,valboxs);
		ncvarget(st->hd.fid,st->hd.fce_vid,startfce,countfce,valfaceid);
		ncvarget(st->hd.fid,st->hd.e_vid,start,count,valexchange[0][0][0]);

		/** Transfer to rawdata arrays **/
		/* Transfer boxes and faces */
		for(b=0; b<bx; b++){
			for(k=0; k<fce; k++){
				tempfaceinfo[valboxs[b]][k][fceid_id] = valfaceid[k];
			}
		}

		/* Transfer timestamps */
		dtmult = 3600; // As Al Herman data in hours not seconds
		if(st->reset_time){
			mintime = MAXDOUBLE;
			for(k=0; k<nt; k++){
				if(valtime[k] < mintime)
					mintime = valtime[k] * dtmult;
			}
			st->time_reset = mintime;
		}

		for(k=0; k<nt; k++){
			exchangedate[k+lastk][time_id] = valtime[k] * dtmult - st->time_reset;

			/* Multiple by dtmult / dt so in days roughly */
			time_of_year = valtime[k] * dtmult / st->dt;
			totnum = 365.0 * (86400 / st->dt);

			if( time_of_year/totnum < 1.0 )
				this_tofy = floor(time_of_year/(86400 / st->dt) + 0.5);
			else
				this_tofy = floor((time_of_year/totnum - floor(time_of_year/totnum))*totnum + 0.5);

			exchangedate[k+lastk][TofY_id] = this_tofy;

			if(st->verbose_pad)
				printf("k: %d, time: %e, tofy: %e valtime: %e\n", k, exchangedate[k+lastk][time_id],  exchangedate[k+lastk][TofY_id] , valtime[k]);
		}

		/* Map netcdf box assignment to bgm assignment */
		map_ROMS_netcdf_faces(st);

		/* Transfer exchanges */
		entry_i = 0;
		for(k=0; k<nt; k++){
			for(b=0; b<bx; b++){
				for(nz=0; nz<lvl; nz++){
					for(ij=0; ij<fce; ij++){
						final_fce = tempfaceinfo[b][ij][rbox_id]; // Actually this is a face ID (overloading rbox_id as just a placeholder index
						rawexchanges[k+lastk][nz][final_fce][entry_i] = valexchange[k][b][nz][ij];

						if(!k && st->verbose_exchange){
							printf("rawexchange[%d][%d][%d][%d] = %e, (val[%d][%d][%d][%d])\n",
								(k+lastk), nz, final_fce, entry_i, rawexchanges[k+lastk][nz][final_fce][entry_i], k, b, nz, ij);
						}
					}
				}
			}
		}
		/* Set lastk */
		lastk += nt;

		/* Free temporary storage */
		switch( sizeof(FPTYPE) ) {
		case sizeof(float):
			f_free4d((float ****)valexchange);
			f_free1d((float *)valtime);
			break;
		case sizeof(double):
			d_free4d((double ****)valexchange);
			d_free1d((double *)valtime);
			break;
		default:
			quit("get_cdf_hydro: Unknown size for FPTYPE\n");
		}
		i_free1d((int *)valboxs);
		i_free1d((int *)valfaceid);

		/* Close hydro file */
		close_hydro(st);
	}

	if(verbose)
		printf("Finished reading raw exchange files\n");

	return;
}

/****************************************************************************************************
Routine to map box and face ids from netcdf file to matching ones in bgm file
*****************************************************************************************************/
void map_ROMS_netcdf_faces(Stuff *st)
{
	int b, k;

	/* Read and load lookup table */
	if(verbose)
		printf("Setting up ROMS-BGM face matching\n");

	/* Read table file */
	for(b=0; b<st->nbox; b++){
		for(k=0; k<st->boxnconn[b]; k++){
			tempfaceinfo[b][k][rbox_id] = st->iface[b][k]; // Actually this is a face ID (overloading rbox_id as just a placeholder index

			if(verbose > 3)
				printf("stored box-%d, face-%d, bgmface-%d\n", b, k, st->iface[b][k]);
		}
	}

	for(b=0; b<st->nface; b++){
		facemap[b] = b;
	}
	return;
}

/****************************************************************************************************
Routine to open netcdf file containing raw hydrodynamic fluxes
*****************************************************************************************************/
void open_ROMS_hydro(Stuff *st, char *name, int filenum)
{
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

	if(verbose > 1)
		printf("open_hydro: opening the file %s\n", name);

    /* Open the file */
    if( (st->hd.fid=ncopen(name,NC_NOWRITE)) < 0 )
	   quit("open_ROMS_hydro: Can't open hydrodynamic model input data file %s\n",name);

    /* Inquire about this file */
    ncopts = NC_VERBOSE | NC_FATAL;
    ncinquire(st->hd.fid,&ndims,&nvars,&natts,&recdim);
    if( ndims < 3 ){
        quit("open_hydro: not enough dimensions in %s\n",name);
    }

	if(verbose > 1)
		printf("open_hydro: check dimensions\n");

    /* Check dimensions are as expected */
	printf("check for time in ROMS cdf\n");

    if( (st->hd.t_did = ncdimid(st->hd.fid,"T12"))  == -1 )
        quit("open_hydro: no t dimension in %s\n",name);
    if( st->hd.t_did != recdim )
        quit("open_hydro: t dimension not unlimited in %s\n",name);

    printf("check for boxes\n");

	if( (st->hd.b_did = ncdimid(st->hd.fid,"BOXID"))  == -1 ){
			quit("open_hydro: no boxes dimension in %s\n",name);
	}

    printf("check for faces\n");

	if( (st->hd.f_did = ncdimid(st->hd.fid,"SECID"))  == -1 ){
			quit("open_hydro: no faces dimension in %s\n",name);
	}

	printf("check for level\n");

    if( (st->hd.z_did = ncdimid(st->hd.fid,"ZBIN"))  == -1 )
        quit("open_hydro: no levels dimension in %s\n",name);

    /* Get dimension sizes. Note check against geometry now happens in get_cdf_hydro() */

	// TO DO CHECK THAT 'check vs geometry' IS SENSIBLE

	if(verbose > 1)
		printf("open_hydro: get dimnensionn\n");

    ncdiminq(st->hd.fid,st->hd.b_did,NULL,&n);
	ncdiminq(st->hd.fid,st->hd.z_did,NULL,&n);

    /* Check that time units and steps match this model */
	if(verbose > 1)
		printf("open_hydro: get units\n");

    st->hd.t_vid = ncvarid(st->hd.fid,"T12");
    if( st->hd.t_vid < 0 )
		quit("open_hydro: no t variable in %s\n",name);
    memset(st->hd.t_units,0,STSLEN);
    ncattget(st->hd.fid,st->hd.t_vid,"units",st->hd.t_units);
    sscanf(st->hd.t_units,"%s",hdu);
	if(!filenum)
		sprintf(st->t_units,"%s",st->hd.t_units);
    sscanf(st->t_units,"%s",stu);
    if( strcmp(hdu,stu) != 0 ){
        quit("open_hydro: Time units (%s) in %s don't match model time units (%s) loaded from %s\n",
			hdu, name, stu, st->hd.fname[0]);
	}

	st->hd.dt = 43200.0;

    /* Find out how many time steps are in the file */
    ncdiminq(st->hd.fid,st->hd.t_did,NULL,&st->hd.nstep);

    /* Get other variable ids */
    st->hd.srcb_vid = ncvarid(st->hd.fid,"BOXID");
    if( st->hd.srcb_vid < 0 )
		quit("open_hydro: no boxid variable in %s\n",name);
 	st->hd.e_vid = ncvarid(st->hd.fid,"TRANS");
    if( st->hd.e_vid < 0 )
		quit("open_hydro: no transport variable in %s\n",name);

    /* Check variable types and dimensions for the transport data */
	if(verbose > 1)
		printf("open_hydro: get exchange data ids\n");

    ncvarinq(st->hd.fid,st->hd.e_vid,NULL,&daty,&ndims,dims,&natts);
    if( nctypelen(daty) != sizeof(FPTYPE) )
        warn("open_hydro: Type of exchange variable doesn't match model (ok if float in file and double FPTYPE)\n");
    if( ndims != 4 || dims[0] != st->hd.t_did || dims[1] != st->hd.b_did || dims[2] != st->hd.z_did || dims[3] != st->hd.f_did)
        quit("open_hydro: transport variable has incorrect dimensions\n");

    ncattget(st->hd.fid,st->hd.e_vid,"units",st->hd.e_units);
	if(verbose > 1)
		printf("open_hydro: get exchange units\n");

	sscanf(st->hd.e_units,"%s",stu);
	if(strstr("Sverdrup",stu) != NULL )
		st->unit_type = 0;
	else
		st->unit_type = 1;

    /* Reset netCDF error handling */
    ncopts = NC_VERBOSE | NC_FATAL;

	return;
}

/***************************************************************************************************
Routine to calculate number of data points, faces and levels to be loaded from netcdf file
****************************************************************************************************/
void get_ROMS_netcdf_dimensions(Stuff *st)
{
	int i, totndata, maxfaces, minlevel, maxlevel, totnvdata, maxboxes;
	long nt, fce, lvl, bx, nvt, vlvl, nvb;

	/* Initialise counters */
	totndata = 0;
	maxfaces = 0;
	minlevel = MAXINT;
	maxlevel = 0;
	maxboxes = 0;

	for(i=0; i<st->hd.nfiles; i++) {
		if(verbose)
			printf("Reading raw exchange file to get dimensions %s\n",st->hd.fname[i]);

		/* Open the exchange file */
	    open_ROMS_hydro(st,st->hd.fname[i],i);

	    /* Check for dimensions in file */
		ncdiminq(st->hd.fid,ncdimid(st->hd.fid,"T12"),NULL,&nt);
		ncdiminq(st->hd.fid,ncdimid(st->hd.fid,"SECID"),NULL,&fce);
		ncdiminq(st->hd.fid,ncdimid(st->hd.fid,"ZBIN"),NULL,&lvl);
		ncdiminq(st->hd.fid,ncdimid(st->hd.fid,"BOXID"),NULL,&bx);

		/* Check against counters */
		totndata += nt;
		if(fce > maxfaces)
			maxfaces = fce;
		if(bx > maxboxes)
			maxboxes = bx;
		if(lvl < minlevel)
			minlevel = lvl;
		if(lvl > maxlevel)
			maxlevel = lvl;

		/* Close hydro file */
		close_hydro(st);

	}

	/* Check aganist stored values */
	st->ndata = totndata;
	st->nhdsteps = st->ndata;

	printf("ndata: %d, nface: %d vs maxface: %d, nbox: %d vs maxbox: %d\n", st->ndata, st->nface, maxfaces, st->nbox, maxboxes);

	if(st->nface < maxfaces){
		warn("One of the raw hydro files contains more faces than the bgm file - reset nface\n");
		st->nface = maxfaces;
	}
	if(st->nbox < maxboxes){
		warn("One of the raw hydro files contains more boxes than the bgm file - reset nbox\n");
		st->nbox = maxboxes;
	}
	if(st->wcnz > minlevel){
		warn("One of the raw hydro files contains less water column layers (%d) than the parameter file specifies (%d)\n", minlevel, st->wcnz);
	}
	if(st->wcnz < maxlevel){
		warn("One of the raw hydro files contains less water column layers (%d) than the parameter file specifies (%d)\n", maxlevel, st->wcnz);
	}

	/* Get vertical exchange dimensions */
	totnvdata = 0;
	minlevel = MAXINT;
	maxlevel = 0;
	maxboxes = 0;
	for(i=0; i<st->ts.ntsfiles; i++) {
		if(verbose)
			printf("Reading raw exchange file to get dimensions %s\n",st->ts.tsfname[i]);

		/* Open the exchange file */
	    open_tsfile(st,st->ts.tsfname[i],i);

	    /* Check for dimensions in file */
		ncdiminq(st->ts.tsfid,ncdimid(st->ts.tsfid,"T12"),NULL,&nvt);
		ncdiminq(st->ts.tsfid,ncdimid(st->ts.tsfid,"ZBIN"),NULL,&vlvl);
		ncdiminq(st->ts.tsfid,ncdimid(st->ts.tsfid,"BOXID"),NULL,&nvb);

		/* Check against counters */
		totnvdata += nvt;
		if(nvb > maxboxes)
			maxboxes = nvb;
		if(vlvl < minlevel)
			minlevel = vlvl;
		if(vlvl > maxlevel)
			maxlevel = vlvl;

		/* Close hydro file */
		close_tsfile(st);
	}

	/* Check aganist stored values */
	st->nvhdsteps = totnvdata;

	printf("set nhd: %d, nvhd: %d\n", st->nhdsteps, st->nvhdsteps);

	if(st->nbox < maxboxes){
		warn("One of the raw vert hydro files contains more boxes than the bgm file - reset nbox\n");
		st->nbox = maxboxes;
	}
	if(st->wcnz > minlevel){
		warn("One of the raw vert hydro files contains less water column layers (%d) than the parameter file specifies (%d)\n", minlevel, st->wcnz);
	}
	if(st->wcnz < maxlevel){
		warn("One of the raw vert hydro files contains less water column layers (%d) than the parameter file specifies (%d)\n", maxlevel, st->wcnz);
	}

	return;
}
