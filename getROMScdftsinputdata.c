/**********************************************************************

    getROMScdftsinputdata.c

	by:	Beth Fulton
	date:	11/3/2009

	comments: This file is the code that reads in NOAA netcdf file version of
			raw hydrodynamics flows, and temperature and salinity profiles

    assumptions: assumes raw input is a series of data entries in netcdf file
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
Routine to get the temperature, salinity and vertical flux data out of the raw data file - Al
Herman ROMS based output. Involves reading data from netcdf file to buffer
****************************************************************************************************/
void get_ROMS_cdf_temp_and_salt(Stuff *st)
{
	int i, k, lastk, entry_i, nz, ij, got_ve;
	double time_of_year, totnum, thisdate, this_tofy, dtmult;
    FPTYPE ***valexchange;
    FPTYPE ***valtemp;
	FPTYPE ***valsalt;
    FPTYPE *valtime;
    int *valbox;
    long start[3];
    long count[3];
	long startt[1];
	long countt[1];
	long nt, lvl, nb;

	st->tsflagflux = 0;
    ncopts = NC_VERBOSE | NC_FATAL;

	if(verbose)
		printf("Read in raw temperature and salinity profiles \n");

	/* Loop through the data files */
	lastk = 0;
	for(i=0; i<st->ts.ntsfiles; i++) {
		if(verbose)
			printf("Reading raw exchange file %s\n",st->ts.tsfname[i]);

		/* Open the exchange file */
	    open_ROMStsfile(st,st->ts.tsfname[i],i);

	    /* Check for dimensions in file */
		ncdiminq(st->ts.tsfid,ncdimid(st->ts.tsfid,"T12"),NULL,&nt);
		ncdiminq(st->ts.tsfid,ncdimid(st->ts.tsfid,"ZBIN"),NULL,&lvl);
		ncdiminq(st->ts.tsfid,ncdimid(st->ts.tsfid,"BOXID"),NULL,&nb);

		/* Allocate temporary storage for one tracer */
		switch( sizeof(FPTYPE) ) {
		case sizeof(float):
			valexchange = (FPTYPE ***)f_alloc3d(lvl,nb,nt);
			valtemp = (FPTYPE ***)f_alloc3d(lvl,nb,nt);
			valsalt = (FPTYPE ***)f_alloc3d(lvl,nb,nt);
			valtime = (FPTYPE *)f_alloc1d(nt);
			break;
		case sizeof(double):
			valexchange = (FPTYPE ***)d_alloc3d(lvl,nb,nt);
			valtemp = (FPTYPE ***)d_alloc3d(lvl,nb,nt);
			valsalt = (FPTYPE ***)d_alloc3d(lvl,nb,nt);
			valtime = (FPTYPE *)d_alloc1d(nt);
			break;
		default:
			quit("get_cdf_hydro: Unknown size for FPTYPE\n");
		}
		valbox = (int *)i_alloc1d(nb);

		/* Set indices for reading tracers */
		start[0] = 0;
		start[1] = 0;
		start[2] = 0;
		count[0] = nt;
		count[1] = nb;
		count[2] = lvl;

		startt[0] = 0;
		countt[0] = nt;

		/* Read data - exchanges, dates and box-ids associated with the faces */
		ncvarget(st->ts.tsfid,st->ts.t_vid,startt,countt,valtime);
		//ncvarget(st->ts.tsfid,st->ts.b_vid,start,count,valbox); Not read in as assumed to match that in horizontal read-in
		ncvarget(st->ts.tsfid,st->ts.temp_vid,start,count,valtemp[0][0]);
		ncvarget(st->ts.tsfid,st->ts.salt_vid,start,count,valsalt[0][0]);
		if( st->ts.ve_vid < 0 ){
			/* Do nothing */
			got_ve = 0;
		} else {
			/* Load the vertical fluxes */
			ncvarget(st->ts.tsfid,st->ts.ve_vid,start,count,valexchange[0][0]);
			got_ve = 1;
		}

		/* Check timestamps */
		dtmult = 3600; // As Al Herman data in hours not seconds
		for(k=0; k<nt; k++){
			thisdate = valtime[k] * dtmult - st->time_reset;

			if(thisdate != exchangedate[k+lastk][time_id]){
				quit("Vertical and horizontal dates don't match (hdate: %e, vdate: %e) - k+lastk: %d (k: %d), valtime: %e, dt: %d - current cdf file %s\n",
					exchangedate[k+lastk][time_id], thisdate, k+lastk , k, valtime[k], st->dt, st->ts.tsfname[i]);
			}

			/* Multiple by 86400 / dt as may be tidal entries so must expand
			list to fit these extra entries in */
			time_of_year = valtime[k] * (dtmult / st->dt);
			totnum = 365.0 * (86400 / st->dt);

			if( time_of_year/totnum < 1.0 )
				this_tofy = floor(time_of_year/(86400 / st->dt) + 0.5);
			else
				this_tofy = floor((time_of_year/totnum - floor(time_of_year/totnum))*totnum + 0.5);

			if(this_tofy != exchangedate[k+lastk][TofY_id]){
				quit("Vertical and horizontal dates don't match (hdate %e, vdate: %d)\n",
					exchangedate[k+lastk][TofY_id], this_tofy);
			}

			/* Transfer results */
			tempsaltdate[k+lastk][TofY_id] = exchangedate[k+lastk][TofY_id];
			tempsaltdate[k+lastk][time_id] = exchangedate[k+lastk][time_id];
		}

		/* Box information assumed to be dealt with in base hydrodynamics cdf readin */

		/* Transfer exchanges, temperature and salinity */
		entry_i = 0;
		for(k=0; k<nt; k++){
			for(nz=0; nz<lvl; nz++){
				for(ij=0; ij<nb; ij++){
					if(got_ve)
						rawvertexchanges[k+lastk][nz][ij][entry_i] = valexchange[k][ij][nz];
					else
						rawvertexchanges[k+lastk][nz][ij][entry_i] = 0;
					rawtempsalt[k+lastk][nz][ij][temp_id] = valtemp[k][ij][nz];
					rawtempsalt[k+lastk][nz][ij][salt_id] = valsalt[k][ij][nz];

					if(verbose > 1)
						printf("Storing rawtemp[%d][%d][%d] = %e\n", k+lastk, nz, ij, rawtempsalt[k+lastk][nz][ij][temp_id]);
				}
			}
		}

		/* Set lastk */
		lastk += nt;

		/* Free temporary storage */
		switch( sizeof(FPTYPE) ) {
		case sizeof(float):
			f_free3d((float ***)valexchange);
			f_free3d((float ***)valtemp);
			f_free3d((float ***)valsalt);
			f_free1d((float *)valtime);
			break;
		case sizeof(double):
			d_free3d((double ***)valexchange);
			d_free3d((double ***)valtemp);
			d_free3d((double ***)valsalt);
			d_free1d((double *)valtime);
			break;
		default:
			quit("get_cdf_hydro: Unknown size for FPTYPE\n");
		}
		i_free1d((int *)valbox);

		/* Close hydro file */
		close_tsfile(st);
	}

	if(verbose)
		printf("Finished reading raw temperature and salinity files\n");

	return;
}

/****************************************************************************************************
Routine to open netcdf file containing raw temperature, salinity and vertical flux data
*****************************************************************************************************/
void open_ROMStsfile(Stuff *st, char *name, int filenum)
{
    int ndims;
    int nvars;
    int natts;
    int recdim;
    int len;
    long n;
    nc_type daty;
    char hdu[STSLEN];
    char stu[STSLEN];

    /* Set netCDF library error handling */
    ncopts = NC_VERBOSE;

    if(verbose > 1)
   		printf("open_ROMStsfile: opening the file %s\n", name);

    /* Open the file */
    if( (st->ts.tsfid=ncopen(name,NC_NOWRITE)) < 0 )
	   quit("open_tsfile: Can't open temperature and salinity model input data file %s\n",name);

    /* Inquire about this file */
    ncopts = NC_VERBOSE | NC_FATAL;
    ncinquire(st->ts.tsfid,&ndims,&nvars,&natts,&recdim);
    if( ndims < 3 )
        quit("open_tsfile: not enough dimensions in %s\n",name);

    /* Check dimensions are as expected */
    if( (st->ts.t_did = ncdimid(st->ts.tsfid,"T12"))  == -1 )
        quit("open_tsfile: no t dimension in %s\n",name);
    if( st->ts.t_did != recdim )
        quit("open_tsfile: t dimension not unlimited in %s\n",name);
    if( (st->ts.b_did = ncdimid(st->ts.tsfid,"BOXID"))  == -1 )
        quit("open_tsfile: no boxes dimension in %s\n",name);
    if( (st->ts.z_did = ncdimid(st->ts.tsfid,"ZBIN"))  == -1 )
        quit("open_tsfile: no levels dimension in %s\n",name);

    /* Get dimension sizes. Note check against geometry now happens in get_cdf_hydro() */
    ncdiminq(st->ts.tsfid,st->ts.b_did,NULL,&n);
	ncdiminq(st->ts.tsfid,st->ts.z_did,NULL,&n);

    /* Check that time units and steps match this model */
    st->ts.t_vid = ncvarid(st->ts.tsfid,"T12");
    if( st->ts.t_vid < 0 )
		quit("open_tsfile: no t variable in %s\n",name);
    memset(st->ts.t_units,0,STSLEN);
    ncattget(st->ts.tsfid,st->ts.t_vid,"units",st->ts.t_units);
    sscanf(st->ts.t_units,"%s",hdu);
	if(!filenum)
		sprintf(st->t_units,"%s",st->ts.t_units);
    sscanf(st->t_units,"%s",stu);
    if( strcmp(hdu,stu) != 0 ){
        quit("open_tsfile: Time units (%s) in %s don't match model time units (%s) loaded from %s\n",
			hdu, name, stu, st->ts.tsfname[0]);
	}

	/* Check for timestep size */
    ncattinq(st->ts.tsfid,st->ts.t_vid,"dt",&daty,&len);
    if( nctypelen(daty) != sizeof(st->ts.tsdt) )
        quit("open_tsfile: dt attribute wrong type\n");
    ncattget(st->ts.tsfid,st->ts.t_vid,"dt",&st->ts.tsdt);
    if( st->ts.tsdt <= 0.0 )
	quit("open_tsfile: dt attribute must have positive value\n");

    /* Find out how many time steps are in the file */
    ncdiminq(st->ts.tsfid,st->ts.t_did,NULL,&st->ts.tsnstep);

    /* Get other variable ids */
    st->ts.b_vid = ncvarid(st->ts.tsfid,"BOXID");
    if( st->ts.b_vid < 0 )
		quit("open_tsfile: no boxes variable in %s\n",name);

	st->ts.temp_vid = ncvarid(st->ts.tsfid,"TEMPMEAN");
    if( st->ts.temp_vid < 0 )
		quit("open_tsfile: no temperature variable in %s\n",name);
    st->ts.salt_vid = ncvarid(st->ts.tsfid,"SALTMEAN");
    if( st->ts.salt_vid < 0 )
		quit("open_tsfile: no salinity variable in %s\n",name);
  	st->ts.ve_vid = ncvarid(st->ts.tsfid,"verticalflux");
    if( st->ts.ve_vid < 0 )
		warn("open_tsfile: no verticalflux variable in %s\n",name);

    /* Check variable types and dimensions for the transport data *
    ncvarinq(st->ts.tsfid,st->ts.ve_vid,NULL,&daty,&ndims,dims,&natts);
    if( nctypelen(daty) != sizeof(FPTYPE) )
        warn("open_tsfile: Type of verticalflux variable doesn't match model (OK if float in file and double FPTYPE)\n");
    if( ndims != 3 || dims[0] != st->ts.t_did || dims[1] != st->ts.b_did || dims[2] != st->ts.z_did)
        quit("open_tsfile: verticalflux variable has incorrect dimensions\n");
    ncattget(st->ts.tsfid,st->ts.ve_vid,"units",st->ts.ve_units);
	sscanf(st->ts.ve_units,"%s",stu);
	if(strstr("Sverdrup",stu) != NULL ){
		if(st->unit_type)
			quit("Flux data mismatch - horizontal fluxes in m3 and vertical in Sverdrup\n");
	} else {
		if(!st->unit_type)
			quit("Flux data mismatch - horizontal fluxes in Sverdrup and vertical in m3\n");
	}
*/
    /* Reset netCDF error handling */
    ncopts = NC_VERBOSE | NC_FATAL;

	return;
}

