/**********************************************************************

    File containing code for writing output file of exchanges
	by:	Beth Fulton
	date:	21/10/2003

	comments: Writes cdf file that is later packed in to the hydodynamics exchange file.
			  These routines write text version of netcdf file so can just pack it
			  and then it is complete

    assumptions:

    Revisions:

			18-03-2009 Bec Gorton
			Changed the exchange output to allow for more than 5 output flow files.
			Also changes the code that creates the output to create a well formed cdf file
			so that no additional changes are required before calling ncgen.

			09-07-2009 Bec Gorton
			Fixed the writeTimeSteps function to print the start dates correctly when more than one
			flow file is created.

			14-09-2009 Bec Gorton
			Change the temp/salt output functions to a single function instead of
			two seperate functions.

*************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sjwlib.h>
#include "mainexchange.h"

FILE * initExchangeOutFile(Stuff *st, int index);
FILE * initTempOutFile(Stuff *st);
FILE * initSaltOutFile(Stuff *st);
FILE * initpHOutFile(Stuff *st);
FILE * initSWROutFile(Stuff *st);


void writeTempSaltOut(FILE *fid, Stuff *st, char *name, int index);


/**
 *
 * \brief Print out the timesteps to the output file.
 *
 */
void writeTimeSteps(Stuff *st, FILE *fid, int countStart, int countend)
{
	int lasti, returnnow;
	int i;
	int flag = 0;

	/* Write time stamps */
	lasti = 0;
	returnnow = 0;
	lasti = (countStart ) * st->dt;

	for(i=countStart; i<countend; i++){

		if(flag == 1)
			fprintf(fid,", %d",lasti);
		else
			fprintf(fid,"%d",lasti);
		flag = 1;

		returnnow++;
		if(returnnow == 7){
			fprintf(fid,"\n  ");
			returnnow = 0;
		}
		lasti += (int)(st->dt);
	}
}
void WriteNewOutput(Stuff *st)
{

	/* Write out exchanges */
	write_exchange_out(st);

	/* Write out temperature and salinity profiles - if loaded in the first place */
	if(st->input_style != col_face_type){
		write_temp_salt_out(st);
	}

	return;
}

/*************************************************************************************************
These routines write text version of exchange (hydrodynamic) netcdf file
**************************************************************************************************/
void write_exchange_out(Stuff *st)
{
	int i;
    FILE **outfp = NULL;

    outfp = (FILE **)malloc(sizeof(FILE *) * st->numoutfile);

    for(i = 0; i < st->numoutfile; i++){
    	outfp[i] = initExchangeOutFile(st, i);
    }
    for(i = 0; i < st->numoutfile; i++){
    	writeExchangeOut(outfp[i],st,i + 1);
	}

    free(outfp);
	return;

}

/*************************************************************************************************
File creation routines
**************************************************************************************************/
FILE * initExchangeOutFile(Stuff *st, int index)
/* Routine to initialise output file */
{
    FILE *fid;
    char fname[2*STSLEN];
	//int tout = (int)(((st->tstop - st->tstart) * (86400 / st->dt)) / st->numoutfile) - 1;
	int tout = (int)(st->num_entries / st->numoutfile);
	// -1 as will stop in output file writing because of < in for statement

	if(st->numoutfile == 1){
		sprintf(fname, "%s", st->transportOfname);
	}
	else{
		sprintf(fname, "%s%d", st->transportOfname, index + 1);
	}
    if(verbose)
		printf("Create hyrdo output file\n");

	/* Create file */
    if( (fid=fopen(fname,"w")) == NULL )
		quit("initOutFile: Can't open %s\n",fname);

    /* File header info */
    fprintf(fid,"netcdf exchange_data { \ndimensions:\n	t = UNLIMITED ; // (%d currently)\n",tout + 1);
	fprintf(fid,"	b = %d ;\n	z = %d ;\n	dest = %d ;\n",st->nbox, st->wcnz, st->fndest);
	fprintf(fid,"variables:\n");
	fprintf(fid, "	double t(t) ;\n		t:units = \"seconds since %d-01-01 00:00:00 +10\" ;\n		t:dt = %d. ;", st->REFyear, (int)(st->dt));
	fprintf(fid,"\n	double exchange(t, b, z, dest) ;\n		exchange:_FillValue = 0. ; \n\t\texchange:units = \"m^3\"; \n\t\texchange:long_name= \"Change in volume in this time step\";");
	fprintf(fid,"\n	int dest_b(t, b, z, dest) ;\n		dest_b:_FillValue = -1 ;");
	fprintf(fid,"\n	int dest_k(t, b, z, dest) ;\n		dest_k:_FillValue = -1 ;");
	fprintf(fid,"\n// global attributes:\n		:title = \"trivial\";\n		:geometry = \"%s\" ;",st->geomIfname);
	fprintf(fid,"\n		:parameters = \"\" ;\ndata:\n\n t = ");


    /* Return file pointer */
    return(fid);
}

/*************************************************************************************************
Routines to write out hydro files
**************************************************************************************************/
void writeExchangeOut(FILE *fid, Stuff *st, int numid)
{
    int dayj, i, j, k, this_id, counteach, endeach, countend,
		indexed_countend, indexed_start, tstart = 0;
	FILE *logfp;
    char logfname[STSLEN];
    int flag = 0;

	counteach = (int)((st->tstop * (86400 / st->dt)) / st->numoutfile);
	endeach = (int)(((st->tstop - st->tstart) * (86400 / st->dt)) / st->numoutfile);

	tstart = st->tstart + (numid - 1) * counteach;
	countend = numid * counteach;
	indexed_countend = numid * endeach;
	indexed_start = (numid - 1) * endeach;

	printf("st->tstop = %d, st->tstart= %d\n", st->tstop, st->tstart);

	printf("st->numoutfile = %d, 86400 / st->dt= %e\n", st->numoutfile, 86400 / st->dt);
	printf("indexed_countend = %d, counteach= %d, endeach = %d\n",indexed_countend, counteach, endeach);


	/* Check the transport will be stable over time */
	checkTransport(st, indexed_countend, indexed_start, numid);
	sprintf(logfname,"%s%d","hydrolog.txt",numid);

	/* Create file */
    if( (logfp=fopen(logfname,"w")) == NULL )
		quit("Can't open log file%s\n",logfname);

	if(verbose){
		printf("Write to hydro output file %s%d\n", st->transportOfname, numid);
	}
	warn("Write to hydro output file %s%d\n", st->transportOfname, numid);

	/* Write time stamps */
	writeTimeSteps(st, fid, tstart, countend);

	printf("indexed_start = %d, indexed_countend= %d\n", indexed_start, indexed_countend);

	/* Write time stamps */
	for(i=indexed_start; i<indexed_countend; i++){
		if(st->pad_time)
			this_id = paddeddates[exchange_id][i];
		else
			this_id = i;

		fprintf(logfp,"i: %d, id: %d, date %e, time: %e\n",i, this_id, exchangedate[this_id][TofY_id], exchangedate[this_id][time_id]);
	}
	/* Write exchanges */
	flag = 0;
	fprintf(fid,";\n\n exchange =\n  ");
	for(dayj=indexed_start; dayj<indexed_countend; dayj++){

		if(st->pad_time)
			this_id = paddeddates[exchange_id][dayj];
		else
			this_id = dayj;

		for(i=0; i<st->nbox; i++){
			for(j=0; j<st->wcnz; j++){
				if(verbose > 3)
					printf("boxndest[box%d]: %d, ndest: %d\n", i, boxndest[i], st->fndest);

				for(k=0; k<boxndest[i]; k++){
					if(flag == 1)
						fprintf(fid,", %.12f", finalexchanges[this_id][k][j][i]);
					else
						fprintf(fid,"%.12f", finalexchanges[this_id][k][j][i]);
					flag = 1;
					if((verbose > 1) && finalexchanges[this_id][k][j][i])
						printf("day: %d, exchange: %e (%.12f) from box%d-%d to %d-%d\n",
							dayj, finalexchanges[this_id][k][j][i], finalexchanges[this_id][k][j][i], i, j, boxdest[k][i], layerdest[k][i][j]);
				}

				for(k=boxndest[i]; k<st->fndest; k++){
					if(flag == 1)
						fprintf(fid,",_ ");
					else
						fprintf(fid,"_ ");
					flag = 1;
				}
				fprintf(fid,"\n  ");
			}
		}
	}

	flag = 0;
	/* Write destination boxes */
	fprintf(fid,";\n\n dest_b =\n  ");
	for(dayj=indexed_start; dayj<indexed_countend; dayj++){
		for(i=0; i<st->nbox; i++){
			for(j=0; j<st->wcnz; j++){
				for(k=0; k<boxndest[i]; k++){
					if(flag == 1)
						fprintf(fid,", %d", boxdest[k][i]);
					else
						fprintf(fid,"%d", boxdest[k][i]);
					flag = 1;
				}
				for(k=boxndest[i]; k<st->fndest; k++){
					if(flag == 1)
						fprintf(fid,",_ ");
					else
						fprintf(fid,"_ ");
					flag = 1;
				}
				fprintf(fid,"\n  ");
			}
		}
	}

	flag = 0;
	/* Write destination layers */
	fprintf(fid,";\n\n dest_k =\n  ");
	for(dayj=indexed_start; dayj<indexed_countend; dayj++){
		for(i=0; i<st->nbox; i++){
			for(j=0; j<st->wcnz; j++){
				for(k=0; k<boxndest[i]; k++){
					if(flag == 1)
						fprintf(fid,", %d", layerdest[k][i][j]);
					else
						fprintf(fid,"%d", layerdest[k][i][j]);
					flag = 1;
				}
				for(k=boxndest[i]; k<st->fndest; k++){
					if(flag == 1)
						fprintf(fid,",_ ");
					else
						fprintf(fid,"_ ");
					flag = 1;
				}
				fprintf(fid,"\n  ");
			}
		}
	}

	fprintf(fid,";\n\n}");

	/* Close the out file */
	fclose(fid);
	fclose(logfp);

	return;
}

/*************************************************************************************************
These routines write text version of temperature and salinity netcdf file
**************************************************************************************************/

void write_temp_salt_out(Stuff *st)
{
    FILE *outfp1 = NULL;
	FILE *outfp2 = NULL;
	FILE *outfp3 = NULL;
	FILE *outfp4 = NULL;


	/* Only do this if have temp and salt data */
	if( !st->ts.ntsfiles )
		return;

	/* Initialise files if necessary */
    if( !outfp1 )
		outfp1 = initTempOutFile(st);
	if( !outfp2 )
		outfp2 = initSaltOutFile(st);
	if(st->doPH){
		if(!outfp3)
			outfp3 = initpHOutFile(st);
	}

	if(st->doSWR){
		if(!outfp4)
			outfp4 = initSWROutFile(st);
	}

	/* Write files */

	writeTempSaltOut(outfp1, st, "temperature", temp_id);
	writeTempSaltOut(outfp2, st, "salinity", salt_id);
	if(st->doPH){
		writeTempSaltOut(outfp3, st, "pH", ph_id);
	}

	if(st->doSWR){
		writeTempSaltOut(outfp4, st, "swr", swr_id);
	}


	return;
}

FILE * initTempOutFile(Stuff *st)
/* Routine to initialise output file */
{
    FILE *fid;
    char *fname = st->temperatureOfname;


	// -1 as will stop in output file writing because of < in for statement

    if(verbose)
		printf("Create temperature output file\n");

	/* Create file */
    if( (fid=fopen(fname,"w")) == NULL )
		quit("initOutFile: Can't open %s\n",fname);

    /* File header info */
    fprintf(fid,"netcdf temp_data { \ndimensions:\n	t = UNLIMITED ; // (%d currently)\n", st->num_entries + 1);
	fprintf(fid,"	b = %d ;\n	z = %d ;\n",st->nbox, st->wcnz+1);
	fprintf(fid,"variables:\n");
	fprintf(fid, "	double t(t) ;\n		t:units = \"seconds since %d-01-01 00:00:00 +10\" ;\n		t:dt = %d. ;", st->REFyear, (int)(st->dt));
	fprintf(fid,"\n	double temperature(t, b, z) ;\n		temperature:_FillValue = 0. ;");
	fprintf(fid,"\n// global attributes:\n		:title = \"trivial\";\n		:geometry = \"%s\" ;",st->geomIfname);
	fprintf(fid,"\n		:parameters = \"\" ;\ndata:\n\n t = ");

    /* Return file pointer */
    return(fid);
}


FILE * initSaltOutFile(Stuff *st)
/* Routine to initialise output file */
{
    FILE *fid;
    char *fname = st->salinityOfname;
	// -1 as will stop in output file writing because of < in for statement

    if(verbose)
		printf("Create salinity output file\n");

	/* Create file */
    if( (fid=fopen(fname,"w")) == NULL )
		quit("initOutFile: Can't open %s\n",fname);

    /* File header info */
    fprintf(fid,"netcdf salt_data { \ndimensions:\n	t = UNLIMITED ; // (%d currently)\n",st->num_entries + 1);
	fprintf(fid,"	b = %d ;\n	z = %d ;\n",st->nbox, st->wcnz+1);
	fprintf(fid,"variables:\n");
	fprintf(fid, "	double t(t) ;\n		t:units = \"seconds since %d-01-01 00:00:00 +10\" ;\n		t:dt = %d. ;", st->REFyear, (int)(st->dt));
	fprintf(fid,"\n	double salinity(t, b, z) ;\n		salinity:_FillValue = 35. ;");
	fprintf(fid,"\n// global attributes:\n		:title = \"trivial\";\n		:geometry = \"%s\" ;",st->geomIfname);
	fprintf(fid,"\n		:parameters = \"\" ;\ndata:\n\n t = ");

    /* Return file pointer */
    return(fid);
}

FILE * initpHOutFile(Stuff *st)
/* Routine to initialise output file */
{
    FILE *fid;
    char *fname = st->pHOfname;


	// -1 as will stop in output file writing because of < in for statement

    if(verbose)
		printf("Create pH output file\n");

	/* Create file */
    if( (fid=fopen(fname,"w")) == NULL )
		quit("initOutFile: Can't open %s\n",fname);

    /* File header info */
    fprintf(fid,"netcdf pH_series { \ndimensions:\n	t = UNLIMITED ; // (%d currently)\n", st->num_entries + 1);
	fprintf(fid,"	b = %d ;\n	z = %d ;\n",st->nbox, st->wcnz+1);
	fprintf(fid,"variables:\n");
	fprintf(fid, "	double t(t) ;\n		t:units = \"seconds since %d-01-01 00:00:00 +10\" ;\n		t:dt = %d. ;", st->REFyear, (int)(st->dt));
	fprintf(fid,"\n	double pH(t, b, z) ;\n		pH:_FillValue = 0. ;");
	fprintf(fid,"\n// global attributes:\n		:title = \"trivial\";\n		:geometry = \"%s\" ;",st->geomIfname);
	fprintf(fid,"\n		:parameters = \"\" ;\ndata:\n\n t = ");

    /* Return file pointer */
    return(fid);
}


FILE * initSWROutFile(Stuff *st)
/* Routine to initialise output file */
{
    FILE *fid;
    char *fname = st->swrOfname;


	// -1 as will stop in output file writing because of < in for statement

    if(verbose)
		printf("Create SWR output file\n");

	/* Create file */
    if( (fid=fopen(fname,"w")) == NULL )
		quit("initOutFile: Can't open %s\n",fname);

    /* File header info */
    fprintf(fid,"netcdf SWR_series { \ndimensions:\n	t = UNLIMITED ; // (%d currently)\n", st->num_entries + 1);
	fprintf(fid,"	b = %d ;\n",st->nbox);
	fprintf(fid,"variables:\n");
	fprintf(fid, "	double t(t) ;\n		t:units = \"seconds since %d-01-01 00:00:00 +10\" ;\n		t:dt = %d. ;", st->REFyear, (int)(st->dt));
	fprintf(fid,"\n	double swr(t, b) ;\n		swr:_FillValue = 0. ;");
	fprintf(fid,"\n// global attributes:\n		:title = \"trivial\";\n		:geometry = \"%s\" ;",st->geomIfname);
	fprintf(fid,"\n		:parameters = \"\" ;\ndata:\n\n t = ");

    /* Return file pointer */
    return(fid);
}


void writeTempSaltOut(FILE *fid, Stuff *st, char *name, int index){
    int dayj, i, j, this_id = 0;

	int indexed_start = 0;
	int flag = 0;
	int countend, indexed_countend;

	if(st->dt > 86400){
		countend = (int)(st->tstop / (st->dt/86400));
		indexed_countend = (int)((st->tstop - st->tstart) / (st->dt/86400));
	}else{
		countend = (int)(st->tstop * (86400 / st->dt));
		indexed_countend = (int)((st->tstop - st->tstart) * (86400 / st->dt));
	}


	if(verbose){
		printf("Write to %s output file\n", name);
	}
	warn("Write to %s output file\n", name);

	/* Write time stamps */
	writeTimeSteps(st, fid, st->tstart, countend);

	flag = 0;
	/* Write temperatures */
	fprintf(fid,";\n\n %s =\n  ", name);

	printf("indexed_start = %d\n", indexed_start);
	printf("indexed_countend = %d\n", indexed_countend);


	for(dayj=indexed_start; dayj<indexed_countend; dayj++){

		if(st->pad_time)
			this_id = paddeddates[index][dayj];
		else
			this_id = dayj;

		//printf("dayj = %d, this_id = %d\n", dayj, this_id);
		for(i=0; i<st->nbox; i++){
			if( index == swr_id){
				if(flag == 1)
					fprintf(fid,", %f", finaltempsalt[this_id][0][i][index]);
				else
					fprintf(fid,"%f", finaltempsalt[this_id][0][i][index]);

				flag = 1;
			}else{
				/* Print water column values */
				for(j=0; j<st->wcnz; j++){
					if(flag == 1)
						fprintf(fid,", %f", finaltempsalt[this_id][j][i][index]);
					else
						fprintf(fid,"%f", finaltempsalt[this_id][j][i][index]);
					flag = 1;


				}
				/* Print sediment value (= bottom water column value) */
				fprintf(fid,",%f ", finaltempsalt[this_id][0][i][index]);
			}

			fprintf(fid,"\n  ");
		}
	}

	fprintf(fid,";\n\n}");
    fclose(fid);

	return;

}
