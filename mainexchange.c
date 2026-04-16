/**********************************************************************

    Exchange calculating program to accompany box model appropach
	by:	Beth Fulton
	date: 20/10/2003

	comments: This file is the code that converts hydrodynamic data
		files into boxmodel compatible netcdf files

    revisions:	28/11/2004 updated to include multiple read-in files and
				temperature and salinity data

				11/12/2006 added the code to read in excel estimated flows
				(under col_face_type input_type)

				14-09-2009 Bec Gorton
				Fixed a bug in the code that set the initial values
				for the padDate array. This array was being indexed in the
				incorrect order.

				01-10-2009 Bec Gorton
				Added additional code to free up arrays.
*************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <math.h>
#include <string.h>
#include <sjwlib.h>
#include "mainexchange.h"

/* Global variables */
int killed = 0;
int verbose = 1;
double reference_time;

/* Prototypes for routines used in this file */
void    setupStuff(int argc, char *argv[], Stuff *st);
void    Usage(void);
void    setkill(void);

/*******************************************************************************
Main run routine
*/
int main(int argc, char *argv[])
{
	Stuff st;

	/* Set up array space */
	setupStuff(argc,argv,&st);

#if 0
	/* Set signal to allow clean up on SIGTERM */
	signal(SIGTERM,setkill);
#endif

	/* Do conversion */
	ExchangeStuff(&st);

	/* Determine Temperatures and Salinities */
	CalcTempSalt(&st);

	/* Pad out timeseries */
	if(st.pad_time == 1)
		PadTimeSeries(&st);

	/* WriteOutput */
	WriteNewOutput(&st);

	freeStuff(&st);
	/* End program */
	exit(0);
}


void freeStuff(Stuff *st){

	/* Free the geometry */
	freeGeom(st);


	i_free3d(col_faceid);
	d_free3d(col_faceProportion);
	d_free2d(depth_layer);


	d_free3d(TempSaltInit );
	d_free2d(lookupdz);

	d_free1d(st->nominal_dz);


	d_free4d(rawexchanges);
	i_free1d(checkedthis);
	i_free1d(facemap);
	i_free1d(test_ok);
	i_free1d(next_test);

	d_free4d(rawvertexchanges);
	i_free1d(checkedthisbox);
	i_free1d(boxmap);

	d_free2d(st->boxdz);
	d_free2d(st->boxvol);
	i_free1d(st->boxnz);
	i_free1d(st->boundaries);
	d_free1d(st->boxscale);

	d_free1d(st->estinflux);
	d_free1d(st->estoutflux);
	i_free1d(st->estnum);

	d_free4d(finalexchanges);

	i_free2d(leftrightcheck);
	i_free2d(topbotcheck);

	i_free1d(boxndest);

	i_free2d(boxdest);
	i_free3d(layerdest);

	i_free2d(paddeddates);
	d_free2d(exchangedate);

	d_free4d(rawtempsalt);
	d_free4d(finaltempsalt);

	d_free3d(boxstats);

	d_free2d(tempsaltdate);
	d_free4d(tottempsaltflux);

	d_free3d(faceinfo);
	i_free3d(tempfaceinfo);
	if(boxLookup)
		i_free1d(boxLookup);

	if (st->hd.fname != NULL)
		c_free2d(st->hd.fname);

	if (st->vhd.nfiles > 0)
		c_free2d(st->vhd.fname);

	if (st->ts.tsfname != NULL)
		c_free2d(st->ts.tsfname);

}

/*****************************************************************************************************
Setup box structure and read in raw data
******************************************************************************************************/

void setupStuff(int argc, char *argv[], Stuff *st)
{
	int i, j, b, k, nentries;
	double current_time, time_correction;
	char key[STSLEN];
	char *progname = argv[0];

	/* Clear the stuff space */
	memset(st, 0, sizeof(Stuff));
	st->doPH = 0;
	st->doSWR = 0;

	/* Process arguments */
	if( argc < 8 ) Usage();
	while( --argc > 0 ) {
		if( (*++argv)[0] == '-' ) {
			switch( (*argv)[1] ) {
			case 'f':
				if( argc < 2 ) Usage();
				strncpy(st->transportOfname,*++argv, STSLEN - strlen(st->transportOfname) - 1);
				argc --;
				break;
			case 't':
				if( argc < 2 ) Usage();
				strncpy(st->temperatureOfname,*++argv, STSLEN - strlen(st->temperatureOfname) - 1);
				argc --;
				break;
			case 's':
				if( argc < 2 ) Usage();
				strncpy(st->salinityOfname,*++argv, STSLEN - strlen(st->salinityOfname) - 1);
				argc --;
				break;
			case 'r':
				if( argc < 2 ) Usage();
				strncpy(st->runprmIfname,*++argv, STSLEN - strlen(st->runprmIfname) - 1);
				argc--;
				break;
			case 'p':
				if( argc < 2 ) Usage();
				strncpy(st->pHOfname,*++argv, STSLEN - strlen(st->pHOfname) - 1);
				st->doPH = 1;
				argc--;
				break;
			case 'w':
				if( argc < 2 ) Usage();
				strncpy(st->swrOfname,*++argv, STSLEN - strlen(st->swrOfname) - 1);
				st->doSWR = 1;
				argc--;
				break;
			default: Usage(); break;
			}
		} else
			Usage();
	}

	/* Check that all required arguments are present */
	if( !st->transportOfname[0] || !st->runprmIfname[0] || !st->temperatureOfname[0] || !st->salinityOfname[0]){
		printf("One of the required filenames is missing\n");
		Usage();
	}

	/* Write arguments to parameter string */
	if(st->doPH){
		if(st->doSWR)
			sprintf(st->params,"%s -f %s -t %s -s %s -p %s -w %s -r %s",progname,st->transportOfname,st->temperatureOfname,st->salinityOfname,st->pHOfname, st->swrOfname, st->runprmIfname);
		else
			sprintf(st->params,"%s -f %s -t %s -s %s -p %s -r %s",progname,st->transportOfname,st->temperatureOfname,st->salinityOfname,st->pHOfname, st->runprmIfname);
	}else{
		if(st->doSWR)
			sprintf(st->params,"%s -f %s -t %s -s %s -w %s -r %s",progname,st->transportOfname,st->temperatureOfname,st->salinityOfname, st->swrOfname, st->runprmIfname);
		else
			sprintf(st->params,"%s -f %s -t %s -s %s -r %s",progname,st->transportOfname,st->temperatureOfname,st->salinityOfname,st->runprmIfname);
	}

	/* Set error routine to quit if parameters not found */
    set_keyprm_errfn(quit);

	/* Read the model run parameters */
    readRunParams(st->runprmIfname,st);

    if(verbose)
		printf("Initialising exchange storage matrices\n");

	/* Initiate exchange matrices */
	if(!st->nhdsteps)
		warn("How can have no hdsteps - may want to check the code!\n");

	st->nfacegeo = st->nface;

	if(!st->input_style)
		st->nface = st->altfaces + 1;

	if(verbose)
		printf("Initialising rawexchanges[nstep=%d][wcnz=%d][nface=%d][nflow=2]\n",st->nhdsteps,st->wcnz,st->nface);

	rawexchanges = (double ****)d_alloc4d(2,st->nface,st->wcnz,st->nhdsteps); // two entries per face as can flow two ways
	checkedthis = (int *)i_alloc1d(st->nface);
	facemap = (int *)i_alloc1d(st->nface);

	st->nitcount = 6;
	test_ok = (int *)i_alloc1d(st->nitcount);
	next_test = (int *)i_alloc1d(st->nitcount);

	for(j=0; j<st->nitcount; j++){
		test_ok[j] = 1;
		next_test[j] = (int)(pow(10,j));
	}

	for(j=0; j<st->nhdsteps; j++){
		for(b=0; b<st->wcnz; b++){
			for(i=0; i<st->nface; i++){
				for(k=0; k<2; k++){
					rawexchanges[j][b][i][k] = 0.0;
				}
			}
		}
	}

	for(i=0; i<st->nface; i++){
		checkedthis[i] = 0;
		facemap[i] = -1;
	}

	if(verbose)
		printf("Initialising rawvertexchanges[nstep=%d][wcnz=%d][nbox=%d][nflow=2]\n",st->nvhdsteps,st->wcnz,st->nbox);


	if(st->input_style != cdf_type_id)
		st->nvhdsteps = st->nhdsteps;

	rawvertexchanges = (double ****)d_alloc4d(2,st->nbox,st->wcnz,st->nvhdsteps); // two entries per face as can flow two ways
	checkedthisbox = (int *)i_alloc1d(st->nbox);
	boxmap = (int *)i_alloc1d(st->nbox);

	for(j=0; j<st->nvhdsteps; j++){
		for(b=0; b<st->wcnz; b++){
			for(i=0; i<st->nbox; i++){
				for(k=0; k<2; k++){
					rawvertexchanges[j][b][i][k] = 0.0;
				}
			}
		}
	}

	for(i=0; i<st->nbox; i++){
		checkedthisbox[i] = 0;
		boxmap[i] = -1;
	}

	if(verbose)
		printf("Initialising boxdz[nbox=%d][wcnz=%d]\n",st->nbox,st->wcnz);

	st->boxdz = (double **)d_alloc2d(st->wcnz,st->nbox);

	for(i=0; i<st->nbox; i++){
		for(k=0; k<st->wcnz; k++){
			st->boxdz[i][k] = 0;
		}
	}

	if(verbose)
		printf("Initialising finalexchanges[nstep=%d][ndest=%d][wcnz=%d][nbox=%d]\n",st->nhdsteps,st->ndest,st->wcnz,st->nbox);

	finalexchanges = (double ****)d_alloc4d(st->nbox,st->wcnz,st->ndest,st->nhdsteps);

	for(k=0; k<st->nhdsteps; k++){
		for(j=0; j<st->ndest; j++){
			for(b=0; b<st->wcnz; b++){
				for(i=0; i<st->nbox; i++){
					finalexchanges[k][j][b][i] = 0.0;
				}
			}
		}
	}

	if(verbose)
		printf("Initialising leftrightcheck and topbotcheck\n");

	leftrightcheck = (int **)i_alloc2d(st->wcnz+1,st->nface);
	for(i=0; i<st->nface; i++){
		for(b=0; b<st->wcnz+1; b++){
			leftrightcheck[i][b] = 0;
		}
	}

	topbotcheck = (int **)i_alloc2d(st->wcnz+1,st->nbox);
	for(i=0; i<st->nbox; i++){
		for(b=0; b<st->wcnz+1; b++){
			leftrightcheck[i][b] = 0;
		}
	}

	if(verbose > 3)
		printf("Initialising boxndest[nbox=%d]\n",st->nbox);

	boxndest = (int *)i_alloc1d(st->nbox);

	if(verbose)
		printf("Initialising boxdest[ndest=%d][nbox=%d]\n",st->ndest,st->nbox);

	boxdest = (int **)i_alloc2d(st->nbox,st->ndest);

	for(i=0; i<st->nbox; i++){
		boxndest[i] = 0;
		for(b=0; b<st->ndest; b++){
			boxdest[b][i] = -1;
		}
	}

	if(verbose)
		printf("Initialising layerdest[ndest=%d][nbox=%d][nlayer=%d]\n",st->ndest,st->nbox,st->wcnz);

	layerdest = (int ***)i_alloc3d(st->wcnz,st->nbox,st->ndest);

	//if(verbose)
		printf("Initialising exchangedate[nstep=%d][2]\n",st->nhdsteps);

	exchangedate = (double **)d_alloc2d(3,st->nhdsteps);

	for(i=0; i<st->nhdsteps; i++){
		exchangedate[i][time_id] = 0;
		exchangedate[i][TofY_id] = 0;
		exchangedate[i][done_id] = 0;
	}

	/* Initiate temperature & salinity matrices */
	switch(st->input_style){
		case col_data_type:
		case new_col_data_type:
			st->ntssteps = (int)(st->ts.ntsfiles * 31 * 86400 / st->dt);
			break;
		case flat_type_id:
			if(st->ts.ntsfiles > 0)
				st->ntssteps = st->ts.ntsfiles;
			else
				st->ntssteps = 1;
			break;
		default:
			st->ntssteps = st->ndata;
	}


	if(st->tsflagflux && (st->ntssteps != st->nhdsteps))
		quit("For case where have fluxes of temperature and salinity not absolute values then must have as many water flux data points as salinity and temperature fluxes\n");


	if(!st->input_style)
		nentries = st->altfaces;
	else
		nentries = st->nbox;

	if(verbose)
		printf("Initialising rawtempsalt[nsteps=%d][wcnz=%d][nbox=%d][2]\n",st->ntssteps,st->wcnz,nentries);

	rawtempsalt = (double ****)d_alloc4d(num_TSP_values,nentries,st->wcnz,st->ntssteps);
	finaltempsalt = (double ****)d_alloc4d(num_TSP_values,st->nbox,st->wcnz,st->ntssteps);

	for(i=0; i<nentries; i++){
		for(b=0; b<st->wcnz; b++){
			for(j=0; j<st->ntssteps; j++){
				for(k = 0; k < num_TSP_values; k++){
					rawtempsalt[j][b][i][k] = 0.0;
				}
			}
		}
	}

	for(i=0; i<st->nbox; i++){
		for(b=0; b<st->wcnz; b++){
			for(j=0; j<st->ntssteps; j++){
				for(k=0; k<num_TSP_values; k++){
					finaltempsalt[j][b][i][k] = 0;
				}
			}
		}
	}

	if(verbose)
		printf("Initialising boxstats[nbox=%d][wcnz=%d][2]\n",st->nbox,st->wcnz);
	boxstats = (double ***)d_alloc3d(2,st->wcnz,st->nbox);

	for(i=0; i<st->nbox; i++){
		for(b=0; b<st->wcnz; b++){
			for(k=0; k<2; k++){
				boxstats[i][b][k] = 0.0;
			}
		}
	}

	if(verbose)
		printf("Initialising tempsaltdate[nsteps=%d][2]\n",st->ntssteps);

	printf("st->ntssteps = %d\n", st->ntssteps);
	tempsaltdate = (double **)d_alloc2d(2,st->ntssteps);

	for(i=0; i<st->ntssteps; i++){
		tempsaltdate[i][0] = 0.0;
		tempsaltdate[i][1] = 0.0;
	}

	if(verbose)
		printf("Initialising tottempsaltflux[nsteps=%d][wcnz=%d][nbox=%d][ndata=%d]\n",st->ntssteps,st->wcnz,st->nbox,2);

	tottempsaltflux = (double ****)d_alloc4d(num_TSP_values,st->nbox,st->wcnz,st->ntssteps);

	for(k=0; k<st->ntssteps; k++){
		for(b=0; b<st->wcnz; b++){
			for(i=0; i<st->nbox; i++){
				for(j=0; j<num_TSP_values; j++){
					tottempsaltflux[k][b][i][j] = 0.0;
				}
			}
		}
	}

	/* Padded date id array */
	if(verbose)
		printf("Initialising paddeddates[nt=%d][3]\n",st->nt);

//	paddeddates = (int **)i_alloc2d(3,st->nt);
//
//	for(i=0; i<st->nt; i++){
//		paddeddates[i][0] = 0;
//		paddeddates[i][1] = 0;
//		paddeddates[i][2] = 0;
//	}

	paddeddates = (int **)i_alloc2d(st->nt, num_TSP_values + 1);

	for(i=0; i<st->nt; i++){
		paddeddates[0][i] = 0;
		paddeddates[1][i] = 0;
		paddeddates[2][i] = 0;
	}


	if(verbose)
		printf("Initialising faceinfo[%d][st->nface=%d][st->hd.nfiles=%d]\n",num_faceinfo, st->nface, st->hd.nfiles);

	faceinfo = (double ***)d_alloc3d(num_faceinfo,st->nface,st->hd.nfiles);
	for(i=0; i<st->hd.nfiles; i++){
		for(k=0; k<st->nface; k++){
			for(b=0; b<num_faceinfo; b++){
				faceinfo[i][k][b] = -1;
			}
		}
	}

	tempfaceinfo = (int ***)i_alloc3d(num_faceinfo,st->nface,st->nbox);
	for(i=0; i<st->nbox; i++){
		for(k=0; k<st->nface; k++){
			for(b=0; b<num_faceinfo; b++){
				tempfaceinfo[i][k][b] = -1;
			}
		}
	}
	/* Read exchange and temperature/salinity data */
    get_hydro(st);

	st->ntssteps = st->nhdsteps;

	/* Reset times with respect to earliest date */
	if(st->timecorrect_on == 1){
		/* Exchanges */
		if(st->monthEx < 10){
			if(st->dayEx<10){
				sprintf(key,"%d-0%d-0%d 00:00:00", st->minExyear, st->monthEx, st->dayEx);
			} else {
				sprintf(key,"%d-0%d-%d 00:00:00", st->minExyear, st->monthEx, st->dayEx);
			}
		} else {
			if(st->dayEx<10){
				sprintf(key,"%d-%d-0%d 00:00:00", st->minExyear, st->monthEx, st->dayEx);
			} else {
				sprintf(key,"%d-%d-%d 00:00:00", st->minExyear, st->monthEx, st->dayEx);
			}
		}
		current_time = GetTime(key);
		time_correction = current_time;

		warn("Correcting exchange dates for %s by %e\n", key, time_correction);

		for(i=0; i<st->nhdsteps; i++) {
			exchangedate[i][time_id] -= time_correction;
		}

		/* Temperature and salinity */
		if(st->monthTS < 10){
			if(st->dayTS<10){
				sprintf(key,"%d-0%d-0%d 00:00:00", st->minTSyear, st->monthTS, st->dayTS);
			} else {
				sprintf(key,"%d-0%d-%d 00:00:00", st->minTSyear, st->monthTS, st->dayTS);
			}
		} else {
			if(st->dayTS<10){
				sprintf(key,"%d-%d-0%d 00:00:00", st->minTSyear, st->monthTS, st->dayTS);
			} else {
				sprintf(key,"%d-%d-%d 00:00:00", st->minTSyear, st->monthTS, st->dayTS);
			}
		}
		current_time = GetTime(key);
		time_correction = current_time;

		warn("Correcting temp and salt dates for %s by %e\n", key, time_correction);

		for(i=0; i<st->ntssteps; i++) {
			tempsaltdate[i][time_id] -= time_correction;
		}
	} else if(st->timecorrect_on == 2){
		time_correction = MAXDOUBLE;
		for(i=0; i<st->nhdsteps; i++) {
			if(exchangedate[i][time_id] < time_correction)
				time_correction = exchangedate[i][time_id];
		}
		for(i=0; i<st->nhdsteps; i++) {
			exchangedate[i][time_id] -= time_correction;
		}
		time_correction = MAXDOUBLE;
		for(i=0; i<st->ntssteps; i++) {
			if(tempsaltdate[i][time_id] < time_correction)
				time_correction = tempsaltdate[i][time_id];
		}
		for(i=0; i<st->ntssteps; i++) {
			tempsaltdate[i][time_id] -= time_correction;
		}
	}

	return;

}

/**************************************************************************************************
Usage and termination routines
***************************************************************************************************/
void Usage()
{
    quit("Usage: hydroconstruct -f output_flow.txt -t output_temp.txt -s output_salt.txt [-p ph_output.txt] [-w swr_output.txt] -r run.prm\n");
}

void setkill()
{
    killed = 1;
}
