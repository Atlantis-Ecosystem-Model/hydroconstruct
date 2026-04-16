/**********************************************************************

    getaltinputdata.c

	by:	Beth Fulton
	date:	22/6/2005

	comments: This file is the code that reads in run parameters and 
			raw hydrodynamics flows, and temperature and salinity profiles
			using data from advective model (of the format type of Al Herman)

    assumptions: assumes raw input in one text files with entry values
	        per depth per face in sverdups for flow, flux deg C for temperature
			and flux PSU for salinity

*************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <math.h>
#include <string.h>
#include <sjwlib.h>
#include "mainexchange.h"

/***************************************************************************************************
Routine to get the exchange data out of the alt raw data file
****************************************************************************************************/
void get_alt_data(Stuff *st)
{
	char key[STSLEN];
	FILE *fp;
    int i, fce, dp, line, numline, ti, b, newtime, newupdate, lastupdate, ii, 
		entry_i, lastdp, valf;
	double time_of_year, valtime, valface, valbox, vallayer, valflow, valtemp, valsalt, 
		maxtime, counter, dirscale, valfloworig, totnum, this_tofy, len;


	time_of_year = 0;
	if(verbose)
		printf("Read in raw exchanges\n");

	/* Initialise face mapping - assumed to use bgm set as data not coming from netcdf file format */
	for(i=0; i<st->nface; i++)
		facemap[i] = -1;

	/** Al Herman style input **/
	/* Load lookup table so can convert Al Herman coordinates to bgm ids */
	Load_Lookup_Table(st);
	
	/* Loop through the exchange file */
	ti = 0;
	newtime = 0;
	maxtime = -MAXDOUBLE;	
	for(i=0; i<st->hd.nfiles; i++) {
		if(verbose){
			printf("Reading raw exchange file %s\n",st->hd.fname[i]);
			printf("0         100\n");
		}

		/* Open the exchange file */
		if( (fp=fopen(st->hd.fname[i],"r")) == NULL )
			quit("hydro_exchange_dataread: Can't open %s\n",st->hd.fname[i]);

		/* rewind file */
		fseek(fp, 0L, SEEK_SET);

		/* Read time */
		sprintf(key,"file_info");
		skipToKeyEnd(fp,key);
		if( fscanf(fp,"%d %d",&ti, &numline) != 2 )
			quit("the start line of the flow file must be give as 'file_info start_day number_of_lines' (e.g. file_info 1 255846) with the day given as an integer value\n");
		else if(verbose > 2)
			printf("looking at file %s with start time %d and %d lines\n", st->hd.fname[i], ti, numline);

		/* Reduce ti by 1 due to c arrays being base 0 */
		ti--;

		ti *= (int)(86400 / st->dt); 

		/* Reduce again as getting ti++ below when
		   valtime > maxtime check done */
		ti--;

		if(ti > st->ndata)
			quit("Not enough entries in the data array (i.e. number data points > array size - check ndata in prm file)\n");
		
		/* Initialise tracker ids */
		lastupdate = -1;
		entry_i = 1; // So when do first ti identification endup starting with 0 for entry_i
		lastdp = -1;
		/* Read raw exchanges */
		for(line=1; line<numline; line++){
			valtime = -1.0;
			valface = 0.0;
			valbox = 0.0;
			vallayer = 0.0;
			valflow = 0.0;
			valtemp = 0.0;
			valsalt = 0.0;

			/* rewind file - if appropriate */
			if(st->rewind)
				fseek(fp, 0L, SEEK_SET);

			sprintf(key,"Line%d",line);
			skipToKeyEnd(fp,key);

			if(verbose > 4)
				printf("reading line %d of %s ", line, st->hd.fname[i]);
			else if(verbose > 3)
				printf("%d %d\n",line,i);

			counter = 10 * line / numline;
			newupdate = (int)(counter);
			if(verbose && (newupdate > lastupdate)){
				printf("*, %d, %s", line / numline, st->hd.fname[i]);
				lastupdate = newupdate;

				if(lastupdate >= 10)
					printf("*\n");
			}

			/* Read data line */
			if( fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf",&valtime,&valface,&valbox,&vallayer,&valflow,&valtemp,&valsalt) != 7)
				quit("not enough entries for %s - with %lf %lf %lf %lf %lf %lf %lf\n", key, valtime, valface, valbox, vallayer, valflow, valtemp, valsalt);
			else if(verbose > 3)
				printf("read %s - with %lf %lf %lf %lf %lf %lf %lf\n", key, valtime, valface, valbox, vallayer, valflow, valtemp, valsalt);

			valfloworig = valflow;

			/* Up date valtime so chain year to year */
			valtime = i*365 + valtime;

			if(valtime > maxtime){
				if(verbose > 3 && maxtime > -MAXDOUBLE)
					printf("valtime %f, maxtime: %f - so ti up to %d\n", valtime, maxtime, ti+1);

				maxtime = valtime;
				for(ii=0; ii<st->nface; ii++){
					checkedthis[ii] = 0;
				}
				ti++;
				newtime++;

				if(ti > st->ndata)
					quit("Not enough entries in the data array (i.e. number data points > array size - check ndata in prm file)\n");
		
			} else
				newtime = 0;

			/* Convert time */
			if(newtime){
				exchangedate[ti][time_id] = (st->start_day + valtime) * st->dt;

				/* Multiple by 86400 / dt as may be tidal entries so must expand
				list to fit these extra entries in */
				time_of_year = (st->start_day + valtime) * 86400 / st->dt; 
				totnum = 365.0 * (86400 / st->dt); 

				if( time_of_year/totnum < 1.0 )
					this_tofy = floor(time_of_year/(86400 / st->dt) + 0.5);
				else
					this_tofy = floor((time_of_year/totnum - floor(time_of_year/totnum))*totnum + 0.5);

				exchangedate[ti][TofY_id] = this_tofy;
			}
			
			if(verbose > 3)
				printf("which matches %e (TofY = %e)\n", exchangedate[ti][time_id], exchangedate[ti][TofY_id]);

			/* Convert layer and face */
			Get_Info(st, valface, valbox, vallayer, &fce, &b, &dp, &dirscale, &len);

			valf = (int)valface;
			facemap[valf] = fce;
			faceinfo[i][valf][lbox_id] = b;
			faceinfo[i][valf][rbox_id] = valbox;
			faceinfo[i][valf][length_id] = len;

			if(dp != lastdp){
				for(ii=0; ii<st->nface; ii++){
					checkedthis[ii] = 0;
				}
				lastdp = dp;
			}

			if(!checkedthis[fce]){
				entry_i = 0;
				checkedthis[fce]++;
			} else
				entry_i = 1;

			if(verbose > 3)
				printf("calced box %d face %d, layer: %d\n", b, fce, dp);

			if(valflow < st->hd.missing_data){
				if(!st->sdiff)
					rawexchanges[ti][dp][valf][entry_i] += 0;
				else
					rawexchanges[ti][dp][valf][entry_i] += 0.00000001;
			}
			else {
				if(verbose > 3 && ti < 2)
					printf("setting actual value from rawexchanges[%d][%d][%d][%d] = %e to ",
						ti, dp, (int)(valface), entry_i, rawexchanges[ti][dp][valf][entry_i]);

				valflow *= (st->boxscale[b] * dirscale);
				rawexchanges[ti][dp][valf][entry_i] += valflow;

				if(verbose > 3 && ti < 2)
					printf("%e (dirscale %f, boxscale[%d] %f)\n", rawexchanges[ti][dp][valf][entry_i], dirscale, b, st->boxscale[b]);

			}

			if(verbose > 3 && ((ti < 2) || (ti > st->nhdsteps - 2)))
				printf("set rawexchanges[%d][%d][%d][%d] = %e (from l%f-b%f-dp%f, t=%f flow:%e (%e))\n",
					ti, dp, (int)(valface), entry_i, rawexchanges[ti][dp][valf][entry_i], valface, valbox, vallayer, valtime, valflow, valfloworig);

			/* Get temperature and salt data */
			tempsaltdate[ti][TofY_id] = (int)(time_of_year);

			if(valtemp < st->ts.temp_missing_data)
				rawtempsalt[ti][dp][valf][temp_id] += 0;
			else {
				valtemp *= (st->boxscale[b] * dirscale);
				rawtempsalt[ti][dp][valf][temp_id] += valtemp;
			}
			if(valsalt < st->ts.salt_missing_data)
				rawtempsalt[ti][dp][valf][salt_id] += 0;
			else { 
				valsalt *= (st->boxscale[b] * dirscale);
				rawtempsalt[ti][dp][valf][salt_id] += valsalt;
			}

			if(verbose > 2 && ((ti < 2) || (ti > st->nhdsteps - 2)))
				printf("stored rawexchanges[%d][%d][%d][%i]: %e, rawtemp[%d][%d][%d][%d]: %e (%e), rawsalt[%d][%d][%d][%d]: %e (%e)\n",
					ti, dp, valf, entry_i, rawexchanges[ti][dp][valf][entry_i], ti, dp, valf, temp_id, rawtempsalt[ti][dp][valf][temp_id], valtemp, ti, dp, valf, salt_id, rawtempsalt[ti][dp][valf][salt_id], valsalt);
		}
		/* Close the exchange file */
		fclose(fp);
    }

	if(verbose)
		printf("Finished reading raw exchange files\n");
	
	return;
}

/***************************************************************************************************
Routine to get the vertical exchange data out of the alt raw vert data file
****************************************************************************************************/
void get_alt_vert_data(Stuff *st)
{
	char key[STSLEN];
	FILE *fp;
    int i, fce, dp, line, numline, ti, b, newtime, newupdate, lastupdate,
		ii, entry_i, lastdp;
	double time_of_year, valtime, valface, valbox, vallayer, totnum, this_TofY,
		valflow, valvel, maxtime, counter, dirscale, valfloworig, thisdate, len;

	if(verbose)
		printf("Read in raw vertical exchanges\n");

	/** Al Herman style vertical input **/
	
	/* Loop through the vertical exchange file */
	ti = 0;
	newtime = 0;
	maxtime = -MAXDOUBLE;	
	for(i=0; i<st->vhd.nfiles; i++) {
		if(verbose){
			printf("Reading raw vertical exchange file %s\n",st->vhd.fname[i]);
			printf("0         100\n");
		}

		/* Open the exchange file */
		if( (fp=fopen(st->vhd.fname[i],"r")) == NULL )
			quit("hydro_exchange_dataread: Can't open %s\n",st->vhd.fname[i]);

		/* rewind file */
		fseek(fp, 0L, SEEK_SET);

		/* Read time */
		sprintf(key,"file_info");
		skipToKeyEnd(fp,key);
		if( fscanf(fp,"%d %d",&ti, &numline) != 2 )
			quit("the start line of the vertical flow file must be give as 'file_info start_day number_of_lines' (e.g. file_info 1 255846) with the day given as an integer value\n");
		else if(verbose > 2)
			printf("looking at file %s with start time %d and %d lines\n", st->vhd.fname[i], ti, numline);

		/* Reduce ti by 1 due to c arrays being base 0 */
		ti--;

		ti *= (int)(86400 / st->dt); 

		/* Reduce again as getting ti++ below when
		   valtime > maxtime check done */
		ti--;

		if(ti > st->ndata)
			quit("Not enough entries in the data array (i.e. number data points > array size - check ndata in prm file)\n");

		/* Initialise tracker ids */
		lastupdate = -1;
		entry_i = 1; // So when do first ti identification endup starting with 0 for entry_i
		lastdp = -1;
		/* Read raw exchanges */
		for(line=1; line<numline; line++){
			valtime = -1.0;
			valbox = 0.0;
			vallayer = 0.0;
			valvel = 0.0;
			valflow = 0.0;
			valface = 0.0;

			/* rewind file - if appropriate */
			if(st->rewind)
				fseek(fp, 0L, SEEK_SET);

			sprintf(key,"Line%d",line);
			skipToKeyEnd(fp,key);

			if(verbose > 4)
				printf("reading line %d of %s ", line, st->vhd.fname[i]);
			else if(verbose > 3)
				printf("%d %d\n",line,i);

			counter = 10 * line / numline;
			newupdate = (int)(counter);
			if(verbose && (newupdate > lastupdate)){
				printf("*, %d, %s", line / numline, st->vhd.fname[i]);
				lastupdate = newupdate;

				if(lastupdate >= 10)
					printf("*\n");
			}

			/* Read data line */
			if( fscanf(fp,"%lf %lf %lf %lf %lf",&valtime,&valbox,&vallayer,&valvel,&valflow) != 5)
				quit("not enough entries for %s - with %lf %lf %lf %lf %lf\n", key, valtime, valbox, vallayer, valvel, valflow);
			else if(verbose > 3)
				printf("read %s - with %lf %lf %lf %lf %lf\n", key, valtime, valbox, vallayer, valvel, valflow);

			valfloworig = valflow;

			/* Up date valtime so chain year to year */
			valtime = i*365 + valtime;

			if(valtime > maxtime){
				if(verbose > 3 && maxtime > -MAXDOUBLE)
					printf("valtime %f, maxtime: %f - so ti up to %d\n", valtime, maxtime, ti+1);

				maxtime = valtime;
				for(ii=0; ii<st->nbox; ii++){
					checkedthisbox[ii] = 0;
				}
				ti++;
				newtime++;

				if(ti > st->ndata)
					quit("Not enough entries in the data array (i.e. number data points > array size - check ndata in prm file)\n");
		
			} else
				newtime = 0;

			/* Convert time */
			if(newtime){
				thisdate = (st->start_day + valtime) * st->dt;

				if(thisdate != exchangedate[ti][time_id]){
					quit("line%d Vertical and horizontal dates don't match (hdate: %e, vdate: %e) - ti: %d, startday: %d, valtime: %e, dt: %d\n",
						line, exchangedate[ti][time_id], thisdate, ti, st->start_day, valtime, st->dt);
				}

				/* Multiple by 86400 / dt as may be tidal entries so must expand
				list to fit these extra entries in */
				time_of_year = (st->start_day + valtime) * (86400 / st->dt); 
				totnum = 365.0 * (86400 / st->dt); 

				if( time_of_year/totnum < 1.0 )
					this_TofY = floor(time_of_year/(86400 / st->dt) + 0.5);
				else
					this_TofY = floor((time_of_year/totnum - floor(time_of_year/totnum))*totnum + 0.5);

				if(this_TofY != exchangedate[ti][TofY_id]){
					quit("line%d Vertical and horizontal dates don't match (hdate %e, vdate: %d)\n",
						line, exchangedate[ti][TofY_id], this_TofY);
				}
			}
			
			if(verbose > 3)
				printf("which matches %e (TofY = %e)\n", exchangedate[ti][time_id], exchangedate[ti][TofY_id]);

			/* Convert layer and box */
			Get_Info(st, valface, valbox, vallayer, &fce, &b, &dp, &dirscale, &len);
			b = (int)(valbox);

			if(dp != lastdp){
				for(ii=0; ii<st->nbox; ii++){
					checkedthisbox[ii] = 0;
				}
				lastdp = dp;
			}

			if(!checkedthisbox[b]){
				entry_i = 0;
				checkedthis[b]++;
			} else
				entry_i = 1;

			if(verbose > 3)
				printf("calced box %d, layer: %d\n", b, dp);

			if(valflow < st->hd.missing_data){
				if(!st->sdiff)
					rawvertexchanges[ti][dp][b][entry_i] += 0;
				else
					rawvertexchanges[ti][dp][b][entry_i] += 0.00000001;
			}
			else {
				if(verbose > 3 && ti < 2)
					printf("setting actual value from rawvertexchanges[%d][%d][%d][%d] = %e to ",
						ti, dp, b, entry_i, rawexchanges[ti][dp][b][entry_i]);

				valflow *= (st->boxscale[b]);
				rawvertexchanges[ti][dp][b][entry_i] += valflow;

				if(verbose > 3 && ti < 2)
					printf("%e (vertdirscale %f, boxscale[%d] %f)\n", rawvertexchanges[ti][dp][b][entry_i], dirscale, b, st->boxscale[b]);

			}

			if(verbose > 3 && ((ti < 2) || (ti > st->nvhdsteps - 2)))
				printf("set rawvertexchanges[%d][%d][%d][%d] = %e (from b%f-dp%f, t=%f flow:%e (%e))\n",
					ti, dp, b, entry_i, rawvertexchanges[ti][dp][b][entry_i], valbox, vallayer, valtime, valflow, valfloworig);
		}
		/* Close the exchange file */
		fclose(fp);
    }

	if(verbose)
		printf("Finished reading raw vertical exchange files\n");
	
	return;
}

/***************************************************************************************************
Routine to get the temperature and salinity data out of the raw data files
****************************************************************************************************/
void get_alt_temp_salt(Stuff *st)
{
	char key[STSLEN];
	FILE *fp;
	int b, i, nk;
	double *val;

	/* Create local array */
	val = (double *)d_alloc1d(st->wcnz);

	if(verbose)
		printf("Read in initial temperatures and salinities\n");

	/** Get initial temperature and saliinity values to use as starting points
	so can derive final temperature and salinity values from the fluxes **/

	/* Open the temperature initial conditions file */
	if( (fp=fopen(st->tsTEMPinitfname,"r")) == NULL )
		quit("hydro_exchange_dataread: Can't open %s\n",st->tsTEMPinitfname);

	/* Read temperature initial conditions data */
	for(b=0; b<st->nbox; b++){
		for(i=0; i<st->wcnz; i++){
			val[i] = 0;
		}
	
		/* rewind file */
		fseek(fp, 0L, SEEK_SET);

		/* Read data */
		sprintf(key,"Box%d",b);
		skipToKeyEnd(fp,key);

		for(i=0; i<st->wcnz; i++){
			if( (nk = fscanf(fp,"%lf ", &val[i])) != 1){
				quit("not enough entries for %s - can't load layer %d\n", key, i);
			}

			TempSaltInit[b][i][temp_id] = val[i];
		}
	}

	/* Close the temperature file */
	fclose(fp);

	/* Open the salinity initial conditions file */
	if( (fp=fopen(st->tsSALTinitfname,"r")) == NULL )
		quit("hydro_exchange_dataread: Can't open %s\n",st->tsSALTinitfname);

	/* Read salinity initial conditions data */
	for(b=0; b<st->nbox; b++){
		for(i=0; i<st->wcnz; i++){
			val[i] = 0;
		}
	
		/* rewind file */
		fseek(fp, 0L, SEEK_SET);

		/* Read data */
		sprintf(key,"Box%d",b);
		skipToKeyEnd(fp,key);

		for(i=0; i<st->wcnz; i++){
			if( (nk = fscanf(fp,"%lf ", &val[i])) != 1){
				quit("not enough entries for %s - can't load layer %d\n", key, i);
			}

			TempSaltInit[b][i][salt_id] = val[i];
		}
	}

	/* Close the salinity file */
	fclose(fp);

	/** Loading of salinity and temperature fluxes done in get_alt_data() above **/
	
	/* Free local array */
	d_free1d(val);

	return;
	
}

/****************************************************************************************************
Routine to check look-up table for conversion of Al Herman line coordinates to bgm face ids
*****************************************************************************************************/
void Load_Lookup_Table(Stuff *st)
{
	char key[STSLEN];
	FILE *fp;
	int line, maxdeep, layer, dp;
	double valface, valaltbox, valbox, valscale, vallength, depth;
    
	maxdeep = -MAXINT;
	for(layer = 0; layer<st->wcnz; layer++){
		if(st->nominal_dz[layer] > maxdeep)
			maxdeep = (int)(st->nominal_dz[layer]);
	}

	/* Read and load lookup table */
	if(verbose)
		printf("Reading lookup table file %s\n",st->lookupIfname);

	/* Initialise look-up table array */
	lookuptable = (int **)i_alloc2d(nalt_data_id,st->altfaces);
	lookuplength = (double *)d_alloc1d(st->altfaces);
	lookupdir = (double *)d_alloc1d(st->altfaces);
	lookuplayer = (int *)i_alloc1d(maxdeep);
	
	/* Open the lookup table file */
	if( (fp=fopen(st->lookupIfname,"r")) == NULL )
		quit("Lookup table file: Can't open %s\n",st->lookupIfname);

	/* rewind file */
	fseek(fp, 0L, SEEK_SET);

	/* Read table file */
	for(line=0; line<st->altfaces; line++){
		valface = 0.0;
		valaltbox = 0.0;
		valbox = 0.0;
		valscale = 1.0;

		/* rewind file */
		fseek(fp, 0L, SEEK_SET);

		sprintf(key,"AlFace%d",line);
		skipToKeyEnd(fp,key);

		/* Read data line */
		if( fscanf(fp,"%lf %lf %lf %lf %lf",&valface,&valaltbox,&valbox,&valscale,&vallength) != 5)				
			quit("not enough entries for %s - with %lf %lf %lf %lf %lf\n", key, valface, valaltbox, valbox, valscale, vallength);

		if(verbose > 3)
			printf("readin altface-%d, face-%lf altbox-%lf box-%lf scale-%lf length-%lf\n", line, valface, valaltbox, valbox, valscale, vallength);

		/* Store table */
		lookuptable[line][altface_id] = line;
		lookuptable[line][bgmface_id] = (int)(valface);
		lookuptable[line][altbox_id] = (int)(valaltbox);
		lookuptable[line][bgmbox_id] = (int)(valbox);
		lookupdir[line] = valscale;
		//lookuplength[line] = metres_in_degree * vallength;
		lookuplength[line] = vallength;
		
		if(verbose > 3)
			printf("stored altface-%d, face-%d altbox-%d box-%d scale-%lf\n", 
				lookuptable[line][altface_id], lookuptable[line][bgmface_id], lookuptable[line][altbox_id], lookuptable[line][bgmbox_id], lookupdir[line]);
	}

	/* Load-lookup depths */
	for(layer = 0; layer<maxdeep; layer++){
		lookuplayer[layer] = 0;
		depth = 0;
		for(dp=0; dp<st->wcnz; dp++){
			depth += st->nominal_dz[dp];
			if(layer < depth){
				lookuplayer[layer] = dp;
				if(verbose > 3)
					printf("layer %d matches lookup %d\n", layer, lookuplayer[layer]);
				break;
			}
		}
	}

	/* Close the exchange file */
	fclose(fp);

	return;
}

/****************************************************************************************************
Routine to check look-up table for conversion of Al Herman coordinates to bgm ids
*****************************************************************************************************/
void Get_Info(Stuff *st, double valface, double valbox, double vallayer, int *bgm_face, int *bgm_box, 
			  int *bgm_layer, double *dirscale, double *len) 
{
	int ans;
	int checkface = (int)(valface);
	int checkbox = (int)(valbox);
	int checklayer = (int)(vallayer);
	int totbgmface = st->nfacegeo - 1;
	int totbgmbox = st->nbox - 1;
	
	/* Initialise responses */
	*bgm_face = 0;
	*bgm_box = 0;
	*bgm_layer = 0;

	/* Get bgm face id */
	ans = lookuptable[checkface][bgmface_id];
	if((ans > totbgmface) && (ans != st->notbgmface))
		quit("Trying to access a face that doesn't exist (altface_id: %d, bgmfaceid: %d, nbgmface: %d)\n",
			checkface, bgm_face, totbgmface);
	*bgm_face = ans;

	/* Get bgm box id */
	ans = lookuptable[checkface][bgmbox_id];
	if(ans > totbgmbox)
		warn("Trying to access a box that doesn't exist (altbox_id: %d, bgmboxid: %d, nbgmbox: %d)\n",
			checkbox, bgm_box, totbgmbox);
	*bgm_box = ans;

	/* Get bgm layer id */
	*bgm_layer = lookuplayer[checklayer];

	/* Get bgm directional scalar */
	*dirscale = lookupdir[checkface];

	/* Get face length */
	*len = lookuplength[checkface];

}

