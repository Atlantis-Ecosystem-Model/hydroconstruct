/**********************************************************************

    Code for padding out time series
	by:	Beth Fulton
	date:	28/11/2004

	comments: This file is the code for padding out timeseries by reusing values
		from later in series to fill gaps at same calendar date earlier in the series

    revisions:

*************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sjwlib.h>
#include "mainexchange.h"

/***************************************************************************************************
Main body of calculations and output of exchanges
****************************************************************************************************/
void PadTimeSeries(Stuff *st)
{
	int i, k, min_id, min_tofy, max_tofy, max_id, max_count, update_id, current_id,
		this_id, step1, dayj, which_i,	min_tstofy, max_tstofy, checked, end_yr,
		yr_correct, mth;
	double min_time, max_time, min_tstime, max_tstime, max_tsday, min_day, max_day,
		min_tsday, last_day, next_day, end_yr_day;
	FILE *logfp = NULL;
	char logfname[STSLEN];

	min_id = 0;
	max_id = 0;
	min_day = 0;
	end_yr_day = 365.0 * (86400 / st->dt);
	end_yr = 0;

	if(st->verbose_pad){
		sprintf(logfname,"%s","padlog.txt");

		/* Create verbose log file */
		if( (logfp=fopen(logfname,"w")) == NULL )
			quit("Can't open log file%s\n",logfname);
		else
			printf("Created %s\n", logfname);
	}

	if(verbose)
		printf("Begin padding timeseries\n");

	/****** Padding temperature and salinity profile time series (identify which entry ids
	to use for the time series padding) ****/
	if(verbose)
		printf("Begin padding temperature and salinity timeseries\n");


	min_time = MAXINT;
	max_time = 0;
	if((st->input_style != flat_type_id) && (st->input_style != col_face_type)){
		min_id = 0;
		min_time = MAXINT;
		for(i=0; i<st->ntssteps; i++) {
			if(st->verbose_pad){
				fprintf(logfp,"i: %d, tsmin_time: %e, tsTofY: %e, tstime: %e\n", i, min_time, tempsaltdate[i][TofY_id], tempsaltdate[i][time_id]);
			}
			if(tempsaltdate[i][time_id] < min_time){
				min_id = i;
				min_day = tempsaltdate[i][TofY_id];
				min_time = tempsaltdate[i][time_id];
			}
		}
		max_time = 0;
		for(i=0; i<st->ntssteps; i++) {
			if(tempsaltdate[i][time_id] > max_time){
				/* Divide by dt not 86400 as may be tidal entries not daily */
				max_id = i;
				max_day = tempsaltdate[i][TofY_id];
				max_time = tempsaltdate[i][time_id];
			}
		}
	}
	min_tofy = MAXINT;
	max_tofy = -1;
	for(i=0; i<st->ntssteps; i++) {
		if(tempsaltdate[i][TofY_id] < min_tofy){
			min_tofy = (int)(tempsaltdate[i][TofY_id]);
		}
		if(tempsaltdate[i][TofY_id] > max_tofy){
			max_tofy = (int)(tempsaltdate[i][TofY_id]);
		}
	}

printf("min_tofy = %d,max_tofy= %d\n", min_tofy, max_tofy);
	//fprintf(logfp, "min_tofy = %d, max_tofy = %d, min_time = %e, max_time = %e\n", min_tofy, max_tofy, min_time, max_time);
	if(min_tofy > (st->tstop * (86400 / st->dt)))
		quit("No valid entries (earliest data calendar day (%d) > latest required calendar day (%d))\n", min_tofy, st->tstop);

	current_id = 0;
	yr_correct = 0;
	last_day = -1;
	next_day = min_day;
	update_id = 0;
	max_count = st->ntssteps;

	if(st->input_style == flat_type_id){
		min_id = st->TSmin_id;
		min_time =tempsaltdate[min_id][time_id];
		min_day = 0;
		for(mth=0; mth<st->monthTS-1; mth++){
			switch(mth){
			case 0:
			case 2:
			case 4:
			case 6:
			case 7:
			case 9:
			case 11:
				min_day += 31;
				break;
			case 1:
				min_day += 28;
				break;
			case 3:
			case 5:
			case 8:
			case 10:
				min_day += 30;
				break;
			}
		}
		min_day += st->dayTS;

		max_day = 0;
		for(mth=0; mth<st->mmonthTS-1; mth++){
			switch(mth){
			case 0:
			case 2:
			case 4:
			case 6:
			case 7:
			case 9:
			case 11:
				max_day += 31;
				break;
			case 1:
				max_day += 28;
				break;
			case 3:
			case 5:
			case 8:
			case 10:
				max_day += 30;
				break;
			}
		}
		max_day += st->mdayTS;
		max_id = st->TSmax_id;
		max_time =tempsaltdate[max_id][time_id];
	}

	/* Loop over each data point time */
	for(i=0; i<st->nt; i++){
		step1 = (int)(i/365);
		dayj = i - step1 * 365;
		if(i < min_day){
			/* Padding before */
			if(min_tofy > min_day){
				/* Back fill with smallest available */
				warn("Had to use data for day %d for day %d as back fill = only option available\n", min_day, i);
				current_id = min_id;
				paddeddates[temp_id][i] = min_id;
				paddeddates[salt_id][i] = min_id;
				paddeddates[ph_id][i] = min_id;

			} else {
				/* Fill with entries with appropriate TofY */
				for(k=0; k<st->ntssteps; k++) {
					if(tempsaltdate[k][TofY_id] == i) {
						last_day = tempsaltdate[k][TofY_id];
						current_id = k;
						break;
					} else if((tempsaltdate[k][TofY_id] > last_day) &&  (tempsaltdate[k][TofY_id] < i)){
						last_day = tempsaltdate[k][TofY_id];
						current_id = k;
						break;
					}
				}
				paddeddates[temp_id][i] = current_id;
				paddeddates[salt_id][i] = current_id;
				paddeddates[ph_id][i] = current_id;


				if(st->verbose_pad)
					fprintf(logfp,"tsi: %d, using id %d last_day: %e\n", i, current_id, last_day);

			}
		} else if (update_id < max_count-1){
			if(st->verbose_pad)
				fprintf(logfp,"tsi: %d (%d), last_day: %e, next_day: %e", i, (int)(i - end_yr_day * yr_correct), last_day, next_day);

			/* Straight read-off with interpolation */
			if(i == min_day){
				last_day = min_day;
				current_id = min_id;
				update_id = min_id + 1;
				next_day = tempsaltdate[update_id][TofY_id];
			} else {
				if((!end_yr && ((i - end_yr_day * yr_correct) < next_day)) || (end_yr && ((i - end_yr_day * yr_correct) < end_yr_day)));
					// Nothing to update
				else {
					last_day = tempsaltdate[current_id][TofY_id];
					current_id++;
					update_id++;
					next_day = tempsaltdate[update_id][TofY_id];
					if(last_day > next_day)
						end_yr = 1;
					else{
						if(end_yr)
							yr_correct++;
						end_yr = 0;
					}
				}
			}
			paddeddates[temp_id][i] = current_id;
			paddeddates[salt_id][i] = current_id;
			paddeddates[ph_id][i] = current_id;


			if(st->verbose_pad)
				fprintf(logfp," tsdecided on id: %d, current_day: %e, new_next_day: %e\n", current_id, tempsaltdate[current_id][TofY_id], tempsaltdate[update_id][TofY_id]);
		} else {
			/* Padding after */
			checked = -1;
			for(k=0; k<st->ntssteps; k++) {
				if(tempsaltdate[k][TofY_id] == i) {
					last_day = tempsaltdate[k][TofY_id];
					current_id = k;
					checked++;
					break;
				} else if((fabs(tempsaltdate[k][TofY_id] - dayj) < fabs(last_day - dayj)) && (tempsaltdate[k][TofY_id] < dayj)){
					last_day = tempsaltdate[k][TofY_id];
					current_id = k;
					checked++;
					break;
				}
			}

			/* Double check any match found */
			if(checked < 0){
				k = max_id;
				last_day = tempsaltdate[k][TofY_id];
				current_id = k;
			}
			paddeddates[temp_id][i] = current_id;
			paddeddates[salt_id][i] = current_id;
			paddeddates[ph_id][i] = current_id;

			if(st->verbose_pad)
				fprintf(logfp,"tsi: %d, using id %d last_day: %e\n", i, current_id, last_day);
		}
	}

	min_tstime = min_time;
	max_tstime = max_time;
	min_tsday = min_day;
	max_tsday = max_day;
	min_tstofy = min_tofy;
	max_tstofy = max_tofy;

	/****** Padding exchange time series (identify which entry ids
	to use for the time series padding) ****/
	if(verbose)
		printf("Begin padding exchange timeseries\n");

	if(st->input_style != flat_type_id){
		min_id = 0;
		min_day = MAXINT;
		min_time = MAXINT;
		for(i=0; i<st->nhdsteps; i++) {
			if(st->verbose_pad)
				fprintf(logfp,"i: %d, min_time: %e, TofY: %e, time: %e\n", i, min_time, exchangedate[i][TofY_id], exchangedate[i][time_id]);
			if(exchangedate[i][time_id] < min_time){
				min_id = i;
				min_day = exchangedate[i][TofY_id];
				min_time = exchangedate[i][time_id];
			}
		}
		max_time = 0;
		for(i=0; i<st->nhdsteps; i++) {
			if(exchangedate[i][time_id] > max_time){
				/* Divide by dt not 86400 as may be tidal entries not daily */
				max_id = i;
				max_day = exchangedate[i][TofY_id];
				max_time = exchangedate[i][time_id];
			}
		}
	}
	min_tofy = MAXINT;
	max_tofy = -1;
	max_count = st->nhdsteps;
	for(i=0; i<st->nhdsteps; i++) {
		if(exchangedate[i][TofY_id] < min_tofy){
			min_tofy = (int)(exchangedate[i][TofY_id]);
		}
		if(exchangedate[i][TofY_id] > max_tofy){
			max_tofy = (int)(exchangedate[i][TofY_id]);
		}
	}

	if(st->input_style == flat_type_id){
		min_id = st->Exmin_id;
		min_time =exchangedate[min_id][time_id];
		min_day = 0;
		for(mth=0; mth<st->monthEx-1; mth++){
			switch(mth){
			case 0:
			case 2:
			case 4:
			case 6:
			case 7:
			case 9:
			case 11:
				min_day += 31;
				break;
			case 1:
				min_day += 28;
				break;
			case 3:
			case 5:
			case 8:
			case 10:
				min_day += 30;
				break;
			}
		}
		min_day += st->dayEx;

		max_day = 0;
		for(mth=0; mth<st->mmonthEx-1; mth++){
			switch(mth){
			case 0:
			case 2:
			case 4:
			case 6:
			case 7:
			case 9:
			case 11:
				max_day += 31;
				break;
			case 1:
				max_day += 28;
				break;
			case 3:
			case 5:
			case 8:
			case 10:
				max_day += 30;
				break;
			}
		}
		max_day += st->mdayEx;
		max_id = st->Exmax_id;
		max_time = exchangedate[max_id][time_id];
	}

	if(min_tofy > (st->tstop * (86400 / st->dt)))
		quit("No valid entries (earliest data calendar day (%d) > latest required calendar day (%d))\n", min_tofy, st->tstop);

	if(!st->recycle_flow){
		current_id = 0;
		last_day = -1;
		next_day = min_day;
		update_id = 0;
		end_yr = 0;
		yr_correct = 0;
		for(i=0; i<st->nt; i++){
			step1 = (int)(i/365);
			dayj = i - step1 * 365;
			if(i < min_day){
				/* Padding before - use appropriate time of year */
				for(k=0; k<st->nhdsteps; k++) {
					if(exchangedate[k][TofY_id] == i) {
						last_day = exchangedate[k][TofY_id];
						current_id = k;
						break;
					} else if((fabs(exchangedate[k][TofY_id] - dayj) < fabs(last_day - dayj)) &&  (exchangedate[k][TofY_id] < dayj)){
						last_day = exchangedate[k][TofY_id];
						current_id = k;
						break;
					}
				}

				/* If none found as match just use first entry */
				if(last_day < 0){
					k = min_id;
					last_day = exchangedate[k][TofY_id];
					current_id = k;
				}

				paddeddates[exchange_id][i] = current_id;

				if(st->verbose_pad)
					fprintf(logfp,"i: %d, using id %d last_day: %e\n", i, current_id, last_day);

			} else if (update_id < max_count-1){
				if(st->verbose_pad)
					fprintf(logfp,"i: %d (%d), last_day: %e, next_day: %e", i, (int)(i - end_yr_day * yr_correct), last_day, next_day);

				/* Straight read-off with interpolation */
				if(i == min_day){
					last_day = min_day;
					current_id = min_id;
					update_id = min_id + 1;
					next_day = exchangedate[update_id][TofY_id];
				} else {
					if((!end_yr && ((i - end_yr_day * yr_correct) < next_day)) || (end_yr && ((i - end_yr_day * yr_correct) < end_yr_day)));  // Nothing to update
					else {
						last_day = exchangedate[current_id][TofY_id];
						current_id++;
						update_id++;
						next_day = exchangedate[update_id][TofY_id];
						if(last_day > next_day)
							end_yr = 1;
						else{
							if(end_yr)
								yr_correct++;
							end_yr = 0;
						}
					}
				}
				paddeddates[exchange_id][i] = current_id;

				if(st->verbose_pad)
					fprintf(logfp," decided on id: %d, current_day: %e, new_next_day: %e\n", current_id, tempsaltdate[current_id][TofY_id], tempsaltdate[update_id][TofY_id]);

			} else {
				/* Padding after */
				checked = -1;
				for(k=0; k<st->nhdsteps; k++) {
					if(exchangedate[k][TofY_id] == i) {
						last_day = exchangedate[k][TofY_id];
						current_id = k;
						checked++;
						break;
					} else if((fabs(exchangedate[k][TofY_id] - dayj) < fabs(last_day - dayj)) && (exchangedate[k][TofY_id] < dayj)){
						last_day = exchangedate[k][TofY_id];
						current_id = k;
						checked++;
						break;
					}
				}

				/* Double check any match found */
				if(checked < 0){
					k = max_id;
					last_day = exchangedate[k][TofY_id];
					current_id = k;
				}
				paddeddates[exchange_id][i] = current_id;

				if(st->verbose_pad)
					fprintf(logfp,"i: %d (dayj %d), using id %d last_day: %e\n", i, dayj, current_id, last_day);
			}
		}
	} else {
		/* Just cycle the flows repeatedly and recycle them through time  */
		which_i = 0;
		for(i=0; i<st->nt; i++){
			paddeddates[exchange_id][i] = which_i;
			which_i++;
			if(which_i >= st->nhdsteps)
				which_i = 0;
		}
	}

	if(st->verbose_pad){
		fprintf(logfp,"nt: %d, ntex: %d, exchange stats: min_time: %e, min_day: %e, min_tofy: %d, max_time: %e, max_day: %e, max_tofy: %d\n",
			st->nt, st->nhdsteps, min_time, min_day, min_tofy, max_time, max_day, max_tofy);

		/* Write exchange time list */
		for(i=0; i<st->nt; i++){
			this_id = paddeddates[exchange_id][i];
			fprintf(logfp,"i-%d paddeddates[exchange_id][%d]: %d, exchangedate[%d][TofY_id]: %e, exchangedate[%d][time_id]: %e\n",
				i, i, paddeddates[exchange_id][i], this_id, exchangedate[this_id][TofY_id], this_id, exchangedate[this_id][time_id]);
		}

		fprintf(logfp,"nt: %d, ntts: %d, tempsalt stats: min_tstime: %e, min_tsday: %e, min_tstofy: %d, max_tstime: %e, max_tsday: %e, max_tstofy: %d\n",
			st->nt, st->ntssteps, min_tstime, min_tsday, min_tstofy, max_tstime, max_tsday, max_tstofy);

		/* Write temperature and salt time list */
		for(i=0; i<st->nt; i++){
			this_id = paddeddates[temp_id][i];
			fprintf(logfp,"i-%d paddeddates[temp_id][%d]: %d, paddeddates[salt_id][%d]: %d, tempsaltdate[%d][TofY_id]: %e, tempsaltdate[%d][time_id]: %e\n",
				i, i, paddeddates[temp_id][i], i, paddeddates[salt_id][i], i, tempsaltdate[this_id][TofY_id], i, tempsaltdate[this_id][time_id]);
		}

		/* Close pad log file */
		fclose(logfp);

	}

	return;

}
