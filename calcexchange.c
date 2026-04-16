/**********************************************************************

 Code for calculating exchanges per box
 by:	Beth Fulton
 date:	20/10/2003

 comments: This file is the code that converts hydrodynamic transports per face
 into transports from each box. Also pads out timeseries by reusing values
 from later in series to fill gaps at same calndar date earlier in the serious

 revisions: 28/11/2004 Updated for new raw data format (from many data files).

 18-03-2009 Bec Gorton
 Added time correction of the vertical flux.
 Change the code that checks the layer of arrival and leaving.
 leave layer check changed from:
 if(dpleave >= st->boxnz[leave_box])

 to if(dpleave > st->boxnz[leave_box])

 and the arrival layer check is now the same.
 This was forcing the flux in the top water column layer (surface) to the bottom
 water column layer.

 28-04-2009 Bec Gorton
 Added more debugging output.

 07-06-2009 Bec Gorton
 Added debugging output for the temp/salt calculations.

 14-09-2009 Bec Gorton
 Changed the HydroConstruct code to store the st->boxnz as the actual number of
 layers in a box not the max index of the layers.

 24-09-2009 Bec Gorton
 Added code to spit a warning if there is not horizontal or vertical flow
 associated with a box.

01-10-2009 Bec Gorton
Added code to check that exchange values are finite and quit if they are not.
This is because some of the matlab code is generating infinite vertical flux values.
 *************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sjwlib.h>
#include "mainexchange.h"
#include "unixUtils.h"

int maxdest = -1;

int getDestinationID(int arrive_box, int dparrive, int leave_box, int dpleave) {
	int dest_already_noted;
	int i, dest_id;

	/* Get destination ids */
	dest_already_noted = 0;
	for (i = 0; i < boxndest[leave_box]; i++) {
		if ((boxdest[i][leave_box] == arrive_box) && (layerdest[i][leave_box][dpleave] == dparrive)) {
			dest_already_noted = 1;
			dest_id = i;
			break;
		}
	}

	/* If destination box not already listed update the list */
	if (!dest_already_noted) {
		if (boxndest[leave_box] < 0)
			boxndest[leave_box] = 1;
		else
			boxndest[leave_box]++;
		dest_id = boxndest[leave_box] - 1;
		boxdest[dest_id][leave_box] = arrive_box;
		layerdest[dest_id][leave_box][dpleave] = dparrive;

		if (dest_id > maxdest)
			maxdest = dest_id;
	}

	return dest_id;

}
/***************************************************************************************************
 Main body of calculations and output of exchanges
 ****************************************************************************************************/
void ExchangeStuff(Stuff *st) {
	int actfce, i, ii, k, leave_box, arrive_box, dest_already_noted, dpleave, nzleave, fce, dp, dest_id,
			nzarrive, dparrive, tl, ta, b, dparriveSTART, dpleaveSTART, lengthcheck, altbox, filei, ndat, ns_orient, we_orient;
	double exchange_amt, units_scalar, temp_flux, salt_flux, face_length, dirscale, temp_amt, salt_amt, tot_facelength,
			step1fluxT, step1fluxS, ph_flux, step1fluxPH, ph_amt;

	int n;
	maxdest = 0;
	filei = 0;
	ndat = st->nhdsteps / st->hd.nfiles;

	if (verbose)
		printf("Begin calculating exchanges\n");

	lengthcheck = 2;
	dirscale = 1.0;

	if (st->verbose_exchange) {
		for (k = 0; k < lengthcheck; k++) {
			for (fce = 0; fce < st->nface; fce++) {
				for (dp = 0; dp < st->wcnz; dp++) {
					for (i = 0; i < 2; i++) {
						printf("k = %d, dp = %d, fce = %d, i = %d\n", k, dp, fce, i);
						printf("rawexchanges[%d][%d][%d][%d]: %e\n", k, dp, fce, i, rawexchanges[k][dp][fce][i]);
					}
				}
			}
		}
		//			for(k=st->nhdsteps-lengthcheck; k<st->nhdsteps; k++){
		//				for(fce=0; fce<st->nface; fce++){
		//					for(dp=0; dp<st->wcnz; dp++){
		//						for(i=0; i<2; i++){
		//							printf("rawexchanges[%d][%d][%d][%d]: %e\n", k, dp, fce, i, rawexchanges[k][dp][fce][i]);
		//						}
		//					}
		//				}
		//			}
	}



	/* Get water column layer volumes incase re-scaling from sverdrups required */
	Calc_Wc_Vols(st);

	/* Assign exchanges to correct destination cells and layers */
	for (k = 0; k < st->nhdsteps; k++) {
		for (fce = 0; fce < st->nface; fce++) {
			if (st->verbose_exchange > 4)
				printf("\n");
			for (dp = 0; dp < st->wcnz; dp++) {
			  leftrightcheck[fce][dp] = 0;
			  leftrightcheck[fce][st->wcnz] = 0;

			  actfce = facemap[fce];
			  /* If face unassigned (i.e. no equivalent bgm face) then skip to next face */
			  if (actfce < 0 || actfce == st->notbgmface){
			    continue;
			  }

			  if ((st->faces[actfce].ibl < 0) || (st->faces[actfce].ibr < 0)){
			  printf("\nTEST 2");
			  continue;
			}

				for (ii = 0; ii < 2; ii++) {

					/* Get raw exchange per second and multiply by total timestep to give final exchanges */
					//printf("rawexchanges[%d][%d][%d][%d] = %e\n", k, dp, fce, ii, rawexchanges[k][dp][fce][ii]);
					exchange_amt = rawexchanges[k][dp][fce][ii] * st->dt;

					/* If Al Herman style data check for winding consistency */
					if (st->input_style == NOAA_type_id) {
						filei = (int) (floor(k / ndat));

						altbox = (int) (faceinfo[filei][fce][lbox_id]);

						if (altbox == st->faces[actfce].ibl)
							dirscale = 1.0;
						else if (altbox == st->faces[actfce].ibr)
							dirscale = -1.0;
						else
							quit("Face (fce %d or actfce %d) box mismatch -> altbox: %d, bgmbox-l: %d, bgmbox-r: %d\n", fce, actfce, altbox,
									st->faces[actfce].ibl, st->faces[actfce].ibr);

						/* Correct flow from "positive if flowing into box" to
						 "positive if flowing to left across the face"
						 */
						exchange_amt *= dirscale;
					}

					/* check sign of raw exchanges, positive means flowing right-to-left across face */
					if (st->uniflow) {
						// Flows assumed always positive so adjust leave and arrive box accordingly
						if (exchange_amt < 0) {
							leave_box = st->faces[actfce].ibl;
							arrive_box = st->faces[actfce].ibr;
							exchange_amt *= -1.0;
						} else {
							leave_box = st->faces[actfce].ibr;
							arrive_box = st->faces[actfce].ibl;
						}
					} else {
						// Allow for sign of flow to signify direction
						arrive_box = st->faces[actfce].ibl;
						leave_box = st->faces[actfce].ibr;

						// Left-right check to see if getting duplicate flows
						if (exchange_amt > 0) {
							leftrightcheck[fce][dp]++;
							leftrightcheck[fce][st->wcnz]++;
						} else if (exchange_amt < 0) {
							leftrightcheck[fce][dp]--;
							leftrightcheck[fce][st->wcnz]++;
						}
					}

					dpleave = dp;
					nzleave = st->boxnz[leave_box] - 1;
					if (dpleave > nzleave)
						dpleave = 0;
					else
						dpleave = (int) (abs(dp - nzleave));

					/* If there is no exchange at this layer don't bother doing anything */
					if (exchange_amt == 0) {
						//if (st->verbose_exchange > 4)
							//printf("Exchange is nil\n");
						continue;
					}

					//if (st->verbose_exchange > 4)
					//	printf("nzleave = %d, dpleave = %d\n", nzleave, dpleave);

					if ((st->boundaries[arrive_box] > 0) && (st->boundaries[leave_box] > 0)) {
						/* If coming from and going to boundary then skip */
						if (st->verbose_exchange > 4)
							printf("Continue - boundary box\n");
						continue;

					} else if (st->boundaries[arrive_box] > 2) {
						/* For reflective boundaries reverse arrive and leave boxes */
						tl = leave_box;
						ta = arrive_box;
						leave_box = ta;
						arrive_box = tl;
					} else if (st->boundaries[leave_box] == 2) {
						/* If fenced absorptive boundary that doesn't feed into the model domain then skip */
						if (st->verbose_exchange > 4)
							printf("Continue - boundary box\n");
						continue;
					}

					/* Update exchange for that time step */
					nzarrive = st->boxnz[arrive_box] - 1;

					/* Invert water columns as flow files typically work from surface down,
					 not sediment up like in bgm structure - extra layers are all crammed
					 into the benthic surface (lower-most water column) layer */

					// BEC COMMENTED OUT - NEEDS TO BE CHECKED.
					//if (dp >= nzarrive)
					if (dp > nzarrive)
						dparrive = 0;
					else
						dparrive = (int) (abs(dp - nzarrive));

					//printf("dparrive = %d\n", dparrive);

					/* Wind back exchange to depth that is appropriate
					 for the departure as well as arrival box */
					if (((st->verbose_exchange && (st->verbose_box == -1 || dparrive == st->verbose_box || leave_box == st->verbose_box))) && dparrive
							> st->boxnz[leave_box])
						printf("dparrive (%d) > nz in dest box (%d): %d\n", leave_box, dparrive, st->boxnz[leave_box]);

//					while (dparrive >= st->boxnz[leave_box]) {
//						dparrive--;
//						if (st->verbose_exchange > 4)
//							printf("dparrive now: %d\n", dparrive);
//					}


					if(nzarrive == -1 || dparrive == -1)
						continue;

					/* Make sure exchanges in correct units (m3/s) */
					if (st->unit_type == 0)
						/* Sverdrups (10^6 m3/s) - multiply by 10^6 to get (m3/s) */
						units_scalar = 1000000.0;
					else
						/* m3/s */
						units_scalar = 1.0;

					/* Area correct flows to try and avoid hyperdiffusion */
					if (st->area_correct_flow == 1) {
						/* Simple area correction - as if square */
						units_scalar /= st->boxarea[arrive_box];
					} else if (st->area_correct_flow == 2) {
						/* x-y decomposition of flows so correct differently in x and y */
						for (i = 0; i < st->boxnconn[leave_box]; i++) {
							if (arrive_box == st->ibox[leave_box][i]) {
								/* Found correct box now see its orientation */
								ns_orient = st->boxns[leave_box][i];
								we_orient = st->boxwe[leave_box][i];


								if (ns_orient != same_lat) {
									/* North or south */
									units_scalar /= (st->boxnscoefft[arrive_box] * st->boxarea[arrive_box]);

									if (verbose >= 1)
										printf("ns: %d, abox: %d, lbox: %d, units_scalar: %e (coefft: %e, area: %e)\n", ns_orient, arrive_box, leave_box,
												units_scalar, st->boxnscoefft[arrive_box], st->boxarea[arrive_box]);

								} else {
									/* West or east */
									if (we_orient != same_long) {
										units_scalar /= (st->boxwecoefft[arrive_box] * st->boxarea[arrive_box]);

										if (verbose >= 1)
											printf("we: %d, abox: %d, lbox: %d, units_scalar: %e (coefft: %e, area: %e)\n", we_orient, arrive_box, leave_box,
													units_scalar, st->boxwecoefft[arrive_box], st->boxarea[arrive_box]);
									} else{
										/* Screw up and quit */

										quit("How can box %d not be east, west, north or south of box %d?\n", arrive_box, leave_box);
									}
								}
								break;
							}
						}
					}

					if(!_finite(units_scalar)){
						quit("scaler is not finite. st->boxnscoefft[arrive_box] = %e, st->boxarea[arrive_box] = %e\n", st->boxnscoefft[arrive_box], st->boxarea[arrive_box]);
					}


					//exchange_amt = rawexchanges[k][dp][fce][ii];
					//units_scalar = 1.0;
					exchange_amt *= units_scalar;

					/* Get destination ids */
					dest_id = getDestinationID(arrive_box, dparrive, leave_box, dpleave);
					//
					//					dest_already_noted = 0;
					//					for (i = 0; i < boxndest[leave_box]; i++) {
					//						if ((boxdest[i][leave_box] == arrive_box) && (layerdest[i][leave_box][dpleave] == dparrive)) {
					//							dest_already_noted = 1;
					//							dest_id = i;
					//							break;
					//						}
					//					}
					//
					//					/* If destination box not already listed update the list */
					//					if (!dest_already_noted) {
					//						if (boxndest[leave_box] < 0)
					//							boxndest[leave_box] = 1;
					//						else
					//							boxndest[leave_box]++;
					//						dest_id = boxndest[leave_box] - 1;
					//						boxdest[dest_id][leave_box] = arrive_box;
					//						layerdest[dest_id][leave_box][dpleave] = dparrive;
					//
					//						if (dest_id > maxdest)
					//							maxdest = dest_id;
					//					}

					/* Get final exchanges */
					finalexchanges[k][dest_id][dpleave][leave_box] += exchange_amt;

					//if(leave_box == 4 || arrive_box == 4){
						//printf("leave_box = %d, arrive_box = %d, finalexchanges[k][dest_id][dpleave][leave_box] = %e\n", leave_box, arrive_box, finalexchanges[k][dest_id][dpleave][leave_box]);
					//}

					if(finalexchanges[k][dest_id][dpleave][leave_box] > 1e15){
						printf("\n\nExchange is too large - check input data. \n\n");

						printf(
								"k:%d, face: %d, box%d-%d, t: %d, exchange: %e from leave_box: %d, layer %d, arrive_box: %d, in layer: %d (vs %d), finalexchanges: %e, vol: %e, area: %e, dpa: %d, dpl:%d, nza: %d, nzl: %d\n",
								k, fce, leave_box, dp, k, exchange_amt, leave_box, dpleave, arrive_box, dparrive, dp,
								finalexchanges[k][dest_id][dpleave][leave_box], st->boxvol[leave_box][dpleave], st->boxarea[arrive_box], dparrive, dpleave,
								st->boxnz[arrive_box], st->boxnz[leave_box]);

						printf("day: %d, dest_id: %d, exchange: %e (%.12f) from box%d-%d to box%d-%d\n", k, dest_id,
								finalexchanges[k][dest_id][dpleave][leave_box], finalexchanges[k][dest_id][dpleave][leave_box], leave_box, dpleave,
								boxdest[dest_id][leave_box], dparrive);

						quit("");
					}


					if(!_finite(finalexchanges[k][dest_id][dpleave][leave_box])){
						quit("horizontal flux - Exchange is not finite k = %d, leave_box = %d", k, leave_box);
					}


					/* Get final temperature and salinity fluxes if appropriate */
					if (st->tsflagflux == 1) {
						/* Method 1: multiply temperature or salinity flux by water flux get absolute
						 change in that property)
						 */
						temp_flux = exchange_amt * rawtempsalt[k][dpleave][leave_box][temp_id];
						salt_flux = exchange_amt * rawtempsalt[k][dpleave][leave_box][salt_id];
						ph_flux = exchange_amt * rawtempsalt[k][dpleave][leave_box][ph_id];

						tottempsaltflux[k][dpleave][leave_box][temp_id] += temp_flux;
						tottempsaltflux[k][dpleave][leave_box][salt_id] += salt_flux;
						tottempsaltflux[k][dpleave][leave_box][ph_id] += ph_flux;

					}

					if (st->verbose_exchange > 3 && (st->verbose_exchange && (st->verbose_box == -1 || arrive_box == st->verbose_box || leave_box
							== st->verbose_box))) {
						printf("rawexchanges[%d][%d][%d][%d] = %e\n", k, dp, fce, ii, rawexchanges[k][dp][fce][ii]);
						printf(
								"k:%d, face: %d, box%d-%d, t: %d, exchange: %e from leave_box: %d, layer %d, arrive_box: %d, in layer: %d (vs %d), finalexchanges: %e, vol: %e, area: %e, dpa: %d, dpl:%d, nza: %d, nzl: %d\n",
								k, fce, leave_box, dp, k, exchange_amt, leave_box, dpleave, arrive_box, dparrive, dp,
								finalexchanges[k][dest_id][dpleave][leave_box], st->boxvol[leave_box][dpleave], st->boxarea[arrive_box], dparrive, dpleave,
								st->boxnz[arrive_box], st->boxnz[leave_box]);

						printf("day: %d, dest_id: %d, exchange: %e (%.12f) from box%d-%d to box%d-%d\n", k, dest_id,
								finalexchanges[k][dest_id][dpleave][leave_box], finalexchanges[k][dest_id][dpleave][leave_box], leave_box, dpleave,
								boxdest[dest_id][leave_box], dparrive);

					}
					/* if (fce == 0){ */
					/*   printf("rawexchanges[%d][%d][%d][%d] = %e\n", k, dp, fce, ii, rawexchanges[k][dp][fce][ii]); */
					/*   printf("day: %d, dest_id: %d, exchange: %e (%.12f) from box%d-%d to box%d-%d\n", k, dest_id, */
					/* 	 finalexchanges[k][dest_id][dpleave][leave_box], finalexchanges[k][dest_id][dpleave][leave_box], leave_box, dpleave, */
					/* 	 boxdest[dest_id][leave_box], dparrive); */
					/* } */

				}
				if (st->verbose_exchange && leftrightcheck[fce][st->wcnz] && verbose > 3)
					printf("leftrightcheck[%d][%d] = %d\n", fce, dp, leftrightcheck[fce][dp]);
			}

		}
		/*
		 * The following code will subsitute flow in boxes where there is no flow with selected flow between other boxes.
		 * This was written for the ECCAL model where there were a large number of boxes that had no horizontal flow.
		 *
		 *
		 */

		for (n = 0; n < st->numSubBoxes; n++) {
			int originalArriveDp;

			for (i = 0; i < st->ndest; i++) {
				if (boxdest[i][st->subBoxInfo[n].originalLeaveBox] == st->subBoxInfo[n].originalDestBox) {

					arrive_box = st->subBoxInfo[n].targetDestBox;
					leave_box = st->subBoxInfo[n].targetLeaveBox;
					/* dp is the depth of the leaving box */
					for (dp = 0; dp < st->wcnz; dp++) {

						/* Get the original destination depth so the area correction is using the correct layer */
						originalArriveDp = layerdest[i][st->subBoxInfo[n].originalLeaveBox][dp];

						exchange_amt = finalexchanges[k][i][dp][st->subBoxInfo[n].originalLeaveBox];
						if (exchange_amt == 0)
							continue;

						/* Find the leaving depth */
						dpleave = dp;
						nzleave = st->boxnz[leave_box] - 1;
						if (dpleave > nzleave)
							dpleave = 0;

						/* Find the arrival depth */
						dparrive = dp;
						nzarrive = st->boxnz[arrive_box] - 1;
						if (dparrive > nzarrive)
							dparrive = 0;

						/* Get the new dest_id for this transport value */
						dest_id = getDestinationID(arrive_box, dparrive, leave_box, dpleave);

						/*
						 * calculated the final exchange between the new boxes.
						 * Do an area correction - divide by the original volume and mult by the volume of the new box.
						 */
						finalexchanges[k][dest_id][dpleave][leave_box] += (exchange_amt * st->boxvol[arrive_box][dparrive]
								/ st->boxvol[st->subBoxInfo[n].originalDestBox][originalArriveDp]);

						/* Check to see that there are no strange values that will result in an invalid cdf file.*/
						if (!_finite(finalexchanges[k][dest_id][dpleave][leave_box])) {
							quit("st->boxvol[arrive_box][dparrive] = %e\n", st->boxvol[arrive_box][dparrive]);
						}

						if (st->verbose_exchange > 5) {

							printf("day: %d, dest_id: %d, exchange: %e (%.12f) from box%d-%d to box%d-%d\n", k, dest_id,
									finalexchanges[k][dest_id][dpleave][leave_box], finalexchanges[k][dest_id][dpleave][leave_box], leave_box, dpleave,
									boxdest[dest_id][leave_box], dparrive);

							printf("st->boxvol[arrive_box][dparrive] = %e\n", st->boxvol[arrive_box][dparrive]);
							printf("st->boxvol[st->subBoxInfo[n].originalDestBox][originalArriveDp] = %e\n",
									st->boxvol[st->subBoxInfo[n].originalDestBox][originalArriveDp]);

						}
					}
				}
			}
		}

		/* Calculate temperatures and salinities using method 2: Sum over
		 temperature (or salinity) flux through the faces in a box,
		 weighting by the length of the face, and then divide by the water
		 exchange to get final temperature (or salinity). Note initially told the
		 other way round (commented out code), but simple dimensional check shows that
		 gives 1/deg C not deg C so reversed the calculation to get correct units.
		 */
		if (st->tsflagflux == 2) {
			/* Get total circumference of the box - already initialised so no need to do it here */
			if (!k) {
				for (b = 0; b < st->nbox; b++) {
					for (fce = 0; fce < st->nface; fce++) {
						face_length = faceinfo[filei][fce][length_id];
						if (face_length < 0)
							continue;
						altbox = (int) (faceinfo[filei][fce][lbox_id]);
						if (altbox == b) {
							for (dp = 0; dp < st->wcnz; dp++) {
								boxstats[altbox][dp][circum_id] += face_length;
							}
						}
					}
				}
			}
			/* Get temperature and salinity flux contribution */
			for (fce = 0; fce < st->nface; fce++) {
				face_length = faceinfo[filei][fce][length_id];
				if (face_length < 0)
					continue;
				altbox = (int) (faceinfo[filei][fce][lbox_id]);

				for (dp = 0; dp < st->wcnz; dp++) {
					exchange_amt = 0;
					for (ii = 0; ii < 2; ii++) {
						exchange_amt += rawexchanges[k][dp][fce][ii];
					}
					temp_amt = rawtempsalt[k][dp][fce][temp_id];
					salt_amt = rawtempsalt[k][dp][fce][salt_id];
					ph_amt = rawtempsalt[k][dp][fce][ph_id];

					tot_facelength = boxstats[altbox][dp][circum_id];

					/* Step 1: Temperature of the face = exchange_amt / temp_amt */
					step1fluxT = temp_amt / (exchange_amt + tiny);
					step1fluxS = salt_amt / (exchange_amt + tiny);
					step1fluxPH = ph_amt / (exchange_amt + tiny);

					/* Step 2: Final temperature = Sum_allfaces(FaceSurfaceArea * step1flux) / Sum_allfaces(FaceSurfaceArea)
					 where FaceSurfaceArea = face_length * st->nominal_dz[dp]
					 */

					temp_flux = (face_length * st->nominal_dz[dp]) * step1fluxT;
					salt_flux = (face_length * st->nominal_dz[dp]) * step1fluxS;
					ph_flux = (face_length * st->nominal_dz[dp]) * step1fluxPH;


					tottempsaltflux[k][dp][altbox][temp_id] += temp_flux / (tot_facelength * st->nominal_dz[dp] + tiny);
					tottempsaltflux[k][dp][altbox][salt_id] += salt_flux / (tot_facelength * st->nominal_dz[dp] + tiny);
					tottempsaltflux[k][dp][altbox][ph_id] += ph_flux / (tot_facelength * st->nominal_dz[dp] + tiny);

					if (st->verbose_exchange > 3 && (st->verbose_exchange && (st->verbose_box == -1 || altbox == st->verbose_box )) && ((k < lengthcheck) || (k > (st->nhdsteps - lengthcheck)))) {
						printf("box: %d-%d fce %d, added tempflux: %e (area: %e, exchange: %e, rawtemp flux: %e, totfcelength: %e, dz: %e) finaltemp: %e\n",
								altbox, dp, fce, temp_flux, face_length * st->nominal_dz[dp], exchange_amt, temp_amt, tot_facelength, st->nominal_dz[k],
								tottempsaltflux[k][dp][altbox][temp_id]);
						printf("box: %d-%d fce %d, added saltflux: %e (area: %e, exchange: %e, rawsalt flux: %e, totfcelength: %e, dz: %e) finalsalt: %e\n",
								altbox, dp, fce, salt_flux, face_length * st->nominal_dz[dp], exchange_amt, salt_amt, tot_facelength, st->nominal_dz[k],
								tottempsaltflux[k][dp][altbox][salt_id]);
						printf("box: %d-%d fce %d, added phflux: %e (area: %e, exchange: %e, rawph flux: %e, totfcelength: %e, dz: %e) finalPH: %e\n",
								altbox, dp, fce, ph_flux, face_length * st->nominal_dz[dp], exchange_amt, ph_amt, tot_facelength, st->nominal_dz[k],
								tottempsaltflux[k][dp][altbox][ph_id]);
					}
				}
			}

			if ((verbose > 2) && ((k < lengthcheck) || (k > (st->nhdsteps - lengthcheck)))) {
				for (fce = 0; fce < st->nbox; fce++) {
					for (dp = 0; dp < st->wcnz; dp++) {
						printf("Box: %d-%d final temp: %e, salt: %e\n", fce, dp, tottempsaltflux[k][dp][fce][temp_id], tottempsaltflux[k][dp][fce][salt_id]);
					}
				}
			}
		}

		/* Do estuary flows */
		for (fce = 0; fce < st->num_estuaries; fce++) {

			/* Face id */
			actfce = st->estnum[fce];

			/* Flow out of the estuary */
			arrive_box = st->faces[actfce].ibl;
			leave_box = st->faces[actfce].ibr;
			exchange_amt = st->estoutflux[fce];

			/* Make sure exchanges in correct units (m3/s) */
			if (!st->unit_type)
				/* Sverdrups (10^6 m3/s) - multiply by 10^6 to get (m3/s) */
				units_scalar = 1000000.0;
			else
				/* m3/s */
				units_scalar = 1.0;

			/* Area correct flows to try and avoid hyperdiffusion */
			if (st->area_correct_flow)
				units_scalar /= st->boxarea[arrive_box];

			exchange_amt *= units_scalar;

			nzarrive = st->boxnz[arrive_box] - 1;
			nzleave = st->boxnz[leave_box] - 1;

			ii = min(nzarrive,nzleave);
			for (dp = 0; dp < ii; dp++) {
				/* Get destination ids */
				dest_already_noted = 0;
				for (i = 0; i < boxndest[leave_box]; i++) {
					if ((boxdest[i][leave_box] == arrive_box) && (layerdest[i][leave_box][dp] == dp)) {
						dest_already_noted = 1;
						dest_id = i;
						break;
					}
				}

				/* If destination box not already listed update the list */
				if (!dest_already_noted) {
					if (boxndest[leave_box] < 0)
						boxndest[leave_box] = 1;
					else
						boxndest[leave_box]++;
					dest_id = boxndest[leave_box] - 1;
					boxdest[dest_id][leave_box] = arrive_box;
					layerdest[dest_id][leave_box][dp] = dp;

					if (dest_id > maxdest)
						maxdest = dest_id;
				}

				/* Get final exchanges */
				finalexchanges[k][dest_id][dp][leave_box] += exchange_amt;
			}

			/* Flow into the estuary */
			arrive_box = st->faces[actfce].ibr;
			leave_box = st->faces[actfce].ibl;
			exchange_amt = st->estinflux[fce];

			/* Make sure exchanges in correct units (m3/s) */
			if (!st->unit_type)
				/* Sverdrups (10^6 m3/s) - multiply by 10^6 to get (m3/s) */
				units_scalar = 1000000.0;
			else
				/* m3/s */
				units_scalar = 1.0;

			/* Area correct flows to try and avoid hyperdiffusion */
			if (st->area_correct_flow)
				units_scalar /= st->boxarea[arrive_box];

			exchange_amt *= units_scalar;

			nzarrive = st->boxnz[arrive_box] - 1;
			nzleave = st->boxnz[leave_box] - 1;

			ii = min(nzarrive,nzleave);
			for (dp = 0; dp < ii; dp++) {
				/* Get destination ids */
				dest_already_noted = 0;
				for (i = 0; i < boxndest[leave_box]; i++) {
					if ((boxdest[i][leave_box] == arrive_box) && (layerdest[i][leave_box][dp] == dp)) {
						dest_already_noted = 1;
						dest_id = i;
						break;
					}
				}

				/* If destination box not already listed update the list */
				if (!dest_already_noted) {
					if (boxndest[leave_box] < 0)
						boxndest[leave_box] = 1;
					else
						boxndest[leave_box]++;
					dest_id = boxndest[leave_box] - 1;
					boxdest[dest_id][leave_box] = arrive_box;
					layerdest[dest_id][leave_box][dp] = dp;

					if (dest_id > maxdest)
						maxdest = dest_id;
				}

				/* Get final exchanges */
				finalexchanges[k][dest_id][dp][leave_box] += exchange_amt;

			}
		}
	}

	/* Warn about boxes that have no flow - do this before the vertical flux
	 * values are added.
	 */

	fprintf(stderr, "\n\n");
	for (b = 0; b < st->nbox; b++) {
		int found = -1;

		if (st->boundaries[b] == 2 || st->boundaries[b] == 1) {
			/* If fenced absorptive boundary that doesn't feed into the model domain then skip */
			continue;
		}

		for (i = 0; i < st->ndest; i++) {
			if (boxdest[i][b] > 0)
				found++;
		}
		if (found == -1) {
			int destb;
			fprintf(stderr, "There is no flow from box %d\n", b);

			/* Check to see if there is any flow going 'to' this box */

			for (destb = 0; destb < st->nbox; destb++) {
				for (i = 0; i < st->ndest; i++) {
					if (boxdest[i][destb] == b) {

						//fprintf(stderr, "There is flow from box %d to %d\n", destb, b);
						found++;
					}
				}
			}

			if (found == -1) {
				fprintf(stderr, "There is no horizontal flow from or to box %d\n", b);
			}
		}

	}

	/* Vertical exchanges */
	printf("nvhdsteps: %d\n", st->nvhdsteps);

	if (st->verbose_exchange) {
		for (k = 0; k < lengthcheck; k++) {
			for (b = 0; b < st->nbox; b++) {
				for (dp = 0; dp < st->wcnz; dp++) {
					for (i = 0; i < 2; i++) {
						printf("rawvertexchanges[%d][%d][%d][%d]: %e\n", k, dp, b, i, rawvertexchanges[k][dp][b][i]);
					}
				}
			}
		}
		//		for (k = st->nvhdsteps - lengthcheck; k < st->nhdsteps; k++) {
		//			for (fce = 0; fce < st->nbox; fce++) {
		//				for (dp = 0; dp < st->wcnz; dp++) {
		//					for (i = 0; i < 2; i++) {
		//						printf("rawvertexchanges[%d][%d][%d][%d]: %e\n", k, dp, b, i, rawvertexchanges[k][dp][b][i]);
		//					}
		//				}
		//			}
		//		}

	}

	for (k = 0; k < st->nvhdsteps; k++) {
		for (b = 0; b < st->nbox; b++) {

			/* ignore boxes that are land - get into problems if we are correcting for box depth when it is 0 */
			if(st->boxtotdepth[b] < 0){
				continue;
			}

			for (dp = 0; dp < st->wcnz; dp++) {
				topbotcheck[b][dp] = 0;
				topbotcheck[b][st->wcnz] = 0;
				for (ii = 0; ii < 2; ii++) {

					/* Allow for sign of flow to signify direction */
					exchange_amt = rawvertexchanges[k][dp][b][ii] * st->dt;
					arrive_box = b;
					leave_box = b;
					dparriveSTART = dp;
					dpleaveSTART = dp + 1;

					//exchange_amt = 0;
					//printf("Box:%d:%d exchange = %e\n", b, dp, exchange_amt);
					if (!exchange_amt)
						continue;

					/* If coming up from "underneath the world", rather than allow depletion
					 of the system eventually, suck in the excess from static boundary box (as
					 ultimately in reality flow along bottom from outside model domain - given
					 all boundary boxes are static it doesn't really matter which one is feeding it
					 */
					if (dpleaveSTART >= st->wcnz) {
						dpleaveSTART = dp;
						leave_box = 0;
					}

					/* Up-down check to see if getting duplicate flows */
					if (rawvertexchanges[k][dp][b][ii] > 0) {
						topbotcheck[b][dp]++;
						topbotcheck[b][st->wcnz]++;
					} else if (rawvertexchanges[k][dp][b][ii] < 0) {
						topbotcheck[b][dp]--;
						topbotcheck[b][st->wcnz]++;
					}

					/* Update exchange for that time step */
					nzarrive = st->boxnz[arrive_box] - 1;
					nzleave = st->boxnz[leave_box] - 1;

					/* Invert water columns as flow files typically work from surface down,
					 not sediment up like in bgm structure - extra layers are all crammed
					 into the benthic surface (lower-most water column) layer */
					if (dparriveSTART >= nzarrive)
						dparrive = 0;
					else
						dparrive = (int) (abs(dparriveSTART - nzarrive));

					if (dpleaveSTART >= nzleave)
						dpleave = 0;
					else
						dpleave = (int) (abs(dpleaveSTART - nzleave));

					/* Make sure exchanges in correct unitsst->dt (m3/s) */
					if (!st->unit_type)
						/* Sverdrups (10^6 m3/s) - multiply by 10^6 to get (m3/s) */
						units_scalar = 1000000.0;
					else
						/* m3/s */
						units_scalar = 1.0;

					/* Area correct flows to try and avoid hyperdiffusion */
					if (st->area_correct_vflow) {

						/* Divide by the volume of the layer
						 */
						units_scalar /= (st->boxarea[arrive_box] * st->boxdz[arrive_box][dparrive]);
					}

					if(!_finite(exchange_amt)){
						quit("exchange_amt = %e\n", exchange_amt);
					}

					exchange_amt *= units_scalar;

					/* Get destination ids */
					dest_already_noted = 0;
					for (i = 0; i < boxndest[leave_box]; i++) {
						if ((boxdest[i][leave_box] == arrive_box) && (layerdest[i][leave_box][dpleave] == dparrive)) {
							dest_already_noted = 1;
							dest_id = i;
							break;
						}
					}

					/* If destination box not already listed update the list */
					if (!dest_already_noted) {
						if (boxndest[leave_box] < 0)
							boxndest[leave_box] = 1;
						else
							boxndest[leave_box]++;
						dest_id = boxndest[leave_box] - 1;
						boxdest[dest_id][leave_box] = arrive_box;
						layerdest[dest_id][leave_box][dpleave] = dparrive;

						if (dest_id > maxdest)
							maxdest = dest_id;
					}

					if(!_finite(exchange_amt)){
						quit(" vertical exchange - k = %d, box = %d, dparrive = %d, exchange_amt = %e, area = %e\n", k, leave_box, dparrive, exchange_amt, st->boxdz[arrive_box][dparrive]);
					}

					/* Apply final exchanges */
					finalexchanges[k][dest_id][dpleave][leave_box] += exchange_amt;


					if(!_finite(finalexchanges[k][dest_id][dpleave][leave_box])){
						printf("units_scalar = %e\n",units_scalar);
						printf("finalexchanges[k][dest_id][dpleave][leave_box] = %e\n",finalexchanges[k][dest_id][dpleave][leave_box]);

						quit("vertical flux - Exchange is not finite k = %d, leave_box = %d. arrive_box = %d, dparrive = %d, arriveBox depth = %e\n",
								k, leave_box, arrive_box, dparrive, st->boxdz[arrive_box][dparrive]);
					}

					if ((verbose > 3) &&
						(st->verbose_exchange  && (st->verbose_box == -1 || dparrive == st->verbose_box || leave_box == st->verbose_box))
					&&	st->verbose_exchange && exchange_amt)
						printf(
							"vertical exchange box%d-%d, t: %d, exchange: %e from leave_box: %d, arrive_box: %d, in layer: %d (leaving %d), finalexchanges: %e, vol: %e, area: %e, dpa: %d, dpl:%d, nza: %d, nzl: %d\n",
							leave_box, dp, k, exchange_amt, leave_box, arrive_box, dparrive, dpleave, finalexchanges[k][dest_id][dpleave][leave_box],
							st->boxvol[leave_box][dpleave], st->boxarea[arrive_box], dparrive, dpleave, st->boxnz[arrive_box] - 1, st->boxnz[leave_box]);
				}
				if ((verbose > 3) && st->verbose_exchange && topbotcheck[b][st->wcnz])
					printf("topbotcheck[%d][%d] = %d\n", b, dp, topbotcheck[b][dp]);
			}
		}
		//quit("");
	}

	/* Warn about boxes that have no flow - now do this after vertical flux is added
	 * in case there are any boxes that have no flow at all!
	 */
	fprintf(stderr, "\n\n");
	for (b = 0; b < st->nbox; b++) {
		int found = -1;

		if (st->boundaries[b] == 2 || st->boundaries[b] == 1) {
			/* If fenced absorptive boundary that doesn't feed into the model domain then skip */
			continue;

		}

		for (i = 0; i < st->ndest; i++) {
			if (boxdest[i][b] > 0)
				found++;
		}

		if (found == -1 && st->boxnz[b] > 1)
			fprintf(stderr, "There is no flow from box %d. This includes vertical and horizontal flux.\n", b);

	}
	fprintf(stderr, "\n\n");

	/* Calculate actual maximum number of destination boxes visited */
	st->fndest = -1;
	for (i = 0; i < st->nbox; i++) {
		if (boxndest[i] > st->fndest)
			st->fndest = boxndest[i];
	}

	if (st->verbose_exchange)
		printf("fndest set to %d\n", st->fndest);

	if (verbose)
		printf("maxdest: %d\n", maxdest);



	return;
}

/*************************************************************************************************************************
 Routine to calculate volume of each water column layer in each box
 **************************************************************************************************************************/
void Calc_Wc_Vols(Stuff *st) {
	int b, z, inverted_z, k;
	double tot_elapsed_dz, depth;

	/* Initialise volume array */
	st->boxvol = (double **) d_alloc2d(st->wcnz, st->nbox);

	/* Loop through watercolumn layers and calculate realised depth and multiple by box area to get layer volume */
	for (b = 0; b < st->nbox; b++) {
		tot_elapsed_dz = 0;

		depth = 0;
		k = 0;

		//printf("b = %d, st->boxtotdepth[b] = %e\n", b, st->boxtotdepth[b]);

		/* Loop over the nominal_dz values until we reach the surface */
		while (depth < st->boxtotdepth[b] && k < st->wcnz) {
			//printf("depth = %e,  st->nominal_dz[%d] = %e\n", depth,  k, st->nominal_dz[k]);
			depth += st->nominal_dz[k];
			k++;
		}

		if (k != st->boxnz[b])
			warn("The numLayers value for box %d (%e) is %d but the number of layers calculated based on nominal_dz is %d - check and retry.\n", b,
					st->boxtotdepth[b], st->boxnz[b], k);

		st->boxnz[b] = k;


		for (z = 0; z < st->boxnz[b]; z++) {
			/* Invert z so working sediment up not surface down*/
			inverted_z = (int) (abs(st->boxnz[b] - z - 1));

			tot_elapsed_dz += st->nominal_dz[z];

			if (tot_elapsed_dz > st->boxtotdepth[b])
				depth = st->boxtotdepth[b] - (tot_elapsed_dz - st->nominal_dz[z]);
			else
				depth = st->nominal_dz[z];

			if (depth <= 0) { /* Flows off land */
				depth = 0.01;
				warn("According to bathymetry box %d is land (depth = %e)\n", b, depth);
			}

			st->boxvol[b][inverted_z] = depth * st->boxarea[b];

			st->boxdz[b][z] = depth;

			if (verbose > 3)
				printf("box: %d, hydrolayer: %d, bgmlayer: %d, depth: %e, area: %e\n", b, z, inverted_z, depth, st->boxarea[b]);
		}

		//if (tot_elapsed_dz != st->boxtotdepth[b])
			//printf("bpx = %d, numLayers = %d\n", b, st->boxnz[b]);

		/* Set all other layers to volume zero - boxnz already corrected back to zero
		 so must "uncorrect" here */
		for (z = st->boxnz[b]; z < st->wcnz; z++) {
			/* Invert z so working sediment up not surface down */
			st->boxvol[b][z] = 0.0;

			if (verbose > 3)
				printf("box: %d, layer: %d, vol set to zero\n", b, z);
		}

		if (verbose > 3) {
			for (z = 0; z < st->wcnz; z++)
				printf("box: %d, layer: %d, volume: %e\n", b, z, st->boxvol[b][z]);
		}

	}
	return;
}

/***************************************************************************************************
 Calculations of temperatures and salinities
 ****************************************************************************************************/
void CalcTempSalt(Stuff *st) {
	int j, b, i, k, bact, botdp, bi, ti, ij, itcount, test_needed, maxnum;

	if (!st->tsflagflux) {
		/* Absolute values read in so no extra calculations required - just need to reverse layer order */
		for (j = 0; j < st->ntssteps; j++) {
			for (i = 0; i < st->nbox; i++) {

				/* deal with land boxes */
				if(st->boxtotdepth[i] < 0){
					for (k = 0; k < num_TSP_values; k++) {
						finaltempsalt[j][0][i][k] = rawtempsalt[j][0][i][k];
					}
					//finaltempsalt[j][0][i][swr_id] = rawtempsalt[j][0][i][swr_id];

				}else{


					for (b = 0; b < st->boxnz[i]; b++) {
						bact = st->boxnz[i] - 1 - b;
						for (k = 0; k < num_TSP_values; k++) {
							finaltempsalt[j][bact][i][k] = rawtempsalt[j][b][i][k];
						}
						finaltempsalt[j][0][i][swr_id] = rawtempsalt[j][0][i][swr_id];
						if (verbose >= 1)
							printf("step: %d box%d-%d temp: %e salt: %e, swr = %e\n", j, i, bact, finaltempsalt[j][bact][i][temp_id], finaltempsalt[j][bact][i][salt_id], finaltempsalt[j][bact][i][swr_id]);
					}
					/* For the "ghost layers" set it to bottom values */
					botdp = st->boxnz[i];
					for (b = st->boxnz[i]; b < st->wcnz; b++) {
						for (k = 0; k < num_TSP_values; k++) {
							finaltempsalt[j][b][i][k] = finaltempsalt[j][botdp][i][k];
						}
						if (verbose >= 1)
							printf("ghost step: %d box%d-%d temp: %e salt: %e, swr = %e\n", j, i, b, finaltempsalt[j][b][i][temp_id], finaltempsalt[j][b][i][salt_id], finaltempsalt[j][b][i][swr_id]);
							//printf("ghost step: %d box%d-%d temp: %e salt: %e\n", j, i, b, finaltempsalt[j][b][i][temp_id], finaltempsalt[j][b][i][salt_id]);
					}
				}
			}
		}
	} else if (st->tsflagflux == 1) {
		/* Fluxes read in so need to calculate final values */
		for (k = 0; k < st->ntssteps; k++) {
			for (b = 0; b < st->wcnz; b++) {
				for (i = 0; i < st->nbox; i++) {
					/* If first step start from initial conditions */
					if (!k) {
						finaltempsalt[k][b][i][temp_id] = TempSaltInit[i][b][temp_id] + tottempsaltflux[k][b][i][temp_id];
						finaltempsalt[k][b][i][salt_id] = TempSaltInit[i][b][salt_id] + tottempsaltflux[k][b][i][salt_id];
						finaltempsalt[k][b][i][ph_id] = TempSaltInit[i][b][ph_id] + tottempsaltflux[k][b][i][ph_id];
						if(st->doSWR == TRUE){
							finaltempsalt[k][b][i][swr_id] = TempSaltInit[i][b][ph_id] + tottempsaltflux[k][b][i][swr_id];
						}

					} else {
						/* Update last value */
						finaltempsalt[k][b][i][temp_id] = finaltempsalt[k - 1][b][i][temp_id] + tottempsaltflux[k][b][i][temp_id];
						finaltempsalt[k][b][i][salt_id] = finaltempsalt[k - 1][b][i][salt_id] + tottempsaltflux[k][b][i][salt_id];
						finaltempsalt[k][b][i][ph_id] = finaltempsalt[k - 1][b][i][ph_id] + tottempsaltflux[k][b][i][ph_id];

						if(st->doSWR == TRUE){
							finaltempsalt[k][b][i][swr_id] = finaltempsalt[k - 1][b][i][swr_id] + tottempsaltflux[k][b][i][swr_id];
						}

					}

					if (verbose > 3)
						printf("step: %d box%d-%d temp: %e salt: %e\n", k, i, b, finaltempsalt[k][b][i][temp_id], finaltempsalt[k][b][i][salt_id]);
				}
			}
		}
	} else if (st->tsflagflux == 2) {
		/* Fluxes read in so need to calculate final values */
		for (k = 0; k < st->ntssteps; k++) {
			for (b = 0; b < st->wcnz; b++) {
				for (i = 0; i < st->nbox; i++) {
					/* Sum over faces to get final value for the box - sum should already
					 be done as did += in ExchangeStuff calc. Thus only simple assignment
					 to the final array remains to be done
					 */
					finaltempsalt[k][b][i][temp_id] = tottempsaltflux[k][b][i][temp_id];
					finaltempsalt[k][b][i][salt_id] = tottempsaltflux[k][b][i][salt_id];
					finaltempsalt[k][b][i][ph_id] = tottempsaltflux[k][b][i][ph_id];

					if (verbose > 3)
						printf("step: %d box%d-%d temp: %e salt: %e, pH: %e\n", k, i, b, finaltempsalt[k][b][i][temp_id],
								finaltempsalt[k][b][i][salt_id],finaltempsalt[k][b][i][ph_id]);
				}
			}
		}
	} else {
		/* No temperature or salinity data being read-in or stored in this case */
	}

	/* Do one final set of checks

	 At data set 1, 10, 100, 1000 - if they exist - check temperature decreases with depth
	 */
	test_needed = 1;
	maxnum = 0;
	for (i = 0; i < st->nbox; i++) {
		/* Just test the first box with maximum number of layers */
		if (maxnum < st->boxnz[i])
			maxnum = st->boxnz[i] - 1;
		if (test_needed && (st->boxnz[i] == st->wcnz)) {
			itcount = 0;
			ij = 1;
			while (ij < st->ntssteps) {
				for (b = 0; b < st->boxnz[i] - 1; b++) {
					if (finaltempsalt[ij][b][i][temp_id] > finaltempsalt[ij][b + 1][i][temp_id])
						test_ok[itcount] = -i;
				}
				itcount++;
				ij = next_test[itcount];
			}
			test_needed = 0;
			break;
		}
	}

	/* Check if test happened */
	if (test_needed) {
		quit("Didn't find a box to test temperatures in - wcnz = %d, max_num layers found = %d\n", st->wcnz, maxnum);
	}

	/* Otherwise output results of test */
	for (ij = 0; ij < st->nitcount; ij++) {
		if (test_ok[ij] < 1) {
			bi = (-1 * test_ok[ij]);
			printf("For box %d temperatures are (from bottom up)", bi);
			for (b = 0; b < st->boxnz[i]; b++) {
				ti = next_test[ij];
				printf("%e ", finaltempsalt[ti][b][bi][temp_id]);
			}
			printf("\n");
		}
	}

	return;
}
