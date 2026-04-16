/**********************************************************************


 *************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <math.h>
#include <string.h>
#include <sjwlib.h>
#include "mainexchange.h"

void checkTransport(Stuff *st, int indexed_countend, int indexed_start, int numid) {
	FPTYPE **startVolume;
	FPTYPE **newVolume;
	FPTYPE exchange;
	int box, layer, dest;
	int dayj, this_id, destBox, dest_k;
	FPTYPE percentChange;
	FILE *fid;
	int countend = (int) (st->tstop * (86400 / st->dt));
	char fname[100];
	int flag = 0;
	startVolume = (double **) d_alloc2d(st->wcnz, st->nbox);
	newVolume = (double **) d_alloc2d(st->wcnz, st->nbox);

	warn("Checking volume \n");
	// -1 as will stop in output file writing because of < in for statement

	sprintf(fname, "volume%d.cdf", numid);
	if (verbose)
		printf("Create volume output file\n");

	/* Create file */
	if ((fid = fopen(fname, "w")) == NULL)
		quit("initOutFile: Can't open %s\n", fname);

	/* File header info */
	fprintf(fid, "netcdf temp_series { \ndimensions:\n	t = UNLIMITED ; // (%d currently)\n", st->num_entries + 1);
	fprintf(fid, "	b = %d ;\n	z = %d ;\n", st->nbox, st->wcnz + 1);
	fprintf(fid, "variables:\n");
	fprintf(fid, "	double t(t) ;\n		t:units = \"seconds since %d-01-01 00:00:00 +10\" ;\n		t:dt = %d. ;", st->REFyear, (int) (st->dt));
	fprintf(fid, "\n	double volume(t, b, z) ;\n		volume:_FillValue = 0. ;");
	fprintf(fid, "\n	double volume_percent(b, z) ;\n		volume_percent:_FillValue = 0. ;");
	fprintf(fid, "\n// global attributes:\n		:title = \"trivial\";\n		:geometry = \"%s\" ;", st->geomIfname);
	fprintf(fid, "\n		:parameters = \"\" ;\ndata:\n\n t = ");

	/* Write time stamps */
	writeTimeSteps(st, fid, st->tstart, countend);

	for (box = 0; box < st->nbox; box++) {
		for (layer = 0; layer < st->wcnz; layer++) {
			if (layer < st->boxnz[box]) {
				startVolume[box][layer] = st->boxvol[box][layer];
				newVolume[box][layer] = 0;
			} else {
				startVolume[box][layer] = 0;
				newVolume[box][layer] = 0;
			}
		}
	}

	fprintf(fid, ";\n\n volume =\n  ");
	for (dayj = indexed_start; dayj < indexed_countend; dayj++) {
		this_id = paddeddates[exchange_id][dayj];

		for (box = 0; box < st->nbox; box++) {
			for (layer = 0; layer < st->wcnz; layer++) {
				if (verbose > 3)
					printf("boxndest[box%d]: %d, ndest: %d\n", box, boxndest[box], st->fndest);

				for (dest = 0; dest < boxndest[box]; dest++) {
					/* Get the destination box */
					destBox = boxdest[dest][box];

					/* Get the destination layer */
					dest_k = layerdest[dest][box][layer];

					if (!isnan(finalexchanges[this_id][dest][layer][box])) {
						exchange = finalexchanges[this_id][dest][layer][box];

						newVolume[box][layer] += -exchange;
						newVolume[destBox][dest_k] +=  exchange;
					}

				}
				if (flag == 1)
					fprintf(fid, ", %.12f", newVolume[box][layer] + startVolume[box][layer]);
				else
					fprintf(fid, "%.12f", newVolume[box][layer] + startVolume[box][layer]);
				flag = 1;

			}
			fprintf(fid, ",0 ");
			fprintf(fid, "\n  ");
		}
	}

	fprintf(fid, ";\n\n volume_percent =\n  ");
	flag = 0;
	for (box = 0; box < st->nbox; box++) {
		for (layer = 0; layer < st->wcnz; layer++) {
			percentChange = 0.0;

			if (layer < st->boxnz[box]) {
				percentChange = 100 * newVolume[box][layer] / startVolume[box][layer];

				if (fabs(percentChange) == 0 || fabs(percentChange) > 0.01)
					printf("Box %d, layer %d percent change %e, startVolume = %e, newVolume = %e\n",
							box, layer, percentChange, startVolume[box][layer], newVolume[box][layer]);
			}
			if (flag == 1)
				fprintf(fid, ", %.12f", percentChange);
			else
				fprintf(fid, "%.12f", percentChange);
			flag = 1;
		}
		fprintf(fid, ",0 ");
		fprintf(fid, "\n  ");
	}

	fprintf(fid, ";\n\n}");
	d_free2d(startVolume);
	d_free2d(newVolume);


	fclose(fid);
}
