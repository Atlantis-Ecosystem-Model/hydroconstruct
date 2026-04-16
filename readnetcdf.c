/**********************************************************************

 File reads in model geometry
 by:	Beth Fulton
 date:	21/10/2003

 comments: This file is the code that reads in the model geometry

 Changes:

 14-09-2009 Bec Gorton
 Changed the code to work out which face each vertices belongs to - this is
 required for the new ROMS data. The new ROMS input file has data for all faces - regardless
 of if they are into the model or out of the model.


 01-10-2009 Bec Gorton
 Added code to free up the geometry structure.

 *************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <malloc.h>
#include <sjwlib.h>
#include <string.h>
#include "mainexchange.h"

/**//**
 /brief Free the memory associated with the geometry.
 ****/
void freeGeom(Stuff *st) {

	if (st->nbgmface > 0) {
		free(st->faces);
		i_free1d(st->faceVertIndex);
	}

	/* Get inside points */
	free(st->boxinsidey);
	free(st->boxinsidex);

	free(st->boxnconn);

	/* Allocate memory for face and neighbour box lists */
	i_free2d(st->ibox);
	i_free2d(st->iface);
	i_free2d(st->boxns);
	i_free2d(st->boxwe);

	/** Get box areas and layer volumes **/
	/* Create area and depth arrays */
	if (st->nbox > 0) {
		free(st->boxarea);
		free(st->boxtotdepth);
	}

	/* Get box heights and widths */
	free(st->boxht);
	free(st->boxwd);
	free(st->boxnscoefft);
	free(st->boxwecoefft);

}

/************************************************************************************************
 code is of SJW's readBoxModelGeom in geomIO.c
 This routine reads a box model geometry from an ascii file. It does not return if an error is encountered.
 Input file format for this part is shown below, where N is an integer and X and Y are floating point
 numbers.

 # Number of boxes in horizontal plane
 nbox        N

 # Number of faces in horizontal plane
 nface       N

 # Box Model polyline (should be closed)
 bnd_vert    X    Y
 bnd_vert    X    Y
 bnd_vert    X    Y
 .
 .
 .

 The boxes and faces are read using readBoxGeom() and
 readFaceGeom() below. This routine assumes that the
 boxes and faces are each in order in the input file
 (ie, data for box0 comes before data for box1 etc). This
 minimises the number of times which the file has to be rewound
 and re-read.

 Arguments:     name   input file name
 bm      box model pointer
 ***************************************************************************************************/
void readBoxGeom(char *name, Stuff *st) {
	char key[STSLEN];
	FILE *fp;
	int i, nconn, nndest, b, faceID, vertIndex;
	double barea, botzz, maxx, maxy, minx, miny,
			px, py;
	int buflen = 2000;
	char buf[2000];
	double prevPx = 0.0, prevPy = 0.0;

	printf("Reading in bgm file %s\n", name);
	if (verbose)
		printf("Read basic box geometry\n");

	/* open the file */
	if ((fp = fopen(name, "r")) == NULL)
		quit("readBoxGeom: Can't open %s\n", name);

	/* rewind file */
	fseek(fp, 0L, SEEK_SET);

	/* Read number of boxes and faces */
	skipToKeyEnd(fp, "nbox");
	if (fscanf(fp, "%d", &st->nbox) != 1 || st->nbox < 1 || st->nbox > 10000)
		quit("readBoxModelGeom: Can't read nbox\n");
	skipToKeyEnd(fp, "nface");
	if (fscanf(fp, "%d", &st->nbgmface) != 1 || st->nbgmface < 0
			|| st->nbgmface > 20000)
		quit("readBoxGeom: Can't read nface\n");
	else
		printf("Found %d faces\n", st->nbgmface);
	st->nface = st->nbgmface;

	/* allocate face arrays */
	if (st->nbgmface > 0) {
		st->faces = (Face *) malloc(st->nbgmface * sizeof(Face));
		st->faceVertIndex = i_alloc1d(st->nbgmface);
		if (st->faces == NULL)
			quit("readBoxGeom: No memory for face array\n");
	} else
		st->faces = NULL;

	/* Set up faces */
	if (verbose)
		printf("Read face definitions\n");

	fseek(fp, 0L, SEEK_SET);
	for (i = 0; i < st->nbgmface; i++) {
		st->faces[i].n = i;
		st->faceVertIndex[i] = -1;
		readFaceGeom(fp, &st->faces[i]);
	}

	/* Get inside points */
	st->boxinsidey = (double *) malloc(st->nbox * sizeof(double));
	st->boxinsidex = (double *) malloc(st->nbox * sizeof(double));

	/* Read station location */
	for (i = 0; i < st->nbox; i++) {
		sprintf(key, "box%d.inside", i);
		skipToKeyEnd(fp, key);
		if (fscanf(fp, "%lf %lf", &st->boxinsidex[i], &st->boxinsidey[i]) != 2){
			warn("readBoxGeom: Can't read %s\n", key);
			st->boxinsidex[i] = 0;
			st->boxinsidey[i] = 0;
		}
	}

	/* Determine ndest */
	st->boxnconn = (int *) malloc(st->nbox * sizeof(int));
	nndest = 0;
	for (i = 0; i < st->nbox; i++) {
		/* Read number of connections */
		sprintf(key, "box%d.nconn", i);
		skipToKeyEnd(fp, key);
		if (fscanf(fp, "%d", &nconn) != 1 || nconn < 0)
			quit("readBoxGeom: Can't read %s\n", key);
		st->boxnconn[i] = nconn;

		/* Compare against ndest so that save the largest value of connections as the
		 number of destinations to worry about*/
		if (nconn > nndest)
			nndest = nconn;
	}
	if ((nndest > st->ndest) || (!st->ndest))
		st->ndest = nndest;

	/* Allocate memory for face and neighbour box lists */
	st->ibox = (int **) i_alloc2d(nndest, st->nbox);
	st->iface = (int **) i_alloc2d(nndest, st->nbox);
	st->boxns = (int **) i_alloc2d(nndest, st->nbox);
	st->boxwe = (int **) i_alloc2d(nndest, st->nbox);

	for (b = 0; b < st->nbox; b++) {
		for (i = 0; i < nndest; i++) {
			/* Initialise box index list */
			st->ibox[b][i] = -1;

			/* Initialise boxns and boxwe arrays */
			st->boxns[b][i] = -1;
			st->boxwe[b][i] = -1;
		}
	}

	/*Read the face index list */
	for (b = 0; b < st->nbox; b++) {
		sprintf(key, "box%d.iface", b);
		skipToKeyEnd(fp, key);
		for (i = 0; i < st->boxnconn[b]; i++)
			if (fscanf(fp, "%d", &st->iface[b][i]) != 1)
				quit("readBoxGeom: Can't read %s\n", key);
	}

	for (b = 0; b < st->nbox; b++) {
		/* Read the box index list */
		sprintf(key, "box%d.ibox", b);
		skipToKeyEnd(fp, key);
		for (i = 0; i < st->boxnconn[b]; i++){
			if (fscanf(fp, "%d", &st->ibox[b][i]) != 1)
				quit("readBoxGeom: Can't read %s\n", key);
		}
	}

	/** Get box areas and layer volumes **/
	/* Create area and depth arrays */
	if (st->nbox > 0) {
		st->boxarea = (double *) malloc(st->nbox * sizeof(double));
		if (st->boxarea == NULL)
			quit("readBoxGeom: No memory for box area\n");
	} else
		st->boxarea = NULL;
	if (st->nbox > 0) {
		st->boxtotdepth = (double *) malloc(st->nbox * sizeof(double));
		if (st->boxtotdepth == NULL)
			quit("readBoxGeom: No memory for total depth of the box\n");
	} else
		st->boxtotdepth = NULL;

	for (i = 0; i < st->nbox; i++) {
		/* Read Area */
		sprintf(key, "box%d.area", i);
		skipToKeyEnd(fp, key);
		if (fscanf(fp, "%lf", &barea) != 1 || barea < 0)
			quit("readBoxGeom: Can't read valid area for %s (area = %d)\n",
					key, barea);
		st->boxarea[i] = barea;

		/* Read tot water column volume in the cell */
		sprintf(key, "box%d.botz", i);
		skipToKeyEnd(fp, key);
		if (fscanf(fp, "%lf", &botzz) != 1)
			quit(
					"readBoxGeom: Can't read valid total depth for %s (depth = %d)\n",
					key, botzz);
		st->boxtotdepth[i] = -botzz;
	}

	/* Get box heights and widths */
	st->boxht = (double *) malloc(st->nbox * sizeof(double));
	st->boxwd = (double *) malloc(st->nbox * sizeof(double));
	st->boxnscoefft = (double *) malloc(st->nbox * sizeof(double));
	st->boxwecoefft = (double *) malloc(st->nbox * sizeof(double));
	for (i = 0; i < st->nbox; i++) {
		//printf("box = %d\n", i);
		maxx = -MAXDOUBLE;
		maxy = -MAXDOUBLE;
		minx = MAXDOUBLE;
		miny = MAXDOUBLE;
		sprintf(key, "box%d.vert ", i);
		skipToKeyStart(fp, key);
		sprintf(key, "box%d.vert %%lf %%lf ", i);
		// sprintf(key,"%%lf %%lf ");
		vertIndex = 0;
		//printf("fscanf(fp,key,&px,&py) = %d\n", fscanf(fp,key,&px,&py));

		prevPx = -1;
		while (fgets(buf, buflen, fp) != NULL) {

			//	printf("buf = %s\n", buf);
			if (strstr(buf, "# Data") != NULL)
				break;

			if (sscanf(buf, key, &px, &py) != 2)
				break;

			//while( fscanf(fp,key,&px,&py) == 2 ){
			if (py > maxy)
				maxy = py;
			if (py < miny)
				miny = py;
			if (px > maxx)
				maxx = px;
			if (px < minx)
				minx = px;

			/* Find which face this corresponds to */
			if (vertIndex > 0) {

				//printf("vertIndex = %d, px = %e, py = %e, prevPx = %e, prevPy = %e\n", vertIndex - 1, px, py, prevPx, prevPy);
				for (faceID = 0; faceID < st->nbgmface; faceID++) {
					//printf("st->faces[%d].p2.x = %e, st->faces[%d].p2.y = %e\n", faceID, st->faces[faceID].p2.x, faceID, st->faces[faceID].p2.y);

					if ((st->faces[faceID].p2.x == px && st->faces[faceID].p2.y
							== py && st->faces[faceID].p1.x == prevPx
							&& st->faces[faceID].p1.y == prevPy)
							|| (st->faces[faceID].p1.x == px
									&& st->faces[faceID].p1.y == py
									&& st->faces[faceID].p2.x == prevPx
									&& st->faces[faceID].p2.y == prevPy)) {

						/* Check to see that this face is to do with this box */
						if (st->faces[faceID].ibl == i || st->faces[faceID].ibr
								== i) {
							st->faceVertIndex[faceID] = vertIndex - 1;
							//printf("FOUND!!. faceID = %d, vertIndex = %d\n", faceID, vertIndex - 1);
							break;
						}
					}
				}
			}

			prevPx = px;
			prevPy = py;
			vertIndex++;
		}

		/* Assign width and height of the box - note that for diagonal boxes these will be overstated */
		st->boxht[i] = maxy - miny;
		st->boxwd[i] = maxx - minx;

		/* Initialise coefficients */
		st->boxnscoefft[i] = 1;
		st->boxwecoefft[i] = 1;

		/* Rescale coefficents based on shape */
		if (st->boxht[i] < st->boxwd[i]) // Short and fat
			st->boxwecoefft[i] = st->boxht[i] / (st->boxwd[i] + tiny);
		if (st->boxwd[i] < st->boxht[i]) // Tall and thin
			st->boxnscoefft[i] = st->boxwd[i] / (st->boxht[i] + tiny);

		/* Check for null or negative sized boxes */
		if (st->boxnscoefft[i] <= 0) {
			warn("Had to reset box %d height\n", i);
			st->boxnscoefft[i] = tiny;
		}
		if (st->boxwecoefft[i] <= 0) {
			warn("Had to reset box %d width\n", i);
			st->boxwecoefft[i] = tiny;
		}

	}
	/* Close the file */
	fclose(fp);

	if (verbose)
		printf("Geometry file closed\n");
}

/*************************************************************************************************
 Code is of SJW's readFaceGeom in geomIO.c
 This routine reads a face geometry from an ascii file. It assumes that the face already contains
 valid values for n. It also assumes that the file is positioned at or before the data for the face
 requested. No file positioning is done here, apart from reading data. This means that my parameter
 library (in sjwlib) should not be used, as its routines do explicit file positioning.
 Input file format for data for 1 face is shown below. Here X, Y, L, C and S are all floating point values.

 # Geometry for face nn
 # Start point of face
 facenn.p1	X    Y

 # End point of face
 facenn.p2	X    Y

 # Length of face
 facenn.length   L

 # Cos and Sin of angle from +ve x axis to normal to face
 facenn.cs   C    S

 # Indices of boxes to left and right of this face
 facenn.lr   N    N


 Arguments:     fp    input file pointer
 f     face pointer
 ***************************************************************************************************/
void readFaceGeom(FILE *fp, Face *f) {
	char key[STSLEN];

	/* Sanity checks */
	if (f == NULL)
		quit("readFaceGeom: NULL face pointer!\n");
	if (f->n < 0)
		quit("readFaceGeom: Face number < 0!\n");

	/* Read face endpoints */
	sprintf(key, "face%d.p1", f->n);
	skipToKeyEnd(fp, key);
	if (fscanf(fp, "%lf %lf", &f->p1.x, &f->p1.y) != 2)
		quit("readFaceGeom: Can't read %s\n", key);
	sprintf(key, "face%d.p2", f->n);
	skipToKeyEnd(fp, key);
	if (fscanf(fp, "%lf %lf", &f->p2.x, &f->p2.y) != 2)
		quit("readFaceGeom: Can't read %s\n", key);

	/* Read face length */
	sprintf(key, "face%d.length", f->n);
	skipToKeyEnd(fp, key);
	if (fscanf(fp, "%lf", &f->len) != 1 || f->len <= 0.0)
		quit("readFaceGeom: Can't read %s\n", key);

	/* Read face orientation */
	sprintf(key, "face%d.cs", f->n);
	skipToKeyEnd(fp, key);
	if (fscanf(fp, "%lf %lf", &f->cos, &f->sin) != 2)
		quit("readFaceGeom: Can't read %s\n", key);

	/* Read indices of boxes to left and right */
	sprintf(key, "face%d.lr", f->n);
	skipToKeyEnd(fp, key);
	if (fscanf(fp, "%d %d", &f->ibl, &f->ibr) != 2 || f->ibl < 0 || f->ibr < 0)
		quit("readFaceGeom: Can't read %s\n", key);
}

/************************************************************************************************
 This routine reads only the face info in lat-long to supplement existing info in x-y from basic
 readBoxGeom() routine -

 This routine reads in the face geometry in lat-long coords from an ascii file (to supplement
 x-y coordinate info from basic bgm file). This is so that code can check Jeff Dunn style
 netcdf input files are mapped to correct box and face ID numbers in bgm file

 It assumes that the face already contains valid values for n. It also assumes that the
 file is positioned at or before the data for the face requested. No file positioning is
 done here, apart from reading data. Input file format for data for 1 face is shown below.

 # Geometry for face nn
 # Start point of face
 facenn.p1	X    Y

 # End point of face
 facenn.p2	X    Y

 Arguments:     name   input file name
 bm      box model pointer
 ***************************************************************************************************/
void readBoxGeomll(char *name, Stuff *st) {
	char key[STSLEN];
	FILE *fp;
	int i;

	if (verbose)
		printf("Read lat-long face geometry\n");

	/* open the file */
	if ((fp = fopen(name, "r")) == NULL)
		quit("readBoxGeomll: Can't open %s\n", name);

	/* rewind file */
	fseek(fp, 0L, SEEK_SET);

	for (i = 0; i < st->nbgmface; i++) {
		/* Sanity checks */
		if (st->faces[i].n < 0)
			quit("readllFaceGeom: Face number < 0!\n");

		/* Read face endpoints in latitude and longitude */
		sprintf(key, "face%d.p1", st->faces[i].n);

		skipToKeyEnd(fp, key);
		if (fscanf(fp, "%lf %lf", &st->faces[i].pll1.x, &st->faces[i].pll1.y)
				!= 2)
			quit("read11FaceGeom: Can't read %s in latlong faces file\n", key);
		sprintf(key, "face%d.p2", st->faces[i].n);
		skipToKeyEnd(fp, key);
		if (fscanf(fp, "%lf %lf", &st->faces[i].pll2.x, &st->faces[i].pll2.y)
				!= 2)
			quit("read11FaceGeom: Can't read %s in latlong faces file\n", key);

		if (verbose > 1)
			printf("bgmfcell: %d, p1(%e, %e) p2(%e, %e)\n", i,
					st->faces[i].pll1.x, st->faces[i].pll1.y,
					st->faces[i].pll2.x, st->faces[i].pll2.y);
	}

	/* Close the file */
	fclose(fp);

	if (verbose)
		printf("Geometry file closed\n");

	return;
}

/**************************************************************************************
 Routine to identify which surrounding boxes are north-south and east-west for each box,
 so can do u-v decomposition of currents and scale to account for hyperdiffusion
 differently based on whether flow is in x or y - so can cope with long thin boxes, rather
 than almost square boxes (i.e. wide as tall).
 **************************************************************************************/
void setup_NorthSouth(Stuff *st) {
	int b, ij, checkbox;
	double currentx, currenty, checkx, checky;

	for (b = 0; b < st->nbox; b++) {
		currentx = st->boxinsidex[b];
		currenty = st->boxinsidey[b];
		for (ij = 0; ij < st->boxnconn[b]; ij++) {
			checkbox = st->ibox[b][ij];
			if (checkbox >= 0) {
				/* Only check valid boxes */
				checkx = st->boxinsidex[checkbox];
				checky = st->boxinsidey[checkbox];

				/* Check if north or south */
				if ((checky - currenty) > (st->boxht[b] * 0.5))
					st->boxns[b][ij] = to_north; // Box to north
				else if ((currenty - checky) > (st->boxht[b] * 0.5))
					st->boxns[b][ij] = to_south; // Box to south
				else
					st->boxns[b][ij] = same_lat; // In same latitudinal band

				/* Check of east or west */
				if (checkx < currentx)
					st->boxwe[b][ij] = to_west; // Box to west
				else if (checkx > currentx)
					st->boxwe[b][ij] = to_east; // Box to east
				else
					st->boxwe[b][ij] = same_long; // In same longitudinal band

			}
		}
	}
	return;
}
