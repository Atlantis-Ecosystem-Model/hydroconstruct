/*
 * unixUtils.c
 *
 *
 * Stuff possibly needed for conversion to unix.
 *
 *  Created on: Feb 19, 2013
 *      Author: bec
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

int _isnan( double v ) { return isnan( v ); }
int _finite( double v ) { return isfinite( v ); }

//#define _finite isfinite




