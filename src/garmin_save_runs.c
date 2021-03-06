/*
  Garmintools software package
  Copyright (C) 2006-2008 Dave Bailey
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include "config.h"
#include <stdio.h>
#include <unistd.h>
#include "garmin.h"


int
main ( int argc, char ** argv )
{
  garmin_unit garmin;
  int         verbose;

  /* Set the verbosity if the -v option was provided. */

  verbose = (getopt(argc,argv,"v") != -1);

  if ( garmin_init(&garmin,verbose) != 0 ) {
    /* Read and save the runs. */
    garmin_save_runs(&garmin);
  } else {
    printf("garmin unit could not be opened!\n");
  }

  return 0;
}
