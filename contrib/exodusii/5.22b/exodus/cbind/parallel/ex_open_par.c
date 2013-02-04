/*
 * Copyright (c) 2005 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Governement
 * retains certain rights in this software.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 * 
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 * 
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.  
 * 
 *     * Neither the name of Sandia Corporation nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 */
/*****************************************************************************
*
* exopen - ex_open
*
* entry conditions - 
*   input parameters:
*       char*   path                    exodus filename path
*       int     mode                    access mode w/r
*
* exit conditions - 
*       int     exoid                   exodus file id
*       int*    comp_ws                 computer word size
*       int*    io_ws                   storage word size
*       float*  version                 EXODUSII interface version number
*
* revision history - 
*
*
*****************************************************************************/

#include <stdio.h>
#include <mpi.h>

#include "netcdf_par.h"
#include "exodusII.h"
#include "exodusII_int.h"

/*!  

The function ex_open() opens an existing exodus file and returns
an ID that can subsequently be used to refer to the file, the word
size of the floating point values stored in the file, and the version
of the exodus database (returned as a ``float'', regardless of the
compute or I/O word size). Multiple files may be ``open'' simultaneously.

\return In case of an error, ex_open() returns a negative
number. Possible causes of errors include:
  -  The specified file does not exist.
  -  The mode specified is something other than the predefined constant \fparam{EX_READ} or \fparam{EX_WRITE}.
  -  Database version is earlier than 2.0.

\param path The file name of the exodus file. This can be given as either an
            absolute path name (from the root of the file system) or a relative
            path name (from the current directory).

\param mode Access mode. Use one of the following predefined constants:
        -  \fparam{EX_READ} To open the file just for reading.
        -  \fparam{EX_WRITE} To open the file for writing and reading.

\param[in,out] comp_ws The word size in bytes (0, 4 or 8) of the floating point variables
               used in the application program. If 0 (zero) is passed, the default
               size of floating point values for the machine will be used and
               returned in this variable. WARNING: all exodus functions requiring
               reals must be passed reals declared with this passed in or returned
               compute word size (4 or 8).


\param[in,out] io_ws The word size in bytes (0, 4 or 8) of the floating 
                    point data as they are stored in the exodus file. If the word 
                    size does not match the word size of data stored in the file, 
                    a fatal error is returned. If this argument is 0, the word size 
                    of the floating point data already stored in the file is returned.

\param[out] version  Returned exodus database version number.

The following opens an exodus file named \file{test.exo} for read
only, using default settings for compute and I/O word sizes:

\code
#include "exodusII.h"
int CPU_word_size,IO_word_size, exoid;
float version;

CPU_word_size = sizeof(float);   \co{float or double}
IO_word_size = 0;                \co{use what is stored in file}

\comment{open exodus files}
exoid = ex_open ("test.exo",     \co{filename path}
                 EX_READ,        \co{access mode = READ}
		 &CPU_word_size, \co{CPU word size}
		 &IO_word_size,  \co{IO word size}
	         &version);      \co{ExodusII library version}
\endcode
 */

static int warning_output = 0;

int ex_open_par_int (const char  *path,
		     int    mode,
		     int   *comp_ws,
		     int   *io_ws,
		     float *version,
		     MPI_Comm comm,
		     MPI_Info info,
		     int    run_version)
{
  int exoid;
  int status, stat_att, stat_dim;
  nc_type att_type = NC_NAT;
  size_t att_len = 0;
  int old_fill;
  int file_wordsize;
  int dim_str_name;
  int int64_status = 0;
  int pariomode = NC_MPIPOSIX;
  
  char errmsg[MAX_ERR_LENGTH];

  exerrval = 0; /* clear error code */
 
  /* set error handling mode to no messages, non-fatal errors */
  ex_opts(exoptval);    /* call required to set ncopts first time through */

  if (run_version != EX_API_VERS_NODOT && warning_output == 0) {
    int run_version_major = run_version / 100;
    int run_version_minor = run_version % 100;
    int lib_version_major = EX_API_VERS_NODOT / 100;
    int lib_version_minor = EX_API_VERS_NODOT % 100;
    fprintf(stderr, "EXODUS: Warning: This code was compiled with exodus version %d.%02d,\n          but was linked with exodus library version %d.%02d\n          This is probably an error in the build process of this code.\n",
	    run_version_major, run_version_minor, lib_version_major, lib_version_minor);
    warning_output = 1;
  }
  

  if ((mode & EX_READ) && (mode & EX_WRITE)) {
    exerrval = EX_BADFILEMODE;
    sprintf(errmsg,"Error: Cannot specify both EX_READ and EX_WRITE");
    ex_err("ex_open",errmsg,exerrval); 
    return (EX_FATAL);
  }

  /* Check parallel io mode.  Valid is NC_MPIPOSIX or NC_MPIIO or NC_PNETCDF
   * Exodus uses different flag values; map to netcdf values
   */
  if (mode & EX_MPIPOSIX)
    pariomode = NC_MPIPOSIX;
  else if (mode & EX_MPIIO)
    pariomode = NC_MPIIO;
  else if (mode & EX_PNETCDF)
    pariomode = NC_PNETCDF;
  
  
  /* The EX_READ mode is the default if EX_WRITE is not specified... */
  if (!(mode & EX_WRITE)) { /* READ ONLY */
      if ((status = nc_open_par (path, NC_NOWRITE|NC_SHARE|pariomode, comm, info, &exoid)) != NC_NOERR)
	{
	  /* NOTE: netCDF returns an id of -1 on an error - but no error code! */
	  if (status == 0) {
	    exerrval = EX_FATAL;
	  }
	  else {
	    /* It is possible that the user is trying to open a netcdf4
	       file, but the netcdf4 capabilities aren't available in the
	       netcdf linked to this library. Note that we can't just use a
	       compile-time define since we could be using a shareable
	       netcdf library, so the netcdf4 capabilities aren't known
	       until runtime...
	  
	       Netcdf-4.X does not (yet?) have a function that can be
	       queried to determine whether the library being used was
	       compiled with --enable-netcdf4, so that isn't very
	       helpful.. 

	       At this time, query the beginning of the file and see if it
	       is an HDF-5 file and if it is assume that the open failure
	       is due to the netcdf library not enabling netcdf4 features...
	    */
	    int type = 0;
	    ex_check_file_type(path, &type);
	  
	    if (type == 5) {
	      /* This is an hdf5 (netcdf4) file. Since the nc_open failed,
		 the assumption is that the netcdf doesn't have netcdf4
		 capabilities enabled.  Tell the user...
	      */
	      fprintf(stderr,
		      "EXODUS: Error: Attempting to open the netcdf-4 file:\n\t'%s'\n\twith a netcdf library that does not support netcdf-4\n",
		      path);
	    }
	    exerrval = status;
	  }
	  sprintf(errmsg,"Error: failed to open %s read only",path);
	  ex_err("ex_open",errmsg,exerrval); 
	  return(EX_FATAL);
	} 
  }
  else /* (mode & EX_WRITE) READ/WRITE */
    {
	if ((status = nc_open_par (path, NC_WRITE|NC_SHARE|pariomode, comm, info, &exoid)) != NC_NOERR)
	  {
	    /* NOTE: netCDF returns an id of -1 on an error - but no error code! */
	    if (status == 0)
	      exerrval = EX_FATAL;
	    else
	      exerrval = status;
	    sprintf(errmsg,"Error: failed to open %s write only",path);
	    ex_err("ex_open",errmsg,exerrval); 
	    return(EX_FATAL);
	  } 

      /* turn off automatic filling of netCDF variables */
      if ((status = nc_set_fill (exoid, NC_NOFILL, &old_fill)) != NC_NOERR) {
	exerrval = status;
	sprintf(errmsg,
		"Error: failed to set nofill mode in file id %d",
		exoid);
	ex_err("ex_open", errmsg, exerrval);
	return (EX_FATAL);
      }

      stat_att = nc_inq_att(exoid, NC_GLOBAL, ATT_MAX_NAME_LENGTH, &att_type, &att_len);
      stat_dim = nc_inq_dimid(exoid, DIM_STR_NAME, &dim_str_name);
      if(stat_att != NC_NOERR || stat_dim != NC_NOERR) {
	nc_redef(exoid);
	if (stat_att != NC_NOERR) {
	  int max_so_far = 32;
	  nc_put_att_int(exoid, NC_GLOBAL, ATT_MAX_NAME_LENGTH, NC_INT, 1, &max_so_far);
	}

	/* If the DIM_STR_NAME variable does not exist on the database, we need to add it now. */
	if(stat_dim != NC_NOERR) {
	  /* Not found; set to default value of 32+1. */
	  int max_name = ex_default_max_name_length < 32 ? 32 : ex_default_max_name_length;
	  nc_def_dim(exoid, DIM_STR_NAME, max_name+1, &dim_str_name);
	}
	nc_enddef (exoid);
      }
    }

  /* determine version of EXODUS II file, and the word size of
   * floating point and integer values stored in the file
   */

  if ((status = nc_get_att_float(exoid, NC_GLOBAL, ATT_VERSION, version)) != NC_NOERR) {
    exerrval  = status;
    sprintf(errmsg,"Error: failed to get database version for file id: %d",
	    exoid);
    ex_err("ex_open",errmsg,exerrval);
    return(EX_FATAL);
  }
   
  /* check ExodusII file version - old version 1.x files are not supported */
  if (*version < 2.0) {
    exerrval  = EX_FATAL;
    sprintf(errmsg,"Error: Unsupported file version %.2f in file id: %d",
	    *version, exoid);
    ex_err("ex_open",errmsg,exerrval);
    return(EX_FATAL);
  }
   
  if (nc_get_att_int (exoid, NC_GLOBAL, ATT_FLT_WORDSIZE, &file_wordsize) != NC_NOERR)
    {  /* try old (prior to db version 2.02) attribute name */
      if (nc_get_att_int (exoid,NC_GLOBAL,ATT_FLT_WORDSIZE_BLANK,&file_wordsize) != NC_NOERR)
	{
	  exerrval  = EX_FATAL;
	  sprintf(errmsg,"Error: failed to get file wordsize from file id: %d",
		  exoid);
	  ex_err("ex_open",errmsg,exerrval);
	  return(exerrval);
	}
    }

  /* See if int64 status attribute exists and if so, what data is stored as int64 
   * Older files don't have the attribute, so it is not an error if it is missing
   */
  if (nc_get_att_int (exoid, NC_GLOBAL, ATT_INT64_STATUS, &int64_status) != NC_NOERR) {
    int64_status = 0; /* Just in case it gets munged by a failed get_att_int call */
  }
  
  /* Merge in API int64 status flags as specified by caller of function... */
  int64_status |= (mode & EX_ALL_INT64_API);
  
  /* initialize floating point and integer size conversion. */
  if (ex_conv_ini( exoid, comp_ws, io_ws, file_wordsize, int64_status ) != EX_NOERR ) {
    exerrval = EX_FATAL;
    sprintf(errmsg,
	    "Error: failed to initialize conversion routines in file id %d",
            exoid);
    ex_err("ex_open", errmsg, exerrval);
    return (EX_FATAL);
  }

  return (exoid);
}
