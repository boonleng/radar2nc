//
// radar2wdssii.c
// converts UF radar data into WDSS-II friendly netcdf file
//
// Created by Boon Leng Cheong on 11/2/2011
// Copyright 2011 Boon Leng Cheong. All rights reserved
//
#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <xlocale.h>
#include <math.h>
#include <dirent.h>

#include <rsl.h>

#include <netcdf.h>

/* Math macros, these are in GNU C++ but not in GNU C */
#ifndef NAN
#define NAN          (0.0/0.0)
#endif
#ifndef MAX
#define MAX(X,Y)     ((X)>(Y)?(X):(Y))
#endif
#ifndef MIN
#define MIN(X,Y)     ((X)>(Y)?(Y):(X))
#endif
#ifndef SQ
#define SQ(X)        ((X)*(X))
#endif
#ifndef ABS
#define ABS(X)       ((X)<0?-(X):(X))
#endif


#define W2_MISSING_DATA       -99900.0
#define W2_RANGE_FOLDED       -99901.0
#define BL_MAX_STRING         256

enum {
	OPTION_NONE                   = 0,
	OPTION_WDSS_FOLDER_STRUCTURE  = 1,
	OPTION_GLOBAL_TIMESTAMP       = 1 << 1,
	OPTION_FIRST_RAY_TIMESTAMP    = 1 << 2
};

int verbose = 0, data_tz = 0;
//int option_flag = OPTION_WDSS_FOLDER_STRUCTURE;
int option_flag = OPTION_NONE;
char output_dir[BL_MAX_STRING] = "netcdf/";

// Some convenient functions to minimize copy & paste ...
static int isDirExist(char *path_name) {
	DIR   *dirptr;
	if ((dirptr = opendir(path_name)) != NULL) {
    	closedir(dirptr);
		return 1;
    }
	return 0;
}

static int put_global_text_att(const int ncid, const char *att, const char *text) {
	return nc_put_att_text(ncid, NC_GLOBAL, att, strlen(text), text);
}

/////////////////////////////////////////////////////////
//
//   w r i t e _ n e t c d f
//
/////////////////////////////////////////////////////////
static int write_netcdf(const Radar *radar) {
	
	int r, iv, is, ir, ig;
	int ncid;	

	char symbol;
	char mystr[BL_MAX_STRING], dest_dir[BL_MAX_STRING], time_str[BL_MAX_STRING];
	char prod_name[BL_MAX_STRING], prod_colormap[BL_MAX_STRING], prod_unit[BL_MAX_STRING];
	
	int az_dimid, gate_dimid, dimids[2];
	int az_varid, el_varid, bw_varid, gw_varid, dat_varid;
	float *array_naz, *array2D;
	float tmpf;

	char site_name[5]; site_name[4] = '\0';
	strncpy(site_name, radar->h.radar_name, 4);

	struct tm *t;
	long tnum_long;

	time(&tnum_long);
	t = gmtime(&tnum_long); // initialize struct tm to get the proper tm_isdst

	Volume *volume = NULL;
	Sweep *sweep = NULL;
	Ray *ray = NULL, *zray = NULL;

#define BL_MAX_NR     1440
#define BL_MAX_NG     4096
	
	if ((array_naz = (float *)malloc(BL_MAX_NR*sizeof(float))) == NULL) {
		fprintf(stderr, "ERROR: Unable to allocate memory.\n");
		return 1;
	}
	if ((array2D = (float *)malloc(BL_MAX_NR*BL_MAX_NG*sizeof(float))) == NULL) {
		fprintf(stderr, "ERROR: Unable to allocate memory for 2D array.\n");
		return 1;
	}

	// Go through all the volumes
	for (iv=0; iv<18; iv++) {
		if (radar->v[iv]) {
			switch (iv) {
				case DZ_INDEX:
					if (radar->v[CZ_INDEX] == NULL) {
						symbol = 'Z';
						snprintf(prod_name, BL_MAX_STRING, "Reflectivity");
						snprintf(prod_colormap, BL_MAX_STRING, "Reflectivity");
						snprintf(prod_unit, BL_MAX_STRING, "dBZ");
					} else {
						symbol = 'X';
					}
					break;
				case CZ_INDEX:
					symbol = 'Z';
					snprintf(prod_name, BL_MAX_STRING, "Reflectivity");
					snprintf(prod_colormap, BL_MAX_STRING, "Reflectivity");
					snprintf(prod_unit, BL_MAX_STRING, "dBZ");
					break;
				case VR_INDEX:
					symbol = 'V';
					snprintf(prod_name, BL_MAX_STRING, "Velocity");
					snprintf(prod_colormap, BL_MAX_STRING, "Velocity");
					snprintf(prod_unit, BL_MAX_STRING, "MetersPerSecond");
					break;
				case SW_INDEX:
					symbol = 'W';
					snprintf(prod_name, BL_MAX_STRING, "SpectrumWidth");
					snprintf(prod_colormap, BL_MAX_STRING, "SpectrumWidth");
					snprintf(prod_unit, BL_MAX_STRING, "MetersPerSecond");
					break;
				case DR_INDEX:
					symbol = 'D';
					snprintf(prod_name, BL_MAX_STRING, "Differential_Reflectivity");
					snprintf(prod_colormap, BL_MAX_STRING, "Differential_Reflectivity");
					snprintf(prod_unit, BL_MAX_STRING, "dB");
					break;
				case PH_INDEX:
					symbol = 'P';
					snprintf(prod_name, BL_MAX_STRING, "PhiDP");
					snprintf(prod_colormap, BL_MAX_STRING, "PhiDP");
					snprintf(prod_unit, BL_MAX_STRING, "Degrees");
					break;
				case RH_INDEX:
					symbol = 'R';
					snprintf(prod_name, BL_MAX_STRING, "RhoHV");
					snprintf(prod_colormap, BL_MAX_STRING, "RhoHV");
					snprintf(prod_unit, BL_MAX_STRING, "Unitless");
					break;
				default:
					symbol = 'X';
					break;
			}
			if (symbol == 'X') {
				continue;
			}
			// Current volume
			volume = radar->v[iv];
			//printf("number of tilts = %d\n", radar->v[iv] );

			//printf("number of tilts = %d\n", volume);

			for (is=0; is<volume->h.nsweeps; is++) {
			  //printf("is = %d\n", is);
				if ((sweep = volume->sweep[is]) == NULL) {
					continue;
				}
				float this_elevation = roundf(sweep->h.elev*10.0f)/10.0f;
				if (option_flag & OPTION_WDSS_FOLDER_STRUCTURE) {
					snprintf(dest_dir, BL_MAX_STRING, "%s%s/%05.2f/", output_dir, prod_name, this_elevation);
				} else {
					snprintf(dest_dir, BL_MAX_STRING, "%s", output_dir);
				}
				// Create the folder if needed
				if (!isDirExist(dest_dir)) {
					snprintf(mystr, BL_MAX_STRING, "mkdir -p %s", dest_dir);
					if (verbose)
						printf("%s\n", mystr);
					system(mystr);
				}
				
//printf("okay variable = %i  %d %d\n", is, option_flag, OPTION_GLOBAL_TIMESTAMP);

				if (option_flag & OPTION_GLOBAL_TIMESTAMP) {
					t->tm_year = radar->h.year - 1900;
					t->tm_mon  = radar->h.month - 1;
					t->tm_mday = radar->h.day;
					t->tm_hour = radar->h.hour;
					t->tm_min  = radar->h.minute;
					t->tm_sec  = radar->h.sec;

				} else {
					if (option_flag & OPTION_FIRST_RAY_TIMESTAMP) {

						if ((ray = sweep->ray[0]) == NULL)
							continue;
					} else {

						int serv=0;
						//printf("here b %p %ld  is:%d\n", t, (long)mktime(t), is);
						//printf(" --> %s\n", radar == NULL?"YES":"NO");
						//printf(" --> %s\n", (radar->v[serv]) == NULL?"YES":"NO");
						//printf("ray %p\n", (radar->v[ser_v]->sweep[is]));
						//printf(" --> %s\n", (radar->v[0]->sweep[is]) == NULL?"YES":"NO");
						//printf(" --> %s\n", (radar->v[0]->sweep[is]->ray[0]) == NULL?"YES":"NO");
						//printf("number of tilts = %d\n", volume->h.nsweeps);

						//if ((radar->v[0] == NULL) || (radar->v[0]->sweep[is] == NULL) || ((ray = radar->v[0]->sweep[is]->ray[0]) == NULL)) { 
					if (radar->v[serv] == NULL)
					{ 	
						serv =serv+1; 
//printf(" serv %i\n", serv);

					}
					//else{
						// printf("ivol=%d  is=%d\n", serv, is);
						if ((ray = radar->v[serv]->sweep[is]->ray[0]) == NULL) {

							continue;
						}
					}


					t->tm_year = ray->h.year - 1900;
					t->tm_mon  = ray->h.month - 1;
					t->tm_mday = ray->h.day;
					t->tm_hour = ray->h.hour;
					t->tm_min  = ray->h.minute;
					t->tm_sec  = ray->h.sec;

				}


				// local time, this machine
				tnum_long = (long)mktime(t);
				// adjust to be UTC (relative to this machine)
				tnum_long -= timezone;        
				strftime(time_str, 32, "%Y%m%d-%H%M%S", gmtime((time_t *)&tnum_long));
   
				// Pseudo timestamp filename based on the sweep's time
				tnum_long += data_tz;
				strftime(time_str, 32, "%Y%m%d-%H%M%S", gmtime((time_t *)&tnum_long));
				if (option_flag & OPTION_WDSS_FOLDER_STRUCTURE) {
					snprintf(mystr, BL_MAX_STRING, "%s%s.netcdf", dest_dir, time_str);
				} else {
					snprintf(mystr, BL_MAX_STRING, "%s%s-%s-E%.1f-%c.nc", dest_dir, site_name, time_str, this_elevation, symbol);					
				}
				tnum_long -= data_tz;


				// Start the file
				if ((r = nc_create(mystr, NC_CLOBBER, &ncid)) > 0) {
					fprintf(stderr, "  NetCDF Error: %s\n", nc_strerror(r));
					return 1;
				}

				// Define the main two dimensions
				const int nr = MIN(BL_MAX_NR, sweep->h.nrays);
				const int ng = MIN(BL_MAX_NG, ray->h.nbins);
				if (sweep->h.nrays >= BL_MAX_NR) {
					fprintf(stderr, "WARNING: %c, Number of radials = %d, truncated to %d.\n", symbol, sweep->h.nrays, nr);
				}
				if (sweep->h.nrays >= BL_MAX_NG) {
					fprintf(stderr, "WARNING: %c, Number of radials = %d, truncated to %d.\n", symbol, sweep->h.nrays, ng);
				}
				nc_def_dim(ncid, "Azimuth", nr, &az_dimid);
				nc_def_dim(ncid, "Gate", ng, &gate_dimid);
				dimids[0] = az_dimid;
				dimids[1] = gate_dimid;

				if (verbose)
					printf("Output filename: %s  %d x %d\n", mystr, nr, ng);
				
				// Define variables
				nc_def_var(ncid, "Azimuth", NC_FLOAT, 1, &az_dimid, &az_varid);
				nc_put_att_text(ncid, az_varid, "Units", 7, "Degrees");
				nc_def_var(ncid, "Elevation", NC_FLOAT, 1, &az_dimid, &el_varid);
				nc_put_att_text(ncid, el_varid, "Units", 7, "Degrees");
				nc_def_var(ncid, "BeamWidth", NC_FLOAT, 1, &az_dimid, &bw_varid);
				nc_put_att_text(ncid, bw_varid, "Units", 7, "Degrees");
				nc_def_var(ncid, "GateWidth", NC_FLOAT, 1, &az_dimid, &gw_varid);
				nc_put_att_text(ncid, gw_varid, "Units", 6, "Meters");
				nc_def_var(ncid, prod_name, NC_FLOAT, 2, dimids, &dat_varid);
				nc_put_att_text(ncid, dat_varid, "Units", strlen(prod_unit), prod_unit);

				// Global attributes - WDSS-II required
				nc_put_att_float(ncid, NC_GLOBAL, "Elevation", NC_FLOAT, 1, &this_elevation);
				put_global_text_att(ncid, "ElevationUnits", "Degrees");
				tmpf = (float)ray->h.range_bin1;
				nc_put_att_float(ncid, NC_GLOBAL, "RangeToFirstGate", NC_FLOAT, 1, &tmpf);
				put_global_text_att(ncid, "RangeToFirstGateUnits", "Meters");
				tmpf = W2_MISSING_DATA;
				nc_put_att_float(ncid, NC_GLOBAL, "MissingData", NC_FLOAT, 1, &tmpf);
				tmpf = W2_RANGE_FOLDED;
				nc_put_att_float(ncid, NC_GLOBAL, "RangeFolded", NC_FLOAT, 1, &tmpf);
				nc_put_att_text(ncid, NC_GLOBAL, "TypeName", strlen(prod_name), prod_name);
				nc_put_att_text(ncid, NC_GLOBAL, "DataType", 9, "RadialSet");
				tmpf = radar->h.latd + radar->h.latm/60.0f + radar->h.lats/3600.0f;
				nc_put_att_float(ncid, NC_GLOBAL, "Latitude", NC_FLOAT, 1, &tmpf);
				tmpf = radar->h.lond + radar->h.lonm/60.0f + radar->h.lons/3600.0f;
				nc_put_att_float(ncid, NC_GLOBAL, "Longitude", NC_FLOAT, 1, &tmpf);
				//tmpf = radar->h.height - ray->h.alt;
				tmpf = radar->h.height;
				nc_put_att_float(ncid, NC_GLOBAL, "Height", NC_FLOAT, 1, &tmpf);
				//int t = (int)tnum_long;
				//nc_put_att_int(ncid, NC_GLOBAL, "Time", NC_INT, 1, &t);
				nc_put_att_long(ncid, NC_GLOBAL, "Time", NC_LONG, 1, &tnum_long);
				tmpf = ray->h.sec - floorf(ray->h.sec);
				nc_put_att_float(ncid, NC_GLOBAL, "FractionalTime", NC_FLOAT, 1, &tmpf);
				put_global_text_att(ncid, "attributes", "Nyquist_Vel Unit radarName vcp ColorMap");
				put_global_text_att(ncid, "Nyquist_Vel-unit", "MetersPerSecond");
				nc_put_att_float(ncid, NC_GLOBAL, "Nyquist_Vel-value", NC_FLOAT, 1, &ray->h.nyq_vel);
				put_global_text_att(ncid, "Unit-unit", "dimensionless");
				put_global_text_att(ncid, "Unit-value", prod_unit);
				put_global_text_att(ncid, "radarName-unit", "dimensionless");
				put_global_text_att(ncid, "radarName-value", site_name);
				put_global_text_att(ncid, "vcp-unit", "dimensionless");
				sprintf(mystr, "%d", radar->h.vcp);
				put_global_text_att(ncid, "vcp-value", mystr);
				put_global_text_att(ncid, "ColorMap-unit", "dimensionless");
				put_global_text_att(ncid, "ColorMap-value", prod_colormap);

				// Global attributes - more for house keeping
				put_global_text_att(ncid, "RadarParameters", "PRF PulseWidth MaximumRange");
				put_global_text_att(ncid, "PRF-unit", "Hertz");
				nc_put_att_int(ncid, NC_GLOBAL, "PRF-value", NC_INT, 1, &ray->h.prf);
				put_global_text_att(ncid, "PulseWidth-unit", "MicroSeconds");
				nc_put_att_float(ncid, NC_GLOBAL, "PulseWidth-value", NC_FLOAT, 1, &ray->h.pulse_width);
				put_global_text_att(ncid, "MaximumRange-unit", "KiloMeters");
				nc_put_att_float(ncid, NC_GLOBAL, "MaximumRange-value", NC_FLOAT, 1, &ray->h.unam_rng);
				
				// Some additional info
				put_global_text_att(ncid, "CoversionSoftware", "radar2nc v1.1");
				put_global_text_att(ncid, "ContactInformation", "http://arrc.ou.edu");

				nc_enddef(ncid);

				// printf(" serv %i\n", nr);

				// Data
				for (ir=0; ir<nr; ir++)
					array_naz[ir] = sweep->ray[ir]->h.azimuth;
				nc_put_var_float(ncid, az_varid, array_naz);

				for (ir=0; ir<nr; ir++)
					array_naz[ir] = sweep->ray[ir]->h.elev;
				nc_put_var_float(ncid, el_varid, array_naz);

				for (ir=0; ir<nr; ir++)
					array_naz[ir] = sweep->ray[ir]->h.beam_width;
				nc_put_var_float(ncid, bw_varid, array_naz);

				for (ir=0; ir<nr; ir++)
					array_naz[ir] = sweep->ray[ir]->h.gate_size;
				nc_put_var_float(ncid, gw_varid, array_naz);

				float *b = array2D;
				for (ir=0; ir<nr; ir++) {
					if ((ray = sweep->ray[ir]) == NULL)
						continue;
					zray =  radar->v[DZ_INDEX]->sweep[is]->ray[ir];
					for (ig=0; ig<ng; ig++) {
						//if (ray->range[ig] == 0) {
						if (ray->range[ig] == 0 || zray->range[ig] == 0) {
							*b++ = W2_MISSING_DATA;
						} else if (ray->range[ig] == 1) {
							*b++ = W2_RANGE_FOLDED;
						} else {
							*b++ = ray->h.f(ray->range[ig]);
						}
					}
				}
				nc_put_var_float(ncid, dat_varid, array2D);

				nc_close(ncid);
			}
		}
	}

	free(array_naz);
	
	return 0;
}


/////////////////////////////////////////////////////////
//
//   s h o w _ h e l p
//
/////////////////////////////////////////////////////////
static void show_help(const char *name) {
	fprintf(stderr, "\n"
			"USAGE\n"
			"     %s [OPTIONS] [FILENAME(s)]\n\n"
			"OPTIONS\n"
			"     Addtional parameters may be specified to change the default behaviors\n\n"
			"     -v     Produce verbose output.\n"
			"            %s will list each file as it is written.\n"
			"            Additional -v options will provide additional detail.\n\n"
			"     -t \033[04moffset\033[0m\n"
			"            Introduce a timezone offset in hours\n\n"
			"     -d \033[04mdirectory\033[0m\n"
			"            Specified the prefix of the output directory.\n"
			"            If this option is not provided, the default output directory\n"
			"            netcdf will be used.\n\n"
			"     -s     Generate all output in a single directory.\n\n"
			"     -g     Use a global timestamp for each volume.\n\n"
			"     -r     Use the 1-st ray's timestamp for each sweep.\n"
			"\n",
			name, name);
	fprintf(stderr, 
			"EXAMPLES\n"
			"     The following creates netcdf files in a folder data\n"
			"            %s -d \033[04mdata\033[0m \033[04m*.uf\033[0m\n\n"
			"     The following creates netcdf files with timezone offset of +8 hours\n"
			"            %s -t \033[04m8\033[0m \033[04m*.uf\033[0m\n\n",
			name, name);
}

/////////////////////////////////////////////////////////
//
//   m a i n
//
/////////////////////////////////////////////////////////
int main(int argc, char **argv) {

	Radar *radar;
	char a;
	char site[5]; site[4] = '\0';
	// Default output directory
	//sprintf(output_dir, "netcdf/");
	
	while ((a = getopt(argc, argv, "v?ht:d:sgr")) != -1) {
		switch (a) {
			case 'v':
				verbose++;
				break;
			case 'h':
				show_help(argv[0]);
				return EXIT_SUCCESS;
				break;
			case 't':
				printf("Using timezone offset = %d --> %d\n", (int)timezone, (int)atoi(optarg));
				data_tz = (int)atoi(optarg) * 3600;
				break;
			case 'd':
				strncpy(output_dir, optarg, BL_MAX_STRING);
				if (strlen(output_dir) < BL_MAX_STRING-1 && output_dir[strlen(output_dir)-1] != '/') {
					strcat(output_dir, "/");
				}
				break;
			case 's':
				option_flag &= ~OPTION_WDSS_FOLDER_STRUCTURE;
				break;
			case 'g':
				option_flag |= OPTION_GLOBAL_TIMESTAMP;
				break;
			case 'r':
				option_flag |= OPTION_FIRST_RAY_TIMESTAMP;
				break;
			case '?':
			default:
				exit(EXIT_FAILURE);
				break;
		}
	}
	// Check output destination
	char cmd[BL_MAX_STRING + 10];
	if (!isDirExist(output_dir)) {
		sprintf(cmd, "mkdir -p %s", output_dir);
		system(cmd);
	}

	if (verbose) {
		printf("%s\n", cmd);
		printf("Output directory: %s\n", output_dir);
		if (verbose > 1)
			RSL_radar_verbose_on();
	}

	RSL_read_these_sweeps("all", NULL);

	// If no input, show_help
	if (optind == argc) {
		show_help(argv[0]);
		exit(EXIT_FAILURE);
	}
	// Treat the rests as filenames
	while (optind < argc) {
		printf("%s\n", argv[optind]);
		// Read the data
		//if (strstr(argv[optind], ".uf") == NULL)
		char *c = argv[optind];
		if (*c++ == 'K' && (*c >= 'A' && *c++ <= 'Z') && (*c >= 'A' && *c++ <= 'Z') && (*c >= 'A' && *c++ <= 'Z')) {
			strncpy(site, argv[optind], 4);
			printf("Assume a WSR-88D\n");
			radar = RSL_anyformat_to_radar(argv[optind], site);
		} else {
			printf("Assume not a WSR-88D - %s\n", argv[optind]);
			radar = RSL_anyformat_to_radar(argv[optind]);
		}
		if (radar == NULL) {
			continue;
		}
		if (verbose) {
			printf("radar @ 0x%p\n", radar);
			printf("Radar->h.date       : %2.2d/%2.2d/%2.2d\n", radar->h.month, radar->h.day, radar->h.year);
			printf("Radar->h.time       : %2.2d:%2.2d:%f\n", radar->h.hour, radar->h.minute, radar->h.sec);
			// Should be one of these: "wsr88d", "lassen", "uf","nsig", "mcgill","kwajalein", "rsl", "toga".
			printf("Radar->h.radar_type : %s\n", radar->h.radar_type);
			printf("Radar->h.nvolumes   : %d\n", radar->h.nvolumes);
			printf("Radar->h.radar_name : %s\n", radar->h.radar_name);
		}
		
		// Write the data... here's the work
		write_netcdf(radar);

		optind++;
	}
	return EXIT_SUCCESS;
}
