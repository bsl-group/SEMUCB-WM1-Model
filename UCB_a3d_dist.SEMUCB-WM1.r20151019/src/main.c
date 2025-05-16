#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <unistd.h>
#include <stdarg.h>
#include <libgen.h>
#include <math.h>

#include "ucb_A3d.h"

#define OUTPUT_FNAME "model-samples.out"

#define MIN_RADIUS 3480.0
#define MAX_RADIUS 6341.0
#define MAX_RADIUS_WARN 6321.0

const char version[] = "Version 0.4 (2015/10/19)";

const char default_output_file[] = OUTPUT_FNAME;

void usage(const char *name) {
	printf("\nUsage: %s -r radius [-r radius -r radius ...] (-d discretization | -l lon,lat | -i lon_lat_file) [-o output_file]\n\n", name);
	printf("\tAll argument units are km and degrees\n");
	printf("\tDefault output file: %s\n", default_output_file);
	printf("\n");
}

void error_exit(const char *name, const char *reason_fmt, ...)
{
	if (reason_fmt != NULL) {
		va_list args;
		va_start(args, reason_fmt);
		vprintf(reason_fmt, args);
		va_end(args);
	}
	usage(name);
	exit(EXIT_FAILURE);
}

int main(int argc, char **argv)
{
	int rflag, dflag, lflag, iflag;
	int i, nrad, nradalloc, npts;
	double dsc, lon, lat, rho,
	       vpv, vsv, qk, qmu, vph, vsh, eta,
	       dvs, dxi, vs, xi, vsv_1, vsh_1, vs_1, xi_1;
	double *r, r_moho, r_m;
	char ch, *name;
	const char *outfile, *infile;
	FILE *f_out, *f_in;

	printf("\nUCB A3d model evaluation tool %s\n\n", version);

	name = basename(argv[0]);

	outfile = default_output_file;

	r = NULL;
	nrad = nradalloc = 0;

	rflag = 0;
	dflag = 0;
	lflag = 0;
	iflag = 0;
	while ((ch = getopt(argc, argv, "r:d:l:o:i:")) != -1) {
		switch (ch) {
			case 'r':
				rflag = 1;
				if (nrad == nradalloc) {
					nradalloc += 10;
					r = realloc(r, nradalloc * sizeof(double));
				}
				if (sscanf(optarg, "%lf", &r[nrad++]) != 1)
					error_exit(name, "Error: nonsense radius argument \"%s\"\n", optarg);
				break;
			case 'd':
				dflag = 1;
				if (sscanf(optarg, "%lf", &dsc) != 1)
					error_exit(name, "Error: nonsense discretization argument \"%s\"\n", optarg);
				break;
			case 'l':
				lflag = 1;
				if (sscanf(optarg, "%lf,%lf", &lon, &lat) != 2)
					error_exit(name, "Error: nonsense location argument \"%s\"\n", optarg);
				break;
			case 'i':
				iflag = 1;
				infile = optarg;
				break;
			case 'o':
				outfile = optarg;
				break;
			case '?':
			default:
				error_exit(name, NULL);
		}
	}

	/* check radius set */
	if (!rflag)
		error_exit(name, "Error: no radius supplied\n");

	/* check radius values */
	for (i = 0; i < nrad; i++) {
		if (r[i] < MIN_RADIUS)
			error_exit(name, "Error: model is undefined below %.1f km\n", MIN_RADIUS);
		if (r[i] > MAX_RADIUS)
			error_exit(name, "Error: model is undefined above %.1f km\n", MAX_RADIUS);
		if (r[i] > MAX_RADIUS_WARN)
			printf("Warning: structure above ~ %.1fkm is less robust\n\n", MAX_RADIUS_WARN);
	}

	/* check lateral sampling set */
	if (dflag || lflag || iflag) {
		if (dflag + lflag + iflag > 1)
			error_exit(name, "Error: more than one of discretization, location, or input file flags set, but only one needed\n");
	} else
		error_exit(name, "Error: neither discretization, location, nor input file flags set\n");

	/* initialize model data */
	init_ucb_ref_();
	init_ucb_a3d_();
	init_moho_();

	if ((f_out = fopen(outfile, "w")) == NULL)
		error_exit(name, "Error: cannot open output file \"%s\" - %s\n", outfile, strerror(errno));

	/* loop over radii */
	for (i = 0; i < nrad; i++) {

		/* km -> m for get_ucb_ref_() */
		r_m = 1000 * r[i];

		/* get reference structure */
		get_ucb_ref_(&r_m, &rho, &vpv, &vsv, &qk, &qmu, &vph, &vsh, &eta);

		/* m -> km */
		vsv /= 1000.0;
		vsh /= 1000.0;

		/* (vsv, vsh) -> voigt, xi */
		vs = sqrt((2.0 * vsv * vsv + vsh * vsh) / 3.0);
		xi = vsh * vsh / (vsv * vsv);

		/* sample the model ... */

		if (lflag) {
			/* evaluate model at a single location */
			get_moho_radius_(&lat, &lon, &r_moho);
			if (r_moho <= r[i]) {
				/* in crust */
				fprintf(f_out, "%5.1f %8.3f %8.3f    nan    nan    nan    nan\n", r[i], lon, lat);
			} else {
				get_a3d_perturbation_(&lat, &lon, &r[i], "S", &dvs);
				get_a3d_perturbation_(&lat, &lon, &r[i], "X", &dxi);
				vs_1 = (1.0 + dvs) * vs;
				xi_1 = (1.0 + dxi) * xi;
				vsv_1 = sqrt(3.0 / (2.0 + xi_1)) * vs_1;
				vsh_1 = sqrt(xi_1) * vsv_1;
				fprintf(f_out, "%5.1f %8.3f %8.3f %6.2f %6.2f %6.3f %6.3f\n", r[i], lon, lat, 100 * dvs, 100 * dxi, vsv_1, vsh_1);
			}
		} else if (iflag) {
			/* evaluate the model at lon, lat points in infile
			 * TODO: repeatedly mapping over this file if multiple
			 * radii are specified. change to caching them on the
			 * first traversal */
			f_in = fopen(infile, "r");
			if (f_in == NULL)
				error_exit(name, "Error: cannot open input file \"%s\" - %s\n", infile, strerror(errno));
			npts = 0;
			while (fscanf(f_in, "%lf %lf", &lon, &lat) == 2) {
				get_moho_radius_(&lat, &lon, &r_moho);
				if (r_moho <= r[i]) {
					/* in crust */
					fprintf(f_out, "%5.1f %8.3f %8.3f    nan    nan    nan    nan\n", r[i], lon, lat);
				} else {
					get_a3d_perturbation_(&lat, &lon, &r[i], "S", &dvs);
					get_a3d_perturbation_(&lat, &lon, &r[i], "X", &dxi);
					vs_1 = (1.0 + dvs) * vs;
					xi_1 = (1.0 + dxi) * xi;
					vsv_1 = sqrt(3.0 / (2.0 + xi_1)) * vs_1;
					vsh_1 = sqrt(xi_1) * vsv_1;
					fprintf(f_out, "%5.1f %8.3f %8.3f %6.2f %6.2f %6.3f %6.3f\n", r[i], lon, lat, 100 * dvs, 100 * dxi, vsv_1, vsh_1);
				}
				npts += 1;
			}
			if (i == 0)
				printf("Read %i locations from %s\n", npts, infile);
			fclose(f_in);
		} else {
			/* evaluate model on a regular lat / lon grid */
			lat = -90.0;
			while (lat <= 90.0) {
				lon = -180.0;
				while (lon <= +180.0) {
					get_moho_radius_(&lat, &lon, &r_moho);
					if (r_moho <= r[i]) {
						/* in crust */
						fprintf(f_out, "%5.1f %8.3f %8.3f    nan    nan    nan    nan\n", r[i], lon, lat);
					} else {
						get_a3d_perturbation_(&lat, &lon, &r[i], "S", &dvs);
						get_a3d_perturbation_(&lat, &lon, &r[i], "X", &dxi);
						vs_1 = (1.0 + dvs) * vs;
						xi_1 = (1.0 + dxi) * xi;
						vsv_1 = sqrt(3.0 / (2.0 + xi_1)) * vs_1;
						vsh_1 = sqrt(xi_1) * vsv_1;
						fprintf(f_out, "%5.1f %8.3f %8.3f %6.2f %6.2f %6.3f %6.3f\n", r[i], lon, lat, 100 * dvs, 100 * dxi, vsv_1, vsh_1);
					}
					lon += dsc;
				}
				lat += dsc;
			}
		}
	}

	printf("Wrote model samples to %s\n\n", outfile);

	fclose(f_out);

	return 0;
}
