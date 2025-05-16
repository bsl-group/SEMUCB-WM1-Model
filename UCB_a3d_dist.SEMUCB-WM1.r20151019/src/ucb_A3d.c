#include <stdio.h>
#include <stdlib.h>
#include <stdint.h> /* uint8_t */
#include <errno.h>
#include <string.h>
#include <math.h>

#include "constants.h"

/******************************************************************************
 * 1D model *******************************************************************
 ******************************************************************************/

#define REF_MODEL_FNAME "data/model.ref"
#define REF_MODEL_NPTS_MAX 2000
#define REF_MODEL_LINE_LEN 256

struct ref_model
{
	int ifanis, ifdeck, npts, icb, cmb, noc;
	double tref;
	double r[REF_MODEL_NPTS_MAX],  rho[REF_MODEL_NPTS_MAX], vpv[REF_MODEL_NPTS_MAX], vsv[REF_MODEL_NPTS_MAX], 
	       qk[REF_MODEL_NPTS_MAX], qmu[REF_MODEL_NPTS_MAX], vph[REF_MODEL_NPTS_MAX], vsh[REF_MODEL_NPTS_MAX], eta[REF_MODEL_NPTS_MAX];
	char model_name[REF_MODEL_LINE_LEN];
};
typedef struct ref_model ref_t;

static int ref_initialized = 0;
static ref_t cached_ref;

static void load_ref_model(char *model_file_name, ref_t *ref)
{
	int ipt;
	FILE *f_in;
	char line[REF_MODEL_LINE_LEN];

	if ((f_in = fopen(model_file_name, "r" )) == NULL){
		printf("Error [%s]: cannot open model file %s - %s\n", __func__, model_file_name, strerror(errno));
		exit(EXIT_FAILURE);
	}

	fgets(line, REF_MODEL_LINE_LEN, f_in);
	sscanf(line, "%s", ref->model_name);
	fgets(line, REF_MODEL_LINE_LEN, f_in);
	sscanf(line, "%i %lf %i", &ref->ifanis, &ref->tref, &ref->ifdeck);
	fgets(line, REF_MODEL_LINE_LEN, f_in);
	sscanf(line, "%i %i %i %i", &ref->npts, &ref->icb, &ref->cmb, &ref->noc);

	if (ref->npts > REF_MODEL_NPTS_MAX) {
		printf("Error [%s]: REF_MODEL_NPTS_MAX is too small\n", __func__);
		exit(EXIT_FAILURE);
	}

	for (ipt = 0; ipt < ref->npts; ipt++)
		fscanf(f_in, "%lf %lf %lf %lf %lf %lf %lf %lf %lf", 
		       &ref->r[ipt],  &ref->rho[ipt], &ref->vpv[ipt], &ref->vsv[ipt],
		       &ref->qk[ipt], &ref->qmu[ipt], &ref->vph[ipt], &ref->vsh[ipt],
		       &ref->eta[ipt]);

	fclose(f_in);
}

static void evaluate_ref(double *r, double *rho, double *vpv, double *vsv,
                         double *qk, double *qmu, double *vph, double *vsh, 
                         double *eta, ref_t *ref)
{
	int k, found = 0;
	for (k = 0; k < ref->npts - 1; k++) {
		/* interpolation - discon underside always wins */
		if (ref->r[k] <= *r && ref->r[k+1] >= *r) {
			double beta = (*r - ref->r[k]) / (ref->r[k+1] - ref->r[k]);
			*rho = (1.0 - beta) * ref->rho[k] + beta * ref->rho[k+1];
			*vpv = (1.0 - beta) * ref->vpv[k] + beta * ref->vpv[k+1];
			*vsv = (1.0 - beta) * ref->vsv[k] + beta * ref->vsv[k+1];
			*qk  = (1.0 - beta) * ref->qk[k]  + beta * ref->qk[k+1];
			*qmu = (1.0 - beta) * ref->qmu[k] + beta * ref->qmu[k+1];
			*vph = (1.0 - beta) * ref->vph[k] + beta * ref->vph[k+1];
			*vsh = (1.0 - beta) * ref->vsh[k] + beta * ref->vsh[k+1];
			*eta = (1.0 - beta) * ref->eta[k] + beta * ref->eta[k+1];
			found = 1;
			break;
		}
	}
	if (found == 0) {
		printf("Error [%s]: could not find model radius %f m\n", __func__, *r);
		exit(EXIT_FAILURE);
	}
}

void init_ucb_ref_()
{
	load_ref_model(REF_MODEL_FNAME, &cached_ref);
	ref_initialized = 1;
}

void get_ucb_ref_(double *r, double *rho, double *vpv, double *vsv,
                  double *qk, double *qmu, double *vph, double *vsh, 
                  double *eta)
{
	if (ref_initialized == 0) {
		printf("Error [%s]: called without initialization of model structure\n", __func__);
		exit(EXIT_FAILURE);
	}
	evaluate_ref(r, rho, vpv, vsv, qk, qmu, vph, vsh, eta, &cached_ref);
}

/******************************************************************************
 * 3D model *******************************************************************
 ******************************************************************************/

#define A3D_MODEL_FNAME "data/model.A3d"
#define GRID_FNAME_BASE "data/grid."
#define GRID_FNAME_LEN 256
#define PARAM 2
#define NBKNOT 20
#define NSKNOT 10242

/* for b-spline calculation */
#define TOL 1.5

/* constants */
#define RADIANS 0.017453292519943295

struct a3d_model {
	char desc[PARAM];
	int n, nd, nbspl, nsspl[PARAM], sspl_knot_levels[PARAM][NSKNOT];
	double bspl_knots[NBKNOT];
	double sspl_knot_theta[PARAM][NSKNOT];
	double c_sspl_knot_theta[PARAM][NSKNOT];
	double s_sspl_knot_theta[PARAM][NSKNOT];
	double sspl_knot_phi[PARAM][NSKNOT];
	double c_sspl_knot_phi[PARAM][NSKNOT];
	double s_sspl_knot_phi[PARAM][NSKNOT];
	double values[PARAM][NBKNOT][NSKNOT];
};
typedef struct a3d_model a3d_t;

static int mdl_initialized = 0;
static a3d_t cached_mdl;

/* fast parameter indexing */
#define IDX_CHAR_MIN 65 /* ascii A */
#define IDX_CHAR_MAX 90 /* ascii Z */
#define IDX_MAP_SIZE (IDX_CHAR_MAX - IDX_CHAR_MIN + 1)
#define IDX_NO_MATCH (-1)
static int idx_map[IDX_MAP_SIZE];

static void init_idx_map(a3d_t *mdl)
{
	int i;
	for (i = 0; i < IDX_MAP_SIZE; i++)
		idx_map[i] = IDX_NO_MATCH;
	for (i = 0; i < mdl->n; i++)
		idx_map[(uint8_t)(mdl->desc[i]) - IDX_CHAR_MIN] = i;
}

static int get_idx(char *desc)
{
	return idx_map[(uint8_t)(*desc) - IDX_CHAR_MIN];
}

static void load_a3d_model(const char *fname, a3d_t *mdl)
{
	FILE *f;
	int k, j, p;
	if ((f = fopen(fname, "r")) == NULL) {
		printf("Error [%s]: cannot open %s - %s\n", __func__, fname, strerror(errno));
		exit(EXIT_FAILURE);
	}
	if (fscanf(f, "%i", &mdl->n) != 1) {
		printf("Error [%s]: malformed A3d header (number of parameters)\n", __func__);
		exit(EXIT_FAILURE);
	}
	if (mdl->n > PARAM) {
		printf("Error [%s]: PARAM is too small\n", __func__);
		exit(EXIT_FAILURE);
	}
	for (k = 0; k < mdl->n; k++)
		if (fscanf(f, "%i %i %c", &mdl->nsspl[k], &mdl->nbspl, &mdl->desc[k]) != 3) {
			printf("Error [%s]: malformed A3d header (parameter spec)\n", __func__);
			exit(EXIT_FAILURE);
		} else if (mdl->nsspl[k] > NSKNOT) {
			printf("Error [%s]: NSKNOT is too small\n", __func__);
			exit(EXIT_FAILURE);
		} else if (mdl->nbspl > NBKNOT) {
			printf("Error [%s]: NBKNOT is too small\n", __func__);
			exit(EXIT_FAILURE);
		}
	if (fscanf(f, "%i", &mdl->nd) != 1) {
		printf("Error [%s]: malformed A3d header (number of discons)\n", __func__);
		exit(EXIT_FAILURE);
	}
	if (mdl->nd > 0) {
		printf("Error [%s]: discontinuities not supported\n", __func__);
		exit(EXIT_FAILURE);
	}
	for (k = 0; k < mdl->nbspl; k++)
		if (fscanf(f, "%lf", &mdl->bspl_knots[k]) != 1) {
			printf("Error [%s]: malformed A3d header (b-spline knots)\n", __func__);
			exit(EXIT_FAILURE);
		}
	for (p = 0; p < mdl->n; p++)
		for (j = 0; j < mdl->nbspl; j++)
			for (k = 0; k < mdl->nsspl[p]; k++)
				if (fscanf(f, "%lf", &mdl->values[p][j][k]) != 1) {
					printf("Error [%s]: malformed A3d values field (parameter %c)\n", __func__, mdl->desc[p]);
					exit(EXIT_FAILURE);
				}
	fclose(f);
	for (p = 0; p < mdl->n; p++) {
		int nsspl;
		char grid_fname[GRID_FNAME_LEN+1];
		sprintf(grid_fname, "%s%c", GRID_FNAME_BASE, mdl->desc[p]);
		if ((f = fopen(grid_fname, "r")) == NULL) {
			printf("Error [%s]: cannot open %s - %s\n", __func__, grid_fname, strerror(errno));
			exit(EXIT_FAILURE);
		}
		if (fscanf(f, "%i", &nsspl) != 1) {
			printf("Error [%s]: malformed grid file header\n", __func__);
			exit(EXIT_FAILURE);
		}
		if (mdl->nsspl[p] != nsspl) {
			printf("Error [%s]: number of knots in header of %s does not match %s\n", __func__, grid_fname, fname);
			exit(EXIT_FAILURE);
		}
		for (k = 0; k < mdl->nsspl[p]; k++)
			if (fscanf(f, "%lg %lg %i", &mdl->sspl_knot_phi[p][k], &mdl->sspl_knot_theta[p][k], &mdl->sspl_knot_levels[p][k]) != 3) {
				printf("Error [%s]: malformed grid knot entry in %s\n", __func__, grid_fname);
				exit(EXIT_FAILURE);
			} else {
				mdl->sspl_knot_theta[p][k] = RADIANS * (90 - mdl->sspl_knot_theta[p][k]);
				mdl->sspl_knot_phi[p][k] = RADIANS * mdl->sspl_knot_phi[p][k];
				/* precompute trigonometric functions */
				mdl->c_sspl_knot_theta[p][k] = cos(mdl->sspl_knot_theta[p][k]);
				mdl->s_sspl_knot_theta[p][k] = sin(mdl->sspl_knot_theta[p][k]);
				mdl->c_sspl_knot_phi[p][k] = cos(mdl->sspl_knot_phi[p][k]);
				mdl->s_sspl_knot_phi[p][k] = sin(mdl->sspl_knot_phi[p][k]);
			}
		fclose(f);
	}
}

static void fill_hh(double *hh, double *knot, int Nx)
{
	int ii;
	for (ii = 0; ii < Nx; ii++)
		hh[ii] = knot[ii + 1] - knot[ii];
}

static double fspl(int ord, int nknots, double *knot, double xi)
/* ord: number of rho(x)
   nknots : # of knkots = Nx+1 (Nx=index of highest spline)
   xi: point of interest
   splines defined as :
   f_i(x) = a_i(x-x_i)^3 + b_i(x-x_i)^2 + c_i(x-x_i) + d_i
*/
{
	int Nx;
	double rho_x;
	double coefa, coefb, coefc, coefd;
	static int init_hh = 0;
	static double hh[NBKNOT];

	Nx = nknots - 1;
	/* Compute vector hh of spacings */
	if (init_hh == 0) {
		fill_hh(hh, knot, Nx);
		init_hh = 1;
	}

	/* Consistency checks */
	if ((xi - (double)TOL) > knot[Nx])
		return 0.0;
	else if ((xi + (double)TOL) < knot[0])
		return 0.0;
	else if (ord > Nx)
		return 0.0;

	if (ord == 0) {		/* LHS */
		double denom;
		denom = 3. * hh[ord] * hh[ord] + 3. * hh[ord] * hh[ord + 1] + hh[ord + 1] * hh[ord + 1];
		if (xi >= knot[ord] && xi <= knot[ord + 1]) {	/* x0<=x<=x1 */
			coefa = 4. / (hh[ord] * (hh[ord] + hh[ord + 1]) * denom);
			coefb = 0.0;
			coefc = -12 / denom;
			coefd = 4 * (2 * hh[ord] + hh[ord + 1]) / denom;

			rho_x = coefa * (xi - knot[ord]) * (xi - knot[ord]) * (xi - knot[ord]);
			rho_x += coefb * (xi - knot[ord]) * (xi - knot[ord]);
			rho_x += coefc * (xi - knot[ord]);
			rho_x += coefd;
		} else if (xi > knot[ord + 1] && xi <= knot[ord + 2]) {	/* x1<=x<=x2 */
			coefa = -4. / (hh[ord + 1] * (hh[ord] + hh[ord + 1]) * denom);
			coefb = 12 / ((hh[ord] + hh[ord + 1]) * denom);
			coefc = -12. * hh[ord + 1] / ((hh[ord] + hh[ord + 1]) * denom);
			coefd = 4. * hh[ord + 1] * hh[ord + 1] / ((hh[ord] + hh[ord + 1]) * denom);

			rho_x = coefa * (xi - knot[ord + 1]) * (xi - knot[ord + 1]) * (xi - knot[ord + 1]);
			rho_x += coefb * (xi - knot[ord + 1]) * (xi - knot[ord + 1]);
			rho_x += coefc * (xi - knot[ord + 1]);
			rho_x += coefd;
		} else		/* x>x2 */
			rho_x = 0.0;
	}

	else if (ord == 1) {	/* LHS+1 */
		double denom, denomsum, dd;
		denom = (3. * hh[ord - 1] * hh[ord - 1] + 4. * hh[ord - 1] * hh[ord] + hh[ord] * hh[ord] + 2. * hh[ord - 1] * hh[ord + 1] + hh[ord] * hh[ord + 1]);
		denomsum = hh[ord - 1] + hh[ord] + hh[ord + 1];
		dd = denomsum * denom;
		if (xi >= knot[ord - 1] && xi <= knot[ord]) {
/* x0<=x<=x1 */
			coefa = -4. * (3. * hh[ord - 1] + 2. * hh[ord] + hh[ord + 1]) / (hh[ord - 1] * (hh[ord - 1] + hh[ord]) * dd);
			coefb = 0.;
			coefc = 12. / denom;
			coefd = 0.;

			rho_x = coefa * (xi - knot[ord - 1]) * (xi - knot[ord - 1]) * (xi - knot[ord - 1]);
			rho_x += coefb * (xi - knot[ord - 1]) * (xi - knot[ord - 1]);
			rho_x += coefc * (xi - knot[ord - 1]);
			rho_x += coefd;
		}

		else if (xi >= knot[ord] && xi <= knot[ord + 1]) {	/* x1<=x<=x2 */
			coefa = 4. * (2. * hh[ord - 1] * hh[ord - 1] + 6. * hh[ord - 1] * hh[ord] + 3. * hh[ord] * hh[ord] + 3. * hh[ord - 1] * hh[ord + 1] + 3. * hh[ord] * hh[ord + 1] + hh[ord + 1] * hh[ord + 1]) / (hh[ord] * (hh[ord - 1] + hh[ord]) * (hh[ord] + hh[ord + 1]) * dd);
			coefb = -12. * (3. * hh[ord - 1] + 2. * hh[ord] + hh[ord + 1]) / ((hh[ord - 1] + hh[ord]) * dd);
			coefc = 12. * (-2. * hh[ord - 1] * hh[ord - 1] + hh[ord] * hh[ord] + hh[ord] * hh[ord + 1]) / ((hh[ord - 1] + hh[ord]) * dd);
			coefd = 4. * hh[ord - 1] * (4. * hh[ord - 1] * hh[ord] + 3. * hh[ord] * hh[ord] + 2. * hh[ord - 1] * hh[ord + 1] + 3. * hh[ord] * hh[ord + 1]) / ((hh[ord - 1] + hh[ord]) * dd);

			rho_x = coefa * (xi - knot[ord]) * (xi - knot[ord]) * (xi - knot[ord]);
			rho_x += coefb * (xi - knot[ord]) * (xi - knot[ord]);
			rho_x += coefc * (xi - knot[ord]);
			rho_x += coefd;
		}

		else if (xi >= knot[ord + 1] && xi <= knot[ord + 2]) {	/* x2<=x<=x3 */
			dd *= (hh[ord] + hh[ord + 1]);
			coefa = -4. * (2. * hh[ord - 1] + hh[ord]) / (hh[ord + 1] * dd);
			coefb = 12. * (2. * hh[ord - 1] + hh[ord]) / dd;
			coefc = -12. * (2. * hh[ord - 1] + hh[ord]) * hh[ord + 1] / dd;
			coefd = 4. * (2. * hh[ord - 1] + hh[ord]) * hh[ord + 1] * hh[ord + 1] / dd;

			rho_x = coefa * (xi - knot[ord + 1]) * (xi - knot[ord + 1]) * (xi - knot[ord + 1]);
			rho_x += coefb * (xi - knot[ord + 1]) * (xi - knot[ord + 1]);
			rho_x += coefc * (xi - knot[ord + 1]);
			rho_x += coefd;
		} else		/* x>x3 */
			rho_x = 0.0;
	}

	else if (ord == Nx - 1) {	/* RHS-1 */
		double denom, denomsum, dd;
		denom = hh[ord - 2] * hh[ord - 1] + hh[ord - 1] * hh[ord - 1] + 2. * hh[ord - 2] * hh[ord] + 4. * hh[ord - 1] * hh[ord] + 3. * hh[ord] * hh[ord];
		denomsum = hh[ord - 2] + hh[ord - 1] + hh[ord];
		dd = denomsum * denom;
		if (xi >= knot[ord - 2] && xi <= knot[ord - 1]) {	/* x0<=x<=x1 */
			coefa = 4. * (hh[ord - 1] + 2. * hh[ord]) / (hh[ord - 2] * (hh[ord - 2] + hh[ord - 1]) * dd);
			coefb = coefc = coefd = 0.0;

			rho_x = coefa * (xi - knot[ord - 2]) * (xi - knot[ord - 2]) * (xi - knot[ord - 2]);
			rho_x += coefb * (xi - knot[ord - 2]) * (xi - knot[ord - 2]);
			rho_x += coefc * (xi - knot[ord - 2]);
			rho_x += coefd;
		}

		else if (xi >= knot[ord - 1] && xi <= knot[ord]) {	/* x1<=x<=x2 */
			coefa = -4. * (hh[ord - 2] * hh[ord - 2] + 3. * hh[ord - 2] * hh[ord - 1] + 3. * hh[ord - 1] * hh[ord - 1] + 3. * hh[ord - 2] * hh[ord] + 6. * hh[ord - 1] * hh[ord] + 2. * hh[ord] * hh[ord]) / (hh[ord - 1] * (hh[ord - 2] + hh[ord - 1]) * (hh[ord - 1] + hh[ord]) * dd);
			coefb = 12. * (hh[ord - 1] + 2. * hh[ord]) / ((hh[ord - 2] + hh[ord - 1]) * dd);
			coefc = 12. * hh[ord - 2] * (hh[ord - 1] + 2. * hh[ord]) / ((hh[ord - 2] + hh[ord - 1]) * dd);
			coefd = 4. * hh[ord - 2] * hh[ord - 2] * (hh[ord - 1] + 2. * hh[ord]) / ((hh[ord - 2] + hh[ord - 1]) * dd);

			rho_x = coefa * (xi - knot[ord - 1]) * (xi - knot[ord - 1]) * (xi - knot[ord - 1]);
			rho_x += coefb * (xi - knot[ord - 1]) * (xi - knot[ord - 1]);
			rho_x += coefc * (xi - knot[ord - 1]);
			rho_x += coefd;
		}

		else if (xi >= knot[ord] && xi <= knot[ord + 1]) {	/* x2<=x<=x3 */
			dd *= (hh[ord - 1] + hh[ord]);
			coefa = 4. * (hh[ord - 2] + 2. * hh[ord - 1] + 3. * hh[ord]) / (hh[ord] * dd);
			coefb = -12. * (hh[ord - 2] + 2. * hh[ord - 1] + 3. * hh[ord]) / dd;
			coefc = 12. * (-hh[ord - 2] * hh[ord - 1] - hh[ord - 1] * hh[ord - 1] + 2. * hh[ord] * hh[ord]) / dd;
			coefd = 4. * hh[ord] * (3. * hh[ord - 2] * hh[ord - 1] + 3. * hh[ord - 1] * hh[ord - 1] + 2. * hh[ord - 2] * hh[ord] + 4. * hh[ord - 1] * hh[ord]) / dd;

			rho_x = coefa * (xi - knot[ord]) * (xi - knot[ord]) * (xi - knot[ord]);
			rho_x += coefb * (xi - knot[ord]) * (xi - knot[ord]);
			rho_x += coefc * (xi - knot[ord]);
			rho_x += coefd;
		} else		/* x>x4 */
			rho_x = 0.0;
	}

	else if (ord == Nx) {	/* RHS */
		double denom;
		denom = (hh[ord - 2] + hh[ord - 1]) * (hh[ord - 2] * hh[ord - 2] + 3. * hh[ord - 2] * hh[ord - 1] + 3. * hh[ord - 1] * hh[ord - 1]);
		if (xi >= knot[ord - 2] && xi <= knot[ord - 1]) {	/* x0<=x<=x1 */
			coefa = 4. / (hh[ord - 2] * denom);
			coefb = coefc = coefd = 0.0;

			rho_x = coefa * (xi - knot[ord - 2]) * (xi - knot[ord - 2]) * (xi - knot[ord - 2]);
			rho_x += coefb * (xi - knot[ord - 2]) * (xi - knot[ord - 2]);
			rho_x += coefc * (xi - knot[ord - 2]);
			rho_x += coefd;
		}

		else if (xi >= knot[ord - 1] && xi <= knot[ord]) {	/* x1<=x<=x2 */
			coefa = -4. / (hh[ord - 1] * denom);
			coefb = 12 / denom;
			coefc = 12 * hh[ord - 2] / denom;
			coefd = 4. * hh[ord - 2] * hh[ord - 2] / denom;

			rho_x = coefa * (xi - knot[ord - 1]) * (xi - knot[ord - 1]) * (xi - knot[ord - 1]);
			rho_x += coefb * (xi - knot[ord - 1]) * (xi - knot[ord - 1]);
			rho_x += coefc * (xi - knot[ord - 1]);
			rho_x += coefd;
		}

		else		/* x>x2 */
			rho_x = 0.0;
	}

	else {			/* Away from borders */
		double denom1, denom2, denom;
		denom1 = hh[ord - 2] + hh[ord - 1] + hh[ord] + hh[ord + 1];
		if (xi >= knot[ord - 2] && xi <= knot[ord - 1]) {	/* x0<=x<=x1 */
			coefa = 4. / (hh[ord - 2] * (hh[ord - 2] + hh[ord - 1]) * (hh[ord - 2] + hh[ord - 1] + hh[ord]) * denom1);
			coefb = coefc = coefd = 0.;

			rho_x = coefa * (xi - knot[ord - 2]) * (xi - knot[ord - 2]) * (xi - knot[ord - 2]);
			rho_x += coefb * (xi - knot[ord - 2]) * (xi - knot[ord - 2]);
			rho_x += coefc * (xi - knot[ord - 2]);
			rho_x += coefd;
		} else if (xi >= knot[ord - 1] && xi <= knot[ord]) {	/* x1<=x<=x2 */
			denom2 = (hh[ord - 2] + hh[ord - 1]) * (hh[ord - 2] + hh[ord - 1] + hh[ord]);
			denom = denom1 * denom2;

			coefa = -4. * (hh[ord - 2] * hh[ord - 2] + 3. * hh[ord - 2] * hh[ord - 1] + 3. * hh[ord - 1] * hh[ord - 1] + 2. * hh[ord - 2] * hh[ord] + 4. * hh[ord - 1] * hh[ord] + hh[ord] * hh[ord] + hh[ord - 2] * hh[ord + 1] + 2. * hh[ord - 1] * hh[ord + 1] + hh[ord] * hh[ord + 1]) / (hh[ord - 1] * (hh[ord - 1] + hh[ord]) * (hh[ord - 1] + hh[ord] + hh[ord + 1]) * denom);
			coefb = 12. / denom;
			coefc = 12. * hh[ord - 2] / denom;
			coefd = 4. * hh[ord - 2] * hh[ord - 2] / denom;

			rho_x = coefa * (xi - knot[ord - 1]) * (xi - knot[ord - 1]) * (xi - knot[ord - 1]);
			rho_x += coefb * (xi - knot[ord - 1]) * (xi - knot[ord - 1]);
			rho_x += coefc * (xi - knot[ord - 1]);
			rho_x += coefd;
		}

		else if (xi >= knot[ord] && xi <= knot[ord + 1]) {	/* x2<=x<=x3 */
			denom2 = (hh[ord - 1] + hh[ord]) * (hh[ord - 2] + hh[ord - 1] + hh[ord]) * (hh[ord - 1] + hh[ord] + hh[ord + 1]);
			denom = denom1 * denom2;

			coefa = 4. * (hh[ord - 2] * hh[ord - 1] + hh[ord - 1] * hh[ord - 1] + 2. * hh[ord - 2] * hh[ord] + 4. * hh[ord - 1] * hh[ord] + 3. * hh[ord] * hh[ord] + hh[ord - 2] * hh[ord + 1] + 2. * hh[ord - 1] * hh[ord + 1] + 3. * hh[ord] * hh[ord + 1] + hh[ord + 1] * hh[ord + 1]) / (hh[ord] * (hh[ord] + hh[ord + 1]) * denom);
			coefb = -12. * (hh[ord - 2] + 2. * hh[ord - 1] + 2. * hh[ord] + hh[ord + 1]) / denom;
			coefc = 12. * (-hh[ord - 2] * hh[ord - 1] - hh[ord - 1] * hh[ord - 1] + hh[ord] * hh[ord] + hh[ord] * hh[ord + 1]) / denom;
			coefd = 4. * (2. * hh[ord - 2] * hh[ord - 1] * hh[ord] + 2. * hh[ord - 1] * hh[ord - 1] * hh[ord] + hh[ord - 2] * hh[ord] * hh[ord] + 2. * hh[ord - 1] * hh[ord] * hh[ord] + hh[ord - 2] * hh[ord - 1] * hh[ord + 1] + hh[ord - 1] * hh[ord - 1] * hh[ord + 1] + hh[ord - 2] * hh[ord] * hh[ord + 1] + 2. * hh[ord - 1] * hh[ord] * hh[ord + 1]) / denom;

			rho_x = coefa * (xi - knot[ord]) * (xi - knot[ord]) * (xi - knot[ord]);
			rho_x += coefb * (xi - knot[ord]) * (xi - knot[ord]);
			rho_x += coefc * (xi - knot[ord]);
			rho_x += coefd;
		}

		else if (xi >= knot[ord + 1] && xi <= knot[ord + 2]) {	/* x3<=x<=x4 */
			denom2 = (hh[ord] + hh[ord + 1]) * (hh[ord - 1] + hh[ord] + hh[ord + 1]);
			denom = denom1 * denom2;

			coefa = -4. / (hh[ord + 1] * denom);
			coefb = 12 / denom;
			coefc = -12 * hh[ord + 1] / denom;
			coefd = 4. * hh[ord + 1] * hh[ord + 1] / denom;

			rho_x = coefa * (xi - knot[ord + 1]) * (xi - knot[ord + 1]) * (xi - knot[ord + 1]);
			rho_x += coefb * (xi - knot[ord + 1]) * (xi - knot[ord + 1]);
			rho_x += coefc * (xi - knot[ord + 1]);
			rho_x += coefd;
		}

		else		/* x>x4 */
			rho_x = 0.0;
	}
	return rho_x;
}

static double evaluate_a3d(double *theta, double *phi, double *r, char *desc, a3d_t *mdl)
{
	/* standard inter-knot distances (radians) by level */
	const double adel[] = {
		1.106538745764405,
		0.5532693728822025,
		0.27576202181510406,
		0.13788101090755203,
		0.06894050545377602,
		0.034557519189487726,
		0.017278759594743863};
	/* cosine of twice the standard inter-knot distances (radians) by level */
	const double cos2adel[] = {
		-0.5990235985155857,
		 0.4477590878387697,
		 0.8517269341430476,
		 0.9622179935292854,
		 0.9905094632383088,
		 0.9976125063612252,
		 0.9994029483549729};
	/* local vars */
	double val, ct, st, cp, sp, v[NBKNOT], h[NSKNOT];
	int j, k, l, nix, ix[NSKNOT], p_ix;
	/* fetch parameter index 
	 * default behavior: return zero on match failure */
	p_ix = get_idx(desc);
	if (p_ix == IDX_NO_MATCH) {
#ifdef ERROR_ON_MATCH_FAILURE
		printf("Error [%s]: cannot find requested parameter '%c' (IDX_NO_MATCH)\n",
		       __func__, *desc);
		exit(EXIT_FAILURE);
#else
		return 0.0;
#endif
	}
	/* evaluate radial b-splines */
	for (k = 0; k < mdl->nbspl; k++)
		v[k] = fspl(k, mdl->nbspl, mdl->bspl_knots, *r);
	/* precompute trigonometric functions for this sampling point */
	ct = cos(*theta);
	st = sin(*theta);
	cp = cos(*phi);
	sp = sin(*phi);
	/* loop over all spherical-spline knots
	 * general idea: use precomputed values and identities to avoid _any_
	 * evaluation of trig functions unless absolutely necessary */
	for (nix = 0, k = 0; k < mdl->nsspl[p_ix]; k++) {
		double arg, delta;
		arg = mdl->c_sspl_knot_theta[p_ix][k] * ct +
		      mdl->s_sspl_knot_theta[p_ix][k] * st *
		      (cp * mdl->c_sspl_knot_phi[p_ix][k] +
		       sp * mdl->s_sspl_knot_phi[p_ix][k]); /* cos(a - b) = ca * cb + sa * sb */
		if (arg >  1.0) arg =  1.0;
		if (arg < -1.0) arg = -1.0;
		if (arg > cos2adel[mdl->sspl_knot_levels[p_ix][k] - 1]) {
			delta = acos(arg) / adel[mdl->sspl_knot_levels[p_ix][k] - 1];
			if (delta < 1.0)
				h[nix] = ( 0.75 * delta * delta * delta - 1.5 * delta * delta + 1.0 );
			else
				h[nix] = 0.25 * ( 2.0 - delta ) * ( 2.0 - delta ) * ( 2.0 - delta );
			ix[nix++] = k;
		}
	}
	/* finally, summation */
	val = 0.0;
	for (j = 0; j < mdl->nbspl; j++)
		if (v[j] > 0.0)
			for (l = 0; l < nix; l++)
				val += mdl->values[p_ix][j][ix[l]] * v[j] * h[l];
	return val;
}

void init_ucb_a3d_()
{
	load_a3d_model(A3D_MODEL_FNAME, &cached_mdl);
	init_idx_map(&cached_mdl);
	mdl_initialized = 1;
}

void get_a3d_perturbation_(double *lat, double *lon, double *r, char *desc, double *dp)
{
	double theta = RADIANS * (90 - *lat), phi = RADIANS * *lon;
	if (mdl_initialized == 0) {
		printf("Error [%s]: called without initialization of model structure\n", __func__);
		exit(EXIT_FAILURE);
	}
	*dp = evaluate_a3d(&theta, &phi, r, desc, &cached_mdl);
}

/******************************************************************************
 * Moho surface ***************************************************************
 ******************************************************************************/

#define MOHO_FNAME "data/moho.dat"
#define MOHO_NLAT 181
#define MOHO_NLON 361

struct moho_surface
{
	int nlat, nlon;
	double dlon, dlat;
	double depth[MOHO_NLAT][MOHO_NLON];
};
typedef struct moho_surface moho_t;

static int moho_initialized = 0;
static moho_t cached_moho;

static void load_moho_surface(const char *fname, moho_t *moho)
{
	FILE *f;
	int i, j;
	if ((f = fopen(fname, "r")) == NULL) {
		printf("Error [%s]: cannot open %s - %s\n", __func__, fname, strerror(errno));
		exit(EXIT_FAILURE);
	}
	if (fscanf(f, "%i %i", &moho->nlat, &moho->nlon) != 2) {
		printf("Error [%s]: malformed moho header (number of lats / lons)\n", __func__);
		exit(EXIT_FAILURE);
	}
	if (moho->nlat != MOHO_NLAT) {
		printf("Error [%s]: unexpected number of lats specified in %s\n", __func__, fname);
		exit(EXIT_FAILURE);
	}
	if (moho->nlon != MOHO_NLON) {
		printf("Error [%s]: unexpected number of lons specified in %s\n", __func__, fname);
		exit(EXIT_FAILURE);
	}
	for (i = 0; i < moho->nlat; i++)
		for (j = 0; j < moho->nlon; j++)
			if (fscanf(f, "%lf", &moho->depth[i][j]) != 1) {
				printf("Error [%s]: could not read moho depth # %i\n", __func__, i * moho->nlon + j);
				exit(EXIT_FAILURE);
			}
	moho->dlat = 180.0 / (moho->nlat - 1);
	moho->dlon = 360.0 / (moho->nlon - 1);
	fclose(f);
}

static double evaluate_moho_depth(double *lat, double *lon)
{
	int ix, iy;
	double plon, clat, ay, ax;
	/* bound the longitude to the interval [0,360) */
	plon = *lon;
	if (plon < 0)
		plon += 360.0;
	if (plon >= 360.0)
		plon -= 360.0;
	/* define colatitude */
	clat = 90.0 - *lat;
	/* derive indices for bilinear interp */
	ix = plon / cached_moho.dlon;
	if (ix < 0)
		ix = 0;
	if (ix >= cached_moho.nlon - 1)
		ix = cached_moho.nlon - 2;
	iy = (180.0 - clat) / cached_moho.dlat;
	if (iy < 0)
		iy = 0;
	if (iy >= cached_moho.nlat - 1)
		iy = cached_moho.nlat - 2;
	/* bilinear interp */
	ay = ((180.0 - clat) - iy * cached_moho.dlat) / cached_moho.dlat;
	ax = (plon - ix * cached_moho.dlon) / cached_moho.dlon;
	return 
		(1.0 - ay) * (1.0 - ax) * cached_moho.depth[iy][ix] +
		       ay  * (1.0 - ax) * cached_moho.depth[iy+1][ix] +
		(1.0 - ay) *        ax  * cached_moho.depth[iy][ix+1] +
		       ay  *        ax  * cached_moho.depth[iy+1][ix+1];
}

void init_moho_()
{
	load_moho_surface(MOHO_FNAME, &cached_moho);
	moho_initialized = 1;
}

void get_moho_radius_(double *lat, double *lon, double *radius)
{
	if (moho_initialized == 0) {
		printf("Error [%s]: called without initialization of moho structure\n", __func__);
		exit(EXIT_FAILURE);
	}
	*radius = R_EARTH_KM - evaluate_moho_depth(lat, lon);
}
