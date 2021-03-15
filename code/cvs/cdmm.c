#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

SEXP cdmm(SEXP y, SEXP x, SEXP fac, SEXP lams, SEXP relax, SEXP maxvar, SEXP maxit, SEXP tol) {
	R_len_t	n = length(y), p = ncols(x), nlam = length(lams), i, j, k;
	int		maxit_cd = INTEGER(maxit)[0], maxit_mm = INTEGER(maxit)[1], itcd, itmm, nvar;
	double	*r_y = REAL(y), *r_x = REAL(x), *r_fac = REAL(fac), *r_lams = REAL(lams),
			mu = REAL(relax)[0], maxv = REAL(maxvar)[0], tol_cd = REAL(tol)[0], tol_mm = REAL(tol)[1],
			*sol, *bet, *res, alp, t0, lam, lam1, mu1, sum, del, old, dif, norm, norm_old;
	SEXP	ans;
	
	sol = (double *) R_alloc(p*nlam, sizeof(double));
	bet = (double *) R_alloc(p, sizeof(double));
	res = (double *) R_alloc(n, sizeof(double));
	
	for (i = 0; i < p; i++) bet[i] = 0.0;
	for (i = 0; i < n; i++) res[i] = r_y[i];
	mu1 = 1.0/(1.0 + mu); sum = 0.0; norm = 0.0;
	for (j = 0; j < nlam; j++) {
		lam = r_lams[j]; lam1 = lam*mu1; alp = 0.0;
		for (itmm = 0; itmm < maxit_mm; itmm++) {
			for (itcd = 0; itcd < maxit_cd; itcd++) {
				del = 0.0; norm_old = norm; norm = 0.0;
				for (i = 0; i < p; i++) {
					old = bet[i]; t0 = 0.0;
					for (k = 0; k < n; k++) t0 += res[k]*r_x[k + i*n];
					t0 = (t0 - mu*(alp + sum))*mu1 + bet[i];
					bet[i] = copysign(fdim(fabs(t0), lam1), t0);
					if (bet[i] != old) {
						dif = bet[i] - old; sum += r_fac[i]*dif;
						for (k = 0; k < n; k++) res[k] -= r_x[k + i*n]*dif;
						del += fabs(dif);
					}
					norm += fabs(bet[i]);
				}
				if (del < tol_cd*norm_old) break;
			}
			if (fabs(sum) < tol_mm) break;
			alp += sum;
		}
		nvar = 0;
		for (i = 0; i < p; i++) {
			sol[i + j*p] = bet[i];
			if (bet[i] != 0.0) nvar++;
		}
		if (R_FINITE(maxv) && nvar > maxv) break;
	}
	
	PROTECT(ans = allocVector(VECSXP, 2));
	SET_VECTOR_ELT(ans, 0, allocMatrix(REALSXP, p, j));
	for (i = 0; i < p*j; i++) REAL(VECTOR_ELT(ans, 0))[i] = sol[i];
	SET_VECTOR_ELT(ans, 1, allocVector(REALSXP, j));
	for (i = 0; i < j; i++) REAL(VECTOR_ELT(ans, 1))[i] = r_lams[i];
	
	UNPROTECT(1);
	return ans;
}
