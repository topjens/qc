\def\frac#1#2{{\begingroup #1\endgroup\over#2}}

@* Introduction. We have to go over the following steps:

\item{(1)} $\displaystyle R^2=(\omega+\mu)^2+\nu$
\item{(2)} $\displaystyle\frac{\theta^2}{2\Delta}=\frac{1}{2\Delta}\left/(\sigma_P+\sigma_Q)\right.$
\item{(3)} $\displaystyle\frac{T}{2\Delta}=\frac{\theta^2}{2\Delta}R^2$
\item{(4)} \quad\quad\quad IF $\displaystyle\frac{T}{2\Delta} < \frac{T_{crit}}{2\Delta}$
\item{(5)} $\displaystyle j = INT\left(\frac{T}{2\Delta}\right)$
\item{(6)} $\displaystyle G_L(T) = g_0(j,L) + \frac{T}{2\Delta}\left(g_1(j,L)+\frac{T}{2\Delta}\left(g_2(j,L)+\frac{T}{2\Delta}g_3(j,L)\right)\right)$
\item{(7)} $\displaystyle \exp(-T) = e_0(j,L) + \frac{T}{2\Delta}\left(e_1(j,L)+\frac{T}{2\Delta}\left(e_2(j,L)+\frac{T}{2\Delta}e_3(j,L)\right)\right)$
\item{(8)} $\displaystyle G_m(T)$
\item{(9a)} $\displaystyle [0]^{(m)}$
\item{} \quad\quad\quad ELSE
\item{(9b)} $\displaystyle [0]^{(m)}$

@d VERBOSITY_LEVEL 5

@c
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <math.h>

void vvprintf(int level, const char *s, ...)
{
	if(level <= VERBOSITY_LEVEL) {
		va_list args;
		va_start(args, s);
		vprintf(s, args);
		va_end(args);
	}
}

int main(int argc, char *argv[])
{
	@<Prepare Chebyshev@>@;

	return 0;
}

@ Chebyshev error theory learns us that $$\Delta = \left(2^n(n+1)!\epsilon\right)^{1/(n+1)}.$$ Since we elect to use cubic interpollation ($n=3$) this simplefies to $$\Delta = \left(192\epsilon\right)^{(1/4)}.$$

@d eps 1.0e-12

@<Prepare Chebyshev@>=
double Delta = sqrt(sqrt(192*eps));
double epslog = log(eps * 2.0);

@ We use Newton-Raphson to approximate $T_{crit}$, with an initial guess $T_{crit}=-\log2\epsilon$.

@<Prepare Chebyshev@>=
double Tcrit = -epslog;
double change; /* to store the change between steps */
register int i; /* to count the number of steps */

vvprintf(5, "Delta = %.15f\n", Delta);
vvprintf(5, "Calculating Tcrit\n");
for(i = 0; i < 100; i++) {
	change = (Tcrit + log(Tcrit) + epslog)/(1.0 + 1.0/Tcrit);
	@\
	vvprintf(5, "%d -- %.15e\n", i, Tcrit);
	@\
	Tcrit -= change;
	if(fabs(change) < 1.0e-6)
		break;
}
@\
vvprintf(5, "%d -- %.15e [FINAL]\n", i + 1, Tcrit);

@ Next we will generate our interpolation tables $g$ and $e$.

@<Prepare Chebyshev@>=
@<Generate Chebyshev interpolation tables@>@;

@ @<Generate Cheb...@>=
int table_size;
double *X;
double *e;
double *g;

table_size = (int)(Tcrit/2.0/Delta+0.5);

vvprintf(5, "Table size = %d (from Delta to %dDelta)\n", table_size, 2*table_size-1);

X = malloc(sizeof *X * table_size);

@<Populate $X_i$@>@;

free(X);

@ Now we populate $X_i$, using $X_i = (2i+1)\Delta$.

@<Populate $X_i$@>=
for(i = 0; i < table_size; i++)
	X[i] = (2*i+1)*Delta;

@* Index.
