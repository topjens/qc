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
@d NaN 0x7fffffff

@c
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <math.h>

double FACT[20] = {1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0, 5040.0,
40320.0, 362880.0, 3628800.0, 39916800.0, 479001600.0, 6227020800.0,
87178291200.0, 1307674368000.0, 20922789888000.0, 355687428096000.0,
6402373705728000.0, 121645100408832000.0};

double fact(int n)
{
	if(n < 20)
		return FACT[n];
	else
		return NaN;
}

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

@f Delta TeX

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
double a[4];
double *e[4];
double *g[4];

table_size = (int)(Tcrit/2.0/Delta+0.5);

vvprintf(5, "Table size = %d (from Delta to %dDelta)\n", table_size, 2*table_size-1);

X = malloc(sizeof *X * table_size);

@<Populate $X_i$@>@;
@<Populate $e_k(i)$@>@;

free(e[0]);
free(e[1]);
free(e[2]);
free(e[3]);
free(X);

@ Now we populate $X_i$, using $X_i = (2i+1)\Delta$.

@<Populate $X_i$@>=
for(i = 0; i < table_size; i++)
	X[i] = (2*i+1)*Delta;

@ To make matters easier, we shall calculate $a$ separatly for each value of $k$. $$a_k = (2-\delta_{0k}).$$ If we wait with adding $\exp{-X_i}$ until the end, we only need to calculate four $a$ values.

@<Populate $e_k(i)$@>=
int m;

a[0] = 1.0; a[1] = 1.0; a[2] = 1.0/2.0; a[3] = 1.0/6.0;

double fact = 1.0;
double term;

double hafdel = -Delta / 2.0;
double hafde2 = hafdel * hafdel;
double r4Delt = 1.0/(4.0*Delta);

for(m = 1; m < 100; m++) {
	vvprintf(5, "*");
	fact *= hafde2 / (m * m);
	term = fact;
	
	if(a[0] == a[0] + term)
		break;
	a[0] += term;
	term /= (m + 1.0);
	a[1] += term;
	term /= (m + 2.0);
	a[2] += term;
	term /= (m + 3.0);
	a[3] += term;
}

vvprintf(5,"\na[0] = %20.15f\ta[1] = %20.15f\ta[2] = %20.15f\ta[3] = %20.15f\n", a[0], a[1], a[2], a[3]);

fact = -Delta;
a[1] *= fact;
fact *= hafdel;
a[2] *= fact;
fact *= hafdel;
a[3] *= fact;

vvprintf(5,"\na[0] = %20.15f\ta[1] = %20.15f\ta[2] = %20.15f\ta[3] = %20.15f\n", a[0], a[1], a[2], a[3]);

@ Constructing the $e_k(i)$ from the $a_k(i)$ can be done in one loop.

@<Populate $e_k(i)$@>=
e[0] = malloc(sizeof *e[0] * table_size);
e[1] = malloc(sizeof *e[1] * table_size);
e[2] = malloc(sizeof *e[2] * table_size);
e[3] = malloc(sizeof *e[3] * table_size);

double a00 = a[0] - a[2];
double a01 = a[1] - 3.0*a[3];
double a02 = 2.0*a[2];
double a03 = 4.0*a[3];
double a10 = (2.0*a[1]-6.0*a[3])*r4Delt;
double a11 = (8.0*a[2])*r4Delt;
double a12 = (24.0*a[3])*r4Delt;
double a20 = (8.0*a[2])*r4Delt*r4Delt;
double a21 = (48.0*a[3])*r4Delt*r4Delt;
double a30 = (32.0*a[3])*r4Delt*r4Delt*r4Delt;

vvprintf(5, "e[1][200] = %14.12f\n", e[1][200]);

for(i = 0; i < table_size; i++) {
	double X = -(i+i+1);
	double TNeg = X * Delta;
	double ExpT = exp(TNeg);
	e[0][i] = ExpT*(a00+X*(a01+X*(a02+X*a03)));
	e[1][i] = ExpT*(a10+X*(a11+X*a12));
	e[2][i] = ExpT*(a20+X*a21);
	e[3][i] = ExpT*a30;
	//printf("e[3][%d] = %20.10f\n", i, e[3][i]);
}

double test = 20;
double t2 = 2.0*test;
int j = t2*r4Delt;
printf("j = %d\n", j);
printf("exp(-%f)=%e\te(-%f)=%e\n", test, exp(-test), test, e[0][j]+t2*(e[1][j]+t2*(e[2][j]+t2*e[3][j])));

@* Index.
