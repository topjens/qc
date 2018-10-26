@* Introduction.

@c
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void die(const char *s)
{
    fprintf(stderr, "%s", s);

    exit(EXIT_FAILURE);
}

@<Functions@>@;

int main(int argc, char *argv[])
{
    printf("%15.7e\n", boys(atof(argv[1]), atof(argv[2])));
    return 0;
}

@ @<Functions@>=
double gammln(double xx)
{
    double x, y ,tmp, ser;
    static double cof[6] = {76.18009172947146, -86.50532032941677,
           24.01409824083091, -1.231739572450155,
           0.1208650973866179e-2, -0.5395239384953e-5};
    int j;

    y = x = xx;
    tmp = x + 5.5;
    tmp -= (x + 0.5) * log(tmp);
    ser = 1.000000000190015;
    for(j = 0; j <= 5; j++) ser += cof[j] / ++y;
    
    return -tmp + log(2.5066282746310005*ser/x);
}

@ The function |gser| is implemented as follows:

@d ITMAX 1000
@d EPS 3.0e-7
@d FPMIN 1.0e-30

@<Functions@>=
void gser(double *gamser, double a, double x, double *gln)
{
	int n;
    double sum, del, ap;

    *gln = gammln(a);
    if(x <= 0.0) {
         if(x < 0.0) die("x less than 0 in routine gser");
         *gamser = 0.0;
         return;
    } else {
    	ap = a;
        del = sum = 1.0/a;
        for(n = 1; n <= ITMAX; n++) {
        	ap++;
            del *= x/ap;
            sum += del;
            if(fabs(del) < fabs(sum)*EPS) {
            	*gamser = sum * exp(-x + a * log(x) - (*gln));
                return;
            }
        }
        die("a too large, ITMAX too small in routine gser");
        return;
    }
}



void gcf(double *gammcf, double a, double x, double *gln)
{
    int i;
	double an,b,c,d,del,h;

	*gln=gammln(a);
    b=x+1.0-a;
    c=1.0/FPMIN;
    d=1.0/b;
    h=d;
    for (i = 1; i <= ITMAX; i++) {
		an = -i * (i - a);
        b += 2.0;
        d = an * d + b;
        if (fabs(d) < FPMIN) d=FPMIN;
        c = b + an / c;
        if (fabs(c) < FPMIN) c=FPMIN;
        d = 1.0 / d;
        del = d * c;
        h *= del;
        if (fabs(del - 1.0) < EPS) break;
    }
    if (i > ITMAX) die("a too large, ITMAX too small in gcf");
    *gammcf = exp(-x + a * log(x) - (*gln)) * h;
}




    
double gammp(double a, double x)
{
	double gamser, gammcf, gln;

    if(x < 0.0 || a <= 0.0) die("Invalid arguments in routine gammp");
    if(x < (a + 1.0)) {
    	gser(&gamser, a, x, &gln);
        return gamser * exp(gln);
    } else {
    	gcf(&gammcf, a, x, &gln);
        return ((1.0 - gammcf)*exp(gln));
    }
}

double boys(double m, double x)
{
    return gammp(m + 0.5, x)/(2.0*pow(x, m + 0.5));
}

@* Auxiliary primitive integrals $[00^^7c00]^{(m)}$. These auxiliary primitive integrals can be calculated from the Boys function as
$$
$[00^^7c00]^{(m)}$=
$$

@d PI_5_2 1e530098423926407

@<Functions@>=
double primitives(double alpha, coord A,
                  double beta, coord B,
                  double gamma, coord C,
                  double delta, coord D,
                  double m)
{
    double p = alpha + beta;
    double q = gamma + delta;
    double rho = p*q/(p+q);
    
	coord P = (alpha*A + beta*B)/p;
    coord Q = (gamma*C + delta*D)/q;

    double R = abs(P-Q);
    R *= R;

	return ((2*PI_5_2)/(p*q*sqrt(p+q))*boys(m, rho*R);
}
                  

@* Index.
