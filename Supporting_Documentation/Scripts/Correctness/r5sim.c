#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

// Sum up the coefficients of a l-long vector v
double sumcoeff(double *v, int l){
    double sum = 0.0;
    unsigned i;
    for(i=0; i<l; i++){
        sum += v[i];
    }
    return sum;
}

double mean(double *v, int l){
    double mean=sumcoeff(v, l) / 1.0 / l;
    return mean;
}

double variance(double *v, int l){
    double sumsq = 0.0;
    unsigned i;
    for(i=0;i<l;i++){
        sumsq += v[i] * v[i];
    }
    double sum = sumcoeff(v,l);
    double var = (sumsq - sum*sum/l) / (l-1);
    return var;
}

double choose(int n, int k)
{
    double res=1;
    for(int i=0; i<k; i++) res *= (n-i);
    for(int i=1; i<=k; i++) res /= i;
    return res;
}

int randominteger(int x)
{
    /* Replace arc4random() to work on non-BSD systems */
    return (int)(x*( (double)rand()/(double)RAND_MAX ));
    // return (int)(x*(arc4random()/(2.0+0xffffffff)));
}

// Sample a real following a Gaussian distribution
// of 0 mean and variance var
void gaussian(double *x, double *y, double var)
{
    double u1 = (double)rand()/(double)RAND_MAX;
    double u2 = (double)rand()/(double)RAND_MAX;

    double sqrt2 = sqrt(2.0);
    double r = sqrt(-log( 1 - u1 ));
    double c = cos( 2*M_PI*u2 );
    double s = sin( 2*M_PI*u2 );
    *x = sqrt(var) * sqrt2 * r * c;
    *y = sqrt(var) * sqrt2 * r * s;
    return;
}

// Sample a n-long vector of real Gaussians
void gaussian_err(double *err, int n, double var){
    unsigned i;
    double x, y;
    for(i=0; i<n; i++){
        if(i%2 == 0){
            gaussian(&x, &y, var);
            err[i] = x;
        } else{
            err[i] = y;
        }
    }
    return;
}

// Debug function, to print means and variances
static void meanvariance(double *v, int l){ 
    printf("\n\nEstimation of mean: \t\t%f\n", mean(v, l));    
    printf("Estimation of variance: \t\t%f\n", variance(v, l));
    
    return;
}

// Polynomial multiplication modulo the polynomial
// x^(n+1)-1, i.e., a convolution
void modN(int *as, int *extra, int *a,
        int *s, int mask, int n)
{
    for(int i=0;i<n;i++){
        as[i]=0;
        for(int j=0;j<=i;j++){
            as[i] += (a[i-j]*s[j]);
        }
        for(int j=i+2;j<n;j++){
            as[i] += (a[n+1+i-j]*s[j]);
        }
        if (as[i]>=0) as[i] &= mask;
        else as[i] = (1+mask) - ((-as[i])&mask);
    }
    *extra=0;
    for(int j=1;j<n;j++)
        *extra += (a[n-j]*s[j])&mask;
    if (*extra >= 0) *extra &= mask;
    else *extra = (1+mask) - ((-*extra)&mask);
}

// Polynomial multiplication modulo the prime-order
// cyclotomic polynomial x^n+...+x+1 = (x^(n+1)-1)/(x-1)
void modPhi(int *as, int *a, int *s, int mask, int n)
{
    int extra;
    modN(as,&extra,a,s,mask,n);
    for(int i=0;i<n;i++) {as[i] -= extra;
        if (as[i]>=0) as[i] &= mask;
        else as[i] = (1+mask) - ((-as[i])&mask);
    }
}

void print(char *x, int *a, int n)
{
    printf("print: n=%d\n",n);
    printf("%s=",x);
    for(int i=0;i<n;i++) printf(" %x",a[i]);
    printf("\n");
}

void shuffle(int *x,int n)
{
    int r,tmp;
    for(int k=n-1;k>=1;k--){
        r = randominteger(k);
        tmp = x[r];
        x[r]=x[k];
        x[k]=tmp;
    }
}

void modN_double(double *re, double *extra, int *r, 
        double *e, int n)
{
    for(int i=0;i<n;i++){
        re[i]=0;
        for(int j=0;j<=i;j++){
            re[i] += (r[i-j]*e[j]);
        }
        for(int j=i+2;j<n;j++){
            re[i] += (r[n+1+i-j]*e[j]);
        }
    }
    *extra=0;
    for(int j=1;j<n;j++) *extra += (r[n-j]*e[j]);
}

void modPhi_double(double *re, int *r, double *e, int n)
{
    double extra;
    modN_double(re,&extra,r,e,n);
    for(int i=0;i<n;i++) re[i] -= extra; 
}

// Rounding function, where the moduli are powers of 2
int rnd(int x, int shift)
{
    if (shift == 0) return x;
    int b = (x>>(shift-1))&1;
    return (x>>shift)+b;
}

void rndarray(int *out, int *in, int outmask,
        int shift, int n)
{
    for(int i=0; i<n; i++){
        out[i] = (rnd(in[i],shift)&outmask);
    }
}

// Simulate the ring variant of Round5
int main(int argc, char *argv[])
{
    // Whether XEf forward error correction is used
    int xef_used=(int)atoi(argv[1]);

    // Ring dimension, Round5 parameter n
    int n=(int)atoi(argv[2]);

    // Hamming-weight of Round5 secret-key vectors
    int h=(int)atoi(argv[3]);

    // Exponent of the Round5 large modulus q
    int qexp=(int)atoi(argv[4]);

    // Exponent of the Round5 rounding modulus p
    int pexp=(int)atoi(argv[5]);

    // Exponent of the Round5 ciphertext compression 
    // modulus t
    int texp=(int)atoi(argv[6]);

    // q/p
    int qp = 1<<(qexp-pexp);

    // q/t
    int qt = 1<<(qexp-texp);

    // Variances for rounding noises,
    // used to sample *equivalent* 
    // additive Gaussian and uniform noises
    double varqp = (qp*qp - 1.0) / 12.0;
    double varqt = (qt*qt - 1.0) / 12.0;

    int nruns = (int)atoi(argv[7]);

    // File to write simulation results to
    FILE *f;
    if (!(f=fopen(argv[8], "a"))) {
        printf("\nError opening file!\n");
        exit(1);
    }

    // Noise distribution, for pretty printing
    int n_gaussian = (int)atoi(argv[9]) & 1;              // 0 if false, 1 if true
    int n_unifadditive = (int)atoi(argv[10]) & 1;
    int n_rounding = (int)atoi(argv[11]) & 1;

    // number of runs after which to flush output to disk
    int intermediate_runs = (int)atoi(argv[12]);

    int qmask = (1<<qexp)-1;
    int pmask = (1<<pexp)-1;
    int tmask = (1<<texp)-1;

    // Round5 public parameter, the polynomial a
    int *a=(int *)malloc(n*sizeof(int));

    // Polynomial product of a and the 
    // initiator's/Alice's/decryptor's secret-key s
    int *as=(int *)malloc(n*sizeof(int));

    // Polynomial product of a and the 
    // initiator's/Alice's/decryptor's secret-key s
    int *ar=(int *)malloc(n*sizeof(int));

    // Public-key polynomial of the initiator/Alice/decryptor
    int *b=(int *)malloc(n*sizeof(int));

    // ``Raw key'' of the responder/Bob/encryptor
    int *br=(int *)malloc(n*sizeof(int));

    // First ciphertext component of the responder/Bob/encryptor
    int *u=(int *)malloc(n*sizeof(int));

    // Second ciphertext component of the responder/Bob/encryptor
    int *v=(int *)malloc(n*sizeof(int));

    // ``Raw key'' of the initiator/Alice/decryptor
    int *us=(int *)malloc(n*sizeof(int));
    int *x=(int *)malloc(n*sizeof(int));

    // When Round5 variants *without* ronding noise are
    // simulated (i.e., with additive Gaussian or uniform noise),
    // these are used to sample the independent noise that is 
    // then used in Round5
    double *eb = (double *)malloc(n*sizeof(double));
    double *eu = (double *)malloc(n*sizeof(double));
    double *ev = (double *)malloc(n*sizeof(double));
    double *ex = (double *)malloc(n*sizeof(double));
    double *reb = (double *)malloc(n*sizeof(double));
    double *seu = (double *)malloc(n*sizeof(double)); 

    // Secret-key polynomial of responder/Bob/decryptor
    int *r = (int *)malloc(n*sizeof(int));

    // Secret-key polynomial of initiator/Alice/decryptor
    int *s = (int *)malloc(n*sizeof(int));
    double *z=(double *)malloc(n*sizeof(double));

    // Decryption errors that occur in a Round5 exchange
    int *errors=(int *)calloc(n+1,sizeof(int));
    double *summederrors=(double *)calloc(n+1,sizeof(double));

    // Init the secret keys
    for(int i=0;i<n;i++) { s[i]=0; r[i]=0; }
    for(int i=0;i<h/2;i++){
        r[i]=1; s[i]=1;
        r[i+h/2]=-1; s[i+h/2]=-1;
    }

    // Output
    fprintf(f, "\nRound5 simulations for analyzing error behavior ");

    if(xef_used==0)         fprintf(f, "in the Prime cyclotomic polynomial ring,\n");
    else                    fprintf(f, "in the cyclic x^(n+1) - 1 polynomial ring,\n");

    if(n_gaussian==1)       fprintf(f, "with Gaussian errors sampled with variance ((q/p)^2 - 1)/12.");
    else if(n_unifadditive) fprintf(f, "with uniform additive noise in [0,q/p) and [0,q/t).");
    else if(n_rounding)     fprintf(f, "with actual rounding noise from Z_q->Z_p and Z_q->Z_t.");

    fprintf(f, "\n\nRunning %d simulations...",nruns);

    // First, parameters:
    fprintf(f, "\n\nParameters:\nn:\t\t%d, \nh:\t\t%d, \nlog2(q):\t%d, \nlog2(p):\t%d, \nlog2(t):\t%d,",n,h,qexp,pexp,texp);

    fprintf(f, "\n\n========Simulation results (Failure rates)=======\n\n(Running the ");

    fflush(f);

    int run;
    for(run=1; run<=nruns; run++){

        // Number of errors, i.e., positions in which the
        // final decrypted message differs from the 
        // original, occurring in this run
        int nerr=0;
        
        // Sample secret-keys
        shuffle(s,n);
        shuffle(r,n);
        
        if(n_rounding){
            // Simulate a Round5 cpapke variant with rounding noise

            int extra;
	    // We need the following only in a full 
	    // simulation, i.e., in actual rounding.
	    // For the other cases it is enough to 
	    // sample the error polynomials/vectors
	    // directly to estimate the error.
            // FULL simulation in case of rounding errors

            // Sample the public polynomial a
            for(int i=0;i<n;i++) a[i]=randominteger(1<<qexp);

            // Compute public-key b, using polynomial 
            // multiplication modulo \Phi_{n+1}(x)
            modPhi(as,a,s,qmask,n);
            // Coefficient wise rounding from ZZ_q to ZZ_p
            rndarray(b,as,pmask,qexp-pexp,n);

            // Compute ciphertext u, using polynomial 
            // multiplication modulo \Phi_{n+1}(x)
            modPhi(ar,a,r,qmask,n);
            // Coefficient wise rounding from ZZ_q to ZZ_p
            rndarray(u,ar,pmask,qexp-pexp,n);

            // If using error correction, ciphertext v
            // computation involves polynomial multiplication
            // modulo N_{n+1}(x), otherwise modulo \Phi_{n+1}(x)
	    if (xef_used==0) modPhi(br,b,r,pmask,n); else modN(br,&extra,b,r,pmask,n);
            // Coefficient wise rounding from ZZ_p to ZZ_t
            rndarray(v,br,tmask,pexp-texp,n);

            // If using error correction, decryption involves 
            // polynomial multiplication modulo N_{n+1}(x), 
            // otherwise modulo \Phi_{n+1}(x)
            if (xef_used==0) modPhi(us,u,s,pmask,n); else modN(us,&extra,u,s,pmask,n);
            rndarray(x,us,tmask,pexp-texp,n);

            // Compute differences between encrypted and decrypted
            // message, actually differences between the ``shared
            // secrets'' computed by Alice and Bob
            for(int i=0;i<n;i++) {
                z[i] = (v[i]-x[i])&tmask;
                // Compute number of positions where 
                // the error crosses the threshold, i.e.,
                // number of errors nerr in this run
                if (z[i] >= (1.0+tmask)/4 && z[i] < 3*(1.0+tmask)/4){ nerr++; }
            }
        }else{
            double extra;
            if(n_gaussian){
		// Sample Gaussian errors.
		gaussian_err(eb, n, varqp);
		gaussian_err(eu, n, varqp);
		gaussian_err(ev, n, varqt);
		gaussian_err(ex, n, varqt);
            }else if(n_unifadditive){
                // Sample noises from uniform intervals 
		for(int i=0;i<n;i++) eb[i]= ( (1<<(qexp-1)) - randominteger(1<<qexp) ) / 1.0 / (1<<pexp);
		for(int i=0;i<n;i++) eu[i]= ( (1<<(qexp-1)) - randominteger(1<<qexp) ) / 1.0 / (1<<pexp);
		for(int i=0;i<n;i++) ev[i]= ( (1<<(qexp-1)) - randominteger(1<<qexp) ) / 1.0 / (1<<texp);
		for(int i=0;i<n;i++) ex[i]= ( (1<<(qexp-1)) - randominteger(1<<qexp) ) / 1.0 / (1<<texp);
            }
            // Add independent noises and compute public keys
	    if (xef_used==0) modPhi_double(reb,r,eb,n); else modN_double(reb,&extra,r,eb,n);
	    if (xef_used==0) modPhi_double(seu,s,eu,n); else modN_double(seu,&extra,s,eu,n);

            // Compute differences between encrypted and decrypted
            // message, actually differences between the ``shared
            // secrets'' computed by Alice and Bob
	    for(int i=0;i<n;i++) {
		z[i] = reb[i] - seu[i] + ev[i] - ex[i];
                // Compute number of positions where the
                // error crosses the threshold, i.e.,
                // number of errors nerr in this run
		if (fabs(z[i]) >= (1.0+qmask)/4){ nerr++; }
	    }
        }

        // Histogram of the decryption error
        errors[nerr]++;

        if(run%intermediate_runs == 0){        
	    for(int i=1; i<=n; i++){
		for(int j=i; j<=n; j++){
                    // An error event in which there are j errors 
                    // contributes to all summederrors[i] for i<=j, 
                    // choose(j,i) is there since that's the number 
                    // of ways the i chosen bits can be chosen from 
                    // the j bits with an error. 
		    if (errors[j]) summederrors[i] += choose(j,i)*errors[j];
		}
	    }

            // Bit failure probability
	    double bfp = log(summederrors[1]/1.0/n/nruns)/log(2);
            // Bit failure probability, on the condition
            // that already one *error* has occurred.
	    double bfp_cond1err = log(2*summederrors[2]/1.0/summederrors[1]/(n-1))/log(2);
            // Bit failure probability, on the condition
            // that already two *errors* have occurred.
	    double bfp_cond2err = log(3*summederrors[3]/1.0/summederrors[2]/(n-2))/log(2);

	    fprintf(f, "%d'th run.)\n\nBit failure probability:\t\t\t%lf \n\nConditional bit failure probability,\nassuming that ONE error has already occurred:\t%lf \n\nConditional bit failure probability, \nassuming that TWO errors have already occurred: %lf\n\n=================================================\n",run,bfp,bfp_cond1err,bfp_cond2err);
            fflush(f);
        }
    }
    // Free everything
    free(a);
    free(as); free(ar);
    free(b);
    free(br);
    free(u);
    free(v);
    free(us);
    free(x);
    free(r); free(s);
    free(eb); free(eu); free(ev); free(ex);
    free(reb); free(seu);
    free(z);
    free(errors);
    free(summederrors);
    if(f!=NULL) fclose(f);
    return 0;
}
