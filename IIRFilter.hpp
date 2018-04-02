#ifndef __IIRFILTER__
#define __IIRFILTER__
#define MAXPZ 10
#include <vector>
#include <string>
#define _USE_MATH_DEFINES
#include <math.h>
#include <complex>

#define PI	    M_PI
#define TWOPI	    (2.0 * PI)
#define MAXORDER    10
#define MAXPOLES    (2*MAXORDER)    /* to allow for doubling of poles in BP filter */


#define opt_be 0x00001	/* -Be		Bessel characteristic	       */
#define opt_bu 0x00002	/* -Bu		Butterworth characteristic     */
#define opt_ch 0x00004	/* -Ch		Chebyshev characteristic       */
#define opt_re 0x00008	/* -Re		Resonator		       */
#define opt_pi 0x00010	/* -Pi		proportional-integral	       */

#define opt_lp 0x00020	/* -Lp		lowpass			       */
#define opt_hp 0x00040	/* -Hp		highpass		       */
#define opt_bp 0x00080	/* -Bp		bandpass		       */
#define opt_bs 0x00100	/* -Bs		bandstop		       */
#define opt_ap 0x00200	/* -Ap		allpass			       */

#define opt_a  0x00400	/* -a		alpha value		       */
#define opt_l  0x00800	/* -l		just list filter parameters    */
#define opt_o  0x01000	/* -o		order of filter		       */
#define opt_p  0x02000	/* -p		specified poles only	       */
#define opt_w  0x04000	/* -w		don't pre-warp		       */
#define opt_z  0x08000	/* -z		use matched z-transform	       */
#define opt_Z  0x10000	/* -Z		additional zero		       */

struct pzrep
{
    std::complex<double> poles[MAXPZ], zeros[MAXPZ];
    int numpoles, numzeros;
};

    /* H(z) =\frac{Y(z)}{X(z)}= \frac{\sum_{k=0}^{N}b_k z^{-k}}{\sum_{k=0}^{M}a_k z^{-k}} */
    /* y(n)=b_0x(n)+...b_N x(n-N)-a_1 y(n-1)-...-a_M y(n-M) */
class IIRFilter 
{
public :
    pzrep splane, zplane;    
    int order;
    double raw_alpha1, raw_alpha2, raw_alphaz;
    std::complex<double>  dc_gain, fc_gain, hf_gain;
    unsigned int opts;	/* option flag bits */

    double warped_alpha1, warped_alpha2, chebrip, qfactor;
    unsigned int polemask;
    bool optsok;
    bool infq;

//    std::complex<double>  poles[MAXPOLES], zeros[MAXPOLES];
    double xcoeffs[MAXPOLES+1], ycoeffs[MAXPOLES+1];
    std::vector <double> a; // denominator  
    std::vector <double> b; // Numerator 
    int n;                  // Filter order y= bx-ay Y=B/A

    double fLow,fHigh;

    IIRFilter(std::string filterType,int order,double fs,std::vector<double> frequency);

std::complex<double>   eval(std::complex<double>   coeffs[], int np, std::complex<double>   z);

std::complex<double>  evaluate(std::complex<double>  topco[],std::complex<double>  botco[], int np, std::complex<double>  z);
void compute_s();
void choosepole(std::complex<double>);
void prewarp();
void normalize();
void compute_z_blt();
std::complex<double> blt(std::complex<double>);
void compute_z_mzt();
void compute_notch();
void compute_apres();
std::complex<double> reflect(std::complex<double>);
void compute_bpres();
void add_extra_zero();
void expandpoly();
void expand(std::complex<double>[], int, std::complex<double>[]);
void multin(std::complex<double>, int, std::complex<double>[]);
std::complex<double> evaluate(std::complex<double> topco[], int nz, std::complex<double> botco[], int np, std::complex<double> z);
void printrecurrence(); /* given (real) Z-plane poles & zeros, compute & print recurrence relation */
void printfilter();


};


#endif
