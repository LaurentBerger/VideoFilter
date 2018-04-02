#include "IIRFilter.hpp"
#include <opencv2/core/core.hpp>

/* mkfilter -- given n, compute recurrence relation
   to implement Butterworth, Bessel or Chebyshev filter of order n
   A.J. Fisher, University of York   <fisher@minster.york.ac.uk>
   September 1992 */

/* Main module */

#include <stdio.h>
#include <cmath>
/* mkfilter -- given n, compute recurrence relation
   to implement Butterworth or Bessel filter of order n
   A.J. Fisher, University of York   <fisher@minster.york.ac.uk>
   September 1992 */

/* Header file */


/*
#define uint	    unsigned int
#define bool	    int
#define true	    1
#define false	    0
#define word	    int		    
*/
#define EPS	    1e-10

#define seq(s1,s2)  (strcmp(s1,s2) == 0)
#define unless(x)   if(!(x))
#define until(x)    while(!(x))
#define ifix(x)	    (int) (((x) >= 0.0) ? (x) + 0.5 : (x) - 0.5)


using namespace std;



IIRFilter::IIRFilter(string ftype, int o, double fs, vector<double> frequency)
{
    // Frequency Bande passante où le gain du filtre est de 1 : 
    // 0 f1 passebas
    // f1 f2 passebande
    // f1 fs/2 passe haut
    if (frequency.size() != 2)
    {
        cv::Exception e;
        e.code = -1;
        e.msg = "you must give fLow and fHigh";

        throw e;

    }
    polemask = ~0;
    order = o;
    opts = 0;
    raw_alpha1 = frequency[0]/fs;
    raw_alpha2 = frequency[1]/fs;;
    if (abs(raw_alpha2-0.5)<0.000001)
        opts |=  opt_hp; 
    else if (raw_alpha1 == 0)
    {
        opts |=  opt_lp; 
        raw_alpha1 = raw_alpha2;
    }
    else
        opts |=  opt_bp; 
    if (ftype == "butt")
        opts = opts | opt_bu;
    else if (ftype == "cheb")
        opts = opts | opt_ch;
    else
    {
        cv::Exception e;
        e.code = -2;
        e.msg = "only butterworth and chebychef filter are implemented";
    }

    if (opts & opt_re)
    {
        if (opts & opt_bp) compute_bpres();	   /* bandpass resonator	 */
        if (opts & opt_bs) compute_notch();	   /* bandstop resonator (notch) */
        if (opts & opt_ap) compute_apres();	   /* allpass resonator		 */
    }
    else
    {
        if (opts & opt_pi)
        {
            prewarp();
            splane.poles[0] = 0.0;
            splane.zeros[0] = -TWOPI * warped_alpha1;
            splane.numpoles = splane.numzeros = 1;
        }
        else
        {
            compute_s();
            prewarp();
            normalize();
        }
        if (opts & opt_z) compute_z_mzt(); else compute_z_blt();
    }
    if (opts & opt_Z) add_extra_zero();
        expandpoly();
//    printresults(argv);
    b.resize(zplane.numzeros+1);
    a.resize(zplane.numpoles+1);
    for (int i=0; i < zplane.numzeros +1; i++)
    { 
        if (i > 0) printf("     + ");
        printf("(%3g * x[n-%2d])\n", xcoeffs[i], zplane.numzeros -i);
        b[zplane.numzeros-i]= xcoeffs[i]/abs(fc_gain);
    }
    putchar('\n');
    for (int i=0; i < zplane.numpoles; i++)
    { 
        printf("     + (%14.10f * y[n-%2d])\n", ycoeffs[i], zplane.numpoles-i);
        a[zplane.numpoles-i]=ycoeffs[i];
    }
    a[0]=0;
    putchar('\n');


}




complex<double>   IIRFilter::eval(complex<double>   coeffs[], int np, complex<double>   z)
  { /* evaluate polynomial in z, substituting for z */
    complex<double>   sum(0,0); int i;
    for (i=np; i >= 0; i--) sum = sum* z+ coeffs[i];
    return sum;
  }

complex<double>  IIRFilter::evaluate(complex<double>  topco[],complex<double>  botco[], int np, complex<double>  z)
  { /* evaluate response, substituting for z */
    return eval(topco, np, z)/ eval(botco, np, z);
  }



complex<double>  bessel_poles[] = {
complex<double>{ -1.00000000000e+00, 0.00000000000e+00 }, complex<double>{ -1.10160133059e+00, 6.36009824757e-01 },
complex<double>{ -1.32267579991e+00, 0.00000000000e+00 }, complex<double>{ -1.04740916101e+00, 9.99264436281e-01 },
complex<double>{ -1.37006783055e+00, 4.10249717494e-01 }, complex<double>{ -9.95208764350e-01, 1.25710573945e+00 },
complex<double>{ -1.50231627145e+00, 0.00000000000e+00 }, complex<double>{ -1.38087732586e+00, 7.17909587627e-01 },
complex<double>{ -9.57676548563e-01, 1.47112432073e+00 }, complex<double>{ -1.57149040362e+00, 3.20896374221e-01 },
complex<double>{ -1.38185809760e+00, 9.71471890712e-01 }, complex<double>{ -9.30656522947e-01, 1.66186326894e+00 },
complex<double>{ -1.68436817927e+00, 0.00000000000e+00 }, complex<double>{ -1.61203876622e+00, 5.89244506931e-01 },
complex<double>{ -1.37890321680e+00, 1.19156677780e+00 }, complex<double>{ -9.09867780623e-01, 1.83645135304e+00 },
complex<double>{ -1.75740840040e+00, 2.72867575103e-01 }, complex<double>{ -1.63693941813e+00, 8.22795625139e-01 },
complex<double>{ -1.37384121764e+00, 1.38835657588e+00 }, complex<double>{ -8.92869718847e-01, 1.99832584364e+00 },
complex<double>{ -1.85660050123e+00, 0.00000000000e+00 }, complex<double>{ -1.80717053496e+00, 5.12383730575e-01 },
complex<double>{ -1.65239648458e+00, 1.03138956698e+00 }, complex<double>{ -1.36758830979e+00, 1.56773371224e+00 },
complex<double>{ -8.78399276161e-01, 2.14980052431e+00 }, complex<double>{ -1.92761969145e+00, 2.41623471082e-01 },
complex<double>{ -1.84219624443e+00, 7.27257597722e-01 }, complex<double>{ -1.66181024140e+00, 1.22110021857e+00 },
complex<double>{ -1.36069227838e+00, 1.73350574267e+00 }, complex<double>{ -8.65756901707e-01, 2.29260483098e+00 } };

complex<double>  cmone( -1.0, 0.0 );
complex<double>  czero = {	 0.0, 0.0 };
complex<double>  cone  = {	 1.0, 0.0 };
complex<double>  ctwo  = {	 2.0, 0.0 };
complex<double>  chalf = {	 0.5, 0.0 };




#define cneg(z) csub(czero, z)



void IIRFilter::compute_s() /* compute S-plane poles for prototype LP filter */
{
    splane.numpoles = 0;
    if (opts & opt_be)
    { /* Bessel filter */
        int p = (order*order) / 4; /* ptr into table */
        if (order & 1) choosepole(bessel_poles[p++]);
        for (int i = 0; i < order / 2; i++)
        {
            choosepole(bessel_poles[p]);
            choosepole(conj(bessel_poles[p]));
            p++;
        }
    }
    if (opts & (opt_bu | opt_ch))
    { /* Butterworth filter */
        for (int i = 0; i < 2 * order; i++)
        {
            double theta = (order & 1) ? (i*PI) / order : ((i + 0.5)*PI) / order;
            choosepole(complex<double>(cos(theta), sin(theta)));
        }
    }
    if (opts & opt_ch)
    { /* modify for Chebyshev (p. 136 DeFatta et al.) */
        if (chebrip >= 0.0)
        {
            fprintf(stderr, "mkfilter: Chebyshev ripple is %g dB; must be .lt. 0.0\n", chebrip);
            exit(1);
        }
        double rip = pow(10.0, -chebrip / 10.0);
        double eps = sqrt(rip - 1.0);
        double y = asinh(1.0 / eps) / (double)order;
        if (y <= 0.0)
        {
            fprintf(stderr, "mkfilter: bug: Chebyshev y=%g; must be .gt. 0.0\n", y);
            exit(1);
        }
        for (int i = 0; i < splane.numpoles; i++)
        {
            splane.poles[i] = complex<double>(splane.poles[i].real() * sinh(y),splane.poles[i].imag() * cosh(y));
        }
    }
}

void IIRFilter::choosepole(complex<double> z)
{
    if (z.real() < 0.0)
    {
        if (polemask & 1) splane.poles[splane.numpoles++] = z;
        polemask >>= 1;
    }
}

void IIRFilter::prewarp()
{ /* for bilinear transform, perform pre-warp on alpha values */
    if (opts & (opt_w | opt_z))
    {
        warped_alpha1 = raw_alpha1;
        warped_alpha2 = raw_alpha2;
    }
    else
    {
        warped_alpha1 = tan(PI * raw_alpha1) / PI;
        warped_alpha2 = tan(PI * raw_alpha2) / PI;
    }
}

void IIRFilter::normalize()		/* called for trad, not for -Re or -Pi */
{
    double w1 = TWOPI * warped_alpha1;
    double w2 = TWOPI * warped_alpha2;
    /* transform prototype into appropriate filter type (lp/hp/bp/bs) */
    switch (opts & (opt_lp | opt_hp | opt_bp | opt_bs))
    {
    case opt_lp:
    { for (int i = 0; i < splane.numpoles; i++) splane.poles[i] = splane.poles[i] * w1;
    splane.numzeros = 0;
    break;
    }

    case opt_hp:
    { int i;
    for (i = 0; i < splane.numpoles; i++) splane.poles[i] = w1 / splane.poles[i];
    for (i = 0; i < splane.numpoles; i++) splane.zeros[i] = 0.0;	 /* also N zeros at (0,0) */
    splane.numzeros = splane.numpoles;
    break;
    }

    case opt_bp:
    { double w0 = sqrt(w1*w2), bw = w2 - w1; int i;
    for (i = 0; i < splane.numpoles; i++)
    {
        complex<double> hba = 0.5 * (splane.poles[i] * bw);
        complex<double> temp = sqrt(1.0 - (w0*w0 / (hba*hba)));
        splane.poles[i] = hba * (1.0 + temp);
        splane.poles[splane.numpoles + i] = hba * (1.0 - temp);
    }
    for (i = 0; i < splane.numpoles; i++) splane.zeros[i] = 0.0;	 /* also N zeros at (0,0) */
    splane.numzeros = splane.numpoles;
    splane.numpoles *= 2;
    break;
    }

    case opt_bs:
    { double w0 = sqrt(w1*w2), bw = w2 - w1; int i;
    for (i = 0; i < splane.numpoles; i++)
    {
        complex<double> hba = 0.5 * (bw / splane.poles[i]);
        complex<double> temp = sqrt(1.0 - (w0*w0 / (hba*hba)));
        splane.poles[i] = hba * (1.0 + temp);
        splane.poles[splane.numpoles + i] = hba * (1.0 - temp);
    }
    for (i = 0; i < splane.numpoles; i++)	   /* also 2N zeros at (0, +-w0) */
    {
        splane.zeros[i] = complex<double>(0.0, +w0);
        splane.zeros[splane.numpoles + i] = complex<double>(0.0, -w0);
    }
    splane.numpoles *= 2;
    splane.numzeros = splane.numpoles;
    break;
    }
    }
}

void IIRFilter::compute_z_blt() /* given S-plane poles & zeros, compute Z-plane poles & zeros, by bilinear transform */
{
    int i;
    zplane.numpoles = splane.numpoles;
    zplane.numzeros = splane.numzeros;
    for (i = 0; i < zplane.numpoles; i++) 
        zplane.poles[i] = blt(splane.poles[i]);
    for (i = 0; i < zplane.numzeros; i++) 
        zplane.zeros[i] = blt(splane.zeros[i]);
    while (zplane.numzeros < zplane.numpoles) 
        zplane.zeros[zplane.numzeros++] = -1.0;
}

complex<double> IIRFilter::blt(complex<double> pz)
{
    return (2.0 + pz) / (2.0 - pz);
}

void IIRFilter::compute_z_mzt() /* given S-plane poles & zeros, compute Z-plane poles & zeros, by matched z-transform */
{
    int i;
    zplane.numpoles = splane.numpoles;
    zplane.numzeros = splane.numzeros;
    for (i = 0; i < zplane.numpoles; i++) zplane.poles[i] = exp(splane.poles[i]);
    for (i = 0; i < zplane.numzeros; i++) zplane.zeros[i] = exp(splane.zeros[i]);
}

void IIRFilter::compute_notch()
{ /* compute Z-plane pole & zero positions for bandstop resonator (notch filter) */
    compute_bpres();		/* iterate to place poles */
    double theta = TWOPI * raw_alpha1;
    complex<double> zz = complex<double>(cos(theta),sin(theta));	/* place zeros exactly */
    zplane.zeros[0] = zz; zplane.zeros[1] = conj(zz);
}

void IIRFilter::compute_apres()
{ /* compute Z-plane pole & zero positions for allpass resonator */
    compute_bpres();		/* iterate to place poles */
    zplane.zeros[0] = reflect(zplane.poles[0]);
    zplane.zeros[1] = reflect(zplane.poles[1]);
}

complex<double> IIRFilter::reflect(complex<double> z)
{
    double r = hypot(z.real(),z.imag());
    return z / (r*r);
}

void IIRFilter::compute_bpres()
{ /* compute Z-plane pole & zero positions for bandpass resonator */
    zplane.numpoles = zplane.numzeros = 2;
    zplane.zeros[0] = 1.0; zplane.zeros[1] = -1.0;
    double theta = TWOPI * raw_alpha1; /* where we want the peak to be */
    if (infq)
    { /* oscillator */
        complex<double> zp = complex<double>(cos(theta), sin(theta));
        zplane.poles[0] = zp; zplane.poles[1] = conj(zp);
    }
    else
    { /* must iterate to find exact pole positions */
        complex<double> topcoeffs[MAXPZ + 1]; expand(zplane.zeros, zplane.numzeros, topcoeffs);
        double r = exp(-theta / (2.0 * qfactor));
        double thm = theta, th1 = 0.0, th2 = PI;
        bool cvg = false;
        for (int i = 0; i < 50 && !cvg; i++)
        {
            complex<double> zp = r * complex<double>(cos(thm), sin(thm));
            zplane.poles[0] = zp; zplane.poles[1] = conj(zp);
            complex<double> botcoeffs[MAXPZ + 1]; expand(zplane.poles, zplane.numpoles, botcoeffs);
            complex<double> g = evaluate(topcoeffs, zplane.numzeros, botcoeffs, zplane.numpoles, complex<double>(cos(theta), sin(theta)));
            double phi = g.imag() / g.real(); /* approx to atan2 */
            if (phi > 0.0) th2 = thm; else th1 = thm;
            if (fabs(phi) < EPS) cvg = true;
            thm = 0.5 * (th1 + th2);
        }
        unless(cvg) fprintf(stderr, "mkfilter: warning: failed to converge\n");
    }
}

void IIRFilter::add_extra_zero()
{
    if (zplane.numzeros + 2 > MAXPZ)
    {
        fprintf(stderr, "mkfilter: too many zeros; can't do -Z\n");
        exit(1);
    }
    double theta = TWOPI * raw_alphaz;
    complex<double> zz(cos(theta), sin(theta));
    zplane.zeros[zplane.numzeros++] = zz;
    zplane.zeros[zplane.numzeros++] = conj(zz);
    while (zplane.numpoles < zplane.numzeros) zplane.poles[zplane.numpoles++] = 0.0;	 /* ensure causality */
}

void IIRFilter::expandpoly() /* given Z-plane poles & zeros, compute top & bot polynomials in Z, and then recurrence relation */
{
    complex<double> topcoeffs[MAXPZ + 1], botcoeffs[MAXPZ + 1]; int i;
    expand(zplane.zeros, zplane.numzeros, topcoeffs);
    expand(zplane.poles, zplane.numpoles, botcoeffs);
    dc_gain = evaluate(topcoeffs, zplane.numzeros, botcoeffs, zplane.numpoles, 1.0);
    double theta = TWOPI * 0.5 * (raw_alpha1 + raw_alpha2); /* "jwT" for centre freq. */
    fc_gain = evaluate(topcoeffs, zplane.numzeros, botcoeffs, zplane.numpoles, complex<double>(cos(theta), sin(theta)));
    hf_gain = evaluate(topcoeffs, zplane.numzeros, botcoeffs, zplane.numpoles, -1.0);
    for (i = 0; i <= zplane.numzeros; i++) xcoeffs[i] = +(topcoeffs[i].real() / botcoeffs[zplane.numpoles].real());
    for (i = 0; i <= zplane.numpoles; i++) ycoeffs[i] = -(botcoeffs[i].real() / botcoeffs[zplane.numpoles].real());
}

void IIRFilter::expand(complex<double> pz[], int npz, complex<double> coeffs[])
{ /* compute product of poles or zeros as a polynomial of z */
    int i;
    coeffs[0] = 1.0;
    for (i = 0; i < npz; i++) coeffs[i + 1] = 0.0;
    for (i = 0; i < npz; i++) multin(pz[i], npz, coeffs);
    /* check computed coeffs of z^k are all real */
    for (i = 0; i < npz + 1; i++)
    {
        if (fabs(coeffs[i].imag()) > EPS)
        {
            fprintf(stderr, "mkfilter: coeff of z^%d is not real; poles/zeros are not complex conjugates\n", i);
            exit(1);
        }
    }
}

void IIRFilter::multin(complex<double> w, int npz, complex<double> coeffs[])
{ /* multiply factor (z-w) into coeffs */
    complex<double> nw = -w;
    for (int i = npz; i >= 1; i--) coeffs[i] = (nw * coeffs[i]) + coeffs[i - 1];
    coeffs[0] = nw * coeffs[0];
}

complex<double> IIRFilter::evaluate(complex<double> topco[], int nz, complex<double> botco[], int np, complex<double> z)
{ /* evaluate response, substituting for z */
    return eval(topco, nz, z) / eval(botco, np, z);
}




void IIRFilter::printrecurrence() /* given (real) Z-plane poles & zeros, compute & print recurrence relation */
  { int i;
    printf("Recurrence relation:\n");
    printf("y[n] = ");
    for (i=0; i < splane.numpoles+1; i++)
      { if (i > 0) printf("     + ");
	printf("(%3g * x[n-%2d])\n", xcoeffs[i], splane.numpoles-i);
      }
    putchar('\n');
    for (i=0; i < splane.numpoles; i++)
      { printf("     + (%14.10f * y[n-%2d])\n", ycoeffs[i], splane.numpoles-i);
      }
    putchar('\n');
  }


void IIRFilter::printfilter()
{ /*
    printf("raw alpha1    = %14.10f\n", raw_alpha1);
    printf("warped alpha1 = %14.10f\n", warped_alpha1);
    printf("raw alpha2    = %14.10f\n", raw_alpha2);
    printf("warped alpha2 = %14.10f\n", warped_alpha2);
    printgain("dc    ", dc_gain);
    printgain("centre", fc_gain);
    printgain("hf    ", hf_gain);
    putchar('\n');
    printrat_s();
    printrat_z();
    printrecurrence();*/
  }









  /*int main(int argc, char **argv) 
  { 
    compute_s();
    normalize();
    compute_z();
    expandpoly();
    printrecurrence();


    return(0);
  }*/
