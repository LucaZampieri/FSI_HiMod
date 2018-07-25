#include <vector>
#include <cmath>
#include <iostream>
#include <lifev/navier_stokes/function/bessel/bessel.hpp>

namespace bessel
{
void bjndd(unsigned int n,double x, std::vector<double>& bj, std::vector<double>& dj, std::vector<double>& fj)
//     =====================================================
//     Purpose: Compute Bessel functions Jn(x)and their
//     first and second derivatives(n= 0,1,úúú)
//     Input:   x ---  Argument of Jn(x,x ò 0)
//     n ---  Order of Jn(x)
//     Output:  BJ(n+1)---  Jn(x)
//     DJ(n+1)---  Jn'(x)
//     FJ(n+1)---  Jn"(x)
//     =====================================================
{
    int mt;
    double bs, f, f0, f1;
    unsigned int nt, m;

    for(nt=1; nt<=900; ++nt)
    {
        mt=static_cast<int>(0.5*log10(6.28*nt)-nt*log10(1.36*fabs(x)/nt));
        if(mt > 20)
            break;
    }

    m=nt;

    bs=0.0E0;
    f0=0.0E0;
    f1=1.0E-35;
    for (int k=m; k>=0; --k)
    {
        f=2.0E0*(k+1.0E0)*f1/x-f0;
        if(k <= n)bj[k]=f;
        if(k == 2*static_cast<int>(k/2)) bs+=2.0E0*f;
        f0=f1;
        f1=f;
    }
//    k=0-1;
    for (unsigned int k=0; k<=n; ++k)
    {
        bj[k]=bj[k]/(bs-f);
    }
//    k=static_cast<int>(n)+1;
    dj[0]=-bj[1];
    fj[0]=-1.0E0*bj[0]-dj[0]/x;
    for (unsigned int k=1; k<=n; ++k)
    {
        dj[k]=bj[k-1]-k*bj[k]/x;
        fj[k]=(k*k/(x*x)-1.0E0)*bj[k]-dj[k]/x;
    }
//    k=static_cast<int>(n)+1;
    return;
}

//       ===========================================================
//       Purpose: Compute the zeros of Bessel functions Jn(x) and
//                Jn'(x), and arrange them in the order of their
//                magnitudes
//       Input :  mbessel - positive integer describing how many bessel functions to evaluate
//                (always start with 'J0') - number of columns
//          nzeros - how many zeros to calculate - number of rows
//       Output:  zjm - this is the matrix of zeros of Jm(x)
//          zjmprime - this is the matrix of zeros of Jm'(x)
//       =============================================================
// function mjdzo
// This program is a direct conversion of the corresponding Fortran program in
// S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).

int besselzeros(unsigned int const mbessel, unsigned int const nzeros, std::vector<std::vector<double> >& zjm, std::vector<std::vector<double> >& zjmprime)
{
    unsigned int dim = 101;
    std::vector<double> bj(dim,0), dj(dim,0), fj(dim,0);
    double x=0.0, x0, x1, x2;

    unsigned int nm = mbessel; // columns
    unsigned int mm = nzeros;  // rows

    zjm.resize(mm);
    zjm.assign(mm, std::vector<double>(nm,0)); // each row is a vector containing nm 0s

    zjmprime.resize(mm);
    zjmprime.assign(mm, std::vector<double>(nm,0)); // each row is a vector containing nm 0s

    for(unsigned int  i=1; i<=nm; ++i)
    {
        x1=.407658+.4795504*sqrt(i-1)+.983618*(i-1);
        x2=1.99535+.8333883*sqrt(i-1)+.984584*(i-1);
        for  (unsigned int j=1; j<=mm; ++j)
        {
            unsigned int ifoo=1;
            if(!(i == 1 & j == 1))
            {
                x=x1;
                while (1)
                {
                    bjndd(i,x,bj,dj,fj);
                    x0=x;
                    x-=dj[i-1]/fj[i-1];
                    if(!(fabs(x-x0)> 1.0E-10))
                        break;
                }
            }
            if(ifoo == 1)
            {            
                zjmprime[j-1][i-1] = x;
                
                if(i <= 15)
                    x1=x+3.057+.0122*(i-1)+(1.555+.41575*(i-1))/((j+1)*(j+1));
                else
                    x1=x+2.918+.01924*(i-1)+(6.26+.13205*(i-1))/((j+1)*(j+1));
            }

            x=x2;
            ifoo=1;
            while (1)
            {
                bjndd(i,x,bj,dj,fj);
                x0=x;
                x-=bj[i-1]/dj[i-1];
                if(!(fabs(x-x0)> 1.0E-10))
                    break;
            }
            if(ifoo == 1)
            {
                zjm[j-1][i-1] = x;
                
                if(i <= 15)
                    x2=x+3.11+.0138*(i-1)+(.04832+.2804*(i-1))/((j+1)*(j+1));
                else
                    x2=x+3.001+.0105*(i-1)+(11.52+.48525*(i-1))/((j+3)*(j+3));
            }
        } // end j-for
    } // end i-for

    return 0;
}

} // end namespace bessel
