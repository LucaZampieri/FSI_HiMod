#include <lifev/himod/util/CheckModalBasis1D.hpp>

namespace LifeV
{

CheckModalBasis1D::
CheckModalBasis1D (const modalbasis_ptrType& modalbasis ) : M_modalbasis ( modalbasis )
{
}

void CheckModalBasis1D::
VerifyOrthonormality() const
{
    std::cout << std::endl;
    std::cout << "Checking orthonormality... " << std::endl;
    std::cout << std::endl;
    std::cout << "In theory, the basis should be orthonormal. So when you try to compute " << std::endl;
    std::cout << "the integral on the slice of Phi_j * Phi_i. You should obtain a Kronecker Delta as a result." << std::endl;
    std::cout << "The result depend a lot on the order of the quadrature rule, if the number of modes increases you will have problems." << std::endl;
    std::cout << "In addition, consider that problems arise when the maximum number between all the subindeces p,q is big." << std::endl;
    std::cout << "So this can happen when you have a big value of mtot or when you have a domain with Ly>>Lz o or" << std::endl;
    std::cout << "when you use HiMod to solve a 2D problem with big values of mtot." << std::endl;
    std::cout << std::endl;

    UInt m = M_modalbasis->mtot();
    Real err (0);
    std::vector<Real> delta (m * m, 0.0);
    for (UInt j (0); j < m; ++j)
    {
        for (UInt i (0); i < m; ++i)
        {
            delta[i + m * j] = M_modalbasis->compute_PhiPhi (i, j);
            //std::cout<<i<<'\t'<<j<<'\t'<<delta[i+m*j]<<"  p_i "<<M_modalbasis->eigenvalues(i).wp<<" q_i "<<M_modalbasis->eigenvalues(i).wq<<std::endl;

            if (i == j)
            {
                err += std::abs (delta[i + m * j] - 1);
            }
            else
            {
                err += std::abs (delta[i + m * j]);
            }
        }
    }
    std::cout << "Now we compute the integral and compare it with the kronecker delta, obtaining a matrix that should be " << std::endl;
    std::cout << "equal to an identity matrix, we compute the difference and, then, the sum of the absolut value of the result" << std::endl;
    std::cout << "The result should be close to zero." << std::endl;
    std::cout << "We also provide the maximum value between all the p and the q." << std::endl;
    std::cout << std::endl;
    std::cout << "The sum of the absolute value is: " << err << " maximum sub-frequency: " << M_modalbasis->subMMax() << std::endl;
    std::cout << std::endl;

}

void CheckModalBasis1D::
VerifyBC (const Real& mu, const Real& chi) const
{

    std::cout << "We want to verify that the basis satisfies the required boundary conditions so we report the value on the boundary." << std::endl;
    UInt nz = M_modalbasis->qrZ().nbQuadPt() - 1;
    UInt ny = M_modalbasis->qrY().nbQuadPt() - 1;

    std::cout << "Here we put the coefficient, pay attention that they may not make sense for a Dirichlet or Neumann BC." << std::endl;
    std::cout << "mu <- " << mu << std::endl;
    std::cout << "chi <- " << chi << std::endl;

    std::cout << "If you use a quadrature rule that doesn't have a node in the boundary" << std::endl;
    std::cout << "the results may be inaccurate, so we report here the position of the point where the basis is actually evaluated" << std::endl;

    std::cout << "First quadrature node in y <- ";
    std::cout << M_modalbasis->qrY().quadPointCoor (0, 0) << std::endl;
    std::cout << "First quadrature node in z <- ";
    std::cout << M_modalbasis->qrZ().quadPointCoor (0, 0) << std::endl;
    std::cout << "Last quadrature node in y <- ";
    std::cout << M_modalbasis->qrY().quadPointCoor (ny, 0) << std::endl;
    std::cout << "Last quadrature node in z <- ";
    std::cout << M_modalbasis->qrZ().quadPointCoor (nz, 0) << std::endl;

    for (UInt k (0); k < M_modalbasis->mtot(); ++k)
    {
        UInt p = M_modalbasis->eigenvalues (k).p - 1;
        UInt q = M_modalbasis->eigenvalues (k).q - 1;

        std::cout << "###########################" << std::endl;
        std::cout << " m <- " << k + 1 << " p <- " << p << " q <- " << q << std::endl;
        std::cout << " Dirichlet  " << std::endl;
        std::cout << "1 (down)" << '\t' << "<- ";
        std::cout << chi* M_modalbasis->phiz (q, 0) << std::endl;
        std::cout << "2 (right)" << '\t' << "<- ";
        std::cout << chi* M_modalbasis->phiy (p, ny) << std::endl;
        std::cout << "3 (up)   " << '\t' << "<- ";
        std::cout << chi* M_modalbasis->phiz (q, nz) << std::endl;
        std::cout << "4 (left)" << '\t' << "<- ";
        std::cout << chi* M_modalbasis->phiy (p, 0) << std::endl;
        std::cout << " Neumann " << std::endl;
        std::cout << "1 (down)" << '\t' << "<- ";
        std::cout << -mu* M_modalbasis->dphiz (q, 0) << std::endl;
        std::cout << "2 (right)" << '\t' << "<- ";
        std::cout << mu* M_modalbasis->dphiy (p, ny) << std::endl;
        std::cout << "3 (up)   " << '\t' << "<- ";
        std::cout << mu* M_modalbasis->dphiz (q, nz) << std::endl;
        std::cout << "4 (left)" << '\t' << "<- ";
        std::cout << -mu* M_modalbasis->dphiy (p, 0) << std::endl;
        std::cout << " Robin " << std::endl;
        std::cout << "1 (down)" << '\t' << "<- ";
        std::cout << chi* M_modalbasis->phiz (q, 0) - mu* M_modalbasis->dphiz (q, 0) << std::endl;
        std::cout << "2 (right)" << '\t' << "<- ";
        std::cout << chi* M_modalbasis->phiy (p, ny) + mu* M_modalbasis->dphiy (p, ny) << std::endl;
        std::cout << "3 (up)   " << '\t' << "<- ";
        std::cout << chi* M_modalbasis->phiz (q, nz) + mu* M_modalbasis->dphiz (q, nz) << std::endl;
        std::cout << "4 (left)" << '\t' << "<- ";
        std::cout << chi* M_modalbasis->phiy (p, 0) - mu* M_modalbasis->dphiy (p, 0) << std::endl;
    }

}
} //End lifev namespace
