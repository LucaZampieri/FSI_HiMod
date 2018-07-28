#include <lifev/himod/modalbasis/ModalSpaceRectangular.hpp>

namespace LifeV
{

bool Comparison::
operator() (const EigenMap& a, const EigenMap& b) const
{
    if ( ( pow (a.wp, 2) + pow (a.wq, 2) ) != ( pow (b.wp, 2) + pow (b.wq, 2) ) )
    {
        return (pow (a.wp, 2) + pow (a.wq, 2) ) < (pow (b.wp, 2) + pow (b.wq, 2) );
    }
    else
    {
        return a.wp < b.wp ;
    }
}


//Fourier Coefficients in the case g is x-independent
std::vector<Real> ModalSpaceRectangular::
fourierCoefficients( const function_Type& g ) const
{
    // Initialization of FourCoeff as 0, length M_mtot
    std::vector<Real> FourCoeff (M_mtot, 0.0);
    // Loop to evaluate g on the quadrature nodes
    // a matrix 32x32 (NnodesXNnodes)
    std::vector<std::vector<Real> > evaluate_g;

    // first resize
    evaluate_g.resize (M_quadruleY->nbQuadPt() );

    // second resize
    for (UInt k = 0; k < evaluate_g.size(); ++k)
    {
        evaluate_g[k].resize (M_quadruleZ->nbQuadPt() );
    }
    // evaluation of g
    for (UInt n = 0; n < M_quadruleY->nbQuadPt(); ++n)
        for (UInt m = 0; m < M_quadruleZ->nbQuadPt(); ++m)
        {
            evaluate_g[n][m] = g ( 0 , 0 , M_quadruleY->quadPointCoor (n, 0) * M_Ly , M_quadruleZ->quadPointCoor (m, 0) * M_Lz , 0);
        }
    // Loop on the modes, for each mode compute the associated fourier coefficients
    // \f[ \int_{[0,Ly]\times[0,Lz]}  g(y,z) \phi_k dy dz \f]
    for (UInt k = 0; k < M_mtot; ++k)
    {
        // extraction of the sub-indeces
        UInt p_k = (M_eigenvalues[k].p - 1); //y
        UInt q_k = (M_eigenvalues[k].q - 1); //z

        Real normy = 1.0 / std::sqrt (M_Ly);
        Real normz = 1.0 / std::sqrt (M_Lz);

        //loop over quadrature nodes
        for (UInt n = 0; n < M_quadruleY->nbQuadPt(); ++n) //y
            for (UInt m = 0; m < M_quadruleZ->nbQuadPt(); ++m) //z
                FourCoeff[k] += evaluate_g[n][m]*                                           // function evaluated in the right nodes
                                M_phiy [p_k][ n] * normy *
                                M_phiz [q_k][m] * normz*                  // \f$ \phi_k\f$ evaluated on the nodes
                                M_Ly * M_Lz*                                                // Jacobian
                                M_quadruleY->weight (n) * M_quadruleZ->weight (m);            // weights

    }
    return FourCoeff;
}

//Single Fourier coefficient in the case f is x-dependent
Real ModalSpaceRectangular::
fourierCoeffPointWise (const Real& x, const function_Type& f, const UInt& k) const
{
    Real coeff = 0.0;
    UInt p_k = (M_eigenvalues[k].p - 1);
    UInt q_k = (M_eigenvalues[k].q - 1);

    Real normy = 1.0 / std::sqrt (M_Ly);
    Real normz = 1.0 / std::sqrt (M_Lz);

    for (UInt n = 0; n < M_quadruleY->nbQuadPt(); ++n)
        for (UInt m = 0; m < M_quadruleZ->nbQuadPt(); ++m)
            coeff += f ( 0 , x , M_quadruleY->quadPointCoor ( n , 0 ) * M_Ly , M_quadruleZ->quadPointCoor ( m , 0 ) * M_Lz , 0 ) *
                     M_phiy[p_k][n] * normy * M_phiz[q_k][m] * normz * M_Ly * M_Lz * M_quadruleY->weight (n) * M_quadruleZ->weight (m);
    return coeff;
}

//FindSubMMax
void ModalSpaceRectangular::
findSubMMax()
{
    if (M_subMMax == 0)
    {
        M_subMMax = 1;
        for (UInt i = 0; i < M_mtot; ++i)
        {
            M_subMMax = std::max (M_subMMax, std::max (M_eigenvalues[i].p, M_eigenvalues[i].q) );
        }
    }
    return;
}

//Compute the integral on modal basis
//PHI PHI
// Should return d(i,j) we should add this improvement
Real ModalSpaceRectangular::
compute_PhiPhi (const UInt& j, const UInt& k) const
{
    Real coeff_y = 0.0;
    Real coeff_z = 0.0;

    UInt p_j = (M_eigenvalues[j].p - 1);
    UInt p_k = (M_eigenvalues[k].p - 1);
    UInt q_j = (M_eigenvalues[j].q - 1);
    UInt q_k = (M_eigenvalues[k].q - 1);

    // If you have a function which is normalized in (0,1) in respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normy = 1.0 / std::sqrt (M_Ly);
    Real normz = 1.0 / std::sqrt (M_Lz);

    for (UInt n = 0; n < M_quadruleY->nbQuadPt(); ++n)
    {
        coeff_y += M_phiy [p_j][n] * normy *
                   M_phiy [p_k][n] * normy *
                   M_Ly * M_quadruleY->weight (n);
    }

    for (UInt n = 0; n < M_quadruleZ->nbQuadPt(); ++n)
    {
        coeff_z += M_phiz[q_j][n] * normz *
                   M_phiz[q_k][n] * normz *
                   M_Lz * M_quadruleZ->weight (n);
    }

    return coeff_y * coeff_z;
}

//DYPHI PHI
Real ModalSpaceRectangular::
compute_DyPhiPhi (const UInt& j, const UInt& k) const
{
    Real coeff_y = 0.0;
    Real coeff_z = 0.0;

    UInt p_j = (M_eigenvalues[j].p - 1);
    UInt p_k = (M_eigenvalues[k].p - 1);
    UInt q_j = (M_eigenvalues[j].q - 1);
    UInt q_k = (M_eigenvalues[k].q - 1);
    // If you have a function which is normalized in (0,1) in respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normy = 1.0 / std::sqrt (M_Ly);
    Real normz = 1.0 / std::sqrt (M_Lz);
    for (UInt n = 0; n < M_quadruleY->nbQuadPt(); ++n)
    {
        coeff_y += M_dphiy[p_j][n] * normy * M_phiy[p_k][n] * normy * M_quadruleY->weight (n);
    }

    for (UInt n = 0; n < M_quadruleZ->nbQuadPt(); ++n)
    {
        coeff_z += M_phiz[q_j][n] * normz * M_phiz[q_k][n] * normz * M_Lz * M_quadruleZ->weight (n);
    }

    return coeff_y * coeff_z;
}


//DYPHI DYPHI
Real ModalSpaceRectangular::
compute_DyPhiDyPhi (const UInt& j, const UInt& k) const
{
    Real coeff_y = 0;
    Real coeff_z = 0;

    UInt p_j = (M_eigenvalues[j].p - 1);
    UInt p_k = (M_eigenvalues[k].p - 1);
    UInt q_j = (M_eigenvalues[j].q - 1);
    UInt q_k = (M_eigenvalues[k].q - 1);
    // If you have a function which is normalized in (0,1) in respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normy = 1.0 / sqrt (M_Ly);
    Real normz = 1.0 / sqrt (M_Lz);

    for (UInt n = 0; n < M_quadruleY->nbQuadPt(); ++n)
        coeff_y +=  M_dphiy[p_j][n] * normy *
                    M_dphiy[p_k][n] * normy / M_Ly *
                    M_quadruleY->weight (n);

    for (UInt n = 0; n < M_quadruleZ->nbQuadPt(); ++n)
        coeff_z +=  M_phiz [q_j][n] * normz *
                    M_phiz [q_k][n] * normz * M_Lz *
                    M_quadruleZ->weight (n);

    return coeff_y * coeff_z;
}

//DZPHI PHI
Real ModalSpaceRectangular::
compute_DzPhiPhi (const UInt& j, const UInt& k) const
{
    Real coeff_y = 0;
    Real coeff_z = 0;

    UInt p_j = (M_eigenvalues[j].p - 1);
    UInt p_k = (M_eigenvalues[k].p - 1);
    UInt q_j = (M_eigenvalues[j].q - 1);
    UInt q_k = (M_eigenvalues[k].q - 1);
    // If you have a function which is normalized in (0,1) in respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normy = 1.0 / sqrt (M_Ly);
    Real normz = 1.0 / sqrt (M_Lz);
    for (UInt n = 0; n < M_quadruleY->nbQuadPt(); ++n)
    {
        coeff_y += M_phiy [p_j][n] * normy * M_phiy [p_k][n] * normy * M_Ly * M_quadruleY->weight (n);
    }

    for (UInt n = 0; n < M_quadruleZ->nbQuadPt(); ++n)
    {
        coeff_z += M_dphiz [q_j][n] * normz * M_phiz [q_k][n] * normz * M_quadruleZ->weight (n);
    }
    return coeff_y * coeff_z;
}


//DZPHI DZPHI
Real ModalSpaceRectangular::
compute_DzPhiDzPhi (const UInt& j, const UInt& k) const
{
    Real coeff_y = 0;
    Real coeff_z = 0;

    UInt p_j = (M_eigenvalues[j].p - 1);
    UInt p_k = (M_eigenvalues[k].p - 1);
    UInt q_j = (M_eigenvalues[j].q - 1);
    UInt q_k = (M_eigenvalues[k].q - 1);
    // If you have a function which is normalized in (0,1) in respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normy = 1.0 / std::sqrt (M_Ly);
    Real normz = 1.0 / std::sqrt (M_Lz);

    for (UInt n = 0; n < M_quadruleY->nbQuadPt(); ++n)
        coeff_y +=  M_phiy [p_j][n] * normy *
                    M_phiy [p_k][n] * normy *
                    M_Ly * M_quadruleY->weight (n);

    for (UInt n = 0; n < M_quadruleZ->nbQuadPt(); ++n)
        coeff_z +=  M_dphiz [q_j][n] * normz *
                    M_dphiz [q_k][n] * normz / M_Lz *
                    M_quadruleZ->weight (n);

    return coeff_y * coeff_z;
}

//PHI
Real ModalSpaceRectangular::
compute_Phi (const UInt& k) const
{
    Real coeff_y = 0;
    Real coeff_z = 0;

    UInt p_k = (M_eigenvalues[k].p - 1);
    UInt q_k = (M_eigenvalues[k].q - 1);
    // If you have a function which is normalized in (0,1) in respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normy = 1.0 / std::sqrt (M_Ly);
    Real normz = 1.0 / std::sqrt (M_Lz);
    for (UInt n = 0; n < M_quadruleY->nbQuadPt(); ++n)
    {
        coeff_y += M_phiy [p_k][n] * normy * M_Ly * M_quadruleY->weight (n);
    }

    for (UInt n = 0; n < M_quadruleZ->nbQuadPt(); ++n)
    {
        coeff_z += M_phiz [q_k][n] * normz * M_Lz * M_quadruleZ->weight (n);
    }
    return coeff_y * coeff_z;
}

//Boundary contribute
// Y
Real ModalSpaceRectangular::
compute_y_PhiPhi (const UInt& j, const UInt& k) const
{
    Real coeff_y = 0.0;

    UInt p_j = (M_eigenvalues[j].p - 1);
    UInt p_k = (M_eigenvalues[k].p - 1);

    Real normy = 1.0 / std::sqrt (M_Ly);

    for (UInt n = 0; n < M_quadruleY->nbQuadPt(); ++n)
    {
        coeff_y += M_phiy [p_j][n] * normy *
                   M_phiy [p_k][n] * normy *
                   M_Ly * M_quadruleY->weight (n);
    }

    return coeff_y;
}

// Z
Real ModalSpaceRectangular::
compute_z_PhiPhi (const UInt& j, const UInt& k) const
{
    Real coeff_z = 0.0;

    UInt q_j = (M_eigenvalues[j].q - 1);
    UInt q_k = (M_eigenvalues[k].q - 1);

    Real normz = 1.0 / std::sqrt (M_Lz);

    for (UInt n = 0; n < M_quadruleZ->nbQuadPt(); ++n)
    {
        coeff_z += M_phiz[q_j][n] * normz *
                   M_phiz[q_k][n] * normz *
                   M_Lz * M_quadruleZ->weight (n);
    }

    return coeff_z;
}




void ModalSpaceRectangular::
addSliceBCY (const std::string& left, const std::string& right, const Real& mu, const Real& chi)
{

    M_genbasisY = Basis1DFactory::instance().createObject (left + right);
    M_genbasisY->setL (M_Ly);
    M_genbasisY->setMu (mu);
    M_genbasisY->setChi (chi);

    return;
}


void ModalSpaceRectangular::
addSliceBCZ (const std::string& down, const std::string& up, const Real& mu, const Real& Chi)
{

    M_genbasisZ = Basis1DFactory::instance().createObject (down + up);
    M_genbasisZ->setL (M_Lz);
    M_genbasisZ->setMu (mu);
    M_genbasisZ->setChi (Chi);

    return;
}

void ModalSpaceRectangular::
eigensProvider()
{
    Real wpold;
    Real wqold;
    UInt p = 1;
    UInt q = 1;

    wpold = M_genbasisY->Next();
    wqold = M_genbasisZ->Next();

    //wp and wq are the square roots of the sub-eigenvalues
    M_eigenvalues.resize (M_mtot);
    M_eigenvalues[0].wp = wpold;
    M_eigenvalues[0].wq = wqold;
    M_eigenvalues[0].p = p;
    M_eigenvalues[0].q = q;
    M_eigenvaluesY.push_back (wpold);
    M_eigenvaluesZ.push_back (wqold);

    /*
         Choosing a set instead of a multiset prevent from creating twice or even more times the same eigenvalues and p q indeces.
         Note that it's impossibile to insert an undesired copy after having erased its mate.
     */
    std::set< EigenMap, Comparison > subeigen;

    Real wp = M_genbasisY->Next();
    Real wq = M_genbasisZ->Next();
    M_eigenvaluesY.push_back (wp);
    M_eigenvaluesZ.push_back (wq);
    subeigen.insert ( EigenMap::make_eigenmap (wp, wqold, p + 1, q ) );
    subeigen.insert ( EigenMap::make_eigenmap (wpold, wq, p, q + 1 ) );

    for (UInt j = 1; j < M_mtot; ++j)
    {
        // Estraggo quello da inserire nel vettore di uscita
        M_eigenvalues[j].wp = subeigen.begin()->wp;
        M_eigenvalues[j].wq = subeigen.begin()->wq;
        M_eigenvalues[j].p  = subeigen.begin()->p;
        M_eigenvalues[j].q  = subeigen.begin()->q;
        /*
        std::cout<<"inserted "<<M_eigenvalues[j].wp<<'\t';
        std::cout<<M_eigenvalues[j].wq<<'\t';
        std::cout<<M_eigenvalues[j].p<<'\t';
        std::cout<<M_eigenvalues[j].q<<'\t'<<std::endl;
        */
        //cancello il primo elemento della lista di supporto
        subeigen.erase (subeigen.begin() );

        Real wpnext = -1.;
        Real wqnext = -1.;
        UInt pnext = M_eigenvalues[j].p + 1;
        UInt qnext = M_eigenvalues[j].q + 1;


        // Cerco se p attuale c'è già stato prima nel vettore
        for ( Int k = j - 1; k >= 0; --k)
        {
            if (M_eigenvalues[k].p == pnext)
            {
                wpnext = M_eigenvalues[k].wp;
            }
            if (M_eigenvalues[k].q == qnext)
            {
                wqnext = M_eigenvalues[k].wq;
            }
        }

        // Cerco nel set
        std::set< EigenMap, Comparison >::iterator iter;
        if (wpnext + wqnext == -2.0)
        {
            for ( iter = subeigen.begin(); iter != subeigen.end(); ++iter)
            {
                ;
            }
            {
                if (iter->p == pnext)
                {
                    wpnext = iter->wp;
                }
                if (iter->q == qnext)
                {
                    wqnext = iter->wq;
                }
            }
        }

        bool check = true;

        if (wpnext == -1.)
            for ( iter = subeigen.begin(); iter != subeigen.end() && check; ++iter)
            {
                if (iter->p == pnext)
                {
                    wpnext = iter->wp;
                    check = false;
                }
            }

        check = true;
        if (wqnext == -1.)
            for ( iter = subeigen.begin(); iter != subeigen.end() && check; ++iter)
            {
                if (iter->q == qnext)
                {
                    wqnext = iter->wq;
                    check = false;
                }
            }


        // Infine calcolo se non gli ho trovati
        if (wpnext == -1.)
        {
            wpnext = M_genbasisY->Next();
            M_eigenvaluesY.push_back (wpnext);
        }
        if (wqnext == -1.)
        {
            wqnext = M_genbasisZ->Next();
            M_eigenvaluesZ.push_back (wqnext);
        }

        //Inserisco il nuovo branch nell'albero
        subeigen.insert ( EigenMap::make_eigenmap (wpnext, M_eigenvalues[j].wq, M_eigenvalues[j].p + 1, M_eigenvalues[j].q ) );
        subeigen.insert ( EigenMap::make_eigenmap (M_eigenvalues[j].wp, wqnext, M_eigenvalues[j].p, M_eigenvalues[j].q + 1 ) );


        /*
                std::cout<<"Set di supporto -----------------"<<std::endl;
                for(iter = subeigen.begin();iter!=subeigen.end();++iter)
                    std::cout<<iter->p<<'\t'<<iter->q<<std::endl;
        */
    }
    //Stampo il vettore e il set di supporto
    /* std::cout << "Vettore di uscita ---------------" << std::endl;
     for (UInt l = 0; l < M_eigenvaluesY.size(); ++l)
     {
         std::cout << M_eigenvaluesY[l] << std::endl;
     }*/

}
void ModalSpaceRectangular::
evaluateBasis()
{
    eigensProvider();
    findSubMMax();
    M_genbasisY->evaluateBasis (M_phiy, M_dphiy, M_eigenvaluesY, M_quadruleY);
    M_genbasisZ->evaluateBasis (M_phiz, M_dphiz, M_eigenvaluesZ, M_quadruleZ);
}
void ModalSpaceRectangular::
showMe() const
{
    std::cout << "---- MODAL SPACE SHOWME ---" << std::endl;
    /*
        std::cout << "Eigenvalues vector:" << std::endl;
        std::cout << "Lambda" << '\t' << '\t';
        std::cout << "p" << '\t' << "q" << std::endl;

        for (UInt j = 0; j < M_mtot; ++j)
        {
            std::cout << std::pow (M_eigenvalues[j].wp, 2) + std::pow (M_eigenvalues[j].wq, 2)  << '\t' << '\t';
            std::cout << M_eigenvalues[j].p << '\t' << M_eigenvalues[j].q << std::endl;
        }
    */
    std::cout << "Ly = "             << M_Ly      << std::endl;
    std::cout << "Lz = "             << M_Lz      << std::endl;
    std::cout << "M = "              << M_mtot    << std::endl;
    std::cout << "MaxSubfrequency = " << M_subMMax << std::endl;
    std::cout << "---------------------------" << std::endl;
}


Real
ModalSpaceRectangular::
compute_R11 (const UInt& j, const UInt& k, const function_Type& mu, const Real& x) const
{
    Real coeff = 0;


    UInt p_j = (M_eigenvalues[j].p - 1);
    UInt p_k = (M_eigenvalues[k].p - 1);
    UInt q_j = (M_eigenvalues[j].q - 1);
    UInt q_k = (M_eigenvalues[k].q - 1);

    Real normy = 1.0 / std::sqrt (M_Ly);
    Real normz = 1.0 / std::sqrt (M_Lz);
    for (UInt n = 0; n < M_quadruleY->nbQuadPt(); ++n)
    {
        for (UInt m = 0; m < M_quadruleZ->nbQuadPt(); ++m)
        {

            coeff += mu (0, x, M_quadruleY->quadPointCoor (n, 0) * M_Ly,
                         M_quadruleZ->quadPointCoor (m, 0) * M_Lz, 0) *
                     M_Ly * M_Lz *
                     M_phiy [p_j][n] * normy * M_phiy [p_k][n] * normy * M_quadruleY->weight (n) *
                     M_phiz [q_j][m] * normz * M_phiz [q_k][m] * normz * M_quadruleZ->weight (m);
        }
    }
    /*
    std::cout<<"X "<<x<<std::endl;
    std::cout<<"J "<<j<<std::endl;
    std::cout<<"K "<<k<<std::endl;
    std::cout<<"C "<<coeff<<std::endl;
    */
    return coeff;
}

Real
ModalSpaceRectangular::
compute_R10 (const UInt& j, const UInt& k, const function_Type& beta, const Real& x) const
{
    Real coeff = 0;


    UInt p_j = (M_eigenvalues[j].p - 1);
    UInt p_k = (M_eigenvalues[k].p - 1);
    UInt q_j = (M_eigenvalues[j].q - 1);
    UInt q_k = (M_eigenvalues[k].q - 1);
    // If you have a function which is normalized in (0,1) in respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normy = 1.0 / std::sqrt (M_Ly);
    Real normz = 1.0 / std::sqrt (M_Lz);
    for (UInt n = 0; n < M_quadruleY->nbQuadPt(); ++n)
    {
        for (UInt m = 0; m < M_quadruleZ->nbQuadPt(); ++m)
        {
            coeff += beta (0, x, M_quadruleY->quadPointCoor (n, 0) * M_Ly,
                           M_quadruleZ->quadPointCoor (m, 0) * M_Lz, 0) *
                     M_Ly * M_Lz *
                     M_phiy [p_j][n] * normy * M_phiy [p_k][n] * normy * M_quadruleY->weight (n) *
                     M_phiz [q_j][m] * normz * M_phiz [q_k][m] * normz * M_quadruleZ->weight (m);
        }
    }

    return coeff;
}


Real
ModalSpaceRectangular::
compute_R00 (const UInt& j, const UInt& k, const function_Type& mu, const function_Type& beta, const function_Type& sigma, const Real& x) const
{
    Real coeff = 0;

    UInt p_j = (M_eigenvalues[j].p - 1);
    UInt p_k = (M_eigenvalues[k].p - 1);
    UInt q_j = (M_eigenvalues[j].q - 1);
    UInt q_k = (M_eigenvalues[k].q - 1);

    Real normy = 1.0 / std::sqrt (M_Ly);
    Real normz = 1.0 / std::sqrt (M_Lz);
    for (UInt n = 0; n < M_quadruleY->nbQuadPt(); ++n)
    {
        for (UInt m = 0; m < M_quadruleZ->nbQuadPt(); ++m)
        {
            coeff += mu (0, x, M_quadruleY->quadPointCoor (n, 0) * M_Ly,
                         M_quadruleZ->quadPointCoor (m, 0) * M_Lz, 0) *
                     M_Ly * M_Lz * M_quadruleY->weight (n) * M_quadruleZ->weight (m) *
                     (
                         M_dphiy [p_j][n] * normy / M_Ly * M_dphiy [p_k][n] * normy / M_Ly *
                         M_phiz [q_j][m] * normz * M_phiz [q_k][m] * normz  +
                         M_phiy [p_j][n] * normy * M_phiy [p_k][n] * normy *
                         M_dphiz [q_j][m] * normz / M_Lz * M_dphiz [q_k][m] * normz / M_Lz
                     );
            coeff += beta (0, x, M_quadruleY->quadPointCoor (n, 0) * M_Ly,
                           M_quadruleZ->quadPointCoor (m, 0) * M_Lz, 1) *
                     M_Lz *
                     M_dphiy [p_j][n] * normy * M_phiy [p_k][n] * normy * M_quadruleY->weight (n) *
                     M_phiz [q_j][m] * normz * M_phiz [q_k][m] * normz * M_quadruleZ->weight (m);

            coeff += beta (0, x, M_quadruleY->quadPointCoor (n, 0) * M_Ly,
                           M_quadruleZ->quadPointCoor (m, 0) * M_Lz, 2) *
                     M_Ly *
                     M_phiy [p_j][n] * normy * M_phiy [p_k][n] * normy * M_quadruleY->weight (n) *
                     M_dphiz [q_j][m] * normz * M_phiz [q_k][m] * normz * M_quadruleZ->weight (m);

            coeff += sigma (0, x, M_quadruleY->quadPointCoor (n, 0) * M_Ly,
                            M_quadruleZ->quadPointCoor (m, 0) * M_Lz, 0) *
                     M_Ly * M_Lz *
                     M_phiy [p_j][n] * normy * M_phiy [p_k][n] * normy * M_quadruleY->weight (n) *
                     M_phiz [q_j][m] * normz * M_phiz [q_k][m] * normz * M_quadruleZ->weight (m);
        }
    }

    return coeff;
}
}   //end LifeV namespace
