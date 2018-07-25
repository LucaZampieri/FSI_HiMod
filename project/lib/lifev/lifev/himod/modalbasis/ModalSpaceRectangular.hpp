//@HEADER
/*
*******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University

    This file is part of LifeV.

    LifeV is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    LifeV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with LifeV.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************
*/
//@HEADER

/*!
    @file ModalSpaceAbstract.hpp
    @brief This file contains the description of the class that
    handle all the utilities related to the modal expansion
    on the transversal fiber. At first this file was named ModalSpace.hpp
    in the HiMod project, then it has been renamed due to a 
    generalization of the code to consider different geometries
    for the transversal section of the cylinder.

    @date 06/2013
    @author M. Aletti <teo.aletti@gmail.com>
    @author A. Bortolossi <andrea.bortolossi@gmail.com>

    @date 11/2013
    @contributor S. Guzzetti <sofia.guzzetti@gmail.com>
    @contributor M. Lupo Pasini <massimiliano.lupo.pasini@gmail.com>

 */
#ifndef __MODALSPACERECTANGULAR_HPP__
#define __MODALSPACERECTANGULAR_HPP__

#include <lifev/himod/basis/Basis1DAbstract.hpp>
#include <lifev/himod/modalbasis/ModalSpaceAbstract.hpp>

namespace LifeV
{
//! Eigenmap - simple struct for storing eigenvalues
/*!
 *  @author M. Aletti <teo.aletti@gmail.com>
 *  @author A. Bortolossi <andrea.bortolossi89@gmail.com>
 *
 *   \f$ \lambda=K_y+K_y \f$
 */

struct EigenMap
{
    //! \f$ \sqrt(K_y) \f$
    Real wp;
    //! \f$ \sqrt(K_z) \f$
    Real wq;

    UInt p;
    UInt q;

    static EigenMap make_eigenmap (const Real& _wp, const Real& _wq, const UInt& _p, const UInt& _q)
    {
        EigenMap a;

        a.wp = _wp;
        a.wq = _wq;
        a.p = _p;
        a.q = _q;

        return a;
    }
};

//! Comparison - a struct to compare sub-frequency
/*!
    @author A. Bortolossi
    @author M. Aletti <teo.aletti@gmail.com>

    This structure defines the order relationship between two pairs that represents
    the subfrequency of the problem
*/
struct Comparison
{
    //! @name Operators
    //@{
    /*! Comparison operator
        return a boolean value, true if the eigenvalues associated to a is lower then the eigenvlalue associated to b,
        if both a and b are associated to the same eigenvalue it returns true if a.first<b.first
        consider that this comparison generate an order relationship where, for example,
        (1,2)<(2,1) even if they generate the same eigenvalue.
        In this way the set used in the code can consider both (1,2) and (2,1) if we used only the comparison
        trough \f$\lambda(a)<\lambda(b)\f$ then it would have considered them equally and would have discarded one of them.
    */
    bool operator() ( EigenMap const& a, EigenMap const& b) const;
    //@}
};

//! ModalSpace - Class for Hi-Mod modal basis
/*!
 *  @author M. Aletti <teo.aletti@gmail.com>
 *  @author A. Bortolossi
 *  @contributor S. Guzzetti <sofia.guzzetti@gmail.com>
 *  @contributor M. Lupo Pasini <massimiliano.lupo.pasini@gmail.com>
 *  Class representing the modalbasis in a square.
 *
 */
class ModalSpaceRectangular: public ModalSpaceAbstract
{

public:

    //! @name Public Types
    //@{
    //! typedef for the utility map for the eigenvalues
    typedef std::vector<EigenMap>                        	EigenContainer;

    //! typedef for the educatedbasis
    typedef Basis1DAbstract                         		basis1d_type;

    //! typedef for the pointer to the educatedbasis
    typedef basis1d_type*                                   basis1d_ptrType;

    //@}

    //! @name Constructor & Destructor
    //@{
    /*!
      The only constructor available
      @param Ly length of the domain in y-direction
      @param Lz length of the domain in z-direction
      @param M dimension of the modal basis
      @param quady,quadz the quadrature rules
    */
    ModalSpaceRectangular
    (const Real& Ly, const Real& Lz, const UInt& M,  const QuadratureRule* quady = &quadRuleLobSeg32pt, const QuadratureRule* quadz = &quadRuleLobSeg32pt) :
    ModalSpaceAbstract(M), M_quadruleY (quady), M_quadruleZ (quadz), M_Ly (Ly), M_Lz (Lz), M_subMMax (0) {};

    //ModalSpace() : M_quadrule (&quadRuleSeg32pt) {};

    //! Destroys the monodimensional basis generators
    virtual ~ModalSpaceRectangular()
    {
        delete M_genbasisY;
        delete M_genbasisZ;
    };
    //@}

    //! @name Methods
    //@{

    //! Add the BC on the sides refers to the Y-axe. Set the basis generator on the Y direction.
    void addSliceBCY (const std::string& left, const std::string& right, const Real& mu = 1, const Real& Chi = 1);

    //! Add the BC on the sides refers to the Z-axe. Set the basis generator on the Z direction.
    void addSliceBCZ (const std::string& down,  const std::string& up, const Real& mu = 1, const Real& Chi = 1);

    //! Calls EvaluateBasis() method own by the basis generators.
    void evaluateBasis();

    //! Compute The Fourier Coefficients of the function \f$g=g(y,z)\f$
    virtual    std::vector<Real> fourierCoefficients (const function_Type& g) const;

    //! Compute The k-th Fourier coefficients of a function \f$f=f(x,y,z)\f$ in the point x
    virtual    Real fourierCoeffPointWise (const Real& x, const function_Type& f, const UInt& k) const;

    //! This methods compute the \f$\int_{[0,Ly]x[0,Lz]} \psi_j \psi_k\f$
    virtual    Real compute_PhiPhi (const UInt& j, const UInt& k) const;

    //! This methods compute the \f$\int_{[0,Ly]x[0,Lz]} Dy(\psi_j) \psi_k \f$
    virtual    Real compute_DyPhiPhi (const UInt& j, const UInt& k) const;

    //! This methods compute the \f$\int_{[0,Ly]x[0,Lz]} Dy(\psi_j) Dy(\psi_k)\f$
    virtual    Real compute_DyPhiDyPhi (const UInt& j, const UInt& k) const;

    //! This methods compute the \f$\int_{[0,Ly]x[0,Lz]} Dz(\psi_j) \psi_k\f$
    virtual    Real compute_DzPhiPhi (const UInt& j, const UInt& k) const;

    //! This methods compute the \f$\int_{[0,Ly]x[0,Lz]} Dz(\psi_j) Dz(\psi_k)\f$
    virtual    Real compute_DzPhiDzPhi (const UInt& j, const UInt& k) const;

    //! This methods compute the \f$\int_{[0,Ly]x[0,Lz]} \psi_k\f$
    virtual    Real compute_Phi (const UInt& k) const;

    //! This methods compute the \f$\int_{[0,Ly]} \eta_j \eta_k\f$
    virtual Real compute_y_PhiPhi (const UInt& j, const UInt& k) const;

    //! This methods compute the \f$\int_{[0,Lz]} \xi_j \xi_k\f$
    virtual Real compute_z_PhiPhi (const UInt& j, const UInt& k) const;


    virtual Real compute_R11 (const UInt& j, const UInt& k,
                              const function_Type& mu, const Real& x) const;

    virtual Real compute_R10 (const UInt& j, const UInt& k,
                              const function_Type& beta, const Real& x) const;

    virtual Real compute_R00 (const UInt& j, const UInt& k,
                              const function_Type& mu, const function_Type& beta, const function_Type& sigma, const Real& x) const;


    //! ShowMe method
    virtual    void showMe() const;

    //@}

    //! @name Get Methods
    //@{
    //! return the length on the y-direction
    Real Ly() const
    {
        return M_Ly;
    }

    //! return the length on the z-direction
    Real Lz() const
    {
        return M_Lz;
    }

    //! return the dimension of the modal basis
    UInt mtot() const
    {
        return M_mtot;
    }

    //! return max(dimension of sub-modal basis in y-direction, dimension of sub-modal basis in z-direction)
    UInt subMMax()
    {
        return M_subMMax;
    }

    //! return the value of the sub-base in the node n at the mode j
    virtual Real phiy (const UInt& j, const UInt& n) const
    {
        return M_phiy[j][n];
    }

    //! return the value of the sub-base derivate in the node n at the mode j
    virtual Real dphiy (const UInt& j, const UInt& n) const
    {
        return M_dphiy[j][n];
    }

    //! return the value of the sub-base in the node n at the mode j
    virtual Real phiz (const UInt& j, const UInt& n) const
    {
        return M_phiz[j][n];
    }

    //! return the value of the sub-base derivate in the node n at the mode j
    virtual Real dphiz (const UInt& j, const UInt& n) const
    {
        return M_dphiz[j][n];
    }

    //! return the Eigen triple of the mode j
    EigenMap eigenvalues (UInt j) const
    {
        return M_eigenvalues[j];
    }

    //! return the Quadrature rule
    const QuadratureRule&  qrY() const
    {
        return *M_quadruleY;
    }

    //! return the Quadrature rule
    const QuadratureRule&  qrZ() const
    {
        return *M_quadruleZ;
    }

    //! return the EBY
    const basis1d_ptrType  gbY() const
    {
        return M_genbasisY;
    }

    //! return the EBZ
    const basis1d_ptrType  gbZ() const
    {
        return M_genbasisZ;
    }
    //@}

protected:

    //! @name Private methods
    //@{
    //! This compute the max(dimension of sub-modal basis in y-direction, dimension of sub-modal basis in z-direction)
    void findSubMMax();

    //! This methods compute the list of the eigenvalues along with the subindeces (p,q), it fills the M_eigenvalues
    virtual void eigensProvider();

    //@}

    //! This variable it's usefull to link every frequency (j) with the sub-frequency (p,q)
    EigenContainer M_eigenvalues;

    //! Contains only the subfrequency along y direction
    MBVector_type M_eigenvaluesY;

    //! Contains only the subfrequency along y direction
    MBVector_type M_eigenvaluesZ;

    //! Pointer to the quadrature rule (Is it ok that it is a real pointer
    const QuadratureRule* M_quadruleY;

    //! Pointer to the quadrature rule (Is it ok that it is a real pointer
    const QuadratureRule* M_quadruleZ;

    //! Basis generator along y direction
    basis1d_ptrType M_genbasisY;

    //! Basis generator along z direction
    basis1d_ptrType M_genbasisZ;

    //! length of the domain in the y-direction
    Real M_Ly;

    //! length of the domain in the z-direction
    Real M_Lz;

    //! max(dimension of sub-modal basis in y-direction, dimension of sub-modal basis in z-direction)
    UInt M_subMMax;

    //! Evaluations of  \f$\tilde{\eta}(\tilde{y})_j)\f$
    MBMatrix_type M_phiy;

    //! Evaluations of  \f$\tilde{\xi}(\tilde{z})_j)\f$
    MBMatrix_type M_phiz;

    //! Evaluations of  \f$\partial_y\tilde{\eta}(\tilde{y})_j)\f$
    MBMatrix_type M_dphiy;

    //! Evaluations of  \f$\partial_z\tilde{\xi}(\tilde{z})_j)\f$
    MBMatrix_type M_dphiz;

};
} // namespace LifeV
#endif // __MODALSPACERECTANGULAR_HPP_
