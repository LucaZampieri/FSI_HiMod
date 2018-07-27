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
    @file
    @brief File containing the definition of the FE, quadrature rules, ...

    @contributor Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @mantainer Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */


#include <lifev/core/fem/QuadratureRule.hpp>
#include <lifev/core/fem/ReferenceFEScalar.hpp>
#include <lifev/core/fem/ReferenceFEHdiv.hpp>
#include <lifev/core/fem/ReferenceFEHybrid.hpp>
#include <lifev/core/fem/GeometricMap.hpp>

namespace LifeV
{

/*======================================================================
 *
 *                          Dummy Quadrature Rules
 *
 *=======================================================================*/

//! id of the quadrature rules on nodes
const int QUAD_RULE_DUMMY = 1;


static const QuadraturePoint pt_node_0pt[ 1 ];

const QuadratureRule quadRuleDummy ( pt_node_0pt,
                                     QUAD_RULE_DUMMY,
                                     "Dummy quadrature rule", NONE, 0, 0 );

/*======================================================================
 *
 *                          Quadrature Rules on Nodes
 *
 *=======================================================================*/
//! total number of quadrature rules on segments

const size_t NB_QUAD_RULE_NODE = 3;
//! id of the quadrature rules on nodes
const int QUAD_RULE_NODE_1PT = 1;


static const QuadraturePoint pt_node_1pt[ 1 ] =
{
    QuadraturePoint ( 0., 1. )
};
const QuadratureRule quadRuleNode1pt ( pt_node_1pt,
                                       QUAD_RULE_NODE_1PT,
                                       "Gauss Legendre 1 point on a node", POINT, 1, 1 );


/*======================================================================
 *
 *                          Quadrature Rules on segments
 *
 *=======================================================================*/
//! total number of quadrature rules on segments
const size_t NB_QUAD_RULE_SEG = 7;
//! id of the quadrature rules on segments
const size_t QUAD_RULE_SEG_1PT = 1;
const size_t QUAD_RULE_SEG_2PT = 2;
const size_t QUAD_RULE_SEG_3PT = 3;
const size_t QUAD_RULE_SEG_4PT = 4;
const size_t QUAD_RULE_SEG_32PT = 32;
const size_t QUAD_RULE_LOB_SEG_32PT = 32;
const size_t QUAD_RULE_LOB_SEG_64PT = 64;
const size_t QUAD_RULE_FQ_SEG_32PT = 32;
const size_t QUAD_RULE_CHEB_SEG_32PT = 32;
const size_t QUAD_RULE_QUADLEG_SEG_32PT = 32;

//----------------------------------------------------------------------

static const QuadraturePoint pt_seg_1pt[ 1 ] =
{
    QuadraturePoint ( 0.5, 1. )
};
const QuadratureRule quadRuleSeg1pt ( pt_seg_1pt,
                                      QUAD_RULE_SEG_1PT,
                                      "Gauss Legendre 1 point on a segment", LINE, 1, 1 );
//
//----------------------------------------------------------------------
const Real q2ptx1 = ( 1 - std::sqrt ( 1. / 3. ) ) / 2., q2ptx2 = ( 1 + std::sqrt ( 1. / 3. ) ) / 2.;
const Real q2ptw1 = 0.5, q2ptw2 = 0.5;

static const QuadraturePoint pt_seg_2pt[ 2 ] =
{
    QuadraturePoint ( q2ptx1 , q2ptw1 ),
    QuadraturePoint ( q2ptx2 , q2ptw2 )
};
const QuadratureRule quadRuleSeg2pt ( pt_seg_2pt,
                                      QUAD_RULE_SEG_2PT,
                                      "Gauss Legendre 2 points on a segment", LINE, 2, 3 );
//----------------------------------------------------------------------
const Real q3ptx1 = 0.5, q3ptx2 = ( 1 - std::sqrt ( 3. / 5. ) ) / 2., q3ptx3 = ( 1 + std::sqrt ( 3. / 5. ) ) / 2.;
const Real q3ptw1 = 8. / 18., q3ptw2 = 5. / 18., q3ptw3 = 5. / 18.;

static const QuadraturePoint pt_seg_3pt[ 3 ] =
{
    QuadraturePoint ( q3ptx1, q3ptw1 ),
    QuadraturePoint ( q3ptx2, q3ptw2 ),
    QuadraturePoint ( q3ptx3, q3ptw3 )
};

const QuadratureRule quadRuleSeg3pt ( pt_seg_3pt,
                                      QUAD_RULE_SEG_3PT,
                                      "Gauss Legendre 3 points on a segment", LINE, 3, 5 );
//----------------------------------------------------------------------
const Real q4ptx1 = (1. - sqrt ( (3. - 2.*sqrt (6. / 5.) ) / 7.) ) / 2., q4ptw1 = 0.5 * (18. + sqrt (30) ) / 36.;
const Real q4ptx2 = (1. + sqrt ( (3. - 2.*sqrt (6. / 5.) ) / 7.) ) / 2., q4ptw2 = 0.5 * (18. + sqrt (30) ) / 36.;
const Real q4ptx3 = (1. - sqrt ( (3. + 2.*sqrt (6. / 5.) ) / 7.) ) / 2., q4ptw3 = 0.5 * (18. - sqrt (30) ) / 36.;
const Real q4ptx4 = (1. + sqrt ( (3. + 2.*sqrt (6. / 5.) ) / 7.) ) / 2., q4ptw4 = 0.5 * (18. - sqrt (30) ) / 36.;

static const QuadraturePoint pt_seg_4pt[ 4 ] =
{
    QuadraturePoint ( q4ptx1, q4ptw1 ),
    QuadraturePoint ( q4ptx2, q4ptw2 ),
    QuadraturePoint ( q4ptx3, q4ptw3 ),
    QuadraturePoint ( q4ptx4, q4ptw4 )
};

const QuadratureRule quadRuleSeg4pt ( pt_seg_4pt,
                                      QUAD_RULE_SEG_4PT,
                                      "Gauss Legendre 4 points on a segment", LINE, 4, 7 );
//-------------------------------------------------------------------------
const Real ql32ptx1 = 0.0000000000000000000000e+00, ql32ptw1 = 1.0080645161290322300851e-03;
const Real ql32ptx2 = 3.6955330136193009771262e-03, ql32ptw2 = 6.1990532506869221174295e-03;
const Real ql32ptx3 = 1.2352654758645387200744e-02, ql32ptw3 = 1.1099776444645981926529e-02;
const Real ql32ptx4 = 2.5857580791383838469955e-02, ql32ptw4 = 1.5887567705457732170071e-02;
const Real ql32ptx5 = 4.4075030468134046568451e-02, ql32ptw5 = 2.0517100793031362720997e-02;
const Real ql32ptx6 = 6.6823761993662245117775e-02, ql32ptw6 = 2.4942635668110602759073e-02;
const Real ql32ptx7 = 9.3877634111278807083067e-02, ql32ptw7 = 2.9120248624027934847280e-02;
const Real ql32ptx8 = 1.2496775303166257620191e-01, ql32ptw8 = 3.3008438628577274620568e-02;
const Real ql32ptx9 = 1.5978512219222457124701e-01, ql32ptw9 = 3.6568569801339517733041e-02;
const Real ql32ptx10 = 1.9798370642578944078593e-01, ql32ptw10 = 3.9765262846053126000179e-02;
const Real ql32ptx11 = 2.3918386855921736078301e-01, ql32ptw11 = 4.2566748974834117680288e-02;
const Real ql32ptx12 = 2.8297614139907656394257e-01, ql32ptw12 = 4.4945186478678914032692e-02;
const Real ql32ptx13 = 3.2892529673055925787395e-01, ql32ptw13 = 4.6876937773406908427365e-02;
const Real ql32ptx14 = 3.7657467057489735218212e-01, ql32ptw14 = 4.8342804474001299652741e-02;
const Real ql32ptx15 = 4.2545070159317627256357e-01, ql32ptw15 = 4.9328218270380891352644e-02;
const Real ql32ptx16 = 4.7506763747670338604578e-01, ql32ptw16 = 4.9823385750638389535272e-02;
const Real ql32ptx17 = 5.2493236252329666946537e-01, ql32ptw17 = 4.9823385750638389535272e-02;
const Real ql32ptx18 = 5.7454929840682378294758e-01, ql32ptw18 = 4.9328218270380891352644e-02;
const Real ql32ptx19 = 6.2342532942510264781788e-01, ql32ptw19 = 4.8342804474001299652741e-02;
const Real ql32ptx20 = 6.7107470326944074212605e-01, ql32ptw20 = 4.6876937773406908427365e-02;
const Real ql32ptx21 = 7.1702385860092343605743e-01, ql32ptw21 = 4.4945186478678914032692e-02;
const Real ql32ptx22 = 7.6081613144078263921699e-01, ql32ptw22 = 4.2566748974834117680288e-02;
const Real ql32ptx23 = 8.0201629357421055921407e-01, ql32ptw23 = 3.9765262846053126000179e-02;
const Real ql32ptx24 = 8.4021487780777537324184e-01, ql32ptw24 = 3.6568569801339517733041e-02;
const Real ql32ptx25 = 8.7503224696833736828694e-01, ql32ptw25 = 3.3008438628577274620568e-02;
const Real ql32ptx26 = 9.0612236588872119291693e-01, ql32ptw26 = 2.9120248624027934847280e-02;
const Real ql32ptx27 = 9.3317623800633775488222e-01, ql32ptw27 = 2.4942635668110602759073e-02;
const Real ql32ptx28 = 9.5592496953186589792040e-01, ql32ptw28 = 2.0517100793031362720997e-02;
const Real ql32ptx29 = 9.7414241920861610601889e-01, ql32ptw29 = 1.5887567705457732170071e-02;
const Real ql32ptx30 = 9.8764734524135455728810e-01, ql32ptw30 = 1.1099776444645981926529e-02;
const Real ql32ptx31 = 9.9630446698638075453403e-01, ql32ptw31 = 6.1990532506869221174295e-03;
const Real ql32ptx32 = 1.0000000000000000000000e+00, ql32ptw32 = 1.0080645161290322300851e-03;

static const QuadraturePoint pt_lob_seg_32pt[32] =
{
    QuadraturePoint ( ql32ptx1, ql32ptw1 ),
    QuadraturePoint ( ql32ptx2, ql32ptw2 ),
    QuadraturePoint ( ql32ptx3, ql32ptw3 ),
    QuadraturePoint ( ql32ptx4, ql32ptw4 ),
    QuadraturePoint ( ql32ptx5, ql32ptw5 ),
    QuadraturePoint ( ql32ptx6, ql32ptw6 ),
    QuadraturePoint ( ql32ptx7, ql32ptw7 ),
    QuadraturePoint ( ql32ptx8, ql32ptw8 ),
    QuadraturePoint ( ql32ptx9, ql32ptw9 ),
    QuadraturePoint ( ql32ptx10, ql32ptw10 ),
    QuadraturePoint ( ql32ptx11, ql32ptw11 ),
    QuadraturePoint ( ql32ptx12, ql32ptw12 ),
    QuadraturePoint ( ql32ptx13, ql32ptw13 ),
    QuadraturePoint ( ql32ptx14, ql32ptw14 ),
    QuadraturePoint ( ql32ptx15, ql32ptw15 ),
    QuadraturePoint ( ql32ptx16, ql32ptw16 ),
    QuadraturePoint ( ql32ptx17, ql32ptw17 ),
    QuadraturePoint ( ql32ptx18, ql32ptw18 ),
    QuadraturePoint ( ql32ptx19, ql32ptw19 ),
    QuadraturePoint ( ql32ptx20, ql32ptw20 ),
    QuadraturePoint ( ql32ptx21, ql32ptw21 ),
    QuadraturePoint ( ql32ptx22, ql32ptw22 ),
    QuadraturePoint ( ql32ptx23, ql32ptw23 ),
    QuadraturePoint ( ql32ptx24, ql32ptw24 ),
    QuadraturePoint ( ql32ptx25, ql32ptw25 ),
    QuadraturePoint ( ql32ptx26, ql32ptw26 ),
    QuadraturePoint ( ql32ptx27, ql32ptw27 ),
    QuadraturePoint ( ql32ptx28, ql32ptw28 ),
    QuadraturePoint ( ql32ptx29, ql32ptw29 ),
    QuadraturePoint ( ql32ptx30, ql32ptw30 ),
    QuadraturePoint ( ql32ptx31, ql32ptw31 ),
    QuadraturePoint ( ql32ptx32, ql32ptw32 )
};

const QuadratureRule quadRuleLobSeg32pt ( pt_lob_seg_32pt,
                                      QUAD_RULE_LOB_SEG_32PT,
                                      "Gauss Lobatto 32 points on a segment", LINE, 32, 61 );

//-------------------------------------------------------------------------
// Legendre quadrature 32 points with quadratic shift on a segment for radius
//-------------------------------------------------------------------------
const Real q32ptx1QuadLeg  = 3.698741779658611e-02, q32ptw1QuadLeg  = 4.743917274834677e-02;
const Real q32ptx2QuadLeg  = 8.481889074590523e-02, q32ptw2QuadLeg  = 4.796807228845817e-02;
const Real q32ptx3QuadLeg  = 1.327360998607643e-01, q32ptw3QuadLeg  = 4.782433967831185e-02;
const Real q32ptx4QuadLeg  = 1.804077659945108e-01, q32ptw4QuadLeg  = 4.749499380484579e-02;
const Real q32ptx5QuadLeg  = 2.276827224823481e-01, q32ptw5QuadLeg  = 4.703463832828563e-02;
const Real q32ptx6QuadLeg  = 2.744379586240122e-01, q32ptw6QuadLeg  = 4.645681989298442e-02;
const Real q32ptx7QuadLeg  = 3.205590460680042e-01, q32ptw7QuadLeg  = 4.576699222682844e-02;
const Real q32ptx8QuadLeg  = 3.659357055957442e-01, q32ptw8QuadLeg  = 4.496843418791501e-02;
const Real q32ptx9QuadLeg  = 4.104605541765157e-01, q32ptw9QuadLeg  = 4.406379210664461e-02;
const Real q32ptx10QuadLeg = 4.540287671278317e-01, q32ptw10QuadLeg = 4.305558449617734e-02;
const Real q32ptx11QuadLeg = 4.965380605088449e-01, q32ptw11QuadLeg = 4.194639386836224e-02;
const Real q32ptx12QuadLeg = 5.378888007150571e-01, q32ptw12QuadLeg = 4.073894682687255e-02;
const Real q32ptx13QuadLeg = 5.779841683462759e-01, q32ptw13QuadLeg = 3.943614881209893e-02;
const Real q32ptx14QuadLeg = 6.167303453487038e-01, q32ptw14QuadLeg = 3.804109842679473e-02;
const Real q32ptx15QuadLeg = 6.540367109028374e-01, q32ptw15QuadLeg = 3.655709170638695e-02;
const Real q32ptx16QuadLeg = 6.898160386335844e-01, q32ptw16QuadLeg = 3.498762101340755e-02;
const Real q32ptx17QuadLeg = 7.239846910286634e-01, q32ptw17QuadLeg = 3.333637082075596e-02;
const Real q32ptx18QuadLeg = 7.564628085976192e-01, q32ptw18QuadLeg = 3.160721154836953e-02;
const Real q32ptx19QuadLeg = 7.871744921718872e-01, q32ptw19QuadLeg = 2.980419208639472e-02;
const Real q32ptx20QuadLeg = 8.160479772299322e-01, q32ptw20QuadLeg = 2.793153136818403e-02;
const Real q32ptx21QuadLeg = 8.430157994161899e-01, q32ptw21QuadLeg = 2.599360921382053e-02;
const Real q32ptx22QuadLeg = 8.680149506005727e-01, q32ptw22QuadLeg = 2.399495658724078e-02;
const Real q32ptx23QuadLeg = 8.909870249450219e-01, q32ptw23QuadLeg = 2.194024536774126e-02;
const Real q32ptx24QuadLeg = 9.118783545326140e-01, q32ptw24QuadLeg = 1.983427771622279e-02;
const Real q32ptx25QuadLeg = 9.306401341926668e-01, q32ptw25QuadLeg = 1.768197511529600e-02;
const Real q32ptx26QuadLeg = 9.472285352458356e-01, q32ptw26QuadLeg = 1.548836719306212e-02;
const Real q32ptx27QuadLeg = 9.616048080507319e-01, q32ptw27QuadLeg = 1.325858056121680e-02;
const Real q32ptx28QuadLeg = 9.737353736426678e-01, q32ptw28QuadLeg = 1.099782835812489e-02;
const Real q32ptx29QuadLeg = 9.835919062135829e-01, q32ptw29QuadLeg = 8.711403249789198e-03;
const Real q32ptx30QuadLeg = 9.911514151701309e-01, q32ptw30QuadLeg = 6.404688759109404e-03;
const Real q32ptx31QuadLeg = 9.963963848652976e-01, q32ptw31QuadLeg = 4.083313372595637e-03;
const Real q32ptx32QuadLeg = 9.993157313505781e-01, q32ptw32QuadLeg = 1.755853978197768e-03;

static const QuadraturePoint pt_QLeg_seg_32pt[32] =
{
    QuadraturePoint ( q32ptx1QuadLeg, q32ptw1QuadLeg ),
    QuadraturePoint ( q32ptx2QuadLeg, q32ptw2QuadLeg ),
    QuadraturePoint ( q32ptx3QuadLeg, q32ptw3QuadLeg ),
    QuadraturePoint ( q32ptx4QuadLeg, q32ptw4QuadLeg ),
    QuadraturePoint ( q32ptx5QuadLeg, q32ptw5QuadLeg ),
    QuadraturePoint ( q32ptx6QuadLeg, q32ptw6QuadLeg ),
    QuadraturePoint ( q32ptx7QuadLeg, q32ptw7QuadLeg ),
    QuadraturePoint ( q32ptx8QuadLeg, q32ptw8QuadLeg ),
    QuadraturePoint ( q32ptx9QuadLeg, q32ptw9QuadLeg ),
    QuadraturePoint ( q32ptx10QuadLeg, q32ptw10QuadLeg ),
    QuadraturePoint ( q32ptx11QuadLeg, q32ptw11QuadLeg ),
    QuadraturePoint ( q32ptx12QuadLeg, q32ptw12QuadLeg ),
    QuadraturePoint ( q32ptx13QuadLeg, q32ptw13QuadLeg ),
    QuadraturePoint ( q32ptx14QuadLeg, q32ptw14QuadLeg ),
    QuadraturePoint ( q32ptx15QuadLeg, q32ptw15QuadLeg ),
    QuadraturePoint ( q32ptx16QuadLeg, q32ptw16QuadLeg ),
    QuadraturePoint ( q32ptx17QuadLeg, q32ptw17QuadLeg ),
    QuadraturePoint ( q32ptx18QuadLeg, q32ptw18QuadLeg ),
    QuadraturePoint ( q32ptx19QuadLeg, q32ptw19QuadLeg ),
    QuadraturePoint ( q32ptx20QuadLeg, q32ptw20QuadLeg ),
    QuadraturePoint ( q32ptx21QuadLeg, q32ptw21QuadLeg ),
    QuadraturePoint ( q32ptx22QuadLeg, q32ptw22QuadLeg ),
    QuadraturePoint ( q32ptx23QuadLeg, q32ptw23QuadLeg ),
    QuadraturePoint ( q32ptx24QuadLeg, q32ptw24QuadLeg ),
    QuadraturePoint ( q32ptx25QuadLeg, q32ptw25QuadLeg ),
    QuadraturePoint ( q32ptx26QuadLeg, q32ptw26QuadLeg ),
    QuadraturePoint ( q32ptx27QuadLeg, q32ptw27QuadLeg ),
    QuadraturePoint ( q32ptx28QuadLeg, q32ptw28QuadLeg ),
    QuadraturePoint ( q32ptx29QuadLeg, q32ptw29QuadLeg ),
    QuadraturePoint ( q32ptx30QuadLeg, q32ptw30QuadLeg ),
    QuadraturePoint ( q32ptx31QuadLeg, q32ptw31QuadLeg ),
    QuadraturePoint ( q32ptx32QuadLeg, q32ptw32QuadLeg )
};

const QuadratureRule quadRuleQuadLegSeg32pt ( pt_QLeg_seg_32pt,
                                      QUAD_RULE_QUADLEG_SEG_32PT,
                                      "Legendre quadrature 32 points with quadratic shift on a segment", LINE, 32, 31 );

//-------------------------------------------------------------------------
// Chebyshev quadrature 32 points on a segment for radius
//-------------------------------------------------------------------------
const Real q32ptx1ch   = 2.454122852291216e-02, q32ptw1ch   = 4.595064493652836e-02;
const Real q32ptx2ch   = 7.356456359966729e-02, q32ptw2ch   = 4.953708679212175e-02;
const Real q32ptx3ch   = 1.224106751992162e-01, q32ptw3ch   = 4.848919620679813e-02;
const Real q32ptx4ch   = 1.709618887603014e-01, q32ptw4ch   = 4.848325136654019e-02;
const Real q32ptx5ch   = 2.191012401568698e-01, q32ptw5ch   = 4.782399359270016e-02;
const Real q32ptx6ch   = 2.667127574748984e-01, q32ptw6ch   = 4.735496641915435e-02;
const Real q32ptx7ch   = 3.136817403988915e-01, q32ptw7ch   = 4.657880656205771e-02;
const Real q32ptx8ch   = 3.598950365349883e-01, q32ptw8ch   = 4.581982545153646e-02;
const Real q32ptx9ch   = 4.052413140049899e-01, q32ptw9ch   = 4.486095221785756e-02;
const Real q32ptx10ch  = 4.496113296546066e-01, q32ptw10ch  = 4.385661345131218e-02;
const Real q32ptx11ch  = 4.928981922297842e-01, q32ptw11ch  = 4.270323874878678e-02;
const Real q32ptx12ch  = 5.349976198870973e-01, q32ptw12ch  = 4.147594478363177e-02;
const Real q32ptx13ch  = 5.758081914178453e-01, q32ptw13ch  = 4.013106751762176e-02;
const Real q32ptx14ch  = 6.152315905806269e-01, q32ptw14ch  = 3.869798415925106e-02;
const Real q32ptx15ch  = 6.531728429537768e-01, q32ptw15ch  = 3.717096130018498e-02;
const Real q32ptx16ch  = 6.895405447370669e-01, q32ptw16ch  = 3.554831377465587e-02;
const Real q32ptx17ch  = 7.242470829514670e-01, q32ptw17ch  = 3.385225277977551e-02;
const Real q32ptx18ch  = 7.572088465064846e-01, q32ptw18ch  = 3.205665144181624e-02;
const Real q32ptx19ch  = 7.883464276266063e-01, q32ptw19ch  = 3.020739221298389e-02;
const Real q32ptx20ch  = 8.175848131515837e-01, q32ptw20ch  = 2.825620400084681e-02;
const Real q32ptx21ch  = 8.448535652497071e-01, q32ptw21ch  = 2.627187444562263e-02;
const Real q32ptx22ch  = 8.700869911087114e-01, q32ptw22ch  = 2.418316819833016e-02;
const Real q32ptx23ch  = 8.932243011955153e-01, q32ptw23ch  = 2.208405563832186e-02;
const Real q32ptx24ch  = 9.142097557035307e-01, q32ptw24ch  = 1.987621092739672e-02;
const Real q32ptx25ch  = 9.329927988347390e-01, q32ptw25ch  = 1.768500734083568e-02;
const Real q32ptx26ch  = 9.495281805930367e-01, q32ptw26ch  = 1.537575790061648e-02;
const Real q32ptx27ch  = 9.637760657954398e-01, q32ptw27ch  = 1.311869657068933e-02;
const Real q32ptx28ch  = 9.757021300385286e-01, q32ptw28ch  = 1.072252810304426e-02;
const Real q32ptx29ch  = 9.852776423889412e-01, q32ptw29ch  = 8.433784617185145e-03;
const Real q32ptx30ch  = 9.924795345987100e-01, q32ptw30ch  = 5.951968489589378e-03;
const Real q32ptx31ch  = 9.972904566786902e-01, q32ptw31ch  = 3.697241011080058e-03;
const Real q32ptx32ch  = 9.996988186962043e-01, q32ptw32ch  = 1.051333888160896e-03;

static const QuadraturePoint pt_cheb_seg_32pt[32] =
{
    QuadraturePoint ( q32ptx1ch, q32ptw1ch ),
    QuadraturePoint ( q32ptx2ch, q32ptw2ch ),
    QuadraturePoint ( q32ptx3ch, q32ptw3ch ),
    QuadraturePoint ( q32ptx4ch, q32ptw4ch ),
    QuadraturePoint ( q32ptx5ch, q32ptw5ch ),
    QuadraturePoint ( q32ptx6ch, q32ptw6ch ),
    QuadraturePoint ( q32ptx7ch, q32ptw7ch ),
    QuadraturePoint ( q32ptx8ch, q32ptw8ch ),
    QuadraturePoint ( q32ptx9ch, q32ptw9ch ),
    QuadraturePoint ( q32ptx10ch, q32ptw10ch ),
    QuadraturePoint ( q32ptx11ch, q32ptw11ch ),
    QuadraturePoint ( q32ptx12ch, q32ptw12ch ),
    QuadraturePoint ( q32ptx13ch, q32ptw13ch ),
    QuadraturePoint ( q32ptx14ch, q32ptw14ch ),
    QuadraturePoint ( q32ptx15ch, q32ptw15ch ),
    QuadraturePoint ( q32ptx16ch, q32ptw16ch ),
    QuadraturePoint ( q32ptx17ch, q32ptw17ch ),
    QuadraturePoint ( q32ptx18ch, q32ptw18ch ),
    QuadraturePoint ( q32ptx19ch, q32ptw19ch ),
    QuadraturePoint ( q32ptx20ch, q32ptw20ch ),
    QuadraturePoint ( q32ptx21ch, q32ptw21ch ),
    QuadraturePoint ( q32ptx22ch, q32ptw22ch ),
    QuadraturePoint ( q32ptx23ch, q32ptw23ch ),
    QuadraturePoint ( q32ptx24ch, q32ptw24ch ),
    QuadraturePoint ( q32ptx25ch, q32ptw25ch ),
    QuadraturePoint ( q32ptx26ch, q32ptw26ch ),
    QuadraturePoint ( q32ptx27ch, q32ptw27ch ),
    QuadraturePoint ( q32ptx28ch, q32ptw28ch ),
    QuadraturePoint ( q32ptx29ch, q32ptw29ch ),
    QuadraturePoint ( q32ptx30ch, q32ptw30ch ),
    QuadraturePoint ( q32ptx31ch, q32ptw31ch ),
    QuadraturePoint ( q32ptx32ch, q32ptw32ch )
//    QuadraturePoint ( q32ptx33nc, q32ptw33nc )
};

const QuadratureRule quadRuleChebSeg32pt ( pt_cheb_seg_32pt,
                                      QUAD_RULE_CHEB_SEG_32PT,
                                      "Even-Chebyshev quadrature 32 points on a segment", LINE, 32, 31 );
                                      
//-------------------------------------------------------------------------
// Fourier quadrature 32 points on a segment for theta
//-------------------------------------------------------------------------
const Real q32ptx1nc  = 0.000000000000000e+00, q32ptw1nc  = 3.125e-02;
const Real q32ptx2nc  = 3.125000000000000e-02, q32ptw2nc  = 3.125e-02;
const Real q32ptx3nc  = 6.250000000000000e-02, q32ptw3nc  = 3.125e-02;
const Real q32ptx4nc  = 9.375000000000000e-02, q32ptw4nc  = 3.125e-02;
const Real q32ptx5nc  = 1.250000000000000e-01, q32ptw5nc  = 3.125e-02;
const Real q32ptx6nc  = 1.562500000000000e-01, q32ptw6nc  = 3.125e-02;
const Real q32ptx7nc  = 1.875000000000000e-01, q32ptw7nc  = 3.125e-02;
const Real q32ptx8nc  = 2.187500000000000e-01, q32ptw8nc  = 3.125e-02;
const Real q32ptx9nc  = 2.500000000000000e-01, q32ptw9nc  = 3.125e-02;
const Real q32ptx10nc = 2.812500000000000e-01, q32ptw10nc = 3.125e-02;
const Real q32ptx11nc = 3.125000000000000e-01, q32ptw11nc = 3.125e-02;
const Real q32ptx12nc = 3.437500000000000e-01, q32ptw12nc = 3.125e-02;
const Real q32ptx13nc = 3.750000000000000e-01, q32ptw13nc = 3.125e-02;
const Real q32ptx14nc = 4.062500000000000e-01, q32ptw14nc = 3.125e-02;
const Real q32ptx15nc = 4.375000000000000e-01, q32ptw15nc = 3.125e-02;
const Real q32ptx16nc = 4.687500000000000e-01, q32ptw16nc = 3.125e-02;
const Real q32ptx17nc = 5.000000000000000e-01, q32ptw17nc = 3.125e-02;
const Real q32ptx18nc = 5.312500000000000e-01, q32ptw18nc = 3.125e-02;
const Real q32ptx19nc = 5.625000000000000e-01, q32ptw19nc = 3.125e-02;
const Real q32ptx20nc = 5.937500000000000e-01, q32ptw20nc = 3.125e-02;
const Real q32ptx21nc = 6.250000000000000e-01, q32ptw21nc = 3.125e-02;
const Real q32ptx22nc = 6.562500000000000e-01, q32ptw22nc = 3.125e-02;
const Real q32ptx23nc = 6.875000000000000e-01, q32ptw23nc = 3.125e-02;
const Real q32ptx24nc = 7.187500000000000e-01, q32ptw24nc = 3.125e-02;
const Real q32ptx25nc = 7.500000000000000e-01, q32ptw25nc = 3.125e-02;
const Real q32ptx26nc = 7.812500000000000e-01, q32ptw26nc = 3.125e-02;
const Real q32ptx27nc = 8.125000000000000e-01, q32ptw27nc = 3.125e-02;
const Real q32ptx28nc = 8.437500000000000e-01, q32ptw28nc = 3.125e-02;
const Real q32ptx29nc = 8.750000000000000e-01, q32ptw29nc = 3.125e-02;
const Real q32ptx30nc = 9.062500000000000e-01, q32ptw30nc = 3.125e-02;
const Real q32ptx31nc = 9.375000000000000e-01, q32ptw31nc = 3.125e-02;
const Real q32ptx32nc = 9.687500000000000e-01, q32ptw32nc = 3.125e-02;
//const Real q32ptx33nc = 1.000000000000000e+00, q32ptw33nc = 6.786386919729440e-03;

static const QuadraturePoint pt_fq_seg_32pt[32] =
{
    QuadraturePoint ( q32ptx1nc, q32ptw1nc ),
    QuadraturePoint ( q32ptx2nc, q32ptw2nc ),
    QuadraturePoint ( q32ptx3nc, q32ptw3nc ),
    QuadraturePoint ( q32ptx4nc, q32ptw4nc ),
    QuadraturePoint ( q32ptx5nc, q32ptw5nc ),
    QuadraturePoint ( q32ptx6nc, q32ptw6nc ),
    QuadraturePoint ( q32ptx7nc, q32ptw7nc ),
    QuadraturePoint ( q32ptx8nc, q32ptw8nc ),
    QuadraturePoint ( q32ptx9nc, q32ptw9nc ),
    QuadraturePoint ( q32ptx10nc, q32ptw10nc ),
    QuadraturePoint ( q32ptx11nc, q32ptw11nc ),
    QuadraturePoint ( q32ptx12nc, q32ptw12nc ),
    QuadraturePoint ( q32ptx13nc, q32ptw13nc ),
    QuadraturePoint ( q32ptx14nc, q32ptw14nc ),
    QuadraturePoint ( q32ptx15nc, q32ptw15nc ),
    QuadraturePoint ( q32ptx16nc, q32ptw16nc ),
    QuadraturePoint ( q32ptx17nc, q32ptw17nc ),
    QuadraturePoint ( q32ptx18nc, q32ptw18nc ),
    QuadraturePoint ( q32ptx19nc, q32ptw19nc ),
    QuadraturePoint ( q32ptx20nc, q32ptw20nc ),
    QuadraturePoint ( q32ptx21nc, q32ptw21nc ),
    QuadraturePoint ( q32ptx22nc, q32ptw22nc ),
    QuadraturePoint ( q32ptx23nc, q32ptw23nc ),
    QuadraturePoint ( q32ptx24nc, q32ptw24nc ),
    QuadraturePoint ( q32ptx25nc, q32ptw25nc ),
    QuadraturePoint ( q32ptx26nc, q32ptw26nc ),
    QuadraturePoint ( q32ptx27nc, q32ptw27nc ),
    QuadraturePoint ( q32ptx28nc, q32ptw28nc ),
    QuadraturePoint ( q32ptx29nc, q32ptw29nc ),
    QuadraturePoint ( q32ptx30nc, q32ptw30nc ),
    QuadraturePoint ( q32ptx31nc, q32ptw31nc ),
    QuadraturePoint ( q32ptx32nc, q32ptw32nc )
};

const QuadratureRule quadRuleFQSeg32pt ( pt_fq_seg_32pt,
                                      QUAD_RULE_FQ_SEG_32PT,
                                      "Fourier quadrature 32 points on a segment", LINE, 32, 31 );
                                      
//-------------------------------------------------------------------------
const Real q32ptx1 = 1.36806907525910403933e-03, q32ptw1 = 3.50930500473500828207e-03;
const Real q32ptx2 = 7.19424422736580915227e-03, q32ptw2 = 8.13719736545264089866e-03;
const Real q32ptx3 = 1.76188722062468050567e-02, q32ptw3 = 1.26960326546312583795e-02;
const Real q32ptx4 = 3.25469620311301111037e-02, q32ptw4 = 1.71369314565104902126e-02;
const Real q32ptx5 = 5.18394221169737878796e-02, q32ptw5 = 2.14179490111136676400e-02;
const Real q32ptx6 = 7.53161931337150702959e-02, q32ptw6 = 2.54990296311881220470e-02;
const Real q32ptx7 = 1.02758102016028862735e-01, q32ptw7 = 2.93420467392677165874e-02;
const Real q32ptx8 = 1.33908940629855199855e-01, q32ptw8 = 3.29111113881807859638e-02;
const Real q32ptx9 = 1.68477866534892439798e-01, q32ptw9 = 3.61728970544242731111e-02;
const Real q32ptx10 = 2.06142121379618625809e-01, q32ptw10 = 3.90969478935349543103e-02;
const Real q32ptx11 = 2.46550045533885375804e-01, q32ptw11 = 4.16559621134734853198e-02;
const Real q32ptx12 = 2.89324361934682361408e-01, q32ptw12 = 4.38260465022020442860e-02;
const Real q32ptx13 = 3.34065698858936110938e-01, q32ptw13 = 4.55869393478821952059e-02;
const Real q32ptx14 = 3.80356318873931620317e-01, q32ptw14 = 4.69221995404019223685e-02;
const Real q32ptx15 = 4.27764019208601742328e-01, q32ptw15 = 4.78193600396371598649e-02;
const Real q32ptx16 = 4.75846167156130650522e-01, q32ptw16 = 4.82700442573640725596e-02;
const Real q32ptx17 = 5.24153832843869071922e-01, q32ptw17 = 4.82700442573634619370e-02;
const Real q32ptx18 = 5.72235980791398368694e-01, q32ptw18 = 4.78193600396370835370e-02;
const Real q32ptx19 = 6.19643681126068490705e-01, q32ptw19 = 4.69221995404022415577e-02;
const Real q32ptx20 = 6.65934301141063778040e-01, q32ptw20 = 4.55869393478817511167e-02;
const Real q32ptx21 = 7.10675638065317860637e-01, q32ptw21 = 4.38260465022019332637e-02;
const Real q32ptx22 = 7.53449954466114735219e-01, q32ptw22 = 4.16559621134732771530e-02;
const Real q32ptx23 = 7.93857878620381374191e-01, q32ptw23 = 3.90969478935357730998e-02;
const Real q32ptx24 = 8.31522133465107504691e-01, q32ptw24 = 3.61728970544242176000e-02;
const Real q32ptx25 = 8.66091059370144966678e-01, q32ptw25 = 3.29111113881806402470e-02;
const Real q32ptx26 = 8.97241897983971248287e-01, q32ptw26 = 2.93420467392676680152e-02;
const Real q32ptx27 = 9.24683806866285151749e-01, q32ptw27 = 2.54990296311881463331e-02;
const Real q32ptx28 = 9.48160577883026212120e-01, q32ptw28 = 2.14179490111134143704e-02;
const Real q32ptx29 = 9.67453037968869833385e-01, q32ptw29 = 1.71369314565111181825e-02;
const Real q32ptx30 = 9.82381127793753194943e-01, q32ptw30 = 1.26960326546313798102e-02;
const Real q32ptx31 = 9.92805755772633968803e-01, q32ptw31 = 8.13719736545309366149e-03;
const Real q32ptx32 = 9.98631930924740673916e-01, q32ptw32 = 3.50930500473497055183e-03;

static const QuadraturePoint pt_seg_32pt[32] =
{
    QuadraturePoint ( q32ptx1, q32ptw1 ),
    QuadraturePoint ( q32ptx2, q32ptw2 ),
    QuadraturePoint ( q32ptx3, q32ptw3 ),
    QuadraturePoint ( q32ptx4, q32ptw4 ),
    QuadraturePoint ( q32ptx5, q32ptw5 ),
    QuadraturePoint ( q32ptx6, q32ptw6 ),
    QuadraturePoint ( q32ptx7, q32ptw7 ),
    QuadraturePoint ( q32ptx8, q32ptw8 ),
    QuadraturePoint ( q32ptx9, q32ptw9 ),
    QuadraturePoint ( q32ptx10, q32ptw10 ),
    QuadraturePoint ( q32ptx11, q32ptw11 ),
    QuadraturePoint ( q32ptx12, q32ptw12 ),
    QuadraturePoint ( q32ptx13, q32ptw13 ),
    QuadraturePoint ( q32ptx14, q32ptw14 ),
    QuadraturePoint ( q32ptx15, q32ptw15 ),
    QuadraturePoint ( q32ptx16, q32ptw16 ),
    QuadraturePoint ( q32ptx17, q32ptw17 ),
    QuadraturePoint ( q32ptx18, q32ptw18 ),
    QuadraturePoint ( q32ptx19, q32ptw19 ),
    QuadraturePoint ( q32ptx20, q32ptw20 ),
    QuadraturePoint ( q32ptx21, q32ptw21 ),
    QuadraturePoint ( q32ptx22, q32ptw22 ),
    QuadraturePoint ( q32ptx23, q32ptw23 ),
    QuadraturePoint ( q32ptx24, q32ptw24 ),
    QuadraturePoint ( q32ptx25, q32ptw25 ),
    QuadraturePoint ( q32ptx26, q32ptw26 ),
    QuadraturePoint ( q32ptx27, q32ptw27 ),
    QuadraturePoint ( q32ptx28, q32ptw28 ),
    QuadraturePoint ( q32ptx29, q32ptw29 ),
    QuadraturePoint ( q32ptx30, q32ptw30 ),
    QuadraturePoint ( q32ptx31, q32ptw31 ),
    QuadraturePoint ( q32ptx32, q32ptw32 )
};

const QuadratureRule quadRuleSeg32pt ( pt_seg_32pt,
                                      QUAD_RULE_SEG_32PT,
                                      "Gauss 32 points on a segment", LINE, 32, 63 );

//-----------------------------------------------------------------------------------------

const Real ql64ptx1 = 0.0000000000000000000000e+00, ql64ptw1 = 2.4801587301587300210537e-04;
const Real ql64ptx2 = 9.1006424891837411905726e-04, ql64ptw2 = 1.5280041224562451382396e-03;
const Real ql64ptx3 = 3.0486366484713811608742e-03, ql64ptw3 = 2.7480081019085784987954e-03;
const Real ql64ptx4 = 6.4036616986298811049494e-03, ql64ptw4 = 3.9606448950233168887491e-03;
const Real ql64ptx5 = 1.0966668584301952904525e-02, ql64ptw5 = 5.1635011834076640843683e-03;
const Real ql64ptx6 = 1.6726444815450380865229e-02, ql64ptw6 = 6.3536995987273680100427e-03;
const Real ql64ptx7 = 2.3668882105668531679044e-02, ql64ptw7 = 7.5283419939807215415239e-03;
const Real ql64ptx8 = 3.1776986262182893572259e-02, ql64ptw8 = 8.6845581922710916172869e-03;
const Real ql64ptx9 = 4.1030913244859190669445e-02, ql64ptw9 = 9.8195203616208599772808e-03;
const Real ql64ptx10 = 5.1408016077074958349868e-02, ql64ptw10 = 1.0930451755759030760262e-02;
const Real ql64ptx11 = 6.2882899671186265333489e-02, ql64ptw11 = 1.2014634072011913523692e-02;
const Real ql64ptx12 = 7.5427482728504524622792e-02, ql64ptw12 = 1.3069414307169218564253e-02;
const Real ql64ptx13 = 8.9011066346241451974208e-02, ql64ptw13 = 1.4092211332924257927157e-02;
const Real ql64ptx14 = 1.0360040908689593930347e-01, ql64ptw14 = 1.5080522249544725879589e-02;
const Real ql64ptx15 = 1.1915980829594002843308e-01, ql64ptw15 = 1.6031928528863512856129e-02;
const Real ql64ptx16 = 1.3565118745558152335562e-01, ql64ptw16 = 1.6944101942062698290536e-02;
const Real ql64ptx17 = 1.5303418935464757622000e-01, ql64ptw17 = 1.7814810262244743649385e-02;
const Real ql64ptx18 = 1.7126627484351175656485e-01, ql64ptw18 = 1.8641922729900586158669e-02;
const Real ql64ptx19 = 1.9030282693078420797050e-01, ql64ptw19 = 1.9423415268903867353378e-02;
const Real ql64ptx20 = 2.1009725996614081466873e-01, ql64ptw20 = 2.0157375440780118946993e-02;
const Real ql64ptx21 = 2.3060113364159978699419e-01, ql64ptw21 = 2.0842007125400976996765e-02;
const Real ql64ptx22 = 2.5176427153197111774574e-01, ql64ptw22 = 2.1475634916800909884893e-02;
const Real ql64ptx23 = 2.7353488388420943433488e-01, ql64ptw23 = 2.2056708223446234995446e-02;
const Real ql64ptx24 = 2.9585969435507297742305e-01, ql64ptw24 = 2.2583805062973850291685e-02;
const Real ql64ptx25 = 3.1868407038686907828406e-01, ql64ptw25 = 2.3055635542144528593589e-02;
const Real ql64ptx26 = 3.4195215690218710991033e-01, ql64ptw26 = 2.3471045013514157950851e-02;
const Real ql64ptx27 = 3.6560701299041498124609e-01, ql64ptw27 = 2.3829016901110319781587e-02;
const Real ql64ptx28 = 3.8959075125152325957956e-01, ql64ptw28 = 2.4128675188207274887597e-02;
const Real ql64ptx29 = 4.1384467945610348138530e-01, ql64ptw29 = 2.4369286561116591111054e-02;
const Real ql64ptx30 = 4.3830944417498601817229e-01, ql64ptw30 = 2.4550262203750655015666e-02;
const Real ql64ptx31 = 4.6292517602694205347547e-01, ql64ptw31 = 2.4671159238569786986695e-02;
const Real ql64ptx32 = 4.8763163618902061191562e-01, ql64ptw32 = 2.4731681810388324305183e-02;
const Real ql64ptx33 = 5.1236836381097938808438e-01, ql64ptw33 = 2.4731681810388324305183e-02;
const Real ql64ptx34 = 5.3707482397305794652453e-01, ql64ptw34 = 2.4671159238569786986695e-02;
const Real ql64ptx35 = 5.6169055582501403733886e-01, ql64ptw35 = 2.4550262203750655015666e-02;
const Real ql64ptx36 = 5.8615532054389651861470e-01, ql64ptw36 = 2.4369286561116591111054e-02;
const Real ql64ptx37 = 6.1040924874847679593159e-01, ql64ptw37 = 2.4128675188207274887597e-02;
const Real ql64ptx38 = 6.3439298700958501875391e-01, ql64ptw38 = 2.3829016901110319781587e-02;
const Real ql64ptx39 = 6.5804784309781294560082e-01, ql64ptw39 = 2.3471045013514157950851e-02;
const Real ql64ptx40 = 6.8131592961313092171594e-01, ql64ptw40 = 2.3055635542144528593589e-02;
const Real ql64ptx41 = 7.0414030564492702257695e-01, ql64ptw41 = 2.2583805062973850291685e-02;
const Real ql64ptx42 = 7.2646511611579056566512e-01, ql64ptw42 = 2.2056708223446234995446e-02;
const Real ql64ptx43 = 7.4823572846802888225426e-01, ql64ptw43 = 2.1475634916800909884893e-02;
const Real ql64ptx44 = 7.6939886635840015749466e-01, ql64ptw44 = 2.0842007125400976996765e-02;
const Real ql64ptx45 = 7.8990274003385918533127e-01, ql64ptw45 = 2.0157375440780118946993e-02;
const Real ql64ptx46 = 8.0969717306921573651834e-01, ql64ptw46 = 1.9423415268903867353378e-02;
const Real ql64ptx47 = 8.2873372515648824343515e-01, ql64ptw47 = 1.8641922729900586158669e-02;
const Real ql64ptx48 = 8.4696581064535236826885e-01, ql64ptw48 = 1.7814810262244743649385e-02;
const Real ql64ptx49 = 8.6434881254441853215553e-01, ql64ptw49 = 1.6944101942062698290536e-02;
const Real ql64ptx50 = 8.8084019170406002707807e-01, ql64ptw50 = 1.6031928528863512856129e-02;
const Real ql64ptx51 = 8.9639959091310406069653e-01, ql64ptw51 = 1.5080522249544725879589e-02;
const Real ql64ptx52 = 9.1098893365375854802579e-01, ql64ptw52 = 1.4092211332924257927157e-02;
const Real ql64ptx53 = 9.2457251727149547537721e-01, ql64ptw53 = 1.3069414307169218564253e-02;
const Real ql64ptx54 = 9.3711710032881367915536e-01, ql64ptw54 = 1.2014634072011913523692e-02;
const Real ql64ptx55 = 9.4859198392292509716128e-01, ql64ptw55 = 1.0930451755759030760262e-02;
const Real ql64ptx56 = 9.5896908675514080933056e-01, ql64ptw56 = 9.8195203616208599772808e-03;
const Real ql64ptx57 = 9.6822301373781716193889e-01, ql64ptw57 = 8.6845581922710916172869e-03;
const Real ql64ptx58 = 9.7633111789433146832096e-01, ql64ptw58 = 7.5283419939807215415239e-03;
const Real ql64ptx59 = 9.8327355518454961913477e-01, ql64ptw59 = 6.3536995987273680100427e-03;
const Real ql64ptx60 = 9.8903333141569804709547e-01, ql64ptw60 = 5.1635011834076640843683e-03;
const Real ql64ptx61 = 9.9359633830137017440620e-01, ql64ptw61 = 3.9606448950233168887491e-03;
const Real ql64ptx62 = 9.9695136335152856332797e-01, ql64ptw62 = 2.7480081019085784987954e-03;
const Real ql64ptx63 = 9.9908993575108162588094e-01, ql64ptw63 = 1.5280041224562451382396e-03;
const Real ql64ptx64 = 1.0000000000000000000000e+00, ql64ptw64 = 2.4801587301587300210537e-04;


static const QuadraturePoint pt_lob_seg_64pt[64] =
{

	QuadraturePoint ( ql64ptx1, ql64ptw1),
	QuadraturePoint ( ql64ptx2, ql64ptw2),
	QuadraturePoint ( ql64ptx3, ql64ptw3),
	QuadraturePoint ( ql64ptx4, ql64ptw4),
	QuadraturePoint ( ql64ptx5, ql64ptw5),
	QuadraturePoint ( ql64ptx6, ql64ptw6),
	QuadraturePoint ( ql64ptx7, ql64ptw7),
	QuadraturePoint ( ql64ptx8, ql64ptw8),
	QuadraturePoint ( ql64ptx9, ql64ptw9),
	QuadraturePoint ( ql64ptx10, ql64ptw10),
	QuadraturePoint ( ql64ptx11, ql64ptw11),
	QuadraturePoint ( ql64ptx12, ql64ptw12),
	QuadraturePoint ( ql64ptx13, ql64ptw13),
	QuadraturePoint ( ql64ptx14, ql64ptw14),
	QuadraturePoint ( ql64ptx15, ql64ptw15),
	QuadraturePoint ( ql64ptx16, ql64ptw16),
	QuadraturePoint ( ql64ptx17, ql64ptw17),
	QuadraturePoint ( ql64ptx18, ql64ptw18),
	QuadraturePoint ( ql64ptx19, ql64ptw19),
	QuadraturePoint ( ql64ptx20, ql64ptw20),
	QuadraturePoint ( ql64ptx21, ql64ptw21),
	QuadraturePoint ( ql64ptx22, ql64ptw22),
	QuadraturePoint ( ql64ptx23, ql64ptw23),
	QuadraturePoint ( ql64ptx24, ql64ptw24),
	QuadraturePoint ( ql64ptx25, ql64ptw25),
	QuadraturePoint ( ql64ptx26, ql64ptw26),
	QuadraturePoint ( ql64ptx27, ql64ptw27),
	QuadraturePoint ( ql64ptx28, ql64ptw28),
	QuadraturePoint ( ql64ptx29, ql64ptw29),
	QuadraturePoint ( ql64ptx30, ql64ptw30),
	QuadraturePoint ( ql64ptx31, ql64ptw31),
	QuadraturePoint ( ql64ptx32, ql64ptw32),
	QuadraturePoint ( ql64ptx33, ql64ptw33),
	QuadraturePoint ( ql64ptx34, ql64ptw34),
	QuadraturePoint ( ql64ptx35, ql64ptw35),
	QuadraturePoint ( ql64ptx36, ql64ptw36),
	QuadraturePoint ( ql64ptx37, ql64ptw37),
	QuadraturePoint ( ql64ptx38, ql64ptw38),
	QuadraturePoint ( ql64ptx39, ql64ptw39),
	QuadraturePoint ( ql64ptx40, ql64ptw40),
	QuadraturePoint ( ql64ptx41, ql64ptw41),
	QuadraturePoint ( ql64ptx42, ql64ptw42),
	QuadraturePoint ( ql64ptx43, ql64ptw43),
	QuadraturePoint ( ql64ptx44, ql64ptw44),
	QuadraturePoint ( ql64ptx45, ql64ptw45),
	QuadraturePoint ( ql64ptx46, ql64ptw46),
	QuadraturePoint ( ql64ptx47, ql64ptw47),
	QuadraturePoint ( ql64ptx48, ql64ptw48),
	QuadraturePoint ( ql64ptx49, ql64ptw49),
	QuadraturePoint ( ql64ptx50, ql64ptw50),
	QuadraturePoint ( ql64ptx51, ql64ptw51),
	QuadraturePoint ( ql64ptx52, ql64ptw52),
	QuadraturePoint ( ql64ptx53, ql64ptw53),
	QuadraturePoint ( ql64ptx54, ql64ptw54),
	QuadraturePoint ( ql64ptx55, ql64ptw55),
	QuadraturePoint ( ql64ptx56, ql64ptw56),
	QuadraturePoint ( ql64ptx57, ql64ptw57),
	QuadraturePoint ( ql64ptx58, ql64ptw58),
	QuadraturePoint ( ql64ptx59, ql64ptw59),
	QuadraturePoint ( ql64ptx60, ql64ptw60),
	QuadraturePoint ( ql64ptx61, ql64ptw61),
	QuadraturePoint ( ql64ptx62, ql64ptw62),
	QuadraturePoint ( ql64ptx63, ql64ptw63),
	QuadraturePoint ( ql64ptx64, ql64ptw64)
};

const QuadratureRule quadRuleLobSeg64pt ( pt_lob_seg_64pt,
                                      QUAD_RULE_LOB_SEG_64PT,
                                      "Gauss 64 points on a segment", LINE, 64, 126 );

/*----------------------------------------------------------------------
  Set of all quadrature rules on segments
  ----------------------------------------------------------------------*/
static const QuadratureRule quad_rule_seg[ NB_QUAD_RULE_SEG ] =
{
    quadRuleSeg1pt,
    quadRuleSeg2pt,
    quadRuleSeg3pt,
    quadRuleSeg4pt,
    quadRuleSeg32pt,
	quadRuleLobSeg32pt,
	quadRuleLobSeg64pt
};
/*======================================================================
 *
 *                     Quadrature Rules 2D on triangles
 *
 *=======================================================================*/
//! total number of quadrature rules in 2D on triangle
#define NB_QUAD_RULE_TRIA 5
//! id of the quadrature rules on triangles
#define QUAD_RULE_TRIA_1PT     1
#define QUAD_RULE_TRIA_3PT     2
#define QUAD_RULE_TRIA_4PT     3
#define QUAD_RULE_TRIA_6PT     4
#define QUAD_RULE_TRIA_7PT     5
//----------------------------------------------------------------------

static const QuadraturePoint pt_tria_1pt[ 1 ] =
{
    QuadraturePoint ( 1. / 3., 1. / 3., 1. / 2. )
};
const QuadratureRule quadRuleTria1pt ( pt_tria_1pt,
                                       QUAD_RULE_TRIA_1PT,
                                       "Quadrature rule 1 point on a triangle", TRIANGLE, 1, 1 );
//----------------------------------------------------------------------
static const QuadraturePoint pt_tria_3pt[ 3 ] =
{
    QuadraturePoint ( 0.5, 0. , 1. / 6. ),
    QuadraturePoint ( 0. , 0.5, 1. / 6. ),
    QuadraturePoint ( 0.5, 0.5, 1. / 6. )
};
const QuadratureRule quadRuleTria3pt ( pt_tria_3pt,
                                       QUAD_RULE_TRIA_3PT,
                                       "Quadrature rule 3 points on a triangle", TRIANGLE, 3, 2 );
//----------------------------------------------------------------------
// 4 points Integration rule for triangle (Ref. e.g. Comincioli pag. 234) D of Ex = 3
const Real t4pt_xb1 = 3. / 5.,
           t4pt_xb2 =  1. / 5.,
           t4pt_w1  = 25. / 96.,
           t4pt_w2  = -9. / 32.,
           t4pt_a   =  1. / 3.;

static const QuadraturePoint pt_tria_4pt[ 4 ] =
{
    QuadraturePoint ( t4pt_xb1, t4pt_xb2, t4pt_w1 ),
    QuadraturePoint ( t4pt_xb2, t4pt_xb1, t4pt_w1 ),
    QuadraturePoint ( t4pt_xb2, t4pt_xb2, t4pt_w1 ),
    QuadraturePoint ( t4pt_a, t4pt_a, t4pt_w2 )
};

const QuadratureRule quadRuleTria4pt ( pt_tria_4pt,
                                       QUAD_RULE_TRIA_4PT,
                                       "Quadrature rule 4 points on a triangle", TRIANGLE, 4, 3 );
//----------------------------------------------------------------------
// 6 points Integration rule for triangle, D of Ex = 4
// Ref: G.R. Cowper,  Gaussian quadrature formulas for triangles,
//      Internat. J. Numer. Methods Engrg.  7 (1973), 405--408.
const Real t6pt_x1 = 0.091576213509770743;
const Real t6pt_x2 = 0.44594849091596488;
const Real t6pt_w1 = 0.054975871827660933;
const Real t6pt_w2 = 0.11169079483900573;

static const QuadraturePoint pt_tria_6pt[ 6 ] =
{
    QuadraturePoint (     t6pt_x1,     t6pt_x1, t6pt_w1 ),
    QuadraturePoint (     t6pt_x1, 1 - 2 * t6pt_x1, t6pt_w1 ),
    QuadraturePoint ( 1 - 2 * t6pt_x1,     t6pt_x1, t6pt_w1 ),
    QuadraturePoint (     t6pt_x2,     t6pt_x2, t6pt_w2 ),
    QuadraturePoint (     t6pt_x2, 1 - 2 * t6pt_x2, t6pt_w2 ),
    QuadraturePoint ( 1 - 2 * t6pt_x2,     t6pt_x2, t6pt_w2 ),
};
const QuadratureRule quadRuleTria6pt ( pt_tria_6pt, QUAD_RULE_TRIA_6PT,
                                       "Quadrature rule 6 points on a triangle",
                                       TRIANGLE, 6, 4 );
//----------------------------------------------------------------------
// 7 points Integration rule for triangle (Ref. Stroud) D of Ex = 5
const Real t7pt_x0 = 1. / 3.;
const Real t7pt_x1 = 0.10128650732345633;
const Real t7pt_x2 = 0.47014206410511508;
const Real t7pt_w0 = 0.1125;
const Real t7pt_w1 = 0.062969590272413576;
const Real t7pt_w2 = 0.066197076394253090;

static const QuadraturePoint pt_tria_7pt[ 7 ] =
{
    QuadraturePoint (     t7pt_x0,     t7pt_x0, t7pt_w0 ),
    QuadraturePoint (     t7pt_x1,     t7pt_x1, t7pt_w1 ),
    QuadraturePoint (     t7pt_x1, 1 - 2 * t7pt_x1, t7pt_w1 ),
    QuadraturePoint ( 1 - 2 * t7pt_x1,     t7pt_x1, t7pt_w1 ),
    QuadraturePoint (     t7pt_x2,     t7pt_x2, t7pt_w2 ),
    QuadraturePoint (     t7pt_x2, 1 - 2 * t7pt_x2, t7pt_w2 ),
    QuadraturePoint ( 1 - 2 * t7pt_x2,     t7pt_x2, t7pt_w2 ),
};
const QuadratureRule quadRuleTria7pt ( pt_tria_7pt, QUAD_RULE_TRIA_7PT,
                                       "Quadrature rule 7 points on a triangle",
                                       TRIANGLE, 7, 5 );
/*----------------------------------------------------------------------
  Set of all quadrature rules on triangle
  ----------------------------------------------------------------------*/
static const QuadratureRule quad_rule_tria[ NB_QUAD_RULE_TRIA ] =
{
    quadRuleTria1pt,
    quadRuleTria3pt,
    quadRuleTria4pt,
    quadRuleTria6pt,
    quadRuleTria7pt
};
//----------------------------------------------------------------------
/*======================================================================
 *
 *                     Quadrature Rules 2D on quadrangles
 *
 *=======================================================================*/
//! total number of quadrature rules in 2D on quadrangle
#define NB_QUAD_RULE_QUAD 4
//! id of the quadrature rules on quadrangles
#define QUAD_RULE_QUAD_1PT     1
#define QUAD_RULE_QUAD_4PT     2
#define QUAD_RULE_QUAD_9PT     3
#define QUAD_RULE_QUAD_16PT    4
//----------------------------------------------------------------------

static const QuadraturePoint pt_quad_1pt[ 1 ] =
{
    QuadraturePoint ( .5, .5, 1. )
};
const QuadratureRule quadRuleQuad1pt ( pt_quad_1pt,
                                       QUAD_RULE_QUAD_1PT,
                                       "Quadrature rule 1 point on a quadrangle", QUAD, 1, 1 );
//----------------------------------------------------------------------
static const QuadraturePoint pt_quad_4pt[ 4 ] =
{
    QuadraturePoint ( q2ptx1, q2ptx1, q2ptw1 * q2ptw1 ),
    QuadraturePoint ( q2ptx1, q2ptx2, q2ptw1 * q2ptw2 ),
    QuadraturePoint ( q2ptx2, q2ptx1, q2ptw2 * q2ptw1 ),
    QuadraturePoint ( q2ptx2, q2ptx2, q2ptw2 * q2ptw2 )
};
const QuadratureRule quadRuleQuad4pt ( pt_quad_4pt,
                                       QUAD_RULE_QUAD_4PT,
                                       "Quadrature rule 4 points on a quadrangle", QUAD, 4, 3 );
//----------------------------------------------------------------------
// 4 points Integration rule for quadrangle

static const QuadraturePoint pt_quad_9pt[ 9 ] =
{
    QuadraturePoint ( q3ptx1, q3ptx1, q3ptw1 * q3ptw1 ),
    QuadraturePoint ( q3ptx2, q3ptx1, q3ptw2 * q3ptw1 ),
    QuadraturePoint ( q3ptx3, q3ptx1, q3ptw3 * q3ptw1 ),
    QuadraturePoint ( q3ptx1, q3ptx2, q3ptw1 * q3ptw2 ),
    QuadraturePoint ( q3ptx2, q3ptx2, q3ptw2 * q3ptw2 ),
    QuadraturePoint ( q3ptx3, q3ptx2, q3ptw3 * q3ptw2 ),
    QuadraturePoint ( q3ptx1, q3ptx3, q3ptw1 * q3ptw3 ),
    QuadraturePoint ( q3ptx2, q3ptx3, q3ptw2 * q3ptw3 ),
    QuadraturePoint ( q3ptx3, q3ptx3, q3ptw3 * q3ptw3 )
};

const QuadratureRule quadRuleQuad9pt ( pt_quad_9pt,
                                       QUAD_RULE_QUAD_9PT,
                                       "Quadrature rule 9 points on a quadrangle", QUAD, 9, 5 );
//----------------------------------------------------------------------
// 4 points Integration rule for quadrangle

static const QuadraturePoint pt_quad_16pt[ 16 ] =
{
    QuadraturePoint ( q4ptx1, q4ptx1, q4ptw1 * q4ptw1 ),
    QuadraturePoint ( q4ptx2, q4ptx1, q4ptw2 * q4ptw1 ),
    QuadraturePoint ( q4ptx3, q4ptx1, q4ptw3 * q4ptw1 ),
    QuadraturePoint ( q4ptx4, q4ptx1, q4ptw4 * q4ptw1 ),
    QuadraturePoint ( q4ptx1, q4ptx2, q4ptw1 * q4ptw2 ),
    QuadraturePoint ( q4ptx2, q4ptx2, q4ptw2 * q4ptw2 ),
    QuadraturePoint ( q4ptx3, q4ptx2, q4ptw3 * q4ptw2 ),
    QuadraturePoint ( q4ptx4, q4ptx2, q4ptw4 * q4ptw2 ),
    QuadraturePoint ( q4ptx1, q4ptx3, q4ptw1 * q4ptw3 ),
    QuadraturePoint ( q4ptx2, q4ptx3, q4ptw2 * q4ptw3 ),
    QuadraturePoint ( q4ptx3, q4ptx3, q4ptw3 * q4ptw3 ),
    QuadraturePoint ( q4ptx4, q4ptx3, q4ptw4 * q4ptw3 ),
    QuadraturePoint ( q4ptx1, q4ptx4, q4ptw1 * q4ptw4 ),
    QuadraturePoint ( q4ptx2, q4ptx4, q4ptw2 * q4ptw4 ),
    QuadraturePoint ( q4ptx3, q4ptx4, q4ptw3 * q4ptw4 ),
    QuadraturePoint ( q4ptx4, q4ptx4, q4ptw4 * q4ptw4 )
};

const QuadratureRule quadRuleQuad16pt ( pt_quad_16pt,
                                        QUAD_RULE_QUAD_16PT,
                                        "Quadrature rule 16 points on a quadrangle", QUAD, 16, 5 );
/*----------------------------------------------------------------------
  Set of all quadrature rules on quadrangle
  ----------------------------------------------------------------------*/
static const QuadratureRule quad_rule_quad[ NB_QUAD_RULE_QUAD ] =
{
    quadRuleQuad1pt,
    quadRuleQuad4pt,
    quadRuleQuad9pt,
    quadRuleQuad16pt
};
//----------------------------------------------------------------------
/*======================================================================
 *
 *                     Quadrature Rules 3D on tetraedras
 *
 *=======================================================================*/
//! total number of quadrature rules in 3D on tetraedra
#define NB_QUAD_RULE_TETRA 5
//! id of the quadrature rules on tetraedra
#define QUAD_RULE_TETRA_1PT     1
#define QUAD_RULE_TETRA_4PT     2
#define QUAD_RULE_TETRA_5PT     3
#define QUAD_RULE_TETRA_15PT    4
#define QUAD_RULE_TETRA_64PT    5
//----------------------------------------------------------------------

static const QuadraturePoint pt_tetra_1pt[ 1 ] =
{
    QuadraturePoint ( 1. / 4., 1. / 4., 1. / 4., 1. / 6. )
};
const QuadratureRule quadRuleTetra1pt ( pt_tetra_1pt,
                                        QUAD_RULE_TETRA_1PT,
                                        "Quadrature rule 1 point on a tetraedra", TETRA, 1, 1 );
//----------------------------------------------------------------------
const Real tet4ptx1 = ( 5. - std::sqrt ( 5. ) ) / 20., tet4ptx2 = ( 5. + 3 * std::sqrt ( 5. ) ) / 20.;

static const QuadraturePoint pt_tetra_4pt[ 4 ] =
{
    QuadraturePoint ( tet4ptx1, tet4ptx1, tet4ptx1, 1. / 24. ),
    QuadraturePoint ( tet4ptx1, tet4ptx1, tet4ptx2, 1. / 24. ),
    QuadraturePoint ( tet4ptx1, tet4ptx2, tet4ptx1, 1. / 24. ),
    QuadraturePoint ( tet4ptx2, tet4ptx1, tet4ptx1, 1. / 24. )
};
const QuadratureRule quadRuleTetra4pt ( pt_tetra_4pt,
                                        QUAD_RULE_TETRA_4PT,
                                        "Quadrature rule 4 points on a tetraedra", TETRA, 4, 2 );
//----------------------------------------------------------------------
// 5 points Integration rule for tetraedra (Ref. e.g. Comincioli pag. 236)
const Real tet5ptx1 = 1. / 6. , tet5ptx2 = 1. / 2., tet5ptx3 = 1. / 4.;

static const QuadraturePoint pt_tetra_5pt[ 5 ] =
{
    QuadraturePoint ( tet5ptx1, tet5ptx1, tet5ptx1, 9. / 120. ),
    QuadraturePoint ( tet5ptx1, tet5ptx1, tet5ptx2, 9. / 120. ),
    QuadraturePoint ( tet5ptx1, tet5ptx2, tet5ptx1, 9. / 120. ),
    QuadraturePoint ( tet5ptx2, tet5ptx1, tet5ptx1, 9. / 120. ),
    QuadraturePoint ( tet5ptx3, tet5ptx3, tet5ptx3, -16. / 120. )
};

const QuadratureRule quadRuleTetra5pt ( pt_tetra_5pt,
                                        QUAD_RULE_TETRA_5PT,
                                        "Quadrature rule 5 points on a tetraedra", TETRA, 5, 3 );
//
//----------------------------------------------------------------------
//                     15 points integration rule for tetra.
//                   D o E = 5 (Stroud, T3:5-1 pag. 315)
// r
const Real r5 = 0.25;
// s
const Real s5[ 4 ] =
{
    0.09197107805272303, 0.3197936278296299
}
; // (7 \mp \sqrt(15))/34
// t
const Real t5[ 4 ] =
{
    0.7240867658418310, 0.04061911651111023
}
; // (13 \pm 3*sqrt(15))/34
// u
const Real u5 = 0.05635083268962915; // (10-2*sqrt(15))/40
// v
const Real v5 = 0.4436491673103708; // (10+2*sqrt(15))/40
// A
const Real A5 = 0.01975308641975309; // 16/135*1/6
// B
const Real B5[ 2 ] =
{
    0.01198951396316977, 0.01151136787104540
}
; // 1/6*(2665 \pm 14*sqrt(15))/37800
// C
const Real C5 = 0.008818342151675485; // 20/378*1/6
//
static const QuadraturePoint pt_tetra_15pt[ 15 ] =
{
    QuadraturePoint ( r5, r5, r5, A5 ),
    QuadraturePoint ( s5[ 0 ], s5[ 0 ], s5[ 0 ], B5[ 0 ] ),
    QuadraturePoint ( t5[ 0 ], s5[ 0 ], s5[ 0 ], B5[ 0 ] ),
    QuadraturePoint ( s5[ 0 ], t5[ 0 ], s5[ 0 ], B5[ 0 ] ),
    QuadraturePoint ( s5[ 0 ], s5[ 0 ], t5[ 0 ], B5[ 0 ] ),
    QuadraturePoint ( s5[ 1 ], s5[ 1 ], s5[ 1 ], B5[ 1 ] ),
    QuadraturePoint ( t5[ 1 ], s5[ 1 ], s5[ 1 ], B5[ 1 ] ),
    QuadraturePoint ( s5[ 1 ], t5[ 1 ], s5[ 1 ], B5[ 1 ] ),
    QuadraturePoint ( s5[ 1 ], s5[ 1 ], t5[ 1 ], B5[ 1 ] ),
    QuadraturePoint ( u5, u5, v5, C5 ),
    QuadraturePoint ( u5, v5, u5, C5 ),
    QuadraturePoint ( v5, u5, u5, C5 ),
    QuadraturePoint ( v5, v5, u5, C5 ),
    QuadraturePoint ( v5, u5, v5, C5 ),
    QuadraturePoint ( u5, v5, v5, C5 )
};
//
const QuadratureRule quadRuleTetra15pt ( pt_tetra_15pt,
                                         QUAD_RULE_TETRA_15PT,
                                         "Quadrature rule 15 points on a tetraedra",
                                         TETRA, 15, 5 );
//----------------------------------------------------------------------
//                     64 points integration rule for tetra.
//                   D o E = 7 (Stroud, T3:7-1 pag. 315)
//
// t
const Real t[ 4 ] =
{
    0.0485005494, 0.2386007376, 0.5170472951, 0.7958514179
};
// s
const Real s[ 4 ] =
{
    0.0571041961, 0.2768430136, 0.5835904324, 0.8602401357
};
// r
const Real r[ 4 ] =
{
    0.0694318422, 0.3300094782, 0.6699905218, 0.9305681558
};
// A
const Real A[ 4 ] =
{
    0.1739274226, 0.3260725774, 0.3260725774, 0.1739274226
};
// B
const Real B[ 4 ] =
{
    0.1355069134, 0.2034645680, 0.1298475476, 0.0311809709
};
// C
const Real C[ 4 ] =
{
    0.1108884156, 0.1434587898, 0.0686338872, 0.0103522407
};
//
/*
for (i=0;i<4;i++){
  for (j=0;j<4;j++){
    for (k=0;k<4;k++){
      pt_tetra_64pt[k+4*j+16*i]=
 QuadraturePoint(t[k], s[j]*(1-t[k]), r[i]*(1-s[j])*(1-t[k]), A[i]*B[j]*C[k]);
    }
  }
}
*/
static const QuadraturePoint pt_tetra_64pt[ 64 ] =
{
    QuadraturePoint ( t[ 0 ], s[ 0 ] * ( 1 - t[ 0 ] ), r[ 0 ] * ( 1 - s[ 0 ] ) * ( 1 - t[ 0 ] ), A[ 0 ] * B[ 0 ] * C[ 0 ] ),
    QuadraturePoint ( t[ 1 ], s[ 0 ] * ( 1 - t[ 1 ] ), r[ 0 ] * ( 1 - s[ 0 ] ) * ( 1 - t[ 1 ] ), A[ 0 ] * B[ 0 ] * C[ 1 ] ),
    QuadraturePoint ( t[ 2 ], s[ 0 ] * ( 1 - t[ 2 ] ), r[ 0 ] * ( 1 - s[ 0 ] ) * ( 1 - t[ 2 ] ), A[ 0 ] * B[ 0 ] * C[ 2 ] ),
    QuadraturePoint ( t[ 3 ], s[ 0 ] * ( 1 - t[ 3 ] ), r[ 0 ] * ( 1 - s[ 0 ] ) * ( 1 - t[ 3 ] ), A[ 0 ] * B[ 0 ] * C[ 3 ] ),
    QuadraturePoint ( t[ 0 ], s[ 1 ] * ( 1 - t[ 0 ] ), r[ 0 ] * ( 1 - s[ 1 ] ) * ( 1 - t[ 0 ] ), A[ 0 ] * B[ 1 ] * C[ 0 ] ),
    QuadraturePoint ( t[ 1 ], s[ 1 ] * ( 1 - t[ 1 ] ), r[ 0 ] * ( 1 - s[ 1 ] ) * ( 1 - t[ 1 ] ), A[ 0 ] * B[ 1 ] * C[ 1 ] ),
    QuadraturePoint ( t[ 2 ], s[ 1 ] * ( 1 - t[ 2 ] ), r[ 0 ] * ( 1 - s[ 1 ] ) * ( 1 - t[ 2 ] ), A[ 0 ] * B[ 1 ] * C[ 2 ] ),
    QuadraturePoint ( t[ 3 ], s[ 1 ] * ( 1 - t[ 3 ] ), r[ 0 ] * ( 1 - s[ 1 ] ) * ( 1 - t[ 3 ] ), A[ 0 ] * B[ 1 ] * C[ 3 ] ),
    QuadraturePoint ( t[ 0 ], s[ 2 ] * ( 1 - t[ 0 ] ), r[ 0 ] * ( 1 - s[ 2 ] ) * ( 1 - t[ 0 ] ), A[ 0 ] * B[ 2 ] * C[ 0 ] ),
    QuadraturePoint ( t[ 1 ], s[ 2 ] * ( 1 - t[ 1 ] ), r[ 0 ] * ( 1 - s[ 2 ] ) * ( 1 - t[ 1 ] ), A[ 0 ] * B[ 2 ] * C[ 1 ] ),
    QuadraturePoint ( t[ 2 ], s[ 2 ] * ( 1 - t[ 2 ] ), r[ 0 ] * ( 1 - s[ 2 ] ) * ( 1 - t[ 2 ] ), A[ 0 ] * B[ 2 ] * C[ 2 ] ),
    QuadraturePoint ( t[ 3 ], s[ 2 ] * ( 1 - t[ 3 ] ), r[ 0 ] * ( 1 - s[ 2 ] ) * ( 1 - t[ 3 ] ), A[ 0 ] * B[ 2 ] * C[ 3 ] ),
    QuadraturePoint ( t[ 0 ], s[ 3 ] * ( 1 - t[ 0 ] ), r[ 0 ] * ( 1 - s[ 3 ] ) * ( 1 - t[ 0 ] ), A[ 0 ] * B[ 3 ] * C[ 0 ] ),
    QuadraturePoint ( t[ 1 ], s[ 3 ] * ( 1 - t[ 1 ] ), r[ 0 ] * ( 1 - s[ 3 ] ) * ( 1 - t[ 1 ] ), A[ 0 ] * B[ 3 ] * C[ 1 ] ),
    QuadraturePoint ( t[ 2 ], s[ 3 ] * ( 1 - t[ 2 ] ), r[ 0 ] * ( 1 - s[ 3 ] ) * ( 1 - t[ 2 ] ), A[ 0 ] * B[ 3 ] * C[ 2 ] ),
    QuadraturePoint ( t[ 3 ], s[ 3 ] * ( 1 - t[ 3 ] ), r[ 0 ] * ( 1 - s[ 3 ] ) * ( 1 - t[ 3 ] ), A[ 0 ] * B[ 3 ] * C[ 3 ] ),
    QuadraturePoint ( t[ 0 ], s[ 0 ] * ( 1 - t[ 0 ] ), r[ 1 ] * ( 1 - s[ 0 ] ) * ( 1 - t[ 0 ] ), A[ 1 ] * B[ 0 ] * C[ 0 ] ),
    QuadraturePoint ( t[ 1 ], s[ 0 ] * ( 1 - t[ 1 ] ), r[ 1 ] * ( 1 - s[ 0 ] ) * ( 1 - t[ 1 ] ), A[ 1 ] * B[ 0 ] * C[ 1 ] ),
    QuadraturePoint ( t[ 2 ], s[ 0 ] * ( 1 - t[ 2 ] ), r[ 1 ] * ( 1 - s[ 0 ] ) * ( 1 - t[ 2 ] ), A[ 1 ] * B[ 0 ] * C[ 2 ] ),
    QuadraturePoint ( t[ 3 ], s[ 0 ] * ( 1 - t[ 3 ] ), r[ 1 ] * ( 1 - s[ 0 ] ) * ( 1 - t[ 3 ] ), A[ 1 ] * B[ 0 ] * C[ 3 ] ),
    QuadraturePoint ( t[ 0 ], s[ 1 ] * ( 1 - t[ 0 ] ), r[ 1 ] * ( 1 - s[ 1 ] ) * ( 1 - t[ 0 ] ), A[ 1 ] * B[ 1 ] * C[ 0 ] ),
    QuadraturePoint ( t[ 1 ], s[ 1 ] * ( 1 - t[ 1 ] ), r[ 1 ] * ( 1 - s[ 1 ] ) * ( 1 - t[ 1 ] ), A[ 1 ] * B[ 1 ] * C[ 1 ] ),
    QuadraturePoint ( t[ 2 ], s[ 1 ] * ( 1 - t[ 2 ] ), r[ 1 ] * ( 1 - s[ 1 ] ) * ( 1 - t[ 2 ] ), A[ 1 ] * B[ 1 ] * C[ 2 ] ),
    QuadraturePoint ( t[ 3 ], s[ 1 ] * ( 1 - t[ 3 ] ), r[ 1 ] * ( 1 - s[ 1 ] ) * ( 1 - t[ 3 ] ), A[ 1 ] * B[ 1 ] * C[ 3 ] ),
    QuadraturePoint ( t[ 0 ], s[ 2 ] * ( 1 - t[ 0 ] ), r[ 1 ] * ( 1 - s[ 2 ] ) * ( 1 - t[ 0 ] ), A[ 1 ] * B[ 2 ] * C[ 0 ] ),
    QuadraturePoint ( t[ 1 ], s[ 2 ] * ( 1 - t[ 1 ] ), r[ 1 ] * ( 1 - s[ 2 ] ) * ( 1 - t[ 1 ] ), A[ 1 ] * B[ 2 ] * C[ 1 ] ),
    QuadraturePoint ( t[ 2 ], s[ 2 ] * ( 1 - t[ 2 ] ), r[ 1 ] * ( 1 - s[ 2 ] ) * ( 1 - t[ 2 ] ), A[ 1 ] * B[ 2 ] * C[ 2 ] ),
    QuadraturePoint ( t[ 3 ], s[ 2 ] * ( 1 - t[ 3 ] ), r[ 1 ] * ( 1 - s[ 2 ] ) * ( 1 - t[ 3 ] ), A[ 1 ] * B[ 2 ] * C[ 3 ] ),
    QuadraturePoint ( t[ 0 ], s[ 3 ] * ( 1 - t[ 0 ] ), r[ 1 ] * ( 1 - s[ 3 ] ) * ( 1 - t[ 0 ] ), A[ 1 ] * B[ 3 ] * C[ 0 ] ),
    QuadraturePoint ( t[ 1 ], s[ 3 ] * ( 1 - t[ 1 ] ), r[ 1 ] * ( 1 - s[ 3 ] ) * ( 1 - t[ 1 ] ), A[ 1 ] * B[ 3 ] * C[ 1 ] ),
    QuadraturePoint ( t[ 2 ], s[ 3 ] * ( 1 - t[ 2 ] ), r[ 1 ] * ( 1 - s[ 3 ] ) * ( 1 - t[ 2 ] ), A[ 1 ] * B[ 3 ] * C[ 2 ] ),
    QuadraturePoint ( t[ 3 ], s[ 3 ] * ( 1 - t[ 3 ] ), r[ 1 ] * ( 1 - s[ 3 ] ) * ( 1 - t[ 3 ] ), A[ 1 ] * B[ 3 ] * C[ 3 ] ),
    QuadraturePoint ( t[ 0 ], s[ 0 ] * ( 1 - t[ 0 ] ), r[ 2 ] * ( 1 - s[ 0 ] ) * ( 1 - t[ 0 ] ), A[ 2 ] * B[ 0 ] * C[ 0 ] ),
    QuadraturePoint ( t[ 1 ], s[ 0 ] * ( 1 - t[ 1 ] ), r[ 2 ] * ( 1 - s[ 0 ] ) * ( 1 - t[ 1 ] ), A[ 2 ] * B[ 0 ] * C[ 1 ] ),
    QuadraturePoint ( t[ 2 ], s[ 0 ] * ( 1 - t[ 2 ] ), r[ 2 ] * ( 1 - s[ 0 ] ) * ( 1 - t[ 2 ] ), A[ 2 ] * B[ 0 ] * C[ 2 ] ),
    QuadraturePoint ( t[ 3 ], s[ 0 ] * ( 1 - t[ 3 ] ), r[ 2 ] * ( 1 - s[ 0 ] ) * ( 1 - t[ 3 ] ), A[ 2 ] * B[ 0 ] * C[ 3 ] ),
    QuadraturePoint ( t[ 0 ], s[ 1 ] * ( 1 - t[ 0 ] ), r[ 2 ] * ( 1 - s[ 1 ] ) * ( 1 - t[ 0 ] ), A[ 2 ] * B[ 1 ] * C[ 0 ] ),
    QuadraturePoint ( t[ 1 ], s[ 1 ] * ( 1 - t[ 1 ] ), r[ 2 ] * ( 1 - s[ 1 ] ) * ( 1 - t[ 1 ] ), A[ 2 ] * B[ 1 ] * C[ 1 ] ),
    QuadraturePoint ( t[ 2 ], s[ 1 ] * ( 1 - t[ 2 ] ), r[ 2 ] * ( 1 - s[ 1 ] ) * ( 1 - t[ 2 ] ), A[ 2 ] * B[ 1 ] * C[ 2 ] ),
    QuadraturePoint ( t[ 3 ], s[ 1 ] * ( 1 - t[ 3 ] ), r[ 2 ] * ( 1 - s[ 1 ] ) * ( 1 - t[ 3 ] ), A[ 2 ] * B[ 1 ] * C[ 3 ] ),
    QuadraturePoint ( t[ 0 ], s[ 2 ] * ( 1 - t[ 0 ] ), r[ 2 ] * ( 1 - s[ 2 ] ) * ( 1 - t[ 0 ] ), A[ 2 ] * B[ 2 ] * C[ 0 ] ),
    QuadraturePoint ( t[ 1 ], s[ 2 ] * ( 1 - t[ 1 ] ), r[ 2 ] * ( 1 - s[ 2 ] ) * ( 1 - t[ 1 ] ), A[ 2 ] * B[ 2 ] * C[ 1 ] ),
    QuadraturePoint ( t[ 2 ], s[ 2 ] * ( 1 - t[ 2 ] ), r[ 2 ] * ( 1 - s[ 2 ] ) * ( 1 - t[ 2 ] ), A[ 2 ] * B[ 2 ] * C[ 2 ] ),
    QuadraturePoint ( t[ 3 ], s[ 2 ] * ( 1 - t[ 3 ] ), r[ 2 ] * ( 1 - s[ 2 ] ) * ( 1 - t[ 3 ] ), A[ 2 ] * B[ 2 ] * C[ 3 ] ),
    QuadraturePoint ( t[ 0 ], s[ 3 ] * ( 1 - t[ 0 ] ), r[ 2 ] * ( 1 - s[ 3 ] ) * ( 1 - t[ 0 ] ), A[ 2 ] * B[ 3 ] * C[ 0 ] ),
    QuadraturePoint ( t[ 1 ], s[ 3 ] * ( 1 - t[ 1 ] ), r[ 2 ] * ( 1 - s[ 3 ] ) * ( 1 - t[ 1 ] ), A[ 2 ] * B[ 3 ] * C[ 1 ] ),
    QuadraturePoint ( t[ 2 ], s[ 3 ] * ( 1 - t[ 2 ] ), r[ 2 ] * ( 1 - s[ 3 ] ) * ( 1 - t[ 2 ] ), A[ 2 ] * B[ 3 ] * C[ 2 ] ),
    QuadraturePoint ( t[ 3 ], s[ 3 ] * ( 1 - t[ 3 ] ), r[ 2 ] * ( 1 - s[ 3 ] ) * ( 1 - t[ 3 ] ), A[ 2 ] * B[ 3 ] * C[ 3 ] ),
    QuadraturePoint ( t[ 0 ], s[ 0 ] * ( 1 - t[ 0 ] ), r[ 3 ] * ( 1 - s[ 0 ] ) * ( 1 - t[ 0 ] ), A[ 3 ] * B[ 0 ] * C[ 0 ] ),
    QuadraturePoint ( t[ 1 ], s[ 0 ] * ( 1 - t[ 1 ] ), r[ 3 ] * ( 1 - s[ 0 ] ) * ( 1 - t[ 1 ] ), A[ 3 ] * B[ 0 ] * C[ 1 ] ),
    QuadraturePoint ( t[ 2 ], s[ 0 ] * ( 1 - t[ 2 ] ), r[ 3 ] * ( 1 - s[ 0 ] ) * ( 1 - t[ 2 ] ), A[ 3 ] * B[ 0 ] * C[ 2 ] ),
    QuadraturePoint ( t[ 3 ], s[ 0 ] * ( 1 - t[ 3 ] ), r[ 3 ] * ( 1 - s[ 0 ] ) * ( 1 - t[ 3 ] ), A[ 3 ] * B[ 0 ] * C[ 3 ] ),
    QuadraturePoint ( t[ 0 ], s[ 1 ] * ( 1 - t[ 0 ] ), r[ 3 ] * ( 1 - s[ 1 ] ) * ( 1 - t[ 0 ] ), A[ 3 ] * B[ 1 ] * C[ 0 ] ),
    QuadraturePoint ( t[ 1 ], s[ 1 ] * ( 1 - t[ 1 ] ), r[ 3 ] * ( 1 - s[ 1 ] ) * ( 1 - t[ 1 ] ), A[ 3 ] * B[ 1 ] * C[ 1 ] ),
    QuadraturePoint ( t[ 2 ], s[ 1 ] * ( 1 - t[ 2 ] ), r[ 3 ] * ( 1 - s[ 1 ] ) * ( 1 - t[ 2 ] ), A[ 3 ] * B[ 1 ] * C[ 2 ] ),
    QuadraturePoint ( t[ 3 ], s[ 1 ] * ( 1 - t[ 3 ] ), r[ 3 ] * ( 1 - s[ 1 ] ) * ( 1 - t[ 3 ] ), A[ 3 ] * B[ 1 ] * C[ 3 ] ),
    QuadraturePoint ( t[ 0 ], s[ 2 ] * ( 1 - t[ 0 ] ), r[ 3 ] * ( 1 - s[ 2 ] ) * ( 1 - t[ 0 ] ), A[ 3 ] * B[ 2 ] * C[ 0 ] ),
    QuadraturePoint ( t[ 1 ], s[ 2 ] * ( 1 - t[ 1 ] ), r[ 3 ] * ( 1 - s[ 2 ] ) * ( 1 - t[ 1 ] ), A[ 3 ] * B[ 2 ] * C[ 1 ] ),
    QuadraturePoint ( t[ 2 ], s[ 2 ] * ( 1 - t[ 2 ] ), r[ 3 ] * ( 1 - s[ 2 ] ) * ( 1 - t[ 2 ] ), A[ 3 ] * B[ 2 ] * C[ 2 ] ),
    QuadraturePoint ( t[ 3 ], s[ 2 ] * ( 1 - t[ 3 ] ), r[ 3 ] * ( 1 - s[ 2 ] ) * ( 1 - t[ 3 ] ), A[ 3 ] * B[ 2 ] * C[ 3 ] ),
    QuadraturePoint ( t[ 0 ], s[ 3 ] * ( 1 - t[ 0 ] ), r[ 3 ] * ( 1 - s[ 3 ] ) * ( 1 - t[ 0 ] ), A[ 3 ] * B[ 3 ] * C[ 0 ] ),
    QuadraturePoint ( t[ 1 ], s[ 3 ] * ( 1 - t[ 1 ] ), r[ 3 ] * ( 1 - s[ 3 ] ) * ( 1 - t[ 1 ] ), A[ 3 ] * B[ 3 ] * C[ 1 ] ),
    QuadraturePoint ( t[ 2 ], s[ 3 ] * ( 1 - t[ 2 ] ), r[ 3 ] * ( 1 - s[ 3 ] ) * ( 1 - t[ 2 ] ), A[ 3 ] * B[ 3 ] * C[ 2 ] ),
    QuadraturePoint ( t[ 3 ], s[ 3 ] * ( 1 - t[ 3 ] ), r[ 3 ] * ( 1 - s[ 3 ] ) * ( 1 - t[ 3 ] ), A[ 3 ] * B[ 3 ] * C[ 3 ] )
};
//
const QuadratureRule quadRuleTetra64pt ( pt_tetra_64pt,
                                         QUAD_RULE_TETRA_64PT,
                                         "Quadrature rule 64 points on a tetraedra",
                                         TETRA, 64, 7 );
/*----------------------------------------------------------------------
  Set of all quadrature rules on tetraedra
  ----------------------------------------------------------------------*/
static const QuadratureRule quad_rule_tetra[ NB_QUAD_RULE_TETRA ] =
{
    quadRuleTetra1pt,
    quadRuleTetra4pt,
    quadRuleTetra5pt,
    quadRuleTetra15pt,
    quadRuleTetra64pt
};
//----------------------------------------------------------------------
/*======================================================================
 *
 *                     Quadrature Rules 3D on hexaedras
 *
 *=======================================================================*/
//! total number of quadrature rules in 3D on hexa
#define NB_QUAD_RULE_HEXA 2
//! id of the quadrature rules on quadrangles
#define QUAD_RULE_HEXA_1PT     1
#define QUAD_RULE_HEXA_8PT     2
//----------------------------------------------------------------------

static const QuadraturePoint pt_hexa_1pt[ 1 ] =
{
    QuadraturePoint ( .5, .5, .5, 1. )
};
const QuadratureRule quadRuleHexa1pt ( pt_hexa_1pt,
                                       QUAD_RULE_HEXA_1PT,
                                       "Quadrature rule 1 point on a hexa", HEXA, 1, 1 );
//----------------------------------------------------------------------
static const QuadraturePoint pt_hexa_8pt[ 8 ] =
{
    QuadraturePoint ( q2ptx1, q2ptx1, q2ptx1, q2ptw1* q2ptw1 * q2ptw1 ),
    QuadraturePoint ( q2ptx1, q2ptx2, q2ptx1, q2ptw1* q2ptw2 * q2ptw1 ),
    QuadraturePoint ( q2ptx2, q2ptx1, q2ptx1, q2ptw2* q2ptw1 * q2ptw1 ),
    QuadraturePoint ( q2ptx2, q2ptx2, q2ptx1, q2ptw2* q2ptw2 * q2ptw1 ),
    QuadraturePoint ( q2ptx1, q2ptx1, q2ptx2, q2ptw1* q2ptw1 * q2ptw2 ),
    QuadraturePoint ( q2ptx1, q2ptx2, q2ptx2, q2ptw1* q2ptw2 * q2ptw2 ),
    QuadraturePoint ( q2ptx2, q2ptx1, q2ptx2, q2ptw2* q2ptw1 * q2ptw2 ),
    QuadraturePoint ( q2ptx2, q2ptx2, q2ptx2, q2ptw2* q2ptw2 * q2ptw2 )
};
const QuadratureRule quadRuleHexa8pt ( pt_hexa_8pt,
                                       QUAD_RULE_HEXA_8PT,
                                       "Quadrature rule 8 points on a hexa", HEXA, 8, 3 );
/*----------------------------------------------------------------------
  Set of all quadrature rules on hexa
  ----------------------------------------------------------------------*/
static const QuadratureRule quad_rule_hexa[ NB_QUAD_RULE_HEXA ] =
{
    quadRuleHexa1pt,
    quadRuleHexa8pt
};



//----------------------------------------------------------------------
//
//                       GeometricMaps
//
//----------------------------------------------------------------------

const GeometricMap geoLinearNode ( "Mapping of a point", POINT,
                                   1, 1,
                                   fct_P0_0D, derfct_P0_0D, der2fct_P0_0D,
                                   refcoor_P0_0D,
                                   ( GeometricMap* ) NULL );

const GeometricMap geoLinearSeg ( "Linear mapping on a segment", LINE,
                                  2, 1,
                                  fct_P1_1D, derfct_P1_1D, der2fct_P1_1D,
                                  refcoor_P1_1D,
                                  &geoLinearNode );

const GeometricMap geoQuadraticSeg ( "Quadratic mapping on a segment", LINE,
                                     3, 1,
                                     fct_P2_1D, derfct_P2_1D, der2fct_P2_1D,
                                     refcoor_P2_1D,
                                     &geoLinearNode );

const GeometricMap geoLinearTria ( "Linear mapping on a triangle", TRIANGLE,
                                   3, 2,
                                   fct_P1_2D, derfct_P1_2D, der2fct_P1_2D,
                                   refcoor_P1_2D,
                                   &geoLinearSeg );

const GeometricMap geoBilinearQuad ( "Bilinear mapping on a quadrangle", QUAD,
                                     4, 2,
                                     fct_Q1_2D, derfct_Q1_2D, der2fct_Q1_2D,
                                     refcoor_Q1_2D,
                                     &geoLinearSeg );

const GeometricMap geoBiquadraticQuad ( "Biquadratic mapping on a quadrangle", QUAD,
                                        9, 2,
                                        fct_Q2_2D, derfct_Q2_2D, der2fct_Q2_2D,
                                        refcoor_Q2_2D,
                                        &geoQuadraticSeg );

const GeometricMap geoLinearTetra ( "Linear mapping on a tetraedra", TETRA,
                                    4, 3,
                                    fct_P1_3D, derfct_P1_3D, der2fct_P1_3D,
                                    refcoor_P1_3D,
                                    &geoLinearTria );

const GeometricMap geoBilinearHexa ( "Bilinear mapping on an hexaedra", HEXA,
                                     8, 3,
                                     fct_Q1_3D, derfct_Q1_3D, der2fct_Q1_3D,
                                     refcoor_Q1_3D,
                                     &geoBilinearQuad );

//======================================================================
//
//                            P0  (0D)
//
//======================================================================
/*
                           1
*/
Real fct1_P0_0D ( const GeoVector& )
{
    return 1.;
}
Real derfct1_P0_0D ( const GeoVector& )
{
    return 0.;
}
Real der2fct1_P0_0D ( const GeoVector& )
{
    return 0.;
}

//======================================================================
//
//                            P0  (1D)
//
//======================================================================
/*
                           --1--
*/
Real fct1_P0_1D ( const GeoVector& /*v*/ )
{
    return 1.;
}

Real derfct1_1_P0_1D ( const GeoVector& /*v*/ )
{
    return 0.;
}

Real der2fct1_P0_1D ( const GeoVector& /*v*/ )
{
    return 0.;
}

//======================================================================
//
//                            P1  (1D)
//
//======================================================================
/*
                           1-----2
*/
Real fct1_P1_1D ( const GeoVector& v )
{
    return 1 - v[0];
}
Real fct2_P1_1D ( const GeoVector& v )
{
    return v[0];
}

Real derfct1_1_P1_1D ( const GeoVector& )
{
    return -1;
}
Real derfct2_1_P1_1D ( const GeoVector& )
{
    return 1;
}

Real der2fct1_P1_1D ( const GeoVector& )
{
    return 0;
}

//======================================================================
//
//                            P2  (1D)
//
//======================================================================
/*
                           1--3--2
*/
Real fct1_P2_1D ( const GeoVector& v )
{
    return 2. * ( v[0] - 1. ) * ( v[0] - 0.5 );
}
Real fct3_P2_1D ( const GeoVector& v )
{
    return 4. * v[0] * ( 1. - v[0] );
}
Real fct2_P2_1D ( const GeoVector& v )
{
    return 2. * v[0] * ( v[0] - 0.5 );
}

Real derfct1_1_P2_1D ( const GeoVector& v )
{
    return 4. * v[0] - 3.;
}
Real derfct3_1_P2_1D ( const GeoVector& v )
{
    return -8. * v[0] + 4.;
}
Real derfct2_1_P2_1D ( const GeoVector& v )
{
    return 4. * v[0] - 1.;
}

Real der2fct1_11_P2_1D ( const GeoVector& )
{
    return 4;
}
Real der2fct3_11_P2_1D (  const GeoVector& )
{
    return -8;
}
Real der2fct2_11_P2_1D (  const GeoVector& )
{
    return 4;
}


//======================================================================
//
//                            P0  (2D)
//
//======================================================================
/*

                           |\
                           | \
                           | 1\
                            ---
*/
Real fct1_P0_2D ( const GeoVector& )
{
    return 1. ;
}
// First and Second derivatives are both equal (to 0).
Real derfct1_P0_2D ( const GeoVector& )
{
    return 0. ;
}
Real der2fct1_P0_2D ( const GeoVector& )
{
    return 0. ;
}

//======================================================================
//
//                            P1  (2D)
//
//======================================================================
/*
                           3
                           |\
                           | \
                           |  \
                           1---2
*/
Real fct1_P1_2D ( const GeoVector& v )
{
    return ( 1. - v[0] - v[1] );
}
Real fct2_P1_2D ( const GeoVector& v )
{
    return v[0] ;
}
Real fct3_P1_2D ( const GeoVector& v )
{
    return v[1] ;
}

Real derfct1_1_P1_2D ( const GeoVector& )
{
    return -1 ;
}
Real derfct1_2_P1_2D ( const GeoVector& )
{
    return -1 ;
}
Real derfct2_1_P1_2D ( const GeoVector& )
{
    return 1 ;
}
Real derfct2_2_P1_2D ( const GeoVector& )
{
    return 0 ;
}
Real derfct3_1_P1_2D ( const GeoVector& )
{
    return 0 ;
}
Real derfct3_2_P1_2D ( const GeoVector& )
{
    return 1 ;
}

// Second derivatives
Real der2fctx_xx_P1_2D ( const GeoVector& )
{
    return 0;
}



//======================================================================
//
//                            P1bubble  (2D)
//
//======================================================================
/*
                           3
                           |\
                           | \
                           |4.\
                           1---2
*/

Real fct1_P1bubble_2D ( const GeoVector& v )
{
    return ( 1. - v[0] - v[1] );
}
Real fct2_P1bubble_2D ( const GeoVector& v )
{
    return v[0] ;
}
Real fct3_P1bubble_2D ( const GeoVector& v )
{
    return v[1] ;
}

Real fct4_P1bubble_2D ( const GeoVector& v )
{
    return ( 1. - v[0] - v[1] ) * v[0] * v[1];
}

Real derfct1_1_P1bubble_2D ( const GeoVector& )
{
    return -1 ;
}
Real derfct1_2_P1bubble_2D ( const GeoVector& )
{
    return -1 ;
}
Real derfct2_1_P1bubble_2D ( const GeoVector& )
{
    return 1 ;
}
Real derfct2_2_P1bubble_2D ( const GeoVector& )
{
    return 0 ;
}
Real derfct3_1_P1bubble_2D ( const GeoVector& )
{
    return 0 ;
}
Real derfct3_2_P1bubble_2D ( const GeoVector& )
{
    return 1 ;
}

Real derfct4_1_P1bubble_2D ( const GeoVector& v )
{
    return ( 1 - 2 * v[0] - v[1] ) * v[1];
}
Real derfct4_2_P1bubble_2D ( const GeoVector& v )
{
    return ( 1 - v[0] - 2 * v[1] ) * v[0];
}

// Second derivatives
Real der2fctx_xx_P1bubble_2D ( const GeoVector& )
{
    return 0;
}

Real der2fct4_11_P1bubble_2D ( const GeoVector& v )
{
    return -2 * v[1];
}
Real der2fct4_12_P1bubble_2D ( const GeoVector& v )
{
    return 1 - 2 * v[0] - 2 * v[1];
}
Real der2fct4_21_P1bubble_2D ( const GeoVector& v )
{
    return 1 - 2 * v[0] - 2 * v[1];
}
Real der2fct4_22_P1bubble_2D ( const GeoVector& v )
{
    return -2 * v[0];
}



//======================================================================
//
//                            P2  (2D)
//
//======================================================================
/*
                           3
                           |\
                           6 5
                           |  \
                           1-4-2
*/
Real fct1_P2_2D ( const GeoVector& v )
{
    return ( 1 - v[0] - v[1] ) * ( 1 - v[0] - v[0] - v[1] - v[1] );
}
Real fct2_P2_2D ( const GeoVector& v )
{
    return -v[0] * ( 1 - v[0] - v[0] );
}
Real fct3_P2_2D ( const GeoVector& v )
{
    return -v[1] * ( 1 - v[1] - v[1] );
}
Real fct4_P2_2D ( const GeoVector& v )
{
    return 4 * v[0] * ( 1 - v[0] - v[1] );
}
Real fct5_P2_2D ( const GeoVector& v )
{
    return 4 * v[0] * v[1];
}
Real fct6_P2_2D ( const GeoVector& v )
{
    return 4 * v[1] * ( 1 - v[0] - v[1] );
}

Real derfct1_1_P2_2D ( const GeoVector& v )
{
    return 4 * ( v[0] + v[1] ) - 3;
}
Real derfct1_2_P2_2D ( const GeoVector& v )
{
    return 4 * ( v[0] + v[1] ) - 3;
}
Real derfct2_1_P2_2D ( const GeoVector& v )
{
    return 4 * v[0] - 1;
}
Real derfct2_2_P2_2D ( const GeoVector& )
{
    return 0;
}
Real derfct3_1_P2_2D ( const GeoVector& )
{
    return 0;
}
Real derfct3_2_P2_2D ( const GeoVector& v )
{
    return 4 * v[1] - 1;
}
Real derfct4_1_P2_2D ( const GeoVector& v )
{
    return 4 * ( 1 - v[0] - v[0] - v[1] );
}
Real derfct4_2_P2_2D ( const GeoVector& v )
{
    return -4 * v[0];
}
Real derfct5_1_P2_2D ( const GeoVector& v )
{
    return 4 * v[1];
}
Real derfct5_2_P2_2D ( const GeoVector& v )
{
    return 4 * v[0];
}
Real derfct6_1_P2_2D ( const GeoVector& v )
{
    return -4 * v[1];
}
Real derfct6_2_P2_2D ( const GeoVector& v )
{
    return 4 * ( 1 - v[0] - v[1] - v[1] );
}

Real der2fct1_11_P2_2D ( const GeoVector& )
{
    return 4;
}
Real der2fct1_12_P2_2D ( const GeoVector& )
{
    return 4;
}
Real der2fct1_21_P2_2D ( const GeoVector& )
{
    return 4;
}
Real der2fct1_22_P2_2D ( const GeoVector& )
{
    return 4;
}

Real der2fct2_11_P2_2D ( const GeoVector& )
{
    return 4;
}
Real der2fct2_12_P2_2D ( const GeoVector& )
{
    return 0;
}
Real der2fct2_21_P2_2D ( const GeoVector& )
{
    return 0;
}
Real der2fct2_22_P2_2D ( const GeoVector& )
{
    return 0;
}

Real der2fct3_11_P2_2D ( const GeoVector& )
{
    return 0;
}
Real der2fct3_12_P2_2D ( const GeoVector& )
{
    return 0;
}
Real der2fct3_21_P2_2D ( const GeoVector& )
{
    return 0;
}
Real der2fct3_22_P2_2D ( const GeoVector& )
{
    return 4;
}

Real der2fct4_11_P2_2D ( const GeoVector& )
{
    return -8;
}
Real der2fct4_12_P2_2D ( const GeoVector& )
{
    return -4;
}
Real der2fct4_21_P2_2D ( const GeoVector& )
{
    return -4;
}
Real der2fct4_22_P2_2D ( const GeoVector& )
{
    return 0;
}

Real der2fct5_11_P2_2D ( const GeoVector& )
{
    return 0;
}
Real der2fct5_12_P2_2D ( const GeoVector& )
{
    return 4;
}
Real der2fct5_21_P2_2D ( const GeoVector& )
{
    return 4;
}
Real der2fct5_22_P2_2D ( const GeoVector& )
{
    return 0;
}

Real der2fct6_11_P2_2D ( const GeoVector& )
{
    return 0;
}
Real der2fct6_12_P2_2D ( const GeoVector& )
{
    return -4;
}
Real der2fct6_21_P2_2D ( const GeoVector& )
{
    return -4;
}
Real der2fct6_22_P2_2D ( const GeoVector& )
{
    return -8;
}

//======================================================================
//
//                            RT0 Triangle  (2D)
//
//======================================================================
/*
                           3
                           |\
                           | \
                           |  \
                           1---2
*/
Real fct1_RT0_1_TRIA_2D ( const GeoVector& v )
{
    return v[0];
}
Real fct1_RT0_2_TRIA_2D ( const GeoVector& v )
{
    return v[1] - 1.;
}

Real fct2_RT0_1_TRIA_2D ( const GeoVector& v )
{
    return v[0];
}
Real fct2_RT0_2_TRIA_2D ( const GeoVector& v )
{
    return v[1];
}

Real fct3_RT0_1_TRIA_2D ( const GeoVector& v )
{
    return v[0] - 1.;
}
Real fct3_RT0_2_TRIA_2D ( const GeoVector& v )
{
    return v[1];
}

Real fct1_DIV_RT0_TRIA_2D ( const GeoVector& /*v*/ )
{
    return 2.;
}
Real fct2_DIV_RT0_TRIA_2D ( const GeoVector& /*v*/ )
{
    return 2.;
}
Real fct3_DIV_RT0_TRIA_2D ( const GeoVector& /*v*/ )
{
    return 2.;
}


//======================================================================
//
//                            Q0  (2D)
//
//======================================================================
/*
                            -------
                           |       |
                           |   1   |
                           |       |
                            -------
*/
Real fct1_Q0_2D ( const GeoVector& )
{
    return 1. ;
}
Real derfct1_Q0_2D ( const GeoVector& )
{
    return 0. ;
}
// The second derivative is equal to the first : both are equal to 0.
Real der2fct1_Q0_2D ( const GeoVector& )
{
    return 0. ;
}

//======================================================================
//
//                            Q1  (2D)
//
//======================================================================
/*
                           4-------3
                           |       |
                           |       |
                           |       |
                           1-------2
*/
Real fct1_Q1_2D ( const GeoVector& v )
{
    return ( 1. - v[0] ) * ( 1. - v[1] );
}
Real fct2_Q1_2D ( const GeoVector& v )
{
    return ( 1. - v[1] ) * v[0];
}
Real fct3_Q1_2D ( const GeoVector& v )
{
    return v[0] * v[1];
}
Real fct4_Q1_2D ( const GeoVector& v )
{
    return v[1] * ( 1. - v[0] );
}

Real derfct1_1_Q1_2D ( const GeoVector& v )
{
    return - ( 1. - v[1] );
}
Real derfct1_2_Q1_2D ( const GeoVector& v )
{
    return - ( 1. - v[0] );
}
Real derfct2_1_Q1_2D ( const GeoVector& v )
{
    return ( 1. - v[1] );
}
Real derfct2_2_Q1_2D ( const GeoVector& v )
{
    return -v[0];
}
Real derfct3_1_Q1_2D ( const GeoVector& v )
{
    return v[1];
}
Real derfct3_2_Q1_2D ( const GeoVector& v )
{
    return v[0];
}
Real derfct4_1_Q1_2D ( const GeoVector& v )
{
    return -v[1];
}
Real derfct4_2_Q1_2D ( const GeoVector& v )
{
    return ( 1. - v[0] );
}

// Second derivatives
Real der2fctx_xx_Q1_2D ( const GeoVector& )
{
    return 0;
}
//======================================================================
//
//                            Q2  (2D)
//
//======================================================================
/*
                           4---7---3
                           |       |
                           8   9   6
                           |       |
                           1---5---2
*/
Real fct1_Q2_2D ( const GeoVector& v )
{
    return 4. * ( 1 - v[0] ) * ( 0.5 - v[0] ) * ( 1 - v[1] ) * ( 0.5 - v[1] );
}
Real fct5_Q2_2D ( const GeoVector& v )
{
    return 8. * v[0] * ( 1 - v[0] ) * ( 1 - v[1] ) * ( 0.5 - v[1] );
}
Real fct2_Q2_2D ( const GeoVector& v )
{
    return 4. * v[0] * ( v[0] - 0.5 ) * ( 1 - v[1] ) * ( 0.5 - v[1] );
}
Real fct6_Q2_2D ( const GeoVector& v )
{
    return 8. * v[0] * ( v[0] - 0.5 ) * v[1] * ( 1 - v[1] );
}
Real fct3_Q2_2D ( const GeoVector& v )
{
    return 4. * v[0] * ( v[0] - 0.5 ) * v[1] * ( v[1] - 0.5 );
}
Real fct7_Q2_2D ( const GeoVector& v )
{
    return 8. * v[0] * ( 1 - v[0] ) * v[1] * ( v[1] - 0.5 );
}
Real fct4_Q2_2D ( const GeoVector& v )
{
    return 4. * ( 1 - v[0] ) * ( 0.5 - v[0] ) * v[1] * ( v[1] - 0.5 );
}
Real fct8_Q2_2D ( const GeoVector& v )
{
    return 8. * ( 0.5 - v[0] ) * ( 1 - v[0] ) * v[1] * ( 1 - v[1] );
}
Real fct9_Q2_2D ( const GeoVector& v )
{
    return 16. * v[0] * ( 1 - v[0] ) * v[1] * ( 1 - v[1] );
}

Real derfct1_1_Q2_2D ( const GeoVector& v )
{
    return ( 2. * v[1] - 1. ) * ( v[1] - 1. ) * ( 4. * v[0] - 3. );
}
Real derfct1_2_Q2_2D ( const GeoVector& v )
{
    return ( 2. * v[0] - 1. ) * ( v[0] - 1. ) * ( 4. * v[1] - 3. );
}
Real derfct5_1_Q2_2D ( const GeoVector& v )
{
    return -4. * ( 2. * v[1] - 1. ) * ( v[1] - 1. ) * ( 2. * v[0] - 1. );
}
Real derfct5_2_Q2_2D ( const GeoVector& v )
{
    return -4. * v[0] * ( v[0] - 1. ) * ( 4. * v[1] - 3. );
}
Real derfct2_1_Q2_2D ( const GeoVector& v )
{
    return ( 2. * v[1] - 1. ) * ( v[1] - 1. ) * ( 4. * v[0] - 1. );
}
Real derfct2_2_Q2_2D ( const GeoVector& v )
{
    return v[0] * ( 2. * v[0] - 1. ) * ( 4. * v[1] - 3. );
}
Real derfct6_1_Q2_2D ( const GeoVector& v )
{
    return -4. * v[1] * ( 4. * v[0] - 1. ) * ( v[1] - 1. );
}
Real derfct6_2_Q2_2D ( const GeoVector& v )
{
    return -4. * v[0] * ( 2. * v[0] - 1. ) * ( 2. * v[1] - 1. );
}
Real derfct3_1_Q2_2D ( const GeoVector& v )
{
    return v[1] * ( 4. * v[0] - 1. ) * ( 2. * v[1] - 1. );
}
Real derfct3_2_Q2_2D ( const GeoVector& v )
{
    return v[0] * ( 2. * v[0] - 1. ) * ( 4. * v[1] - 1. );
}
Real derfct7_1_Q2_2D ( const GeoVector& v )
{
    return -4. * v[1] * ( 2. * v[0] - 1. ) * ( 2. * v[1] - 1. );
}
Real derfct7_2_Q2_2D ( const GeoVector& v )
{
    return -4. * v[0] * ( v[0] - 1. ) * ( 4. * v[1] - 1. );
}
Real derfct4_1_Q2_2D ( const GeoVector& v )
{
    return v[1] * ( 4. * v[0] - 3. ) * ( 2. * v[1] - 1. );
}
Real derfct4_2_Q2_2D ( const GeoVector& v )
{
    return ( 2. * v[0] - 1. ) * ( v[0] - 1. ) * ( 4. * v[1] - 1. );
}
Real derfct8_1_Q2_2D ( const GeoVector& v )
{
    return -4. * v[1] * ( 4. * v[0] - 3. ) * ( v[1] - 1. );
}
Real derfct8_2_Q2_2D ( const GeoVector& v )
{
    return -4. * ( 2. * v[0] - 1. ) * ( v[0] - 1. ) * ( 2. * v[1] - 1. );
}
Real derfct9_1_Q2_2D ( const GeoVector& v )
{
    return 16. * v[1] * ( 2. * v[0] - 1. ) * ( v[1] - 1. );
}
Real derfct9_2_Q2_2D ( const GeoVector& v )
{
    return 16. * v[0] * ( v[0] - 1. ) * ( 2. * v[1] - 1. );
}

Real der2fct1_11_Q2_2D ( const GeoVector& v )
{
    return ( 2. * v[1] - 1. ) * ( v[1] - 1. ) * 4.;
}
Real der2fct1_12_Q2_2D ( const GeoVector& v )
{
    return ( 4. * v[1] - 3. ) * ( 4. * v[0] - 3. );
}
Real der2fct1_21_Q2_2D ( const GeoVector& v )
{
    return ( 4. * v[1] - 3. ) * ( 4. * v[0] - 3. );
}
Real der2fct1_22_Q2_2D ( const GeoVector& v )
{
    return ( 2. * v[0] - 1. ) * ( v[0] - 1. ) * 4.;
}

Real der2fct5_11_Q2_2D ( const GeoVector& v )
{
    return -8. * ( 2. * v[1] - 1. ) * ( v[1] - 1. );
}
Real der2fct5_12_Q2_2D ( const GeoVector& v )
{
    return -4. * ( 2. * v[0] - 1 ) * ( 4. * v[1] - 3 );
}
Real der2fct5_21_Q2_2D ( const GeoVector& v )
{
    return -4. * ( 2. * v[0] - 1 ) * ( 4. * v[1] - 3 );
    ;
}
Real der2fct5_22_Q2_2D ( const GeoVector& v )
{
    return -16. * v[0] * ( v[0] - 1. );
}

Real der2fct2_11_Q2_2D ( const GeoVector& v )
{
    return ( 2. * v[1] - 1. ) * ( v[1] - 1. ) * 4.;
}
Real der2fct2_12_Q2_2D ( const GeoVector& v )
{
    return ( 4. * v[0] - 1 ) * ( 4. * v[1] - 3. );
}
Real der2fct2_21_Q2_2D ( const GeoVector& v )
{
    return ( 4. * v[1] - 3. ) * ( 4. * v[0] - 1. );
}
Real der2fct2_22_Q2_2D ( const GeoVector& v )
{
    return v[0] * ( 2. * v[0] - 1. ) * 4.;
}

Real der2fct6_11_Q2_2D ( const GeoVector& v )
{
    return -16. * v[1] * ( v[1] - 1. );
}
Real der2fct6_12_Q2_2D ( const GeoVector& v )
{
    return -4. * ( 4. * v[0] - 1. ) * ( 2. * v[1] - 1. );
}
Real der2fct6_21_Q2_2D ( const GeoVector& v )
{
    return -4. * ( 4. * v[0] - 1. ) * ( 2. * v[1] - 1. );
}
Real der2fct6_22_Q2_2D ( const GeoVector& v )
{
    return -8. * v[0] * ( 2. * v[0] - 1. );
}

Real der2fct3_11_Q2_2D ( const GeoVector& v )
{
    return 4. * v[1] * ( 2. * v[1] - 1. );
}
Real der2fct3_12_Q2_2D ( const GeoVector& v )
{
    return ( 4. * v[0] - 1. ) * ( 4. * v[1] - 1. );
}
Real der2fct3_21_Q2_2D ( const GeoVector& v )
{
    return ( 4. * v[0] - 1. ) * ( 4. * v[1] - 1. );
}
Real der2fct3_22_Q2_2D ( const GeoVector& v )
{
    return 4. * v[0] * ( 2. * v[0] - 1. );
}

Real der2fct7_11_Q2_2D ( const GeoVector& v )
{
    return -8. * v[1] * ( 2. * v[1] - 1. );
}
Real der2fct7_12_Q2_2D ( const GeoVector& v )
{
    return -4. * ( 2. * v[0] - 1. ) * ( 4. * v[1] - 1. );
}
Real der2fct7_21_Q2_2D ( const GeoVector& v )
{
    return -4. * ( 2. * v[0] - 1. ) * ( 4. * v[1] - 1. );
}
Real der2fct7_22_Q2_2D ( const GeoVector& v )
{
    return -16. * v[0] * ( v[0] - 1. );
}

Real der2fct4_11_Q2_2D ( const GeoVector& v )
{
    return 4. * v[1] * ( 2. * v[1] - 1. );
}
Real der2fct4_12_Q2_2D ( const GeoVector& v )
{
    return ( 4. * v[0] - 3. ) * ( 4. * v[1] - 1. );
}
Real der2fct4_21_Q2_2D ( const GeoVector& v )
{
    return ( 4. * v[0] - 3. ) * ( 4. * v[1] - 1. );
}
Real der2fct4_22_Q2_2D ( const GeoVector& v )
{
    return 4. * ( 2. * v[0] - 1. ) * ( v[0] - 1. );
}

Real der2fct8_11_Q2_2D ( const GeoVector& v )
{
    return -16. * v[1] * ( v[1] - 1. );
}
Real der2fct8_12_Q2_2D ( const GeoVector& v )
{
    return -4. * ( 4. * v[0] - 3. ) * ( 2. * v[1] - 1. );
}
Real der2fct8_21_Q2_2D ( const GeoVector& v )
{
    return -4. * ( 4. * v[0] - 3. ) * ( 2. * v[1] - 1. );
}
Real der2fct8_22_Q2_2D ( const GeoVector& v )
{
    return -8. * ( 2. * v[0] - 1. ) * ( v[0] - 1. );
}

Real der2fct9_11_Q2_2D ( const GeoVector& v )
{
    return 32. * v[1] * ( v[1] - 1. );
}
Real der2fct9_12_Q2_2D ( const GeoVector& v )
{
    return 16. * ( 2. * v[0] - 1. ) * ( 2. * v[1] - 1. );
}
Real der2fct9_21_Q2_2D ( const GeoVector& v )
{
    return 16. * ( 2. * v[0] - 1. ) * ( 2. * v[1] - 1. );
}
Real der2fct9_22_Q2_2D ( const GeoVector& v )
{
    return 32. * v[0] * ( v[0] - 1. );
}

//======================================================================
//
//                            P0  (3D)
//
//======================================================================
/*
                4
               / .
              /  \.3
             /  . \\
            / .    \\
           /.       \!
         1 ----------2
*/
Real fct1_P0_3D ( const GeoVector& )
{
    return 1.;
}

Real derfct1_P0_3D ( const GeoVector& )
{
    return 0.;
}

// Second derivatives
Real der2fct1_P0_3D ( const GeoVector& )
{
    return 0;
}


//======================================================================
//
//                            P1  (3D)
//
//======================================================================
/*
                4
               / .
              /  \.3
             /  . \\
            / .    \\
           /.       \!
         1 ----------2
*/
Real fct1_P1_3D ( const GeoVector& v )
{
    return 1 - v[0] - v[1] - v[2];
}
Real fct2_P1_3D ( const GeoVector& v )
{
    return v[0];
}
Real fct3_P1_3D ( const GeoVector& v )
{
    return v[1];
}
Real fct4_P1_3D ( const GeoVector& v )
{
    return v[2];
}

Real derfct1_1_P1_3D ( const GeoVector& )
{
    return -1;
}
Real derfct1_2_P1_3D ( const GeoVector& )
{
    return -1;
}
Real derfct1_3_P1_3D ( const GeoVector& )
{
    return -1;
}
Real derfct2_1_P1_3D ( const GeoVector& )
{
    return 1;
}
Real derfct2_2_P1_3D ( const GeoVector& )
{
    return 0;
}
Real derfct2_3_P1_3D ( const GeoVector& )
{
    return 0;
}
Real derfct3_1_P1_3D ( const GeoVector& )
{
    return 0;
}
Real derfct3_2_P1_3D ( const GeoVector& )
{
    return 1;
}
Real derfct3_3_P1_3D ( const GeoVector& )
{
    return 0;
}
Real derfct4_1_P1_3D ( const GeoVector& )
{
    return 0;
}
Real derfct4_2_P1_3D ( const GeoVector& )
{
    return 0;
}
Real derfct4_3_P1_3D ( const GeoVector& )
{
    return 1;
}

// Second derivatives
Real der2fctx_xx_P1_3D ( const GeoVector& )
{
    return 0;
}
//======================================================================
//======================================================================
//
//                            P1bubble  (3D)
//
//======================================================================
/*
                4
               / .
              /  \.3
             /  . \\
            / . .5 \\
           /.       \!
         1 ----------2
*/
Real fct1_P1bubble_3D ( const GeoVector& v )
{
    return 1 - v[0] - v[1] - v[2];
}
Real fct2_P1bubble_3D ( const GeoVector& v )
{
    return v[0];
}
Real fct3_P1bubble_3D ( const GeoVector& v )
{
    return v[1];
}
Real fct4_P1bubble_3D ( const GeoVector& v )
{
    return v[2];
}
Real fct5_P1bubble_3D ( const GeoVector& v )
{
    return ( 1 - v[0] - v[1] - v[2] ) * v[0] * v[1] * v[2];
}

Real derfct1_1_P1bubble_3D ( const GeoVector& )
{
    return -1;
}
Real derfct1_2_P1bubble_3D ( const GeoVector& )
{
    return -1;
}
Real derfct1_3_P1bubble_3D ( const GeoVector& )
{
    return -1;
}
Real derfct2_1_P1bubble_3D ( const GeoVector& )
{
    return 1;
}
Real derfct2_2_P1bubble_3D ( const GeoVector& )
{
    return 0;
}
Real derfct2_3_P1bubble_3D ( const GeoVector& )
{
    return 0;
}
Real derfct3_1_P1bubble_3D ( const GeoVector& )
{
    return 0;
}
Real derfct3_2_P1bubble_3D ( const GeoVector& )
{
    return 1;
}
Real derfct3_3_P1bubble_3D ( const GeoVector& )
{
    return 0;
}
Real derfct4_1_P1bubble_3D ( const GeoVector& )
{
    return 0;
}
Real derfct4_2_P1bubble_3D ( const GeoVector& )
{
    return 0;
}
Real derfct4_3_P1bubble_3D ( const GeoVector& )
{
    return 1;
}
Real derfct5_1_P1bubble_3D ( const GeoVector& v )
{
    return ( 1 - 2 * v[0] - v[1] - v[2] ) * v[1] * v[2];
}
Real derfct5_2_P1bubble_3D ( const GeoVector& v )
{
    return ( 1 - v[0] - 2 * v[1] - v[2] ) * v[0] * v[2];
}
Real derfct5_3_P1bubble_3D ( const GeoVector& v )
{
    return ( 1 - v[0] - v[1] - 2 * v[2] ) * v[0] * v[1];
}

// Second derivatives
Real der2fctx_xx_P1bubble_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct5_11_P1bubble_3D ( const GeoVector& v )
{
    return -2 * v[1] * v[2];
}
Real der2fct5_12_P1bubble_3D ( const GeoVector& v )
{
    return ( 1 - 2 * v[0] - 2 * v[1] - v[2] ) * v[2];
}
Real der2fct5_13_P1bubble_3D ( const GeoVector& v )
{
    return ( 1 - 2 * v[0] - v[1] - 2 * v[2] ) * v[1];
}
Real der2fct5_21_P1bubble_3D ( const GeoVector& v )
{
    return ( 1 - 2 * v[0] - 2 * v[1] - v[2] ) * v[2];
}
Real der2fct5_22_P1bubble_3D ( const GeoVector& v )
{
    return -2 * v[0] * v[2];
}
Real der2fct5_23_P1bubble_3D ( const GeoVector& v )
{
    return ( 1 - v[0] - 2 * v[1] - 2 * v[2] ) * v[0];
}
Real der2fct5_31_P1bubble_3D ( const GeoVector& v )
{
    return ( 1 - 2 * v[0] - v[1] - 2 * v[2] ) * v[1];
}
Real der2fct5_32_P1bubble_3D ( const GeoVector& v )
{
    return ( 1 - v[0] - 2 * v[1] - 2 * v[2] ) * v[0];
}
Real der2fct5_33_P1bubble_3D ( const GeoVector& v )
{
    return -2 * v[0] * v[1];
}

//======================================================================
//
//                            P2  (3D)
//
//======================================================================
/*
                4
               / .10
              /  \.3
             8  . 9\
            / 7    \6
           /.       \!
         1 -----5----2
*/
Real fct1_P2_3D ( const GeoVector& v )
{
    return - ( 1 - v[0] - v[1] - v[2] ) * ( 1 - 2 * ( 1 - v[0] - v[1] - v[2] ) );
}
Real fct2_P2_3D ( const GeoVector& v )
{
    return -v[0] * ( 1 - 2 * v[0] );
}
Real fct3_P2_3D ( const GeoVector& v )
{
    return -v[1] * ( 1 - 2 * v[1] );
}
Real fct4_P2_3D ( const GeoVector& v )
{
    return -v[2] * ( 1 - 2 * v[2] );
}
Real fct5_P2_3D ( const GeoVector& v )
{
    return 4 * v[0] * ( 1 - v[0] - v[1] - v[2] );
}
Real fct6_P2_3D ( const GeoVector& v )
{
    return 4 * v[0] * v[1];
}
Real fct7_P2_3D ( const GeoVector& v )
{
    return 4 * v[1] * ( 1 - v[0] - v[1] - v[2] );
}
Real fct8_P2_3D ( const GeoVector& v )
{
    return 4 * v[2] * ( 1 - v[0] - v[1] - v[2] );
}
Real fct9_P2_3D ( const GeoVector& v )
{
    return 4 * v[0] * v[2];
}
Real fct10_P2_3D ( const GeoVector& v )
{
    return 4 * v[1] * v[2];
}


Real derfct1_1_P2_3D ( const GeoVector& v )
{
    return -3 + 4 * v[0] + 4 * v[1] + 4 * v[2];
}
Real derfct1_2_P2_3D ( const GeoVector& v )
{
    return -3 + 4 * v[0] + 4 * v[1] + 4 * v[2];
}
Real derfct1_3_P2_3D ( const GeoVector& v )
{
    return -3 + 4 * v[0] + 4 * v[1] + 4 * v[2];
}

Real derfct2_1_P2_3D ( const GeoVector& v )
{
    return -1 + 4 * v[0];
}
Real derfct2_2_P2_3D ( const GeoVector& )
{
    return 0.;
}
Real derfct2_3_P2_3D ( const GeoVector& )
{
    return 0.;
}

Real derfct3_1_P2_3D ( const GeoVector& )
{
    return 0.;
}
Real derfct3_2_P2_3D ( const GeoVector& v )
{
    return -1 + 4 * v[1];
}
Real derfct3_3_P2_3D ( const GeoVector& )
{
    return 0.;
}

Real derfct4_1_P2_3D ( const GeoVector& )
{
    return 0.;
}
Real derfct4_2_P2_3D ( const GeoVector& )
{
    return 0.;
}
Real derfct4_3_P2_3D ( const GeoVector& v )
{
    return -1 + 4 * v[2];
}

Real derfct5_1_P2_3D ( const GeoVector& v )
{
    return 4 - 8 * v[0] - 4 * v[1] - 4 * v[2];
}
Real derfct5_2_P2_3D ( const GeoVector& v )
{
    return -4 * v[0];
}
Real derfct5_3_P2_3D ( const GeoVector& v )
{
    return -4 * v[0];
}

Real derfct6_1_P2_3D ( const GeoVector& v )
{
    return 4 * v[1];
}
Real derfct6_2_P2_3D ( const GeoVector& v )
{
    return 4 * v[0];
}
Real derfct6_3_P2_3D ( const GeoVector& )
{
    return 0.;
}

Real derfct7_1_P2_3D ( const GeoVector& v )
{
    return -4 * v[1];
}
Real derfct7_2_P2_3D ( const GeoVector& v )
{
    return 4 - 4 * v[0] - 8 * v[1] - 4 * v[2];
}
Real derfct7_3_P2_3D ( const GeoVector& v )
{
    return -4 * v[1];
}

Real derfct8_1_P2_3D ( const GeoVector& v )
{
    return -4 * v[2];
}
Real derfct8_2_P2_3D ( const GeoVector& v )
{
    return -4 * v[2];
}
Real derfct8_3_P2_3D ( const GeoVector& v )
{
    return 4 - 4 * v[0] - 4 * v[1] - 8 * v[2];
}

Real derfct9_1_P2_3D ( const GeoVector& v )
{
    return 4 * v[2];
}
Real derfct9_2_P2_3D ( const GeoVector& )
{
    return 0.;
}
Real derfct9_3_P2_3D ( const GeoVector& v )
{
    return 4 * v[0];
}

Real derfct10_1_P2_3D ( const GeoVector& )
{
    return 0.;
}
Real derfct10_2_P2_3D ( const GeoVector& v )
{
    return 4 * v[2];
}
Real derfct10_3_P2_3D ( const GeoVector& v )
{
    return 4 * v[1];
}


Real der2fct1_11_P2_3D ( const GeoVector& )
{
    return 4;
}
Real der2fct1_12_P2_3D ( const GeoVector& )
{
    return 4;
}
Real der2fct1_13_P2_3D ( const GeoVector& )
{
    return 4;
}
Real der2fct1_21_P2_3D ( const GeoVector& )
{
    return 4;
}
Real der2fct1_22_P2_3D ( const GeoVector& )
{
    return 4;
}
Real der2fct1_23_P2_3D ( const GeoVector& )
{
    return 4;
}
Real der2fct1_31_P2_3D ( const GeoVector& )
{
    return 4;
}
Real der2fct1_32_P2_3D ( const GeoVector& )
{
    return 4;
}
Real der2fct1_33_P2_3D ( const GeoVector& )
{
    return 4;
}

Real der2fct2_11_P2_3D ( const GeoVector& )
{
    return 4;
}
Real der2fct2_12_P2_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct2_13_P2_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct2_21_P2_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct2_22_P2_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct2_23_P2_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct2_31_P2_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct2_32_P2_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct2_33_P2_3D ( const GeoVector& )
{
    return 0;
}

Real der2fct3_11_P2_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct3_12_P2_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct3_13_P2_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct3_21_P2_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct3_22_P2_3D ( const GeoVector& )
{
    return 4;
}
Real der2fct3_23_P2_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct3_31_P2_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct3_32_P2_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct3_33_P2_3D ( const GeoVector& )
{
    return 0;
}

Real der2fct4_11_P2_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct4_12_P2_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct4_13_P2_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct4_21_P2_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct4_22_P2_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct4_23_P2_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct4_31_P2_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct4_32_P2_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct4_33_P2_3D ( const GeoVector& )
{
    return 4;
}

Real der2fct5_11_P2_3D ( const GeoVector& )
{
    return -8;
}
Real der2fct5_12_P2_3D ( const GeoVector& )
{
    return -4;
}
Real der2fct5_13_P2_3D ( const GeoVector& )
{
    return -4;
}
Real der2fct5_21_P2_3D ( const GeoVector& )
{
    return -4;
}
Real der2fct5_22_P2_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct5_23_P2_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct5_31_P2_3D ( const GeoVector& )
{
    return -4;
}
Real der2fct5_32_P2_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct5_33_P2_3D ( const GeoVector& )
{
    return 0;
}

Real der2fct6_11_P2_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct6_12_P2_3D ( const GeoVector& )
{
    return 4;
}
Real der2fct6_13_P2_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct6_21_P2_3D ( const GeoVector& )
{
    return 4;
}
Real der2fct6_22_P2_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct6_23_P2_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct6_31_P2_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct6_32_P2_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct6_33_P2_3D ( const GeoVector& )
{
    return 0;
}

Real der2fct7_11_P2_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct7_12_P2_3D ( const GeoVector& )
{
    return -4;
}
Real der2fct7_13_P2_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct7_21_P2_3D ( const GeoVector& )
{
    return -4;
}
Real der2fct7_22_P2_3D ( const GeoVector& )
{
    return -8;
}
Real der2fct7_23_P2_3D ( const GeoVector& )
{
    return -4;
}
Real der2fct7_31_P2_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct7_32_P2_3D ( const GeoVector& )
{
    return -4;
}
Real der2fct7_33_P2_3D ( const GeoVector& )
{
    return 0;
}

Real der2fct8_11_P2_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct8_12_P2_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct8_13_P2_3D ( const GeoVector& )
{
    return -4;
}
Real der2fct8_21_P2_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct8_22_P2_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct8_23_P2_3D ( const GeoVector& )
{
    return -4;
}
Real der2fct8_31_P2_3D ( const GeoVector& )
{
    return -4;
}
Real der2fct8_32_P2_3D ( const GeoVector& )
{
    return -4;
}
Real der2fct8_33_P2_3D ( const GeoVector& )
{
    return -8;
}

Real der2fct9_11_P2_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct9_12_P2_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct9_13_P2_3D ( const GeoVector& )
{
    return 4;
}
Real der2fct9_21_P2_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct9_22_P2_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct9_23_P2_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct9_31_P2_3D ( const GeoVector& )
{
    return 4;
}
Real der2fct9_32_P2_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct9_33_P2_3D ( const GeoVector& )
{
    return 0;
}

Real der2fct10_11_P2_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct10_12_P2_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct10_13_P2_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct10_21_P2_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct10_22_P2_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct10_23_P2_3D ( const GeoVector& )
{
    return 4;
}
Real der2fct10_31_P2_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct10_32_P2_3D ( const GeoVector& )
{
    return 4;
}
Real der2fct10_33_P2_3D ( const GeoVector& )
{
    return 0;
}
//======================================================================
//
//                            P2tilde  (3D)
// NAVIER-STOKES P2 Basis Oriented to the mass lumping
//======================================================================
/*
                4
               / .10
              /  \.3
             8  . 9\
            / 7 .11 \6
           /.       \!
         1 -----5----2
*/
Real fct1_P2tilde_3D ( const GeoVector& v )
{
    return - ( 1 - v[0] - v[1] - v[2] ) * ( 1 - 2 * ( 1 - v[0] - v[1] - v[2] ) ) + 32 * v[0] * v[1] * v[2] * ( 1 - v[0] - v[1] - v[2] );
}
Real fct2_P2tilde_3D ( const GeoVector& v )
{
    return -v[0] * ( 1 - 2 * v[0] ) + 32 * v[0] * v[1] * v[2] * ( 1 - v[0] - v[1] - v[2] );
}
Real fct3_P2tilde_3D ( const GeoVector& v )
{
    return -v[1] * ( 1 - 2 * v[1] ) + 32 * v[0] * v[1] * v[2] * ( 1 - v[0] - v[1] - v[2] );
}
Real fct4_P2tilde_3D ( const GeoVector& v )
{
    return -v[2] * ( 1 - 2 * v[2] ) + 32 * v[0] * v[1] * v[2] * ( 1 - v[0] - v[1] - v[2] );
}

Real fct5_P2tilde_3D ( const GeoVector& v )
{
    return 4 * v[0] * ( 1 - v[0] - v[1] - v[2] ) - 64 * v[0] * v[1] * v[2] * ( 1 - v[0] - v[1] - v[2] );
}
Real fct6_P2tilde_3D ( const GeoVector& v )
{
    return 4 * v[0] * v[1] - 64 * v[0] * v[1] * v[2] * ( 1 - v[0] - v[1] - v[2] );
}
Real fct7_P2tilde_3D ( const GeoVector& v )
{
    return 4 * v[1] * ( 1 - v[0] - v[1] - v[2] ) - 64 * v[0] * v[1] * v[2] * ( 1 - v[0] - v[1] - v[2] );
}
Real fct8_P2tilde_3D ( const GeoVector& v )
{
    return 4 * v[2] * ( 1 - v[0] - v[1] - v[2] ) - 64 * v[0] * v[1] * v[2] * ( 1 - v[0] - v[1] - v[2] );
}
Real fct9_P2tilde_3D ( const GeoVector& v )
{
    return 4 * v[0] * v[2] - 64 * v[0] * v[1] * v[2] * ( 1 - v[0] - v[1] - v[2] );
}
Real fct10_P2tilde_3D ( const GeoVector& v )
{
    return 4 * v[1] * v[2] - 64 * v[0] * v[1] * v[2] * ( 1 - v[0] - v[1] - v[2] );
}

Real fct11_P2tilde_3D ( const GeoVector& v )
{
    return 256 * v[0] * v[1] * v[2] * ( 1 - v[0] - v[1] - v[2] );
}


Real derfct1_1_P2tilde_3D ( const GeoVector& v )
{
    return -3 + 4 * v[0] + 4 * v[1] + 4 * v[2] + 32 * v[1] * v[2] * ( 1 - 2 * v[0] - v[1] - v[2] );
}
Real derfct1_2_P2tilde_3D ( const GeoVector& v )
{
    return -3 + 4 * v[0] + 4 * v[1] + 4 * v[2] + 32 * v[0] * v[2] * ( 1 - v[0] - 2 * v[1] - v[2] );
}
Real derfct1_3_P2tilde_3D ( const GeoVector& v )
{
    return -3 + 4 * v[0] + 4 * v[1] + 4 * v[2] + 32 * v[0] * v[1] * ( 1 - v[0] - v[1] - 2 * v[2] );
}

Real derfct2_1_P2tilde_3D ( const GeoVector& v )
{
    return -1 + 4 * v[0] + 32 * v[1] * v[2] * ( 1 - 2 * v[0] - v[1] - v[2] );
}
Real derfct2_2_P2tilde_3D ( const GeoVector& v )
{
    return 32 * v[0] * v[2] * ( 1 - v[0] - 2 * v[1] - v[2] );
}
Real derfct2_3_P2tilde_3D ( const GeoVector& v )
{
    return 32 * v[0] * v[1] * ( 1 - v[0] - v[1] - 2 * v[2] );
}

Real derfct3_1_P2tilde_3D ( const GeoVector& v )
{
    return 32 * v[1] * v[2] * ( 1 - 2 * v[0] - v[1] - v[2] );
}
Real derfct3_2_P2tilde_3D ( const GeoVector& v )
{
    return -1 + 4 * v[1] + 32 * v[0] * v[2] * ( 1 - v[0] - 2 * v[1] - v[2] );
}
Real derfct3_3_P2tilde_3D ( const GeoVector& v )
{
    return 32 * v[0] * v[1] * ( 1 - v[0] - v[1] - 2 * v[2] );
}

Real derfct4_1_P2tilde_3D ( const GeoVector& v )
{
    return 32 * v[1] * v[2] * ( 1 - 2 * v[0] - v[1] - v[2] );
}
Real derfct4_2_P2tilde_3D ( const GeoVector& v )
{
    return 32 * v[0] * v[2] * ( 1 - v[0] - 2 * v[1] - v[2] );
}
Real derfct4_3_P2tilde_3D ( const GeoVector& v )
{
    return -1 + 4 * v[2] + 32 * v[0] * v[1] * ( 1 - v[0] - v[1] - 2 * v[2] );
}

Real derfct5_1_P2tilde_3D ( const GeoVector& v )
{
    return 4 - 8 * v[0] - 4 * v[1] - 4 * v[2] - 64 * v[1] * v[2] * ( 1 - 2 * v[0] - v[1] - v[2] );
}
Real derfct5_2_P2tilde_3D ( const GeoVector& v )
{
    return -4 * v[0] - 64 * v[0] * v[2] * ( 1 - v[0] - 2 * v[1] - v[2] );
}
Real derfct5_3_P2tilde_3D ( const GeoVector& v )
{
    return -4 * v[0] - 64 * v[0] * v[1] * ( 1 - v[0] - v[1] - 2 * v[2] );
}

Real derfct6_1_P2tilde_3D ( const GeoVector& v )
{
    return 4 * v[1] - 64 * v[1] * v[2] * ( 1 - 2 * v[0] - v[1] - v[2] );
}
Real derfct6_2_P2tilde_3D ( const GeoVector& v )
{
    return 4 * v[0] - 64 * v[0] * v[2] * ( 1 - v[0] - 2 * v[1] - v[2] );
}
Real derfct6_3_P2tilde_3D ( const GeoVector& v )
{
    return - 64 * v[0] * v[1] * ( 1 - v[0] - v[1] - 2 * v[2] );
}

Real derfct7_1_P2tilde_3D ( const GeoVector& v )
{
    return -4 * v[1] - 64 * v[1] * v[2] * ( 1 - 2 * v[0] - v[1] - v[2] );
}
Real derfct7_2_P2tilde_3D ( const GeoVector& v )
{
    return 4 - 4 * v[0] - 8 * v[1] - 4 * v[2] - 64 * v[0] * v[2] * ( 1 - v[0] - 2 * v[1] - v[2] );
}
Real derfct7_3_P2tilde_3D ( const GeoVector& v )
{
    return -4 * v[1] - 64 * v[0] * v[1] * ( 1 - v[0] - v[1] - 2 * v[2] );
}

Real derfct8_1_P2tilde_3D ( const GeoVector& v )
{
    return -4 * v[2] - 64 * v[1] * v[2] * ( 1 - 2 * v[0] - v[1] - v[2] );
}
Real derfct8_2_P2tilde_3D ( const GeoVector& v )
{
    return -4 * v[2] - 64 * v[0] * v[2] * ( 1 - v[0] - 2 * v[1] - v[2] );
}
Real derfct8_3_P2tilde_3D ( const GeoVector& v )
{
    return 4 - 4 * v[0] - 4 * v[1] - 8 * v[2] - 64 * v[0] * v[1] * ( 1 - v[0] - v[1] - 2 * v[2] );
}

Real derfct9_1_P2tilde_3D ( const GeoVector& v )
{
    return 4 * v[2] - 64 * v[1] * v[2] * ( 1 - 2 * v[0] - v[1] - v[2] );
}
Real derfct9_2_P2tilde_3D ( const GeoVector& v )
{
    return - 64 * v[0] * v[2] * ( 1 - v[0] - 2 * v[1] - v[2] );
}
Real derfct9_3_P2tilde_3D ( const GeoVector& v )
{
    return 4 * v[0] - 64 * v[0] * v[1] * ( 1 - v[0] - v[1] - 2 * v[2] );
}

Real derfct10_1_P2tilde_3D ( const GeoVector& v )
{
    return - 64 * v[1] * v[2] * ( 1 - 2 * v[0] - v[1] - v[2] );
}
Real derfct10_2_P2tilde_3D ( const GeoVector& v )
{
    return 4 * v[2] - 64 * v[0] * v[2] * ( 1 - v[0] - 2 * v[1] - v[2] );
}
Real derfct10_3_P2tilde_3D ( const GeoVector& v )
{
    return 4 * v[1] - 64 * v[0] * v[1] * ( 1 - v[0] - v[1] - 2 * v[2] );
}

Real derfct11_1_P2tilde_3D ( const GeoVector& v )
{
    return 256 * v[1] * v[2] * ( 1 - 2 * v[0] - v[1] - v[2] );
}
Real derfct11_2_P2tilde_3D ( const GeoVector& v )
{
    return 256 * v[0] * v[2] * ( 1 - v[0] - 2 * v[1] - v[2] );
}
Real derfct11_3_P2tilde_3D ( const GeoVector& v )
{
    return 256 * v[0] * v[1] * ( 1 - v[0] - v[1] - 2 * v[2] );
}

Real der2fct1_11_P2tilde_3D ( const GeoVector& v )
{
    return 4 - 64 * v[0] * v[1] * v[2];
}
Real der2fct1_12_P2tilde_3D ( const GeoVector& v )
{
    return 4 + 32 * v[2] * ( 1 - 2 * v[0] - 2 * v[1] - v[2] );
}
Real der2fct1_13_P2tilde_3D ( const GeoVector& v )
{
    return 4 + 32 * v[1] * ( 1 - 2 * v[0] - v[1] - 2 * v[2] );
}
Real der2fct1_21_P2tilde_3D ( const GeoVector& v )
{
    return 4 + 32 * v[2] * ( 1 - 2 * v[0] - 2 * v[1] - v[2] );
}
Real der2fct1_22_P2tilde_3D ( const GeoVector& v )
{
    return 4 - 64 * v[0] * v[1] * v[2];
}
Real der2fct1_23_P2tilde_3D ( const GeoVector& v )
{
    return 4 + 32 * v[0] * ( 1 - v[0] - 2 * v[1] - 2 * v[2] );
}
Real der2fct1_31_P2tilde_3D ( const GeoVector& v )
{
    return 4 + 32 * v[1] * ( 1 - 2 * v[0] - v[1] - 2 * v[2] );
}
Real der2fct1_32_P2tilde_3D ( const GeoVector& v )
{
    return 4 + 32 * v[0] * ( 1 - v[0] - 2 * v[1] - 2 * v[2] );
    ;
}
Real der2fct1_33_P2tilde_3D ( const GeoVector& v )
{
    return 4 - 64 * v[0] * v[1] * v[2];
}

Real der2fct2_11_P2tilde_3D ( const GeoVector& v )
{
    return 4 - 64 * v[0] * v[1] * v[2];
}
Real der2fct2_12_P2tilde_3D ( const GeoVector& v )
{
    return 32 * v[2] * ( 1 - 2 * v[0] - 2 * v[1] - v[2] );
}
Real der2fct2_13_P2tilde_3D ( const GeoVector& v )
{
    return 32 * v[1] * ( 1 - 2 * v[0] - v[1] - 2 * v[2] );
}
Real der2fct2_21_P2tilde_3D ( const GeoVector& v )
{
    return 32 * v[2] * ( 1 - 2 * v[0] - 2 * v[1] - v[2] );
}
Real der2fct2_22_P2tilde_3D ( const GeoVector& v )
{
    return - 64 * v[0] * v[1] * v[2];
}
Real der2fct2_23_P2tilde_3D ( const GeoVector& v )
{
    return 32 * v[0] * ( 1 - v[0] - 2 * v[1] - 2 * v[2] );
}
Real der2fct2_31_P2tilde_3D ( const GeoVector& v )
{
    return 32 * v[1] * ( 1 - 2 * v[0] - v[1] - 2 * v[2] );
}
Real der2fct2_32_P2tilde_3D ( const GeoVector& v )
{
    return 32 * v[0] * ( 1 - v[0] - 2 * v[1] - 2 * v[2] );
}
Real der2fct2_33_P2tilde_3D ( const GeoVector& v )
{
    return - 64 * v[0] * v[1] * v[2];
}

Real der2fct3_11_P2tilde_3D ( const GeoVector& v )
{
    return - 64 * v[0] * v[1] * v[2];
}
Real der2fct3_12_P2tilde_3D ( const GeoVector& v )
{
    return 32 * v[2] * ( 1 - 2 * v[0] - 2 * v[1] - v[2] );
}
Real der2fct3_13_P2tilde_3D ( const GeoVector& v )
{
    return 32 * v[1] * ( 1 - 2 * v[0] - v[1] - 2 * v[2] );
}
Real der2fct3_21_P2tilde_3D ( const GeoVector& v )
{
    return 32 * v[2] * ( 1 - 2 * v[0] - 2 * v[1] - v[2] );
}
Real der2fct3_22_P2tilde_3D ( const GeoVector& v )
{
    return 4 - 64 * v[0] * v[1] * v[2];
}
Real der2fct3_23_P2tilde_3D ( const GeoVector& v )
{
    return 32 * v[0] * ( 1 - v[0] - 2 * v[1] - 2 * v[2] );
}
Real der2fct3_31_P2tilde_3D ( const GeoVector& v )
{
    return 32 * v[1] * ( 1 - 2 * v[0] - v[1] - 2 * v[2] );
}
Real der2fct3_32_P2tilde_3D ( const GeoVector& v )
{
    return 32 * v[0] * ( 1 - v[0] - 2 * v[1] - 2 * v[2] );
}
Real der2fct3_33_P2tilde_3D ( const GeoVector& v )
{
    return - 64 * v[0] * v[1] * v[2];
}

Real der2fct4_11_P2tilde_3D ( const GeoVector& v )
{
    return -64 * v[0] * v[1] * v[2];
}
Real der2fct4_12_P2tilde_3D ( const GeoVector& v )
{
    return 32 * v[2] * ( 1 - 2 * v[0] - 2 * v[1] - v[2] );
}
Real der2fct4_13_P2tilde_3D ( const GeoVector& v )
{
    return 32 * v[1] * ( 1 - 2 * v[0] - v[1] - 2 * v[2] );
}
Real der2fct4_21_P2tilde_3D ( const GeoVector& v )
{
    return 32 * v[2] * ( 1 - 2 * v[0] - 2 * v[1] - v[2] );
}
Real der2fct4_22_P2tilde_3D ( const GeoVector& v )
{
    return - 64 * v[0] * v[1] * v[2];
}
Real der2fct4_23_P2tilde_3D ( const GeoVector& v )
{
    return 32 * v[0] * ( 1 - v[0] - 2 * v[1] - 2 * v[2] );
}
Real der2fct4_31_P2tilde_3D ( const GeoVector& v )
{
    return 32 * v[1] * ( 1 - 2 * v[0] - v[1] - 2 * v[2] );
}
Real der2fct4_32_P2tilde_3D ( const GeoVector& v )
{
    return 32 * v[0] * ( 1 - v[0] - 2 * v[1] - 2 * v[2] );
}
Real der2fct4_33_P2tilde_3D ( const GeoVector& v )
{
    return 4 - 64 * v[0] * v[1] * v[2];
}

Real der2fct5_11_P2tilde_3D ( const GeoVector& v )
{
    return -8 - 128 * v[0] * v[1] * v[2];
}
Real der2fct5_12_P2tilde_3D ( const GeoVector& v )
{
    return -4 + 64 * v[2] * ( 1 - 2 * v[0] - 2 * v[1] - v[2] );
}
Real der2fct5_13_P2tilde_3D ( const GeoVector& v )
{
    return -4 + 64 * v[1] * ( 1 - 2 * v[0] - v[1] - 2 * v[2] );
}
Real der2fct5_21_P2tilde_3D ( const GeoVector& v )
{
    return -4 + 64 * v[2] * ( 1 - 2 * v[0] - 2 * v[1] - v[2] );
}
Real der2fct5_22_P2tilde_3D ( const GeoVector& v )
{
    return - 128 * v[0] * v[1] * v[2];
}
Real der2fct5_23_P2tilde_3D ( const GeoVector& v )
{
    return 64 * v[0] * ( 1 - v[0] - 2 * v[1] - 2 * v[2] );
}
Real der2fct5_31_P2tilde_3D ( const GeoVector& v )
{
    return -4 + 64 * v[1] * ( 1 - 2 * v[0] - v[1] - 2 * v[2] );
}
Real der2fct5_32_P2tilde_3D ( const GeoVector& v )
{
    return 64 * v[0] * ( 1 - v[0] - 2 * v[1] - 2 * v[2] );
}
Real der2fct5_33_P2tilde_3D ( const GeoVector& v )
{
    return - 128 * v[0] * v[1] * v[2];
}

Real der2fct6_11_P2tilde_3D ( const GeoVector& v )
{
    return -128 * v[0] * v[1] * v[2];
}
Real der2fct6_12_P2tilde_3D ( const GeoVector& v )
{
    return 4 + 64 * v[2] * ( 1 - 2 * v[0] - 2 * v[1] - v[2] );
}
Real der2fct6_13_P2tilde_3D ( const GeoVector& v )
{
    return 64 * v[1] * ( 1 - 2 * v[0] - v[1] - 2 * v[2] );
}
Real der2fct6_21_P2tilde_3D ( const GeoVector& v )
{
    return 4 + 64 * v[2] * ( 1 - 2 * v[0] - 2 * v[1] - v[2] );
}
Real der2fct6_22_P2tilde_3D ( const GeoVector& v )
{
    return - 128 * v[0] * v[1] * v[2];
}
Real der2fct6_23_P2tilde_3D ( const GeoVector& v )
{
    return 64 * v[0] * ( 1 - v[0] - 2 * v[1] - 2 * v[2] );
}
Real der2fct6_31_P2tilde_3D ( const GeoVector& v )
{
    return 64 * v[1] * ( 1 - 2 * v[0] - v[1] - 2 * v[2] );
}
Real der2fct6_32_P2tilde_3D ( const GeoVector& v )
{
    return 64 * v[0] * ( 1 - v[0] - 2 * v[1] - 2 * v[2] );
}
Real der2fct6_33_P2tilde_3D ( const GeoVector& v )
{
    return - 128 * v[0] * v[1] * v[2];
}

Real der2fct7_11_P2tilde_3D ( const GeoVector& v )
{
    return -128 * v[0] * v[1] * v[2];
}
Real der2fct7_12_P2tilde_3D ( const GeoVector& v )
{
    return -4 + 64 * v[2] * ( 1 - 2 * v[0] - 2 * v[1] - v[2] );
}
Real der2fct7_13_P2tilde_3D ( const GeoVector& v )
{
    return 64 * v[1] * ( 1 - 2 * v[0] - v[1] - 2 * v[2] );
}
Real der2fct7_21_P2tilde_3D ( const GeoVector& v )
{
    return -4 + 64 * v[2] * ( 1 - 2 * v[0] - 2 * v[1] - v[2] );
}
Real der2fct7_22_P2tilde_3D ( const GeoVector& v )
{
    return -8 - 128 * v[0] * v[1] * v[2];
}
Real der2fct7_23_P2tilde_3D ( const GeoVector& v )
{
    return -4 + 64 * v[0] * ( 1 - v[0] - 2 * v[1] - 2 * v[2] );
}
Real der2fct7_31_P2tilde_3D ( const GeoVector& v )
{
    return 64 * v[1] * ( 1 - 2 * v[0] - v[1] - 2 * v[2] );
}
Real der2fct7_32_P2tilde_3D ( const GeoVector& v )
{
    return -4 + 64 * v[0] * ( 1 - v[0] - 2 * v[1] - 2 * v[2] );
}
Real der2fct7_33_P2tilde_3D ( const GeoVector& v )
{
    return - 128 * v[0] * v[1] * v[2];
}

Real der2fct8_11_P2tilde_3D ( const GeoVector& v )
{
    return -128 * v[0] * v[1] * v[2];
}
Real der2fct8_12_P2tilde_3D ( const GeoVector& v )
{
    return 64 * v[2] * ( 1 - 2 * v[0] - 2 * v[1] - v[2] );
}
Real der2fct8_13_P2tilde_3D ( const GeoVector& v )
{
    return -4 + 64 * v[1] * ( 1 - 2 * v[0] - v[1] - 2 * v[2] );
}
Real der2fct8_21_P2tilde_3D ( const GeoVector& v )
{
    return 64 * v[2] * ( 1 - 2 * v[0] - 2 * v[1] - v[2] );
}
Real der2fct8_22_P2tilde_3D ( const GeoVector& v )
{
    return - 128 * v[0] * v[1] * v[2];
}
Real der2fct8_23_P2tilde_3D ( const GeoVector& v )
{
    return -4 + 64 * v[0] * ( 1 - v[0] - 2 * v[1] - 2 * v[2] );
}
Real der2fct8_31_P2tilde_3D ( const GeoVector& v )
{
    return -4 + 64 * v[1] * ( 1 - 2 * v[0] - v[1] - 2 * v[2] );
}
Real der2fct8_32_P2tilde_3D ( const GeoVector& v )
{
    return -4 + 64 * v[0] * ( 1 - v[0] - 2 * v[1] - 2 * v[2] );
}
Real der2fct8_33_P2tilde_3D ( const GeoVector& v )
{
    return -8 - 128 * v[0] * v[1] * v[2];
}

Real der2fct9_11_P2tilde_3D ( const GeoVector& v )
{
    return -128 * v[0] * v[1] * v[2];
}
Real der2fct9_12_P2tilde_3D ( const GeoVector& v )
{
    return 64 * v[2] * ( 1 - 2 * v[0] - 2 * v[1] - v[2] );
}
Real der2fct9_13_P2tilde_3D ( const GeoVector& v )
{
    return 4 + 64 * v[1] * ( 1 - 2 * v[0] - v[1] - 2 * v[2] );
}
Real der2fct9_21_P2tilde_3D ( const GeoVector& v )
{
    return 64 * v[2] * ( 1 - 2 * v[0] - 2 * v[1] - v[2] );
}
Real der2fct9_22_P2tilde_3D ( const GeoVector& v )
{
    return - 128 * v[0] * v[1] * v[2];
}
Real der2fct9_23_P2tilde_3D ( const GeoVector& v )
{
    return 64 * v[0] * ( 1 - v[0] - 2 * v[1] - 2 * v[2] );
}
Real der2fct9_31_P2tilde_3D ( const GeoVector& v )
{
    return 4 + 64 * v[1] * ( 1 - 2 * v[0] - v[1] - 2 * v[2] );
}
Real der2fct9_32_P2tilde_3D ( const GeoVector& v )
{
    return 64 * v[0] * ( 1 - v[0] - 2 * v[1] - 2 * v[2] );
}
Real der2fct9_33_P2tilde_3D ( const GeoVector& v )
{
    return - 128 * v[0] * v[1] * v[2];
}

Real der2fct10_11_P2tilde_3D ( const GeoVector& v )
{
    return -128 * v[0] * v[1] * v[2];
}
Real der2fct10_12_P2tilde_3D ( const GeoVector& v )
{
    return 64 * v[2] * ( 1 - 2 * v[0] - 2 * v[1] - v[2] );
}
Real der2fct10_13_P2tilde_3D ( const GeoVector& v )
{
    return 64 * v[1] * ( 1 - 2 * v[0] - v[1] - 2 * v[2] );
}
Real der2fct10_21_P2tilde_3D ( const GeoVector& v )
{
    return 64 * v[2] * ( 1 - 2 * v[0] - 2 * v[1] - v[2] );
}
Real der2fct10_22_P2tilde_3D ( const GeoVector& v )
{
    return - 128 * v[0] * v[1] * v[2];
}
Real der2fct10_23_P2tilde_3D ( const GeoVector& v )
{
    return 4 + 64 * v[0] * ( 1 - v[0] - 2 * v[1] - 2 * v[2] );
}
Real der2fct10_31_P2tilde_3D ( const GeoVector& v )
{
    return 64 * v[1] * ( 1 - 2 * v[0] - v[1] - 2 * v[2] );
}
Real der2fct10_32_P2tilde_3D ( const GeoVector& v )
{
    return 4 + 64 * v[0] * ( 1 - v[0] - 2 * v[1] - 2 * v[2] );
}
Real der2fct10_33_P2tilde_3D ( const GeoVector& v )
{
    return - 128 * v[0] * v[1] * v[2];
}

Real der2fct11_11_P2tilde_3D ( const GeoVector& v )
{
    return -512 * v[0] * v[1] * v[2];
}
Real der2fct11_12_P2tilde_3D ( const GeoVector& v )
{
    return 256 * v[2] * ( 1 - 2 * v[0] - 2 * v[1] - v[2] );
}
Real der2fct11_13_P2tilde_3D ( const GeoVector& v )
{
    return 256 * v[1] * ( 1 - 2 * v[0] - v[1] - 2 * v[2] );
}
Real der2fct11_21_P2tilde_3D ( const GeoVector& v )
{
    return 256 * v[2] * ( 1 - 2 * v[0] - 2 * v[1] - v[2] );
}
Real der2fct11_22_P2tilde_3D ( const GeoVector& v )
{
    return -512 * v[0] * v[1] * v[2];
}
Real der2fct11_23_P2tilde_3D ( const GeoVector& v )
{
    return 256 * v[0] * ( 1 - v[0] - 2 * v[1] - 2 * v[2] );
}
Real der2fct11_31_P2tilde_3D ( const GeoVector& v )
{
    return 256 * v[1] * ( 1 - 2 * v[0] - v[1] - 2 * v[2] );
}
Real der2fct11_32_P2tilde_3D ( const GeoVector& v )
{
    return 256 * v[0] * ( 1 - v[0] - 2 * v[1] - 2 * v[2] );
}
Real der2fct11_33_P2tilde_3D ( const GeoVector& v )
{
    return -512 * v[0] * v[1] * v[2];
}

//======================================================================
//
//                            Q0  (3D)
//
//======================================================================
/*
        ________
       /.      /|
      / .     / |
     /_______/  |
     |  .  1 |  |
     |  .....|..|
     | .     | /
     |.      |/
     |_______|

*/
Real fct1_Q0_3D ( const GeoVector& )
{
    return 1.;
}
Real derfct1_Q0_3D ( const GeoVector& )
{
    return 0.;
}
// The second derivative is equal to the first : both are equal to 0.
Real der2fct1_Q0_3D ( const GeoVector& )
{
    return 0.;
}

//======================================================================
//
//                            Q1  (3D)
//
//======================================================================
/*
        8-------7
       /.      /|
      / .     / |
     5_______6  |
     |  .    |  |
     |  4....|..3
     | .     | /
     |.      |/
     1_______2
*/
Real fct1_Q1_3D ( const GeoVector& v )
{
    return ( 1. - v[0] ) * ( 1. - v[1] ) * ( 1. - v[2] );
}
Real fct2_Q1_3D ( const GeoVector& v )
{
    return v[0] * ( 1. - v[1] ) * ( 1. - v[2] );
}
Real fct3_Q1_3D ( const GeoVector& v )
{
    return v[0] * v[1] * ( 1. - v[2] );
}
Real fct4_Q1_3D ( const GeoVector& v )
{
    return ( 1. - v[0] ) * v[1] * ( 1. - v[2] );
}
Real fct5_Q1_3D ( const GeoVector& v )
{
    return ( 1. - v[0] ) * ( 1. - v[1] ) * v[2];
}
Real fct6_Q1_3D ( const GeoVector& v )
{
    return v[0] * ( 1. - v[1] ) * v[2];
}
Real fct7_Q1_3D ( const GeoVector& v )
{
    return v[0] * v[1] * v[2];
}
Real fct8_Q1_3D ( const GeoVector& v )
{
    return ( 1. - v[0] ) * v[1] * v[2];
}

Real derfct1_1_Q1_3D ( const GeoVector& v )
{
    return - ( 1. - v[1] ) * ( 1. - v[2] );
}
Real derfct1_2_Q1_3D ( const GeoVector& v )
{
    return - ( 1. - v[0] ) * ( 1. - v[2] );
}
Real derfct1_3_Q1_3D ( const GeoVector& v )
{
    return - ( 1. - v[0] ) * ( 1. - v[1] );
}
Real derfct2_1_Q1_3D ( const GeoVector& v )
{
    return ( 1. - v[1] ) * ( 1. - v[2] );
}
Real derfct2_2_Q1_3D ( const GeoVector& v )
{
    return -v[0] * ( 1. - v[2] ) ;
}
Real derfct2_3_Q1_3D ( const GeoVector& v )
{
    return -v[0] * ( 1. - v[1] );
}
Real derfct3_1_Q1_3D ( const GeoVector& v )
{
    return v[1] * ( 1. - v[2] );
}
Real derfct3_2_Q1_3D ( const GeoVector& v )
{
    return v[0] * ( 1. - v[2] );
}
Real derfct3_3_Q1_3D ( const GeoVector& v )
{
    return -v[0] * v[1] ;
}
Real derfct4_1_Q1_3D ( const GeoVector& v )
{
    return -v[1] * ( 1. - v[2] );
}
Real derfct4_2_Q1_3D ( const GeoVector& v )
{
    return ( 1. - v[0] ) * ( 1. - v[2] );
}
Real derfct4_3_Q1_3D ( const GeoVector& v )
{
    return - ( 1. - v[0] ) * v[1];
}
Real derfct5_1_Q1_3D ( const GeoVector& v )
{
    return - ( 1. - v[1] ) * v[2];
}
Real derfct5_2_Q1_3D ( const GeoVector& v )
{
    return - ( 1. - v[0] ) * v[2];
}
Real derfct5_3_Q1_3D ( const GeoVector& v )
{
    return ( 1. - v[0] ) * ( 1. - v[1] );
}
Real derfct6_1_Q1_3D ( const GeoVector& v )
{
    return ( 1. - v[1] ) * v[2] ;
}
Real derfct6_2_Q1_3D ( const GeoVector& v )
{
    return -v[0] * v[2];
}
Real derfct6_3_Q1_3D ( const GeoVector& v )
{
    return v[0] * ( 1. - v[1] );
}
Real derfct7_1_Q1_3D ( const GeoVector& v )
{
    return v[1] * v[2];
}
Real derfct7_2_Q1_3D ( const GeoVector& v )
{
    return v[0] * v[2];
}
Real derfct7_3_Q1_3D ( const GeoVector& v )
{
    return v[0] * v[1];
}
Real derfct8_1_Q1_3D ( const GeoVector& v )
{
    return -v[1] * v[2];
}
Real derfct8_2_Q1_3D ( const GeoVector& v )
{
    return ( 1. - v[0] ) * v[2];
}
Real derfct8_3_Q1_3D ( const GeoVector& v )
{
    return ( 1. - v[0] ) * v[1];
}

Real der2fct1_11_Q1_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct1_12_Q1_3D ( const GeoVector& v )
{
    return 1. - v[2];
}
Real der2fct1_13_Q1_3D ( const GeoVector& v )
{
    return 1. - v[1];
}
Real der2fct1_21_Q1_3D ( const GeoVector& v )
{
    return 1. - v[2];
}
Real der2fct1_22_Q1_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct1_23_Q1_3D ( const GeoVector& v )
{
    return 1. - v[0];
}
Real der2fct1_31_Q1_3D ( const GeoVector& v )
{
    return 1. - v[1];
}
Real der2fct1_32_Q1_3D ( const GeoVector& v )
{
    return 1. - v[0];
}
Real der2fct1_33_Q1_3D ( const GeoVector& )
{
    return 0;
}

Real der2fct2_11_Q1_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct2_12_Q1_3D ( const GeoVector& v )
{
    return - ( 1. - v[2] );
}
Real der2fct2_13_Q1_3D ( const GeoVector& v )
{
    return - ( 1. - v[1] );
}
Real der2fct2_21_Q1_3D ( const GeoVector& v )
{
    return - ( 1. - v[2] );
}
Real der2fct2_22_Q1_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct2_23_Q1_3D ( const GeoVector& v )
{
    return v[0];
}
Real der2fct2_31_Q1_3D ( const GeoVector& v )
{
    return - ( 1. - v[1] );
}
Real der2fct2_32_Q1_3D ( const GeoVector& v )
{
    return v[0];
}
Real der2fct2_33_Q1_3D ( const GeoVector& )
{
    return 0;
}

Real der2fct3_11_Q1_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct3_12_Q1_3D ( const GeoVector& v )
{
    return ( 1. - v[2] );
}
Real der2fct3_13_Q1_3D ( const GeoVector& v )
{
    return -v[1];
}
Real der2fct3_21_Q1_3D ( const GeoVector& v )
{
    return ( 1. - v[2] );
}
Real der2fct3_22_Q1_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct3_23_Q1_3D ( const GeoVector& v )
{
    return -v[0];
}
Real der2fct3_31_Q1_3D ( const GeoVector& v )
{
    return -v[1];
}
Real der2fct3_32_Q1_3D ( const GeoVector& v )
{
    return -v[0];
}
Real der2fct3_33_Q1_3D ( const GeoVector& )
{
    return 0;
}

Real der2fct4_11_Q1_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct4_12_Q1_3D ( const GeoVector& v )
{
    return - ( 1. - v[2] );
}
Real der2fct4_13_Q1_3D ( const GeoVector& v )
{
    return v[1];
}
Real der2fct4_21_Q1_3D ( const GeoVector& v )
{
    return - ( 1. - v[2] );
}
Real der2fct4_22_Q1_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct4_23_Q1_3D ( const GeoVector& v )
{
    return - ( 1. - v[0] );
}
Real der2fct4_31_Q1_3D ( const GeoVector& v )
{
    return v[1];
}
Real der2fct4_32_Q1_3D ( const GeoVector& v )
{
    return - ( 1. - v[0] );
}
Real der2fct4_33_Q1_3D ( const GeoVector& )
{
    return 0;
}

Real der2fct5_11_Q1_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct5_12_Q1_3D ( const GeoVector& v )
{
    return v[2];
}
Real der2fct5_13_Q1_3D ( const GeoVector& v )
{
    return - ( 1. - v[1] );
}
Real der2fct5_21_Q1_3D ( const GeoVector& v )
{
    return v[2];
}
Real der2fct5_22_Q1_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct5_23_Q1_3D ( const GeoVector& v )
{
    return - ( 1. - v[0] );
}
Real der2fct5_31_Q1_3D ( const GeoVector& v )
{
    return - ( 1. - v[1] );
}
Real der2fct5_32_Q1_3D ( const GeoVector& v )
{
    return - ( 1. - v[0] );
}
Real der2fct5_33_Q1_3D ( const GeoVector& )
{
    return 0;
}

Real der2fct6_11_Q1_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct6_12_Q1_3D ( const GeoVector& v )
{
    return -v[2];
}
Real der2fct6_13_Q1_3D ( const GeoVector& v )
{
    return 1. - v[1];
}
Real der2fct6_21_Q1_3D ( const GeoVector& v )
{
    return -v[2];
}
Real der2fct6_22_Q1_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct6_23_Q1_3D ( const GeoVector& v )
{
    return -v[0];
}
Real der2fct6_31_Q1_3D ( const GeoVector& v )
{
    return 1. - v[1];
}
Real der2fct6_32_Q1_3D ( const GeoVector& v )
{
    return -v[0];
}
Real der2fct6_33_Q1_3D ( const GeoVector& )
{
    return 0;
}

Real der2fct7_11_Q1_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct7_12_Q1_3D ( const GeoVector& v )
{
    return v[2];
}
Real der2fct7_13_Q1_3D ( const GeoVector& v )
{
    return v[1];
}
Real der2fct7_21_Q1_3D ( const GeoVector& v )
{
    return v[2];
}
Real der2fct7_22_Q1_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct7_23_Q1_3D ( const GeoVector& v )
{
    return v[0];
}
Real der2fct7_31_Q1_3D ( const GeoVector& v )
{
    return v[1];
}
Real der2fct7_32_Q1_3D ( const GeoVector& v )
{
    return v[0];
}
Real der2fct7_33_Q1_3D ( const GeoVector& )
{
    return 0;
}

Real der2fct8_11_Q1_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct8_12_Q1_3D ( const GeoVector& v )
{
    return -v[2];
}
Real der2fct8_13_Q1_3D ( const GeoVector& v )
{
    return -v[1];
}
Real der2fct8_21_Q1_3D ( const GeoVector& v )
{
    return -v[2];
}
Real der2fct8_22_Q1_3D ( const GeoVector& )
{
    return 0;
}
Real der2fct8_23_Q1_3D ( const GeoVector& v )
{
    return 1. - v[0];
}
Real der2fct8_31_Q1_3D ( const GeoVector& v )
{
    return -v[1];
}
Real der2fct8_32_Q1_3D ( const GeoVector& v )
{
    return 1 - v[0];
}
Real der2fct8_33_Q1_3D ( const GeoVector& )
{
    return 0;
}

//======================================================================
//
//                            RT0 Hexa  (3D)
//
//======================================================================
/*

        8-------7
       /.      /|
      / .     / |
     5_______6  |
     |  .    |  |
     |  4....|..3
     | .     | /
     |.      |/
     1_______2

   face 1: 1,4,3,2
   face 2: 1,5,8,4
   face 3: 1,2,6,5
   face 4: 2,3,7,6
   face 5: 3,4,8,7
   face 6: 5,6,7,8

*/

Real fct1_RT0_1_HEXA_3D ( const GeoVector& )
{
    return 0.;
}
Real fct1_RT0_2_HEXA_3D ( const GeoVector& )
{
    return 0.;
}
Real fct1_RT0_3_HEXA_3D ( const GeoVector& v )
{
    return v[2] - 1.;
}

Real fct2_RT0_1_HEXA_3D ( const GeoVector& v )
{
    return v[0] - 1.;
}
Real fct2_RT0_2_HEXA_3D ( const GeoVector& )
{
    return 0.;
}
Real fct2_RT0_3_HEXA_3D ( const GeoVector& )
{
    return 0.;
}

Real fct3_RT0_1_HEXA_3D ( const GeoVector& )
{
    return 0.;
}
Real fct3_RT0_2_HEXA_3D ( const GeoVector& v )
{
    return v[1] - 1.;
}
Real fct3_RT0_3_HEXA_3D ( const GeoVector& )
{
    return 0.;
}

Real fct4_RT0_1_HEXA_3D ( const GeoVector& v )
{
    return v[0];
}
Real fct4_RT0_2_HEXA_3D ( const GeoVector& )
{
    return 0.;
}
Real fct4_RT0_3_HEXA_3D ( const GeoVector& )
{
    return 0.;
}

Real fct5_RT0_1_HEXA_3D ( const GeoVector& )
{
    return 0.;
}
Real fct5_RT0_2_HEXA_3D ( const GeoVector& v )
{
    return v[1];
}
Real fct5_RT0_3_HEXA_3D ( const GeoVector& )
{
    return 0.;
}

Real fct6_RT0_1_HEXA_3D ( const GeoVector& )
{
    return 0.;
}
Real fct6_RT0_2_HEXA_3D ( const GeoVector& )
{
    return 0.;
}
Real fct6_RT0_3_HEXA_3D ( const GeoVector& v )
{
    return v[2];
}

Real fct1_DIV_RT0_HEXA_3D ( const GeoVector& )
{
    return 1.;
}
Real fct2_DIV_RT0_HEXA_3D ( const GeoVector& )
{
    return 1.;
}
Real fct3_DIV_RT0_HEXA_3D ( const GeoVector& )
{
    return 1.;
}
Real fct4_DIV_RT0_HEXA_3D ( const GeoVector& )
{
    return 1.;
}
Real fct5_DIV_RT0_HEXA_3D ( const GeoVector& )
{
    return 1.;
}
Real fct6_DIV_RT0_HEXA_3D ( const GeoVector& )
{
    return 1.;
}

//======================================================================
//
//                            RT0 Tetra  (3D)
//
//======================================================================

/*

                4
               / .
              /  \.3
             /  . \\
            / .    \\
           /.       \!
         1 ----------2

SEE ElementShapes.cc   for the ORIENTATION CONVENTIONS
   point 1: 0, 0, 0
   point 2: 1, 0, 0
   point 3: 0, 1, 0
   point 4: 0, 0, 1

   face 1: 2, 3, 4
   face 2: 1, 4, 3
   face 3: 1, 2, 4
   face 4: 1, 3, 2


*/

Real fct3_RT0_1_TETRA_3D ( const GeoVector& v )
{
    return 2. * v[0];
}
Real fct3_RT0_2_TETRA_3D ( const GeoVector& v )
{
    return 2. * v[1];
}
Real fct3_RT0_3_TETRA_3D ( const GeoVector& v )
{
    return 2. * v[2];
}

Real fct4_RT0_1_TETRA_3D ( const GeoVector& v )
{
    return 2. * v[0] - 2.;
}
Real fct4_RT0_2_TETRA_3D ( const GeoVector& v )
{
    return 2. * v[1];
}
Real fct4_RT0_3_TETRA_3D ( const GeoVector& v )
{
    return 2. * v[2];
}

Real fct2_RT0_1_TETRA_3D ( const GeoVector& v )
{
    return 2. * v[0];
}
Real fct2_RT0_2_TETRA_3D ( const GeoVector& v )
{
    return 2. * v[1] - 2.;
}
Real fct2_RT0_3_TETRA_3D ( const GeoVector& v )
{
    return 2. * v[2];
}

Real fct1_RT0_1_TETRA_3D ( const GeoVector& v )
{
    return 2. * v[0];
}
Real fct1_RT0_2_TETRA_3D ( const GeoVector& v )
{
    return 2. * v[1];
}
Real fct1_RT0_3_TETRA_3D ( const GeoVector& v )
{
    return 2. * v[2] - 2.;
}

Real fct1_DIV_RT0_TETRA_3D ( const GeoVector& )
{
    return 6.;
}
Real fct2_DIV_RT0_TETRA_3D ( const GeoVector& )
{
    return 6.;
}
Real fct3_DIV_RT0_TETRA_3D ( const GeoVector& )
{
    return 6.;
}
Real fct4_DIV_RT0_TETRA_3D ( const GeoVector& )
{
    return 6.;
}

// Transformation functions

std::vector<Real> lagrangianTransform (const std::vector<Real>& values)
{
    return values;
}

std::vector<Real> P1Bubble3DTransform (const std::vector<Real>& nodalValues)
{
    std::vector<Real> FEValues (nodalValues);
    FEValues[4] = 256 * nodalValues[4] - 64 * (nodalValues[0] + nodalValues[1] + nodalValues[2] + nodalValues[3]);
    return FEValues;
}

std::vector<Real> P1Bubble2DTransform (const std::vector<Real>& nodalValues)
{
    std::vector<Real> FEValues (nodalValues);
    FEValues[3] = 27 * nodalValues[3] - 9 * (nodalValues[0] + nodalValues[1] + nodalValues[2]);
    return FEValues;
}


//======================================================================
//
//                            P0  (0D)
//
//======================================================================
/*
                           1
*/

const ReferenceFEScalar fePointP0 ( "Lagrange P0 on a point",
                                    FE_P0_0D,
                                    POINT,
                                    1,                           // nb dof per vertex
                                    0,                           // nb dof per edge
                                    0,                           // nb dof per face
                                    0,                           // nb dof per volume
                                    1,                           // nb dof
                                    1,                           // nb coor
                                    fct_P0_0D,
                                    derfct_P0_0D,
                                    der2fct_P0_0D,
                                    refcoor_P0_0D,
                                    STANDARD_PATTERN,
                                    ( ReferenceFE* ) NULL,
                                    &lagrangianTransform );

//======================================================================
//
//                            P0  (1D)
//
//======================================================================
/*
                           --1--
*/

const ReferenceFEScalar feSegP0 ( "Lagrange P0 on a segment", FE_P0_1D, LINE, 0, 1, 0, 0, 1, 1,
                                  fct_P0_1D, derfct_P0_1D, der2fct_P0_1D, refcoor_P0_1D,
                                  STANDARD_PATTERN, &fePointP0, &lagrangianTransform );

//======================================================================
//
//                            P1  (1D)
//
//======================================================================
/*
                           1-----2
*/

const ReferenceFEScalar feSegP1 ( "Lagrange P1 on a segment", FE_P1_1D, LINE, 1, 0, 0, 0, 2, 1,
                                  fct_P1_1D, derfct_P1_1D, der2fct_P1_1D, refcoor_P1_1D,
                                  STANDARD_PATTERN, &fePointP0, &lagrangianTransform );

//======================================================================
//
//                            P2  (1D)
//
//======================================================================
/*
                           1--3--2
*/

const ReferenceFEScalar feSegP2 ( "Lagrange P2 on a segment", FE_P2_1D, LINE, 1, 1, 0, 0, 3, 1,
                                  fct_P2_1D, derfct_P2_1D, der2fct_P2_1D, refcoor_P2_1D,
                                  STANDARD_PATTERN, &fePointP0, &lagrangianTransform );

//======================================================================
//
//                            P0  (2D)
//
//======================================================================
/*

                           |\
                           | \
                           | 1\
                            ---
*/

const ReferenceFEScalar feTriaP0 ( "Lagrange P0 on a triangle", FE_P0_2D, TRIANGLE, 0, 0, 1, 0, 1, 2,
                                   fct_P0_2D, derfct_P0_2D, der2fct_P0_2D, refcoor_P0_2D,
                                   STANDARD_PATTERN, &feSegP0, &lagrangianTransform );

//======================================================================
//
//                            P1  (2D)
//
//======================================================================
/*
                           3
                           |\
                           | \
                           |  \
                           1---2
*/

const ReferenceFEScalar feTriaP1 ( "Lagrange P1 on a triangle", FE_P1_2D, TRIANGLE, 1, 0, 0, 0, 3, 2,
                                   fct_P1_2D, derfct_P1_2D, der2fct_P1_2D, refcoor_P1_2D,
                                   STANDARD_PATTERN, &feSegP1, &lagrangianTransform );


//======================================================================
//
//                            P1bubble  (2D)
//
//======================================================================
/*
                           3
                           |\
                           | \
                           |4.\
                           1---2
*/

const ReferenceFEScalar feTriaP1bubble ( "P1bubble on a triangle", FE_P1bubble_2D, TRIANGLE, 1, 0, 1, 0, 4, 2,
                                         fct_P1bubble_2D, derfct_P1bubble_2D, der2fct_P1bubble_2D, refcoor_P1bubble_2D,
                                         STANDARD_PATTERN, &feSegP1, &P1Bubble2DTransform );


//======================================================================
//
//                            P2  (2D)
//
//======================================================================
/*
                           3
                           |\
                           6 5
                           |  \
                           1-4-2
*/

const ReferenceFEScalar feTriaP2 ( "Lagrange P2 on a triangle", FE_P2_2D, TRIANGLE, 1, 1, 0, 0, 6, 2,
                                   fct_P2_2D, derfct_P2_2D, der2fct_P2_2D, refcoor_P2_2D,
                                   STANDARD_PATTERN, &feSegP2, &lagrangianTransform );

//======================================================================
//
//                            RT0 (2D)
//
//======================================================================
/*
                           +
                           |\
                           3 2
                           |  \
                           +-1-+
*/

const ReferenceFEHdiv feTriaRT0 ( "Lagrange RT0 on a triangle", FE_RT0_TRIA_2D, TRIANGLE, 0, 1, 0, 0, 3, 2,
                                  fct_RT0_TRIA_2D, fct_DIV_RT0_TRIA_2D, refcoor_RT0_TRIA_2D,
                                  STANDARD_PATTERN, &feSegP0 );


//======================================================================
//
//                            Q0  (2D)
//
//======================================================================
/*
                            -------
                           |       |
                           |   1   |
                           |       |
                            -------
*/

const ReferenceFEScalar feQuadQ0 ( "Lagrange Q0 on a quadrangle", FE_Q0_2D, QUAD, 0, 0, 1, 0, 1, 2,
                                   fct_Q0_2D, derfct_Q0_2D, der2fct_Q0_2D, refcoor_Q0_2D,
                                   STANDARD_PATTERN, &feSegP0, &lagrangianTransform );

//======================================================================
//
//                            Q1  (2D)
//
//======================================================================
/*
                           4-------3
                           |       |
                           |       |
                           |       |
                           1-------2
*/

const ReferenceFEScalar feQuadQ1 ( "Lagrange Q1 on a quadrangle", FE_Q1_2D, QUAD, 1, 0, 0, 0, 4, 2,
                                   fct_Q1_2D, derfct_Q1_2D, der2fct_Q1_2D, refcoor_Q1_2D,
                                   STANDARD_PATTERN, &feSegP1, &lagrangianTransform );


//======================================================================
//
//                            Q2  (2D)
//
//======================================================================
/*
                           4---7---3
                           |       |
                           8   9   6
                           |       |
                           1---5---2
*/

const ReferenceFEScalar feQuadQ2 ( "Lagrange Q2 on a quadrangle", FE_Q2_2D, QUAD, 1, 1, 1, 0, 9, 2,
                                   fct_Q2_2D, derfct_Q2_2D, der2fct_Q2_2D, refcoor_Q2_2D,
                                   STANDARD_PATTERN, &feSegP2, &lagrangianTransform );

//======================================================================
//
//                            P0  (3D)
//
//======================================================================
/*

               / .
              /  \.
             /  . \\
            / . 1  \\
           /.       \!
           ----------
*/
const ReferenceFEScalar feTetraP0 ( "Lagrange P0 on a tetraedra", FE_P0_3D, TETRA, 0, 0, 0, 1, 1, 3,
                                    fct_P0_3D, derfct_P0_3D, der2fct_P0_3D, refcoor_P0_3D,
                                    STANDARD_PATTERN, &feTriaP0, &lagrangianTransform );

//======================================================================
//
//                            P1  (3D)
//
//======================================================================
/*
                4
               / .
              /  \.3
             /  . \\
            / .    \\
           /.       \!
         1 ----------2
*/
const ReferenceFEScalar feTetraP1 ( "Lagrange P1 on a tetraedra", FE_P1_3D, TETRA, 1, 0, 0, 0, 4, 3,
                                    fct_P1_3D, derfct_P1_3D, der2fct_P1_3D, refcoor_P1_3D,
                                    STANDARD_PATTERN, &feTriaP1, &lagrangianTransform );

//======================================================================
//
//                            P1bubble  (3D)
//
//======================================================================
/*
                4
               / .
              /  \.3
             /  . \\
            / . .5 \\
           /.       \!
         1 ----------2
*/
const ReferenceFEScalar feTetraP1bubble ( "Lagrange P1bubble on a tetraedra", FE_P1bubble_3D, TETRA, 1, 0, 0, 1, 5, 3,
                                          fct_P1bubble_3D, derfct_P1bubble_3D, der2fct_P1bubble_3D, refcoor_P1bubble_3D,
                                          STANDARD_PATTERN, &feTriaP1, &P1Bubble3DTransform );


//======================================================================
//
//                            P2  (3D)
//
//======================================================================
/*
                4
               / .10
              /  \.3
             8  . 9\
            / 7    \6
           /.       \!
         1 -----5----2
*/
const ReferenceFEScalar feTetraP2 ( "Lagrange P2 on a tetraedra", FE_P2_3D, TETRA, 1, 1, 0, 0, 10, 3,
                                    fct_P2_3D, derfct_P2_3D, der2fct_P2_3D, refcoor_P2_3D,
                                    STANDARD_PATTERN, &feTriaP2, &lagrangianTransform );
//======================================================================
//
//                            P2tilde  (3D)
// NAVIER-STOKES P2 Basis Oriented to the mass lumping
//======================================================================
/*
                4
               / .10
              /  \.3
             8  . 9\
            / 7 .11 \6
           /.       \!
         1 -----5----2
*/
const ReferenceFEScalar feTetraP2tilde ( "Lagrange P2tilde on a tetraedra", FE_P2tilde_3D,
                                         TETRA, 1, 1, 0, 1, 11, 3, fct_P2tilde_3D,
                                         derfct_P2tilde_3D,
                                         der2fct_P2tilde_3D,
                                         refcoor_P2tilde_3D,
                                         STANDARD_PATTERN, &feTriaP2, &lagrangianTransform );

//======================================================================
//
//                            Q0  (3D)
//
//======================================================================
/*
        ________
       /.      /|
      / .     / |
     /_______/  |
     |  .  1 |  |
     |  .....|..|
     | .     | /
     |.      |/
     |_______|
*/
const ReferenceFEScalar feHexaQ0 ( "Lagrange Q0 on a hexaedra", FE_Q0_3D, HEXA, 0, 0, 0, 1, 1, 3,
                                   fct_Q0_3D, derfct_Q0_3D, der2fct_Q0_3D, refcoor_Q0_3D,
                                   STANDARD_PATTERN, &feQuadQ0, &lagrangianTransform );

//======================================================================
//
//                            Q1  (3D)
//
//======================================================================
/*
        8-------7
       /.      /|
      / .     / |
     5_______6  |
     |  .    |  |
     |  4....|..3
     | .     | /
     |.      |/
     1_______2
*/
const ReferenceFEScalar feHexaQ1 ( "Lagrange Q1 on a hexaedra", FE_Q1_3D, HEXA, 1, 0, 0, 0, 8, 3,
                                   fct_Q1_3D, derfct_Q1_3D, der2fct_Q1_3D, refcoor_Q1_3D,
                                   STANDARD_PATTERN, &feQuadQ1, &lagrangianTransform );

//======================================================================
//
//                            RT0 (3D)
//
//======================================================================
/*

        8-------7
       /.      /|
      / .     / |
     5_______6  |
     |  .    |  |
     |  4....|..3
     | .     | /
     |.      |/
     1_______2

   face 1: 1,4,3,2
   face 2: 1,5,8,4
   face 3: 1,2,6,5
   face 4: 2,3,7,6
   face 5: 3,4,8,7
   face 6: 5,6,7,8

*/
const ReferenceFEHdiv feHexaRT0 ( "Lagrange RT0 on a hexaedra", FE_RT0_HEXA_3D, HEXA, 0, 0, 1, 0, 6, 3,
                                  fct_RT0_HEXA_3D, fct_DIV_RT0_HEXA_3D, refcoor_RT0_HEXA_3D,
                                  STANDARD_PATTERN, &feQuadQ0);

//======================================================================
//
//                            RT0 (3D)
//
//======================================================================
/*
                4
               / .
              /  \.3
             /  . \\
            / .    \\
           /.       \!
         1 ----------2


   face 1: 1, 3, 2
   face 2: 1, 2, 4
   face 3: 2, 3, 4
   face 4: 1, 4, 3
*/
const ReferenceFEHdiv feTetraRT0 ( "Lagrange RT0 on a tetraedra", FE_RT0_TETRA_3D, TETRA, 0, 0, 1, 0, 4, 3,
                                   fct_RT0_TETRA_3D, fct_DIV_RT0_TETRA_3D, refcoor_RT0_TETRA_3D,
                                   STANDARD_PATTERN, &feTriaP0 );


//----------------------------------------------------------------------
//
//                       Mixed Hybrid FE
//
//----------------------------------------------------------------------

//======================================================================
//
//                           RT0 TRIA HYBRID (2D)
//                Element defined on SEG :  P0 on each TRIA face.
//
//======================================================================
/*!


*/
// N.B. : the hybrid classes and arrays depend on the quadrature rules,
//        geometric mappings and other reference elements :
//        thus they must be defined AFTER the definitions of quadrule, geomap, refFE...

//! Total number of Boundary elements for the hybrid MFE for TRIA (= Number of faces. common for RT0,RT1...)
#define NB_BDFE_RT0_HYB_TRIA 3
static const CurrentFEManifold BdFE_RT0_HYB_TRIA_1 ( feSegP0, geoLinearSeg, quadRuleSeg1pt,
                                                     refcoor_HYB_TRIA_SEG_1, 0 );
static const CurrentFEManifold BdFE_RT0_HYB_TRIA_2 ( feSegP0, geoLinearSeg, quadRuleSeg1pt,
                                                     refcoor_HYB_TRIA_SEG_2, 1 );
static const CurrentFEManifold BdFE_RT0_HYB_TRIA_3 ( feSegP0, geoLinearSeg, quadRuleSeg1pt,
                                                     refcoor_HYB_TRIA_SEG_3, 2 );

static const CurrentFEManifold* HybRT0TriaList[ NB_BDFE_RT0_HYB_TRIA ] =
{
    &BdFE_RT0_HYB_TRIA_1, &BdFE_RT0_HYB_TRIA_2, &BdFE_RT0_HYB_TRIA_3
};

static const CurrentFEManifold BdFE_RT0_HYB_TRIA_VdotN_1 ( feSegP0, geoLinearSeg, quadRuleSeg1pt,
                                                           refcoor_HYB_TRIA_SEG_1, 0, 1. );
static const CurrentFEManifold BdFE_RT0_HYB_TRIA_VdotN_2 ( feSegP0, geoLinearSeg, quadRuleSeg1pt,
                                                           refcoor_HYB_TRIA_SEG_2, 1, 1. / std::sqrt ( 2. ) );
static const CurrentFEManifold BdFE_RT0_HYB_TRIA_VdotN_3 ( feSegP0, geoLinearSeg, quadRuleSeg1pt,
                                                           refcoor_HYB_TRIA_SEG_3, 2, 1. );

static const CurrentFEManifold* HybRT0TriaVdotNList[ NB_BDFE_RT0_HYB_TRIA ] =
{
    &BdFE_RT0_HYB_TRIA_VdotN_1, &BdFE_RT0_HYB_TRIA_VdotN_2, &BdFE_RT0_HYB_TRIA_VdotN_3
};

const ReferenceFEHybrid feTriaRT0Hyb ( "Hybrid RT0 elements on a triangle", FE_RT0_HYB_TRIA_2D, TRIANGLE,
                                       0, 1, 0, 0, 3, 2, NB_BDFE_RT0_HYB_TRIA, HybRT0TriaList,
                                       refcoor_RT0HYB_TRIA, STANDARD_PATTERN );

const ReferenceFEHybrid feTriaRT0VdotNHyb ( "Hybrid RT0 elements on a triangle", FE_RT0_HYB_TRIA_2D, TRIANGLE,
                                            0, 1, 0, 0, 3, 2, NB_BDFE_RT0_HYB_TRIA, HybRT0TriaVdotNList,
                                            refcoor_RT0HYB_TRIA, STANDARD_PATTERN );


//======================================================================
//
//                           RT0 HYBRID (3D)
//                Element defined on FACES :  Q0 on each QUAD face.
//
//======================================================================
/*!

        8-------7
       /.      /|
      / .     / |
     5_______6  |
     |  .    |  |
     |  4....|..3
     | .     | /
     |.      |/
     1_______2

   face 1: 1,4,3,2
   face 2: 1,5,8,4
   face 3: 1,2,6,5
   face 4: 2,3,7,6
   face 5: 3,4,8,7
   face 6: 5,6,7,8


*/

// N.B. : the hybrid classes and arrays depend on the quadrature rules,
//        geometrical mappings and other reference elements :
//        thus they must be defined AFTER the definitions of quadrule, geomap, refFE...

//! Total number of Boundary elements for the hybrid MFE for HEXA (= Number of faces, common for RT0,RT1...)
#define NB_BDFE_HYB_HEXA 6
static const CurrentFEManifold BdFE_RT0_HYB_HEXA_1 ( feQuadQ0, geoBilinearQuad, quadRuleQuad1pt,
                                                     refcoor_HYB_HEXA_FACE_1, 0 );
static const CurrentFEManifold BdFE_RT0_HYB_HEXA_2 ( feQuadQ0, geoBilinearQuad, quadRuleQuad1pt,
                                                     refcoor_HYB_HEXA_FACE_2, 1 );
static const CurrentFEManifold BdFE_RT0_HYB_HEXA_3 ( feQuadQ0, geoBilinearQuad, quadRuleQuad1pt,
                                                     refcoor_HYB_HEXA_FACE_3, 2 );
static const CurrentFEManifold BdFE_RT0_HYB_HEXA_4 ( feQuadQ0, geoBilinearQuad, quadRuleQuad1pt,
                                                     refcoor_HYB_HEXA_FACE_4, 3 );
static const CurrentFEManifold BdFE_RT0_HYB_HEXA_5 ( feQuadQ0, geoBilinearQuad, quadRuleQuad1pt,
                                                     refcoor_HYB_HEXA_FACE_5, 4 );
static const CurrentFEManifold BdFE_RT0_HYB_HEXA_6 ( feQuadQ0, geoBilinearQuad, quadRuleQuad1pt,
                                                     refcoor_HYB_HEXA_FACE_6, 5 );

static const CurrentFEManifold* HybRT0HexaList[ NB_BDFE_HYB_HEXA ] =
{
    &BdFE_RT0_HYB_HEXA_1, &BdFE_RT0_HYB_HEXA_2,
    &BdFE_RT0_HYB_HEXA_3, &BdFE_RT0_HYB_HEXA_4,
    &BdFE_RT0_HYB_HEXA_5, &BdFE_RT0_HYB_HEXA_6
};

//const RefHybridFE feHexaRT0Hyb(NB_BDFE_HYB_HEXA,HybRT0HexaList,"Hybrid RT0 elements on a hexaedra",
//         FE_RT0_HYB_HEXA_3D, HEXA, 0,0,1,0,6,3,
//         refcoor_RT0HYB_HEXA,STANDARD_PATTERN);
static const CurrentFEManifold BdFE_RT0_HYB_HEXA_VdotN_1 ( feQuadQ0, geoBilinearQuad, quadRuleQuad1pt,
                                                           refcoor_HYB_HEXA_FACE_1, 0, 1. );
static const CurrentFEManifold BdFE_RT0_HYB_HEXA_VdotN_2 ( feQuadQ0, geoBilinearQuad, quadRuleQuad1pt,
                                                           refcoor_HYB_HEXA_FACE_2, 1, 1. );
static const CurrentFEManifold BdFE_RT0_HYB_HEXA_VdotN_3 ( feQuadQ0, geoBilinearQuad, quadRuleQuad1pt,
                                                           refcoor_HYB_HEXA_FACE_3, 2, 1. );
static const CurrentFEManifold BdFE_RT0_HYB_HEXA_VdotN_4 ( feQuadQ0, geoBilinearQuad, quadRuleQuad1pt,
                                                           refcoor_HYB_HEXA_FACE_4, 3, 1. );
static const CurrentFEManifold BdFE_RT0_HYB_HEXA_VdotN_5 ( feQuadQ0, geoBilinearQuad, quadRuleQuad1pt,
                                                           refcoor_HYB_HEXA_FACE_5, 4, 1. );
static const CurrentFEManifold BdFE_RT0_HYB_HEXA_VdotN_6 ( feQuadQ0, geoBilinearQuad, quadRuleQuad1pt,
                                                           refcoor_HYB_HEXA_FACE_6, 5, 1. );

static const CurrentFEManifold* HybRT0HexaVdotNList[ NB_BDFE_HYB_HEXA ] =
{
    &BdFE_RT0_HYB_HEXA_VdotN_1, &BdFE_RT0_HYB_HEXA_VdotN_2,
    &BdFE_RT0_HYB_HEXA_VdotN_3, &BdFE_RT0_HYB_HEXA_VdotN_4,
    &BdFE_RT0_HYB_HEXA_VdotN_5, &BdFE_RT0_HYB_HEXA_VdotN_6
};

const ReferenceFEHybrid feHexaRT0Hyb ( "Hybrid RT0 elements on a hexaedra", FE_RT0_HYB_HEXA_3D, HEXA,
                                       0, 0, 1, 0, 6, 3, NB_BDFE_HYB_HEXA, HybRT0HexaList,
                                       refcoor_RT0HYB_HEXA, STANDARD_PATTERN );

const ReferenceFEHybrid feHexaRT0VdotNHyb ( "Hybrid RT0 elements on a hexaedra", FE_RT0_HYB_HEXA_3D, HEXA,
                                            0, 0, 1, 0, 6, 3, NB_BDFE_HYB_HEXA, HybRT0HexaVdotNList,
                                            refcoor_RT0HYB_HEXA, STANDARD_PATTERN );


//======================================================================
//
//                           RT0 TETRA HYBRID (3D)
//                Element defined on FACES :  Q0 on each QUAD face.
//
//======================================================================
/*!

                4
               / .
              /  \.3
             /  . \\
            / .    \\
           /.       \!
         1 ----------2

SEE ElementShapes.cc   for the ORIENTATION CONVENTIONS
   point 1: 0, 0, 0
   point 2: 1, 0, 0
   point 3: 0, 1, 0
   point 4: 0, 0, 1

   face 1: 1, 3, 2
   face 2: 1, 2, 4
   face 3: 2, 3, 4
   face 4: 1, 4, 3


*/
// N.B. : the hybrid classes and arrays depend on the quadrature rules,
//        geometric mappings and other reference elements :
//        thus they must be defined AFTER the definitions of quadrule, geomap, refFE...

//! Total number of Boundary elements for the hybrid MFE for TETRA (= Number of faces. common for RT0,RT1...)
#define NB_BDFE_RT0_HYB_TETRA 4
static const CurrentFEManifold BdFE_RT0_HYB_TETRA_1 ( feTriaP0, geoLinearTria, quadRuleTria1pt,
                                                      refcoor_HYB_TETRA_FACE_1, 0 );
static const CurrentFEManifold BdFE_RT0_HYB_TETRA_2 ( feTriaP0, geoLinearTria, quadRuleTria1pt,
                                                      refcoor_HYB_TETRA_FACE_2, 1 );
static const CurrentFEManifold BdFE_RT0_HYB_TETRA_3 ( feTriaP0, geoLinearTria, quadRuleTria1pt,
                                                      refcoor_HYB_TETRA_FACE_3, 2 );
static const CurrentFEManifold BdFE_RT0_HYB_TETRA_4 ( feTriaP0, geoLinearTria, quadRuleTria1pt,
                                                      refcoor_HYB_TETRA_FACE_4, 3 );

static const CurrentFEManifold* HybRT0TetraList[ NB_BDFE_RT0_HYB_TETRA ] =
{
    &BdFE_RT0_HYB_TETRA_1, &BdFE_RT0_HYB_TETRA_2,
    &BdFE_RT0_HYB_TETRA_3, &BdFE_RT0_HYB_TETRA_4
};

/*const RefHybridFE feTetraRT0Hyb (NB_BDFE_RT0_HYB_TETRA,HybRT0TetraList,"Hybrid RT0 elements on a tetrahedron",
    FE_RT0_HYB_TETRA_3D, TETRA, 0,0,1,0,4,3,
    refcoor_RT0HYB_TETRA,STANDARD_PATTERN);*/


static const CurrentFEManifold BdFE_RT0_HYB_TETRA_VdotN_1 ( feTriaP0, geoLinearTria, quadRuleTria1pt,
                                                            refcoor_HYB_TETRA_FACE_1, 0, 2. );
static const CurrentFEManifold BdFE_RT0_HYB_TETRA_VdotN_2 ( feTriaP0, geoLinearTria, quadRuleTria1pt,
                                                            refcoor_HYB_TETRA_FACE_2, 1, 2. );
static const CurrentFEManifold BdFE_RT0_HYB_TETRA_VdotN_3 ( feTriaP0, geoLinearTria, quadRuleTria1pt,
                                                            refcoor_HYB_TETRA_FACE_3, 2, 2. / std::sqrt ( 3. ) );
static const CurrentFEManifold BdFE_RT0_HYB_TETRA_VdotN_4 ( feTriaP0, geoLinearTria, quadRuleTria1pt,
                                                            refcoor_HYB_TETRA_FACE_4, 3, 2. );

static const CurrentFEManifold* HybRT0TetraVdotNList[ NB_BDFE_RT0_HYB_TETRA ] =
{
    &BdFE_RT0_HYB_TETRA_VdotN_1, &BdFE_RT0_HYB_TETRA_VdotN_2,
    &BdFE_RT0_HYB_TETRA_VdotN_3, &BdFE_RT0_HYB_TETRA_VdotN_4
};

const ReferenceFEHybrid feTetraRT0Hyb ( "Hybrid RT0 elements on a tetrahedron", FE_RT0_HYB_TETRA_3D, TETRA,
                                        0, 0, 1, 0, 4, 3, NB_BDFE_RT0_HYB_TETRA, HybRT0TetraList,
                                        refcoor_RT0HYB_TETRA, STANDARD_PATTERN );

const ReferenceFEHybrid feTetraRT0VdotNHyb ( "Hybrid RT0 elements on a tetrahedron", FE_RT0_HYB_TETRA_3D, TETRA,
                                             0, 0, 1, 0, 4, 3, NB_BDFE_RT0_HYB_TETRA, HybRT0TetraVdotNList,
                                             refcoor_RT0HYB_TETRA, STANDARD_PATTERN );


}
