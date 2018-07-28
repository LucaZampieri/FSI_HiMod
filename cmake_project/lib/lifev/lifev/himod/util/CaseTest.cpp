#include <lifev/himod/util/CaseTest.hpp>
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-variable"

namespace LifeV
{

/*
    If you change the coefficients be care, change them also in CaseTest.hpp!!!!!
*/

Real
CaseTest::g (const Real& /*t*/, const Real& /*x*/, const Real& y, const Real& z, const ID& /*i*/)
{
    Real Ly = 0.1;
    Real Lz = 0.1;
    Real Lx = 0.2;
    Real x = 0.0;

    return 1e7 * y * z * (Ly - y) * (Lz - z) * std::pow (x - Lx, 2) * std::exp (std::pow ( (Lx - x), 2) * y * z * 2);
}
Real
CaseTest::f (const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& /*i*/)
{
    Real Lx = 0.2;
    Real Ly = 0.1;
    Real Lz = 0.1;

    using std::exp;
    using std::pow;

    return     (
                   2e7 * y * pow (Lx - x, 2) * (Ly - y)
                   + 2e7 * z * pow (Lx - x, 2) * (Lz - z)
                   + 4e7 * y * y * z * pow (Lx - x, 4) * (Ly - y)
                   + 4e7 * y * z * z * pow (Lx - x, 4) * (Lz - z)
                   - 4e7 * y * y * pow (Lx - x, 4) * (Ly - y) * (Lz - z)
                   - 4e7 * z * z * pow (Lx - x, 4) * (Ly - y) * (Lz - z)
                   - 2e7 * y * z * (Ly - y) * (Lz - z)
                   - 4e7 * y * y * z * z * (Ly - y) * (Lz - z) * pow (2 * Lx - 2 * x, 2)
                   - 4e7 * y * z * z * z * pow (Lx - x, 6) * (Ly - y) * (Lz - z)
                   - 4e7 * y * y * y * z * pow (Lx - x, 6) * (Ly - y) * (Lz - z)
                   - 4e7 * y * y * z * z * pow (Lx - x, 2) * (Ly - y) * (Lz - z)
                   - 4e7 * y * y * y * z * z * z * pow (Lx - x, 2) * (Ly - y) * (Lz - z) * pow (2 * Lx - 2 * x, 2)
               ) * exp (2 * y * z * pow (Lx - x, 2) );

}
Real CaseTest::ues (const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& /*i*/)
{
    Real Lx = 0.2;
    Real Ly = 0.1;
    Real Lz = 0.1;

    return 1e7 * y * z * (Ly - y) * (Lz - z) * std::pow (x - Lx, 2) * std::exp (std::pow ( (Lx - x), 2) * y * z * 2);

}

Real CaseTest::g_adv (const Real& /*t*/, const Real& /*x*/, const Real& y, const Real& z, const ID& /*i*/)
{
    Real Ly = 0.1;
    Real Lz = 0.1;
    Real x = 0.0;
    Real Lx = 0.2;

    return 1e7 * y * z * (Ly - y) * (Lz - z) * std::pow (x - Lx, 2) * std::exp (std::pow ( (Lx - x), 2) * y * z * 2);
}

Real CaseTest::f_adv (const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& /*i*/)
{
    Real Lx = 0.2;
    Real Ly = 0.1;
    Real Lz = 0.1;

    using std::exp;
    using std::pow;

    return 1e7 * (
               2 * y * pow (Lx - x, 2) * (Ly - y) + 2 * z * pow (Lx - x, 2) * (Lz - z)
               + 4 * y * y * z * pow (Lx - x, 4) * (Ly - y) + 4 * y * z * z * pow (Lx - x, 4) * (Lz - z)
               - 4 * y * y * pow (Lx - x, 4) * (Ly - y) * (Lz - z) - 4 * z * z * pow (Lx - x, 4) * (Ly - y) * (Lz - z)
               - 2 * y * z * (Ly - y) * (Lz - z) - y * z * pow (Lx - x, 2) * (Ly - y)
               - y * z * pow (Lx - x, 2) * (Lz - z) + y * pow (Lx - x, 2) * (Ly - y) * (Lz - z)
               + z * pow (Lx - x, 2) * (Ly - y) * (Lz - z) - 5 * y * z * (Ly - y) * (Lz - z) * (2 * Lx - 2 * x)
               - 4 * y * y * z * z * (Ly - y) * (Lz - z) * pow (2 * Lx - 2 * x, 2)
               + 2 * y * z * z * pow (Lx - x, 4) * (Ly - y) * (Lz - z) + 2 * y * y * z * pow (Lx - x, 4) * (Ly - y) * (Lz - z)
               - 4 * y * z * z * z * pow (Lx - x, 6) * (Ly - y) * (Lz - z)
               - 4 * y * y * y * z * pow (Lx - x, 6) * (Ly - y) * (Lz - z)
               - 4 * y * y * z * z * pow (Lx - x, 2) * (Ly - y) * (Lz - z)
               - 4 * y * y * y * z * z * z * pow (Lx - x, 2) * (Ly - y) * (Lz - z) * pow (2 * Lx - 2 * x, 2)
               - 10 * y * y * z * z * pow (Lx - x, 2) * (Ly - y) * (Lz - z) * (2 * Lx - 2 * x)
           ) * exp (2 * y * z * pow (Lx - x, 2) );
}

Real CaseTest::ues_adv (const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& /*i*/)
{
    Real Lx = 0.2;
    Real Ly = 0.1;
    Real Lz = 0.1;

    return 1e7 * y * z * (Ly - y) * (Lz - z) * std::pow (x - Lx, 2) * std::exp (std::pow ( (Lx - x), 2) * y * z * 2);
}

Real CaseTest::g_adr (const Real& /*t*/, const Real& /*x*/, const Real& y, const Real& z, const ID& /*i*/)
{
    Real Ly = 0.1;
    Real Lz = 0.1;
    Real x = 0.0;
    Real Lx = 0.2;

    return 1e7 * y * z * (Ly - y) * (Lz - z) * std::pow (x - Lx, 2) * std::exp (std::pow ( (Lx - x), 2) * y * z * 2);
}



Real CaseTest::f_adr (const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& /*i*/)
{
    Real Lx = 0.2;
    Real Ly = 0.1;
    Real Lz = 0.1;

    using std::exp;
    using std::pow;


    return (
               2e7 * y * pow (Lx - x, 2) * (Ly - y)
               + 2e7 * z * pow (Lx - x, 2) * (Lz - z)
               + 4e7 * y * y * z * pow (Lx - x, 4) * (Ly - y)
               + 4e7 * y * z * z * pow (Lx - x, 4) * (Lz - z)
               - 4e7 * y * y * pow (Lx - x, 4) * (Ly - y) * (Lz - z)
               - 4e7 * z * z * pow (Lx - x, 4) * (Ly - y) * (Lz - z)
               - 2e7 * y * z * (Ly - y) * (Lz - z)
               - 1e7 * y * z * pow (Lx - x, 2) * (Ly - y)
               - 1e7 * y * z * pow (Lx - x, 2) * (Lz - z)
               + 1e7 * y * pow (Lx - x, 2) * (Ly - y) * (Lz - z)
               + 1e7 * z * pow (Lx - x, 2) * (Ly - y) * (Lz - z)
               - 5e7 * y * z * (Ly - y) * (Lz - z) * (2 * Lx - 2 * x)
               + 3e7 * y * z * pow (Lx - x, 2) * (Ly - y) * (Lz - z)
               - 4e7 * y * y * z * z * (Ly - y) * (Lz - z) * pow (2 * Lx - 2 * x, 2)
               + 2e7 * y * z * z * pow (Lx - x, 4) * (Ly - y) * (Lz - z)
               + 2e7 * y * y * z * pow (Lx - x, 4) * (Ly - y) * (Lz - z)
               - 4e7 * y * z * z * z * pow (Lx - x, 6) * (Ly - y) * (Lz - z)
               - 4e7 * y * y * y * z * pow (Lx - x, 6) * (Ly - y) * (Lz - z)
               - 4e7 * y * y * z * z * pow (Lx - x, 2) * (Ly - y) * (Lz - z)
               - 4e7 * y * y * y * z * z * z * pow (Lx - x, 2) * (Ly - y) * (Lz - z) * pow (2 * Lx - 2 * x, 2)
               - 1e8 * y * y * z * z * pow (Lx - x, 2) * (Ly - y) * (Lz - z) * (2 * Lx - 2 * x)
           ) * exp (2 * y * z * pow (Lx - x, 2) );
}

Real CaseTest::ues_adr (const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& /*i*/)
{
    Real Lx = 0.2;
    Real Ly = 0.1;
    Real Lz = 0.1;

    return 1e7 * y * z * (Ly - y) * (Lz - z) * std::pow (x - Lx, 2) * std::exp (std::pow ( (Lx - x), 2) * y * z * 2);
}


Real CaseTest::g_drdr (const Real& /*t*/, const Real& /*x*/, const Real& y, const Real& z, const ID& /*i*/)
{
    Real Ly  = 0.1;
    Real Lz  = 0.1;
    Real x   = 0.0;
    Real Lx  = 0.1;
    Real chi = 3.345;
    Real mu  = 1.0;

    using std::exp;
    using std::pow;

    return 1e5 * pow ( (Lx - x), 2 ) * z * ( Lz - z ) *
           exp (
               70 * y * y / ( x * z + 1 ) - 140 * pow (y, 3) / ( 3 * Ly * ( x * z + 1 ) ) - chi * pow ( Ly / 2. - y , 2 ) / ( Ly * mu )
           );
}

Real CaseTest::f_drdr (const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& /*i*/)
{
    Real Ly  = 0.1;
    Real Lz  = 0.1;
    Real Lx  = 0.1;
    Real chi = 3.345;

    using std::exp;
    using std::pow;

    return (
               - 2e5 * z * (Lz - z) + 2e5 * pow (Lx - x, 2)
               - 2e5 * z * pow (Lx - x, 2) * ( (70 * x * y * y) / pow (x * z + 1, 2)
                                               - (140 * x * y * y * y) / (3 * Ly * pow (x * z + 1, 2) ) )
               + 2e5 * pow (Lx - x, 2) * (Lz - z) * ( (70 * x * y * y) / pow (x * z + 1, 2)
                                                      - (140 * x * y * y * y) / (3 * Ly * pow (x * z + 1, 2) ) )
               - 1e5 * z * pow (Lx - x, 2) * (Lz - z) * pow ( (70 * x * y * y) / pow (x * z + 1, 2)
                                                              - (140 * x * y * y * y) / (3 * Ly * pow (x * z + 1, 2) ), 2) - 1e5 * z * pow (Lx - x, 2) * (Lz - z) * pow ( (70 * y * y * z) / pow (x * z + 1, 2)
                                                                      - (140 * y * y * y * z) / (3 * Ly * pow (x * z + 1, 2) ), 2 ) - 1e5 * z * ( (140 * x * x * y * y) / pow (x * z + 1, 3)
                                                                              - (280 * x * x * y * y * y) / (3 * Ly * pow (x * z + 1, 3) ) ) * pow (Lx - x, 2) * (Lz - z)
               - 1e5 * z * ( (140 * y * y * z * z) / pow (x * z + 1, 3)
                             - (280 * y * y * y * z * z) / (3 * Ly * pow (x * z + 1, 3) ) ) * pow (Lx - x, 2) * (Lz - z)
               - 1e5 * z * pow (Lx - x, 2) * (Lz - z) * pow ( (140 * y) / (x * z + 1) + (chi * (Ly - 2 * y) ) / Ly
                                                              - (140 * y * y) / (Ly * (x * z + 1) ), 2) - 2e5 * z * (Lz - z) * ( (70 * y * y * z) / pow (x * z + 1, 2)
                                                                      - (140 * y * y * y * z) / (3 * Ly * pow (x * z + 1, 2) ) ) * (2 * Lx - 2 * x)
               + 1e5 * z * pow (Lx - x, 2) * (Lz - z) * ( (2 * chi) / Ly - 140 / (x * z + 1)
                                                          + (280 * y) / (Ly * (x * z + 1) ) )
           ) * exp ( (70 * y * y) / (x * z + 1) - (chi * pow (Ly / 2 - y, 2) ) / Ly - (140 * y * y * y) / (3 * Ly * (x * z + 1) ) );
}

Real CaseTest::ues_drdr (const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& /*i*/)
{
    Real Ly  = 0.1;
    Real Lz  = 0.1;
    Real Lx  = 0.1;
    Real chi = 3.345;
    Real mu  = 1.0;

    using std::exp;
    using std::pow;

    return 1e5 * pow ( (Lx - x), 2 ) * z * ( Lz - z ) *
           exp (
               70 * y * y / ( x * z + 1 ) - 140 * pow (y, 3) / ( 3 * Ly * ( x * z + 1 ) ) - chi * pow ( Ly / 2. - y , 2 ) / ( Ly * mu )
           );
}

Real CaseTest::g_rrrr (const Real& /*t*/, const Real& /*x*/, const Real& y, const Real& z, const ID& /*i*/)
{
    Real Ly  = 0.1;
    Real Lz  = 0.1;
    Real x   = 0.0;
    Real Lx  = 0.1;
    Real chi = 4.456;
    Real Mu  = 1.0;

    using std::exp;
    using std::pow;

    return  1e5 * exp ( (70 * y * y) /  (x + 1) - (140 * y * y * y)       / (3 * Ly * (x + 1) ) - (chi * pow (Ly / 2 - y, 2) ) / (Ly * Mu) ) * exp ( (70 * z * z) / (x + 1) - (140 * z * z * z) / (3 * Lz * (x + 1) ) - (chi * pow (Lz / 2 - z, 2) ) / (Lz * Mu) ) * pow (Lx - x, 2);
}

Real CaseTest::f_rrrr (const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& /*i*/)
{
    Real Ly  = 0.1;
    Real Lz  = 0.1;
    Real Lx  = 0.1;
    Real chi = 4.456;

    using std::exp;
    using std::pow;

    return 100000 * exp ( (70 * y * y) / (x + 1) - (140 * y * y * y) / (3 * Ly * (x + 1) ) - (chi * pow (Ly / 2 - y, 2) ) / Ly) * exp ( (70 * z * z) / (x + 1) - (140 * z * z * z) / (3 * Lz * (x + 1) ) - (chi * pow (Lz / 2 - z, 2) ) / Lz) * pow (Lx - x, 2) * ( (2 * chi) / Ly - 140 / (x + 1) + (280 * y) / (Ly * (x + 1) ) ) - 200000 * exp ( (70 * y * y) / (x + 1) - (140 * y * y * y) / (3 * Ly * (x + 1) ) - (chi * pow (Ly / 2 - y, 2) ) / Ly) * exp ( (70 * z * z) / (x + 1) - (140 * z * z * z) / (3 * Lz * (x + 1) ) - (chi * pow (Lz / 2 - z, 2) ) / Lz) + 100000 * exp ( (70 * y * y) / (x + 1) - (140 * y * y * y) / (3 * Ly * (x + 1) ) - (chi * pow (Ly / 2 - y, 2) ) / Ly) * exp ( (70 * z * z) / (x + 1) - (140 * z * z * z) / (3 * Lz * (x + 1) ) - (chi * pow (Lz / 2 - z, 2) ) / Lz) * pow (Lx - x, 2) * ( (2 * chi) / Lz - 140 / (x + 1) + (280 * z) / (Lz * (x + 1) ) ) - 100000 * exp ( (70 * y * y) / (x + 1) - (140 * y * y * y) / (3 * Ly * (x + 1) ) - (chi * pow (Ly / 2 - y, 2) ) / Ly) * exp ( (70 * z * z) / (x + 1) - (140 * z * z * z) / (3 * Lz * (x + 1) ) - (chi * pow (Lz / 2 - z, 2) ) / Lz) * pow (Lx - x, 2) * pow ( (70 * y * y) / pow (x + 1, 2) - (140 * y * y * y) / (3 * Ly * pow (x + 1, 2) ), 2) - 100000 * exp ( (70 * y * y) / (x + 1) - (140 * y * y * y) / (3 * Ly * (x + 1) ) - (chi * pow (Ly / 2 - y, 2) ) / Ly) * exp ( (70 * z * z) / (x + 1) - (140 * z * z * z) / (3 * Lz * (x + 1) ) - (chi * pow (Lz / 2 - z, 2) ) / Lz) * pow (Lx - x, 2) * pow ( (70 * z * z) / pow (x + 1, 2) - (140 * z * z * z) / (3 * Lz * pow (x + 1, 2) ), 2) - 100000 * exp ( (70 * y * y) / (x + 1) - (140 * y * y * y) / (3 * Ly * (x + 1) ) - (chi * pow (Ly / 2 - y, 2) ) / Ly) * exp ( (70 * z * z) / (x + 1) - (140 * z * z * z) / (3 * Lz * (x + 1) ) - (chi * pow (Lz / 2 - z, 2) ) / Lz) * pow (Lx - x, 2) * pow ( (140 * y) / (x + 1) - (140 * y * y) / (Ly * (x + 1) ) + (chi * (Ly - 2 * y) ) / Ly, 2) - 100000 * exp ( (70 * y * y) / (x + 1) - (140 * y * y * y) / (3 * Ly * (x + 1) ) - (chi * pow (Ly / 2 - y, 2) ) / Ly) * exp ( (70 * z * z) / (x + 1) - (140 * z * z * z) / (3 * Lz * (x + 1) ) - (chi * pow (Lz / 2 - z, 2) ) / Lz) * pow (Lx - x, 2) * pow ( (140 * z) / (x + 1) - (140 * z * z) / (Lz * (x + 1) ) + (chi * (Lz - 2 * z) ) / Lz, 2) - 200000 * exp ( (70 * y * y) / (x + 1) - (140 * y * y * y) / (3 * Ly * (x + 1) ) - (chi * pow (Ly / 2 - y, 2) ) / Ly) * exp ( (70 * z * z) / (x + 1) - (140 * z * z * z) / (3 * Lz * (x + 1) ) - (chi * pow (Lz / 2 - z, 2) ) / Lz) * (2 * Lx - 2 * x) * ( (70 * y * y) / pow (x + 1, 2) - (140 * y * y * y) / (3 * Ly * pow (x + 1, 2) ) ) - 100000 * exp ( (70 * y * y) / (x + 1) - (140 * y * y * y) / (3 * Ly * (x + 1) ) - (chi * pow (Ly / 2 - y, 2) ) / Ly) * exp ( (70 * z * z) / (x + 1) - (140 * z * z * z) / (3 * Lz * (x + 1) ) - (chi * pow (Lz / 2 - z, 2) ) / Lz) * pow (Lx - x, 2) * ( (140 * y * y) / pow (x + 1, 3) - (280 * y * y * y) / (3 * Ly * pow (x + 1, 3) ) ) - 200000 * exp ( (70 * y * y) / (x + 1) - (140 * y * y * y) / (3 * Ly * (x + 1) ) - (chi * pow (Ly / 2 - y, 2) ) / Ly) * exp ( (70 * z * z) / (x + 1) - (140 * z * z * z) / (3 * Lz * (x + 1) ) - (chi * pow (Lz / 2 - z, 2) ) / Lz) * (2 * Lx - 2 * x) * ( (70 * z * z) / pow (x + 1, 2) - (140 * z * z * z) / (3 * Lz * pow (x + 1, 2) ) ) - 100000 * exp ( (70 * y * y) / (x + 1) - (140 * y * y * y) / (3 * Ly * (x + 1) ) - (chi * pow (Ly / 2 - y, 2) ) / Ly) * exp ( (70 * z * z) / (x + 1) - (140 * z * z * z) / (3 * Lz * (x + 1) ) - (chi * pow (Lz / 2 - z, 2) ) / Lz) * pow (Lx - x, 2) * ( (140 * z * z) / pow (x + 1, 3) - (280 * z * z * z) / (3 * Lz * pow (x + 1, 3) ) ) - 200000 * exp ( (70 * y * y) / (x + 1) - (140 * y * y * y) / (3 * Ly * (x + 1) ) - (chi * pow (Ly / 2 - y, 2) ) / Ly) * exp ( (70 * z * z) / (x + 1) - (140 * z * z * z) / (3 * Lz * (x + 1) ) - (chi * pow (Lz / 2 - z, 2) ) / Lz) * pow (Lx - x, 2) * ( (70 * y * y) / pow (x + 1, 2) - (140 * y * y * y) / (3 * Ly * pow (x + 1, 2) ) ) * ( (70 * z * z) / pow (x + 1, 2) - (140 * z * z * z) / (3 * Lz * pow (x + 1, 2) ) );

}


Real CaseTest::ues_rrrr (const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& /*i*/)
{

    Real Ly  = 0.1;
    Real Lz  = 0.1;
    Real Lx  = 0.1;
    Real chi = 4.456;
    Real Mu  = 1.0;

    using std::exp;
    using std::pow;

    return  1e5 * exp ( (70 * y * y) /  (x + 1) - (140 * y * y * y)       / (3 * Ly * (x + 1) ) - (chi * pow (Ly / 2 - y, 2) ) / (Ly * Mu) ) * exp ( (70 * z * z) / (x + 1) - (140 * z * z * z) / (3 * Lz * (x + 1) ) - (chi * pow (Lz / 2 - z, 2) ) / (Lz * Mu) ) * pow (Lx - x, 2);

}



Real CaseTest::g_rr2d (const Real& /*t*/, const Real& /*x*/, const Real& y, const Real& /*z*/, const ID& /*i*/)
{
    Real x = 0.0;
    return y + pow (y - 1, 2) - 3 * y * (y - 1) + 4 * y * y * pow (x - 1, 2) * pow (y - 1, 2) * (8 * x * x * y + 8 * x * y * y + 0.75) + 1;
}
Real CaseTest::f_rr2d (const Real& /*t*/, const Real& x, const Real& y, const Real& /*z*/, const ID& /*i*/)
{
    return 2 * pow (y - 1, 2) - 6 * y - 6 * y * (y - 1) - 8 * y * y * pow (x - 1, 2) * (8 * x * x * y + 8 * x * y * y + 3 / 4) - 8 * y * y * pow (y - 1, 2) * (8 * x * x * y + 8 * x * y * y + 3 / 4) - 8 * pow (x - 1, 2) * pow (y - 1, 2) * (8 * x * x * y + 8 * x * y * y + 3 / 4) - 64 * y * y * y * pow (x - 1, 2) * pow (y - 1, 2) - 16 * y * (2 * y - 2) * pow (x - 1, 2) * (8 * x * x * y + 8 * x * y * y + 3 / 4) + 16 * y * pow (x - 1, 2) * pow (y - 1, 2) * (8 * x * x * y + 8 * x * y * y + 3 / 4) - 64 * x * y * y * pow (x - 1, 2) * pow (y - 1, 2) - 16 * y * pow (x - 1, 2) * pow (y - 1, 2) * (8 * x * x + 16 * y * x) + 80 * y * y * (2 * x - 2) * pow (y - 1, 2) * (8 * x * x * y + 8 * x * y * y + 3 / 4) + 8 * y * y * (2 * y - 2) * pow (x - 1, 2) * (8 * x * x * y + 8 * x * y * y + 3 / 4) + 8 * y * y * pow (x - 1, 2) * pow (y - 1, 2) * (8 * x * x * y + 8 * x * y * y + 3 / 4) - 8 * y * y * (2 * y - 2) * pow (x - 1, 2) * (8 * x * x + 16 * y * x) + 8 * y * y * pow (x - 1, 2) * pow (y - 1, 2) * (8 * x * x + 16 * y * x) - 8 * y * y * (2 * x - 2) * pow (y - 1, 2) * (8 * y * y + 16 * x * y) + 80 * y * y * pow (x - 1, 2) * pow (y - 1, 2) * (8 * y * y + 16 * x * y) + 10;

}
Real CaseTest::ues_rr2d (const Real& /*t*/, const Real& x, const Real& y, const Real& /*z*/, const ID& /*i*/)
{
    return y + pow (y - 1, 2) - 3 * y * (y - 1) + 4 * y * y * pow (x - 1, 2) * pow (y - 1, 2) * (8 * x * x * y + 8 * x * y * y + 0.75) + 1;
}




Real CaseTest::g_2drr (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& z, const ID& /*i*/)
{
    Real x = 0.0;
    return z + pow (z - 1, 2) - 3 * z * (z - 1) + 4 * z * z * pow (x - 1, 2) * pow (z - 1, 2) * (8 * x * x * z + 8 * x * z * z + 0.75) + 1;
}
Real CaseTest::f_2drr (const Real& /*t*/, const Real& x, const Real& /*y*/, const Real& z, const ID& /*i*/)
{
    return 2 * pow (z - 1, 2) - 6 * z - 6 * z * (z - 1) - 8 * z * z * pow (x - 1, 2) * (8 * x * x * z + 8 * x * z * z + 3 / 4) - 8 * z * z * pow (z - 1, 2) * (8 * x * x * z + 8 * x * z * z + 3 / 4) - 8 * pow (x - 1, 2) * pow (z - 1, 2) * (8 * x * x * z + 8 * x * z * z + 3 / 4) - 64 * z * z * z * pow (x - 1, 2) * pow (z - 1, 2) - 16 * z * (2 * z - 2) * pow (x - 1, 2) * (8 * x * x * z + 8 * x * z * z + 3 / 4) + 16 * z * pow (x - 1, 2) * pow (z - 1, 2) * (8 * x * x * z + 8 * x * z * z + 3 / 4) - 64 * x * z * z * pow (x - 1, 2) * pow (z - 1, 2) - 16 * z * pow (x - 1, 2) * pow (z - 1, 2) * (8 * x * x + 16 * z * x) + 80 * z * z * (2 * x - 2) * pow (z - 1, 2) * (8 * x * x * z + 8 * x * z * z + 3 / 4) + 8 * z * z * (2 * z - 2) * pow (x - 1, 2) * (8 * x * x * z + 8 * x * z * z + 3 / 4) + 8 * z * z * pow (x - 1, 2) * pow (z - 1, 2) * (8 * x * x * z + 8 * x * z * z + 3 / 4) - 8 * z * z * (2 * z - 2) * pow (x - 1, 2) * (8 * x * x + 16 * z * x) + 8 * z * z * pow (x - 1, 2) * pow (z - 1, 2) * (8 * x * x + 16 * z * x) - 8 * z * z * (2 * x - 2) * pow (z - 1, 2) * (8 * z * z + 16 * x * z) + 80 * z * z * pow (x - 1, 2) * pow (z - 1, 2) * (8 * z * z + 16 * x * z) + 10;

}
Real CaseTest::ues_2drr (const Real& /*t*/, const Real& x, const Real& /*y*/, const Real& z, const ID& /*i*/)
{
    return z + pow (z - 1, 2) - 3 * z * (z - 1) + 4 * z * z * pow (x - 1, 2) * pow (z - 1, 2) * (8 * x * x * z + 8 * x * z * z + 0.75) + 1;
}


Real CaseTest::g_dd2d (const Real& /*t*/, const Real& /*x*/, const Real& y, const Real& /*z*/, const ID& /*i*/)
{
    Real x = 0.0;
    return -y * exp (sin (20 * y * y * y * pow (x - 1, 2) * pow (y - 1, 2) ) ) * (y - 1);
}
Real CaseTest::f_dd2d (const Real& /*t*/, const Real& x, const Real& y, const Real& /*z*/, const ID& /*i*/)
{
    return 2 * exp (sin (20 * y * y * y * pow (x - 1, 2) * pow (y - 1, 2) ) ) - 2 * y * exp (sin (20 * y * y * y * pow (x - 1, 2) * pow (y - 1, 2) ) ) - 2 * exp (sin (20 * y * y * y * pow (x - 1, 2) * pow (y - 1, 2) ) ) * (y - 1) - 2 * y * exp (sin (20 * y * y * y * pow (x - 1, 2) * pow (y - 1, 2) ) ) * (y - 1) + 2 * y * cos (20 * y * y * y * pow (x - 1, 2) * pow (y - 1, 2) ) * exp (sin (20 * y * y * y * pow (x - 1, 2) * pow (y - 1, 2) ) ) * (20 * y * y * y * (2 * y - 2) * pow (x - 1, 2) + 60 * y * y * pow (x - 1, 2) * pow (y - 1, 2) ) + 2 * cos (20 * y * y * y * pow (x - 1, 2) * pow (y - 1, 2) ) * exp (sin (20 * y * y * y * pow (x - 1, 2) * pow (y - 1, 2) ) ) * (20 * y * y * y * (2 * y - 2) * pow (x - 1, 2) + 60 * y * y * pow (x - 1, 2) * pow (y - 1, 2) ) * (y - 1) + 40 * y * y * y * y * cos (20 * y * y * y * pow (x - 1, 2) * pow (y - 1, 2) ) * exp (sin (20 * y * y * y * pow (x - 1, 2) * pow (y - 1, 2) ) ) * pow (y - 1, 3) - 400 * y * y * y * y * cos (20 * y * y * y * pow (x - 1, 2) * pow (y - 1, 2) ) * exp (sin (20 * y * y * y * pow (x - 1, 2) * pow (y - 1, 2) ) ) * (2 * x - 2) * pow (y - 1, 3) - y * sin (20 * y * y * y * pow (x - 1, 2) * pow (y - 1, 2) ) * exp (sin (20 * y * y * y * pow (x - 1, 2) * pow (y - 1, 2) ) ) * pow (20 * y * y * y * (2 * y - 2) * pow (x - 1, 2) + 60 * y * y * pow (x - 1, 2) * pow (y - 1, 2), 2) * (y - 1) - 400 * y * y * y * y * y * y * y * sin (20 * y * y * y * pow (x - 1, 2) * pow (y - 1, 2) ) * exp (sin (20 * y * y * y * pow (x - 1, 2) * pow (y - 1, 2) ) ) * pow (2 * x - 2, 2) * pow (y - 1, 5) + y * pow (cos (20 * y * y * y * pow (x - 1, 2) * pow (y - 1, 2) ), 2) * exp (sin (20 * y * y * y * pow (x - 1, 2) * pow (y - 1, 2) ) ) * pow (20 * y * y * y * (2 * y - 2) * pow (x - 1, 2) + 60 * y * y * pow (x - 1, 2) * pow (y - 1, 2), 2) * (y - 1) - 2 * y * cos (20 * y * y * y * pow (x - 1, 2) * pow (y - 1, 2) ) * exp (sin (20 * y * y * y * pow (x - 1, 2) * pow (y - 1, 2) ) ) * (20 * y * y * y * (2 * y - 2) * pow (x - 1, 2) + 60 * y * y * pow (x - 1, 2) * pow (y - 1, 2) ) * (y - 1) + 400 * y * y * y * y * y * y * y * pow (cos (20 * y * y * y * pow (x - 1, 2) * pow (y - 1, 2) ), 2) * exp (sin (20 * y * y * y * pow (x - 1, 2) * pow (y - 1, 2) ) ) * pow (2 * x - 2, 2) * pow (y - 1, 5) + y * cos (20 * y * y * y * pow (x - 1, 2) * pow (y - 1, 2) ) * exp (sin (20 * y * y * y * pow (x - 1, 2) * pow (y - 1, 2) ) ) * (y - 1) * (40 * y * y * y * pow (x - 1, 2) + 120 * y * pow (x - 1, 2) * pow (y - 1, 2) + 120 * y * y * (2 * y - 2) * pow (x - 1, 2) );
}
Real CaseTest::ues_dd2d (const Real& /*t*/, const Real& x, const Real& y, const Real& /*z*/, const ID& /*i*/)
{
    return -y * exp (sin (20 * y * y * y * pow (x - 1, 2) * pow (y - 1, 2) ) ) * (y - 1);
}


Real CaseTest::g_2ddd (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& z, const ID& /*i*/)
{
    Real x = 0.0;
    return -z * exp (sin (20 * z * z * z * pow (x - 1, 2) * pow (z - 1, 2) ) ) * (z - 1);
}
Real CaseTest::f_2ddd (const Real& /*t*/, const Real& x, const Real& /*y*/, const Real& z, const ID& /*i*/)
{
    return 2 * exp (sin (20 * z * z * z * pow (x - 1, 2) * pow (z - 1, 2) ) ) - 2 * z * exp (sin (20 * z * z * z * pow (x - 1, 2) * pow (z - 1, 2) ) ) - 2 * exp (sin (20 * z * z * z * pow (x - 1, 2) * pow (z - 1, 2) ) ) * (z - 1) - 2 * z * exp (sin (20 * z * z * z * pow (x - 1, 2) * pow (z - 1, 2) ) ) * (z - 1) + 2 * z * cos (20 * z * z * z * pow (x - 1, 2) * pow (z - 1, 2) ) * exp (sin (20 * z * z * z * pow (x - 1, 2) * pow (z - 1, 2) ) ) * (20 * z * z * z * (2 * z - 2) * pow (x - 1, 2) + 60 * z * z * pow (x - 1, 2) * pow (z - 1, 2) ) + 2 * cos (20 * z * z * z * pow (x - 1, 2) * pow (z - 1, 2) ) * exp (sin (20 * z * z * z * pow (x - 1, 2) * pow (z - 1, 2) ) ) * (20 * z * z * z * (2 * z - 2) * pow (x - 1, 2) + 60 * z * z * pow (x - 1, 2) * pow (z - 1, 2) ) * (z - 1) + 40 * z * z * z * z * cos (20 * z * z * z * pow (x - 1, 2) * pow (z - 1, 2) ) * exp (sin (20 * z * z * z * pow (x - 1, 2) * pow (z - 1, 2) ) ) * pow (z - 1, 3) - 400 * z * z * z * z * cos (20 * z * z * z * pow (x - 1, 2) * pow (z - 1, 2) ) * exp (sin (20 * z * z * z * pow (x - 1, 2) * pow (z - 1, 2) ) ) * (2 * x - 2) * pow (z - 1, 3) - z * sin (20 * z * z * z * pow (x - 1, 2) * pow (z - 1, 2) ) * exp (sin (20 * z * z * z * pow (x - 1, 2) * pow (z - 1, 2) ) ) * pow (20 * z * z * z * (2 * z - 2) * pow (x - 1, 2) + 60 * z * z * pow (x - 1, 2) * pow (z - 1, 2), 2) * (z - 1) - 400 * z * z * z * z * z * z * z * sin (20 * z * z * z * pow (x - 1, 2) * pow (z - 1, 2) ) * exp (sin (20 * z * z * z * pow (x - 1, 2) * pow (z - 1, 2) ) ) * pow (2 * x - 2, 2) * pow (z - 1, 5) + z * pow (cos (20 * z * z * z * pow (x - 1, 2) * pow (z - 1, 2) ), 2) * exp (sin (20 * z * z * z * pow (x - 1, 2) * pow (z - 1, 2) ) ) * pow (20 * z * z * z * (2 * z - 2) * pow (x - 1, 2) + 60 * z * z * pow (x - 1, 2) * pow (z - 1, 2), 2) * (z - 1) - 2 * z * cos (20 * z * z * z * pow (x - 1, 2) * pow (z - 1, 2) ) * exp (sin (20 * z * z * z * pow (x - 1, 2) * pow (z - 1, 2) ) ) * (20 * z * z * z * (2 * z - 2) * pow (x - 1, 2) + 60 * z * z * pow (x - 1, 2) * pow (z - 1, 2) ) * (z - 1) + 400 * z * z * z * z * z * z * z * pow (cos (20 * z * z * z * pow (x - 1, 2) * pow (z - 1, 2) ), 2) * exp (sin (20 * z * z * z * pow (x - 1, 2) * pow (z - 1, 2) ) ) * pow (2 * x - 2, 2) * pow (z - 1, 5) + z * cos (20 * z * z * z * pow (x - 1, 2) * pow (z - 1, 2) ) * exp (sin (20 * z * z * z * pow (x - 1, 2) * pow (z - 1, 2) ) ) * (z - 1) * (40 * z * z * z * pow (x - 1, 2) + 120 * z * pow (x - 1, 2) * pow (z - 1, 2) + 120 * z * z * (2 * z - 2) * pow (x - 1, 2) );
}
Real CaseTest::ues_2ddd (const Real& /*t*/, const Real& x, const Real& /*y*/, const Real& z, const ID& /*i*/)
{
    return -z * exp (sin (20 * z * z * z * pow (x - 1, 2) * pow (z - 1, 2) ) ) * (z - 1);
}

//-------------------------------------------------------------------------------------------------------------
//---------------------------------------- SHOW CASE TESTS --------------------------------------------------

Real CaseTest::ues_camini (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return NAN;
}

Real CaseTest::f_camini ( const Real& /* t */, const Real& x, const Real& y, const Real&  z  , const ID& /* i */ )
{

    Real Lx = 0.2;
    Real Ly = 0.1;
    Real Lz = 0.1;

    Real frac1X = 1. / 8.;
    Real frac1Y = 1. / 8.;
    Real frac1Z = 1. / 8.;

    Real frac2X = 2. / 8.;
    Real frac2Y = 7. / 8.;
    Real frac2Z = 7. / 8.;

    Real thick = 900;

    using std::exp;
    using std::pow;

    return  1e2 * ( exp ( -thick * ( pow (x - Lx * frac1X, 2) + pow (y - Ly * frac1Y, 2) + pow (z - Lz * frac1Z, 2) ) )
                    + exp ( -thick * ( pow (x - Lx * frac2X, 2) + pow (y - Ly * frac2Y, 2) + pow (z - Lz * frac2Z, 2) ) ) );

}

Real CaseTest::g_camini (const Real& /*t*/, const Real& /*x*/, const Real& y, const Real& z, const ID& /*i*/)
{
    return 0 * y * z;
}


Real CaseTest::f_circ1 (const Real& /*t*/, const Real& x, const Real& rho, const Real& theta, const ID& /*i*/)
{
    Real Lx = 1;
    Real Rho = 0.1;
    Real mu = 1.;
    Real betax = 3;
    Real betarho = 1.;
    Real betatheta = 1.;
    Real sigma = 0;
    Real Theta = 2 * M_PI;

    using std::exp;
    using std::pow;

    return mu * 2 * ( rho * rho - Rho * Rho + 2 * ( Lx - x ) * ( Lx - x ) )
		    -2 * betarho * rho * ( Lx - x ) * ( Lx - x )
    		- 2 * betax * ( - rho * rho + Rho * Rho ) * ( Lx - x )
    		+ sigma * ( - rho * rho + Rho * Rho) * ( Lx - x ) * ( Lx - x );

}

Real CaseTest::frho_circ1 (const Real& /*t*/, const Real& x, const Real& rho, const Real& theta, const ID& /*i*/)
{
	return f( 0., x, rho, theta, 0 ) * rho;
}

// The exact solution
Real CaseTest::ues_circ1 (const Real& /*t*/, const Real& x, const Real& rho, const Real& theta, const ID& /*i*/)
{
    Real Lx = 1;
    Real Rho = 0.1;
    Real Theta = 2 * M_PI;
    
    return ( Rho * Rho - rho * rho ) * ( Lx - x ) * ( Lx - x );


}

// The dirichlet inflow data
Real CaseTest::g_circ1 (const Real& /*t*/, const Real& /*x*/, const Real& rho, const Real& theta, const ID& /*i*/)
{
    Real Rho = 0.1;
    Real Theta = 2 * M_PI;
    Real x = 0.0;
    Real Lx = 1;
    return ues( 0., 0., rho, theta, 0 );
}



Real CaseTest::ues_stent (const Real& /*t*/, const Real& /*x*/, const Real& /*rho*/, const Real& /*theta*/, const ID& /*i*/)
{
    return NAN;
}

Real CaseTest::f_stent ( const Real& /* t */, const Real& x, const Real& rho, const Real&  theta  , const ID& /* i */ )
{

    Real Lx = 2;
    Real Rho = 0.1;
    Real Theta = 0.1;

	return ( rho == Rho && abs( x - Lx/3 ) < 0.2 );
}

Real CaseTest::g_stent (const Real& /*t*/, const Real& /*x*/, const Real& rho, const Real& theta, const ID& /*i*/)
{
	Real Rho = 0.1;
	
    return Rho * Rho - rho * rho;
}

}
