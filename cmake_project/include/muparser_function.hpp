#ifndef HAVE_MUPARSER_FUNCTION_H
#define HAVE_MUPARSER_FUNCTION_H

#include <muParser.h>
#include <memory>
#include <string>

#include <lifev/core/LifeV.hpp>

namespace LifeV
{

class muparser_function
{
private :

  Real tVar;
  Real xVar;
  Real rVar;
  Real aVar;
  Real iVar;

  mu::Parser p;

public :

  muparser_function (const muparser_function & m)
    : p (m.p)
  { p.DefineVar ("t", &tVar);
    p.DefineVar ("x", &xVar);
    p.DefineVar ("r", &rVar);
    p.DefineVar ("a", &aVar);
    p.DefineVar ("i", &iVar);};

  muparser_function (const std::string & s)
  {
    try
      {
        p.DefineVar ("t", &tVar);
        p.DefineVar ("x", &xVar);
        p.DefineVar ("r", &rVar);
        p.DefineVar ("a", &aVar);
        p.DefineVar ("i", &iVar);

        p.SetExpr (s);
      }
    catch (mu::Parser::exception_type &e)
      { std::cerr << e.GetMsg () << std::endl; }
  };

  Real operator() (Real x, Real r, Real a, Real t, Real i)
  {
    Real y;

    tVar = t;
    xVar = x;
    rVar = r;
    aVar = a;
    iVar = i;

    try
      { y = p.Eval (); }
    catch (mu::Parser::exception_type &e)
      { std::cerr << e.GetMsg () << std::endl; }
    return y;
  };

};

}

#endif
