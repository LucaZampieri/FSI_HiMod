#ifndef HAVE_MUPARSER_TIME_FUNCTION_H
#define HAVE_MUPARSER_TIME_FUNCTION_H

#include <muParser.h>
#include <memory>
#include <string>

#include <lifev/core/LifeV.hpp>

namespace LifeV
{

class muparser_timeFunction
{
private :

  Real var;
  mu::Parser p;

public :

  muparser_timeFunction (const muparser_timeFunction & m)
    : p (m.p)
  { p.DefineVar ("t", &var); };

  muparser_timeFunction (const std::string & s)
  {
    try
      {
        p.DefineVar ("t", &var);
        p.SetExpr (s);
      }
    catch (mu::Parser::exception_type &e)
      { std::cerr << e.GetMsg () << std::endl; }
  };

  Real operator() (Real t)
  {
    Real y;
    var = t;
    try
      { y = p.Eval (); }
    catch (mu::Parser::exception_type &e)
      { std::cerr << e.GetMsg () << std::endl; }
    return (y);
  };

};

}

#endif
