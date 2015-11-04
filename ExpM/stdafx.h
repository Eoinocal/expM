// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#include "targetver.h"

#include <stdio.h>
#include <tchar.h>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/detail/vector_assign.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/operations.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/optional.hpp>

#include <boost/numeric/odeint.hpp>

#include <limits>
#include <random>
#include <vector>
#include <numeric>
#include <fstream>
#include <random>

#pragma warning (push)
#pragma warning (disable : 4003)
#include <vexcl/vexcl.hpp>
#include <vexcl/external/viennacl.hpp>
#pragma warning (push)

// Must be set if you want to use ViennaCL algorithms on ublas objects
#define VIENNACL_WITH_UBLAS 1

// ViennaCL headers
#include <viennacl/scalar.hpp>
#include <viennacl/vector.hpp>
#include <viennacl/matrix.hpp>
#include <viennacl/linalg/prod.hpp>
#include <viennacl/tools/random.hpp>
#include <viennacl/tools/timer.hpp>

namespace ub = boost::numeric::ublas;
namespace ode = boost::numeric::odeint;
namespace vcl = viennacl;
