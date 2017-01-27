/*
 * Code Author: Miroslav Stoyanov
 *
 * Copyright (C) 2015  Miroslav Stoyanov
 *
 * This file is part of
 * Toolkit for Adaptive Stochastic Modeling And Non-Intrusive Approximation
 *              a.k.a. TASMANIAN
 *
 * TASMANIAN is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * TASMANIAN is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with TASMANIAN.  If not, see <http://www.gnu.org/licenses/>
 *
 */

#ifndef __TASMANIAN_SPARSE_GRID_ENUMERATES_HPP
#define __TASMANIAN_SPARSE_GRID_ENUMERATES_HPP

#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <string.h>

namespace TasGrid{

using std::cout; // for debugging purposes
using std::endl; // for debugging purposes

using std::cerr; // for error messages

enum TypeIndexRelation{ // internal for IndexSets, unambiguous set comparison
        type_abeforeb, type_bbeforea, type_asameb
};

enum TypeDepth{
        type_none,
        type_level, type_curved,
        type_iptotal, type_ipcurved,
        type_qptotal, type_qpcurved,
        type_hyperbolic, type_iphyperbolic, type_qphyperbolic,
        type_tensor, type_iptensor, type_qptensor
};

enum TypeOneDRule{ // list of one-d rules
        rule_none, // encountering is an indication of error
        rule_clenshawcurtis,
        rule_clenshawcurtis0, // assumes zero boundary
        rule_chebyshev,
        rule_chebyshevodd,
        rule_gausslegendre,
        rule_gausslegendreodd,
        rule_gausspatterson,
        rule_leja,
        rule_lejaodd,
        rule_rleja,
        rule_rlejadouble2,
        rule_rlejadouble4,
        rule_rlejaodd,
        rule_rlejashifted,
        rule_rlejashiftedeven,
        rule_maxlebesgue,
        rule_maxlebesgueodd,
        rule_minlebesgue,
        rule_minlebesgueodd,
        rule_mindelta,
        rule_mindeltaodd,
        rule_gausschebyshev1,
        rule_gausschebyshev1odd,
        rule_gausschebyshev2,
        rule_gausschebyshev2odd,
        rule_fejer2,
        rule_gaussgegenbauer,
        rule_gaussgegenbauerodd,
        rule_gaussjacobi,
        rule_gaussjacobiodd,
        rule_gausslaguerre,
        rule_gausslaguerreodd,
        rule_gausshermite,
        rule_gausshermiteodd,
        rule_customtabulated,
        // Piece-Wise rules
        rule_localp,
        rule_localp0,
        rule_semilocalp,
        // Wavelet rules
        rule_wavelet
};

enum TypeRefinement{
        refine_classic, refine_parents_first, refine_direction_selective, refine_fds /* FDS = parents_first + direction_selective */
};

/////////////////////////////////////////////////////////
//                                                     //
//     ---  Moved from tsgHardcodedConstants.hpp  ---  //
//                                                     //
/////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////
//
//  On the purpose of this section:
//
//
//  Hardcoding constants is a bad coding practice and should be avoided.
//  On the other hand, numerical algorithms depend on many small tweaking parameters.
//  For example, this code computes the nodes and abscissas of Gauss-Legendre rule on the fly,
//  this is done with an iterative eigenvalue decomposition method that needs stopping criteria,
//  we can add another user specified parameter to the corresponding "makeGlobalGrid()" rule,
//  however, this is an additional tweaking variable that the user has to understand and control.
//  The goal of our code is to provide a most seamless experience to the user,
//  one should think about the desired properties of a quadrature rule or interpolant and not
//  about convergence of an iterative scheme (especially one hidden from top view)
//  Therefore, the tolerance for such convergence criteria needs to be set inside the code to
//  a "reasonable" value, which is a value that would work well for the overwhelming majority
//  of use cases.
//
//  On the other hand, every application is different and the usage of the code will change
//  over time. It is unreasonable to believe that one single value would work for absolutely
//  everyone. Instead of "hiding" hardcoded constant throughout the code, all such constants
//  will be exposed here so the used can make adjustments in compile time.
//
//  Long story short, do not adjust those variables unless you have a good reason.
//
///////////////////////////////////////////////////////////////////////////////////////////////

//  NUM_TOL is used in many places:
// - as a stopping criteria for various iterative schemes (e.g., finding leja points)
// - drop criteria for eigenvalue solver related to Gauss rules
// - comparison between nodes to detect repeated points in non-nested rules (e.g., all odd Chebyshev rules include zero)
// - determining sparse matrix pattern, entries smaller than NUM_TOL will be ignored (for wavelet grids)
// - drop criteria in estimating anisotropic coefficients (refinement or just the coefficients) surpluses or Legendre coefficients below 10^3 times NUM_TOL will be ignored
#define TSG_NUM_TOL 1.E-12 // to move in hardcoded constants



// this defines the maximum number of secant method iterations to be used for finding Leja, Lebesgue, and Delta points
// this is a safeguard criteria to prevent "hanging" in a loop
#define TSG_MAX_SECANT_ITERATIONS 1000 // to move in hardcoded constants

}

#endif
