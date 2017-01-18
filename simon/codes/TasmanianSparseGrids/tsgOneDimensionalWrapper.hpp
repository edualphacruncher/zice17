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

#ifndef __TSG_ONE_DIMENSIONAL_WRAPPER_HPP
#define __TSG_ONE_DIMENSIONAL_WRAPPER_HPP

#include "tsgEnumerates.hpp"
#include "tsgCoreOneDimensional.hpp"
#include "tsgSequenceOptimizer.hpp"
#include "tsgHardCodedTabulatedRules.hpp"

namespace TasGrid{

class OneDimensionalWrapper{
public:
        OneDimensionalWrapper( const OneDimensionalMeta *meta, int max_level, TypeOneDRule crule, double alpha = 0.0, double beta = 0.0 );
        ~OneDimensionalWrapper();

        int getNumPoints( int level ) const;
        int getPointIndex( int level, int j ) const;

        double getNode( int j ) const;
        double getWeight( int level, int j ) const;

        const double* getNodes( int level ) const;
        const double* getCoefficients( int level ) const;

        int getPointsCount( int level ) const;

        TypeOneDRule getType() const;

private:
        bool isNonNested;
        int num_levels;
        TypeOneDRule rule;

        int *num_points; // give the number of points per level
        int *pntr; // gives the cumulative offset of each level
        int *indx; // gives a list of the points associated with each level

        double *weights; // contains the weight associated with each level
        double *nodes; // contains all nodes for each level
        double *unique; // contains the x-coordinate of each sample point

        double *coeff; // the coefficients of the Lagrange
};

}

#endif
