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

#ifndef __TSG_CACHE_LAGRANGE_CPP
#define __TSG_CACHE_LAGRANGE_CPP

#include "tsgCacheLagrange.hpp"

namespace TasGrid{

CacheLagrange::CacheLagrange( int num_dimensions, const int max_levels[], const OneDimensionalWrapper *crule, const double x[] ) : rule(crule){
        cache = new double*[num_dimensions];
        int full_cache = rule->getPointsCount( max_levels[0] + 1 );
        for( int j=1; j<num_dimensions; j++ ){
                full_cache += rule->getPointsCount( max_levels[j]+1 );
        }
        cache[0] = new double[ full_cache ];
        for( int j=0; j<num_dimensions-1; j++ ){
                cache[j+1] = &( cache[j][ rule->getPointsCount( max_levels[j]+1 ) ] );
        }

        for( int dim=0; dim<num_dimensions; dim++ ){
                for( int level=0; level <= max_levels[dim]; level++ ){
                        const double *nodes = rule->getNodes( level );
                        const double *coeff = rule->getCoefficients( level );
                        int num_points = rule->getNumPoints( level );

                        double *c = &( cache[dim][ rule->getPointsCount( level ) ] );
                        c[0] = 1.0;
                        for( int j=0; j<num_points-1; j++ ){
                                c[j+1] = ( x[dim] - nodes[j] ) * c[j];
                        }
                        double w = ( rule->getType() == rule_clenshawcurtis0 ) ? ( x[dim] - 1.0 ) * ( x[dim] + 1.0 ) : 1.0;
                        c[num_points-1] *= w * coeff[num_points-1];
                        for( int j=num_points-2; j>=0; j-- ){
                                w *= ( x[dim] - nodes[j+1] );
                                c[j] *= w * coeff[j];
                        }
                }
        }
}

CacheLagrange::~CacheLagrange(){
        if ( cache != 0 ){
                delete[] cache[0];
                delete[] cache;
        }
        rule = 0;
}

double CacheLagrange::getLagrange( int dimension, int level, int local ) const{
        return cache[dimension][ rule->getPointsCount( level ) + local ];
}

}

#endif
