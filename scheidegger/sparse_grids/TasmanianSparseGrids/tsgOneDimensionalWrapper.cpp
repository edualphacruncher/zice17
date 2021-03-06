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

#ifndef __TSG_ONE_DIMENSIONAL_WRAPPER_CPP
#define __TSG_ONE_DIMENSIONAL_WRAPPER_CPP

#include "tsgOneDimensionalWrapper.hpp"

namespace TasGrid{

OneDimensionalWrapper::OneDimensionalWrapper( const OneDimensionalMeta *meta, int max_level, TypeOneDRule crule, double alpha, double beta) :
        num_levels(max_level+1), rule(crule), indx(0), nodes(0)
{
        // find the points per level and the cumulative pointers
        isNonNested = OneDimensionalMeta::isNonNested( rule );
        num_points = new int[num_levels];
        pntr = new int[num_levels+1]; pntr[0] = 0;
        for( int l=0; l<num_levels; l++ ){
                num_points[l] = meta->getNumPoints( l, rule );
                pntr[l+1] = pntr[l] + num_points[l];
        }
        int num_total = pntr[num_levels];

        weights = new double[num_total];
        coeff = new double[num_total];

        OneDimensionalNodes core;

        if ( isNonNested ){
                indx    = new int[num_total];
                nodes   = new double[num_total];

                int num_unique = 0;
                unique = new double[num_total];

                double *x = 0, *w = 0;

                if ( rule == rule_customtabulated ){
                        if ( num_levels > meta->getCustom()->getNumLevels() ){
                                cerr << "ERROR: custom-tabulated rule needed with levels " << num_levels << ", but only " << meta->getCustom()->getNumLevels() << " are provided." << endl;
                                exit(1);
                        }
                }

                for( int l=0; l<num_levels; l++ ){
                        int n = num_points[l];
                        if ( (rule == rule_chebyshev) || (rule == rule_chebyshevodd) ){
                                core.getChebyshev( n, w, x );
                        }else if ( (rule == rule_gausslegendre) || (rule == rule_gausslegendreodd) ){
                                core.getGaussLegendre( n, w, x );
                        }else if ( (rule == rule_gausschebyshev1) || (rule == rule_gausschebyshev1odd) ){
                                core.getGaussChebyshev1( n, w, x );
                        }else if ( (rule == rule_gausschebyshev2) || (rule == rule_gausschebyshev2odd) ){
                                core.getGaussChebyshev2( n, w, x );
                        }else if ( (rule == rule_gaussgegenbauer) || (rule == rule_gaussgegenbauerodd) ){
                                core.getGaussJacobi( n, w, x, alpha, alpha );
                        }else if ( (rule == rule_gausshermite) || (rule == rule_gausshermiteodd) ){
                                core.getGaussHermite( n, w, x, alpha );
                        }else if ( (rule == rule_gaussjacobi) || (rule == rule_gaussjacobiodd) ){
                                core.getGaussJacobi( n, w, x, alpha, beta );
                        }else if ( (rule == rule_gausslaguerre) || (rule == rule_gausslaguerreodd) ){
                                core.getGaussLaguerre( n, w, x, alpha );
                        }else if ( rule == rule_customtabulated ){
                                meta->getCustom()->getWeightsNodes( l, w, x );
                        }

                        for( int i=0; i<n; i++ ){
                                weights[ pntr[l] + i ] = w[i];
                                nodes[pntr[l] + i] = x[i];

                                int point = -1;
                                for( int j=0; j<num_unique; j++ ){
                                        if ( fabs( x[i] - unique[j] ) < TSG_NUM_TOL ){
                                                point = j;
                                                break;
                                        }
                                }
                                if ( point == - 1){ // new point found
                                        unique[num_unique] = x[i];
                                        point = num_unique++;
                                }
                                indx[pntr[l] + i] = point;
                        }

                }

                double *t = unique;
                unique = new double[num_unique];
                std::copy( t, t + num_unique, unique );

                delete[] t;
                delete[] x;
                delete[] w;

                // make coefficients
                for( int l=0; l<num_levels; l++ ){
                        int n = num_points[l];
                        x = &( nodes[ pntr[l] ] );
                        w = &( coeff[ pntr[l] ] );
                        for( int i=0; i<n; i++ ){
                                w[i] = 1.0;
                                for( int j=0; j<i; j++ ){
                                        w[i] /= ( x[i] - x[j] );
                                }
                                for( int j=i+1; j<n; j++ ){
                                        w[i] /= ( x[i] - x[j] );
                                }
                        }
                }
        }else{
                TableGaussPatterson *gp;
                if ( rule == rule_clenshawcurtis ){
                        unique = core.getClenshawCurtisNodes( max_level );
                }else if ( rule == rule_clenshawcurtis0 ){
                        unique = core.getClenshawCurtisNodesZero( max_level );
                }else if ( rule == rule_fejer2 ){
                        unique = core.getFejer2Nodes( max_level );
                }else if ( rule == rule_gausspatterson ){
                        gp = new TableGaussPatterson();
                        if ( max_level > gp->getMaxLeve() ){
                                cerr << "ERROR: gauss-patterson rule needed with level " << max_level << ", but only " << gp->getMaxLeve() << " are hardcoded." << endl;
                                exit(1);
                        }
                        unique = gp->getNodes( max_level );

                        // move here to avoid a warning
                        for( int l=0; l<num_levels; l++ ){
                                for( int i=0; i<num_points[l]; i++ ){
                                        weights[ pntr[l] + i ] = gp->getWeight( l, i );
                                }
                        }
                        delete gp;
                }else if ( rule == rule_rleja ){
                        unique = core.getRLeja( meta->getNumPoints(max_level,rule) );
                }else if ( (rule == rule_rlejaodd) || (rule == rule_rlejadouble2) || (rule == rule_rlejadouble4) ){
                        unique = core.getRLejaCentered( meta->getNumPoints(max_level,rule) );
                }else if ( (rule == rule_rlejashifted) || (rule == rule_rlejashiftedeven) ){
                        unique = core.getRLejaShifted( meta->getNumPoints(max_level,rule) );
                }else if ( (rule == rule_leja) || (rule == rule_lejaodd) ){
                        GreedySequences greedy;
                        unique = greedy.getLejaNodes( meta->getNumPoints(max_level, rule) );
                }else if ( (rule == rule_maxlebesgue) || (rule == rule_maxlebesgueodd) ){
                        GreedySequences greedy;
                        unique = greedy.getMaxLebesgueNodes( meta->getNumPoints(max_level, rule) );
                }else if ( (rule == rule_minlebesgue) || (rule == rule_minlebesgueodd) ){
                        GreedySequences greedy;
                        unique = greedy.getMinLebesgueNodes( meta->getNumPoints(max_level, rule) );
                }else if ( (rule == rule_mindelta) || (rule == rule_mindeltaodd) ){
                        GreedySequences greedy;
                        unique = greedy.getMinDeltaNodes( meta->getNumPoints(max_level, rule) );
                }

                for( int l=0; l<num_levels; l++ ){
                        int n = num_points[l];
                        double *c = &( coeff[ pntr[l] ] );
                        for( int i=0; i<n; i++ ){
                                c[i] = 1.0;
                                for( int j=0; j<i; j++ ){
                                        c[i] /= ( unique[i] - unique[j] );
                                }
                                for( int j=i+1; j<n; j++ ){
                                        c[i] /= ( unique[i] - unique[j] );
                                }
                                if ( rule == rule_clenshawcurtis0 ){
                                        c[i] /= ( unique[i] - 1.0 ) * ( unique[i] + 1.0 );
                                }
                        }
                }

                if ( rule == rule_clenshawcurtis ){
                        for( int l=0; l<num_levels; l++ ){
                                for( int i=0; i<num_points[l]; i++ ){
                                        weights[ pntr[l] + i ] = core.getClenshawCurtisWeight( l, i );
                                }
                        }
                }else if ( rule == rule_clenshawcurtis0 ){
                        for( int l=0; l<num_levels; l++ ){
                                for( int i=0; i<num_points[l]; i++ ){
                                        weights[ pntr[l] + i ] = core.getClenshawCurtisWeightZero( l, i );
                                }
                        }
                }else if ( rule == rule_fejer2 ){
                        for( int l=0; l<num_levels; l++ ){
                                for( int i=0; i<num_points[l]; i++ ){
                                        weights[ pntr[l] + i ] = core.getFejer2Weight( l, i );
                                }
                        }
                }else /*if ( rule == rule_gausspatterson ){
                        // weights are normally computed here, but I put the GP weights above
                }else*/ if ( (rule == rule_rleja) || (rule == rule_rlejaodd) || (rule == rule_rlejadouble2) || (rule == rule_rlejadouble4) || (rule == rule_leja) || (rule == rule_lejaodd) ||
                            (rule == rule_rlejashifted) || (rule == rule_rlejashiftedeven) ||
                           (rule == rule_maxlebesgue) || (rule == rule_maxlebesgueodd) || (rule == rule_minlebesgue) || (rule == rule_minlebesgueodd) || (rule == rule_mindelta) || (rule == rule_mindeltaodd) ){
                        int n = 1 + num_points[max_level] / 2; // number of Gauss-Legendre points needed to integrate the basis functions
                        double *lag_x = 0, *lag_w = 0;
                        core.getGaussLegendre( n, lag_w, lag_x );
                        std::fill( weights, weights+num_total, 0.0 );
                        for( int l=0; l<num_levels; l++ ){
                                int npl = num_points[l];
                                double *v = new double[ npl ];
                                double *c = &( coeff[ pntr[l] ] );
                                double *w = &( weights[ pntr[l] ] );
                                for( int i=0; i<n; i++ ){
                                        v[0] = 1.0;
                                        for( int j=0; j<npl-1; j++ ){
                                                v[j+1] = ( lag_x[i] - unique[j] ) * v[j];
                                        }
                                        v[npl-1] *= c[npl-1];
                                        double s = 1.0;
                                        for( int j=npl-2; j>=0; j-- ){
                                                s *= ( lag_x[i] - unique[j+1] );
                                                v[j] *= s * c[j];
                                        }
                                        for( int j=0; j<npl; j++ ){
                                                w[j] += lag_w[i] * v[j];
                                        }
                                }
                                delete[] v;
                        }
                        delete[] lag_w;
                        delete[] lag_x;
                }

        }
}

OneDimensionalWrapper::~OneDimensionalWrapper(){
        delete[] num_points;
        delete[] pntr;
        delete[] indx;
        delete[] weights;
        delete[] nodes;
        delete[] unique;
        delete[] coeff;
}

int OneDimensionalWrapper::getNumPoints( int level ) const{  return num_points[level];  }
int OneDimensionalWrapper::getPointIndex( int level, int j ) const{  return indx[ pntr[level] + j ];  }

double OneDimensionalWrapper::getNode( int j ) const{  return unique[j];  }
double OneDimensionalWrapper::getWeight( int level, int j ) const{  return weights[ pntr[level] + j ];  }

const double* OneDimensionalWrapper::getNodes( int level ) const{  return ( isNonNested ) ? &( nodes[pntr[level]] ) : unique;  }
const double* OneDimensionalWrapper::getCoefficients( int level ) const{  return &( coeff[pntr[level]] );  }

int OneDimensionalWrapper::getPointsCount( int level ) const{  return pntr[level]; }

TypeOneDRule OneDimensionalWrapper::getType() const{ return rule; }

}

#endif
