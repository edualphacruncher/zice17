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

#ifndef __TASMANIAN_SPARSE_GRID_GLOBAL_NESTED_HPP
#define __TASMANIAN_SPARSE_GRID_GLOBAL_NESTED_HPP

#include "tsgEnumerates.hpp"
#include "tsgIndexSets.hpp"
#include "tsgCoreOneDimensional.hpp"
#include "tsgIndexManipulator.hpp"
#include "tsgLinearSolvers.hpp"
#include "tsgCacheLagrange.hpp"
#include "tsgOneDimensionalWrapper.hpp"
#include "tsgGridCore.hpp"

namespace TasGrid{

class GridSequence : public BaseCanonicalGrid{
public:
        GridSequence();
        ~GridSequence();

        void write( std::ofstream &ofs ) const;
        void read( std::ifstream &ifs );

        void makeGrid( int cnum_dimensions, int cnum_outputs, int depth, TypeDepth type, TypeOneDRule crule, const int *anisotropic_weights = 0 );
        void setPoints( IndexSet* &pset, int cnum_outputs, TypeOneDRule crule );

        void updateGrid( int depth, TypeDepth type, const int *anisotropic_weights = 0 );
        void updateGrid( IndexSet* &update );

        int getNumDimensions() const;
        int getNumOutputs() const;
        TypeOneDRule getRule() const;

        int getNumLoaded() const;
        int getNumNeeded() const;
        int getNumPoints() const; // returns the number of loaded points unless no points are loaded, then returns the number of needed points

        double* getLoadedPoints() const;
        double* getNeededPoints() const;
        double* getPoints() const; // returns the loaded points unless no points are loaded, then returns the needed points

        double* getQuadratureWeights() const;
        double* getInterpolationWeights( const double x[] ) const;

        void loadNeededPoints( const double *vals );

        void evaluate( const double x[], double y[] ) const;
        void integrate( double q[] ) const;

        int* estimateAnisotropicCoefficients( TypeDepth type ) const;

        void setAnisotropicRefinement( TypeDepth type, int min_growth = 1 );
        void setSurplusRefinement( double tolerance );
        void clearRefinement();

        int* getPolynomialSpace( bool interpolation, int &n ) const;

protected:
        void reset();

        void prepareSequence( int n );
        double** cacheBasisValues( const double x[] ) const;
        double* cacheBasisIntegrals() const;

        void recomputeSurpluses();
        void applyTransformationTransposed( double weights[] ) const;

        double evalBasis( const int f[], const int p[] ) const; // evaluate function corresponding to f at p

private:
        int num_dimensions, num_outputs;
        TypeOneDRule rule;

        IndexSet *points;
        IndexSet *needed;
        int *parents; // NOTE: this is needed only for computing surpluses, maybe there is no need to store it

        double *surpluses;
        double *nodes;
        double *coeff;

        StorageSet *values;

        int *max_levels;
};

}

#endif
