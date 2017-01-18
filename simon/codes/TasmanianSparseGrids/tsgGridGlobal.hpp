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

#ifndef __TASMANIAN_SPARSE_GRID_GLOBAL_HPP
#define __TASMANIAN_SPARSE_GRID_GLOBAL_HPP

#include "tsgEnumerates.hpp"
#include "tsgIndexSets.hpp"
#include "tsgCoreOneDimensional.hpp"
#include "tsgIndexManipulator.hpp"
#include "tsgLinearSolvers.hpp"
#include "tsgCacheLagrange.hpp"
#include "tsgOneDimensionalWrapper.hpp"
#include "tsgGridCore.hpp"

namespace TasGrid{

class GridGlobal : public BaseCanonicalGrid{
public:
        GridGlobal();
        ~GridGlobal();

        void write( std::ofstream &ofs ) const;
        void read( std::ifstream &ifs );

        void makeGrid( int cnum_dimensions, int cnum_outputs, int depth, TypeDepth type, TypeOneDRule crule, const int *anisotropic_weights = 0, double calpha = 0.0, double cbeta = 0.0, const char* custom_filename = 0 );

        void setTensors( IndexSet* &tset, int cnum_outputs, TypeOneDRule crule, double calpha, double cbeta );

        void updateGrid(  int depth, TypeDepth type, const int *anisotropic_weights = 0 );

        int getNumDimensions() const;
        int getNumOutputs() const;
        TypeOneDRule getRule() const;
        const char* getCustomRuleDescription() const; // returns the description of the custom rule (only if rule_customtabulated is used)

        double getAlpha() const;
        double getBeta() const;

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

        int* estimateAnisotropicCoefficients( TypeDepth type, int output ) const;

        void setAnisotropicRefinement( TypeDepth type, int min_growth = 1, int output = 0 );
        void setSurplusRefinement( double tolerance, int output = 0 );
        void clearRefinement();

        int* getPolynomialSpace( bool interpolation, int &n ) const;

protected:
        void reset( bool includeCustom );

        double* computeSurpluses( int output, bool normalize ) const; // only for sequence rules, select the output to compute the surpluses

        static double legendre( int n, double x );
        static double* multi_legendre( int n, double x );

private:
        int num_dimensions, num_outputs;
        TypeOneDRule rule;
        double alpha, beta;

        OneDimensionalWrapper *wrapper;

        IndexSet *tensors;
        IndexSet *active_tensors;
        int *active_w;
        IndexSet *points;
        IndexSet *needed;

        int **tensor_refs;

        int *max_levels; // for evaluation purposes, counts the maximum level in each direction (only counts tensors)

        StorageSet *values;

        IndexSet *updated_tensors;
        IndexSet *updated_active_tensors;
        int *updated_active_w;

        CustomTabulated *custom;
};

}

#endif
