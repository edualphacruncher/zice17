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

#ifndef __TASMANIAN_SPARSE_GRID_HPP
#define __TASMANIAN_SPARSE_GRID_HPP

#include "tsgEnumerates.hpp"

#include "tsgGridGlobal.hpp"
#include "tsgGridSequence.hpp"
#include "tsgGridLocalPolynomial.hpp"
#include "tsgGridWavelet.hpp"

#include <iomanip> // only needed for printStats()

namespace TasGrid{

class TasmanianSparseGrid{
public:
        TasmanianSparseGrid();
        ~TasmanianSparseGrid();

        static const char* getVersion();
        static const char* getLicense();

        void write( const char *filename ) const;
        bool read( const char *filename );

        void write( std::ofstream &ofs ) const;
        bool read( std::ifstream &ifs );

        void makeGlobalGrid( int dimensions, int outputs, int depth, TypeDepth type, TypeOneDRule rule, const int *anisotropic_weights = 0, double alpha = 0.0, double beta = 0.0, const char* custom_filename = 0 );
        void makeSequenceGrid( int dimensions, int outputs, int depth, TypeDepth type, TypeOneDRule rule, const int *anisotropic_weights = 0 );
        void makeLocalPolynomialGrid( int dimensions, int outputs, int depth, int order = 1, TypeOneDRule rule = rule_localp );
        void makeWaveletGrid( int dimensions, int outputs, int depth, int order = 1 );

        void updateGlobalGrid( int depth, TypeDepth type, const int *anisotropic_weights = 0 );
        void updateSequenceGrid( int depth, TypeDepth type, const int *anisotropic_weights = 0 );

        double getAlpha() const;
        double getBeta() const;
        int getOrder() const;

        int getNumDimensions() const;
        int getNumOutputs() const;
        TypeOneDRule getRule() const;
        const char* getCustomRuleDescription() const; // used only for Global Grids with rule_customtabulated

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

        bool isGlobal() const;
        bool isSequence() const;
        bool isLocalPolynomial() const;
        bool isWavelet() const;

        void setDomainTransform( const double a[], const double b[] ); // set the ranges of the box
        bool isSetDomainTransfrom() const;
        void clearDomainTransform();
        void getDomainTransform( double a[], double b[] ) const;

        void setAnisotropicRefinement( TypeDepth type, int min_growth = 1, int output = 0 );
        int* estimateAnisotropicCoefficients( TypeDepth type, int output = 0 );
        void setSurplusRefinement( double tolerance, int output = 0 );
        void setSurplusRefinement( double tolerance, TypeRefinement criteria );
        void clearRefinement();

        int* getGlobalPolynomialSpace( bool interpolation, int &num_indexes ) const;

        void printStats() const;

protected:
        void clear();

        void mapCanonicalToTransformed( int num_dimensions, int num_points, TypeOneDRule rule, double x[] ) const;
        void mapTransformedToCanonical( int num_dimensions, TypeOneDRule rule, double x[] ) const;
        double getQuadratureScale( int num_dimensions, TypeOneDRule rule ) const;

private:
        BaseCanonicalGrid *base;

        GridGlobal *global;
        GridSequence *sequence;
        GridLocalPolynomial *pwpoly;
        GridWavelet *wavelet;

        double *domain_transform_a, *domain_transform_b;
};

}

#endif
