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

#ifndef __TASMANIAN_SPARSE_GRID_LPOLY_HPP
#define __TASMANIAN_SPARSE_GRID_LPOLY_HPP

#include "tsgEnumerates.hpp"
#include "tsgIndexSets.hpp"
#include "tsgIndexManipulator.hpp"
#include "tsgGridCore.hpp"
#include "tsgRuleLocalPolynomial.hpp"

namespace TasGrid{

class GridLocalPolynomial : public BaseCanonicalGrid{
public:
        GridLocalPolynomial();
        ~GridLocalPolynomial();

        void write( std::ofstream &ofs ) const;
        void read( std::ifstream &ifs );

        void makeGrid( int cnum_dimensions, int cnum_outputs, int depth, int corder, TypeOneDRule crule );

        int getNumDimensions() const;
        int getNumOutputs() const;
        TypeOneDRule getRule() const;
        int getOrder() const;

        int getNumLoaded() const;
        int getNumNeeded() const;
        int getNumPoints() const;

        double* getLoadedPoints() const;
        double* getNeededPoints() const;
        double* getPoints() const;

        double* getQuadratureWeights() const;
        double* getInterpolationWeights( const double x[] ) const;

        void loadNeededPoints( const double *vals );

        void evaluate( const double x[], double y[] ) const;
        void integrate( double q[] ) const;

        void setSurplusRefinement( double tolerance, TypeRefinement criteria );
        void clearRefinement();

protected:
        void reset();

        void buildTree();

        void recomputeSurpluses();

        double evalBasisRaw( const int point[], const double x[] ) const;
        double evalBasisSupported( const int point[], const double x[], bool &isSupported ) const;

        double* getBasisIntegrals() const;

        double* getNormalization() const;

        int* buildUpdateMap( double tolerance, TypeRefinement criteria ) const;

        bool addParent( const int point[], int direction, GranulatedIndexSet *destination, IndexSet *exclude ) const;
        void addChild( const int point[], int direction, GranulatedIndexSet *destination, IndexSet *exclude )const;

private:
        int num_dimensions, num_outputs, order, top_level;

        double *surpluses;

        IndexSet *points;
        IndexSet *needed;

        StorageSet *values;
        int *parents;

        // three for evaluation
        int num_roots, *roots;
        int *pntr, *indx;

        BaseRuleLocalPolynomial *rule;
        RuleLocalPolynomial rpoly;
        RuleSemiLocalPolynomial rsemipoly;
        RuleLocalPolynomialZero rpoly0;
};

}

#endif
