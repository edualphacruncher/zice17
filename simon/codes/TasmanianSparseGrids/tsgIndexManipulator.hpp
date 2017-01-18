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

#ifndef __TSG_INDEX_MANIPULATOR_HPP
#define __TSG_INDEX_MANIPULATOR_HPP

#include "tsgEnumerates.hpp"
#include "tsgIndexSets.hpp"
#include "tsgCoreOneDimensional.hpp"
#include "tsgOneDimensionalWrapper.hpp"
#include "tsgRuleLocalPolynomial.hpp"

namespace TasGrid{

class IndexManipulator{
public:
        IndexManipulator( int cnum_dimensions, const CustomTabulated* custom = 0 );
        ~IndexManipulator();

        int getIndexWeight( const int index[], TypeDepth type, const int weights[], TypeOneDRule rule ) const;

        IndexSet* selectTensors( int offset, TypeDepth type, const int *anisotropic_weights, TypeOneDRule rule ) const;
        IndexSet* getLowerCompletion( const IndexSet *set ) const;

        int* computeLevels( const IndexSet* set ) const; // returns a vector of its corresponding to the sum of entries of set

        int* makeTensorWeights( const IndexSet* set ) const;
        IndexSet* nonzeroSubset( const IndexSet* set, const int weights[] ) const; // returns a subset corresponding to the non-zeros weights

        UnsortedIndexSet* tensorGenericPoints( const int levels[], const OneDimensionalWrapper *rule ) const;
        IndexSet* generateGenericPoints( const IndexSet *tensors, const OneDimensionalWrapper *rule ) const;
        int* referenceGenericPoints( const int levels[], const OneDimensionalWrapper *rule, const IndexSet *points ) const;

        int* computeDAGdown( const IndexSet *set ) const;
        int* computeDAGup( const IndexSet *set ) const;

        IndexSet* tensorNestedPoints( const int levels[], const OneDimensionalWrapper *rule ) const;
        IndexSet* generateNestedPoints( const IndexSet *tensors, const OneDimensionalWrapper *rule ) const;
        int* referenceNestedPoints( const int levels[], const OneDimensionalWrapper *rule, const IndexSet *points ) const;

        IndexSet* getPolynomialSpace( const IndexSet *tensors, TypeOneDRule rule, bool iexact ) const;

        int getMinChildLevel( const IndexSet *set, TypeDepth type, const int weights[], TypeOneDRule rule );

        IndexSet* selectFlaggedChildren( const IndexSet *set, const bool flagged[] ) const;

        void getMaxLevels( const IndexSet *set, int max_levels[], int &total_max ) const;

        // use by Local Grids
        UnsortedIndexSet* getToalDegreeDeltas( int level ) const;
        IndexSet* generatePointsFromDeltas( const UnsortedIndexSet* deltas, const BaseRuleLocalPolynomial *rule ) const;

        int* computeDAGupLocal( const IndexSet *set, const BaseRuleLocalPolynomial *rule ) const;

protected:
        int getLevel( const int index[], const int weights[] ) const;
        int getCurved( const int index[], const int weights[] ) const;
        int getIPTotal( const int index[], const int weights[], TypeOneDRule rule ) const;
        int getIPCurved( const int index[], const int weights[], TypeOneDRule rule ) const;
        int getQPTotal( const int index[], const int weights[], TypeOneDRule rule ) const;
        int getQPCurved( const int index[], const int weights[], TypeOneDRule rule ) const;
        int getHyperbolic( const int index[], const int weights[] ) const;
        int getIPHyperbolic( const int index[], const int weights[], TypeOneDRule rule ) const;
        int getQPHyperbolic( const int index[], const int weights[], TypeOneDRule rule ) const;

        int* getProperWeights( TypeDepth type, const int *anisotropic_weights ) const;

private:
        int num_dimensions;

        OneDimensionalMeta *meta;
};

}

#endif
