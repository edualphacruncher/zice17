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

#ifndef __TASMANIAN_SPARSE_GRID_SEQUENCE_OPTIMIZER_HPP
#define __TASMANIAN_SPARSE_GRID_SEQUENCE_OPTIMIZER_HPP

#include "tsgEnumerates.hpp"

namespace TasGrid{

class GreedySequences{
public:
        GreedySequences();
        ~GreedySequences();

        double* getLejaNodes( int n ) const;
        double* getMaxLebesgueNodes( int n ) const;
        double* getMinLebesgueNodes( int n ) const;
        double* getMinDeltaNodes( int n ) const;

        double findLebesgueConstant( int n, const double nodes[] ) const;

        int getNumMinLebesgueStored() const;
        double getMinLebesgueStored( int i ) const;
        int getNumMinDeltaStored() const;
        double getMinDeltaStored( int i ) const;
};

class Functional{
public:
        Functional();
        ~Functional();

        virtual double getValue( double x ) const;

        virtual bool hasDerivative() const;
        virtual double getDiff( double x ) const;

        virtual int getNumIntervals() const;
        virtual double* getIntervals() const;

        static double* makeCoeff( int num_nodes, const double nodes[] ); // may not be the best place to put this
        static double* evalLag( int num_nodes, const double nodes[], const double coeff[], double x );
        static double basisDx( int num_nodes, const double nodes[], const double coeff[], int inode, double x );
        static double* sortIntervals( int num_nodes, const double nodes[] );
        static int nodeCompar( const void * a, const void * b );
};

class Residual : public Functional{
public:
        Residual( int cnum_nodes, const double cnodes[] );
        ~Residual();

        double getValue( double x ) const;

        bool hasDerivative() const;
        double getDiff( double x ) const;

        int getNumIntervals() const;
        double* getIntervals() const;

private:
        int num_nodes;
        double *nodes;
};

class MaxLebesgue : public Functional{
public:
        MaxLebesgue( int cnum_nodes, const double cnodes[], double new_node );
        MaxLebesgue( int cnum_nodes, const double cnodes[] );
        ~MaxLebesgue();

        double getValue( double x ) const;

        bool hasDerivative() const;
        double getDiff( double x ) const;

        int getNumIntervals() const;
        double* getIntervals() const;

private:
        int num_nodes;
        double *nodes;
        double *coeff;
};

class MinLebesgue : public Functional{
public:
        MinLebesgue( int cnum_nodes, const double cnodes[] );
        ~MinLebesgue();

        double getValue( double x ) const;
        bool hasDerivative() const;

        int getNumIntervals() const;
        double* getIntervals() const;

private:
        int num_nodes;
        double *nodes;
};

class MaxDelta : public Functional{
public:
        MaxDelta( int cnum_nodes, const double cnodes[], double new_node );
        MaxDelta( int cnum_nodes, const double cnodes[] );
        ~MaxDelta();

        double getValue( double x ) const;

        bool hasDerivative() const;
        double getDiff( double x ) const;

        int getNumIntervals() const;
        double* getIntervals() const;

private:
        int num_nodes;
        double *nodes;
        double *coeff_plus, *coeff_minus;
};

class MinDelta : public Functional{
public:
        MinDelta( int cnum_nodes, const double cnodes[] );
        ~MinDelta();

        double getValue( double x ) const;
        bool hasDerivative() const;

        int getNumIntervals() const;
        double* getIntervals() const;

private:
        int num_nodes;
        double *nodes;
};

struct OptimizerResult{
        double xmax, fmax;
};

class Optimizer{
public:
        Optimizer();
        ~Optimizer();

        static OptimizerResult argMaxGlobal( const Functional *F ); // assumes range is [-1,1] and both points are included
        static OptimizerResult argMaxLocalPattern( const Functional *F, double left, double right );
        static double argMaxLocalSecant( const Functional *F, double left, double right );
};

}

#endif
