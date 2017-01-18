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

#ifndef __TSG_CORE_ONE_DIMENSIONAL_HPP
#define __TSG_CORE_ONE_DIMENSIONAL_HPP

#include "tsgEnumerates.hpp"
#include "tsgLinearSolvers.hpp"

namespace TasGrid{

class CustomTabulated{
public:
        CustomTabulated();
        CustomTabulated( const char* filename );
        ~CustomTabulated();

        // I/O subroutines
        void write( std::ofstream &ofs ) const;
        bool read( std::ifstream &ifs );

        int getNumLevels() const;
        int getNumPoints( int level ) const;
        int getIExact( int level ) const;
        int getQExact( int level ) const;

        void getWeightsNodes( int level, double* &w, double* &x ) const;
        const char* getDescription() const;

protected:
        void reset();

private:
        int num_levels;
        int *num_nodes;
        int *precision;
        int *offsets;
        double *nodes;
        double *weights;
        std::string *description;
};

class OneDimensionalMeta{
public:
        OneDimensionalMeta();
        OneDimensionalMeta( const CustomTabulated *ccustom );
        ~OneDimensionalMeta();

        int getNumPoints( int level, TypeOneDRule rule ) const;
        int getIExact( int level, TypeOneDRule rule ) const;
        int getQExact( int level, TypeOneDRule rule ) const;

        const CustomTabulated *getCustom() const;

        static bool isNonNested( TypeOneDRule rule );
        static bool isSequence( TypeOneDRule rule );
        static bool isGlobal( TypeOneDRule rule );
        static bool isSingleNodeGrowth( TypeOneDRule rule );
        static bool isLocalPolynomial( TypeOneDRule rule );
        static bool isWavelet( TypeOneDRule rule );
        static TypeOneDRule getIORuleString( const char *name );
        static const char* getIORuleString( TypeOneDRule rule );
        static const char* getHumanString( TypeOneDRule rule );

private:
        // add the custom class here, but alias to something created by GlobalGrid
        const CustomTabulated *custom;
};

class OneDimensionalNodes{
public:
        OneDimensionalNodes();
        ~OneDimensionalNodes();

        // non-nested rules
        void getGaussLegendre( int m, double* &w, double* &x ) const;
        void getChebyshev( int m, double* &w, double* &x ) const;
        void getGaussChebyshev1( int m, double* &w, double* &x ) const;
        void getGaussChebyshev2( int m, double* &w, double* &x ) const;
        void getGaussJacobi( int m, double* &w, double* &x, double alpha, double beta ) const;
        void getGaussHermite( int m, double* &w, double* &x, double alpha ) const;
        void getGaussLaguerre( int m, double* &w, double* &x, double alpha ) const;

        // nested rules
        double* getClenshawCurtisNodes( int level ) const;
        double getClenshawCurtisWeight( int level, int point ) const;

        double* getClenshawCurtisNodesZero( int level ) const; // assuming zero boundary
        double getClenshawCurtisWeightZero( int level, int point ) const; // assuming zero boundary

        double* getFejer2Nodes( int level ) const;
        double getFejer2Weight( int level, int point ) const;

        double* getRLeja( int n ) const;
        double* getRLejaCentered( int n ) const;
        double* getRLejaShifted( int n ) const;
};

}

#endif
