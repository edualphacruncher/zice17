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

#ifndef __TASGRID_TESTER_HPP
#define __TASGRID_TESTER_HPP

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <string.h>
#include <math.h>

#include "TasmanianSparseGrid.hpp"

#include "tasgridTestFunctions.hpp"

using std::cout;
using std::endl;
using std::setw;

using namespace TasGrid;

struct TestResults{
        double error;
        int num_points;
};

enum TestType{
        type_integration, type_nodal_interpolation, type_internal_interpolation
};

class ExternalTester{
public:
        ExternalTester(int in_num_mc = 1);
        ~ExternalTester();

        void Test() const;

        bool testGlobalRule( const BaseFunction *f, TasGrid::TypeOneDRule rule, const int *anisotropic, double alpha, double beta, bool interpolation, const int depths[], const double tols[], const char *custom_rule_filename = 0 ) const;
        bool performGLobalTest( const TasGrid::TypeOneDRule rule ) const;

        bool testLocalPolynomialRule( const BaseFunction *f, TasGrid::TypeOneDRule rule, const int depths[], const double tols[] ) const;
        bool testLocalWaveletRule( const BaseFunction *f, const int depths[], const double tols[] ) const;
        bool testSurplusRefinement( const BaseFunction *f, TasmanianSparseGrid *grid, double tol, TypeRefinement rtype, const int np[], const double errs[], int max_iter  ) const;
        bool testAnisotropicRefinement( const BaseFunction *f, TasmanianSparseGrid *grid, TypeDepth type, int min_growth, const int np[], const double errs[], int max_iter  ) const;

        bool testSurplusRefinementMake( const BaseFunction *f, TasmanianSparseGrid *grid, double tol, TypeRefinement rtype, int max_iter  ) const;
        bool testAnisotropicRefinementMake( const BaseFunction *f, TasmanianSparseGrid *grid, TypeDepth type, int min_growth, int max_iter  ) const;

        TestResults getError( const BaseFunction *f, TasGrid::TasmanianSparseGrid *grid, TestType type, const double *x = 0 ) const;

        bool testAllGlobal() const;
        bool testAllPWLocal() const;
        bool testAllWavelet() const;
        bool testAllRefinement() const;

        void debugTest(); // call this with -test debug
        void debugTestII(); // call this with -test debug

protected:

        void setRandomX( int n, double x[] ) const;

private:
        int num_mc;

        TwoOneExpNX2 f21nx2;
        TwoOneCos f21cos;
        TwoOneSinSin f21sinsin;
        TwoOneCosCos f21coscos;
        TwoOneDivisionAnisotropic f21aniso;
        TwoOne1DCurved f21curved;

        TwoOneConstGC1 f21constGC1;
        TwoOneConstGC2 f21constGC2;
        TwoOneConstGG f21constGG;
        TwoOneConstGJ f21constGJ;
        TwoOneConstGGL f21constGGL;
        TwoOneConstGH f21constGH;

        TwoOneENX2aniso f21nx2aniso;
        SixteenOneActive3 f16active3;
};

#endif
