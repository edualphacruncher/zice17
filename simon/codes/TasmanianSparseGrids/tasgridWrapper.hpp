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

#ifndef __TASGRID_WRAPPER_HPP
#define __TASGRID_WRAPPER_HPP

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <string.h>
#include <math.h>

#include "TasmanianSparseGrid.hpp"

using std::cout;
using std::endl;
using std::setw;

using namespace TasGrid;

enum TypeCommand{
        command_none,
        command_makeglobal,
        command_makesequence,
        command_makelocalp,
        command_makewavelet,
        command_makequadrature,
        command_update,

        command_getquadrature,
        command_getinterweights,
        command_getpoints,
        command_getneeded,
        command_getrefcoeff,

        command_loadvalues,

        command_evaluate,
        command_integrate,

        command_getanisocoeff,
        command_refine_surp,
        command_refine_aniso,
        command_refine,
        command_refine_clear,

        command_getpoly,

        command_summary
};

class TasgridWrapper{
public:
        TasgridWrapper();
        ~TasgridWrapper();

        void setCommand( TypeCommand com );
        TypeCommand getCommand() const;

        void setNumDimensions( int dim );
        void setNumOutputs( int out );
        void setNumDepth( int d );
        void setOrder( int o );
        void setDepthType( TypeDepth dt );
        void setRule( TypeOneDRule r );
        TypeOneDRule getRule() const;
        void setAlpha( double a );
        void setBeta( double b );

        void setTolerance( double tol );
        void setRefOutput( int out );
        void setMinGrowth( int mg );
        void setTypeRefinement( TypeRefinement rt );

        void setGridFilename( const char *filename );
        void setOutFilename( const char *filename );
        void setValsFilename( const char *filename );
        void setXFilename( const char *filename );
        void setAnisoFilename( const char *filename );
        void setTransformFilename( const char *filename );
        void setCustomFilename( const char *filename );

        void setPrintPoints( bool pp );

        bool executeCommand();

        static bool isCreateCommand( TypeCommand com );

protected:
        bool checkSane() const;

        void createGlobalGird();
        void createSequenceGird();
        void createLocalPolynomialGird();
        void createWaveletGird();
        void createQuadrature();
        bool updateGrid();
        void writeGrid() const;
        bool readGrid();

        void outputPoints( bool useNeeded ) const;
        void outputQuadrature() const;

        bool loadValues();

        bool getInterWeights();
        bool getEvaluate();
        bool getIntegrate();
        bool getAnisoCoeff();

        bool refineGrid();
        bool cancelRefine();

        bool getPoly();

        bool getSummary();

        int* readAnisotropicFile( int num_weights ) const;
        double* readTransform() const;

        static void readMatrix( const char *filename, int &rows, int &cols, double* &mat );
        static void writeMatrix( const char *filename, int rows, int cols, const double mat[] );
        static void printMatrix( int rows, int cols, const double mat[] );

private:
        TasmanianSparseGrid *grid;

        TypeCommand command;

        int num_dimensions, num_outputs, depth, order;
        TypeDepth depth_type;
        TypeOneDRule rule;

        double alpha, beta;
        bool set_alpha, set_beta;

        double tolerance;
        bool set_tolerance;

        int ref_output, min_growth;
        TypeRefinement tref;
        bool set_tref;

        const char *gridfilename;
        const char *outfilename;
        const char *valsfilename;
        const char *xfilename;
        const char *anisofilename;
        const char *transformfilename;
        const char *customfilename;

        bool printCout;

};

#endif
