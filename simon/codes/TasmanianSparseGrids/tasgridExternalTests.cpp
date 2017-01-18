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

#ifndef __TASGRID_TESTER_CPP
#define __TASGRID_TESTER_CPP

#include "tasgridExternalTests.hpp"

ExternalTester::ExternalTester(int in_num_mc ) : num_mc(in_num_mc){ /*srand( time(0) );*/ srand(10); }
ExternalTester::~ExternalTester(){}

void ExternalTester::setRandomX( int n, double x[] ) const{
        for( int i=0; i<n; i++ ){
                x[i] = 2.0 * ((double) rand()) / ( (double) RAND_MAX ) -1.0;
        };
}

void ExternalTester::Test() const{
        bool pass = true;

        cout << endl << endl;
        cout << "---------------------------------------------------------------------" << endl;
        cout << "       Tasmanian Sparse Grids Module: Functionality Test" << endl;
        cout << "---------------------------------------------------------------------" << endl << endl;


        pass = pass && testAllGlobal();
        pass = pass && testAllPWLocal();
        pass = pass && testAllWavelet();
        pass = pass && testAllRefinement();

        cout << endl;
        if ( pass ){
                cout << "---------------------------------------------------------------------" << endl;
                cout << "               All Tests Completed Successfully" << endl;
                cout << "---------------------------------------------------------------------" << endl << endl;
        }else{
                cout << "FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL" << endl;
                cout << "                 Some Tests Have Failed" << endl;
                cout << "FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL FAIL" << endl << endl;
        }
}

TestResults ExternalTester::getError( const BaseFunction *f, TasGrid::TasmanianSparseGrid *grid, TestType type, const double *x ) const{
        TestResults R;
        int num_dimensions = f->getNumInputs();
        int num_outputs = f->getNumOutputs();
        int num_points = grid->getNumPoints();
        if ( (type == type_integration) || (type == type_nodal_interpolation) ){
                double *points = grid->getPoints(), *weights = 0;
                if ( type == type_integration ){
                        weights = grid->getQuadratureWeights();
                }else{
                        weights = grid->getInterpolationWeights( x );
                }

                double *y = new double[num_outputs];
                double *r = new double[num_outputs];  std::fill( r, r + num_outputs, 0.0 );

                #pragma omp parallel
                {
                        double *y_local = new double[num_outputs];
                        double *r_local = new double[num_outputs];  std::fill( r_local, r_local + num_outputs, 0.0 );
                        #pragma omp for schedule(dynamic)
                        for( int i=0; i<num_points; i++ ){
                                f->eval( &(points[i*num_dimensions]), y_local );
                                for( int j=0; j<num_outputs; j++ ){
                                        r_local[j] += weights[i] * y_local[j];
                                };
                        };

                        #pragma omp critical (sum)
                        {
                                for( int j=0; j<num_outputs; j++ ){
                                        r[j] += r_local[j];
                                }
                        }
                        delete[] y_local;
                        delete[] r_local;
                }

                double err = 0.0;
                if ( type == type_integration ){
                        f->getIntegral( y );
                }else{
                        f->eval( x, y );
                }
                for( int j=0; j<num_outputs; j++ ){
                        err += fabs( y[j] - r[j] );
                };
                R.error = err;

                delete[] r;
                delete[] y;
                delete[] points;
                delete[] weights;
        }else if ( type == type_internal_interpolation ){
                // load needed points
                int num_needed_points = grid->getNumNeeded();
                if ( num_needed_points > 0 ){
                        double *needed_points = grid->getNeededPoints();
                        double *values = new double[num_outputs * num_needed_points];

                        for( int i=0; i<num_needed_points; i++ ){
                                f->eval( &(needed_points[i*num_dimensions]), &(values[i*num_outputs]) );
                        }

                        grid->loadNeededPoints( values );

                        delete[] values;
                        delete[] needed_points;
                }

                double *r = new double[num_outputs];  std::fill( r, r + num_outputs, 0.0 );
                double *n = new double[num_outputs];  std::fill( n, n + num_outputs, 0.0 );

                #pragma omp parallel default(none) shared(f,r,n,num_outputs,num_points,num_dimensions,grid)
                {
                	double *x_local = new double[num_dimensions];
                	double *n_local = new double[num_dimensions];
                	double *r_local = new double[num_dimensions];
                	double *y       = new double[num_dimensions];
                	double *s       = new double[num_dimensions];

                	for( int i=0; i<num_outputs; i++ ){
                		y[i] = s[i] = r_local[i] = n_local[i] = 0.0;
                	}
                        #pragma omp for schedule(dynamic)
                	for( int k=0; k<num_mc; k++ ){
                		setRandomX( num_dimensions, x_local );
                		f->eval( x_local, y );
                		grid->evaluate( x_local, s );
                		for( int j=0; j<num_outputs; j++ ){
                                        double e = fabs( y[j] - s[j] );
                                        if ( r_local[j] < e ) r_local[j] = e;
                                        e = fabs( y[j] );
                                        if ( n_local[j] < e ) n_local[j] = e;
                		}
                	}

                        #pragma omp critical (assignment)
                	{
                		for(int i = 0; i < num_outputs; i++){
                			if ( r[i] < r_local[i] ) r[i] = r_local[i];
                			if ( n[i] < n_local[i] ) n[i] = n_local[i];
                		}
                	}

                	delete[] x_local;
                	delete[] n_local;
                	delete[] r_local;
                	delete[] y;
                	delete[] s;
                }

                double err = r[0] / n[0];
                for( int j=1; j<num_outputs; j++ ){
                        err = ( err > r[j] / n[j] ) ? err : r[j] / n[j];
                }

                delete[] r;
                delete[] n;
                R.error = err;
        }
        R.num_points = grid->getNumPoints();
        return R;
}

bool ExternalTester::testGlobalRule( const BaseFunction *f, TasGrid::TypeOneDRule rule, const int *anisotropic, double alpha, double beta, const bool interpolation, const int depths[], const double tols[], const char *custom_rule_filename ) const{
        TasGrid::TasmanianSparseGrid grid;
        TestResults R;
        int num_global_tests = ( interpolation ) ? 3 : 1;
        TestType tests[3] = { type_integration, type_nodal_interpolation, type_internal_interpolation };
        TasGrid::TypeDepth type = TasGrid::type_iptotal;
        double *x = new double[f->getNumInputs()]; setRandomX(f->getNumInputs(),x);
        bool bPass = true;
        const char *custom_filename = ( rule == rule_customtabulated ) ? "GaussPattersonRule.table" : 0;
        for( int i=0; i<num_global_tests; i++ ){
                if ( interpolation ){
                        grid.makeGlobalGrid( f->getNumInputs(), f->getNumOutputs(), depths[i], type, rule, anisotropic, alpha, beta, custom_filename );
                        R = getError( f, &grid, tests[i], x );
                }else{
                        grid.makeGlobalGrid( f->getNumInputs(), 0, depths[i], type, rule, anisotropic, alpha, beta, custom_filename );
                        R = getError( f, &grid, type_integration );
                }
                if ( R.error > tols[i] ){
                        bPass = false;
                        cout << setw(18) << "ERROR: FAILED" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( rule );
                        if ( interpolation ){
                                if ( tests[i%3] == type_integration ){
                                        cout << setw(25) << "integration test";
                                }else if ( tests[i%3] == type_nodal_interpolation ){
                                        cout << setw(25) << "w-interpolation";
                                }else{
                                        cout << setw(25) << "interpolation";
                                }
                        }else{
                                cout << setw(25) << "integration test";
                        }
                        cout << "   failed function: " << f->getDescription();
                        cout << setw(10) << "observed: " << R.error << "  expected: " << tols[i] << endl;
                }
        }
        delete[] x;
        return bPass;
}

bool ExternalTester::performGLobalTest( TasGrid::TypeOneDRule rule ) const{
        double alpha = 0.3, beta = 0.7;
        bool pass = true;
        if ( rule == TasGrid::rule_clenshawcurtis ){
                { TasGrid::TypeOneDRule oned = TasGrid::rule_clenshawcurtis;
                const int depths1[3] = { 25, 25, 25 };
                const double tols1[3] = { 1.E-12, 1.E-12, 1.E-11 };
                const int depths2[3] = { 25, 27, 27 };
                const double tols2[3] = { 1.E-12, 1.E-10, 1.E-11 };
                if ( testGlobalRule( &f21nx2, oned, 0, alpha, beta, true, depths1, tols1 ) && testGlobalRule( &f21cos, oned, 0, alpha, beta, true, depths2, tols2 ) ){
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "Pass" << endl;
                }else{
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "FAIL" << endl;  pass = false;
                }}
        }else if ( rule == TasGrid::rule_clenshawcurtis0 ){
                { TasGrid::TypeOneDRule oned = TasGrid::rule_clenshawcurtis0;
                const int depths1[3] = { 25, 25, 25 };
                const double tols1[3] = { 1.E-12, 1.E-12, 1.E-11 };
                if ( testGlobalRule( &f21sinsin, oned, 0, alpha, beta, true, depths1, tols1 ) ){
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "Pass" << endl;
                }else{
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "FAIL" << endl;  pass = false;
                }}
        }else if ( rule == TasGrid::rule_chebyshev ){
                { TasGrid::TypeOneDRule oned = TasGrid::rule_chebyshev;
                const int depths1[3] = { 22, 22, 22 };
                const double tols1[3] = { 1.E-12, 1.E-10, 1.E-10 };
                const int depths2[3] = { 22, 22, 22 };
                const double tols2[3] = { 1.E-12, 1.E-09, 1.E-09 };
                if ( testGlobalRule( &f21nx2, oned, 0, alpha, beta, true, depths1, tols1 ) && testGlobalRule( &f21cos, oned, 0, alpha, beta, true, depths2, tols2 ) ){
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "Pass" << endl;
                }else{
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "FAIL" << endl;  pass = false;
                }}
        }else if ( rule == TasGrid::rule_chebyshevodd ){
                { TasGrid::TypeOneDRule oned = TasGrid::rule_chebyshevodd;
                const int depths1[3] = { 22, 22, 22 };
                const double tols1[3] = { 1.E-12, 1.E-10, 1.E-10 };
                const int depths2[3] = { 22, 22, 22 };
                const double tols2[3] = { 1.E-12, 1.E-09, 1.E-09 };
                if ( testGlobalRule( &f21nx2, oned, 0, alpha, beta, true, depths1, tols1 ) && testGlobalRule( &f21cos, oned, 0, alpha, beta, true, depths2, tols2 ) ){
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "Pass" << endl;
                }else{
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "FAIL" << endl;  pass = false;
                }}
        }else if ( rule == TasGrid::rule_leja ){
                { TasGrid::TypeOneDRule oned = TasGrid::rule_leja;
                const int depths1[3] = { 20, 20, 20 };
                const double tols1[3] = { 3.E-10, 5.E-09, 5.E-09 };
                const int depths2[3] = { 20, 20, 20 };
                const double tols2[3] = { 3.E-09, 5.E-08, 5.E-08 };
                if ( testGlobalRule( &f21nx2, oned, 0, alpha, beta, true, depths1, tols1 ) && testGlobalRule( &f21cos, oned, 0, alpha, beta, true, depths2, tols2 ) ){
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "Pass" << endl;
                }else{
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "FAIL" << endl;  pass = false;
                }}
        }else if ( rule == TasGrid::rule_lejaodd ){
                { TasGrid::TypeOneDRule oned = TasGrid::rule_lejaodd;
                const int depths1[3] = { 20, 20, 20 };
                const double tols1[3] = { 3.E-10, 5.E-09, 5.E-09 };
                const int depths2[3] = { 20, 20, 20 };
                const double tols2[3] = { 3.E-09, 5.E-08, 5.E-08 };
                if ( testGlobalRule( &f21nx2, oned, 0, alpha, beta, true, depths1, tols1 ) && testGlobalRule( &f21cos, oned, 0, alpha, beta, true, depths2, tols2 ) ){
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "Pass" << endl;
                }else{
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "FAIL" << endl;  pass = false;
                }}
        }else if ( rule == TasGrid::rule_rleja ){
                { TasGrid::TypeOneDRule oned = TasGrid::rule_rleja;
                const int depths1[3] = { 20, 20, 20 };
                const double tols1[3] = { 3.E-10, 1.E-08, 1.E-08 };
                const int depths2[3] = { 20, 20, 20 };
                const double tols2[3] = { 3.E-09, 5.E-08, 5.E-08 };
                if ( testGlobalRule( &f21nx2, oned, 0, alpha, beta, true, depths1, tols1 ) && testGlobalRule( &f21cos, oned, 0, alpha, beta, true, depths2, tols2 ) ){
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "Pass" << endl;
                }else{
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "FAIL" << endl;  pass = false;
                }}
        }else if ( rule == TasGrid::rule_rlejadouble2 ){
                { TasGrid::TypeOneDRule oned = TasGrid::rule_rlejadouble2;
                const int depths1[3] = { 25, 25, 25 };
                const double tols1[3] = { 1.E-12, 1.E-11, 1.E-11 };
                const int depths2[3] = { 25, 27, 27 };
                const double tols2[3] = { 1.E-12, 1.E-10, 1.E-10 };
                if ( testGlobalRule( &f21nx2, oned, 0, alpha, beta, true, depths1, tols1 ) && testGlobalRule( &f21cos, oned, 0, alpha, beta, true, depths2, tols2 ) ){
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "Pass" << endl;
                }else{
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "FAIL" << endl;  pass = false;
                }}
        }else if ( rule == TasGrid::rule_rlejadouble4 ){
                { TasGrid::TypeOneDRule oned = TasGrid::rule_rlejadouble4;
                const int depths1[3] = { 25, 25, 25 };
                const double tols1[3] = { 1.E-12, 1.E-11, 1.E-11 };
                const int depths2[3] = { 25, 27, 27 };
                const double tols2[3] = { 1.E-12, 1.E-10, 1.E-10 };
                if ( testGlobalRule( &f21nx2, oned, 0, alpha, beta, true, depths1, tols1 ) && testGlobalRule( &f21cos, oned, 0, alpha, beta, true, depths2, tols2 ) ){
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "Pass" << endl;
                }else{
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "FAIL" << endl;  pass = false;
                }}
        }else if ( rule == TasGrid::rule_rlejaodd ){
                { TasGrid::TypeOneDRule oned = TasGrid::rule_rlejaodd;
                const int depths1[3] = { 20, 20, 20 };
                const double tols1[3] = { 3.E-10, 5.E-09, 5.E-09 };
                const int depths2[3] = { 20, 20, 20 };
                const double tols2[3] = { 3.E-09, 5.E-08, 5.E-08 };
                if ( testGlobalRule( &f21nx2, oned, 0, alpha, beta, true, depths1, tols1 ) && testGlobalRule( &f21cos, oned, 0, alpha, beta, true, depths2, tols2 ) ){
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "Pass" << endl;
                }else{
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "FAIL" << endl;  pass = false;
                }}
        }else if ( rule == TasGrid::rule_rlejashifted ){
                { TasGrid::TypeOneDRule oned = TasGrid::rule_rlejashifted;
                const int depths1[3] = { 20, 20, 20 };
                const double tols1[3] = { 3.E-10, 1.E-08, 1.E-08 };
                const int depths2[3] = { 20, 20, 20 };
                const double tols2[3] = { 3.E-09, 5.E-08, 5.E-08 };
                if ( testGlobalRule( &f21nx2, oned, 0, alpha, beta, true, depths1, tols1 ) && testGlobalRule( &f21cos, oned, 0, alpha, beta, true, depths2, tols2 ) ){
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "Pass" << endl;
                }else{
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "FAIL" << endl;  pass = false;
                }}
        }else if ( rule == TasGrid::rule_rlejashiftedeven ){
                { TasGrid::TypeOneDRule oned = TasGrid::rule_rlejashiftedeven;
                const int depths1[3] = { 20, 20, 20 };
                const double tols1[3] = { 3.E-10, 5.E-09, 5.E-09 };
                const int depths2[3] = { 20, 20, 20 };
                const double tols2[3] = { 6.E-09, 5.E-08, 5.E-08 };
                if ( testGlobalRule( &f21nx2, oned, 0, alpha, beta, true, depths1, tols1 ) && testGlobalRule( &f21cos, oned, 0, alpha, beta, true, depths2, tols2 ) ){
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "Pass" << endl;
                }else{
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "FAIL" << endl;  pass = false;
                }}
        }else if ( rule == TasGrid::rule_maxlebesgue ){
                { TasGrid::TypeOneDRule oned = TasGrid::rule_maxlebesgue;
                const int depths1[3] = { 20, 20, 20 };
                const double tols1[3] = { 3.E-10, 5.E-09, 5.E-09 };
                const int depths2[3] = { 20, 20, 20 };
                const double tols2[3] = { 3.E-09, 5.E-08, 5.E-08 };
                if ( testGlobalRule( &f21nx2, oned, 0, alpha, beta, true, depths1, tols1 ) && testGlobalRule( &f21cos, oned, 0, alpha, beta, true, depths2, tols2 ) ){
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "Pass" << endl;
                }else{
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "FAIL" << endl;  pass = false;
                }}
        }else if ( rule == TasGrid::rule_maxlebesgueodd ){
                { TasGrid::TypeOneDRule oned = TasGrid::rule_maxlebesgueodd;
                const int depths1[3] = { 20, 20, 20 };
                const double tols1[3] = { 3.E-10, 5.E-09, 5.E-09 };
                const int depths2[3] = { 20, 20, 20 };
                const double tols2[3] = { 3.E-09, 5.E-08, 5.E-08 };
                if ( testGlobalRule( &f21nx2, oned, 0, alpha, beta, true, depths1, tols1 ) && testGlobalRule( &f21cos, oned, 0, alpha, beta, true, depths2, tols2 ) ){
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "Pass" << endl;
                }else{
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "FAIL" << endl;  pass = false;
                }}
        }else if ( rule == TasGrid::rule_minlebesgue ){
                { TasGrid::TypeOneDRule oned = TasGrid::rule_minlebesgue;
                const int depths1[3] = { 20, 20, 20 };
                const double tols1[3] = { 3.E-10, 5.E-09, 5.E-09 };
                const int depths2[3] = { 20, 20, 20 };
                const double tols2[3] = { 3.E-09, 5.E-08, 5.E-08 };
                if ( testGlobalRule( &f21nx2, oned, 0, alpha, beta, true, depths1, tols1 ) && testGlobalRule( &f21cos, oned, 0, alpha, beta, true, depths2, tols2 ) ){
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "Pass" << endl;
                }else{
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "FAIL" << endl;  pass = false;
                }}
        }else if ( rule == TasGrid::rule_minlebesgueodd ){
                { TasGrid::TypeOneDRule oned = TasGrid::rule_minlebesgueodd;
                const int depths1[3] = { 20, 20, 20 };
                const double tols1[3] = { 3.E-10, 5.E-09, 5.E-09 };
                const int depths2[3] = { 20, 20, 20 };
                const double tols2[3] = { 3.E-09, 5.E-08, 5.E-08 };
                if ( testGlobalRule( &f21nx2, oned, 0, alpha, beta, true, depths1, tols1 ) && testGlobalRule( &f21cos, oned, 0, alpha, beta, true, depths2, tols2 ) ){
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "Pass" << endl;
                }else{
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "FAIL" << endl;  pass = false;
                }}
        }else if ( rule == TasGrid::rule_mindelta ){
                { TasGrid::TypeOneDRule oned = TasGrid::rule_mindelta;
                const int depths1[3] = { 20, 20, 20 };
                const double tols1[3] = { 3.E-10, 5.E-09, 5.E-09 };
                const int depths2[3] = { 20, 20, 20 };
                const double tols2[3] = { 3.E-09, 5.E-08, 5.E-08 };
                if ( testGlobalRule( &f21nx2, oned, 0, alpha, beta, true, depths1, tols1 ) && testGlobalRule( &f21cos, oned, 0, alpha, beta, true, depths2, tols2 ) ){
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "Pass" << endl;
                }else{
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "FAIL" << endl;  pass = false;
                }}
        }else if ( rule == TasGrid::rule_mindeltaodd ){
                { TasGrid::TypeOneDRule oned = TasGrid::rule_mindeltaodd;
                const int depths1[3] = { 20, 20, 20 };
                const double tols1[3] = { 3.E-10, 5.E-09, 5.E-09 };
                const int depths2[3] = { 20, 20, 20 };
                const double tols2[3] = { 3.E-09, 5.E-08, 5.E-08 };
                if ( testGlobalRule( &f21nx2, oned, 0, alpha, beta, true, depths1, tols1 ) && testGlobalRule( &f21cos, oned, 0, alpha, beta, true, depths2, tols2 ) ){
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "Pass" << endl;
                }else{
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "FAIL" << endl;  pass = false;
                }}
        }else if ( rule == TasGrid::rule_gausslegendre ){
                { TasGrid::TypeOneDRule oned = TasGrid::rule_gausslegendre;
                const int depths1[3] = { 20, 36, 38 };
                const double tols1[3] = { 1.E-10, 1.E-07, 1.E-07 };
                const int depths2[3] = { 24, 36, 36 };
                const double tols2[3] = { 1.E-10, 1.E-07, 1.E-07 };
                if ( testGlobalRule( &f21nx2, oned, 0, alpha, beta, true, depths1, tols1 ) && testGlobalRule( &f21cos, oned, 0, alpha, beta, true, depths2, tols2 ) ){
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "Pass" << endl;
                }else{
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "FAIL" << endl;  pass = false;
                }}
        }else if ( rule == TasGrid::rule_gausslegendreodd ){
                { TasGrid::TypeOneDRule oned = TasGrid::rule_gausslegendreodd;
                const int depths1[3] = { 20, 36, 38 };
                const double tols1[3] = { 1.E-10, 1.E-07, 1.E-07 };
                const int depths2[3] = { 24, 36, 36 };
                const double tols2[3] = { 1.E-10, 1.E-07, 1.E-07 };
                if ( testGlobalRule( &f21nx2, oned, 0, alpha, beta, true, depths1, tols1 ) && testGlobalRule( &f21cos, oned, 0, alpha, beta, true, depths2, tols2 ) ){
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "Pass" << endl;
                }else{
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "FAIL" << endl;  pass = false;
                }}
        }else if ( rule == TasGrid::rule_gausspatterson ){
                { TasGrid::TypeOneDRule oned = TasGrid::rule_gausspatterson;
                const int depths1[3] = { 20, 36, 38 };
                const double tols1[3] = { 1.E-10, 1.E-07, 1.E-07 };
                const int depths2[3] = { 24, 36, 36 };
                const double tols2[3] = { 1.E-10, 1.E-07, 1.E-07 };
                if ( testGlobalRule( &f21nx2, oned, 0, alpha, beta, true, depths1, tols1 ) && testGlobalRule( &f21cos, oned, 0, alpha, beta, true, depths2, tols2 ) ){
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "Pass" << endl;
                }else{
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "FAIL" << endl;  pass = false;
                }}
        }else if ( rule == TasGrid::rule_customtabulated ){
                {
                std::ifstream ftest( "GaussPattersonRule.table" );
                if ( !ftest.good() ){
                        ftest.close();
                        cout << "WARNING: cannot find GaussPattersonRule.table file and cannot test the custom rule!" << endl;
                }else{
                        ftest.close();
                        TasGrid::TypeOneDRule oned = TasGrid::rule_customtabulated;
                        const int depths1[3] = { 20, 36, 38 };
                        const double tols1[3] = { 1.E-10, 1.E-07, 1.E-07 };
                        const int depths2[3] = { 24, 36, 36 };
                        const double tols2[3] = { 1.E-10, 1.E-07, 1.E-07 };
                        if ( testGlobalRule( &f21nx2, oned, 0, alpha, beta, true, depths1, tols1 ) && testGlobalRule( &f21cos, oned, 0, alpha, beta, true, depths2, tols2 ) ){
                                cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "Pass" << endl;
                        }else{
                                cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "FAIL" << endl;  pass = false;
                        }
                }}
        }else if ( rule == TasGrid::rule_fejer2 ){
                { TasGrid::TypeOneDRule oned = TasGrid::rule_fejer2;
                const int depths1[3] = { 20, 40, 40 };
                const double tols1[3] = { 1.E-14, 1.E-12, 1.E-12 };
                if ( testGlobalRule( &f21coscos, oned, 0, alpha, beta, true, depths1, tols1 ) ){
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "Pass" << endl;
                }else{
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "FAIL" << endl;  pass = false;
                }}
        }else if ( rule == TasGrid::rule_gausschebyshev1 ){
                { TasGrid::TypeOneDRule oned = TasGrid::rule_gausschebyshev1;
                const int depths1[3] = { 20, 20, 20 };
                const double tols1[3] = { 5.E-14, 1.E-05, 1.E-05 };
                if ( testGlobalRule( &f21constGC1, oned, 0, alpha, beta, true, depths1, tols1 ) ){
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "Pass" << endl;
                }else{
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "FAIL" << endl;  pass = false;
                }}
        }else if ( rule == TasGrid::rule_gausschebyshev1odd ){
                { TasGrid::TypeOneDRule oned = TasGrid::rule_gausschebyshev1odd;
                const int depths1[3] = { 20, 20, 20 };
                const double tols1[3] = { 1.E-14, 1.E-05, 1.E-05 };
                if ( testGlobalRule( &f21constGC1, oned, 0, alpha, beta, true, depths1, tols1 ) ){
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "Pass" << endl;
                }else{
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "FAIL" << endl;  pass = false;
                }}
        }else if ( rule == TasGrid::rule_gausschebyshev2 ){
                { TasGrid::TypeOneDRule oned = TasGrid::rule_gausschebyshev2;
                const int depths1[3] = { 20, 20, 20 };
                const double tols1[3] = { 1.E-14, 1.E-05, 1.E-05 };
                if ( testGlobalRule( &f21constGC2, oned, 0, alpha, beta, true, depths1, tols1 ) ){
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "Pass" << endl;
                }else{
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "FAIL" << endl;  pass = false;
                }}
        }else if ( rule == TasGrid::rule_gausschebyshev2odd ){
                { TasGrid::TypeOneDRule oned = TasGrid::rule_gausschebyshev2odd;
                const int depths1[3] = { 20, 20, 20 };
                const double tols1[3] = { 1.E-14, 1.E-05, 1.E-05 };
                if ( testGlobalRule( &f21constGC2, oned, 0, alpha, beta, true, depths1, tols1 ) ){
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "Pass" << endl;
                }else{
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "FAIL" << endl;  pass = false;
                }}
        }else if ( rule == TasGrid::rule_gaussgegenbauer ){
                { TasGrid::TypeOneDRule oned = TasGrid::rule_gaussgegenbauer;
                const int depths1[3] = { 20, 20, 20 };
                const double tols1[3] = { 1.E-11, 1.E-05, 1.E-05 };
                if ( testGlobalRule( &f21constGG, oned, 0, alpha, beta, true, depths1, tols1 ) ){
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "Pass" << endl;
                }else{
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "FAIL" << endl;  pass = false;
                }}
        }else if ( rule == TasGrid::rule_gaussgegenbauerodd ){
                { TasGrid::TypeOneDRule oned = TasGrid::rule_gaussgegenbauerodd;
                const int depths1[3] = { 20, 20, 20 };
                const double tols1[3] = { 1.E-11, 1.E-05, 1.E-05 };
                if ( testGlobalRule( &f21constGG, oned, 0, alpha, beta, true, depths1, tols1 ) ){
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "Pass" << endl;
                }else{
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "FAIL" << endl;  pass = false;
                }}
        }else if ( rule == TasGrid::rule_gaussjacobi ){
                { TasGrid::TypeOneDRule oned = TasGrid::rule_gaussjacobi;
                const int depths1[3] = { 20, 20, 20 };
                const double tols1[3] = { 1.E-08, 1.E-05, 1.E-05 };
                if ( testGlobalRule( &f21constGJ, oned, 0, alpha, beta, true, depths1, tols1 ) ){
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "Pass" << endl;
                }else{
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "FAIL" << endl;  pass = false;
                }}
        }else if ( rule == TasGrid::rule_gaussjacobiodd ){
                { TasGrid::TypeOneDRule oned = TasGrid::rule_gaussjacobiodd;
                const int depths1[3] = { 20, 20, 20 };
                const double tols1[3] = { 1.E-08, 1.E-05, 1.E-05 };
                if ( testGlobalRule( &f21constGJ, oned, 0, alpha, beta, true, depths1, tols1 ) ){
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "Pass" << endl;
                }else{
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "FAIL" << endl;  pass = false;
                }}
        }else if ( rule == TasGrid::rule_gausslaguerre ){
                { TasGrid::TypeOneDRule oned = TasGrid::rule_gausslaguerre;
                const int depths1[1] = { 20 };
                const double tols1[1] = { 1.E-08 };
                if ( testGlobalRule( &f21constGGL, oned, 0, alpha, beta, false, depths1, tols1 ) ){
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "Pass" << endl;
                }else{
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "FAIL" << endl;  pass = false;
                }}
        }else if ( rule == TasGrid::rule_gausslaguerreodd ){
                { TasGrid::TypeOneDRule oned = TasGrid::rule_gausslaguerreodd;
                const int depths1[1] = { 20 };
                const double tols1[1] = { 1.E-08 };
                if ( testGlobalRule( &f21constGGL, oned, 0, alpha, beta, false, depths1, tols1 ) ){
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "Pass" << endl;
                }else{
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "FAIL" << endl;  pass = false;
                }}
        }else if ( rule == TasGrid::rule_gausshermite ){
                { TasGrid::TypeOneDRule oned = TasGrid::rule_gausshermite;
                const int depths1[1] = { 20 };
                const double tols1[1] = { 1.E-09 };
                if ( testGlobalRule( &f21constGH, oned, 0, alpha, beta, false, depths1, tols1 ) ){
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "Pass" << endl;
                }else{
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "FAIL" << endl;  pass = false;
                }}
        }else if ( rule == TasGrid::rule_gausshermiteodd ){
                { TasGrid::TypeOneDRule oned = TasGrid::rule_gausshermiteodd;
                const int depths1[1] = { 20 };
                const double tols1[1] = { 1.E-09 };
                if ( testGlobalRule( &f21constGH, oned, 0, alpha, beta, false, depths1, tols1 ) ){
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "Pass" << endl;
                }else{
                        cout << setw(10) << "Rule" << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned ) << setw(25) << "FAIL" << endl;  pass = false;
                }}
        }
        return pass;
}

bool ExternalTester::testAllGlobal() const{
        bool pass = true;
        const int nrules = 35; // sync with below
        TasGrid::TypeOneDRule rules[nrules] = {
                                        TasGrid::rule_chebyshev,
                                        TasGrid::rule_chebyshevodd,
                                        TasGrid::rule_clenshawcurtis,
                                        TasGrid::rule_clenshawcurtis0,
                                        TasGrid::rule_fejer2,
                                        TasGrid::rule_leja,
                                        TasGrid::rule_lejaodd,
                                        TasGrid::rule_rleja,
                                        TasGrid::rule_rlejadouble2,
                                        TasGrid::rule_rlejadouble4,
                                        TasGrid::rule_rlejaodd,
                                        TasGrid::rule_rlejashifted,
                                        TasGrid::rule_rlejashiftedeven,
                                        TasGrid::rule_maxlebesgue,
                                        TasGrid::rule_maxlebesgueodd,
                                        TasGrid::rule_minlebesgue,
                                        TasGrid::rule_minlebesgueodd,
                                        TasGrid::rule_mindelta,
                                        TasGrid::rule_mindeltaodd,
                                        TasGrid::rule_gausslegendre,
                                        TasGrid::rule_gausslegendreodd,
                                        TasGrid::rule_gausspatterson,
                                        TasGrid::rule_gausschebyshev1,
                                        TasGrid::rule_gausschebyshev1odd,
                                        TasGrid::rule_gausschebyshev2,
                                        TasGrid::rule_gausschebyshev2odd,
                                        TasGrid::rule_gaussgegenbauer,
                                        TasGrid::rule_gaussgegenbauerodd,
                                        TasGrid::rule_gaussjacobi,
                                        TasGrid::rule_gaussjacobiodd,
                                        TasGrid::rule_gausslaguerre,
                                        TasGrid::rule_gausslaguerreodd,
                                        TasGrid::rule_gausshermite,
                                        TasGrid::rule_gausshermiteodd,
                                        TasGrid::rule_customtabulated   };

        for( int i=0; i<nrules; i++ ){
                if ( !performGLobalTest( rules[i] ) ){
                        pass = false;
                }
        }
        return pass;
}

bool ExternalTester::testLocalPolynomialRule( const BaseFunction *f, TasGrid::TypeOneDRule rule, const int depths[], const double tols[] ) const{
        TasGrid::TasmanianSparseGrid grid;
        TestResults R;
        TestType tests[3] = { type_integration, type_nodal_interpolation, type_internal_interpolation };
        int orders[5] = { 1, 2, 3, 4, -1 };
        double *x = new double[f->getNumInputs()]; setRandomX(f->getNumInputs(),x);
        bool bPass = true;
        for( int i=0; i<15; i++ ){
                grid.makeLocalPolynomialGrid( f->getNumInputs(), f->getNumOutputs(), depths[i], orders[i/3], rule );
                R = getError( f, &grid, tests[i%3], x );
                if ( R.error > tols[i] ){
                        bPass = false;
                        cout << setw(18) << "ERROR: FAILED";
                        cout << setw(6) << TasGrid::OneDimensionalMeta::getIORuleString( rule );
                        cout << " order: " << orders[i%3];

                        if ( tests[i%3] == type_integration ){
                                cout << setw(25) << "integration test";
                        }else if ( tests[i%3] == type_nodal_interpolation ){
                                cout << setw(25) << "w-interpolation";
                        }else{
                                cout << setw(25) << "interpolation";
                        }

                        cout << "   failed function: " << f->getDescription();
                        cout << setw(10) << "observed: " << R.error << "  expected: " << tols[i] << endl;
                }
        }
        delete[] x;
        return bPass;
}

bool ExternalTester::testSurplusRefinement( const BaseFunction *f, TasmanianSparseGrid *grid, double tol, TypeRefinement rtype, const int np[], const double errs[], int max_iter  ) const{
        for( int itr=0; itr<max_iter; itr++ ){
                TestResults R = getError( f, grid, type_internal_interpolation );
                if (  (R.num_points != np[itr]) || (R.error > errs[itr])  ){
                        cout << setw(18) << "ERROR: FAILED refinement test at iteration: " << itr << endl;
                        cout << " expected: " << np[itr] << "  " << errs[itr] << "   computed: " << R.num_points << "  " << R.error << endl;
                        return false;
                }
                if ( (grid->isGlobal()) || (grid->isSequence())  ){
                        grid->setSurplusRefinement( tol );
                }else{
                        grid->setSurplusRefinement( tol, rtype );
                }
        }
        return true;
}
bool ExternalTester::testSurplusRefinementMake( const BaseFunction *f, TasmanianSparseGrid *grid, double tol, TypeRefinement rtype, int max_iter  ) const{
        for( int itr=0; itr<max_iter; itr++ ){
                TestResults R = getError( f, grid, type_internal_interpolation );
                cout << std::scientific;  cout.precision(16);
                cout << setw(10) << R.num_points << setw(30) << R.error << endl;
                if ( (grid->isGlobal()) || (grid->isSequence())  ){
                        grid->setSurplusRefinement( tol );
                }else{
                        grid->setSurplusRefinement( tol, rtype );
                }
        }
        return true;
}
bool ExternalTester::testAnisotropicRefinement( const BaseFunction *f, TasmanianSparseGrid *grid, TypeDepth type, int min_growth, const int np[], const double errs[], int max_iter  ) const{
        for( int itr=0; itr<max_iter; itr++ ){
                TestResults R = getError( f, grid, type_internal_interpolation );
                if (  (R.num_points != np[itr]) || (R.error > errs[itr])  ){
                        cout << setw(18) << "ERROR: FAILED refinement test at iteration: " << itr << endl;
                        cout << " expected: " << np[itr] << "  " << errs[itr] << "   computed: " << R.num_points << "  " << R.error << endl;
                        return false;
                }
                grid->setAnisotropicRefinement( type, min_growth );
        }
        return true;
}
bool ExternalTester::testAnisotropicRefinementMake( const BaseFunction *f, TasmanianSparseGrid *grid, TypeDepth type, int min_growth, int max_iter  ) const{
        TestResults *R = new TestResults[max_iter];
        for( int itr=0; itr<max_iter; itr++ ){
                R[itr] = getError( f, grid, type_internal_interpolation );
                grid->setAnisotropicRefinement( type, min_growth );
        }
        cout << std::scientific;  cout.precision(2);
        cout << "int np[" << max_iter << "] = { " << R[0].num_points;
        for( int i=1; i<max_iter; i++ ){
                cout << ", " << R[i].num_points;
        }
        cout << " };" << endl;
        cout << "int err[" << max_iter << "] = { " << R[0].error;
        for( int i=1; i<max_iter; i++ ){
                cout << ", " << R[i].error;
        }
        cout << " };" << endl;

        delete[] R;
        return true;
}

bool ExternalTester::testAllPWLocal() const{
        bool pass = true;
        { TasGrid::TypeOneDRule oned = TasGrid::rule_semilocalp;
        const int depths1[15] = { 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8 };
        const double tols1[15] = { 1.E-03, 1.E-03, 1.E-03, 1.E-07, 1.E-04, 1.E-04, 1.E-07, 1.E-05, 1.E-05, 1.E-07, 4.E-06, 4.E-06, 1.E-07, 4.E-06, 4.E-06 };
        if ( testLocalPolynomialRule( &f21nx2, oned, depths1, tols1 ) ){
                cout << setw(10) << "Rule";
                cout << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned );
                cout << setw(25) << "Pass" << endl;
        }else{
                cout << setw(10) << "Rule";
                cout << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned );
                cout << setw(25) << "FAIL" << endl; pass = false;
        }}
        { TasGrid::TypeOneDRule oned = TasGrid::rule_localp;
        const int depths1[15] = { 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8 };
        const double tols1[15] = { 1.E-03, 1.E-03, 1.E-03, 1.E-07, 1.E-04, 1.E-04, 1.E-07, 1.E-05, 1.E-05, 1.E-07, 4.E-06, 4.E-06, 1.E-07, 4.E-06, 4.E-06 };
        if ( testLocalPolynomialRule( &f21nx2, oned, depths1, tols1 ) ){
                cout << setw(10) << "Rule";
                cout << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned );
                cout << setw(25) << "Pass" << endl;
        }else{
                cout << setw(10) << "Rule";
                cout << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned );
                cout << setw(25) << "FAIL" << endl; pass = false;
        }}

        { TasGrid::TypeOneDRule oned = TasGrid::rule_localp0;
        const int depths1[15] = { 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8 };
        const double tols1[15] = { 1.E-03, 2.E-04, 2.E-04, 1.E-09, 1.E-06, 1.E-06, 1.E-09, 3.E-08, 3.E-08, 1.E-09, 4.E-09, 4.E-09, 1.E-09, 4.E-09, 4.E-09 };
        if ( testLocalPolynomialRule( &f21coscos, oned, depths1, tols1 ) ){
                cout << setw(10) << "Rule";
                cout << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned );
                cout << setw(25) << "Pass" << endl;
        }else{
                cout << setw(10) << "Rule";
                cout << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( oned );
                cout << setw(25) << "FAIL" << endl; pass = false;
        }}
        return pass;
}

bool ExternalTester::testLocalWaveletRule( const BaseFunction *f, const int depths[], const double tols[] ) const{
        TasGrid::TasmanianSparseGrid grid;
        TestResults R;
        TestType tests[3] = { type_integration, type_nodal_interpolation, type_internal_interpolation };
        int orders[2] = { 1, 3 };
        double *x = new double[f->getNumInputs()]; setRandomX(f->getNumInputs(),x);
        bool bPass = true;
        for( int i=0; i<6; i++ ){
                grid.makeWaveletGrid( f->getNumInputs(), f->getNumOutputs(), depths[i], orders[i/3] );
                R = getError( f, &grid, tests[i%3], x );
                if ( R.error > tols[i] ){
                        bPass = false;
                        cout << setw(18) << "ERROR: FAILED";
                        cout << setw(6) << TasGrid::OneDimensionalMeta::getIORuleString( rule_wavelet );
                        cout << " order: " << orders[i%3];

                        if ( tests[i%3] == type_integration ){
                                cout << setw(25) << "integration test";
                        }else if ( tests[i%3] == type_nodal_interpolation ){
                                cout << setw(25) << "w-interpolation";
                        }else{
                                cout << setw(25) << "interpolation";
                        }

                        cout << "   failed function: " << f->getDescription();
                        cout << setw(10) << "observed: " << R.error << "  expected: " << tols[i] << endl;
                }
        }
        delete[] x;
        return bPass;
}
bool ExternalTester::testAllWavelet() const{
        bool pass = true;
        const int depths1[6] = { 7, 7, 7, 5, 5, 5 };
        const double tols1[6] = { 5.E-05, 1.E-04, 1.E-04, 1.E-08, 1.E-07, 1.E-07 };
        if ( testLocalWaveletRule( &f21nx2, depths1, tols1 ) ){
                cout << setw(10) << "Rule";
                cout << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( rule_wavelet );
                cout << setw(25) << "Pass" << endl;
        }else{
                cout << setw(10) << "Rule";
                cout << setw(25) << TasGrid::OneDimensionalMeta::getIORuleString( rule_wavelet );
                cout << setw(25) << "FAIL" << endl; pass = false;
        }
        return pass;
}

bool ExternalTester::testAllRefinement() const{
        TasmanianSparseGrid grid;
        bool pass = true;
        {
                const BaseFunction *f = &f21nx2;
                grid.makeGlobalGrid( f->getNumInputs(), f->getNumOutputs(), 3, type_iptotal, rule_leja );
                int np[13] = { 10, 15, 21, 28, 36, 45, 55, 66, 78, 91, 105, 118, 130 };
                double err[13] = { 2, 2.E-1, 5.E-1, 2.E-2, 4.E-2, 2.E-3, 4.E-3, 2.E-4, 2.E-4, 2.E-5, 2.E-5, 8.E-7, 8.E-7 };
                if ( !testSurplusRefinement( f, &grid, 1.E-6, refine_classic, np, err, 13 ) ){
                        cout << "ERROR: failed leja surplus refinement for " << f->getDescription() << endl;  pass = false;
                }
        }{
                const BaseFunction *f = &f21coscos;
                int np[9] = {        21,    24,    30,    39,    49,    60,    72,    79,    85 };
                double err[9] = { 2.E-1, 7.E-3, 2.E-2, 3.E-4, 6.E-4, 4.E-6, 9.E-6, 5.E-7, 5.E-7 };
                grid.makeGlobalGrid( f->getNumInputs(), f->getNumOutputs(), 5, type_iptotal, rule_rleja );
                if ( !testSurplusRefinement( f, &grid, 1.E-6, refine_classic, np, err, 9 ) ){
                        cout << "ERROR: failed rleja global surplus refinement for " << f->getDescription() << endl;  pass = false;
                }
        }{
                const BaseFunction *f = &f21coscos;
                int np[9] = {        21,    24,    30,    39,    49,    60,    72,    79,    85 };
                double err[9] = { 2.E-1, 7.E-3, 2.E-2, 3.E-4, 6.E-4, 4.E-6, 9.E-6, 5.E-7, 5.E-7 };
                grid.makeSequenceGrid( f->getNumInputs(), f->getNumOutputs(), 5, type_iptotal, rule_rleja );
                if ( !testSurplusRefinement( f, &grid, 1.E-6, refine_classic, np, err, 9 ) ){
                        cout << "ERROR: failed rleja sequence surplus refinement for " << f->getDescription() << endl;  pass = false;
                }
        }{
                const BaseFunction *f = &f21nx2;
                grid.makeLocalPolynomialGrid( f->getNumInputs(), f->getNumOutputs(), 3, 2, rule_semilocalp );
                int np[8] = { 29, 65, 145, 321, 705, 1521, 2753, 3569 };
                double err[8] = { 4.E-2, 1.E-2, 1.E-3, 2.E-4, 4.E-5, 5.E-6, 1.E-6, 5.E-7 };
                if ( !testSurplusRefinement( f, &grid, 1.E-6, refine_classic, np, err, 8 ) ){
                        cout << "ERROR: failed semi-local classic refinement for " << f->getDescription() << endl;  pass = false;
                }
                grid.makeLocalPolynomialGrid( f->getNumInputs(), f->getNumOutputs(), 3, 2, rule_semilocalp );
        }{
                const BaseFunction *f = &f21nx2;
                grid.makeLocalPolynomialGrid( f->getNumInputs(), f->getNumOutputs(), 3, 2, rule_semilocalp );
                int np[8] = { 29, 65, 145, 321, 705, 1521, 2753, 3569 };
                double err[8] = { 4.E-2, 1.E-2, 1.E-3, 2.E-4, 4.E-5, 5.E-6, 1.E-6, 5.E-7 };
                if ( !testSurplusRefinement( f, &grid, 1.E-6, refine_parents_first, np, err, 8 ) ){
                        cout << "ERROR: failed semi-local parents refinement for " << f->getDescription() << endl;  pass = false;
                }
        }{
                const BaseFunction *f = &f21nx2;
                grid.makeLocalPolynomialGrid( f->getNumInputs(), f->getNumOutputs(), 2, 2, rule_semilocalp );
                int np[6] = { 13, 29, 65, 145, 321, 545 };
                double err[6] = { 8.E-02, 5.E-02, 7.E-03, 2.E-03, 3.E-04, 6.E-05 };
                if ( !testSurplusRefinement( f, &grid, 1.E-4, refine_direction_selective, np, err, 6 ) ){
                        cout << "ERROR: failed semi-local direction refinement for " << f->getDescription() << endl;  pass = false;
                }
        }{
                const BaseFunction *f = &f21nx2;
                grid.makeLocalPolynomialGrid( f->getNumInputs(), f->getNumOutputs(), 2, 2, rule_semilocalp );
                int np[6] = { 13, 29, 65, 145, 321, 545 };
                double err[6] = { 8.E-02, 5.E-02, 7.E-03, 2.E-03, 3.E-04, 6.E-05 };
                if ( !testSurplusRefinement( f, &grid, 1.E-4, refine_fds, np, err, 6 ) ){
                        cout << "ERROR: failed semi-local fds refinement for " << f->getDescription() << endl;  pass = false;
                }
        }{
                const BaseFunction *f = &f21nx2;
                grid.makeLocalPolynomialGrid( f->getNumInputs(), f->getNumOutputs(), 3, 1, rule_localp );
                int np[10] = { 29, 65, 145, 321, 705, 1537, 3321, 6981, 13517, 19113 };
                double err[10] = { 4.E-2, 2.E-2, 6.E-3, 2.E-3, 6.E-4, 2.E-4, 6.E-5, 2.E-5, 6.E-6, 2.E-6 };
                if ( !testSurplusRefinement( f, &grid, 1.E-6, refine_classic, np, err, 10 ) ){
                        cout << "ERROR: failed localp classic refinement for " << f->getDescription() << endl;  pass = false;
                }
        }{
                const BaseFunction *f = &f21nx2;
                grid.makeLocalPolynomialGrid( f->getNumInputs(), f->getNumOutputs(), 3, 1, rule_localp );
                int np[10] = { 29, 65, 145, 321, 705, 1537, 3321, 6981, 13517, 19113 };
                double err[10] = { 4.E-2, 2.E-2, 6.E-3, 2.E-3, 6.E-4, 2.E-4, 6.E-5, 2.E-5, 6.E-6, 2.E-6 };
                if ( !testSurplusRefinement( f, &grid, 1.E-6, refine_parents_first, np, err, 10 ) ){
                        cout << "ERROR: failed localp parents refinement for " << f->getDescription() << endl;  pass = false;
                }
        }{
                const BaseFunction *f = &f21nx2;
                grid.makeLocalPolynomialGrid( f->getNumInputs(), f->getNumOutputs(), 2, 1, rule_localp );
                int np[8] = { 13, 29, 65, 145, 321, 673, 1233, 1433 };
                double err[8] = { 1.E-01, 5.E-02, 2.E-02, 5.E-03, 2.E-03, 6.E-04, 2.E-04, 1.E-04 };
                if ( !testSurplusRefinement( f, &grid, 1.E-4, refine_direction_selective, np, err, 8 ) ){
                        cout << "ERROR: failed localp direction refinement for " << f->getDescription() << endl;  pass = false;
                }
        }{
                const BaseFunction *f = &f21nx2;
                grid.makeLocalPolynomialGrid( f->getNumInputs(), f->getNumOutputs(), 2, 1, rule_localp );
                int np[8] = { 13, 29, 65, 145, 321, 673, 1233, 1433 };
                double err[8] = { 1.E-01, 5.E-02, 2.E-02, 5.E-03, 2.E-03, 6.E-04, 2.E-04, 1.E-04 };
                if ( !testSurplusRefinement( f, &grid, 1.E-4, refine_fds, np, err, 8 ) ){
                        cout << "ERROR: failed localp fds refinement for " << f->getDescription() << endl;  pass = false;
                }
        }{
                const BaseFunction *f = &f21coscos;
                grid.makeLocalPolynomialGrid( f->getNumInputs(), f->getNumOutputs(), 3, 2, rule_localp0 );
                int np[6] = { 49, 129, 321, 769, 1761, 2209 };
                double err[6] = { 2.E-3, 3.E-4, 5.E-5, 7.E-6, 8.E-7, 5.E-7 };
                if ( !testSurplusRefinement( f, &grid, 1.E-6, refine_classic, np, err, 6 ) ){
                        cout << "ERROR: failed localp-zero classic refinement for " << f->getDescription() << endl;  pass = false;
                }
        }{
                const BaseFunction *f = &f21coscos;
                grid.makeLocalPolynomialGrid( f->getNumInputs(), f->getNumOutputs(), 3, 2, rule_localp0 );
                int np[6] = { 49, 129, 321, 769, 1761, 2209 };
                double err[6] = { 2.E-3, 3.E-4, 5.E-5, 7.E-6, 8.E-7, 5.E-7 };
                if ( !testSurplusRefinement( f, &grid, 1.E-6, refine_parents_first, np, err, 6 ) ){
                        cout << "ERROR: failed localp-zero parents refinement for " << f->getDescription() << endl;  pass = false;
                }
        }{
                const BaseFunction *f = &f21coscos;
                grid.makeLocalPolynomialGrid( f->getNumInputs(), f->getNumOutputs(), 2, 2, rule_localp0 );
                int np[4] = { 17, 49, 129, 305 };
                double err[4] = { 7.E-03, 2.E-03, 4.E-04, 4.E-05 };
                if ( !testSurplusRefinement( f, &grid, 1.E-4, refine_direction_selective, np, err, 4 ) ){
                        cout << "ERROR: failed localp-zero direction refinement for " << f->getDescription() << endl;  pass = false;
                }
        }{
                const BaseFunction *f = &f21coscos;
                grid.makeLocalPolynomialGrid( f->getNumInputs(), f->getNumOutputs(), 2, 2, rule_localp0 );
                int np[4] = { 17, 49, 129, 305 };
                double err[4] = { 7.E-03, 2.E-03, 4.E-04, 4.E-05 };
                if ( !testSurplusRefinement( f, &grid, 1.E-4, refine_fds, np, err, 4 ) ){
                        cout << "ERROR: failed localp-zero fds refinement for " << f->getDescription() << endl;  pass = false;
                }
        }{
                const BaseFunction *f = &f21coscos;
                grid.makeWaveletGrid( f->getNumInputs(), f->getNumOutputs(), 2, 1 );
                int np[7] = { 49, 81, 193, 449, 993, 1921, 1937 };
                double err[7] = { 6.E-02, 3.E-02, 6.E-03, 3.E-03, 6.E-04, 3.E-04, 2.E-04 };
                if ( !testSurplusRefinement( f, &grid, 1.E-4, refine_classic, np, err, 7 ) ){
                        cout << "ERROR: failed wavelet classic refinement for " << f->getDescription() << endl;  pass = false;
                }
        }{
                const BaseFunction *f = &f21coscos;
                grid.makeWaveletGrid( f->getNumInputs(), f->getNumOutputs(), 2, 1 );
                int np[7] = { 49, 81, 193, 449, 993, 1921, 1937 };
                double err[7] = { 6.E-02, 3.E-02, 6.E-03, 3.E-03, 6.E-04, 3.E-04, 2.E-04 };
                if ( !testSurplusRefinement( f, &grid, 1.E-4, refine_parents_first, np, err, 7 ) ){
                        cout << "ERROR: failed wavelet parents refinement for " << f->getDescription() << endl;  pass = false;
                }
        }{
                const BaseFunction *f = &f21nx2;
                grid.makeWaveletGrid( f->getNumInputs(), f->getNumOutputs(), 2, 1 );
                int np[6] = { 49, 113, 257, 561, 1113, 1481 };
                double err[6] = { 6.E-02, 1.E-02, 5.E-03, 1.E-03, 5.E-04, 1.E-04 };
                if ( !testSurplusRefinement( f, &grid, 1.E-4, refine_direction_selective, np, err, 6 ) ){
                        cout << "ERROR: failed wavelet direction refinement for " << f->getDescription() << endl;  pass = false;
                }
        }{
                const BaseFunction *f = &f21coscos;
                grid.makeWaveletGrid( f->getNumInputs(), f->getNumOutputs(), 2, 1 );
                int np[7] = { 49, 81, 161, 385, 889, 1737, 1769 };
                double err[7] = { 6.E-02, 3.E-02, 6.E-03, 3.E-03, 6.E-04, 3.E-04, 2.E-04 };
                if ( !testSurplusRefinement( f, &grid, 1.E-4, refine_fds, np, err, 7 ) ){
                        cout << "ERROR: failed wavelet fds refinement for " << f->getDescription() << endl;  pass = false;
                }
        }

        if ( pass ){
                cout << "      Refinement            surplus" << setw(25) << "Pass" << endl;
        }else{
                cout << "      Refinement            surplus" << setw(25) << "FAIL" << endl;
        }

        pass = true;
        {
                const BaseFunction *f = &f21aniso;
                grid.makeGlobalGrid( f->getNumInputs(), f->getNumOutputs(), 3, type_iptotal, rule_leja );
                int np[40] = { 10, 15, 21, 28, 29, 30, 31, 32, 34, 35, 37, 40, 41, 45, 49, 54, 59, 64, 70, 77, 84, 92, 100, 108, 117, 126, 135, 145, 155, 165, 176, 187, 198, 210, 212, 224, 237, 248, 251, 263 };
                double errs[40] = { 9.04e-01, 4.24e-01, 5.73e-01, 2.78e-01, 3.15e-01, 2.49e-01, 3.00e-01, 8.85e-02, 9.30e-02, 9.67e-02, 2.06e-01, 3.03e-01, 5.24e-02, 4.63e-02, 5.85e-02, 5.11e-02, 9.80e-03, 2.71e-02, 5.42e-03, 7.85e-03, 6.21e-03, 5.41e-03, 2.56e-03, 3.32e-03, 5.18e-04, 6.14e-04, 3.66e-04, 4.87e-04, 8.19e-05, 2.58e-04, 5.76e-05, 5.54e-05, 5.22e-05, 4.89e-05, 4.68e-05, 8.92e-06, 2.20e-05, 5.56e-06, 5.14e-06, 5.79e-06 };
                if ( !testAnisotropicRefinement( f, &grid, type_iptotal, 1, np, errs, 40 ) ){
                        cout << "ERROR: failed anisotropic refinement using leja iptotal nodes for " << f->getDescription() << endl;  pass = false;
                }
        }{
                const BaseFunction *f = &f21curved;
                grid.makeGlobalGrid( f->getNumInputs(), f->getNumOutputs(), 3, type_iptotal, rule_leja );
                int np[10] = { 10, 12, 17, 24, 32, 34, 41, 42, 57, 59 };
                double errs[10] = { 9.48e-03, 9.50e-03, 6.85e-03, 5.11e-04, 6.26e-05, 7.11e-06, 5.07e-06, 5.19e-06, 1.17e-08, 1.86e-08 };
                if ( !testAnisotropicRefinement( f, &grid, type_ipcurved, 1, np, errs, 10 ) ){
                        cout << "ERROR: failed anisotropic refinement using leja iptotal nodes for " << f->getDescription() << endl;  pass = false;
                }
        }{
                const BaseFunction *f = &f21curved;
                grid.makeGlobalGrid( f->getNumInputs(), f->getNumOutputs(), 3, type_iptotal, rule_clenshawcurtis );
                int np[3] = { 13, 21, 33 };
                double errs[3] = { 6.12e-04, 6.05e-04, 1.33e-08 };
                if ( !testAnisotropicRefinement( f, &grid, type_ipcurved, 1, np, errs, 3 ) ){
                        cout << "ERROR: failed anisotropic refinement using leja iptotal nodes for " << f->getDescription() << endl;  pass = false;
                }
        }{
                const BaseFunction *f = &f21curved;
                grid.makeSequenceGrid( f->getNumInputs(), f->getNumOutputs(), 3, type_iptotal, rule_leja );
                int np[10] = { 10, 12, 17, 24, 32, 34, 41, 42, 57, 59 };
                double errs[10] = { 9.48e-03, 9.50e-03, 6.85e-03, 5.11e-04, 6.26e-05, 7.11e-06, 5.07e-06, 5.19e-06, 1.17e-08, 1.86e-08 };
                if ( !testAnisotropicRefinement( f, &grid, type_ipcurved, 1, np, errs, 10 ) ){
                        cout << "ERROR: failed anisotropic refinement using leja iptotal nodes for " << f->getDescription() << endl;  pass = false;
                }
        }

        if ( pass ){
                cout << "      Refinement        anisotropic" << setw(25) << "Pass" << endl;
        }else{
                cout << "      Refinement        anisotropic" << setw(25) << "FAIL" << endl;
        }

        return true;
}

void ExternalTester::debugTest(){
        cout << "Debug Test" << endl;
        //cout << std::scientific;  cout.precision(17);

        cout << "Put here testing code and call this with ./tasgrid -test debug" << endl;
}

void ExternalTester::debugTestII(){
        cout << "Debug Test II" << endl;
        //cout << std::scientific;  cout.precision(17);

        cout << "Put here testing code and call this with ./tasgrid -test db" << endl;
}

#endif
