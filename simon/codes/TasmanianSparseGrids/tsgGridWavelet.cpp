/*
 * Code Author: Drayton Munster, July 2013
 * Edited: Miroslav Stoyanov, Sep 2015
 *
 * Copyright (C) 2013  Drayton Munster
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

#ifndef __TASMANIAN_SPARSE_GRID_WAVELET_CPP
#define __TASMANIAN_SPARSE_GRID_WAVELET_CPP

#include "tsgGridWavelet.hpp"

namespace TasGrid{

GridWavelet::GridWavelet() : num_dimensions(0), num_outputs(0), order(1), coefficients(0), points(0), needed(0), values(0), inter_matrix(0)
{}
GridWavelet::~GridWavelet(){ reset(); }

void GridWavelet::reset(){
        if ( coefficients != 0 ){  delete[] coefficients; coefficients = 0;  }
        if ( points != 0 ){  delete points; points = 0;  }
        if ( needed != 0 ){  delete needed; needed = 0;  }
        if ( values != 0 ){  delete values; values = 0;  }
        if ( inter_matrix != 0 ){  delete inter_matrix; inter_matrix = 0;  }
}

void GridWavelet::write( std::ofstream &ofs ) const{
        ofs << std::scientific; ofs.precision(17);
        ofs << num_dimensions << " " << num_outputs << " " << order << endl;
        if ( num_dimensions > 0 ){
                if ( points == 0 ){
                        ofs << "0" << endl;
                }else{
                        ofs << "1 ";
                        points->write( ofs );
                }
                if ( coefficients == 0 ){
                        ofs << "0" << endl;
                }else{
                        ofs << "1 ";
                        for( int i=0; i<points->getNumIndexes() * num_outputs; i++ ){  ofs << " " << coefficients[i];  } ofs << endl;
                }
                if ( needed == 0 ){
                        ofs << "0" << endl;
                }else{
                        ofs << "1 ";
                        needed->write( ofs );
                }
                if ( num_outputs > 0 ) values->write( ofs );
        }
}
void GridWavelet::read( std::ifstream &ifs ){
        reset();
        int flag;
        ifs >> num_dimensions >> num_outputs >> order;
        if ( num_dimensions > 0 ){
                rule1D.updateOrder( order );

                ifs >> flag;  if ( flag == 1 ){  points    = new IndexSet( num_dimensions );   points->read( ifs );  }
                ifs >> flag;  if ( flag == 1 ){  coefficients = new double[ points->getNumIndexes() * num_outputs ]; for( int i=0; i<points->getNumIndexes() * num_outputs; i++ ){  ifs >> coefficients[i];  }  }
                ifs >> flag;  if ( flag == 1 ){  needed    = new IndexSet( num_dimensions );   needed->read( ifs );  }

                if ( num_outputs > 0 ){  values = new StorageSet( 0, 0 ); values->read( ifs );  }
        }
}

void GridWavelet::makeGrid( int cnum_dimensions, int cnum_outputs, int depth, int corder ){
        reset();
        num_dimensions = cnum_dimensions;
        num_outputs = cnum_outputs;
        order = corder;

        rule1D.updateOrder( order );

        IndexManipulator IM(num_dimensions);
        UnsortedIndexSet* deltas = IM.getToalDegreeDeltas( depth );

        WaveletLevels wavelevels( order ); // describes the wavelet hierarchy, namely points per delta

        needed = IM.generatePointsFromDeltas( deltas, &wavelevels );
        delete deltas;

        if ( num_outputs == 0 ){
                points = needed;
                needed = 0;
        }else{
                values = new StorageSet( num_outputs, needed->getNumIndexes() );
        }

        buildInterpolationMatrix();
}
void GridWavelet::setNodes( IndexSet* &nodes, int cnum_outputs, int corder ){
        reset();
        num_dimensions = nodes->getNumDimensions();
        num_outputs = cnum_outputs;
        order = corder;

        rule1D.updateOrder( order );

        needed = nodes; nodes = 0;

        if ( num_outputs == 0 ){
                points = needed;
                needed = 0;
        }else{
                values = new StorageSet( num_outputs, needed->getNumIndexes() );
        }

        buildInterpolationMatrix();
}

int GridWavelet::getNumDimensions() const{  return num_dimensions;  }
int GridWavelet::getNumOutputs() const{  return num_outputs;  }
TypeOneDRule GridWavelet::getRule() const{  return rule_wavelet;  }
int GridWavelet::getOrder() const{  return order;  }

int GridWavelet::getNumLoaded() const{  return ((points == 0) ? 0 : points->getNumIndexes());  }
int GridWavelet::getNumNeeded() const{  return ((needed == 0) ? 0 : needed->getNumIndexes());  }
int GridWavelet::getNumPoints() const{  return ((points == 0) ? getNumNeeded() : getNumLoaded());  }

double* GridWavelet::getLoadedPoints() const{
        if ( points == 0 ) return 0;
        int num_points = points->getNumIndexes();
        if ( num_points == 0 ) return 0;
        double *x = new double[ num_dimensions * num_points ];
        #pragma omp parallel for schedule(static)
        for( int i=0; i<num_points; i++ ){
                const int *p = points->getIndex( i );
                for( int j=0; j<num_dimensions; j++ ){
                        x[i*num_dimensions + j] = rule1D.getX( p[j] );
                }
        }
        return x;
}
double* GridWavelet::getNeededPoints() const{
        if ( needed == 0 ) return 0;
        int num_points = needed->getNumIndexes();
        if ( num_points == 0 ) return 0;
        double *x = new double[ num_dimensions * num_points ];
        #pragma omp parallel for schedule(static)
        for( int i=0; i<num_points; i++ ){
                const int *p = needed->getIndex( i );
                for( int j=0; j<num_dimensions; j++ ){
                        x[i*num_dimensions + j] = rule1D.getX( p[j] );
                }
        }
        return x;
}
double* GridWavelet::getPoints() const{
        return (( points == 0 ) ? getNeededPoints() : getLoadedPoints());
}

double* GridWavelet::getQuadratureWeights() const{
        IndexSet *work = ( points == 0 ) ? needed : points;
        int num_points = work->getNumIndexes();
	double *weights = new double[num_points];
	#pragma omp parallel for
	for( int i=0; i<num_points; i++ ){
		weights[i] = evalIntegral( work->getIndex(i) );
	}

	solveTransposed(weights);
	return weights;
}
double* GridWavelet::getInterpolationWeights( const double x[] ) const{
        IndexSet *work = ( points == 0 ) ? needed : points;

        int num_points = work->getNumIndexes();
	double *weights = new double[num_points];
	#pragma omp parallel for
	for( int i=0; i<num_points; i++ ){
                weights[i] = evalBasis( work->getIndex(i), x );
	}

	solveTransposed(weights);
	return weights;
}
void GridWavelet::loadNeededPoints( const double *vals ){
        if ( points == 0 ){
                values->setValues( vals );
                points = needed;
                needed = 0;
        }else{
                values->addValues( points, needed, vals );
                points->addIndexSet( needed );
                delete needed; needed = 0;
                buildInterpolationMatrix();
        }
        recomputeCoefficients();
}
void GridWavelet::evaluate( const double x[], double y[] ) const{
        int num_points = points->getNumIndexes();
	double *basis_values = new double[num_points];
	#pragma omp parallel for
	for( int i=0; i<num_points; i++ ){
                basis_values[i] = evalBasis( points->getIndex(i), x );
	}
	for( int j=0; j<num_outputs; j++ ){
                double sum = 0.0;
                #pragma omp parallel for reduction( + : sum )
                for( int i=0; i<num_points; i++ ){
                                sum += basis_values[i] * coefficients[i*num_outputs + j];
                }
                y[j] = sum;
	}
	delete[] basis_values;
}
void GridWavelet::integrate( double q[] ) const{
        int num_points = points->getNumIndexes();
	double *basis_integrals = new double[num_points];
	#pragma omp parallel for
	for( int i=0; i<num_points; i++ ){
                basis_integrals[i] = evalIntegral( points->getIndex(i) );
	}
	for( int j=0; j<num_outputs; j++ ){
                double sum = 0.0;

                #pragma omp parallel for reduction( + : sum )
                for( int i=0; i<num_points; i++ ){
                        sum += basis_integrals[i] * coefficients[i*num_outputs + j];
                }
                q[j] = sum;
	}
	delete[] basis_integrals;
}

double GridWavelet::evalBasis( const int p[], const double x[] ) const{
	/*
	 * Evaluates the wavelet basis given at point p at the coordinates given by x.
	 */
	double v = 1.;
	for(int i = 0; i < num_dimensions; i++){
		v *= rule1D.eval(p[i], x[i]);
		if ( v == 0 ){ break; }; // MIRO: reduce the expensive wavelet evaluations
	}
	return v;
}
double GridWavelet::evalIntegral( const int p[] ) const{
	/*
	 * For a given node p, evaluate the integral of the associated wavelet.
	 */
	double v = 1.;
	for(int i = 0; i < num_dimensions; i++){
		v *= rule1D.getWeight(p[i]);
	}
	return v;
}

void GridWavelet::buildInterpolationMatrix(){
	/*
	 * Using the IndexSet points, constructs the interpolation matrix needed for methods
	 * such as recomputeCoefficients, solveTransposed, etc.
	 */
	IndexSet *work = ( points == 0 ) ? needed : points;
	if(inter_matrix != 0) { delete inter_matrix; }

	int num_points = work->getNumIndexes();

	TasSparse::TsgSparseCOO coo_mat(num_points, num_points);

#ifdef _OPENMP
	int num_threads = omp_get_max_threads();
	TasSparse::TsgSparseCOO *coos = new TasSparse::TsgSparseCOO[num_threads];
	#pragma omp parallel
#endif
	{
		double *xs = new double[num_dimensions];

#ifdef _OPENMP
		int thread_id = omp_get_thread_num();
		#pragma omp for
#endif
		for(int i = 0; i < num_points; i++){ /* Loop over points */
			const int *point = work->getIndex(i);

			for(int k = 0; k < num_dimensions; k++){ /* Fill in x values of point */
				xs[k] = rule1D.getX(point[k]);
			}

			for(int j = 0; j < num_points; j++){ /* Loop over wavelets */
				const int *wavelet = work->getIndex(j);
				double v = 1.;

				for(int k = 0; k < num_dimensions; k++){ /* Loop over dimensions */
					v *= rule1D.eval(wavelet[k], xs[k]);
					if ( v == 0.0 ){ break; }; // MIRO: evaluating the wavelets is expensive, stop if any one of them is zero
				} /* End for dimensions */

				if(v != 0.0){
				/*
				 * Testing for equal to zero is safe since v*0 is the only way
				 * for this to happen.
				 */
#ifdef _OPENMP
					coos[thread_id].addPoint(i, j, v);
#else
					coo_mat.addPoint(i,j,v);
#endif
				}


			} /* End for wavelets */

		} /* End for points */

		delete[] xs;

#ifdef _OPENMP
		#pragma omp master
		{ /* Master thread combines the sub-thread work */

			for(int i = 0; i < num_threads; i++){
				coo_mat.combine(coos[i]);
			}
			delete[] coos;
		} /* End master section */
#endif
	} /* End parallel section */
	inter_matrix = new TasSparse::SparseMatrix(coo_mat);
}

void GridWavelet::recomputeCoefficients(){
	/*
	 * Recalculates the coefficients to interpolate the values in points.
	 * Make sure buildInterpolationMatrix has been called since the list was updated.
	 */

	int num_points = points->getNumIndexes();
	if ( coefficients != 0 ){ delete[] coefficients; }
	coefficients = new double[num_points * num_outputs];

	if ( (inter_matrix == 0) || (inter_matrix->getNumRows() != num_points) ){
                buildInterpolationMatrix();
	}

	double *workspace = new double[2*num_points];
	double *b = workspace;
	double *x = workspace + num_points;

	for(int output = 0; output < num_outputs; output++){
		// Copy relevant portion
		std::fill( workspace, workspace + 2*num_points, 0.0 );

		// Populate RHS
		for(int i = 0; i < num_points; i++){
			b[i] = values->getValues( i )[output];
		}

		// Solve system
		inter_matrix->solve( b, x );

		// Populate surplus
		for(int i = 0; i < num_points; i++){
			coefficients[num_outputs * i + output] = x[i];
		}
	}

	delete[] workspace;
}

void GridWavelet::solveTransposed(double w[]) const{
	/*
	 * Solves the system A^T * w = y. Used to calculate interpolation and integration
	 * weights. RHS values should be passed in through w. At exit, w will contain the
	 * required weights.
	 */
	int num_points = inter_matrix->getNumRows();

	double *y = new double[num_points];

	std::copy( w, w + num_points, y );

	inter_matrix->solve( y, w, true );

	delete[] y;
}

double* GridWavelet::getNormalization() const{
        double* norm = new double[num_outputs];  std::fill( norm, norm + num_outputs, 0.0 );
        for( int i=0; i<points->getNumIndexes(); i++ ){
                const double *v = values->getValues( i );
                for( int j=0; j<num_outputs; j++ ){
                        if ( norm[j] < fabs( v[j] ) ) norm[j] = fabs( v[j] );
                }
        }
        return norm;
}
int* GridWavelet::buildUpdateMap( double tolerance, TypeRefinement criteria ) const{
        int num_points = points->getNumIndexes();
        int *map = new int[ num_points * num_dimensions ];  std::fill( map, map + num_points * num_dimensions, 0 );

        double *norm = getNormalization();

        if ( (criteria == refine_classic) || (criteria == refine_parents_first) ){
                // classic refinement
                #pragma omp parallel for
                for( int i=0; i<num_points; i++ ){
                        bool small = true;
                        for( int k=0; k<num_outputs; k++ ){
                                if ( small && ( (fabs( coefficients[i*num_outputs + k] ) / norm[k]) > tolerance ) ) small = false;
                        }
                        if ( !small ){
                                std::fill( &(map[i*num_dimensions]), &(map[i*num_dimensions]) + num_dimensions, 1 );
                        }
                }
        }else{
                SplitDirections split( points );

                for( int s=0; s<split.getNumJobs(); s++ ){
                        int d = split.getJobDirection( s );
                        int nump = split.getJobNumPoints( s );
                        const int *pnts = split.getJobPoints( s );

                        int *global_to_pnts = new int[num_points];

                        double *vals = new double[ nump * num_outputs ];
                        int *indexes = new int[ nump * num_dimensions ];

                        for( int i=0; i<nump; i++ ){
                                const double* v = values->getValues( pnts[i] );
                                std::copy( v, v + num_outputs, &( vals[i*num_outputs] ) );
                                const int *p = points->getIndex( pnts[i] );
                                std::copy( p, p + num_dimensions, &( indexes[i*num_dimensions] ) );
                                global_to_pnts[ pnts[i] ] = i;
                        }
                        IndexSet *pointset = new IndexSet( num_dimensions, nump, indexes );

                        GridWavelet direction_grid;
                        direction_grid.setNodes( pointset, num_outputs, order );
                        direction_grid.loadNeededPoints( vals );
                        const double *coeff = direction_grid.coefficients;

                        for( int i=0; i<nump; i++ ){
                                bool small = true;
                                for( int k=0; k<num_outputs; k++ ){
                                        if ( small && ( (fabs( coefficients[pnts[i]*num_outputs + k] ) / norm[k]) > tolerance ) && ( (fabs( coeff[i*num_outputs + k] ) / norm[k] ) > tolerance ) ) small = false;
                                }
                                map[ pnts[i]*num_dimensions + d ] = ( small ) ? 0 : 1;;
                        }

                        delete[] vals;
                        delete[] global_to_pnts;
                }
        }

        delete[] norm;

        return map;
}

bool GridWavelet::addParent( const int point[], int direction, GranulatedIndexSet *destination, IndexSet *exclude ) const{
        int *dad = new int[num_dimensions];  std::copy( point, point + num_dimensions, dad );
        bool added = false;
        dad[direction] = rule1D.getParent( point[direction] );
        if ( dad[direction] == -2 ){
                for( int c=0; c<rule1D.getNumPoints(0); c++ ){
                        dad[direction] = c;
                        if ( exclude->getSlot( dad ) == -1 ){
                                destination->addIndex( dad );
                                added = true;
                        }
                }
        }else if ( dad[direction] >= 0 ){
                if ( exclude->getSlot( dad ) == -1 ){
                        destination->addIndex( dad );
                        added = true;
                }
        }

        delete[] dad;
        return added;
}
void GridWavelet::addChild( const int point[], int direction, GranulatedIndexSet *destination, IndexSet *exclude )const{
        int *kid = new int[num_dimensions];  std::copy( point, point + num_dimensions, kid );
        int L, R; rule1D.getChildren( point[direction], L, R );
        kid[direction] = L;
        if ( (kid[direction] != -1) && (exclude->getSlot(kid) == -1) ){
                destination->addIndex( kid );
        }
        kid[direction] = R;
        if ( (kid[direction] != -1) && (exclude->getSlot(kid) == -1) ){
                destination->addIndex( kid );
        }
        delete[] kid;
}

void GridWavelet::clearRefinement(){
        if ( needed != 0 ){  delete needed;  needed = 0;  }
}
void GridWavelet::setSurplusRefinement( double tolerance, TypeRefinement criteria ){
        clearRefinement();

        int *map = buildUpdateMap( tolerance, criteria );

        bool useParents = (criteria == refine_fds) || (criteria == refine_parents_first);

        GranulatedIndexSet *refined = new GranulatedIndexSet( num_dimensions );

        int num_points = points->getNumIndexes();

        for( int i=0; i<num_points; i++ ){
                for( int j=0; j<num_dimensions; j++ ){
                        if ( map[i*num_dimensions+j] == 1 ){ // if this dimension needs to be refined
                                if ( !(useParents && addParent( points->getIndex(i), j, refined, points )) ){
                                        addChild( points->getIndex(i), j, refined, points );
                                }
                        }
                }
        }

        if ( refined->getNumIndexes() > 0 ){
                needed = new IndexSet( refined );
        }
        delete refined;

        delete[] map;
}

}

#endif
