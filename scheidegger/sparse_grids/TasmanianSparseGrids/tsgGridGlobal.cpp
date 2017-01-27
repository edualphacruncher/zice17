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

#ifndef __TASMANIAN_SPARSE_GRID_GLOBAL_CPP
#define __TASMANIAN_SPARSE_GRID_GLOBAL_CPP

#include "tsgGridGlobal.hpp"

namespace TasGrid{

GridGlobal::GridGlobal() : num_dimensions(0), num_outputs(0), alpha(0.0), beta(0.0), wrapper(0), tensors(0), active_tensors(0), active_w(0),
        points(0), needed(0), tensor_refs(0), max_levels(0), values(0),
        updated_tensors(0), updated_active_tensors(0), updated_active_w(0), custom(0)
{}
GridGlobal::~GridGlobal(){ reset( true ); }

void GridGlobal::write( std::ofstream &ofs ) const{
        ofs << std::scientific; ofs.precision(17);
        ofs << num_dimensions << " " << num_outputs << " " << alpha << " " << beta << endl;
        if ( num_dimensions > 0 ){
                ofs << OneDimensionalMeta::getIORuleString( rule ) << endl;
                if ( rule == rule_customtabulated ){
                        custom->write( ofs );
                }
                tensors->write( ofs );
                active_tensors->write( ofs );
                ofs << active_w[0];
                for( int i=1; i<active_tensors->getNumIndexes(); i++ ){
                        ofs << " " << active_w[i];
                }
                ofs << endl;
                if ( points == 0 ){
                        ofs << "0" << endl;
                }else{
                        ofs << "1 ";
                        points->write( ofs );
                }
                if ( needed == 0 ){
                        ofs << "0" << endl;
                }else{
                        ofs << "1 ";
                        needed->write( ofs );
                }
                ofs << max_levels[0];
                for( int j=1; j<num_dimensions; j++ ){
                        ofs << " " << max_levels[j];
                }
                ofs << endl;
                if ( num_outputs > 0 ) values->write( ofs );
                if ( updated_tensors != 0 ){
                        ofs << "1" << endl;
                        updated_tensors->write( ofs );
                        updated_active_tensors->write( ofs );
                        ofs << updated_active_w[0];
                        for( int i=1; i<updated_active_tensors->getNumIndexes(); i++ ){
                                ofs << " " << updated_active_w[i];
                        }
                }else{
                        ofs << "0";
                }
                ofs << endl;
        }
}
void GridGlobal::read( std::ifstream &ifs ){
        reset( true ); // true deletes any custom rule
        int flag;
        ifs >> num_dimensions >> num_outputs >> alpha >> beta;
        if ( num_dimensions > 0 ){
                std::string T;
                ifs >> T;
                rule = OneDimensionalMeta::getIORuleString( T.c_str() );
                if ( rule == rule_customtabulated ){
                        custom = new CustomTabulated();
                        custom->read( ifs );
                }
                tensors = new IndexSet( num_dimensions );  tensors->read( ifs );
                active_tensors = new IndexSet( num_dimensions );  active_tensors->read( ifs );
                active_w = new int[ active_tensors->getNumIndexes() ];  for( int i=0; i<active_tensors->getNumIndexes(); i++ ){  ifs >> active_w[i]; }
                ifs >> flag; if ( flag == 1 ){  points = new IndexSet(num_dimensions); points->read( ifs );  }
                ifs >> flag; if ( flag == 1 ){  needed = new IndexSet(num_dimensions); needed->read( ifs );  }
                max_levels = new int[ num_dimensions ];  for( int j=0; j<num_dimensions; j++ ){  ifs >> max_levels[j]; }
                if ( num_outputs > 0 ){  values = new StorageSet( 0, 0 ); values->read( ifs ); }
                ifs >> flag;
                IndexManipulator IM(num_dimensions);
                int oned_max_level;
                if ( flag == 1 ){
                        updated_tensors = new IndexSet( num_dimensions );  updated_tensors->read( ifs );
                        updated_active_tensors = new IndexSet( num_dimensions );  updated_active_tensors->read( ifs );
                        updated_active_w = new int[ updated_active_tensors->getNumIndexes() ];  for( int i=0; i<updated_active_tensors->getNumIndexes(); i++ ){  ifs >> updated_active_w[i]; }
                        IM.getMaxLevels( updated_tensors, 0, oned_max_level );
                }else{
                        oned_max_level = max_levels[0];
                        for( int j=1; j<num_dimensions; j++ ) if ( oned_max_level < max_levels[j] ) oned_max_level = max_levels[j];
                }
                OneDimensionalMeta meta( custom );
                wrapper = new OneDimensionalWrapper( &meta, oned_max_level, rule, alpha, beta );

                int nz_weights = active_tensors->getNumIndexes();
                IndexSet *work = (points != 0 ) ? points : needed;
                if ( OneDimensionalMeta::isNonNested( rule ) ){
                        tensor_refs = new int*[ nz_weights ];
                        #pragma omp parallel for schedule(dynamic)
                        for( int i=0; i<nz_weights; i++ ){
                                tensor_refs[i] = IM.referenceGenericPoints( active_tensors->getIndex(i), wrapper, work );
                        }
                }else{
                        tensor_refs = new int*[ nz_weights ];
                        #pragma omp parallel for schedule(dynamic)
                        for( int i=0; i<nz_weights; i++ ){
                                tensor_refs[i] = IM.referenceNestedPoints( active_tensors->getIndex(i), wrapper, work );
                        }
                }
                work = 0;
        }
}

void GridGlobal::reset( bool includeCustom ){
        if ( tensor_refs != 0 ){ for( int i=0; i<active_tensors->getNumIndexes(); i++ ){  delete[] tensor_refs[i];  } delete[] tensor_refs;  tensor_refs = 0; }
        if ( wrapper != 0 ){  delete wrapper;  wrapper = 0; }
        if ( tensors != 0 ){  delete tensors;  tensors = 0; }
        if ( active_tensors != 0 ){  delete active_tensors;  active_tensors = 0; }
        if ( active_w != 0 ){  delete[] active_w;  active_w = 0; }
        if ( points != 0 ){  delete points;  points = 0; }
        if ( needed != 0 ){  delete needed;  needed = 0; }
        if ( max_levels != 0 ){  delete[] max_levels;  max_levels = 0; }
        if ( values != 0 ){  delete values;  values = 0; }
        if ( updated_tensors != 0 ){  delete updated_tensors;  updated_tensors = 0; }
        if ( updated_active_tensors != 0 ){  delete updated_active_tensors;  updated_active_tensors = 0; }
        if ( updated_active_w != 0 ){  delete[] updated_active_w;  updated_active_w = 0; }
        if ( includeCustom && (custom != 0) ){  delete custom;  custom = 0; }
        num_dimensions = num_outputs = 0;
}

void GridGlobal::clearRefinement(){
        if ( needed != 0 ){  delete needed;  needed = 0; }
        if ( updated_tensors != 0 ){  delete updated_tensors;  updated_tensors = 0; }
        if ( updated_active_tensors != 0 ){  delete updated_active_tensors;  updated_active_tensors = 0; }
        if ( updated_active_w != 0 ){  delete[] updated_active_w;  updated_active_w = 0; }
}

void GridGlobal::makeGrid( int cnum_dimensions, int cnum_outputs, int depth, TypeDepth type, TypeOneDRule crule, const int *anisotropic_weights, double calpha, double cbeta, const char* custom_filename ){
        if ( crule == rule_customtabulated ){
                custom = new CustomTabulated( custom_filename );
        }
        IndexManipulator IM(cnum_dimensions, custom);

        IndexSet *tset = IM.selectTensors( depth, type, anisotropic_weights, crule );

        setTensors( tset, cnum_outputs, crule, calpha, cbeta );
}

void GridGlobal::setTensors( IndexSet* &tset, int cnum_outputs, TypeOneDRule crule, double calpha, double cbeta ){
        reset( false );
        num_dimensions = tset->getNumDimensions();
        num_outputs = cnum_outputs;
        rule = crule;
        alpha = calpha;  beta = cbeta;

        tensors = tset;
        tset = 0;

        IndexManipulator IM(num_dimensions, custom);

        OneDimensionalMeta meta( custom );
        max_levels = new int[num_dimensions];
        int max_level; IM.getMaxLevels( tensors, max_levels, max_level );
        wrapper = new OneDimensionalWrapper( &meta, max_level, rule, alpha, beta );


        int* tensors_w = IM.makeTensorWeights( tensors );
        active_tensors = IM.nonzeroSubset( tensors, tensors_w );

        int nz_weights = active_tensors->getNumIndexes();

        active_w = new int[ nz_weights ];
        int count = 0;
        for( int i=0; i<tensors->getNumIndexes(); i++ ){  if ( tensors_w[i] != 0 ) active_w[count++] = tensors_w[i];  }

        delete[] tensors_w;

        if ( OneDimensionalMeta::isNonNested( rule ) ){
                needed = IM.generateGenericPoints( active_tensors, wrapper );

                tensor_refs = new int*[ nz_weights ];
                #pragma omp parallel for schedule(dynamic)
                for( int i=0; i<nz_weights; i++ ){
                        tensor_refs[i] = IM.referenceGenericPoints( active_tensors->getIndex(i), wrapper, needed );
                }
        }else{
                needed = IM.generateNestedPoints( tensors, wrapper ); // nested grids exploit nesting

                tensor_refs = new int*[ nz_weights ];
                #pragma omp parallel for schedule(dynamic)
                for( int i=0; i<nz_weights; i++ ){
                        tensor_refs[i] = IM.referenceNestedPoints( active_tensors->getIndex(i), wrapper, needed );
                }
        }

        if ( num_outputs == 0 ){
                points = needed;
                needed = 0;
        }else{
                values = new StorageSet( num_outputs, needed->getNumIndexes() );
        }
}

void GridGlobal::updateGrid( int depth, TypeDepth type, const int *anisotropic_weights ){
        if ( (num_outputs == 0) || (points == 0) ){
                makeGrid( num_dimensions, num_outputs, depth, type, rule, anisotropic_weights, alpha, beta );
        }else{
                clearRefinement();

                IndexManipulator IM(num_dimensions, custom);

                updated_tensors = IM.selectTensors( depth, type, anisotropic_weights, rule );

                needed = updated_tensors->diffSets( tensors );
                if ( (needed != 0) && (needed->getNumIndexes() > 0) ){
                        delete needed;

                        updated_tensors->addIndexSet( tensors ); // avoids the case where existing points in tensor are not included in the update

                        OneDimensionalMeta meta( custom );
                        int max_level; IM.getMaxLevels( updated_tensors, 0, max_level );
                        delete wrapper;
                        wrapper = new OneDimensionalWrapper( &meta, max_level, rule, alpha, beta );

                        int* updates_tensor_w = IM.makeTensorWeights( updated_tensors );
                        updated_active_tensors = IM.nonzeroSubset( updated_tensors, updates_tensor_w );

                        int nz_weights = updated_active_tensors->getNumIndexes();

                        updated_active_w = new int[ nz_weights ];
                        int count = 0;
                        for( int i=0; i<updated_tensors->getNumIndexes(); i++ ){  if ( updates_tensor_w[i] != 0 ) updated_active_w[count++] = updates_tensor_w[i];  }

                        IndexSet *new_points = 0;
                        if ( OneDimensionalMeta::isNonNested( rule ) ){
                                new_points = IM.generateGenericPoints( updated_active_tensors, wrapper );
                        }else{
                                new_points = IM.generateNestedPoints( updated_tensors, wrapper );
                        }

                        needed = new_points->diffSets( points );

                        delete new_points;
                        delete[] updates_tensor_w;
                }else{
                        clearRefinement();
                }
        }
}

int GridGlobal::getNumDimensions() const{ return num_dimensions; }
int GridGlobal::getNumOutputs() const{ return num_outputs; }
TypeOneDRule GridGlobal::getRule() const{ return rule; }
const char* GridGlobal::getCustomRuleDescription() const{  return ( custom != 0 ) ? custom->getDescription() : "";  }

double GridGlobal::getAlpha() const{ return alpha; }
double GridGlobal::getBeta() const{ return beta; }

int GridGlobal::getNumLoaded() const{ return ((points == 0) ? 0 : points->getNumIndexes()); }
int GridGlobal::getNumNeeded() const{ return ((needed == 0) ? 0 : needed->getNumIndexes()); }
int GridGlobal::getNumPoints() const{ return ((points == 0) ? getNumNeeded() : getNumLoaded()); }

double* GridGlobal::getLoadedPoints() const{
        if ( points == 0 ) return 0;
        int num_points = points->getNumIndexes();
        double *x = new double[ num_dimensions * num_points ];
        #pragma omp parallel for schedule(static)
        for( int i=0; i<num_points; i++ ){
                const int *p = points->getIndex( i );
                for( int j=0; j<num_dimensions; j++ ){
                        x[i*num_dimensions + j] = wrapper->getNode( p[j] );
                }
        }
        return x;
}
double* GridGlobal::getNeededPoints() const{
        if ( needed == 0 ) return 0;
        int num_points = needed->getNumIndexes();
        if ( num_points == 0 ) return 0;
        double *x = new double[ num_dimensions * num_points ];
        #pragma omp parallel for schedule(static)
        for( int i=0; i<num_points; i++ ){
                const int *p = needed->getIndex( i );
                for( int j=0; j<num_dimensions; j++ ){
                        x[i*num_dimensions + j] = wrapper->getNode( p[j] );
                }
        }
        return x;
}

double* GridGlobal::getPoints() const{
        return (( points == 0 ) ? getNeededPoints() : getLoadedPoints());
}

double* GridGlobal::getQuadratureWeights() const{
        IndexSet *work = ( points == 0 ) ? needed : points;

        int num_points = work->getNumIndexes();
        double *weights = new double[ num_points ];
        std::fill( weights, weights + num_points, 0.0 );

        int *num_oned_points = new int[num_dimensions];
        for( int n=0; n<active_tensors->getNumIndexes(); n++ ){
                const int* levels = active_tensors->getIndex( n );
                num_oned_points[0] = wrapper->getNumPoints( levels[0] );
                int num_tensor_points = num_oned_points[0];
                for( int j=1; j<num_dimensions; j++ ){
                        num_oned_points[j] = wrapper->getNumPoints( levels[j] );
                        num_tensor_points *= num_oned_points[j];
                }
                double tensor_weight = (double) active_w[n];

                #pragma omp parallel for schedule(static)
                for( int i=0; i<num_tensor_points; i++ ){
                        int t = i;
                        double w = 1.0;
                        for( int j=num_dimensions-1; j>=0; j-- ){
                                w *= wrapper->getWeight( levels[j], t % num_oned_points[j] );
                                t /= num_oned_points[j];
                        }
                        weights[ tensor_refs[n][i] ] += tensor_weight * w;
                }
        }
        delete[] num_oned_points;
        work = 0;

        return weights;
}

double* GridGlobal::getInterpolationWeights( const double x[] ) const{
        IndexSet *work = ( points == 0 ) ? needed : points;

        CacheLagrange *lcache = new CacheLagrange( num_dimensions, max_levels, wrapper, x );

        int num_points = work->getNumIndexes();
        double *weights = new double[ num_points ];
        std::fill( weights, weights + num_points, 0.0 );

        int *num_oned_points = new int[num_dimensions];
        for( int n=0; n<active_tensors->getNumIndexes(); n++ ){
                const int* levels = active_tensors->getIndex( n );
                num_oned_points[0] = wrapper->getNumPoints( levels[0] );
                int num_tensor_points = num_oned_points[0];
                for( int j=1; j<num_dimensions; j++ ){
                        num_oned_points[j] = wrapper->getNumPoints( levels[j] );
                        num_tensor_points *= num_oned_points[j];
                }
                double tensor_weight = (double) active_w[n];
                for( int i=0; i<num_tensor_points; i++ ){
                        int t = i;
                        double w = 1.0;
                        for( int j=num_dimensions-1; j>=0; j-- ){
                                w *= lcache->getLagrange( j, levels[j], t % num_oned_points[j] );
                                t /= num_oned_points[j];
                        }
                        weights[ tensor_refs[n][i] ] += tensor_weight * w;
                }
        }
        delete[] num_oned_points;
        work = 0;

        delete lcache;

        return weights;
}

void GridGlobal::loadNeededPoints( const double *vals ){
        if ( points == 0 ){
                values->setValues( vals );
                points = needed;
                needed = 0;
        }else{
                values->addValues( points, needed, vals );
                points->addIndexSet( needed );
                delete needed; needed = 0;
                for( int i=0; i<active_tensors->getNumIndexes(); i++ ){ delete[] tensor_refs[i]; } delete[] tensor_refs;

                delete tensors; tensors = updated_tensors; updated_tensors = 0;
                delete active_tensors; active_tensors = updated_active_tensors; updated_active_tensors = 0;
                delete[] active_w; active_w = updated_active_w; updated_active_w = 0;

                int nz_weights = active_tensors->getNumIndexes();
                IndexManipulator IM(num_dimensions, custom);
                int m; IM.getMaxLevels( tensors, max_levels, m );
                if ( OneDimensionalMeta::isNonNested( rule ) ){
                        tensor_refs = new int*[ nz_weights ];
                        #pragma omp parallel for schedule(dynamic)
                        for( int i=0; i<nz_weights; i++ ){
                                tensor_refs[i] = IM.referenceGenericPoints( active_tensors->getIndex(i), wrapper, points );
                        }
                }else{
                        tensor_refs = new int*[ nz_weights ];
                        #pragma omp parallel for schedule(dynamic)
                        for( int i=0; i<nz_weights; i++ ){
                                tensor_refs[i] = IM.referenceNestedPoints( active_tensors->getIndex(i), wrapper, points );
                        }
                }
        }
}

void GridGlobal::evaluate( const double x[], double y[] ) const{
        std::fill( y, y + num_outputs, 0.0 );

        double *w = getInterpolationWeights( x );
        std::fill( y, y+num_outputs, 0.0 );
        for( int k=0; k<num_outputs; k++ ){
                for( int i=0; i<points->getNumIndexes(); i++ ){
                        const double *v = values->getValues( i );
                        y[k] += w[i] * v[k];
                }
        }
        delete[] w;
}

void GridGlobal::integrate( double q[] ) const{
        double *w = getQuadratureWeights();
        std::fill( q, q+num_outputs, 0.0 );
        #pragma omp parallel for schedule(static)
        for( int k=0; k<num_outputs; k++ ){
                for( int i=0; i<points->getNumIndexes(); i++ ){
                        const double *v = values->getValues( i );
                        q[k] += w[i] * v[k];
                }
        }
        delete[] w;
}

double* GridGlobal::computeSurpluses( int output, bool normalize ) const{
        int n = points->getNumIndexes();
        double *surp = new double[ n ];
        double max_surp = 0.0;

        if ( OneDimensionalMeta::isSequence( rule ) ){
                for( int i=0; i<n; i++ ){
                        const double* v = values->getValues( i );
                        surp[i] = v[output];
                        if ( fabs( surp[i] ) > max_surp ) max_surp = fabs( surp[i] );
                }

                IndexManipulator IM(num_dimensions);
                int *level = IM.computeLevels( points );
                int top_level = level[0];  for( int i=1; i<n; i++ ){  if ( top_level < level[i] ) top_level = level[i];  }
                int top_1d = 0; const int* id = points->getIndex( 0 );  for( int i=0; i<n*num_dimensions; i++ ) if ( top_1d < id[i] ) top_1d = id[i];

                int *parents = IM.computeDAGup( points );

                const double* nodes = wrapper->getNodes( 0 );
                double *coeff = new double[ top_1d+1 ];
                coeff[0] = 1.0;
                for( int i=1; i<=top_1d; i++ ){
                        coeff[i] = 1.0;
                        for( int j=0; j<i; j++ ) coeff[i] *= ( nodes[i] - nodes[j] );
                }

                for( int l=1; l<=top_level; l++ ){
                        #pragma omp parallel for schedule(dynamic)
                        for( int i=0; i<n; i++ ){
                                if ( level[i] == l ){
                                        const int* p = points->getIndex( i );

                                        int *monkey_count = new int[top_level + 1 ];
                                        int *monkey_tail = new int[top_level + 1 ];
                                        bool *used = new bool[n];
                                        std::fill( used, used + n, false );

                                        int current = 0;

                                        monkey_count[0] = 0;
                                        monkey_tail[0] = i;

                                        while( monkey_count[0] < num_dimensions ){
                                                if ( monkey_count[current] < num_dimensions ){
                                                        int branch = parents[ monkey_tail[current] * num_dimensions + monkey_count[current] ];
                                                        if ( (branch == -1) || (used[branch]) ){
                                                                monkey_count[current]++;
                                                        }else{
                                                                double basis_value = 1.0;
                                                                const int* f = points->getIndex(branch);
                                                                for( int j=0; j<num_dimensions; j++ ){
                                                                        double x = nodes[ p[j] ];
                                                                        double w = 1.0;
                                                                        for( int k=0; k<f[j]; k++ ){
                                                                                w *= ( x - nodes[k] );
                                                                        }
                                                                        basis_value *= w / coeff[ f[j] ];
                                                                }
                                                                surp[ i ] -= basis_value * surp[ branch ];

                                                                used[branch] = true;

                                                                monkey_count[++current] = 0;
                                                                monkey_tail[current] = branch;
                                                        }
                                                }else{
                                                        monkey_count[--current]++;
                                                }
                                        }

                                        delete[] used;
                                        delete[] monkey_tail;
                                        delete[] monkey_count;
                                }
                        }
                }

                delete[] coeff;
                delete[] parents;
                delete[] level;

                if ( normalize ){
                        #pragma omp parallel for schedule(static)
                        for( int i=0; i<n; i++ ) surp[i] /= max_surp;
                }
        }else{
                IndexManipulator IM(num_dimensions);
                int *max_lvl = new int[num_dimensions];
                int ml; IM.getMaxLevels( points, max_lvl, ml );

                GridGlobal *gg = new GridGlobal();
                gg->makeGrid(num_dimensions, 0, 2*ml, type_qptotal, rule_gausspatterson );

                int qn = gg->getNumPoints();
                double *w = gg->getQuadratureWeights();
                double *x = gg->getPoints();
                double *I = new double[qn];
                delete gg;
                #pragma omp parallel for schedule(static)
                for( int i=0; i<qn; i++ ){
                        double *y = new double[num_outputs];
                        evaluate( &( x[i*num_dimensions]), y );
                        I[i] = w[i] * y[output];
                        delete[] y;
                }

                #pragma omp parallel for schedule(dynamic)
                for( int i=0; i<n; i++ ){
                        // for each surp, do the quadrature
                        const int* p = points->getIndex( i );
                        double c = 0.0;
                        for( int k=0; k<qn; k++ ){
                                double v = legendre( p[0], x[k*num_dimensions] );
                                for( int j=1; j<num_dimensions; j++ ){
                                        v *= legendre( p[j], x[k*num_dimensions+j] );
                                }
                                c += v * I[k];
                        }
                        double nrm = sqrt( (double) p[0] + 0.5 );
                        for( int j=1; j<num_dimensions; j++ ){
                                nrm *= sqrt( (double) p[j] + 0.5 );
                        }
                        surp[i] = c * nrm;
                }

                delete[] max_lvl;
                delete[] I;
                delete[] x;
                delete[] w;
        }

        return surp;
}

int* GridGlobal::estimateAnisotropicCoefficients( TypeDepth type, int output ) const{
        double tol = 1000 * TSG_NUM_TOL;
        double *surp = computeSurpluses( output, false );

        int num_points = points->getNumIndexes();

        int n = 0, m;
        for( int j=0; j<num_points; j++ ){
                surp[j] = fabs( surp[j] );
                if ( surp[j] > tol ) n++;
        }

        double *A, *b;

        if ( (type == type_curved) || (type == type_ipcurved) || (type == type_qpcurved) ){
                m = 2*num_dimensions + 1;
                A = new double[ n * m ];
                b = new double[ n ];

                int count = 0;
                for( int c=0; c<num_points; c++ ){
                        const int *indx = points->getIndex(c);
                        if ( surp[c] > tol ){
                                for( int j=0; j<num_dimensions; j++ ){
                                        A[ j*n + count ] = ( (double) indx[j] );
                                }
                                for( int j=0; j<num_dimensions; j++ ){
                                        A[ (num_dimensions + j)*n + count ] = log( (double) (indx[j] + 1 ) );
                                }
                                A[ 2*num_dimensions*n + count ] = 1.0;
                                b[count++] = -log( surp[c] );
                        }
                }
        }else{
                m = num_dimensions + 1;
                A = new double[ n * m ];
                b = new double[ n ];

                int count = 0;
                for( int c=0; c<num_points; c++ ){
                        const int *indx = points->getIndex(c);
                        if ( surp[c] > tol ){
                                for( int j=0; j<num_dimensions; j++ ){
                                        A[ j*n + count ] = ( (double) indx[j] );
                                }
                                A[ num_dimensions*n + count ] = 1.0;
                                b[count++] = - log( surp[c] );
                        }
                }
        }
        delete[] surp;

        double *x = new double[m];
        TasmanianDenseSolver::solveLeastSquares( n, m, A, b, 1.E-5, x );

        int *weights = new int[m--];
        for( int j=0; j<m; j++ ){
                weights[j] = (int)( x[j] * 1000.0 + 0.5 );
        }

        int min_weight = weights[0];  for( int j=1; j<num_dimensions; j++ ){  if ( min_weight < weights[j] ) min_weight = weights[j];  } // start by finding the largest weight
        if ( min_weight < 0 ){ // all directions are diverging, default to isotropic total degree
                for( int j=0; j<num_dimensions; j++ )  weights[j] = 1;
                if ( m == 2*num_dimensions ) for( int j=num_dimensions; j<2*num_dimensions; j++ )  weights[j] = 0;
        }else{
                for( int j=0; j<num_dimensions; j++ )  if( (weights[j]>0) && (weights[j]<min_weight) ) min_weight = weights[j];  // find the smallest positive weight
                for( int j=0; j<num_dimensions; j++ ){
                        if ( weights[j] <= 0 ){
                                weights[j] = min_weight;
                                if ( m == 2*num_dimensions ){
                                        if ( fabs( weights[num_dimensions + j] ) > weights[j] ){
                                                weights[num_dimensions + j] = ( weights[num_dimensions + j] > 0.0 ) ? weights[j] : -weights[j];
                                        }
                                }
                        }
                }
        }

        delete[] A;
        delete[] b;
        delete[] x;

        return weights;
}

void GridGlobal::setAnisotropicRefinement( TypeDepth type, int min_growth, int output ){
        clearRefinement();
        int *weights = estimateAnisotropicCoefficients( type, output ); //for( int i=0; i<2*num_dimensions; i++ ) cout << weights[i] << "  "; cout << endl;

        IndexManipulator IM(num_dimensions);
        int level = IM.getMinChildLevel( tensors, type, weights, rule );

        updated_tensors = IM.selectTensors( level, type, weights, rule );
        needed = updated_tensors->diffSets( tensors ); // this exploits the 1-1 correspondence between points and tensors for sequence rules (see the correction below)

        while( (needed == 0) || (needed->getNumIndexes() < min_growth) ){ // for CC min_growth is lots of points
                delete updated_tensors;
                if ( needed != 0 ) delete needed;
                updated_tensors = IM.selectTensors( ++level, type, weights, rule );
                needed = updated_tensors->diffSets( tensors );
        }
        delete[] weights;

        updated_tensors->addIndexSet( tensors ); // avoids the case where existing points in tensor are not included in the update

        OneDimensionalMeta meta( custom );
        int max_level; IM.getMaxLevels( updated_tensors, 0, max_level );
        delete wrapper;
        wrapper = new OneDimensionalWrapper( &meta, max_level, rule, alpha, beta );

        int* updates_tensor_w = IM.makeTensorWeights( updated_tensors );
        updated_active_tensors = IM.nonzeroSubset( updated_tensors, updates_tensor_w );

        int nz_weights = updated_active_tensors->getNumIndexes();

        updated_active_w = new int[ nz_weights ];
        int count = 0;
        for( int i=0; i<updated_tensors->getNumIndexes(); i++ ){  if ( updates_tensor_w[i] != 0 ) updated_active_w[count++] = updates_tensor_w[i];  }

        if ( !( OneDimensionalMeta::isSequence( rule ) ) ){ // non-sequence rules need to recompute needed
                delete needed;

                IndexSet *new_points = IM.generateNestedPoints( updated_tensors, wrapper );

                needed = new_points->diffSets( points );

                delete new_points;
        }
        delete[] updates_tensor_w;
}

void GridGlobal::setSurplusRefinement( double tolerance, int output ){
        clearRefinement();
        double *surp = computeSurpluses( output, true );

        int n = points->getNumIndexes();
        bool *flagged = new bool[n];

        #pragma omp parallel for
        for( int i=0; i<n; i++ ){
                flagged[i] = ( fabs( surp[i] ) > tolerance );
        }

        IndexManipulator IM(num_dimensions);
        IndexSet *kids = IM.selectFlaggedChildren( points, flagged );

        if ( (kids != 0) && (kids->getNumIndexes() > 0) ){
                kids->addIndexSet( points );

                updated_tensors = IM.getLowerCompletion( kids );
                if ( updated_tensors == 0 ){
                        updated_tensors = kids;
                        kids = 0;
                }else{
                        updated_tensors->addIndexSet( kids );
                        delete kids;
                }

                OneDimensionalMeta meta( custom );
                int max_level; IM.getMaxLevels( updated_tensors, 0, max_level );
                delete wrapper;
                wrapper = new OneDimensionalWrapper( &meta, max_level, rule, alpha, beta );

                int* updates_tensor_w = IM.makeTensorWeights( updated_tensors );
                updated_active_tensors = IM.nonzeroSubset( updated_tensors, updates_tensor_w );

                int nz_weights = updated_active_tensors->getNumIndexes();

                updated_active_w = new int[ nz_weights ];
                int count = 0;
                for( int i=0; i<updated_tensors->getNumIndexes(); i++ ){  if ( updates_tensor_w[i] != 0 ) updated_active_w[count++] = updates_tensor_w[i];  }

                delete[] updates_tensor_w;

                needed = updated_tensors->diffSets( tensors );
        }

        delete[] flagged;
        delete[] surp;
}

double GridGlobal::legendre( int n, double x ){
        if ( n == 0 ) return 1.0;
        if ( n == 1 ) return x;

        double lm = 1, l = x, lp;
        for( int i=2; i<=n; i++ ){
                lp = ( ( 2*i - 1 ) * x * l ) / ( (double) i ) - ( ( i - 1 ) * lm ) / ( (double) i );
                lm = l; l = lp;
        }
        return l;
}
double* GridGlobal::multi_legendre( int n, double x ){
        double *l = new double[n];
        l[0] = 1.0;
        if ( n > 1 ) l[1] = x;

        for( int i=2; i<n; i++ ){
                l[i] = ( ( 2*i - 1 ) * x * l[i-1] ) / ( (double) i ) - ( ( i - 1 ) * l[i-2] ) / ( (double) i );
        }
        for( int i=0; i<n; i++ ) l[i] *= sqrt( ( (double) i ) + 0.5 );

        return l;
}

int* GridGlobal::getPolynomialSpace( bool interpolation, int &n ) const{
        IndexManipulator IM( num_dimensions, custom );
        IndexSet* set = IM.getPolynomialSpace( active_tensors, rule, interpolation );

        n = set->getNumIndexes();

        int *poly = new int[ n * num_dimensions ];
        const int* p = set->getIndex( 0 );

        std::copy( p, p + n * num_dimensions, poly );

        delete set;

        return poly;
}

}

#endif