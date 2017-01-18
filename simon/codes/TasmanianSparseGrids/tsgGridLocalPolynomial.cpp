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

#ifndef __TASMANIAN_SPARSE_GRID_LPOLY_CPP
#define __TASMANIAN_SPARSE_GRID_LPOLY_CPP

#include "tsgGridLocalPolynomial.hpp"

namespace TasGrid{

GridLocalPolynomial::GridLocalPolynomial() : num_dimensions(0), num_outputs(0), order(1), top_level(0),
                                             surpluses(0), points(0), needed(0), values(0), parents(0), num_roots(0), roots(0), pntr(0), indx(0), rule(0)  {}
GridLocalPolynomial::~GridLocalPolynomial(){
        reset();
}
void GridLocalPolynomial::reset(){
        if ( surpluses != 0 ){  delete[] surpluses; surpluses = 0;  }
        if ( points != 0 ){  delete points; points = 0;  }
        if ( needed != 0 ){  delete needed; needed = 0;  }
        if ( values != 0 ){  delete values; values = 0;  }
        if ( parents != 0 ){  delete[] parents; parents = 0;  }
        if ( roots != 0 ){  delete[] roots;  roots = 0;  }
        if ( pntr != 0 ){  delete[] pntr;  pntr = 0;  }
        if ( indx != 0 ){  delete[] indx;  indx = 0;  }
        rule = 0;
}

void GridLocalPolynomial::write( std::ofstream &ofs ) const{
        ofs << std::scientific; ofs.precision(17);
        ofs << num_dimensions << " " << num_outputs << " " << order << " " << top_level << endl;
        if ( num_dimensions > 0 ){
                ofs << OneDimensionalMeta::getIORuleString( rule->getType() ) << endl;
                if ( points == 0 ){
                        ofs << "0" << endl;
                }else{
                        ofs << "1 ";
                        points->write( ofs );
                }
                if ( surpluses == 0 ){
                        ofs << "0" << endl;
                }else{
                        ofs << "1 ";
                        for( int i=0; i<points->getNumIndexes() * num_outputs; i++ ){  ofs << " " << surpluses[i];  } ofs << endl;
                }
                if ( needed == 0 ){
                        ofs << "0" << endl;
                }else{
                        ofs << "1 ";
                        needed->write( ofs );
                }
                if ( parents == 0 ){
                        ofs << "0" << endl;
                }else{
                        ofs << "1";
                        if ( rule->isSemiLocal() ){
                                for( int i=0; i<2*(points->getNumIndexes()); i++ ){ ofs << " " << parents[i]; } ofs << endl;
                        }else{
                                for( int i=0; i<points->getNumIndexes(); i++ ){ ofs << " " << parents[i]; } ofs << endl;
                        }
                }
                int num_points = ( points == 0 ) ? needed->getNumIndexes() : points->getNumIndexes();
                ofs << num_roots; for( int i=0; i<num_roots; i++ ){  ofs << " " << roots[i];  } ofs << endl;
                ofs << pntr[0]; for( int i=1; i<=num_points; i++ ){  ofs << " " << pntr[i];  }  ofs << endl;
                ofs << indx[0]; for( int i=1; i<pntr[num_points]; i++ ){  ofs << " " << indx[i];  }  ofs << endl;
                if ( num_outputs > 0 ) values->write( ofs );
        }
}
void GridLocalPolynomial::read( std::ifstream &ifs ){
        reset();
        int flag;
        ifs >> num_dimensions >> num_outputs >> order >> top_level;
        if ( num_dimensions > 0 ){
                std::string T;
                ifs >> T;
                TypeOneDRule crule = OneDimensionalMeta::getIORuleString( T.c_str() );
                if ( crule == rule_localp ){
                        rule = &rpoly;
                }else if ( crule == rule_semilocalp ){
                        rule = &rsemipoly;
                }else if ( crule == rule_localp0 ){
                        rule = &rpoly0;
                }
                rule->setMaxOrder( order );

                ifs >> flag;  if ( flag == 1 ){  points    = new IndexSet( num_dimensions );   points->read( ifs );  }
                ifs >> flag;  if ( flag == 1 ){  surpluses = new double[ points->getNumIndexes() * num_outputs ]; for( int i=0; i<points->getNumIndexes() * num_outputs; i++ ){  ifs >> surpluses[i];  }  }
                ifs >> flag;  if ( flag == 1 ){  needed    = new IndexSet( num_dimensions );   needed->read( ifs );  }
                ifs >> flag;  if ( flag == 1 ){  int num_parents = (( rule->isSemiLocal() ) ? 2 : 1 ) * points->getNumIndexes();  parents = new int[num_parents];  for( int i=0; i<num_parents; i++ ){  ifs >> parents[i];  }  }

                int num_points = ( points == 0 ) ? needed->getNumIndexes() : points->getNumIndexes();

                ifs >> num_roots;
                roots = new int[num_roots];         for( int i=0; i<num_roots; i++ ) ifs >> roots[i];
                pntr  = new int[num_points + 1];    for( int i=0; i<=num_points; i++ ) ifs >> pntr[i];
                if ( pntr[num_points] > 0 ){
                        indx  = new int[pntr[num_points]];  for( int i=0; i<pntr[num_points]; i++ ) ifs >> indx[i];
                }else{
                        indx  = new int[1];  ifs >> indx[0]; // there is a special case when the grid has only one point without any children
                }

                if ( num_outputs > 0 ){  values = new StorageSet( 0, 0 ); values->read( ifs );  }
        }
}

void GridLocalPolynomial::makeGrid( int cnum_dimensions, int cnum_outputs, int depth, int corder, TypeOneDRule crule ){
        reset();
        num_dimensions = cnum_dimensions;
        num_outputs = cnum_outputs;
        order = corder;

        TypeOneDRule effective_rule = ( crule == rule_localp0 ) ? rule_localp0 : ( ( (crule == rule_semilocalp) && ( (order == -1) || (order > 1) ) ) ? rule_semilocalp : rule_localp );
        if ( effective_rule == rule_localp ){
                rule = &rpoly;
        }else if ( effective_rule == rule_semilocalp ){
                rule = &rsemipoly;
        }else if ( effective_rule == rule_localp0 ){
                rule = &rpoly0;
        }
        rule->setMaxOrder( order );

        IndexManipulator IM(num_dimensions);
        UnsortedIndexSet* deltas = IM.getToalDegreeDeltas( depth );

        needed = IM.generatePointsFromDeltas( deltas, rule );
        delete deltas;

        buildTree();

        if ( num_outputs == 0 ){
                points = needed;
                needed = 0;
                parents = IM.computeDAGupLocal( points, rule );
        }else{
                values = new StorageSet( num_outputs, needed->getNumIndexes() );
        }
}

int GridLocalPolynomial::getNumDimensions() const{  return num_dimensions;  }
int GridLocalPolynomial::getNumOutputs() const{  return num_outputs;  }
TypeOneDRule GridLocalPolynomial::getRule() const{  return rule->getType();  }
int GridLocalPolynomial::getOrder() const{  return order;  }

int GridLocalPolynomial::getNumLoaded() const{  return ((points == 0) ? 0 : points->getNumIndexes());  }
int GridLocalPolynomial::getNumNeeded() const{  return ((needed == 0) ? 0 : needed->getNumIndexes());  }
int GridLocalPolynomial::getNumPoints() const{  return ((points == 0) ? getNumNeeded() : getNumLoaded());  }

double* GridLocalPolynomial::getLoadedPoints() const{
        if ( points == 0 ) return 0;
        int num_points = points->getNumIndexes();
        if ( num_points == 0 ) return 0;
        double *x = new double[ num_dimensions * num_points ];
        #pragma omp parallel for schedule(static)
        for( int i=0; i<num_points; i++ ){
                const int *p = points->getIndex( i );
                for( int j=0; j<num_dimensions; j++ ){
                        x[i*num_dimensions + j] = rule->getNode( p[j] );
                }
        }
        return x;
}
double* GridLocalPolynomial::getNeededPoints() const{
        if ( needed == 0 ) return 0;
        int num_points = needed->getNumIndexes();
        if ( num_points == 0 ) return 0;
        double *x = new double[ num_dimensions * num_points ];
        #pragma omp parallel for schedule(static)
        for( int i=0; i<num_points; i++ ){
                const int *p = needed->getIndex( i );
                for( int j=0; j<num_dimensions; j++ ){
                        x[i*num_dimensions + j] = rule->getNode( p[j] );
                }
        }
        return x;
}
double* GridLocalPolynomial::getPoints() const{
        return (( points == 0 ) ? getNeededPoints() : getLoadedPoints());
}

void GridLocalPolynomial::evaluate( const double x[], double y[] ) const{
        int *monkey_count = new int[top_level+1];
        int *monkey_tail = new int[top_level+1];

        double basis_value;
        bool isSupported;
        int offset;

        std::fill( y, y + num_outputs, 0.0 );

        for( int r=0; r<num_roots; r++ ){
                basis_value = evalBasisSupported( points->getIndex( roots[r] ), x, isSupported );

                if ( isSupported ){
                        offset = roots[r] * num_outputs;
                        for( int k=0; k<num_outputs; k++ ) y[k] += basis_value * surpluses[ offset + k ];

                        int current = 0;
                        monkey_tail[0] = roots[r];
                        monkey_count[0] = pntr[roots[r]];

                        while( monkey_count[0] < pntr[ monkey_tail[0]+1 ] ){
                                if ( monkey_count[current] < pntr[ monkey_tail[current]+1 ] ){
                                        offset = indx[ monkey_count[current] ];
                                        basis_value = evalBasisSupported( points->getIndex( offset ), x, isSupported );
                                        if ( isSupported ){
                                                offset *= num_outputs;
                                                for( int k=0; k<num_outputs; k++ ) y[k] += basis_value * surpluses[ offset + k ];
                                                offset /= num_outputs;

                                                monkey_tail[++current] = offset;
                                                monkey_count[current] = pntr[ offset ];
                                        }else{
                                                monkey_count[current]++;
                                        }
                                }else{
                                        monkey_count[--current]++;
                                }
                        }
                }
        }

        delete[] monkey_count;
        delete[] monkey_tail;
}

void GridLocalPolynomial::loadNeededPoints( const double *vals ){
        if ( points == 0 ){
                values->setValues( vals );
                points = needed;
                needed = 0;
        }else{
                values->addValues( points, needed, vals );
                points->addIndexSet( needed );
                delete needed; needed = 0;
                buildTree();
        }
        recomputeSurpluses();
}

double* GridLocalPolynomial::getInterpolationWeights( const double x[] ) const{
        IndexSet *work = ( points == 0 ) ? needed : points;

        int *active_points = new int[ work->getNumIndexes() ];
        double *weights = new double[ work->getNumIndexes() ];   std::fill( weights, weights + work->getNumIndexes(), 0.0 );
        int num_active_points = 0;

        int *monkey_count = new int[top_level+1];
        int *monkey_tail = new int[top_level+1];

        double basis_value;
        bool isSupported;
        int offset;

        for( int r=0; r<num_roots; r++ ){

                basis_value = evalBasisSupported( work->getIndex( roots[r] ), x, isSupported );

                if ( isSupported ){
                        active_points[num_active_points++] = roots[r];
                        weights[roots[r]] = basis_value;

                        int current = 0;
                        monkey_tail[0] = roots[r];
                        monkey_count[0] = pntr[roots[r]];

                        while( monkey_count[0] < pntr[ monkey_tail[0]+1 ] ){
                                if ( monkey_count[current] < pntr[ monkey_tail[current]+1 ] ){
                                        offset = indx[ monkey_count[current] ];

                                        basis_value = evalBasisSupported( work->getIndex( offset ), x, isSupported );

                                        if ( isSupported ){
                                                active_points[num_active_points++] = offset;
                                                weights[offset] = basis_value;

                                                monkey_tail[++current] = offset;
                                                monkey_count[current] = pntr[ offset ];
                                        }else{
                                                monkey_count[current]++;
                                        }
                                }else{
                                        monkey_count[--current]++;
                                }
                        }
                }
        }

        // apply the transpose of the surplus transformation
        const int *dagUp;
        if ( num_outputs == 0 ){
                dagUp = parents;
        }else{
                IndexManipulator IM(num_dimensions);
                dagUp = IM.computeDAGupLocal( work, rule );
        }

        int *level = new int[num_active_points];
        int active_top_level = 0;
        for( int i=0; i<num_active_points; i++ ){
                const int *p = work->getIndex( active_points[i] );
                level[i] = rule->getLevel( p[0] );
                for( int j=1; j<num_dimensions; j++ ){
                        level[i] += rule->getLevel( p[j] );
                }
                if ( active_top_level < level[i] ) active_top_level = level[i];
        }

        bool *used = new bool[work->getNumIndexes()];
        double *node = new double[num_dimensions];
        int max_parents = ( rule->isSemiLocal() ) ? 2*num_dimensions : num_dimensions;

        for( int l=active_top_level; l>0; l-- ){
                for( int i=0; i<num_active_points; i++ ){
                        if ( level[i] == l ){
                                const int* p = work->getIndex( active_points[i] );
                                for( int j=0; j<num_dimensions; j++ ) node[j] = rule->getNode( p[j] );

                                std::fill( used, used + work->getNumIndexes(), false );

                                monkey_count[0] = 0;
                                monkey_tail[0] = active_points[i];
                                int current = 0;

                                while( monkey_count[0] < max_parents ){
                                        if ( monkey_count[current] < max_parents ){
                                                int branch = dagUp[  monkey_tail[current] * max_parents + monkey_count[current]  ];
                                                if ( (branch == -1) || used[branch] ){
                                                        monkey_count[current]++;
                                                }else{
                                                        const int *func = work->getIndex( branch );
                                                        basis_value = rule->evalRaw( func[0], node[0] );
                                                        for( int j=1; j<num_dimensions; j++ ) basis_value *= rule->evalRaw( func[j], node[j] );
                                                        weights[branch] -= weights[active_points[i]] * basis_value;
                                                        used[branch] = true;

                                                        monkey_count[++current] = 0;
                                                        monkey_tail[current] = branch;
                                                }
                                        }else{
                                                monkey_count[--current]++;
                                        }
                                }
                        }
                }
        }

        delete[] used;
        delete[] node;

        delete[] level;

        delete[] monkey_count;
        delete[] monkey_tail;

        delete[] active_points;

        if ( num_outputs > 0 ){
                delete[] dagUp;
        }

        return weights;
}

void GridLocalPolynomial::recomputeSurpluses(){
        int num_ponits = points->getNumIndexes();
        if ( surpluses != 0 ) delete[] surpluses;
        surpluses = new double[ num_ponits * num_outputs ];

        if ( num_outputs > 2 ){
                #pragma omp parallel for schedule(static)
                for( int i=0; i<num_ponits; i++ ){
                        const double* v = values->getValues( i );
                        std::copy( v, v + num_outputs, &( surpluses[i*num_outputs] ) );
                }
        }else{
                const double* v = values->getValues( 0 );
                std::copy( v, v + num_ponits * num_outputs, surpluses );
        }

        IndexManipulator IM(num_dimensions);
        int *dagUp = IM.computeDAGupLocal( points, rule );

        int max_parents = ( rule->isSemiLocal() ) ? 2*num_dimensions : num_dimensions;

        int *level = new int[num_ponits];
        #pragma omp parallel for schedule(static)
        for( int i=0; i<num_ponits; i++ ){
                const int *p = points->getIndex( i );
                level[i] = rule->getLevel( p[0] );
                for( int j=1; j<num_dimensions; j++ ){
                        level[i] += rule->getLevel( p[j] );
                }
        }

        for( int l=1; l<=top_level; l++ ){
                #pragma omp parallel for schedule(dynamic)
                for( int i=0; i<num_ponits; i++ ){
                        if ( level[i] == l ){
                                const int* p = points->getIndex( i );
                                double *x = new double[num_dimensions];
                                for( int j=0; j<num_dimensions; j++ ) x[j] = rule->getNode( p[j] );

                                int *monkey_count = new int[top_level + 1 ];
                                int *monkey_tail = new int[top_level + 1 ];
                                bool *used = new bool[num_ponits];
                                std::fill( used, used + num_ponits, false );

                                int current = 0;

                                monkey_count[0] = 0;
                                monkey_tail[0] = i;

                                while( monkey_count[0] < max_parents ){
                                        if ( monkey_count[current] < max_parents ){
                                                int branch = dagUp[ monkey_tail[current] * max_parents + monkey_count[current] ];
                                                if ( (branch == -1) || (used[branch]) ){
                                                        monkey_count[current]++;
                                                }else{
                                                        double basis_value = evalBasisRaw( points->getIndex(branch), x );
                                                        for( int k=0; k<num_outputs; k++ ){
                                                                surpluses[ i*num_outputs + k ] -= basis_value * surpluses[ branch * num_outputs + k ];
                                                        }
                                                        used[branch] = true;

                                                        monkey_count[++current] = 0;
                                                        monkey_tail[current] = branch;
                                                }
                                        }else{
                                                monkey_count[--current]++;
                                        }
                                }

                                delete[] x;
                                delete[] used;
                                delete[] monkey_tail;
                                delete[] monkey_count;
                        }
                }
        }

        delete[] dagUp;
        delete[] level;
}

double GridLocalPolynomial::evalBasisRaw( const int point[], const double x[] ) const{
        double f = rule->evalRaw( point[0], x[0] );
        for( int j=1; j<num_dimensions; j++ ) f *= rule->evalRaw( point[j], x[j] );
        return f;
}
double GridLocalPolynomial::evalBasisSupported( const int point[], const double x[], bool &isSupported ) const{
        double f = rule->evalSupport( point[0], x[0], isSupported );
        if ( !isSupported ) return 0.0;
        for( int j=1; j<num_dimensions; j++ ){
                f *= rule->evalSupport( point[j], x[j], isSupported );
                if ( !isSupported ) return 0.0;
        }
        return f;
}

void GridLocalPolynomial::buildTree(){
        if ( roots != 0 ) delete[] roots;
        if ( pntr != 0 ) delete[] pntr;
        if ( indx != 0 ) delete[] indx;

        IndexSet *work = ( points == 0 ) ? needed : points;
        int num_points = work->getNumIndexes();

        int *level = new int[num_points];
        #pragma omp parallel for schedule(static)
        for( int i=0; i<num_points; i++ ){
                const int *p = work->getIndex( i );
                level[i] = rule->getLevel( p[0] );
                for( int j=1; j<num_dimensions; j++ ){
                        level[i] += rule->getLevel( p[j] );
                }
        }

        top_level = 0;
        for( int i=0; i<num_points; i++ ) if ( top_level < level[i] ) top_level = level[i];

        int max_kids = 2*num_dimensions;
        int* monkey_count = new int[top_level + 1];
        int* monkey_tail  = new int[top_level + 1];

        int  *tree = new int[max_kids * num_points];  std::fill( tree, tree + max_kids * num_points, -1 );
        bool *free = new bool[num_points];            std::fill( free, free + num_points, true );
        int  *kid  = new int[num_dimensions];

        int next_root = 0;

        int crts = 12;
        int *rts = new int[crts];
        int num_rts = 0;

        while( next_root != -1 ){
                if ( num_rts == crts ){
                        int *tmp = rts;
                        crts *= 2;
                        rts = new int[crts];
                        std::copy( tmp, tmp + num_rts, rts );
                        delete[] tmp;
                }
                rts[num_rts++] = next_root;
                free[next_root] = false;

                monkey_tail[0] = next_root;
                monkey_count[0] = 0;
                int current = 0;

                while( monkey_count[0] < max_kids ){
                        if ( monkey_count[current] < max_kids ){
                                const int *p = work->getIndex( monkey_tail[current] );

                                int dir = monkey_count[current] / 2;
                                int ikid = ( monkey_count[current] % 2 == 0 ) ? rule->getKidLeft( p[dir] ) : rule->getKidRight( p[dir] );

                                if ( ikid == -1 ){
                                        monkey_count[current]++;
                                }else{
                                        std::copy( p, p + num_dimensions, kid );
                                        kid[dir] = ikid;
                                        int t = work->getSlot( kid );
                                        if ( (t == -1) || (!free[t]) ){
                                                monkey_count[current]++;
                                        }else{
                                                tree[ monkey_tail[current] * max_kids + monkey_count[current] ] =  t;
                                                monkey_count[++current] = 0;
                                                monkey_tail[current] = t;
                                                free[t] = false;
                                        }
                                }
                        }else{
                                monkey_count[--current]++;
                        }
                }

                next_root = -1;
                int next_level = top_level+1;
                for( int i=0; i<num_points; i++ ){
                        if ( free[i] && (level[i] < next_level) ){
                                next_root = i;
                                next_level = level[i];
                        }
                }
        }

        num_roots = num_rts;
        roots = new int[num_roots];
        std::copy( rts, rts + num_roots, roots );

        pntr = new int[num_points + 1];
        pntr[0] = 0;
        for( int i=0; i<num_points; i++ ){
                pntr[i+1] = pntr[i];
                for( int j=0; j<max_kids; j++ ) if ( tree[ i * max_kids + j ] > -1 ) pntr[i+1]++;
        }

        indx = ( pntr[num_points] > 0 ) ? new int[pntr[num_points]] : new int[1]; indx[0] = 0;
        int count = 0;
        for( int i=0; i<num_points; i++ ){
                for( int j=0; j<max_kids; j++ ){
                        int t = tree[ i * max_kids + j ];
                        if ( t > -1 ) indx[count++] = t;
                }
        }

        delete[] kid;
        delete[] free;
        delete[] tree;
        delete[] monkey_tail;
        delete[] monkey_count;
        delete[] rts;
        delete[] level;
}

double* GridLocalPolynomial::getBasisIntegrals() const{
        IndexSet *work = ( points == 0 ) ? needed : points;

        int n = 0;
        double *w = 0, *x = 0;

        if ( (rule->getMaxOrder() == -1) || (rule->getMaxOrder() > 3)  ){
                OneDimensionalNodes GL;
                n = top_level / 2 + 1;
                GL.getGaussLegendre( n, w, x );
        }

        double *integrals = new double[work->getNumIndexes()];
        for( int i=0; i<work->getNumIndexes(); i++ ){
                const int* p = work->getIndex( i );
                integrals[i] = rule->getArea( p[0], n, w, x );
                for( int j=1; j<num_dimensions; j++ ){
                        integrals[i] *= rule->getArea( p[j], n, w, x );
                }
        }

        if ( (rule->getMaxOrder() == -1) || (rule->getMaxOrder() > 3) ){
                delete[] x;
                delete[] w;
        }

        return integrals;
}

double* GridLocalPolynomial::getQuadratureWeights() const{
        IndexSet *work = ( points == 0 ) ? needed : points;
        double *weights = getBasisIntegrals();

        int *monkey_count = new int[top_level+1];
        int *monkey_tail = new int[top_level+1];

        double basis_value;

        const int *dagUp;
        if ( num_outputs == 0 ){
                dagUp = parents;
        }else{
                IndexManipulator IM(num_dimensions);
                dagUp = IM.computeDAGupLocal( work, rule );
        }

        int num_points = work->getNumIndexes();
        int *level = new int[num_points];
        for( int i=0; i<num_points; i++ ){
                const int *p = work->getIndex( i );
                level[i] = rule->getLevel( p[0] );
                for( int j=1; j<num_dimensions; j++ ){
                        level[i] += rule->getLevel( p[j] );
                }
        }

        bool *used = new bool[work->getNumIndexes()];
        double *node = new double[num_dimensions];
        int max_parents = ( rule->isSemiLocal() ) ? 2*num_dimensions : num_dimensions;

        for( int l=top_level; l>0; l-- ){
                for( int i=0; i<num_points; i++ ){
                        if ( level[i] == l ){
                                const int* p = work->getIndex( i );
                                for( int j=0; j<num_dimensions; j++ ) node[j] = rule->getNode( p[j] );

                                std::fill( used, used + work->getNumIndexes(), false );

                                monkey_count[0] = 0;
                                monkey_tail[0] = i;
                                int current = 0;

                                while( monkey_count[0] < max_parents ){
                                        if ( monkey_count[current] < max_parents ){
                                                int branch = dagUp[  monkey_tail[current] * max_parents + monkey_count[current]  ];
                                                if ( (branch == -1) || used[branch] ){
                                                        monkey_count[current]++;
                                                }else{
                                                        const int *func = work->getIndex( branch );
                                                        basis_value = rule->evalRaw( func[0], node[0] );
                                                        for( int j=1; j<num_dimensions; j++ ) basis_value *= rule->evalRaw( func[j], node[j] );
                                                        weights[branch] -= weights[i] * basis_value;
                                                        used[branch] = true;

                                                        monkey_count[++current] = 0;
                                                        monkey_tail[current] = branch;
                                                }
                                        }else{
                                                monkey_count[--current]++;
                                        }
                                }
                        }
                }
        }


        delete[] used;
        delete[] node;

        delete[] level;

        delete[] monkey_count;
        delete[] monkey_tail;

        if ( num_outputs > 0 )  delete[] dagUp;

        return weights;
}

void GridLocalPolynomial::integrate( double q[] ) const{
        double *integrals = getBasisIntegrals();

        std::fill( q, q + num_outputs, 0.0 );

        int n = points->getNumIndexes();

        for( int i=0; i<n; i++ ){
                for( int k=0; k<num_outputs; k++ ){
                        q[k] += integrals[i] * surpluses[ i*num_outputs + k ];
                }
        }
        delete[] integrals;
}

double* GridLocalPolynomial::getNormalization() const{
        double* norm = new double[num_outputs];  std::fill( norm, norm + num_outputs, 0.0 );
        for( int i=0; i<points->getNumIndexes(); i++ ){
                const double *v = values->getValues( i );
                for( int j=0; j<num_outputs; j++ ){
                        if ( norm[j] < fabs( v[j] ) ) norm[j] = fabs( v[j] );
                }
        }
        return norm;
}

int* GridLocalPolynomial::buildUpdateMap( double tolerance, TypeRefinement criteria ) const{
        int num_points = points->getNumIndexes();
        int *map = new int[ num_points * num_dimensions ];  std::fill( map, map + num_points * num_dimensions, 0 );

        double *norm = getNormalization();

        if ( (criteria == refine_classic) || (criteria == refine_parents_first) ){
                #pragma omp parallel for
                for( int i=0; i<num_points; i++ ){
                        bool small = true;
                        for( int k=0; k<num_outputs; k++ ){
                                if ( small && ( (fabs( surpluses[i*num_outputs + k] ) / norm[k]) > tolerance ) ) small = false;
                        }
                        if ( !small ){
                                std::fill( &(map[i*num_dimensions]), &(map[i*num_dimensions]) + num_dimensions, 1 );
                        }
                }
        }else{
                IndexManipulator IM(num_dimensions);
                int *dagUp = IM.computeDAGupLocal( points, rule );

                int max_parents = ( rule->isSemiLocal() ) ? 2*num_dimensions : num_dimensions;
                int max_1D_parents = ( rule->isSemiLocal() ) ? 2 : 1;

                SplitDirections split( points );

                #pragma omp parallel for
                for( int s=0; s<split.getNumJobs(); s++ ){
                        int d = split.getJobDirection( s );
                        int nump = split.getJobNumPoints( s );
                        const int *pnts = split.getJobPoints( s );

                        int *global_to_pnts = new int[num_points];
                        int *levels = new int[nump];
                        int max_level = 0;

                        double *vals = new double[ nump * num_outputs ];

                        for( int i=0; i<nump; i++ ){
                                const double* v = values->getValues( pnts[i] );
                                const int *p = points->getIndex( pnts[i] );
                                std::copy( v, v + num_outputs, &( vals[i*num_outputs] ) );
                                global_to_pnts[ pnts[i] ] = i;
                                levels[i] = rule->getLevel( p[d] );
                                if ( max_level < levels[i] ) max_level = levels[i];
                        }

                        int *monkey_count = new int[max_level + 1 ];
                        int *monkey_tail = new int[max_level + 1 ];
                        bool *used = new bool[nump];

                        for( int l=1; l<=max_level; l++ ){
                                for( int i=0; i<nump; i++ ){
                                        if ( levels[i] == l ){
                                                const int *p = points->getIndex( pnts[i] );
                                                double x = rule->getNode( p[d] );

                                                int current = 0;
                                                monkey_count[0] = d * max_1D_parents;
                                                monkey_tail[0] = pnts[i]; // uses the global indexes
                                                std::fill( used, used + nump, false );

                                                while( monkey_count[0] < (d+1) * max_1D_parents ){
                                                        if ( monkey_count[current] < (d+1) * max_1D_parents ){
                                                                int branch = dagUp[ monkey_tail[current] * max_parents + monkey_count[current] ];
                                                                if ( (branch == -1) || (used[global_to_pnts[branch]]) ){
                                                                        monkey_count[current]++;
                                                                }else{
                                                                        const int *branch_point = points->getIndex( branch );
                                                                        double basis_value = rule->evalRaw( branch_point[d], x );
                                                                        for( int k=0; k<num_outputs; k++ ){
                                                                                vals[ i * num_outputs + k ] -= basis_value * vals[ global_to_pnts[branch] * num_outputs + k ];
                                                                        }

                                                                        used[global_to_pnts[branch]] = true;
                                                                        monkey_count[++current] = d * max_1D_parents;
                                                                        monkey_tail[current] = branch;
                                                                }
                                                        }else{
                                                                monkey_count[--current]++;
                                                        }
                                                }
                                        }
                                }
                        }

                        delete[] used;
                        delete[] monkey_tail;
                        delete[] monkey_count;

                        // at this point, vals contains the one directional surpluses
                        for( int i=0; i<nump; i++ ){
                                bool small = true;
                                for( int k=0; k<num_outputs; k++ ){
                                        if ( small && ( (fabs( surpluses[pnts[i]*num_outputs + k] ) / norm[k]) > tolerance ) && ( (fabs( vals[i*num_outputs + k] ) / norm[k] ) > tolerance ) ) small = false;
                                }
                                map[ pnts[i]*num_dimensions + d ] = ( small ) ? 0 : 1;;
                        }

                        delete[] vals;
                        delete[] global_to_pnts;
                        delete[] levels;
                }

                delete[] dagUp;
        }

        delete[] norm;

        return map;
}

bool GridLocalPolynomial::addParent( const int point[], int direction, GranulatedIndexSet *destination, IndexSet *exclude ) const{
        int *dad = new int[num_dimensions];  std::copy( point, point + num_dimensions, dad );
        bool added = false;
        dad[direction] = rule->getParent( point[direction] );
        if ( (dad[direction] != -1) && (exclude->getSlot( dad ) == -1) ){
                destination->addIndex( dad );
                added = true;
        }
        dad[direction] = rule->getStepParent( point[direction] );
        if ( (dad[direction] != -1) && (exclude->getSlot( dad ) == -1) ){
                destination->addIndex( dad );
                added = true;
        }
        delete[] dad;
        return added;
}
void GridLocalPolynomial::addChild( const int point[], int direction, GranulatedIndexSet *destination, IndexSet *exclude )const{
        int *kid = new int[num_dimensions];  std::copy( point, point + num_dimensions, kid );
        kid[direction] = rule->getKidLeft( point[direction] );
        if ( (kid[direction] != -1) && (exclude->getSlot(kid) == -1) ){
                destination->addIndex( kid );
        }
        kid[direction] = rule->getKidRight( point[direction] );
        if ( (kid[direction] != -1) && (exclude->getSlot(kid) == -1) ){
                destination->addIndex( kid );
        }
        delete[] kid;
}

void GridLocalPolynomial::clearRefinement(){
        if ( needed != 0 ){  delete needed;  needed = 0;  }
}
void GridLocalPolynomial::setSurplusRefinement( double tolerance, TypeRefinement criteria ){
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
