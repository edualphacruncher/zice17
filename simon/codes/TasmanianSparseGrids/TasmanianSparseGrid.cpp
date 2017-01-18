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

#ifndef __TASMANIAN_SPARSE_GRID_CPP
#define __TASMANIAN_SPARSE_GRID_CPP

#include "TasmanianSparseGrid.hpp"

namespace TasGrid{

const char* TasmanianSparseGrid::getVersion(){
        return "3.0";
}
const char* TasmanianSparseGrid::getLicense(){
        return "GPLv3";
}

TasmanianSparseGrid::TasmanianSparseGrid() : base(0), global(0), sequence(0), pwpoly(0), wavelet(0), domain_transform_a(0), domain_transform_b(0){}
TasmanianSparseGrid::~TasmanianSparseGrid(){
        clear();
}

void TasmanianSparseGrid::clear(){
        if ( global != 0 ){ delete global; global = 0; }
        if ( sequence != 0 ){ delete sequence; sequence = 0; }
        if ( pwpoly != 0 ){ delete pwpoly; pwpoly = 0; }
        if ( wavelet != 0 ){ delete wavelet; wavelet = 0; }
        if ( domain_transform_a != 0 ){ delete[] domain_transform_a; domain_transform_a = 0; }
        if ( domain_transform_b != 0 ){ delete[] domain_transform_b; domain_transform_b = 0; }
        base = 0;
}

void TasmanianSparseGrid::write( const char *filename ) const{
        std::ofstream ofs; ofs.open( filename );
        write( ofs );
        ofs.close();
}
bool TasmanianSparseGrid::read( const char *filename ){
        std::ifstream ifs; ifs.open( filename );
        bool isGood = read( ifs );
        ifs.close();
        return isGood;
}

void TasmanianSparseGrid::write( std::ofstream &ofs ) const{
        ofs << "TASMANIAN SG " << getVersion() << endl;
        ofs << "WARNING: do not edit this manually" << endl;
        if ( global != 0 ){
                ofs << "global" << endl;
                global->write( ofs );
        }else if ( sequence != 0 ){
                ofs << "sequence" << endl;
                sequence->write( ofs );
        }else if ( pwpoly != 0 ){
                ofs << "localpolynomial" << endl;
                pwpoly->write( ofs );
        }else if ( wavelet != 0 ){
                ofs << "wavelet" << endl;
                wavelet->write( ofs );
        }else{
                ofs << "empty" << endl;
        }
        ofs << "TASMANIAN SG end" << endl;
}
bool TasmanianSparseGrid::read( std::ifstream &ifs ){
        std::string T;
        ifs >> T;  if ( !(T.compare("TASMANIAN") == 0) ){  cerr << "ERROR: wrong file format, first word in not 'TASMANIAN'" << endl; return false; }
        ifs >> T;  if ( !(T.compare("SG") == 0) ){  cerr << "ERROR: wrong file format, second word in not 'SG'" << endl; return false; }
        getline( ifs, T ); T.erase(0,1); if ( !(T.compare( getVersion() ) == 0) ){  cerr << "WARNING: Version mismatch, possibly undefined behavior!" << endl; }
        getline( ifs, T ); if ( !(T.compare("WARNING: do not edit this manually") == 0) ){ cerr << "ERROR: wrong file format, did not match 'WARNING: do not edit this manually'" << endl; return false; }
        ifs >> T;
        if ( T.compare( "global" ) == 0 ){
                clear();
                global = new GridGlobal();
                global->read( ifs );
                base = global;
        }else if ( T.compare( "sequence" ) == 0 ){
                clear();
                sequence = new GridSequence();
                sequence->read( ifs );
                base = sequence;
        }else if ( T.compare( "localpolynomial" ) == 0 ){
                clear();
                pwpoly = new GridLocalPolynomial();
                pwpoly->read( ifs );
                base = pwpoly;
        }else if ( T.compare( "wavelet" ) == 0 ){
                clear();
                wavelet = new GridWavelet();
                wavelet->read( ifs );
                base = wavelet;
        }else if ( T.compare( "empty" ) == 0 ){
                clear();
        }else{
                cerr << "ERROR: wrong file format!" << endl; return false;
        }
        getline( ifs, T ); getline( ifs, T );
        if ( !(T.compare("TASMANIAN SG end") == 0) ){ cerr << "WARNING: stream did not end with 'TASMANIAN SG end', this may result in undefined behavior" << endl; }
        return true;
}

void TasmanianSparseGrid::makeGlobalGrid( int dimensions, int outputs, int depth, TypeDepth type, TypeOneDRule rule, const int *anisotropic_weights, double alpha, double beta, const char* custom_filename ){
        clear();
        global = new GridGlobal();
        global->makeGrid( dimensions, outputs, depth, type, rule, anisotropic_weights, alpha, beta, custom_filename );
        base = global;
}
void TasmanianSparseGrid::makeSequenceGrid( int dimensions, int outputs, int depth, TypeDepth type, TypeOneDRule rule, const int *anisotropic_weights ){
        if ( outputs < 1 ){
                cerr << "ERROR: makeSequenceGrid is called with zero outputs, for zero outputs use makeGlobalGrid instead" << endl;
                return;
        }
        if ( OneDimensionalMeta::isSequence( rule ) ){
                clear();
                sequence = new GridSequence();
                sequence->makeGrid( dimensions, outputs, depth, type, rule, anisotropic_weights );
                base = sequence;
        }else{
                cerr << "ERROR: makeSequenceGrid is called with rule " << OneDimensionalMeta::getIORuleString(rule) << " which is not a sequence rule" << endl;
        }
}
void TasmanianSparseGrid::makeLocalPolynomialGrid( int dimensions, int outputs, int depth, int order, TypeOneDRule rule ){
        if ( (rule != rule_localp) && (rule != rule_localp0) && (rule != rule_semilocalp) ){
                cout << "ERROR: makeLocalPolynomialGrid is called with rule " << OneDimensionalMeta::getIORuleString(rule) << " which is not a local polynomial rule" << endl;
                cout << "       use either " << OneDimensionalMeta::getIORuleString(rule_localp) << " or " << OneDimensionalMeta::getIORuleString(rule_semilocalp)
                                 << " or " << OneDimensionalMeta::getIORuleString(rule_localp0)  << endl;
                return;
        }
        if ( order < -1 ){
                cout << "ERROR: makeLocalPolynomialGrid is called with order " << order << ", but the order cannot be less than -1." << endl;
                return;
        }
        clear();
        pwpoly = new GridLocalPolynomial();
        pwpoly->makeGrid( dimensions, outputs, depth, order, rule );
        base = pwpoly;
}
void TasmanianSparseGrid::makeWaveletGrid( int dimensions, int outputs, int depth, int order ){
        if ( (order != 1) && (order != 3) ){
                cout << "ERROR: makeWaveletGrid is called with order " << order << ", but wavelets are implemented only for orders 1 and 3." << endl;
                return;
        }
        clear();
        wavelet = new GridWavelet();
        wavelet->makeGrid( dimensions, outputs, depth, order );
        base = wavelet;
}

void TasmanianSparseGrid::updateGlobalGrid( int depth, TypeDepth type, const int *anisotropic_weights ){
        if ( global != 0 ){
                global->updateGrid( depth, type, anisotropic_weights );
        }else{
                cerr << "ERROR: updateGlobalGrid called, but the grid is not global" << endl;
        }
}
void TasmanianSparseGrid::updateSequenceGrid( int depth, TypeDepth type, const int *anisotropic_weights ){
        if ( sequence != 0 ){
                sequence->updateGrid( depth, type, anisotropic_weights );
        }else{
                cerr << "ERROR: updateSequenceGrid called, but the grid is not sequence" << endl;
        }
}

double TasmanianSparseGrid::getAlpha() const{
        return ( global != 0 ) ? global->getAlpha() : 0.0;
}
double TasmanianSparseGrid::getBeta() const{
        return ( global != 0 ) ? global->getBeta() : 0.0;
}
int TasmanianSparseGrid::getOrder() const{
        return ( pwpoly != 0 ) ? pwpoly->getOrder() : ( ( wavelet != 0 ) ? wavelet->getOrder() : -1 );
}

int TasmanianSparseGrid::getNumDimensions() const{ return base->getNumDimensions(); }
int TasmanianSparseGrid::getNumOutputs() const{ return base->getNumOutputs(); }
TypeOneDRule TasmanianSparseGrid::getRule() const{ return base->getRule(); }
const char* TasmanianSparseGrid::getCustomRuleDescription() const{ return (global != 0) ? global->getCustomRuleDescription() : ""; }

int TasmanianSparseGrid::getNumLoaded() const{ return base->getNumLoaded(); }
int TasmanianSparseGrid::getNumNeeded() const{ return base->getNumNeeded(); }
int TasmanianSparseGrid::getNumPoints() const{ return base->getNumPoints(); }

double* TasmanianSparseGrid::getLoadedPoints() const{
        if ( domain_transform_a == 0 ){
                return base->getLoadedPoints();
        }else{
                double *x = base->getLoadedPoints();
                mapCanonicalToTransformed( base->getNumDimensions(), base->getNumLoaded(), base->getRule(), x );
                return x;
        }
}
double* TasmanianSparseGrid::getNeededPoints() const{
        if ( domain_transform_a == 0 ){
                return base->getNeededPoints();
        }else{
                double *x = base->getNeededPoints();
                mapCanonicalToTransformed( base->getNumDimensions(), base->getNumNeeded(), base->getRule(), x );
                return x;
        }
}
double* TasmanianSparseGrid::getPoints() const{
        if ( domain_transform_a == 0 ){
                return base->getPoints();
        }else{
                double *x = base->getPoints();
                mapCanonicalToTransformed( base->getNumDimensions(), base->getNumPoints(), base->getRule(), x );
                return x;
        }
}

double* TasmanianSparseGrid::getQuadratureWeights() const{
        if ( domain_transform_a == 0 ){
                return base->getQuadratureWeights();
        }else{
                double *w = base->getQuadratureWeights();
                double scale = getQuadratureScale( base->getNumDimensions(), base->getRule() );
                #pragma omp parallel for schedule(static)
                for( int i=0; i<getNumPoints(); i++ ) w[i] *= scale;
                return w;
        }
}
double* TasmanianSparseGrid::getInterpolationWeights( const double x[] ) const{
        if ( domain_transform_a == 0 ){
                return base->getInterpolationWeights( x );
        }else{
                int num_dimensions = base->getNumDimensions();
                double *x_canonical = new double[num_dimensions];  std::copy( x, x + num_dimensions, x_canonical );
                mapTransformedToCanonical( num_dimensions, base->getRule(), x_canonical );
                double *w = base->getInterpolationWeights( x_canonical );
                delete[] x_canonical;
                return w;
        }
}

void TasmanianSparseGrid::loadNeededPoints( const double *vals ){ base->loadNeededPoints( vals ); }

void TasmanianSparseGrid::evaluate( const double x[], double y[] ) const{
        if ( domain_transform_a == 0 ){
                base->evaluate( x, y );
        }else{
                int num_dimensions = base->getNumDimensions();
                double *x_canonical = new double[num_dimensions];  std::copy( x, x + num_dimensions, x_canonical );
                mapTransformedToCanonical( num_dimensions, base->getRule(), x_canonical );
                base->evaluate( x_canonical, y );
                delete[] x_canonical;
        }
}
void TasmanianSparseGrid::integrate( double q[] ) const{
        base->integrate( q );
        if ( domain_transform_a != 0 ){
                double scale = getQuadratureScale( base->getNumDimensions(), base->getRule() );
                for( int k=0; k<getNumOutputs(); k++ ) q[k] *= scale;
        }
}

bool TasmanianSparseGrid::isGlobal() const{
        return (global != 0);
}
bool TasmanianSparseGrid::isSequence() const{
        return (sequence != 0);
}
bool TasmanianSparseGrid::isLocalPolynomial() const{
        return (pwpoly != 0);
}
bool TasmanianSparseGrid::isWavelet() const{
        return (wavelet != 0);
}

void TasmanianSparseGrid::setDomainTransform( const double a[], const double b[] ){
        if ( (base == 0) || (base->getNumDimensions() == 0) ){
                cerr << "ERROR: cannot call setDomainTransform on uninitialized grid!" << endl;
                return;
        }
        clearDomainTransform();
        int num_dimensions = base->getNumDimensions();
        domain_transform_a = new double[num_dimensions];  std::copy( a, a + num_dimensions, domain_transform_a );
        domain_transform_b = new double[num_dimensions];  std::copy( b, b + num_dimensions, domain_transform_b );
}
bool TasmanianSparseGrid::isSetDomainTransfrom() const{
        return ( domain_transform_a != 0 );
}
void TasmanianSparseGrid::clearDomainTransform(){
        if ( domain_transform_a != 0 ){ delete[] domain_transform_a; domain_transform_a = 0; }
        if ( domain_transform_b != 0 ){ delete[] domain_transform_b; domain_transform_b = 0; }
}
void TasmanianSparseGrid::getDomainTransform( double a[], double b[] ) const{
        if ( (base == 0) || (base->getNumDimensions() == 0) || (domain_transform_a == 0) ){
                cerr << "ERROR: cannot call getDomainTransform on uninitialized grid or if no transform has been set!" << endl;
                return;
        }
        int num_dimensions = base->getNumDimensions();
        std::copy( domain_transform_a, domain_transform_a + num_dimensions, a );
        std::copy( domain_transform_b, domain_transform_b + num_dimensions, b );
}

void TasmanianSparseGrid::mapCanonicalToTransformed( int num_dimensions, int num_points, TypeOneDRule rule, double x[] ) const{
        if ( rule == rule_gausslaguerre ){
                for( int i=0; i<num_points; i++ ){
                        for( int j=0; j<num_dimensions; j++ ){
                                x[i*num_dimensions+j] /= domain_transform_b[j];
                                x[i*num_dimensions+j] += domain_transform_a[j];
                        }
                }
        }else if ( (rule == rule_gausshermite) || (rule == rule_gausshermiteodd) ){
                for( int i=0; i<num_points; i++ ){
                        for( int j=0; j<num_dimensions; j++ ){
                                x[i*num_dimensions+j] /= sqrt( domain_transform_b[j] );
                                x[i*num_dimensions+j] += domain_transform_a[j];
                        }
                }
        }else{
                for( int i=0; i<num_points; i++ ){
                        for( int j=0; j<num_dimensions; j++ ){
                                x[i*num_dimensions+j] *= 0.5* ( domain_transform_b[j] - domain_transform_a[j] );
                                x[i*num_dimensions+j] += 0.5* ( domain_transform_b[j] + domain_transform_a[j] );
                        }
                }
        }
}
void TasmanianSparseGrid::mapTransformedToCanonical( int num_dimensions, TypeOneDRule rule, double x[] ) const{
        if ( rule == rule_gausslaguerre ){
                for( int j=0; j<num_dimensions; j++ ){
                        x[j] -= domain_transform_a[j];
                        x[j] *= domain_transform_b[j];
                }
        }else if ( (rule == rule_gausshermite) || (rule == rule_gausshermiteodd) ){
                for( int j=0; j<num_dimensions; j++ ){
                        x[j] -= domain_transform_a[j];
                        x[j] *= sqrt( domain_transform_b[j] );
                }
        }else{
                for( int j=0; j<num_dimensions; j++ ){
                        x[j] *= 2.0;
                        x[j] -= ( domain_transform_b[j] + domain_transform_a[j] );
                        x[j] /= ( domain_transform_b[j] - domain_transform_a[j] );
                }
        }
}
double TasmanianSparseGrid::getQuadratureScale( int num_dimensions, TypeOneDRule rule ) const{
        double scale = 1.0;
        if ( (rule == rule_gausschebyshev1) || (rule == rule_gausschebyshev2) || (rule == rule_gaussgegenbauer) || (rule == rule_gaussjacobi) ){
                double alpha = ( rule == rule_gausschebyshev1 ) ? -0.5 : ( rule == rule_gausschebyshev2 ) ? 0.5 : global->getAlpha();
                double beta = ( rule == rule_gausschebyshev1 ) ? -0.5 : ( rule == rule_gausschebyshev2 ) ? 0.5 : ( rule == rule_gaussgegenbauer ) ? global->getAlpha() : global->getBeta();
                for( int j=0; j<num_dimensions; j++ ) scale *= pow( 0.5*( domain_transform_b[j] - domain_transform_a[j] ), alpha + beta + 1.0 );
        }else if ( rule == rule_gausslaguerre ){
                for( int j=0; j<num_dimensions; j++ ) scale *= pow( domain_transform_b[j], global->getAlpha() + 1.0 );
        }else if ( (rule == rule_gausshermite) || (rule == rule_gausshermiteodd) ){
                double power = -0.5 * ( 1.0 + global->getAlpha() );
                for( int j=0; j<num_dimensions; j++ ) scale *= pow( domain_transform_b[j], power );
        }else{
                for( int j=0; j<num_dimensions; j++ ) scale *= (domain_transform_b[j] - domain_transform_a[j] ) / 2.0;
        }
        return scale;
}

void TasmanianSparseGrid::setAnisotropicRefinement( TypeDepth type, int min_growth, int output ){
        if ( sequence != 0 ){
                sequence->setAnisotropicRefinement( type, min_growth );
        }else if ( global != 0 ){
                if ( OneDimensionalMeta::isNonNested( global->getRule() ) ){
                        cerr << "ERROR: setAnisotropicRefinement called for a global grid with non-nested rule" << endl;
                }else{
                        global->setAnisotropicRefinement( type, min_growth, output );
                }
        }else{
                cerr << "ERROR: setAnisotropicRefinement called for grid that is neither sequence nor Global with sequence rule" << endl;
        }
}
int* TasmanianSparseGrid::estimateAnisotropicCoefficients( TypeDepth type, int output ){
        if ( sequence != 0 ){
                return sequence->estimateAnisotropicCoefficients( type );
        }else if ( global != 0 ){
                if ( OneDimensionalMeta::isNonNested( global->getRule() ) ){
                        cerr << "ERROR: estimateAnisotropicCoefficients called for a global grid with non-nested rule" << endl;
                }else{
                        return global->estimateAnisotropicCoefficients( type, output );
                }
        }else{
                cerr << "ERROR: estimateAnisotropicCoefficients called for grid that is neither sequence nor Global with sequence rule" << endl;
        }
        return 0;
}
void TasmanianSparseGrid::setSurplusRefinement( double tolerance, int output ){
        if ( sequence != 0 ){
                sequence->setSurplusRefinement( tolerance );
        }else if ( global != 0 ){
                if ( OneDimensionalMeta::isSequence( global->getRule() ) ){
                        global->setSurplusRefinement( tolerance, output );
                }else{
                        cerr << "ERROR: setSurplusRefinement called for a global grid with non-sequence rule" << endl;
                }
        }else{
                cerr << "ERROR: setSurplusRefinement( double, int ) called for grid that is neither sequence nor Global with sequence rule" << endl;
        }
}
void TasmanianSparseGrid::setSurplusRefinement( double tolerance, TypeRefinement criteria ){
        if ( pwpoly != 0 ){
                pwpoly->setSurplusRefinement( tolerance, criteria );
        }else if ( wavelet != 0 ){
                wavelet->setSurplusRefinement( tolerance, criteria );
        }else{
                cerr << "ERROR: setSurplusRefinement( double, TypeRefinement ) called for grid that is neither local polynomial nor wavelet" << endl;
        }
}
void TasmanianSparseGrid::clearRefinement(){
        base->clearRefinement();
}

int* TasmanianSparseGrid::getGlobalPolynomialSpace( bool interpolation, int &num_indexes ) const{
        if ( global != 0 ){
                return global->getPolynomialSpace( interpolation, num_indexes );
        }else if ( sequence != 0 ){
                return sequence->getPolynomialSpace( interpolation, num_indexes );
        }else{
                cout << "ERROR: getGlobalPolynomialSpace() called for a grid that is neither Global nor Sequence." << endl;
                num_indexes = 0;
                return 0;
        }
}

void TasmanianSparseGrid::printStats() const{
        using std::setw;

        const int L1 = 20, L2 = 10;
        cout << endl;
        cout << setw(L1) << "Grid Type:" << setw(L2) << " " << " ";
        if ( isGlobal() ) cout << "Global" << endl;
        if ( isSequence() ) cout << "Sequence" << endl;
        if ( isLocalPolynomial() ) cout << "Local Polynomial" << endl;
        if ( isWavelet() ) cout << "Wavelets" << endl;

        cout << setw(L1) << "Dimensions:" << setw(L2) << " " << getNumDimensions() << endl;
        cout << setw(L1) << "Outputs:" << setw(L2) << " " << getNumOutputs() << endl;
        if ( getNumOutputs() == 0 ){
                cout << setw(L1) << "Nodes:" << setw(L2) << " " << getNumPoints() << endl;
        }else{
                cout << setw(L1) << "Loaded nodes:" << setw(L2) << " " << getNumLoaded() << endl;
                cout << setw(L1) << "Needed nodes:" << setw(L2) << " " << getNumNeeded() << endl;
        }
        cout << setw(L1) << "Rule:" << setw(L2) << " " << " " << OneDimensionalMeta::getHumanString( getRule() ) << endl;
        if ( getRule() == rule_customtabulated ){
                cout << setw(L1) << "Description:" << setw(L2) << " " << " " << getCustomRuleDescription() << endl;
        }
        if ( isSetDomainTransfrom() ){
                cout << setw(L1) << "Domain:" << setw(L2) << "Custom" << endl;
        }else{
                cout << setw(L1) << "Domain:" << setw(L2) << "Canonical" << endl;
        }

        if ( isGlobal() ){
                TypeOneDRule rr = getRule();
                if ( (rr == rule_gaussgegenbauer) || (rr == rule_gausslaguerre) || (rr == rule_gausshermite) || (rr == rule_gaussgegenbauerodd) || (rr == rule_gausshermiteodd)  ){
                        cout << setw(L1) << "Alpha:" << setw(L2) << " " << getAlpha() << endl;
                }
                if ( rr == rule_gaussjacobi ){
                        cout << setw(L1) << "Alpha:" << setw(L2) << " " << getAlpha() << endl;
                        cout << setw(L1) << "Beta:" << setw(L2) << " " << getBeta() << endl;
                }
        }else if ( isSequence() ){
                // sequence rules are simple, nothing to specify here
        }else if ( isLocalPolynomial() ){
                int o = getOrder();
                cout << setw(L1) << "Order:" << setw(L2) << " " << o << endl;
        }else if ( isWavelet() ){
                int o = getOrder();
                cout << setw(L1) << "Order:" << setw(L2) << " " << o << endl;
        }else{
                cerr << endl << "ERROR: unknown grid type!" << endl;
        }

        cout << endl;
}

}

#endif
