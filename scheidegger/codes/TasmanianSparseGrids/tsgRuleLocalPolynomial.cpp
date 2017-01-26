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

#ifndef __TSG_RULE_LOCAL_POLYNOMIAL_CPP
#define __TSG_RULE_LOCAL_POLYNOMIAL_CPP

#include "tsgRuleLocalPolynomial.hpp"

namespace TasGrid{

BaseRuleLocalPolynomial::BaseRuleLocalPolynomial(){}
BaseRuleLocalPolynomial::~BaseRuleLocalPolynomial(){}

int BaseRuleLocalPolynomial::getMaxOrder() const{ return max_order; }
void BaseRuleLocalPolynomial::setMaxOrder( int order ){ max_order = order; }
TypeOneDRule BaseRuleLocalPolynomial::getType() const{ return rule_none; }

int BaseRuleLocalPolynomial::getNumPoints( int level ) const{  return 0;  }
double BaseRuleLocalPolynomial::getNode( int point ) const{ return 0.0; }
int BaseRuleLocalPolynomial::getLevel( int point ) const{  return 0;  }
bool BaseRuleLocalPolynomial::isSemiLocal() const{  return false;  }

int BaseRuleLocalPolynomial::getParent( int point ) const{  return 0;  }
int BaseRuleLocalPolynomial::getStepParent( int point ) const{  return 0;  }
int BaseRuleLocalPolynomial::getKidLeft( int point ) const{  return 0;  }
int BaseRuleLocalPolynomial::getKidRight( int point ) const{  return 0;  }

double BaseRuleLocalPolynomial::evalRaw( int point, double x ) const{  return 0.0;  }
double BaseRuleLocalPolynomial::evalSupport( int point, double x, bool &isSupported ) const{  return 0.0;  }
double BaseRuleLocalPolynomial::getArea( int point, int n, const double w[], const double x[] ) const{  return 0.0;  }

int BaseRuleLocalPolynomial::intlog2( int i ){ // this is effectively: floor( log_2( i ) )
        int result = 0;
        while (i >>= 1){ result++; }
        return result;
}
int BaseRuleLocalPolynomial::int2log2( int i ){ // this is effectively: 2^( floor( log_2( i ) ) ), when i == 0, this returns 1
        int result = 1;
        while (i >>= 1){ result <<= 1; }
        return result;
}

RuleLocalPolynomial::RuleLocalPolynomial(){}
RuleLocalPolynomial::~RuleLocalPolynomial(){}

int RuleLocalPolynomial::getMaxOrder() const{ return max_order; }
void RuleLocalPolynomial::setMaxOrder( int order ){ max_order = order; }
TypeOneDRule RuleLocalPolynomial::getType() const{  return rule_localp; }
int RuleLocalPolynomial::getNumPoints( int level ) const{
        return (level == 0 ) ? 1 : ( (1 << level) + 1 );
}
double RuleLocalPolynomial::getNode( int point ) const{
        if ( point == 0 ) return  0.0;
        if ( point == 1 ) return -1.0;
        if ( point == 2 ) return  1.0;
        return ( (double)( 2*point - 1)  ) / ( (double) int2log2( point - 1 )  ) - 3.0;
}
int RuleLocalPolynomial::getLevel( int point ) const{
        return (point == 0) ? 0 : (point == 1) ? 1 : intlog2(point - 1) + 1;
}
bool RuleLocalPolynomial::isSemiLocal() const{  return false;  }

int RuleLocalPolynomial::getParent( int point ) const{
        int dad = ( point + 1 ) / 2;
        if ( point < 4 ) dad--;
        return dad;
}
int RuleLocalPolynomial::getStepParent( int point ) const{  return -1;  }

int RuleLocalPolynomial::getKidLeft( int point ) const{
        if ( point == 0 ) return 1;
        if ( point == 1 ) return 3;
        if ( point == 2 ) return 4;
        return 2*point-1;
}
int RuleLocalPolynomial::getKidRight( int point ) const{
        if ( point == 0 ) return 2;
        if ( (point == 1) || (point == 2) ) return -1;
        return 2*point;
}
double RuleLocalPolynomial::evalRaw( int point, double x ) const{
        if ( (max_order == 0) || (point == 0) ) return 1.0;
        double xn = scaleX( point, x );
        if ( max_order == 1 ) return 1.0 - fabs( xn );
        if ( max_order == 2 ) return evalPWQuadratic( point, xn );
        if ( max_order == 3 ) return evalPWCubic( point, xn );
        return evalPWPower( point, xn );
}

double RuleLocalPolynomial::evalSupport( int point, double x, bool &isSupported ) const{
        isSupported = true;
        if ( point == 0 ) return 1.0;
        double xn = scaleX( point, x );
        if ( fabs( xn ) <= 1.0 ){
                if ( max_order == 0 ) return 1.0;
                if ( max_order == 1 ) return 1.0 - fabs( xn );
                if ( max_order == 2 ) return evalPWQuadratic( point, xn );
                if ( max_order == 3 ) return evalPWCubic( point, xn );
                return evalPWPower( point, xn );
        }else{
                isSupported = false;
                return 0.0;
        }
}
double RuleLocalPolynomial::getArea( int point, int n, const double w[], const double x[] ) const{
        if ( point == 0 ) return 2.0;
        if ( max_order == 0 ){
                if ( (point == 1) || (point == 2) ) return 1.0;
                return 2.0 * getSupport( point );
        }
        if ( (point == 1) || (point == 2) ) return 0.5;
        if ( max_order == 1 ) return getSupport( point );
        if ( (max_order == 2) || (max_order == 3) || (point <= 8) )  return (4.0/3.0) * getSupport( point );
        double sum = 0.0;
        for( int i=0; i<n; i++ ) sum += w[i] * evalPWPower( point, x[i] );
        return sum * getSupport( point );
}
double RuleLocalPolynomial::getSupport( int point ) const{
        return ( point == 0 ) ? 1.0 : 1.0 / ( (double) int2log2( point - 1) );
}
double RuleLocalPolynomial::scaleX( int point, double x ) const{
        if ( point == 0 ) return x;
        if ( point == 1 ) return ( x + 1.0 );
        if ( point == 2 ) return ( x - 1.0 );
        return ( (double) int2log2( point - 1 ) * ( x + 3.0 ) + 1.0 - (double) ( 2*point ) );
}
double RuleLocalPolynomial::evalPWQuadratic( int point, double x ) const{
        if ( point == 1 ) return 1.0 - x;
        if ( point == 2 ) return 1.0 + x;
        return ( 1.0 - x ) * ( 1.0 + x );
}
double RuleLocalPolynomial::evalPWCubic( int point, double x ) const{
        if ( point == 0 ) return 1.0;
        if ( point == 1 ) return 1.0 - x;
        if ( point == 2 ) return 1.0 + x;
        if ( point <= 4 ) return ( 1.0 - x ) * ( 1.0 + x );
        return ( point % 2 == 0 ) ? ( 1.0 - x ) * ( 1.0 + x ) * ( 3.0 + x ) / 3.0 : ( 1.0 - x ) * ( 1.0 + x ) * ( 3.0 - x ) / 3.0;
}
double RuleLocalPolynomial::evalPWPower( int point, double x ) const{
        if ( point <= 8 ) return evalPWCubic( point, x ); //cout << " x = " << x << endl;
        int level = getLevel( point );
        int z, mod = 1;
        double value = ( 1.0 - x )*( 1.0 + x ), offset = 1.0;
        int imax = ( max_order < 0 ) ? level-2 : ( ( max_order < level ) ? max_order -2 : level-2 );
        for( int j=0; j < imax; j++ ){
                mod *= 2;
                offset = 2.0 * offset + 1.0;
                z = (point-1) % mod;
                if ( z < mod / 2 ){
                        value *= ( x - offset + 2.0 * ( (double) z) ) / ( - offset + 2.0 * ( (double) z) );
                }else{
                        value *= ( x + offset - 2.0 * ((double) ( mod - 1 - z )) ) / ( offset - 2.0 * ((double) ( mod - 1 - z )) );
                }
        }
        return value;
}

RuleSemiLocalPolynomial::RuleSemiLocalPolynomial(){}
RuleSemiLocalPolynomial::~RuleSemiLocalPolynomial(){}

int RuleSemiLocalPolynomial::getMaxOrder() const{ return max_order; }
void RuleSemiLocalPolynomial::setMaxOrder( int order ){ max_order = order; }
TypeOneDRule RuleSemiLocalPolynomial::getType() const{  return rule_semilocalp; }
int RuleSemiLocalPolynomial::getNumPoints( int level ) const{
        return (level == 0 ) ? 1 : ( (1 << level) + 1 );
}
double RuleSemiLocalPolynomial::getNode( int point ) const{
        if ( point == 0 ){ return  0.0; }else
        if ( point == 1 ){ return -1.0; }else
        if ( point == 2 ){ return  1.0; }
        return ( (double)( 2*point - 1)  ) / ( (double) int2log2( point - 1 )  ) - 3.0;
}
int RuleSemiLocalPolynomial::getLevel( int point ) const{
        return (point == 0) ? 0 : (point == 1) ? 1 : intlog2(point - 1) + 1;
}
bool RuleSemiLocalPolynomial::isSemiLocal() const{  return true;  }

int RuleSemiLocalPolynomial::getParent( int point ) const{
        int dad = ( point + 1 ) / 2;
        if ( point < 4 ) dad--;
        return dad;
}
int RuleSemiLocalPolynomial::getStepParent( int point ) const{
        if ( point == 3 ) return 2;
        if ( point == 4 ) return 1;
        return -1;
}

int RuleSemiLocalPolynomial::getKidLeft( int point ) const{
        if ( point == 0 ) return 1;
        if ( point == 1 ) return 3;
        if ( point == 2 ) return 4;
        return 2*point-1;
}
int RuleSemiLocalPolynomial::getKidRight( int point ) const{
        if ( point == 0 ) return 2;
        if ( (point == 1) || (point == 2) ) return -1;
        return 2*point;
}
double RuleSemiLocalPolynomial::evalRaw( int point, double x ) const{
        if ( point == 0 ) return 1.0;
        if ( point == 1 ) return 0.5 * x * ( x - 1.0 );
        if ( point == 2 ) return 0.5 * x * ( x + 1.0 );
        double xn = scaleX( point, x );
        if ( max_order == 2 ) return evalPWQuadratic( point, xn );
        if ( max_order == 3 ) return evalPWCubic( point, xn );
        return evalPWPower( point, xn );
}

double RuleSemiLocalPolynomial::evalSupport( int point, double x, bool &isSupported ) const{
        isSupported = true;
        if ( point == 0 ) return 1.0;
        if ( point == 1 ) return 0.5 * x * ( x - 1.0 );
        if ( point == 2 ) return 0.5 * x * ( x + 1.0 );
        double xn = scaleX( point, x );
        if ( fabs( xn ) <= 1.0 ){
                if ( max_order == 2 ) return evalPWQuadratic( point, xn );
                if ( max_order == 3 ) return evalPWCubic( point, xn );
                return evalPWPower( point, xn );
        }else{
                isSupported = false;
                return 0.0;
        }
}
double RuleSemiLocalPolynomial::getArea( int point, int n, const double w[], const double x[] ) const{
        if ( point == 0 ) return 2.0;
        if ( max_order == 0 ){
                if ( (point == 1) || (point == 2) ) return 1.0;
                return 2.0 * getSupport( point );
        }
        if ( (point == 1) || (point == 2) ) return 1.0/3.0;
        if ( (max_order == 2) || (max_order == 3) || (point <= 4) )  return (4.0/3.0) * getSupport( point );
        double sum = 0.0;
        for( int i=0; i<n; i++ ) sum += w[i] * evalPWPower( point, x[i] );
        return sum * getSupport( point );
}
double RuleSemiLocalPolynomial::getSupport( int point ) const{
        return ( point == 0 ) ? 1.0 : 1.0 / ( (double) int2log2( point - 1) );
}
double RuleSemiLocalPolynomial::scaleX( int point, double x ) const{
        return ( (double) int2log2( point - 1 ) * ( x + 3.0 ) + 1.0 - (double) ( 2*point ) );
}
double RuleSemiLocalPolynomial::evalPWQuadratic( int point, double x ) const{
        return ( 1.0 - x ) * ( 1.0 + x );
}
double RuleSemiLocalPolynomial::evalPWCubic( int point, double x ) const{
        return ( point % 2 == 0 ) ? ( 1.0 - x ) * ( 1.0 + x ) * ( 3.0 + x ) / 3.0 : ( 1.0 - x ) * ( 1.0 + x ) * ( 3.0 - x ) / 3.0;
}
double RuleSemiLocalPolynomial::evalPWPower( int point, double x ) const{
        if ( point <= 4 ) return evalPWCubic( point, x ); //cout << " x = " << x << endl;
        int level = getLevel( point );
        int z, mod = 1;
        double value = ( 1.0 - x )*( 1.0 + x ), offset = 1.0;
        int imax = ( max_order < 0 ) ? level-1 : ( ( max_order-1 < level ) ? max_order -2 : level-1 );
        for( int j=0; j < imax; j++ ){
                mod *= 2;
                offset = 2.0 * offset + 1.0;
                z = (point-1) % mod;
                if ( z < mod / 2 ){
                        value *= ( x - offset + 2.0 * ( (double) z) ) / ( - offset + 2.0 * ( (double) z) );
                }else{
                        value *= ( x + offset - 2.0 * ((double) ( mod - 1 - z )) ) / ( offset - 2.0 * ((double) ( mod - 1 - z )) );
                }
        }
        return value;
}

RuleLocalPolynomialZero::RuleLocalPolynomialZero(){}
RuleLocalPolynomialZero::~RuleLocalPolynomialZero(){}

int RuleLocalPolynomialZero::getMaxOrder() const{ return max_order; }
void RuleLocalPolynomialZero::setMaxOrder( int order ){ max_order = order; }
TypeOneDRule RuleLocalPolynomialZero::getType() const{  return rule_localp0; }
int RuleLocalPolynomialZero::getNumPoints( int level ) const{
        return ( 1 << (level+1) ) -1;
}
double RuleLocalPolynomialZero::getNode( int point ) const{
        if ( point == 0 ) return  0.0;
        return ( (double)( 2*point +3)  ) / ( (double) int2log2( point + 1 )  ) - 3.0;
}
int RuleLocalPolynomialZero::getLevel( int point ) const{
        return intlog2( point + 1 );
}
bool RuleLocalPolynomialZero::isSemiLocal() const{  return false;  }

int RuleLocalPolynomialZero::getParent( int point ) const{
        if ( point == 0 ) return -1;
        return ( point - 1 ) / 2;
}
int RuleLocalPolynomialZero::getStepParent( int point ) const{  return -1;  }

int RuleLocalPolynomialZero::getKidLeft( int point ) const{
        return 2*point+1;
}
int RuleLocalPolynomialZero::getKidRight( int point ) const{
        return 2*point+2;
}
double RuleLocalPolynomialZero::evalRaw( int point, double x ) const{
        if ( max_order == 0 ) return 1.0;
        double xn = scaleX( point, x );
        if ( max_order == 1 ) return 1.0 - fabs( xn );
        if ( max_order == 2 ) return ( 1.0 - xn ) * ( 1.0 + xn );
        if ( max_order == 3 ) return evalPWCubic( point, xn );
        return evalPWPower( point, xn );
}

double RuleLocalPolynomialZero::evalSupport( int point, double x, bool &isSupported ) const{
        isSupported = true;
        double xn = scaleX( point, x );
        if ( fabs( xn ) <= 1.0 ){
                if ( max_order == 0 ) return 1.0;
                if ( max_order == 1 ) return 1.0 - fabs( xn );
                if ( max_order == 2 ) return ( 1.0 - xn ) * ( 1.0 + xn );
                if ( max_order == 3 ) return evalPWCubic( point, xn );
                return evalPWPower( point, xn );
        }else{
                isSupported = false;
                return 0.0;
        }
}
double RuleLocalPolynomialZero::getArea( int point, int n, const double w[], const double x[] ) const{
        if ( max_order == 0 ){
                if ( point == 0 ) return 2.0;
                return 2.0 * getSupport( point );
        }
        if ( max_order == 1 ) return getSupport( point );
        if ( (max_order == 2) || (max_order == 3) || (point <= 2) )  return (4.0/3.0) * getSupport( point );
        double sum = 0.0;
        for( int i=0; i<n; i++ ) sum += w[i] * evalPWPower( point, x[i] );
        return sum * getSupport( point );
}
double RuleLocalPolynomialZero::getSupport( int point ) const{
        return 1.0 / ( (double) int2log2( point + 1) );
}
double RuleLocalPolynomialZero::scaleX( int point, double x ) const{
        if ( point == 0 ) return x;
        return ( (double) int2log2( point + 1 ) * ( x + 3.0 ) - 3.0 - (double) ( 2*point ) );
}
double RuleLocalPolynomialZero::evalPWCubic( int point, double x ) const{
        if ( point == 0 ) return ( 1.0 - x ) * ( 1.0 + x );
        return ( point % 2 == 0 ) ? ( 1.0 - x ) * ( 1.0 + x ) * ( 3.0 + x ) / 3.0 : ( 1.0 - x ) * ( 1.0 + x ) * ( 3.0 - x ) / 3.0;
}
double RuleLocalPolynomialZero::evalPWPower( int point, double x ) const{
        if ( point <= 2 ) return evalPWCubic( point, x );
        int level = getLevel( point );
        int z, mod = 1;
        double value = ( 1.0 - x )*( 1.0 + x ), offset = 1.0;
        int imax = ( max_order < 0 ) ? level : ( ( max_order-2 < level ) ? max_order -2 : level );
        for( int j=0; j < imax; j++ ){
                mod *= 2;
                offset = 2.0 * offset + 1.0;
                z = (point+1) % mod;
                if ( z < mod / 2 ){
                        value *= ( x - offset + 2.0 * ( (double) z) ) / ( - offset + 2.0 * ( (double) z) );
                }else{
                        value *= ( x + offset - 2.0 * ((double) ( mod - 1 - z )) ) / ( offset - 2.0 * ((double) ( mod - 1 - z )) );
                }
        }
        return value;
}

}

#endif
