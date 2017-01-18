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

#ifndef __TSG_RULE_LOCAL_POLYNOMIAL_HPP
#define __TSG_RULE_LOCAL_POLYNOMIAL_HPP

#include "tsgEnumerates.hpp"

namespace TasGrid{

class BaseRuleLocalPolynomial{
public:
        BaseRuleLocalPolynomial();
        ~BaseRuleLocalPolynomial();

        virtual int getMaxOrder() const;
        virtual void setMaxOrder( int order );

        virtual TypeOneDRule getType() const;

        virtual int getNumPoints( int level ) const; // for building the initial grid

        virtual double getNode( int point ) const;

        virtual int getLevel( int point ) const;

        virtual bool isSemiLocal() const; // second level for quadratic functions, two ways

        virtual int getParent( int point ) const;
        virtual int getStepParent( int point ) const;
        virtual int getKidLeft( int point ) const;
        virtual int getKidRight( int point ) const;

        virtual double evalRaw( int point, double x ) const; // normalizes x (i.e., (x-node) / support ), but it does not check the support
        virtual double evalSupport( int point, double x, bool &isSupported ) const; // // normalizes x (i.e., (x-node) / support ) and checks if x is within the support

        virtual double getArea( int point, int n, const double w[], const double x[] ) const;
        // integrate the function associated with the point, constant to cubic are known analytically, higher order need a 1-D quadrature rule

        static int intlog2( int i );
        static int int2log2( int i );

protected:
        int max_order;
};

class RuleLocalPolynomial : public BaseRuleLocalPolynomial{
public:
        RuleLocalPolynomial();
        ~RuleLocalPolynomial();

        int getMaxOrder() const;
        void setMaxOrder( int order );

        TypeOneDRule getType() const;

        int getNumPoints( int level ) const;

        double getNode( int point ) const;

        int getLevel( int point ) const;

        bool isSemiLocal() const;

        int getParent( int point ) const;
        int getStepParent( int point ) const;
        int getKidLeft( int point ) const;
        int getKidRight( int point ) const;

        double evalRaw( int point, double x ) const;
        double evalSupport( int point, double x, bool &isSupported ) const;
        double getArea( int point, int n, const double w[], const double x[] ) const;

protected:
        double scaleX( int point, double x ) const;
        double getSupport( int point ) const;
        // those take a normalized x
        double evalPWQuadratic( int point, double x ) const;
        double evalPWCubic( int point, double x ) const;
        double evalPWPower( int point, double x ) const;
};

class RuleSemiLocalPolynomial : public BaseRuleLocalPolynomial{
public:
        RuleSemiLocalPolynomial();
        ~RuleSemiLocalPolynomial();

        int getMaxOrder() const;
        void setMaxOrder( int order );

        TypeOneDRule getType() const;

        int getNumPoints( int level ) const;

        double getNode( int point ) const;

        int getLevel( int point ) const;

        bool isSemiLocal() const;

        int getParent( int point ) const;
        int getStepParent( int point ) const;
        int getKidLeft( int point ) const;
        int getKidRight( int point ) const;

        double evalRaw( int point, double x ) const;
        double evalSupport( int point, double x, bool &isSupported ) const;
        double getArea( int point, int n, const double w[], const double x[] ) const;

protected:
        double scaleX( int point, double x ) const;
        double getSupport( int point ) const;
        // those take a normalized x
        double evalPWQuadratic( int point, double x ) const;
        double evalPWCubic( int point, double x ) const;
        double evalPWPower( int point, double x ) const;
};

class RuleLocalPolynomialZero : public BaseRuleLocalPolynomial{
public:
        RuleLocalPolynomialZero();
        ~RuleLocalPolynomialZero();

        int getMaxOrder() const;
        void setMaxOrder( int order );

        TypeOneDRule getType() const;

        int getNumPoints( int level ) const;

        double getNode( int point ) const;

        int getLevel( int point ) const;

        bool isSemiLocal() const;

        int getParent( int point ) const;
        int getStepParent( int point ) const;
        int getKidLeft( int point ) const;
        int getKidRight( int point ) const;

        double evalRaw( int point, double x ) const;
        double evalSupport( int point, double x, bool &isSupported ) const;
        double getArea( int point, int n, const double w[], const double x[] ) const;

protected:
        double scaleX( int point, double x ) const;
        double getSupport( int point ) const;
        // those take a normalized x
        double evalPWCubic( int point, double x ) const;
        double evalPWPower( int point, double x ) const;
};

}

#endif
