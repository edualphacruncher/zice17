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

#ifndef __TASMANIAN_TASGRID_FUNCTIONS_HPP
#define __TASMANIAN_TASGRID_FUNCTIONS_HPP

#include "math.h"

class BaseFunction{
public:
        BaseFunction();
        ~BaseFunction();

        virtual int getNumInputs() const;
        virtual int getNumOutputs() const;
        virtual const char* getDescription() const;
        virtual void eval( const double x[], double y[] ) const;
        virtual void getIntegral( double y[] ) const;
};

class OneOneP0: public BaseFunction{
public:
        OneOneP0(); ~OneOneP0();

        int getNumInputs() const; int getNumOutputs() const; const char* getDescription() const; void eval( const double x[], double y[] ) const; void getIntegral( double y[] ) const;
};

class OneOneP3: public BaseFunction{
public:
        OneOneP3(); ~OneOneP3();

        int getNumInputs() const; int getNumOutputs() const; const char* getDescription() const; void eval( const double x[], double y[] ) const; void getIntegral( double y[] ) const;
};

class OneOneP4: public BaseFunction{
public:
        OneOneP4(); ~OneOneP4();

        int getNumInputs() const; int getNumOutputs() const; const char* getDescription() const; void eval( const double x[], double y[] ) const; void getIntegral( double y[] ) const;
};

class OneOneExpMX: public BaseFunction{
public:
        OneOneExpMX(); ~OneOneExpMX();

        int getNumInputs() const; int getNumOutputs() const; const char* getDescription() const; void eval( const double x[], double y[] ) const; void getIntegral( double y[] ) const;
};

class TwoOneP4: public BaseFunction{
public:
        TwoOneP4(); ~TwoOneP4();

        int getNumInputs() const; int getNumOutputs() const; const char* getDescription() const; void eval( const double x[], double y[] ) const; void getIntegral( double y[] ) const;
};

class TwoOneP5: public BaseFunction{
public:
        TwoOneP5(); ~TwoOneP5();

        int getNumInputs() const; int getNumOutputs() const; const char* getDescription() const; void eval( const double x[], double y[] ) const; void getIntegral( double y[] ) const;
};

class TwoOneExpNX2: public BaseFunction{
public:
        TwoOneExpNX2(); ~TwoOneExpNX2();

        int getNumInputs() const; int getNumOutputs() const; const char* getDescription() const; void eval( const double x[], double y[] ) const; void getIntegral( double y[] ) const;
};

class ThreeOneExpNX2: public BaseFunction{
public:
        ThreeOneExpNX2(); ~ThreeOneExpNX2();

        int getNumInputs() const; int getNumOutputs() const; const char* getDescription() const; void eval( const double x[], double y[] ) const; void getIntegral( double y[] ) const;
};

class TwoOneCos: public BaseFunction{
public:
        TwoOneCos(); ~TwoOneCos();

        int getNumInputs() const; int getNumOutputs() const; const char* getDescription() const; void eval( const double x[], double y[] ) const; void getIntegral( double y[] ) const;
};

class TwoOneSinSin: public BaseFunction{
public:
        TwoOneSinSin(); ~TwoOneSinSin();

        int getNumInputs() const; int getNumOutputs() const; const char* getDescription() const; void eval( const double x[], double y[] ) const; void getIntegral( double y[] ) const;
};

class TwoOneCosCos: public BaseFunction{
public:
        TwoOneCosCos(); ~TwoOneCosCos();

        int getNumInputs() const; int getNumOutputs() const; const char* getDescription() const; void eval( const double x[], double y[] ) const; void getIntegral( double y[] ) const;
};

class TwoOneExpm40: public BaseFunction{
public:
        TwoOneExpm40(); ~TwoOneExpm40();

        int getNumInputs() const; int getNumOutputs() const; const char* getDescription() const; void eval( const double x[], double y[] ) const; void getIntegral( double y[] ) const;
};

class FiveOneExpSum: public BaseFunction{
public:
        FiveOneExpSum(); ~FiveOneExpSum();

        int getNumInputs() const; int getNumOutputs() const; const char* getDescription() const; void eval( const double x[], double y[] ) const; void getIntegral( double y[] ) const;
};

class SixOneExpSum: public BaseFunction{
public:
        SixOneExpSum(); ~SixOneExpSum();

        int getNumInputs() const; int getNumOutputs() const; const char* getDescription() const; void eval( const double x[], double y[] ) const; void getIntegral( double y[] ) const;
};

class EightOneCosSum: public BaseFunction{
public:
        EightOneCosSum(); ~EightOneCosSum();

        int getNumInputs() const; int getNumOutputs() const; const char* getDescription() const; void eval( const double x[], double y[] ) const; void getIntegral( double y[] ) const;
};

class ThreeOneUnitBall: public BaseFunction{
public:
        ThreeOneUnitBall(); ~ThreeOneUnitBall();

        int getNumInputs() const; int getNumOutputs() const; const char* getDescription() const; void eval( const double x[], double y[] ) const; void getIntegral( double y[] ) const;
};

// Function to test special integration rules, i.e. integrals against weight functions
class TwoOneConstGC1: public BaseFunction{
public:
        TwoOneConstGC1(); ~TwoOneConstGC1();

        int getNumInputs() const; int getNumOutputs() const; const char* getDescription() const; void eval( const double x[], double y[] ) const; void getIntegral( double y[] ) const;
};
class TwoOneConstGC2: public BaseFunction{
public:
        TwoOneConstGC2(); ~TwoOneConstGC2();

        int getNumInputs() const; int getNumOutputs() const; const char* getDescription() const; void eval( const double x[], double y[] ) const; void getIntegral( double y[] ) const;
};
class TwoOneConstGG: public BaseFunction{
public:
        TwoOneConstGG(); ~TwoOneConstGG();

        int getNumInputs() const; int getNumOutputs() const; const char* getDescription() const; void eval( const double x[], double y[] ) const; void getIntegral( double y[] ) const;
};
class TwoOneConstGJ: public BaseFunction{
public:
        TwoOneConstGJ(); ~TwoOneConstGJ();

        int getNumInputs() const; int getNumOutputs() const; const char* getDescription() const; void eval( const double x[], double y[] ) const; void getIntegral( double y[] ) const;
};
class TwoOneConstGGL: public BaseFunction{
public:
        TwoOneConstGGL(); ~TwoOneConstGGL();

        int getNumInputs() const; int getNumOutputs() const; const char* getDescription() const; void eval( const double x[], double y[] ) const; void getIntegral( double y[] ) const;
};
class TwoOneConstGH: public BaseFunction{
public:
        TwoOneConstGH(); ~TwoOneConstGH();

        int getNumInputs() const; int getNumOutputs() const; const char* getDescription() const; void eval( const double x[], double y[] ) const; void getIntegral( double y[] ) const;
};

class TwoOneENX2aniso: public BaseFunction{
public:
        TwoOneENX2aniso(); ~TwoOneENX2aniso();

        int getNumInputs() const; int getNumOutputs() const; const char* getDescription() const; void eval( const double x[], double y[] ) const; void getIntegral( double y[] ) const;
};
class SixteenOneActive3: public BaseFunction{
public:
        SixteenOneActive3(); ~SixteenOneActive3();

        int getNumInputs() const; int getNumOutputs() const; const char* getDescription() const; void eval( const double x[], double y[] ) const; void getIntegral( double y[] ) const;
};

class TwoOneDivisionAnisotropic: public BaseFunction{
public:
        TwoOneDivisionAnisotropic(); ~TwoOneDivisionAnisotropic();

        int getNumInputs() const; int getNumOutputs() const; const char* getDescription() const; void eval( const double x[], double y[] ) const; void getIntegral( double y[] ) const;
};
class TwoOne1DCurved: public BaseFunction{
public:
        TwoOne1DCurved(); ~TwoOne1DCurved();

        int getNumInputs() const; int getNumOutputs() const; const char* getDescription() const; void eval( const double x[], double y[] ) const; void getIntegral( double y[] ) const;
};

#endif
