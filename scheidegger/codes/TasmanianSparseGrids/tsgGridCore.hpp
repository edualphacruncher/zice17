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

#ifndef __TSG_BASE_CLASS_HPP
#define __TSG_BASE_CLASS_HPP

#include "tsgEnumerates.hpp"
#include "tsgIndexSets.hpp"

namespace TasGrid{

class BaseCanonicalGrid{
public:
        BaseCanonicalGrid();
        virtual ~BaseCanonicalGrid();

        virtual int getNumDimensions() const;
        virtual int getNumOutputs() const;
        virtual TypeOneDRule getRule() const;

        virtual int getNumLoaded() const;
        virtual int getNumNeeded() const;
        virtual int getNumPoints() const;

        virtual double* getLoadedPoints() const;
        virtual double* getNeededPoints() const;
        virtual double* getPoints() const;

        virtual double* getQuadratureWeights() const;
        virtual double* getInterpolationWeights( const double x[] ) const;

        virtual void loadNeededPoints( const double *vals );

        virtual void evaluate( const double x[], double y[] ) const;
        virtual void integrate( double q[] ) const;

        virtual void clearRefinement();
};

class SplitDirections{
public:
        SplitDirections( const IndexSet *points );
        ~SplitDirections();

        int getNumJobs() const;
        int getJobDirection( int job ) const;
        int getJobNumPoints( int job ) const;
        const int* getJobPoints( int job ) const;

protected:
        bool doesBelongSameLine( const int a[], const int b[], int direction ) const;

private:
        int num_dimensions, num_allocated_jobs, num_jobs;
        int *job_directions, *num_job_pnts, **job_pnts;
};

}

#endif
