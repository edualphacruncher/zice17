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

#ifndef __TSG_BASE_CLASS_CPP
#define __TSG_BASE_CLASS_CPP

#include "tsgGridCore.hpp"

namespace TasGrid{

BaseCanonicalGrid::BaseCanonicalGrid(){}
BaseCanonicalGrid::~BaseCanonicalGrid(){}

int BaseCanonicalGrid::getNumDimensions() const{ return 0; }
int BaseCanonicalGrid::getNumOutputs() const{ return 0; }
TypeOneDRule BaseCanonicalGrid::getRule() const{ return rule_none; }

int BaseCanonicalGrid::getNumLoaded() const{ return 0; }
int BaseCanonicalGrid::getNumNeeded() const{ return 0; }
int BaseCanonicalGrid::getNumPoints() const{ return 0; }

double* BaseCanonicalGrid::getLoadedPoints() const{ return 0; }
double* BaseCanonicalGrid::getNeededPoints() const{ return 0; }
double* BaseCanonicalGrid::getPoints() const{ return 0; }

double* BaseCanonicalGrid::getQuadratureWeights() const{ return 0; }
double* BaseCanonicalGrid::getInterpolationWeights( const double x[] ) const{ return 0; }

void BaseCanonicalGrid::loadNeededPoints( const double *vals ){}

void BaseCanonicalGrid::evaluate( const double x[], double y[] ) const{}
void BaseCanonicalGrid::integrate( double q[] ) const{}

void BaseCanonicalGrid::clearRefinement(){}

SplitDirections::SplitDirections( const IndexSet *points ) : num_dimensions(0), num_allocated_jobs(0), num_jobs(0), job_directions(0), num_job_pnts(0), job_pnts(0) {
        int num_points = points->getNumIndexes();
        num_dimensions = points->getNumDimensions();
        int *map = new int[ num_points * num_dimensions ];  std::fill( map, map + num_points * num_dimensions, 0 );

        num_allocated_jobs = 256;
        job_directions = new int[num_allocated_jobs];
        num_job_pnts   = new int[num_allocated_jobs];
        job_pnts       = new int*[num_allocated_jobs];

        // This section is taking too long! Find a way to parallelize it!
        for( int s=0; s<num_points; s++ ){
                const int *p = points->getIndex( s );
                for( int d=0; d<num_dimensions; d++ ){
                        if ( map[s*num_dimensions + d] == 0 ){
                                if ( num_jobs == num_allocated_jobs ){
                                        num_allocated_jobs *= 2;
                                        int *tmp = job_directions;
                                        job_directions = new int[num_allocated_jobs];  std::copy( tmp, tmp + num_jobs, job_directions );
                                        delete[] tmp;
                                        tmp = num_job_pnts;
                                        num_job_pnts = new int[num_allocated_jobs];  std::copy( tmp, tmp + num_jobs, num_job_pnts );
                                        delete[] tmp;
                                        int **ttmp = job_pnts;
                                        job_pnts = new int*[num_allocated_jobs];  std::copy( ttmp, ttmp + num_jobs, job_pnts );
                                        delete[] ttmp;
                                }

                                job_directions[num_jobs] = d;

                                int nump_allocated = 128;
                                job_pnts[num_jobs] = new int[nump_allocated];
                                num_job_pnts[num_jobs] = 0;
                                for( int i=0; i<num_points; i++ ){
                                        if ( (map[i*num_dimensions + d] == 0) && ( doesBelongSameLine( p, points->getIndex(i), d ) ) ){
                                                // adding a new point
                                                if ( num_job_pnts[num_jobs] == nump_allocated ){
                                                        nump_allocated *= 2;
                                                        int *tmp = job_pnts[num_jobs];
                                                        job_pnts[num_jobs] = new int[nump_allocated];  std::copy( tmp, tmp + num_job_pnts[num_jobs], job_pnts[num_jobs] );
                                                        delete[] tmp;
                                                }

                                                job_pnts[num_jobs][num_job_pnts[num_jobs]] = i;
                                                num_job_pnts[num_jobs]++;
                                                map[i*num_dimensions + d] = -1;
                                        }
                                }

                                num_jobs++;
                        }
                }
        }

        delete[] map;
}
SplitDirections::~SplitDirections(){
        for( int s=0; s<num_jobs; s++ ){
                delete[] job_pnts[s];
        }
        delete[] job_pnts;
        delete[] num_job_pnts;
        delete[] job_directions;
}

int SplitDirections::getNumJobs() const{  return num_jobs;  }
int SplitDirections::getJobDirection( int job ) const{  return job_directions[job];  }
int SplitDirections::getJobNumPoints( int job ) const{  return num_job_pnts[job];  }
const int* SplitDirections::getJobPoints( int job ) const{  return job_pnts[job];  }

bool SplitDirections::doesBelongSameLine( const int a[], const int b[], int direction ) const{
        for( int i=0; i<num_dimensions; i++ ){
                if ( (i != direction) && (a[i] != b[i]) ) return false;
        }
        return true;
}

}

#endif
