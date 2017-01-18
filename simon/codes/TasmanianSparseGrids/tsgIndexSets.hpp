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

#ifndef __TASMANIAN_SPARSE_GRID_INDEX_SETS_HPP
#define __TASMANIAN_SPARSE_GRID_INDEX_SETS_HPP

#include "tsgEnumerates.hpp"

namespace TasGrid{

class UnsortedIndexSet{ // assumes the elements would be added in a sorted order
public:
        UnsortedIndexSet( int cnum_dimensions, int cnum_slots );
        ~UnsortedIndexSet();

        int getNumDimensions() const;
        int getNumIndexes() const;

        void addIndex( const int p[] );
        const int* getIndex( int i ) const;

        int* getIndexesSorted() const;
        // returns an array of num_dimensions X num_indexes of the sorted indexes

protected:
        TypeIndexRelation compareIndexes( const int a[], const int b[] ) const;
        void merge( const int listA[], int sizeA, const int listB[], int sizeB, int destination[] ) const;

private:
        int num_dimensions;
        int num_indexes;
        int *index;
};

class GranulatedIndexSet{ // optimized to add one index at a time
public:
        GranulatedIndexSet( int cnum_dimensions, int cnum_slots = 64 );
        ~GranulatedIndexSet();

        int getNumDimensions() const;
        int getNumIndexes() const;

        void addIndex( const int p[] );
        int getSlot( const int p[] ) const;

        const int* getIndex( int j ) const;

        void addUnsortedSet( const UnsortedIndexSet *set );
        void addGranulatedSet( const GranulatedIndexSet *set );

        const int* getIndexes() const;
        const int* getMap() const;

protected:
        TypeIndexRelation compareIndexes( const int a[], const int b[] ) const;

        void merge( const int newIndex[], int sizeNew );
        void mergeMapped( const int newIndex[], const int newMap[], int sizeNew );

private:
        int num_dimensions;
        int num_indexes;
        int num_slots;
        int *index;
        int *map;
};

class IndexSet{ // ridgit set but optimal in storage size
public:
        IndexSet( int cnum_dimensions );
        IndexSet( const UnsortedIndexSet *set );
        IndexSet( const GranulatedIndexSet *set );
        IndexSet( const IndexSet *set );
        IndexSet( int cnum_dimensions, const int cindex[] );
        IndexSet( int cnum_dimensions, int cnum_indexes, int* &cindex );
        ~IndexSet();

        int getNumDimensions() const;
        int getNumIndexes() const;

        void write( std::ofstream &ofs ) const;
        void read( std::ifstream &ifs );

        int getSlot( const int p[] ) const;

        const int* getIndex( int i ) const;

        void addUnsortedSet( const UnsortedIndexSet *set );
        void addGranulatedSet( const GranulatedIndexSet *set );
        void addIndexSet( const IndexSet *set );

        IndexSet* diffSets( const IndexSet *set ) const; // returns this set minus the points in set

protected:
        TypeIndexRelation compareIndexes( const int a[], const int b[] ) const;

        void merge( const int newIndex[], int sizeNew );
        void mergeMapped( const int newIndex[], const int map[], int sizeNew );

private:
        int num_dimensions;
        int num_indexes;
        int *index;
};

class StorageSet{ // stores the values of the function
public:
        StorageSet( int cnum_outputs, int cnum_values );
        ~StorageSet();

        void write( std::ofstream &ofs ) const;
        void read( std::ifstream &ifs );

        const double* getValues( int i ) const;
        int getNumOutputs() const;

        void setValues( const double vals[] );
        void setValue( int i, const double val[] );
        void addValues( const IndexSet *old_set, const IndexSet *new_set, const double new_vals[] );

protected:
        TypeIndexRelation compareIndexes( int num_dimensions, const int a[], const int b[] ) const;

private:
        int num_outputs, num_values;
        double *values;
};

}

#endif
