/*
 * Code Author: Miroslav Stoyanov
 *
 * Copyright (C) 2015  Miroslav Stoyanov
 *   Large amount of the sparse matrix code is due to Drayton Munster
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

#ifndef __TASMANIAN_LINEAR_SOLVERS_HPP
#define __TASMANIAN_LINEAR_SOLVERS_HPP

#include <list>
#include <vector>
#include <fstream>
#include <math.h>

#include "tsgEnumerates.hpp"

namespace TasGrid{

class TasmanianDenseSolver{
public:
        TasmanianDenseSolver();
        ~TasmanianDenseSolver();

        static void solveLeastSquares( int n, int m, const double A[], const double b[], double reg, double *x );
};

class TasmanianTridiagonalSolver{
public:
        TasmanianTridiagonalSolver();
        ~TasmanianTridiagonalSolver();

        static void decompose( int n, double d[], double e[], double z[] );
};

}

// Sparse linear algebra is used by the wavelet grids only
// tsgEnumerates.hpp is needed by TasSparse for the hardcoded constants, i.e., TSG_DROP_TOL
namespace TasSparse{

using std::endl;

// the enumerates defined here are used only by the sparse code
enum TsgSparseSorting{
	unsorted,
	csc_sort,
	csr_sort
};

enum TsgSparseType{
	coo,
	csc,
	csr,
	none
};

enum TsgCgStatus{
	converged,
	max_iter_reached
};

// (I,J,V) triplet used to represent a sparse matrix entry where A[i,j] += v.
typedef struct COO_Triplet {
	unsigned int i;
	unsigned int j;
	double v;
	inline bool operator==(COO_Triplet &other) const{
		return (i == other.i) && (j == other.j);
	}
} COO_Triplet;

// Sparse matrix represented by above triplets.
class TsgSparseCOO {
public:
	TsgSparseCOO();
	TsgSparseCOO(int in_m, int in_n);
	~TsgSparseCOO();

	void addPoint(unsigned int i, unsigned int j, double v);
	void clear();
	int getNumRows();
	int getNumCols();
	void combine(TsgSparseCOO &m);
	friend class TsgSparseCSC;
	friend class TsgSparseCSR;
	friend class SparseMatrix;

protected:
	std::list<COO_Triplet> values;
	static bool csrComp(COO_Triplet a, COO_Triplet b);
	static bool cscComp(COO_Triplet a, COO_Triplet b);
	void sort(TsgSparseSorting type);
	void coalesce();
	int m, n, nnz;
private:
	TsgSparseSorting sorted;

}; // End TsgSparseCOO

// Make a new matrix class that would store the ILU preconditioner and solve either the regular or adjoined problem all in one class
class SparseMatrix{
public:
        SparseMatrix(TsgSparseCOO &M);
        SparseMatrix(std::ifstream &ifs);
        ~SparseMatrix();

        void write(std::ofstream &ofs) const;
	bool read(std::ifstream &ifs);

        int getNumRows() const;

	void solve( const double b[], double x[], bool transposed = false ) const;

protected:
        void clear();

        void computeILU();

private:
        const double tol;
        int num_rows;
        int *pntr, *indx, *indxD;
        double *vals, *ilu;
};

}

#endif
