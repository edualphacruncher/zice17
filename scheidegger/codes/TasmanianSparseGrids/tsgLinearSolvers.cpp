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

#ifndef __TASMANIAN_LINEAR_SOLVERS_CPP
#define __TASMANIAN_LINEAR_SOLVERS_CPP

#include "tsgLinearSolvers.hpp"

namespace TasGrid{

TasmanianDenseSolver::TasmanianDenseSolver(){}
TasmanianDenseSolver::~TasmanianDenseSolver(){}

void TasmanianDenseSolver::solveLeastSquares( int n, int m, const double A[], const double b[], double reg, double *x ){
        // form Ar = A' * A
        // form x = A' * b
        double* Ar = new double[ m * m ];
        #pragma omp parallel for
        for( int i=0; i<m; i++ ){
                for( int j=0; j<m; j++ ){
                        Ar[ i*m + j ] = 0.0;
                        for( int k=0; k<n; k++ ){
                                Ar[ i*m + j ] += A[ i*n + k ] * A[ j*n + k ];
                        }
                }
                x[i] = 0.0;
                for( int k=0; k<n; k++ ){
                        x[i] += A[ i*n+k ] * b[k];
                }
        }

        // add regularization
        for( int i=0; i<m; i++ ){  Ar[ i *m + i ] += reg;  };

        // factorize Ar
        for( int i=0; i<m; i++ ){
                Ar[ i*m + i ] = sqrt( Ar[ i*m + i ] );

                for( int j=i+1; j<m; j++ ){
                        Ar[ i*m + j ] /= Ar[ i*m + i ];
                }

                for( int k=i+1; k<m; k++ ){
                        for( int j=i+1; j <= k; j++ ){
                                Ar[ j*m + k ] -= Ar[ i*m + j ] * Ar[ i*m + k ];
                        }
                }
        }

        // solve L^(-1) x
        for( int i=0; i<m; i++ ){
                x[i] /= Ar[ i*m + i ];
                for( int j=i+1; j<m; j++ ){
                        x[j] -= x[i] * Ar[ i*m + j ];
                }
        }

        // solve L^(-T) x
        for( int i=m-1; i>=0; i-- ){
                for( int j=i+1; j<m; j++ ){
                        x[i] -= x[j] * Ar[ i*m + j ];
                }
                x[i] /= Ar[ i*m + i ];
        }

        delete[] Ar;
}

TasmanianTridiagonalSolver::TasmanianTridiagonalSolver(){}
TasmanianTridiagonalSolver::~TasmanianTridiagonalSolver(){}

void TasmanianTridiagonalSolver::decompose( int n, double d[], double e[], double z[] ){
        const double tol = TSG_NUM_TOL;
        if ( n == 1 ){ z[0] = z[0]*z[0]; return; }

        for( int l=0; l<n-1; l++ ){
                int m = l;
                while( (m < n-1) && (fabs(e[m]) > tol) ) m++;

                while ( m != l ){
                        double p = d[l];
                        double g = ( d[l+1] - p ) / ( 2.0 * e[l] );
                        double r = sqrt( g*g + 1.0 );

                        g = d[m] - p + e[l] / ( g + ( (g>=0) ? 1.0 : -1.0 ) * r ); // sign function here may be unstable

                        double s = 1.0;
                        double c = 1.0;
                        p = 0.0;

                        for( int i=m-1; i>=l; i-- ){
                                double f = s * e[i];
                                double b = c * e[i];

                                if ( fabs( f ) >= fabs( g ) ){
                                        c = g / f;
                                        r = sqrt( c*c + 1.0 );
                                        e[i+1] = f*r;
                                        s = 1.0 / r;
                                        c *= s;
                                }else{
                                        s = f / g;
                                        r =  sqrt( s*s + 1.0 );
                                        e[i+1] = g * r;
                                        c = 1.0 / r;
                                        s *= c;
                                }

                                g = d[i+1] - p;
                                r = ( d[i] - g ) * s + 2.0 * c * b;
                                p = s * r;
                                d[i+1] = g + p;
                                g = c * r - b;
                                f = z[i+1];
                                z[i+1] = s * z[i] + c * f;
                                z[i] = c * z[i] - s * f;
                        }

                        d[l] -= p;
                        e[l] = g;
                        e[m] = 0.0;

                        m = l;
                        while( (m < n-1) && (fabs(e[m]) > tol) ) m++;
                }
        }

        for( int i=1; i<n; i++ ){
                for( int j=0; j<n-1; j++ ){
                        if ( d[j] > d[j+1] ){
                                double p = d[j];
                                d[j] = d[j+1];
                                d[j+1] = p;
                                p = z[j];
                                z[j] = z[j+1];
                                z[j+1] = p;
                        }
                }
        }
        for( int i=0; i<n; i++ ){
               z[i] = z[i]*z[i];
        }
}

}

namespace TasSparse{

/* BEGIN TsgSparseCOO */

TsgSparseCOO::TsgSparseCOO() : m(0), n(0), nnz(0), sorted(unsorted){}

TsgSparseCOO::TsgSparseCOO(int in_m, int in_n) : sorted(unsorted) {
	m = in_m;
	n = in_n;
	nnz = 0;
}

TsgSparseCOO::~TsgSparseCOO() {
	clear();
}

void TsgSparseCOO::combine(TsgSparseCOO &mat){
	values.splice(values.end(), mat.values);
	nnz += mat.nnz;
	mat.clear();
}

void TsgSparseCOO::clear(){
	sorted = unsorted;
	values.clear();
	nnz = 0;
}

bool TsgSparseCOO::csrComp(const COO_Triplet a, const COO_Triplet b){
	return (a.i < b.i) || (a.i == b.i && a.j < b.j);
}

bool TsgSparseCOO::cscComp(const COO_Triplet a, const COO_Triplet b){
	return (a.j < b.j) || (a.j == b.j && a.i < b.i);
}

void TsgSparseCOO::addPoint(unsigned int i, unsigned int j, double v){
	/*
	 * Adds the given point to the sparse matrix. Note: assumes the entry is inside the
	 * bounds defined in the constructor (i.e. i < m, j < n). If this assumption is
	 * violated, conversion routines will result in undefined behavior (likely segfault).
	 */
	sorted = unsorted;
	COO_Triplet t;
	t.i = i;
	t.j = j;
	t.v = v;
	nnz++;
	values.push_back(t);
}

void TsgSparseCOO::sort(TsgSparseSorting type){
	/* Already sorted or no sort specified*/
	if (sorted == type || type == unsorted) return;

	values.sort((type == csc_sort) ? cscComp : csrComp);
	sorted = type;
}

void TsgSparseCOO::coalesce(){
	/* Assumes any equivalent entries are adjacent e.g. sorted list */

//	If only list would pass references rather than copies.
//	/* Combines equivalent elements */
//	values.unique(coalesce_helper);
//
//	/* Removes elements below DROP_TOL */
//	values.remove_if(near_zero);


	std::list<COO_Triplet>::iterator it = values.begin();
	it++;
	std::list<COO_Triplet>::iterator old_elem = values.begin();

	while(it != values.end()){
		/* Coalesce duplicate elements */
		if (*old_elem == *it){
			old_elem->v += it->v;
			it = values.erase(it); /* Removes the list element and advances the iter */
			nnz--;
		} else{
			/* Drop elements below DROP_TOL */
			if (fabs(old_elem->v) < TSG_NUM_TOL){
				old_elem = values.erase(old_elem);
				nnz--;
			}else{
				old_elem++;
			}
			it++;
		}
	}

}
/* END TsgSparseCOO */

SparseMatrix::SparseMatrix(TsgSparseCOO &M) : tol(TSG_NUM_TOL), num_rows(0), pntr(0), indx(0), indxD(0), vals(0), ilu(0) {
//
//      Drayton's code for collapse of COO:
//	Constructs a sparse matrix in CSR (compressed sparse row / compressed row storage)
//	format given a sparse matrix in COO (coordinate) format.
//
	M.sort(csr_sort);
	M.coalesce();

	int length = M.nnz; // MIRO: assume the matrix is square
	num_rows = M.n;
	vals = new double[length];
	indx = new int[length];
	pntr = new int[num_rows+1];

	int current_row = -1;

	std::list<COO_Triplet>::const_iterator it = M.values.begin();
	for(int i = 0; i < length; i++, it++){
		COO_Triplet current = *it;
		while ( ((int)current.i) != current_row){
			/* Either denote the end of the row or add empty rows until we've caught up */
			current_row++;
			pntr[current_row] = i;
		}
		indx[i] = current.j;
		vals[i] = current.v;
	}
	pntr[current_row+1] = length;

	computeILU();
}

SparseMatrix::SparseMatrix(std::ifstream &ifs) : tol(TSG_NUM_TOL), num_rows(0), pntr(0), indx(0), indxD(0), vals(0), ilu(0) {
        read( ifs );
}

SparseMatrix::~SparseMatrix(){
        clear();
}

void SparseMatrix::write(std::ofstream &ofs) const{
        ofs << std::scientific; ofs.precision(17);
        ofs << num_rows << endl;
        ofs << pntr[0]; for( int i=1; i<=num_rows; i++ ){  ofs << " " << pntr[i];  } ofs << endl;
        if ( num_rows > 0 ){
                ofs << indx[0];  for( int i=1; i<pntr[num_rows]; i++ ){  ofs << " " << indx[i];  } ofs << endl;
                ofs << vals[0];  for( int i=1; i<pntr[num_rows]; i++ ){  ofs << " " << vals[i];  } ofs << endl;
                ofs << indxD[0]; for( int i=1; i<num_rows; i++ ){        ofs << " " << indxD[i]; } ofs << endl;
                ofs << ilu[0];   for( int i=1; i<pntr[num_rows]; i++ ){  ofs << " " << ilu[i];   } ofs << endl;
        }
}
bool SparseMatrix::read(std::ifstream &ifs){
        clear();
        ifs >> num_rows;
        if ( num_rows > 0 ){
                pntr = new int[num_rows+1];
                for( int i=0; i<=num_rows; i++ ){  ifs >> pntr[i];  }
                indx = new int[pntr[num_rows]];
                for( int i=0; i<pntr[num_rows]; i++ ){  ifs >> indx[i];  }
                vals = new double[pntr[num_rows]];
                for( int i=0; i<pntr[num_rows]; i++ ){  ifs >> vals[i];  }
                indxD = new int[num_rows];
                for( int i=0; i<num_rows; i++ ){  ifs >> indxD[i];  }
                ilu = new double[pntr[num_rows]];
                for( int i=0; i<pntr[num_rows]; i++ ){  ifs >> ilu[i];  }
        }
        return true;
}

void SparseMatrix::clear(){
        if ( pntr != 0 ){ delete[] pntr; pntr=0; }
        if ( indx != 0 ){ delete[] indx; indx=0; }
        if ( indxD != 0 ){ delete[] indxD; indxD=0; }
        if ( vals != 0 ){ delete[] vals; vals=0; }
        if ( ilu != 0 ){ delete[] ilu; ilu=0; }
        num_rows = 0;
}

int SparseMatrix::getNumRows() const{ return num_rows; }

void SparseMatrix::computeILU(){
        if ( indxD != 0 ){ delete[] indxD; }
        if ( ilu != 0 ){ delete[] ilu; }
        indxD = new int[num_rows];
        ilu = new double[ pntr[num_rows] ];
        for( int i=0; i<num_rows; i++ ){
		int j = pntr[i];
		while( indx[j] < i ){ j++; };
		indxD[i] = j;
	}

	std::copy( vals, vals + pntr[num_rows], ilu );

        for( int i=0; i<num_rows-1; i++ ){
                double u = ilu[ indxD[i] ];
                #pragma omp parallel for
                for( int j=i+1; j<num_rows; j++ ){ // update the rest of the matrix, each row can be done in parallel
                        int jc = pntr[j];
                        while( indx[jc] < i ){ jc++; }
                        if ( indx[jc] == i ){
                                ilu[jc] /= u;
                                double l = ilu[jc];
                                int ik = indxD[i]+1;
                                int jk = jc+1;
                                while( (ik<pntr[i+1]) && (jk<pntr[j+1]) ){
                                        if ( indx[ik] == indx[jk] ){
                                                ilu[jk] -= l * ilu[ik];
                                                ik++; jk++;
                                        }else if ( indx[ik] < indx[jk] ){
                                                ik++;
                                        }else{
                                                jk++;
                                        }
                                }
                        }
                }
        }
}

void SparseMatrix::solve( const double b[], double x[], bool transposed ) const{ // using GMRES
        int max_inner = 30;
        int max_outer = 80;
        double *W = new double[(max_inner+1) * num_rows]; // Krylov basis

        double *H = new double[ max_inner * (max_inner+1) ]; // holds the transformation for the normalized basis
        double *S = new double[ max_inner ]; // sin and cos of the Givens rotations
        double *C = new double[ max_inner+1 ];
        double *Z = new double[ max_inner ]; // holds the coefficients of the solution

        double alpha, h_k; // temp variables

        double outer_res = tol + 1.0, inner_res; // outer and inner residual
        int outer_itr = 0, inner_itr; // counts the inner and outer iterations

        double *pb = new double[num_rows];
        if ( !transposed ){
                std::copy( b, b + num_rows, pb );

                // action of the preconditioner
                for( int i=1; i<num_rows; i++ ){
                        for( int j=pntr[i]; j<indxD[i]; j++){
                                pb[i] -= ilu[j] * pb[ indx[j] ];
                        }
                }
                for( int i=num_rows-1; i>=0; i-- ){
                        for( int j=indxD[i]+1; j<pntr[i+1]; j++ ){
                                pb[i] -= ilu[j] * pb[ indx[j] ];
                        }
                        pb[i] /= ilu[ indxD[i] ];
                }
        }

        std::fill( x, x + num_rows, 0.0 ); // zero initial guess, I wonder if we can improve this

        while ( (outer_res > tol ) && ( outer_itr < max_outer ) ){
                std::fill( W, W + num_rows, 0.0 );
                if ( transposed ){
                        std::copy( x, x + num_rows, pb );
                        for( int i=0; i<num_rows; i++ ){
                                pb[i] /= ilu[indxD[i]];
                                for( int j=indxD[i]+1; j<pntr[i+1]; j++ ){
                                        pb[ indx[j] ] -= ilu[j] * pb[i];
                                }
                        }
                        for( int i=num_rows-2; i>=0; i-- ){
                                for( int j=pntr[i]; j<indxD[i]; j++ ){
                                        pb[ indx[j] ] -= ilu[j] * pb[i];
                                }
                        }
                        for( int i=0; i<num_rows; i++ ){
                                for( int j=pntr[i]; j<pntr[i+1]; j++ ){
                                        W[indx[j]] += vals[j] * pb[i];
                                }
                        }
                        for( int i=0; i<num_rows; i++ ){
                                W[i] = b[i] - W[i];
                        }
                }else{
                        for( int i=0; i<num_rows; i++ ){
                                for( int j=pntr[i]; j<pntr[i+1]; j++ ){
                                        W[i] += vals[j] * x[indx[j]];
                                }
                        }
                        for( int i=1; i<num_rows; i++ ){
                                for( int j=pntr[i]; j<indxD[i]; j++){
                                        W[i] -= ilu[j] * W[ indx[j] ];
                                }
                        }
                        for( int i=num_rows-1; i>=0; i-- ){
                                for( int j=indxD[i]+1; j<pntr[i+1]; j++ ){
                                        W[i] -= ilu[j] * W[ indx[j] ];
                                }
                                W[i] /= ilu[ indxD[i] ];
                        }

                        for( int i=0; i<num_rows; i++ ){
                                W[i] = pb[i] - W[i];
                        }
                }

                Z[0] = 0.0;  for( int i=0; i<num_rows; i++ ){  Z[0] += W[i]*W[i];  };  Z[0] = sqrt( Z[0] );
                for( int i=0; i<num_rows; i++ ){  W[i] /= Z[0]; };

                inner_res = Z[0]; // first residual
                inner_itr = 0; // counts the size of the basis

                while ( (inner_res > tol) && (inner_itr < max_inner-1) ){
                        inner_itr++;

                        std::fill( &(W[inner_itr*num_rows]), &(W[inner_itr*num_rows]) + num_rows, 0.0 );
                        if ( transposed ){
                                std::copy( &(W[num_rows*(inner_itr-1)]), &(W[num_rows*(inner_itr-1)]) + num_rows, pb );
                                for( int i=0; i<num_rows; i++ ){
                                        pb[i] /= ilu[indxD[i]];
                                        for( int j=indxD[i]+1; j<pntr[i+1]; j++ ){
                                                pb[ indx[j] ] -= ilu[j] * pb[i];
                                        }
                                }
                                for( int i=num_rows-2; i>=0; i-- ){
                                        for( int j=pntr[i]; j<indxD[i]; j++ ){
                                                pb[ indx[j] ] -= ilu[j] * pb[i];
                                        }
                                }
                                for( int i=0; i<num_rows; i++ ){
                                        for( int j=pntr[i]; j<pntr[i+1]; j++ ){
                                                W[inner_itr*num_rows + indx[j]] += vals[j] * pb[i];
                                        }
                                }
                        }else{
                                for( int i=0; i<num_rows; i++ ){
                                        for( int j=pntr[i]; j<pntr[i+1]; j++ ){
                                                W[inner_itr*num_rows + i] += vals[j] * W[num_rows*(inner_itr-1) + indx[j]];
                                        }
                                }
                                for( int i=1; i<num_rows; i++ ){
                                        for( int j=pntr[i]; j<indxD[i]; j++){
                                                W[inner_itr*num_rows + i] -= ilu[j] * W[inner_itr*num_rows + indx[j] ];
                                        }
                                }
                                for( int i=num_rows-1; i>=0; i-- ){
                                        for( int j=indxD[i]+1; j<pntr[i+1]; j++ ){
                                                W[inner_itr*num_rows + i] -= ilu[j] * W[inner_itr*num_rows + indx[j] ];
                                        }
                                        W[inner_itr*num_rows + i] /= ilu[ indxD[i] ];
                                }
                        }

                        #pragma omp parallel for
                        for( int i=0; i<inner_itr; i++ ){
                                H[i*max_inner + inner_itr-1] = 0.0; for( int j=0; j<num_rows; j++ ){  H[i*max_inner + inner_itr-1] += W[inner_itr*num_rows+j] * W[i*num_rows+j];  };
                        }

                        #pragma omp parallel for
                        for( int j=0; j<num_rows; j++ ){
                                for( int i=0; i<inner_itr; i++ ){
                                        W[inner_itr*num_rows+j] -= H[i*max_inner + inner_itr-1] * W[i*num_rows+j];
                                }
                        };

                        h_k = 0.0;  for( int i=0; i<num_rows; i++ ){  h_k += W[inner_itr*num_rows+i]*W[inner_itr*num_rows+i];  }; h_k = sqrt( h_k ); //cout << "h_k = " << h_k << endl;
                        for( int i=0; i<num_rows; i++ ){  W[inner_itr*num_rows+i] /= h_k;  };

                        for ( int i=0; i<inner_itr-1; i++ ){ // form the next row of the transformation
                                alpha = H[i*max_inner + inner_itr-1];
                                H[    i*max_inner + inner_itr-1] = C[i] * alpha + S[i] * H[(i+1)*max_inner + inner_itr-1];
                                H[(i+1)*max_inner + inner_itr-1] = S[i] * alpha - C[i] * H[(i+1)*max_inner + inner_itr-1];
                        };

                        alpha = sqrt( h_k * h_k  +  H[(inner_itr-1)*max_inner + inner_itr-1] * H[(inner_itr-1)*max_inner + inner_itr-1] );

                        // set the next set of Givens rotations
                        S[inner_itr-1] = h_k / alpha;
                        C[inner_itr-1] = H[(inner_itr-1)*max_inner + inner_itr-1] / alpha;

                        H[(inner_itr-1) * max_inner + inner_itr-1] = alpha;

                        // Z is used to reconstruct the solution in the end
                        Z[inner_itr] = S[inner_itr-1]*Z[inner_itr-1];
			Z[inner_itr-1] = C[inner_itr-1]*Z[inner_itr-1]; // apply it on z

                        inner_res = fabs(Z[inner_itr]);
                }

                inner_itr--;

                if ( inner_itr > -1 ){ // if the first guess was not within TOL of the true solution
			Z[inner_itr] /= H[inner_itr * max_inner + inner_itr];
			for( int i=inner_itr-1; i>-1; i-- ){
				h_k = 0.0;
				for( int j=i+1; j<=inner_itr; j++ ){
					h_k += H[i*max_inner + j] * Z[j];
				};
				Z[i] = ( Z[i] - h_k ) / H[i * max_inner + i];
			}

			for( int i=0; i<=inner_itr; i++ ){
				for( int j=0; j<num_rows; j++ ){
                                        x[j] += Z[i] * W[i*num_rows+j];
				}
			}

			if ( transposed ){
                                for( int i=0; i<num_rows; i++ ){
                                        x[i] /= ilu[indxD[i]];
                                        for( int j=indxD[i]+1; j<pntr[i+1]; j++ ){
                                                x[ indx[j] ] -= ilu[j] * x[i];
                                        }
                                }
                                for( int i=num_rows-2; i>=0; i-- ){
                                        for( int j=pntr[i]; j<indxD[i]; j++ ){
                                                x[ indx[j] ] -= ilu[j] * x[i];
                                        }
                                }
			}

		}

                outer_res = inner_res;
                outer_itr++;
        }

        delete[] pb;
        delete[] W;

        delete[] H;
        delete[] S;
        delete[] C;
        delete[] Z;
}

} /* namespace TasSparse */

#endif
