/* -*- Mode: C++; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ***************************************************************************
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License version 2 as 
 *  published by the Free Software Foundation.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 *
 *  As a special exception, you may use this file as part of a free software
 *  library without restriction.  Specifically, if other files instantiate
 *  templates or use macros or inline functions from this file, or you compile
 *  this file and link it with other files to produce an executable, this
 *  file does not by itself cause the resulting executable to be covered by
 *  the GNU General Public License.  This exception does not however
 *  invalidate any other reasons why the executable file might be covered by
 *  the GNU General Public License.
 *
 ****************************************************************************
 */

/* 
 * Author: Massimo Torquati <torquati@di.unipi.it> 
 * Date:   August 2014
 */

/*
 * NxM matrix multiplication using the Divide and Conquer approach to solve 
 * the Strassen algorithm O(n^2.8074).
 * This algorithm is particularly suitable for shared-cache multi-core since
 * recursion allows to better exploit cache hierarchy.
 *
 * command:  
 *          DCstrassen <M> <N> <P> [base_case_size] [pfworkers:chuncksize] [check]
 *          <-> required argument, [-] optional argument
 *
 * example:               
 *          DCstrassen 2048 2048 2048 $((256*256*256)) 4:0 0
 * 
 *  With the above command, it executes a square N by M matrix multiplication 
 *  with N = 2048, using the ParallelFor pattern with 4 workers (and using the 
 *  default static scheduling policy). No result check is executed (since last
 *  argument is 0). The recursion base case is 256*256*256.
 * 
 */

#include <cassert>
#include <cmath>
#include <cstdio>  
#include <cstdlib> 
#include <cstring>
#include <sys/time.h>

#include <string>

#include <tbb/task_scheduler_init.h>
#include <tbb/task.h>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>

using namespace tbb;

const double THRESHOLD = 0.001;
unsigned long long GRAIN = 512*512*512; // recursion base case size
bool  check = false;                    // set to true for checking the result
bool PFSPINWAIT=true;                   // enabling spinWait
int  PFWORKERS=1;                       // parallel_for parallelism degree
int  PFGRAIN  =0;                       // default static scheduling of iterations


void random_init (long M, long N, long P, double *A, double *B) {
    for (long i = 0; i < M; i++) {  
        for (long j = 0; j < P; j++) {  
            A[i*P+j] = 5.0 - ((double)(rand()%100) / 10.0); 
        }
    }
	
    for (long i = 0; i < P; i++) {  
        for (long j = 0; j < N; j++) {  
            B[i*N+j] = 5.0 - ((double)(rand()%100) / 10.0);
        }     
    }
}

long parse_arg(const char *string) {
    char *endptr;
    long r = strtol(string, &endptr, 10);
    assert(endptr[0]==0); // used up the whole string
    assert(r>0);          // need a positive integer
    return r;
}
//C99 unsigned long long
unsigned long long parse_argull(const char *string) {
    char *endptr;
    unsigned long long r = strtoull(string, &endptr, 10);
    assert(endptr[0]==0); // used up the whole string
    assert(r>0);          // need a positive integer
    return r;
}

void printarray(const double *A, long m, long n, long N) {
    for (long i=0; i<m; i++) {
	for (long j=0; j<n; j++)
	    printf("%f\t", A[i*N+j]);
	printf("\n");
    }
}

// triple nested loop (ijk) implementation 
inline void seqMatMultPF( long m, long n, long p, 
                       const double* A, const long AN, 
                       const double* B, const long BN, 
                       double* C, const long CN)  {   

    tbb::parallel_for(tbb::blocked_range<long>(0,m),
                      [&] (const tbb::blocked_range<long>& r) {
                          for (long i = r.begin(); i != r.end(); ++i) {
                              for (long j = 0; j < n; j++) {
                                  C[i*CN+j] = 0.0;  
                                  for (long k = 0; k < p; k++)  
                                      C[i*CN+j] += A[i*AN+k]*B[k*BN+j];  
                              }                                
                          }
                      }
                      );
} 
inline void seqMatMult(long m, long n, long p, 
                       const double* A, const long AN, 
                       const double* B, const long BN, 
                       double* C, const long CN)  {   
    for (long i = 0; i < m; i++)  
        for (long j = 0; j < n; j++) {
            C[i*CN+j] = 0.0;  
            for (long k = 0; k < p; k++)  
                C[i*CN+j] += A[i*AN+k]*B[k*BN+j];  
        }  
} 

// m by n with row stride XN for X YN for Y and CN for C
void mmsum ( const double *X, long XN, const double *Y, long YN,
            double *C, long CN,  long m, long n) {
    for(long i=0;i<m;++i) {
            for (long j=0; j<n; j++) 
                C[i*CN+j] = X[i*XN+j] + Y[i*YN + j];
        }
   
}

// m by n with row stride XN for X YN for Y and CN for C
void mmsub ( const double *X, long XN,  const double *Y, long YN, 
	    double *C, long CN,  long m, long n) {
    for(long i=0;i<m;++i) {
            for (long j=0; j<n; j++)
                C[i*CN+j] = X[i*XN+j] - Y[i*YN + j];
        }

}


class strassenMMult : public tbb::task {
	long m, n, p;
    const double *A;
    const long AN;
    const double *B;
    const long BN;
    double *C;
    const long CN;
public:
    strassenMMult ( long m, long n, long p,
                    const double *A, const long AN,
                    const double *B, const long BN,
                    double *C, const long CN): m(m),n(n),p(p),
                                               A(A),AN(AN),
                                               B(B),BN(BN),
                                               C(C),CN(CN) { }
	
	tbb::task* execute();		
};



/* 
 * Strassen algorithm: 
 *
 *  S1  = A11 + A22
 *  S2  = B11 + B22
 *  P1  = S1 * S2
 *  S3  = A21 + A22
 *  P2  = S3 * B11
 *  S4  = B12 - B22
 *  P3  = A11 * S4
 *  S5  = B21 - B11
 *  P4  = A22 * S5
 *  S6  = A11 + A12
 *  P5  = S6 * B22
 *  S7  = A21 - A11
 *  S8  = B11 + B12
 *  P6  = S7 * S8
 *  S9  = A12 - A22
 *  S10 = B21 + B22
 *  P7  = S9*S10
 *  C11 = P1 + P4 - P5 + P7
 *  C12 = P3 + P5
 *  C21 = P2 + P4
 *  C22 = P1 - P2 + P3 + P6
 *
 */
// Version with reduced memory usage. 
// Instead of using temporary arrays (P2, P3, P6 and P7)
// it uses the C matrices
tbb::task* strassenMMult::execute() {

    if ( (m==1) || (n==1) || (p==1) ||
         (((unsigned long long)m*n*p) < GRAIN) ) {
        seqMatMultPF(m,n,p, A, AN, B, BN, C, CN);
	} else {
        long m2 = m/2;
        long n2 = n/2;
        long p2 = p/2;    
   
        const double *A11 = &A[0];
        const double *A12 = &A[p2];
        const double *A21 = &A[m2*AN];
        const double *A22 = &A[m2*AN+p2];
        
        const double *B11 = &B[0];
        const double *B12 = &B[n2];
        const double *B21 = &B[p2*BN];
        const double *B22 = &B[p2*BN+n2];
        
        double *C11 = &C[0];
        double *C12 = &C[n2];
        double *C21 = &C[m2*CN];
        double *C22 = &C[m2*CN+n2];
        
        double *P1  = (double*)malloc(m2*n2*sizeof(double));
        double *P2  = (double*)malloc(m2*n2*sizeof(double));
        double *P3  = (double*)malloc(m2*n2*sizeof(double));
        double *P4  = (double*)malloc(m2*n2*sizeof(double));
        double *P5  = (double*)malloc(m2*n2*sizeof(double));
        double *P6  = (double*)malloc(m2*n2*sizeof(double));
        double *P7  = (double*)malloc(m2*n2*sizeof(double));		
        
        double *sumA1= (double*)malloc(m2*p2*sizeof(double));
        double *sumB1= (double*)malloc(p2*n2*sizeof(double));
        double *sumA2= (double*)malloc(m2*p2*sizeof(double));
        double *sumB3= (double*)malloc(p2*n2*sizeof(double));
        double *sumB4= (double*)malloc(p2*n2*sizeof(double));
        double *sumA5= (double*)malloc(m2*p2*sizeof(double));
        double *sumA6= (double*)malloc(m2*p2*sizeof(double));
        double *sumB6= (double*)malloc(p2*n2*sizeof(double));
        double *sumA7= (double*)malloc(m2*p2*sizeof(double));
        double *sumB7= (double*)malloc(p2*n2*sizeof(double));

        int count = 8;
		this->set_ref_count(count);

        mmsum( A11, AN, A22, AN, sumA1, p2, m2, p2);               // S1
        mmsum( B11, BN, B22, BN, sumB1, n2, p2, n2);               // S2
        tbb::task *t1 = new (allocate_child()) strassenMMult( m2,n2,p2, sumA1, p2, sumB1, n2, P1, n2);   // P1
        spawn(*t1);

        mmsum( A21, AN, A22, AN, sumA2, p2, m2, p2);               // S3 
        tbb::task *t2 = new (allocate_child()) strassenMMult( m2,n2,p2, sumA2, p2, B11, BN, P2, n2);    // P2
        spawn(*t2);

        mmsub( B12, BN, B22, BN, sumB3, n2, p2, n2);               // S4
        tbb::task *t3 = new (allocate_child()) strassenMMult( m2,n2,p2, A11, AN, sumB3, n2, P3, n2);    // P3
        spawn(*t3);

        mmsub( B21, BN, B11, BN, sumB4, n2, p2, n2);               // S5
        tbb::task *t4 = new (allocate_child()) strassenMMult( m2,n2,p2, A22, AN, sumB4, n2, P4, n2);    // P4
        spawn(*t4);

        mmsum( A11, AN, A12, AN, sumA5, p2, m2, p2);               // S6
        tbb::task *t5 = new (allocate_child()) strassenMMult( m2,n2,p2, sumA5, p2, B22, BN, P5, n2);    // P5
        spawn(*t5);
        
        mmsub( A21, AN, A11, AN, sumA6, p2, m2, p2);               // S7
        mmsum( B11, BN, B12, BN, sumB6, n2, p2, n2);               // S8
        tbb::task *t6 = new (allocate_child()) strassenMMult( m2, n2, p2, sumA6, p2, sumB6, n2, P6, n2); // P6
        spawn(*t6);
        
        mmsub( A12, AN, A22, AN, sumA7, p2, m2, p2);               // S9
        mmsum( B21, BN, B22, BN, sumB7, n2, p2, n2);               // S10
        tbb::task *t7 = new (allocate_child()) strassenMMult( m2, n2, p2, sumA7, p2, sumB7, n2, P7, n2); // P7
        
        spawn_and_wait_for_all(*t7);

        tbb::parallel_for(tbb::blocked_range<long>(0,m2),
                          [&] (const tbb::blocked_range<long>& r) {
                              for (long i = r.begin(); i != r.end(); ++i) {
                                  for(long j=0; j<n2; ++j) {
                                      C11[i*CN + j] = P1[i*n2 + j] + P4[i*n2 + j] - P5[i*n2 + j] + P7[i*n2 + j];
                                      C12[i*CN + j] = P3[i*n2 + j] + P5[i*n2 + j];
                                      C21[i*CN + j] = P2[i*n2 + j] + P4[i*n2 + j];
                                      C22[i*CN + j] = P1[i*n2 + j] - P2[i*n2 + j] + P3[i*n2 + j] + P6[i*n2 + j];
                                  }                                  
                              }
                                                   }
                          );
        
        free(P1); free(P2); free(P3); free(P4); free(P5); free(P6); free(P7);
        free(sumA1); free(sumB1);
        free(sumA2);
        free(sumB3);
        free(sumB4);
        free(sumA5);
        free(sumA6); free(sumB6);        
        free(sumA7); free(sumB7);
    }
    return nullptr;
}
void startWrapper ( long m, long n, long p,
                    const double *A, const long AN,
                    const double *B, const long BN,
                    double *C, const long CN) {     
    
    strassenMMult &t = *new (tbb::task::allocate_root()) strassenMMult(m,n,p,A,AN,B,BN,C,CN);
    tbb::task::spawn_root_and_wait(t);
}



long CheckResults(long m, long n, const double *C, const double *C1)
{
	for (long i=0; i<m; i++)
		for (long j=0; j<n; j++) {
			long idx = i*n+j;
			if (fabs(C[idx] - C1[idx]) > THRESHOLD) {
				printf("ERROR %ld,%ld %f != %f\n", i, j, C[idx], C1[idx]);
				return 1;
			}
		}
	printf("OK.\n");
	return 0;
}
 
void get_time( void (*F)( long m, long n, long p, const double* A, long AN, const double* B, long BN, double* C, long CN),
                long m, long n, long p, const double* A, long AN, const double* B, long BN, double *C,
	       const char *descr)
{
	printf("Executing %-40s", descr);
	fflush(stdout);
	struct timeval before,after;
	gettimeofday(&before, NULL);

	F( m, n, p, A, AN, B, BN, C, BN);
	gettimeofday(&after,  NULL);

	double tdiff = after.tv_sec - before.tv_sec + (1e-6)*(after.tv_usec - before.tv_usec);
	printf("  Done in secs: %11.6f\n", tdiff);

}

void get_time_and_check(void (*F)( long m, long n, long p, const double* A, long AN, const double* B, long BN, double* C, long CN),
                         long m, long n, long p, const double* A, long AN, const double* B, long BN, const double *C_expected, double *C,
			const char *descr)
{
    get_time(F,  m, n, p, A, AN, B, BN, C, descr);
    CheckResults(m, n, C_expected, C);
}



int main(int argc, char* argv[])  {     
    if (argc < 4) {
        printf("\n\tuse: %s <M> <N> <P> [base_case_size=512*512*512] [pfworkers:pfgrain=1:0] [check=0]\n", argv[0]);
        printf("\t       <-> required argument, [-] optional argument\n");
        printf("\t       A is M by P\n");
        printf("\t       B is P by N\n");
        printf("\t       base_case_size is the base case for the recursion (default 512*512*512)\n");
        printf("\t       pfworkers is the  n. of workers of the ParallelFor pattern (default 1)\n");
        printf("\t       pfgrain is the ParallelFor grain size (0, default static scheduling)\n");
        printf("\t       check!=0 executes also the standard ijk algo for checking the result\n\n");
        printf("\tNOTE: M, N and P should be even.\n\n");
        return -1;
    }
    long M = parse_arg(argv[1]);
    long N = parse_arg(argv[2]);
    long P = parse_arg(argv[3]);
    if (argc >= 5) GRAIN=parse_argull(argv[4]);
    if (argc >= 6) {
        std::string pfarg(argv[5]);
        int n = pfarg.find_first_of(":");
        if (n>0) {
            PFWORKERS = atoi(pfarg.substr(0,n).c_str());
            PFGRAIN   = atoi(pfarg.substr(n+1).c_str());
        } else PFWORKERS = atoi(argv[4]);
    }
    if (argc >= 7) check = (atoi(argv[6])?true:false);

    const double *A = (double*)malloc(M*P*sizeof(double));
    const double *B = (double*)malloc(P*N*sizeof(double));
    assert(A); assert(B);

    random_init(M, N, P, const_cast<double*>(A), const_cast<double*>(B));
    double *C       = (double*)malloc(M*N*sizeof(double));
  
	printf("Executing with %d threads\n",PFWORKERS);
    tbb::task_scheduler_init init(PFWORKERS);

  
    if (check) {
        double *C2   = (double*)malloc(M*N*sizeof(double));
        get_time(seqMatMult,  M,N,P, A, P, B, N, C, "Standard ijk algorithm");
        get_time_and_check(startWrapper,  M, N, P, A, P, B, N, C, C2, "Strassen");
        free(C2);
    } else 
        get_time(startWrapper,  M, N, P, A, P, B, N, C, "Strassen");
    
    free((void*)A); free((void*)B); free(C);

    return 0;  
} 
