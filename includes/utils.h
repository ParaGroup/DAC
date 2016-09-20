/* ***************************************************************************
 *
 *  This program is free software; you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License version 3 as
 *  published by the Free Software Foundation.
 *
 *  This program is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program; if not, write to the Free Software Foundation,
 *  Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 *
 ****************************************************************************


 Author: Tiziano De Matteis <dematteis@di.unipi.it>


 Set of utilities used by the various programs (e.g. allocation of arrays, matrices, ...)
*/
#ifndef UTILS_H
#define UTILS_H
#include <time.h>
#include <sys/time.h>

/**
    Array functions
*/

//maximum value for arrays elements
#define MAX_NUM 99999.9
#define MAX_INT_NUM 999
#define MAX_DBL_NUM 999.9

//threshold for comparing double matrices
const double THRESHOLD = 0.001;

//Array functions

int *generateRandomArray(int n);
void printArray(int *a, int n);

//simple check
bool isArraySorted(int *a,int n);


//Matrix function
double **generateRandomMatrix(int n);
void printMatrix(double **a,int n);
double **matmul(double **a, double **b, int n);

//size n*n
double *generateCompactRandomMatrix(int n);
void printCompactMatrix(double *a, int n, int row_stride);
double *compactMatmul(double *a, int rs_a, double *b, int rs_b, int n);

//C=A+B
void addMatrix(double **a, double **b,double **c,int n);
//C=A-B
void subtMatrix(double **a,double **b,double **c,int n);

bool areMatrixEqual(double **a, double **b, int n);

//Asssuming square matrices, we have to consider the row stride of the operand
//and of the results
//C=A+B
/**
 * @brief addCompactMatrix add two matrices
 * @param a first operand
 * @param rs_a row stride of the first operand
 * @param b second operand
 * @param rs_b row stride of the first operand
 * @param c result
 * @param rs_c row stride
 * @param n matrices size (assuming square)
 */
void addCompactMatrix(double *a, int rs_a, double *b,int rs_b,double *c,int rs_c,int n);
//C=A-B
//using row_strides
void subtCompactMatrix(double *a, int rs_a, double *b, int rs_b, double *c, int rs_c, int n);


//assuming row stride n
bool areCompactMatrixEqual(double *a, double *b, int n);

//Matrix functions
inline double **allocateMatrix(int n) __attribute__((always_inline));
/**
 * @brief allocateMatrix allocate the space for a matrix of size n*n
 * @return
 */
inline double **allocateMatrix(int n)
{
    double **matrix=new double*[n]();
    for(int i=0;i<n;i++)
    {
        matrix[i]=new double[n]();
    }
    return matrix;
}

/**
 * @brief allocateCompactMatrix allocate the matrix on a single continuous memory space of size n*n (row stripe=n)
 * @param n
 * @return a contiguos memory area
 */
inline double *allocateCompactMatrix(int n) __attribute__((always_inline));
inline double *allocateCompactMatrix(int n)
{
    return new double[n*n]();
}

inline void deallocateMatrix(double **m,int n) __attribute__((always_inline));
inline void deallocateMatrix(double **m,int n)
{
    for(int i=0;i<n;i++)
        delete []m[i];
    delete []m;
}

inline void deallocateCompactMatrix(double *m,int n) __attribute__((always_inline));
inline void deallocateCompactMatrix(double *m,int n)
{
    delete []m;
}
/**
    Timing functions
*/
inline long current_time_usecs() __attribute__((always_inline));
inline long current_time_usecs(){
    struct timeval t;
    gettimeofday(&t, NULL);
    return (t.tv_sec)*1000000L + t.tv_usec;

}

inline long current_time_nsecs() __attribute__((always_inline));
inline long current_time_nsecs(){
    struct timespec t;
    clock_gettime(CLOCK_REALTIME, &t);
    return (t.tv_sec)*1000000000L + t.tv_nsec;
}

inline int isPowerOfTwo (unsigned int x)
{
  return ((x != 0) && !(x & (x - 1)));
}

#endif // UTILS_H
