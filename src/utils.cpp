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


 Set of utilities
*/
#include <stdlib.h>
#include <iostream>
#include "../includes/utils.h"

using namespace std;
int *generateRandomArray(int n)
{
    srand ((time(0)));
    int *numbers=new int[n];
    for(int i=0;i<n;i++)
        numbers[i]=(int) (rand()) / ((RAND_MAX/MAX_NUM));
    return numbers;
}
void printArray(int *a, int n)
{
    for(int i=0;i<n;i++)
        cout << a[i]<<endl;
}

//simple check
bool isArraySorted(int *a, int n)
{
    for(int i=1;i<n;i++)
        if(a[i]<a[i-1])
            return false;
    return true;
}





double **generateRandomMatrix(int n)
{
    srand ((time(0)));
    double **matrix=allocateMatrix(n);

    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
                matrix[i][j]=(double) (rand()) / ((double)(RAND_MAX/MAX_INT_NUM));
            //matrix[i][j]=i*n+j;
        }

    }
    return matrix;
}


void printMatrix(double **a,int n)
{
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            std::cout << a[i][j] <<"\t";
        }
        std::cout <<endl;
    }

}
double **matmul(double **a, double **b, int n)
{
    //allocate space for the result
    double **c=allocateMatrix(n);

    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            c[i][j]=0;
            for(int k=0;k<n;k++)
                c[i][j]+=a[i][k]*b[k][j];
        }
    }
    return c;
}

void addMatrix(double **a, double **b,double **c,int n)
{
    for(int i=0;i<n;i++)
        for(int j=0;j<n;j++)
            c[i][j]=a[i][j]+b[i][j];
}

void subtMatrix(double **a,double **b,double **c,int n)
{
    for(int i=0;i<n;i++)
        for(int j=0;j<n;j++)
            c[i][j]=a[i][j]-b[i][j];
}

bool areMatrixEqual(double **a, double **b, int n)
{
    for(int i=0;i<n;i++)
        for(int j=0;j<n;j++)
            if(a[i][j]!=b[i][j])
            {
                std::cout << a[i][j]<<"!="<<b[i][j]<<endl;
                return false;
            }
    return true;
}


// COMPACT VERSION: matrix allocate on a single continuous space, with row stride=n
double *generateCompactRandomMatrix(int n)
{
    srand ((time(0)));
    double *matrix=allocateCompactMatrix(n);

    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
                matrix[i*n+j]=(double) (rand()) / ((double)(RAND_MAX/MAX_DBL_NUM));
        }

    }
    return matrix;
}



void printCompactMatrix(double *a, int n, int row_stride)
{
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            std::cout << a[i*row_stride+j] <<"\t";
        }
        std::cout <<endl;
    }
}

double *compactMatmul(double *a, int rs_a, double *b, int rs_b, int n)
{
    //allocate space for the result
    double *c=allocateCompactMatrix(n);

    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            c[i*n+j]=0;
            for(int k=0;k<n;k++)
                c[i*n+j]+=a[i*rs_a+k]*b[k*rs_b+j];
        }
    }
    return c;
}



void addCompactMatrix(double *a, int rs_a, double *b, int rs_b, double *c, int rs_c, int n)
{
    for(int i=0;i<n;i++)
        for(int j=0;j<n;j++)
            c[i*rs_c+j]=a[i*rs_a+j]+b[i*rs_b+j];
}

//C=A-B
void subtCompactMatrix(double *a, int rs_a, double *b, int rs_b, double *c, int rs_c, int n)
{
    for(int i=0;i<n;i++)
        for(int j=0;j<n;j++)
            c[i*rs_c+j]=a[i*rs_a+j]-b[i*rs_b+j];
}

bool areCompactMatrixEqual(double *a, double *b, int n)
{

    for(int i=0;i<n;i++)
        for(int j=0;j<n;j++)
            if((a[i*n+j]-b[i*n+j])>THRESHOLD)
            {
                std::cout << a[i*n+j]<<"!="<<b[i*n+j]<<endl;
                return false;
            }
    return true;
}
