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

 Handmade version of quicksort parallelized with openmp

*/

#include <iostream>
#include <functional>
#include <algorithm>
#include "../includes/utils.h"
#include <omp.h>
using namespace std;
#define CUTOFF 2000



/*
 * The divide chooses as pivot the middle element and redistributes the elements
 */
int divide(int *a, int left, int right)
{
    int pivot=a[(left+right)/2];
    int i = left-1, j = right+1;
    int tmp;

    while(true)
    {
		do{
			i++;
		}while(a[i]<pivot);
		do{
			j--;
		}while(a[j]>pivot);

		if(i>=j)
		   break;

		//swap
		tmp=a[i];
		a[i]=a[j];
		a[j]=tmp;
    }

    //j is the pivot
    return j;
}



/*
 * Base case: we resort on std::sort
 */
void seq(int *a, int left, int right)
{
    std::sort(&(a[left]),&(a[right+1]));

}


//sort the array from index left to index right (included)
void quicksort(int *a, int left, int right)
{
    if(!(right-left<=CUTOFF))
    {
	int p=divide(a,left,right);
#pragma omp task
	quicksort(a,left,p);
#pragma omp task
	quicksort(a,p+1,right);
#pragma omp taskwait
    }
    else
	seq(a,left,right);
}

int main(int argc, char *argv[])
{
    if(argc<2)
    {
		cerr << "Usage: "<<argv[0]<< " <num_elements> <num_workers>"<<endl;
		exit(-1);
    }

    int num_elem=atoi(argv[1]);
    int nwork=atoi(argv[2]);
    int *numbers=generateRandomArray(num_elem);

    //build the operand
    //Operand *op=new Operand();

    long start_t=current_time_usecs();

    //compute
#pragma omp parallel num_threads(nwork)
#pragma omp single
    quicksort(numbers,0,num_elem-1);

    long end_t=current_time_usecs();

    if(!isArraySorted(numbers,num_elem))
    {
		fprintf(stderr,"Error: array is not sorted!!\n");
		exit(-1);
    }
    printf("Time (usecs): %Ld\n",end_t-start_t);

    return 0;
}
