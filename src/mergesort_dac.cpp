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


 Mergesort: sort an array of N integer in parallel using the DAC pattern
  and C++11 semantics (iterator)


*/
#include <iostream>
#include <functional>
#include <vector>
#include <algorithm>
#include <cstring>
#include "../includes/utils.h"
#if USE_FF
#include <ff/DC.hpp>
using namespace ff;
#endif
#if USE_OPENMP
#include "../includes/dac_openmp.hpp"
#endif
#if USE_TBB
#include "../includes/dac_tbb.hpp"
#endif
using namespace std;
#define CUTOFF 2000


// Operand (i.e. the Problem) and Results share the same format
struct ops{
	vector<int>::iterator left;
	vector<int>::iterator right;
};

typedef struct ops Operand;
typedef struct ops Result;


/*
 * The divide simply 'split' the array in two: the splitting is only logical.
 * The recursion occur on the left and on the right part
 */
void divide(const Operand &op,std::vector<Operand> &subops)
{
	vector<int>::iterator mid=op.left+(op.right-op.left)/2;
	Operand a;
	a.left=op.left;
	a.right=mid;
	subops.push_back(a);

	Operand b;
	b.left=mid;
	b.right=op.right;
	subops.push_back(b);
}


/*
 * For the base case we resort to std::sort
 */
void seq(const Operand &op, Result &ret)
{
	ret=op;
	std::sort(ret.left,ret.right);
}


/*
 * The Merge (Combine) function start from two ordered sub array and construct the original one
 * It uses additional memory
 */
void mergeMS(vector<Result>&ress, Result &ret)
{
	//compute what is needed: array pointer, mid, ...
	vector<int>::iterator i=ress[0].left;
	vector<int>::iterator mid=ress[0].right;
	vector<int>::iterator j=mid;
	int size=ress[1].right-ress[0].left;
	vector<int> tmp(size);

	//merge in order
	for(int k=0;k<size;k++)
	{
		if(i<mid && (j>=ress[1].right || *i<=*j))
		{
			tmp[k]=*i;
			i++;
		}
		else
		{
			tmp[k]=*j;
			j++;
		}
	}

	//copy back

	std::copy(tmp.begin(),tmp.end(),ress[0].left);
	//build the result
	ret.left=ress[0].left;
	ret.right=ress[1].right;
}


/*
 * Base case condition
 */
bool cond(const Operand &op)
{
	return (op.right-op.left<=CUTOFF);
}

//simple check
bool isVectorSorted(vector<int> a, int n)
{
	for(int i=1;i<n;i++)
		if(a[i]<a[i-1])
			return false;
	return true;
}

int main(int argc, char *argv[])
{
	if(argc<2)
	{
		cerr << "Usage: "<<argv[0]<< " <num_elements> <num_workers>"<<endl;
		exit(-1);
	}
	std::function<void(const Operand&,vector<Operand>&)> div(divide);
	std::function <void(const Operand &,Result &)> sq(seq);
	std::function <void(vector<Result>&,Result &)> mergef(mergeMS);
	std::function<bool(const Operand &)> cf(cond);

	int num_elem=atoi(argv[1]);
	int nwork=atoi(argv[2]);
	//generate a random array
	int *numbers=generateRandomArray(num_elem);
	//fill the vector
	vector<int> v(numbers, numbers+num_elem ); // use some utility to avoid hardcoding the size here

	//build the operand
	Operand op;

	op.left=v.begin();
	op.right=v.end();

	Result res;
#if USE_FF
	ff_DC<Operand, Result> dac(div,mergef,sq,cf,op,res,nwork);
#endif
#if USE_OPENMP
	DacOpenmp<Operand, Result> dac(div,mergef,sq,cf,op,res,nwork);
#endif
#if USE_TBB
	DacTBB<Operand, Result> dac(div,mergef,sq,cf,op,res,nwork);
#endif

	long start_t=current_time_usecs();

	//compute
#if USE_FF
	dac.run_and_wait_end();
#else
	dac.compute();
#endif
	long end_t=current_time_usecs();


	//Correctness check
	if(!isVectorSorted(v,num_elem))
	{
		fprintf(stderr,"Error: array is not sorted!!\n");
		exit(-1);
	}
	printf("Time (usecs): %Ld\n",end_t-start_t);

	return 0;
}

