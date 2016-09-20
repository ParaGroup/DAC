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

 Fibonacci: computes the n-th number of the fibonacci sequence
*/

#include <iostream>
#include <functional>
#include <vector>
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

/*
 * Problem and Result are just integers
 */

/*
 * Divide Function: recursively compute n-1 and n-2
 */
void divide(const unsigned int &op,std::vector<unsigned int> &subops)
{
	subops.push_back(op-1);
	subops.push_back(op-2);
}

/*
 * Base Case
 */
void seq(const unsigned int &op, unsigned int &res)
{
	res=1;
}

/*
 * Combine function
 */
void combine(vector<unsigned int>& res, unsigned int &ret)
{
	ret=res[0]+res[1];
}

/*
 * Condition for base case
 */
bool cond(const unsigned int &op)
{
	return (op<=2);
}



int main(int argc, char *argv[])
{

	if(argc<3)
	{
		fprintf(stderr,"Usage: %s <N> <pardegree>\n",argv[0]);
		exit(-1);
	}
	unsigned int start=atoi(argv[1]);
	int nwork=atoi(argv[2]);

	unsigned int res;

	//lambda version just for testing it
#if USE_FF
	ff_DC<unsigned int, unsigned int> dac(
				[](const unsigned int &op,std::vector<unsigned int> &subops){
						subops.push_back(op-1);
						subops.push_back(op-2);
					},
				[](vector<unsigned int>& res, unsigned int &ret){
						ret=res[0]+res[1];
					},
				[](const unsigned int &op, unsigned int &res){
						res=1;
					},
				[](const unsigned int &op){
						return (op<=2);
					},
				start,
				res,
				nwork
				);
#endif
#if USE_OPENMP
	DacOpenmp<unsigned int, unsigned int> dac(divide,combine,seq,cond,start,res,nwork);
#endif
#if USE_TBB
	DacTBB<unsigned int, unsigned int> dac(divide,combine,seq,cond,start,res,nwork);
#endif

    long start_t=current_time_usecs();

	//compute
#if USE_FF
	dac.run_and_wait_end();
#else
	dac.compute();
#endif
    long end_t=current_time_usecs();
	printf("Result: %d\n",res);
	printf("Time (usecs): %Ld\n",end_t-start_t);

}
