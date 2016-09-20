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


 Backend implementation for the DAC pattern in OpenMP
*/

#ifndef DAC_OPENMP_HPP
#define DAC_OPENMP_HPP
#include <vector>
#include <functional>
#include <omp.h>

// This is a first implementation prone to optimizations (especially considering cutoff that for the moment being is statically found)

template<typename OperandType,typename ResultType>
class DacOpenmp{

public:

	DacOpenmp(const std::function<void(const OperandType&,std::vector<OperandType>&)>& divide_fn,
			  const std::function<void(std::vector<ResultType>&,ResultType&)>& combine_fn,
			  const std::function<void(const OperandType&, ResultType&)>& seq_fn,
			  const std::function<bool(const OperandType&)>& cond_fn, const OperandType& op, ResultType& res, int pardegree):
				_divide_fn(divide_fn), _combine_fn(combine_fn), _seq_fn(seq_fn), _condition_fn(cond_fn), _op(&op), _res(&res), _pardegree(pardegree)
	{}

	void compute()
	{
		//call recursive DAC
#pragma omp parallel num_threads(_pardegree)
#pragma omp single
		recursiveDac(_op,_res);
	}


private:

	void recursiveDac(const OperandType *op, ResultType *ret)
	{

		if(!_condition_fn(*op)) //not the base case
		{
			//divide
			std::vector<OperandType> *ops=new std::vector<OperandType>();
			_divide_fn(*op,*ops);
			int branch_factor=ops->size();

			//create the space for the partial results
			std::vector<ResultType> *ress=new std::vector<ResultType>(branch_factor);

			//create recursive tasks
			for(int i=0;i<branch_factor;i++)
			{
#pragma omp task
				{
					recursiveDac(&(*ops)[i],&(*ress)[i]);
				}
			}
#pragma omp taskwait



			//combine results
			_combine_fn(*ress,*ret);

			//cleanup memory

			delete ops;
			delete ress;
		}
		else
		{
			_seq_fn(*op,*ret);
		}

	}

	//function pointers
	const std::function<void(const OperandType&,std::vector<OperandType>&)>& _divide_fn;
	const std::function<void(const OperandType& ,  ResultType&)>& _seq_fn;
	const std::function<void(std::vector<ResultType>&,ResultType&)>& _combine_fn;
	const std::function<bool(const OperandType&)>& _condition_fn;
	const OperandType* _op;
	ResultType* _res;

	int _pardegree;
};

#endif // DAC_OPENMP_HPP
