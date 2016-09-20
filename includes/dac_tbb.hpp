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


 Backend implementation of the DAC pattern for the Intel TBB framework
*/

#ifndef DAC_TBB_HPP
#define DAC_TBB_HPP

#include <vector>
#include <functional>
#include <tbb/task_scheduler_init.h>
#include <tbb/task.h>



/**
   The pattern is provided by  using the same API of the dac_mdf (or very similar)
	but exploiting  tbb::task (as indicated by Intel if we are looking for the best performance and scalability)
	[1] - https://www.threadingbuildingblocks.org/docs/help/reference/task_scheduler/catalog_of_recommended_task_patterns.html
*/



template<typename OperandType,typename ResultType>
class DacTask :public tbb::task{

public:
	DacTask(const std::function<void(const OperandType&,std::vector<OperandType>&)>& divide_fn,
			const std::function<void(std::vector<ResultType>&,ResultType&)>& combine_fn,
			const std::function<void(const OperandType&, ResultType&)>& seq_fn,
			const std::function<bool(const OperandType&)>& cond_fn, const OperandType* op, ResultType* res ):
			  _divide_fn(divide_fn), _combine_fn(combine_fn), _seq_fn(seq_fn), _condition_fn(cond_fn), _op(op), _res(res)
	{
	}


	//execute method required by tbb
	tbb::task* execute()
	{
		if(!_condition_fn(*_op)) //not the base case
		{
			//divide
			std::vector<OperandType> *ops=new std::vector<OperandType>();
			_divide_fn(*_op,*ops);
			int branch_factor=ops->size();

			//create the space for the partial results
			std::vector<ResultType> *ress=new std::vector<ResultType>(branch_factor);

			//create the tasks
			//The call to set_ref_count uses k+1 as its argument. The extra 1 is critical. (source [1])
			this->set_ref_count(branch_factor+1);
			for(int i=0;i<branch_factor-1;i++)
			{
				tbb::task *t =new (allocate_child()) DacTask(_divide_fn,_combine_fn,_seq_fn,_condition_fn,&(*ops)[i],&(*ress)[i]);
				spawn(*t);
			}
			//last one
			tbb::task *t =new (allocate_child()) DacTask(_divide_fn,_combine_fn,_seq_fn,_condition_fn,&(*ops)[branch_factor-1],&(*ress)[branch_factor-1]);
			spawn_and_wait_for_all(*t);


			//combine results
			_combine_fn(*ress,*_res);

			//cleanup memory

			delete ops;
			delete ress;
		}
		else
		{
			_seq_fn(*_op,*_res);
		}
		return nullptr;
	}

private:

	const std::function<void(const OperandType&,std::vector<OperandType>&)>& _divide_fn;
	const std::function<void(const OperandType& ,  ResultType&)>& _seq_fn;
	const std::function<void(std::vector<ResultType>&,ResultType&)>& _combine_fn;
	const std::function<bool(const OperandType&)>& _condition_fn;
	const OperandType* _op;
	ResultType* _res;


};


template<typename OperandType,typename ResultType>
class DacTBB  {

public:

	DacTBB(const std::function<void(const OperandType&,std::vector<OperandType>&)>& divide_fn,
			  const std::function<void(std::vector<ResultType>&,ResultType&)>& combine_fn,
			  const std::function<void(const OperandType&, ResultType&)>& seq_fn,
			  const std::function<bool(const OperandType&)>& cond_fn, const OperandType& op, ResultType& res, int pardegree):
				_divide_fn(divide_fn), _combine_fn(combine_fn), _seq_fn(seq_fn), _condition_fn(cond_fn), _op(&op), _res(&res), _pardegree(pardegree), _task_scheduler(pardegree)
	{


	}

	void compute()
	{
		//create the first task
		DacTask<OperandType,ResultType> *dac=new (tbb::task::allocate_root()) DacTask<OperandType,ResultType>(_divide_fn,_combine_fn,_seq_fn,_condition_fn,_op,_res);
		tbb::task::spawn_root_and_wait(*dac);


	}


private:


	//function pointers
	const std::function<void(const OperandType&,std::vector<OperandType>&)>& _divide_fn;
	const std::function<void(const OperandType& ,  ResultType&)>& _seq_fn;
	const std::function<void(std::vector<ResultType>&,ResultType&)>& _combine_fn;
	const std::function<bool(const OperandType&)>& _condition_fn;
	const OperandType* _op;
	ResultType* _res;
	int _pardegree;
	tbb::task_scheduler_init _task_scheduler;		//needed to set par degree
};


#endif // DAC_TBB_HPP
