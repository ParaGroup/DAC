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
 
 Handmade version of quicksort parallelized with TBB

*/

#include <iostream>
#include <functional>
#include <algorithm>
#include <tbb/task_scheduler_init.h>
#include <tbb/task.h>

#include "../includes/utils.h"
using namespace std;
#define CUTOFF 2000



//definition of the task class that encapsulate the algorithm

class QuickSort: public tbb::task{

public:
    QuickSort(int *array, int left, int right): _array(array),_left(left),_right(right)
    {}



private:


    //sort the array from index left to index right (included)
    tbb::task* execute()
    {
		if(!(_right-_left<=CUTOFF))
		{
			int p=divide();

			//spawn the tasks
			this->set_ref_count(3);
			tbb::task *t =new (allocate_child()) QuickSort(_array, _left,p);
			spawn(*t);

			tbb::task *t2 =new (allocate_child()) QuickSort(_array, p+1,_right);
			spawn_and_wait_for_all(*t2);

		}
		else
			seq();

		return nullptr;
    }


    /*
     * The divide chooses as pivot the middle element and redistributes the elements
     */
	int divide()
    {
		int pivot=_array[(_left+_right)/2];
		int i = _left-1, j = _right+1;
		int tmp;

		while(true)
		{
			do{
				i++;
			}while(_array[i]<pivot);
			do{
				j--;
			}while(_array[j]>pivot);

			if(i>=j)
			   break;

			//swap
			tmp=_array[i];
			_array[i]=_array[j];
			_array[j]=tmp;
		}

		//j is the pivot
		return j;
    }



    /*
     * Base case: we resort on std::sort
     */
	void seq()
    {
		std::sort(&(_array[_left]),&(_array[_right+1]));

    }



    int *_array;
    int _left;
    int _right;
};

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

    tbb::task_scheduler_init _task_scheduler(nwork);		//needed to set par degree
    QuickSort *qs=new (tbb::task::allocate_root())QuickSort(numbers, 0, num_elem-1);

    long start_t=current_time_usecs();

    //create and spawn the first task
    tbb::task::spawn_root_and_wait(*qs);

    long end_t=current_time_usecs();


    if(!isArraySorted(numbers,num_elem))
    {
		fprintf(stderr,"Error: array is not sorted!!\n");
		exit(-1);
    }
    printf("Time (usecs): %Ld\n",end_t-start_t);

    return 0;
}
