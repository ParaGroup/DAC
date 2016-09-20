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



 Implementation of the Intel Parallel Stable Merge Sort using the DAC pattern.
 Source: https://software.intel.com/en-us/articles/a-parallel-stable-sort-using-c11-for-tbb-cilk-plus-and-openmp

 The basic algorithm is the same. We want to show that our solution do not performs poorly wrt intel one
 and requires less coding

*/

#include <iostream>
#include <functional>
#include <vector>
#include <algorithm>
#include <assert.h>
#include <atomic>
#include <iterator>
#include <cstring>
#include <iostream>

#include "../includes/utils.h"

//library taken from intel stable sort implementation
#include <pss_common.h>
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

#define CUTOFF 500	//same value of Intel source code (INTEL)

//---------------------------------------------------------------------
//Types and definition inherithed by test.cpp in parallel stable sort (INTEL)
//---------------------------------------------------------------------


//! Number of extant keys
static std::atomic<int> KeyCount;

//! Highest index in array to be sorted.
static unsigned HighValidIndex;

//! A key to be sorted, with lots of checking.
class Key {
	//! Value used by comparator
	int value;
	//! Original position or special vlaue (Empty or Dead)
	int index;
	//! Special value used to mark constructed object without a comparable value.
	static const int Empty = -1;
	//! Special value used to mark destroyed objects.
	static const int Dead = -2;
public:
	Key() {
		++KeyCount;
		index = Empty;
	}
	Key( const Key& k ) : value(k.value), index(k.index) {
		assert(k.isLive());
		++KeyCount;
	}
	Key( Key&& k ) : value(k.value), index(k.index) {
		assert(k.isConstructed());
		k.index = Empty;
		++KeyCount;
	}
	~Key() {
		assert(isConstructed());

		index = Dead;
		--KeyCount;
	}
	void operator=( Key&& k ) {
		assert(k.isConstructed());
		assert(isConstructed());
		value = k.value;
		index = k.index;
		k.index = Empty;
	}
	void operator=( const Key& k ) {
		assert(k.isLive());
		assert(isConstructed());

		value = k.value;
		index = k.index;
	}
	bool isConstructed() const {
		return isLive() || index==Empty;
	}
	bool isLive() const {
		return unsigned(index)<=HighValidIndex;
	}
	friend class KeyCompare;
	friend void CreateDataSet( size_t n );
	friend void CheckIsSorted( size_t n );
};

class KeyCompare {
public: //edit:  status types are public now
	enum statusType {
		//! Special value used to defined object.
		Live = 0xabcd,
		//! Special value used to mark default-constructed object.
		Empty = -1,
		//! Special value used to mark destroyed objects.
		Dead = -2
	} status;
	KeyCompare( statusType s ) : status(s) {}

	KeyCompare() {status = Empty;}
	~KeyCompare() {status = Dead;}
	bool operator()( const Key& j, const Key& k ) const {
		assert(status==Live);
		assert(j.isLive());
		assert(k.isLive());
		return j.value<k.value;
	}
	friend void Test( size_t n );
};

//! Iterator that points to a key, with some checking.
class Iterator {
	Key* my_ptr;
	//! [my_begin,my_end) is the legal range for my_ptr if this Iterator has defined value.
	Key *my_begin, *my_end;
	static Key* poison(long val) {
		return (Key*)val;
	}
	bool isDefined() const {
		assert(my_begin<=my_ptr);
		assert(my_ptr<=my_end);
		return true;
	}
public:
	Iterator( Key* begin, Key* end, size_t offset ) : my_ptr(begin+offset), my_begin(begin), my_end(end) {}
	//! Construct undefined iterator
	Iterator() : my_ptr(0), my_begin(0), my_end(poison(-1)) {}
	~Iterator() {
		my_begin = poison(-2);
		my_end = poison(-3);
	}
	Iterator& operator+=( std::ptrdiff_t n ) {
		assert(isDefined());
		assert(my_begin<=my_ptr+n);
		assert(my_ptr+n<=my_end);
		my_ptr += n;
		return *this;
	}
	Iterator& operator++() {
		return operator+=(1);
	}
	Iterator& operator--() {
		return operator+=(-1);
	}
	friend Iterator operator+( const Iterator& i, const std::ptrdiff_t n ) {
		Iterator j(i);
		return j += n;
	}
	friend Iterator operator-( const Iterator& i, const std::ptrdiff_t n ) {
		Iterator j(i);
		return j += -n;
	}
	friend std::ptrdiff_t operator-( const Iterator& i, const Iterator& j ) {
		assert(i.isDefined());
		assert(j.isDefined());
		return i.my_ptr-j.my_ptr;
	}
	friend bool operator==( const Iterator& i, const Iterator& j ) {
		assert(i.isDefined());
		assert(j.isDefined());
		return i.my_ptr==j.my_ptr;
	}
	friend bool operator!=( const Iterator& i, const Iterator& j ) {
		assert(i.isDefined());
		assert(j.isDefined());
		return i.my_ptr!=j.my_ptr;
	}
	friend bool operator<( const Iterator& i, const Iterator& j ) {
		assert(i.isDefined());
		assert(j.isDefined());
		return i.my_ptr<j.my_ptr;
	}
	Key& operator*() {
		assert(isDefined());
		return *my_ptr;
	}
};

namespace std {
	template<>
	class iterator_traits<Iterator> {
	public:
		typedef random_access_iterator_tag iterator_category;
		typedef Key value_type;
		typedef value_type& reference;
		typedef std::ptrdiff_t difference_type;
	};
}

const size_t N_MAX = 100000000;

static Key Array[N_MAX];
static char Flag[N_MAX];

// Initialize Array with n elements.
void CreateDataSet( size_t n ) {
	HighValidIndex = n-1;
	// Keys will be in [0..m-1].  The limit almost ensures that some duplicate keys will occur.
	int m = 2*n;
	for( size_t i=0; i<n; ++i ) {
		Array[i].value = rand() % m;
		Array[i].index = i;
	}
}

// Check that Array is sorted, and sort is stable.
void CheckIsSorted( size_t n ) {

	std::memset(Flag,0,sizeof(Flag));
	for( size_t i=0; i<n; ++i ) {
		int k = Array[i].index;
		if( Flag[k] ) {
			printf("ERROR: duplicate!\n");
			abort();
		}
		Flag[k] = 1;
	}
	if( memchr(Flag,0,n) ) {
		printf("ERROR: missing value!\n");
		abort();
	}
	if( n<2 )
		return;
	for( size_t i=0; i<n-2; ++i ) {
		int ai = Array[i].value;
		int aj = Array[i+1].value;
		int bi = Array[i].index;
		int bj = Array[i+1].index;

		if( ai>aj ) {
			printf("ERROR: not sorted! array[%ld].value = %d > %d = array[%ld].value\n", long(i), ai, aj, long(i+1) );
			abort();
		} else if( ai==aj && bi>bj ) {
			printf("ERROR: not stable! array[%ld].index = %d > %d = array[%ld].index\n", long(i), bi, bj, long(i+1) );
			abort();
		}
	}
}
//----------------------------------------------------------------------
// End of code from test.cpp
//---------------------------------------------------------------------


/*
 * Operand (i.e. the Problem) and Result types are the same. We have to encapsulate all the information needed by the algorithms defined in Intel source code
 */
struct ops{
	Key *Array;
	Iterator start;
	Iterator end;
	Key * temp_buff;
	KeyCompare comp;

	int inplace;	//nel programma originale hanno usato come valore di inplace 2
					//(e usano quel numero come confronto quindi non si puo' mettere un bool)
					//ATTENZIONE: tuttavia tolta la prima chiamata poi lo trattano come bool alternando nella fase di merge
					//l'unico problema e' che in questa maniera non ci ricordiamo del primo inplace
					// e quindi non viene liberata la memoria relativa al buffer temporaneo: lo facciamo a mano una volta eseguito il mergesort
					// per passare il controllo sul KeyCount esistente in test.cpp
};

typedef struct ops Operand;
typedef struct ops Result;

void divide(const Operand &op,std::vector<Operand> &subops)
{

	Iterator xs=op.start;
	Iterator xe=op.end;
	Iterator xm=xs+(xe-xs)/2;

	//quello che era Random Access Iterator 2 e' ora Key*
	Key *zs=op.temp_buff;
	Key* zm=zs+(xm-xs);


	Operand left;
	left.Array=op.Array;
	left.start=xs;
	left.end=xm;
	left.temp_buff=zs;
	left.comp=op.comp;
	left.inplace=!(op.inplace); //come fatto in parallel_stable_sort.h. Di fatto solo la prima chiamata e' in place

	Operand right;
	right.Array=op.Array;
	right.start=xm;
	right.end=xe;
	right.temp_buff=zm;
	right.comp=op.comp;
	right.inplace=!(op.inplace);

	subops.push_back(left);
	subops.push_back(right);
}


/*
 * Base case
 */
void seq(const Operand &op, Result &ret)
{
	pss::internal::stable_sort_base_case(op.start,op.end,op.temp_buff,op.inplace,op.comp);
	ret=op;
	ret.inplace=!op.inplace;	//deve essere invertito (verra' utilizzato a livello superiore)
}


/*
 * Merge function: unlike intel source code we perform a serial merge in this phase
 */
void mergeMS(std::vector<Result>&ress, Result &ret)
{
	if(ress[0].inplace)
	{
		//parallel_move_merge( zs, zm, zm, ze, xs, inplace==2, comp );
		//xs=zs
		Key *xs=ress[0].temp_buff;
		//xe=zm
		Key *xe=ress[1].temp_buff;
		//ys=zm
		Key *ys=xe;
		//ye=ze
		//by definition (in the divide) ze=zs+(xe-xs)
		Key *ye=ress[0].temp_buff+(ress[1].end-ress[0].start);
		//zs=xs
		Iterator zs=ress[0].start;
		//destroy=inplace

		pss::internal::serial_move_merge(xs, xe, ys, ye, zs, ress[0].comp);

		//destroy
		if(ress[0].inplace==2) //dall'algoritmo originale. In realta' qui non entrera' mai perche' l'inplace dell'operando di partenza (=2) si perde
		{
			pss::internal::serial_destroy(xs,xe);
			pss::internal::serial_destroy(ys,ye);
		}

	}
	else
	{
		// parallel_move_merge( xs, xm, xm, xe, zs, false, comp );
		//xs=xs
		Iterator xs=ress[0].start;
		//xe=xm
		Iterator xe=ress[0].end;
		//ys=xm
		Iterator ys=xe;
		//ye=xe (l'originale)
		Iterator ye=ress[1].end;
		//zs=zs(originale)
		Key *zs=ress[0].temp_buff;
		pss::internal::serial_move_merge(xs, xe, ys, ye, zs, ress[0].comp);

	}
	//get the final result
	ret.Array=ress[0].Array;
	ret.start=ress[0].start;
	ret.end=ress[1].end;
	ret.comp=ress[0].comp;
	ret.temp_buff=ress[0].temp_buff;

	//it must be inverted
	ret.inplace=!(ress[0].inplace);
}


/*
 * Base Case condition
 */
bool cond(const Operand &op)
{
	return (op.end-op.start<=CUTOFF);


}
int main(int argc, char *argv[])
{
	//EDIT: ask for the size as command line argument
	if(argc<3)
	{
		fprintf(stderr,"Usage: %s <num_elements> <num_workers>\n",argv[0]);
		exit(-1);
	}

	int n=atoi(argv[1]);
	int nwork=atoi(argv[2]);

	if(n>N_MAX)
	{
		fprintf(stderr,"Max array length: %d\n",N_MAX);
		exit(-1);
	}
	//from test.cpp
	CreateDataSet(n);
	Iterator i(Array,Array+n,0);

	Operand op;
	op.Array=Array;
	op.start=i;
	op.end=i+n;
	op.comp=KeyCompare(KeyCompare::Live);
	op.inplace=2;

	//create the internal buffer (starting from code in parallel_stable sort.h
	typedef typename std::iterator_traits<Iterator>::value_type T;
	pss::internal::raw_buffer z = pss::internal::raw_buffer( sizeof(T)*(op.end-op.start));
	if(z)
		op.temp_buff=(T*)z.get();
	else
	{
		fprintf(stderr,"Error in allocating temp buffer\n");
		exit(-1);
	}



	Result res;
	std::function<void(const Operand&,std::vector<Operand>&)> div(divide);
	std::function <void(const Operand &,Result &)> sq(seq);
	std::function <void(std::vector<Result >&,Result &)> mergef(mergeMS);
	std::function<bool(const Operand &)> cf(cond);

	int count0 = KeyCount;

#if USE_FF
	ff_DC<Operand, Result> dac(div,mergef,sq,cf,op,res,nwork);
#endif
#if USE_OPENMP
	DacOpenmp<Operand, Result> dac(div,mergef,sq,cf,op,res,nwork);
#endif
#if USE_TBB
	DacTBB<Operand, Result> dac(div,mergef,sq,cf,op,res,nwork);
#endif
	//cleanup memory
	pss::internal::serial_destroy(op.temp_buff,op.temp_buff+n);

	long start_t=current_time_usecs();

	//compute
#if USE_FF
	dac.run_and_wait_end();
#else
	dac.compute();
#endif
	long end_t=current_time_usecs();


	//Correctness Checks
	int count1 = KeyCount;
	// Check that number of keys constructed are equal to number destroyed
	if( count0!=count1 ) {
		printf("%d %d\n",count0,count1);
		printf("ERROR: count difference = %d\n",count1-count0);
		//abort();
	}
	// Check that keys were sorted
	CheckIsSorted(n);



	printf("Time (usecs): %Ld\n",end_t-start_t);

}
