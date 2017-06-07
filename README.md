# About DAC 
This repository contains the source code used to perform the experimental evaluations in the paper entitled *A Divide-and-Conquer Parallel Pattern Implementation for Multicores*, presented at *SEPS 2016*.


In the repository you can find the applications used for the evaluations and the backend implementation of the Parallel Divide and Conquer pattern in OpenMP, Intel TBB and Fastflow.

The pattern (and related backend implementations) can be used to easily parallelize other Divide and Conquer algorithms. Details on the interface can be found in the paper.

##Applications
To understand how the pattern works and its interface a basic example for the *n-th* fibonacci number computation is provided.

The main applications used for the evaluation are essentially three: the merge- and quick-sort algorithms
and the Strassen algorithm for matrix multiplication. It is important to notice that for the three applications the main program is the same for the different backends (can be found under the `src/` folder). The different backends can be selected by using proper compiler directives (`USE_OPENMP` for OpenMP`USE_TBB` for the Intel TBB version and `USE_FF` for the FastFlow version).

In addition, to compare the pattern based version with third-party algorithms are present hand-made parallelizations of the aforementioned applications (for the merge-sort comparison we used the stable sort implementation provided by Intel [here](https://software.intel.com/en-us/articles/a-parallel-stable-sort-using-c11-for-tbb-cilk-plus-and-openmp)).


## Usage

### Required software
The artifact uses external software. In particular:

* `FastFlow`: a C++ parallel programming framework targeting shared-memory architectures. Website: http://calvados.di.unipi.it/
*  `Intel Stable Sort`:  a C++11 implementation of a Stable Merge Sort provided by Intel. Website: https://software.intel.com/en-us/articles/a-parallel-stable-sort-using-c11-for-tbb-cilk-plus-and-openmp

In the sequel we will provide a brief description for their downloading

#### Fastflow
It is an header-only library. Therefore, it is only required to download it from the website or the SVN. To download the latest version and save it into the fastflow directory, run the following
command in the shell:

    $ svn checkout svn://svn.code.sf.net/p/mc-fastflow/code/ fastflow


#### Intel Stable Sort
The source code can be downloaded at https://software.intel.com/en-us/articles/a-parallel-stable-sort-using-c11-for-tbb-cilk-plus-and-openmp. Once decompressed and compiled it is ready to be used.


### Compilation
Before compiling the user must provide two diffent environment variables: 
`FASTFLOW_DIR` that points to the FastFlow library directory and `INTEL_STABLESORT_DIR` that points to the directory that contains the Intel source code (used for the comparison), After that, the code can be compiled. The set of command is the following:

     $ export FASTFLOW_DIR=<...path to fastflow...>
     $ export INTEL_STABLESORT_DIR=<... path to intel source code...>
     $ make -j

This will produce different executables:

 - `fibonacci_dac_{openmp,tbb,ff}`: are the the parallel pattern based implementations of the fibonacci  problem that use the OpenMP, Intel TBB and FastFlow backends respectively;
 - `mergesort_dac_{openmp,tbb,ff}`: that are the parallel pattern based implementations of the mergesort problem;
 - `quicksort_dac_{openmp,tbb,ff}`: the  implementations for the quicksort problems for the different backends;
 - `strassen_dac_{openmp,tbb,ff}`: implementations for the Strassen matrices multiplication algorithm;
 - `stable_mergesort_dac_{openmp,tbb,ff}`: implementation of the Intel Stable Sort algorithm used for the comparison. It is essentially the same algorithm (with the same classes and data types) provided by Intel whose divide-and-conquer part is parallelized using the proposed pattern;
 -  `quicksort_hm_{openmp,tbb}` and `strassen_hm_{openmp,tbb}`: hand made parallelizations for OpenMP and TBB
 -  `intel_sort_{openmp,tbb}`: the intel version of the program. Can be compiled directly from the source codes provided in the Intel WebSite.

Each of these programs require certain parameters. To see the right sequence it is sufficient to invoke the program without arguments.


For any further information please write at: dematteis <at> di.unipi.it

