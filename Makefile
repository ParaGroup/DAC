#Please specify the path of the following libraries
#the following environment variable that contains the needed library
#FASTFLOW_DIR			= path to the fasflow library
#INTEL_STABLESORT_DIR	= path to the intel stable sort directory. It can be found at https://software.intel.com/en-us/articles/a-parallel-stable-sort-using-c11-for-tbb-cilk-plus-and-openmp	

CXX				= g++
CXXFLAGS		= -O3 --std=c++11
LIBS			= -lpthread -lm -lrt
SRC				= src
INCLUDES		= includes
EXE				= fibonacci_dac_ff fibonacci_dac_openmp fibonacci_dac_tbb mergesort_dac_ff mergesort_dac_openmp\
					mergesort_dac_tbb quicksort_dac_ff quicksort_dac_openmp quicksort_dac_tbb strassen_dac_ff\
					strassen_dac_openmp strassen_dac_tbb stable_mergesort_dac_ff stable_mergesort_dac_openmp\
					stable_mergesort_dac_tbb strassen_hm_omp strassen_hm_tbb intel_sort_tbb intel_sort_openmp\
					quicksort_hm_openmp quicksort_hm_tbb
FF_FLAGS		= -I$(FASTFLOW_DIR) -DUSE_FF -DDONT_USE_FFALLOC
OMP_FLAGS		= -fopenmp -DUSE_OPENMP
TBB_FLAGS		= -ltbb -DUSE_TBB

.PHONY: clean

all: $(EXE)

fibonacci_dac_ff: $(SRC)/fibonacci_dac.cpp utils.o
	$(CXX) $(CXXFLAGS)  -o $@ $^ $(LIBS) $(FF_FLAGS)

fibonacci_dac_openmp: $(SRC)/fibonacci_dac.cpp utils.o
	$(CXX) $(CXXFLAGS)  -o $@ $^ $(LIBS) $(OMP_FLAGS)

fibonacci_dac_tbb: $(SRC)/fibonacci_dac.cpp utils.o
	$(CXX) $(CXXFLAGS)  -o $@ $^ $(LIBS) $(TBB_FLAGS)

mergesort_dac_ff: $(SRC)/mergesort_dac.cpp utils.o
	$(CXX) $(CXXFLAGS)  -o $@ $^ $(LIBS) $(FF_FLAGS)

mergesort_dac_openmp: $(SRC)/mergesort_dac.cpp utils.o
	$(CXX) $(CXXFLAGS)  -o $@ $^ $(LIBS) $(OMP_FLAGS)

mergesort_dac_tbb: $(SRC)/mergesort_dac.cpp utils.o
	$(CXX) $(CXXFLAGS)  -o $@ $^ $(LIBS) $(TBB_FLAGS)

quicksort_dac_ff: $(SRC)/quicksort_dac.cpp utils.o
	$(CXX) $(CXXFLAGS)  -o $@ $^ $(LIBS) $(FF_FLAGS)

quicksort_dac_openmp: $(SRC)/quicksort_dac.cpp utils.o
	$(CXX) $(CXXFLAGS)  -o $@ $^ $(LIBS) $(OMP_FLAGS)

quicksort_dac_tbb: $(SRC)/quicksort_dac.cpp utils.o
	$(CXX) $(CXXFLAGS)  -o $@ $^ $(LIBS) $(TBB_FLAGS)

quicksort_hm_openmp: $(SRC)/quicksort_hm_openmp.cpp utils.o
	$(CXX) $(CXXFLAGS)  -o $@ $^ $(LIBS) $(OMP_FLAGS)

quicksort_hm_tbb: $(SRC)/quicksort_hm_tbb.cpp utils.o
	$(CXX) $(CXXFLAGS)  -o $@ $^ $(LIBS) $(TBB_FLAGS)

strassen_dac_ff: $(SRC)/strassen_dac.cpp utils.o
	$(CXX) $(CXXFLAGS)  -o $@ $^ $(LIBS) $(FF_FLAGS)

strassen_dac_openmp: $(SRC)/strassen_dac.cpp utils.o
	$(CXX) $(CXXFLAGS)  -o $@ $^ $(LIBS) $(OMP_FLAGS)

strassen_dac_tbb: $(SRC)/strassen_dac.cpp utils.o
	$(CXX) $(CXXFLAGS)  -o $@ $^ $(LIBS) $(TBB_FLAGS)

stable_mergesort_dac_ff: $(SRC)/stable_mergesort_dac.cpp utils.o
	$(CXX) $(CXXFLAGS)  -o $@ $^ $(LIBS) $(FF_FLAGS) -I$(INTEL_STABLESORT_DIR)

stable_mergesort_dac_openmp: $(SRC)/stable_mergesort_dac.cpp utils.o
	$(CXX) $(CXXFLAGS)  -o $@ $^ $(LIBS) $(OMP_FLAGS) -I$(INTEL_STABLESORT_DIR)

stable_mergesort_dac_tbb: $(SRC)/stable_mergesort_dac.cpp utils.o
	$(CXX) $(CXXFLAGS)  -o $@ $^ $(LIBS) $(TBB_FLAGS) -I$(INTEL_STABLESORT_DIR)

strassen_hm_omp: src/strassen_hm_omp.cpp
	$(CXX) $(CXXFLAGS)  -o $@ $^ $(LIBS) -fopenmp

strassen_hm_tbb: src/strassen_hm_tbb.cpp
	$(CXX) $(CXXFLAGS)  -o $@ $^ $(LIBS) -ltbb

intel_sort_tbb: $(INTEL_STABLESORT_DIR)/test.cpp
	$(CXX) $(CXXFLAGS)  -o $@ $^ $(LIBS) -ltbb -DUSE_TBB_LOWLEVEL -I$(INTEL_STABLESORT_DIR)

intel_sort_openmp: $(INTEL_STABLESORT_DIR)/test.cpp
	$(CXX) $(CXXFLAGS)  -o $@ $^ $(LIBS) $(OMP_FLAGS) -I$(INTEL_STABLESORT_DIR)

%.o: $(SRC)/%.cpp $(INCLUDES)/*
	$(CXX) $(CXXFLAGS) -c -o $@ $<  $(LIBS)

clean:
	rm -f *.o $(EXE)
