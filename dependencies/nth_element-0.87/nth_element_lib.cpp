/* nth_element_lib.cpp
 * Ver 0.86
 * Peter H. Li 2012 FreeBSD License 
 * See nth_element.m for documentation. 
 */
#ifndef NTH_ELEMENT_LIB_CPP
#define NTH_ELEMENT_LIB_CPP

#include <algorithm>


// This runs on data inplace!
template <typename T> 
void nth_element_cols(T *data, mwIndex rank, mwSize ncols, mwSize nrows) {
  #pragma omp parallel for
  for (mwIndex i = 0; i < ncols; ++i) {
    // Figure out linear indices into this column
    mwIndex start    = i * nrows;
    mwIndex thisRank = start + rank;
    mwIndex end      = start + nrows;

    // Run nth_element to iteratively partition to the specified rank
    std::nth_element(data + start, data + thisRank, data + end);
  }
}


// Determine type of data, run
// This runs on data inplace!
void run_nth_element(mxArray *inarr, mwIndex rank, mwSize ncols, mwSize nrows) {
  void *indata = mxGetData(inarr);

  switch (mxGetClassID(inarr)) {
    case mxDOUBLE_CLASS:
      nth_element_cols((double *) indata, rank, ncols, nrows);
      break;

    case mxSINGLE_CLASS:
      nth_element_cols((float *) indata, rank, ncols, nrows);
      break;

    case mxINT8_CLASS:
      nth_element_cols((signed char *) indata, rank, ncols, nrows);
      break;

    case mxUINT8_CLASS:
      nth_element_cols((unsigned char *) indata, rank, ncols, nrows);
      break;

    case mxINT16_CLASS:
      nth_element_cols((signed short *) indata, rank, ncols, nrows);
      break;

    case mxUINT16_CLASS:
      nth_element_cols((unsigned short *) indata, rank, ncols, nrows);
      break;

    case mxINT32_CLASS:
      nth_element_cols((signed int *) indata, rank, ncols, nrows);
      break;

    case mxUINT32_CLASS:
      nth_element_cols((unsigned int *) indata, rank, ncols, nrows);
      break;

    // Uncomment these if int64 is needed, but note that on some compilers
    // it's called "__int64" instead of "long long"
    //case mxINT64_CLASS:
      //nth_element_cols((signed long long *) indata, rank, ncols, nrows);
      //break;

    //case mxUINT64_CLASS:
      //nth_element_cols((unsigned long long *) indata, rank, ncols, nrows);
      //break;

    default:
      mexErrMsgIdAndTxt("Numerical:nth_element:prhs", "Unrecognized numeric array type.");
  }
}


template <class T> class IndexingComparatorFunct
{
private:
  T*  _data;
  mwIndex _start;
    
public:
  IndexingComparatorFunct(T* data, mwIndex start) : _data(data), _start(start) {}
  
  bool operator()(unsigned int i, unsigned int j) {
      return _data[_start + i] < _data[_start + j];
  }
};


// This runs on data inplace!
template <typename T> void nth_element_cols(T *data, unsigned int *indices, mwIndex rank, mwSize ncols, mwSize nrows) {
  #pragma omp parallel for
  for (mwIndex i = 0; i < ncols; ++i) {
    T *temp = new T[nrows];
    mwIndex start    = i * nrows;
    mwIndex thisRank = start + rank;
    mwIndex end      = start + nrows;
    
    // Create meta-indexed comparator
    IndexingComparatorFunct<T> f(data, start);
    
    // Run nth_element to iteratively partition meta-indices to the specified rank
    std::nth_element(indices + start, indices + thisRank, indices + end, f);
    
    // Now use partitioned meta-indices to arrange actual data
    for (mwIndex j = 0; j < nrows; ++j) temp[j] = data[start+j];
    for (mwIndex j = 0; j < nrows; ++j) data[start+j] = temp[indices[start+j]];
    delete[] temp;
  }
}


// Determine type of data, run
// This runs on data inplace!
void run_nth_element(mxArray *inarr, unsigned int *indices, mwIndex rank, mwSize ncols, mwSize nrows) {
  void *indata = mxGetData(inarr);

  switch (mxGetClassID(inarr)) {
    case mxDOUBLE_CLASS:
      nth_element_cols((double *) indata, indices, rank, ncols, nrows);
      break;

    case mxSINGLE_CLASS:
      nth_element_cols((float *) indata, indices, rank, ncols, nrows);
      break;

    case mxINT8_CLASS:
      nth_element_cols((signed char *) indata, indices, rank, ncols, nrows);
      break;

    case mxUINT8_CLASS:
      nth_element_cols((unsigned char *) indata, indices, rank, ncols, nrows);
      break;

    case mxINT16_CLASS:
      nth_element_cols((signed short *) indata, indices, rank, ncols, nrows);
      break;

    case mxUINT16_CLASS:
      nth_element_cols((unsigned short *) indata, indices, rank, ncols, nrows);
      break;

    case mxINT32_CLASS:
      nth_element_cols((signed int *) indata, indices, rank, ncols, nrows);
      break;

    case mxUINT32_CLASS:
      nth_element_cols((unsigned int *) indata, indices, rank, ncols, nrows);
      break;

    // Uncomment these if int64 is needed, but note that on some compilers
    // it's called "__int64" instead of "long long"
    //case mxINT64_CLASS:
      //nth_element_cols((signed long long *) indata, indices, rank, ncols, nrows);
      //break;

    //case mxUINT64_CLASS:
      //nth_element_cols((unsigned long long *) indata, indices, rank, ncols, nrows);
      //break;

    default:
      mexErrMsgIdAndTxt("Numerical:nth_element:prhs", "Unrecognized numeric array type.");
  }
}

#endif
