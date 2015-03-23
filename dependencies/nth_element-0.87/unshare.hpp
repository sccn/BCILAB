/* unshare.hpp
 * Ver 0.86
 * Peter H. Li 2012 FreeBSD License 
 * 
 * This is an undocumented MathWorks internal  method to make editing an 
 * array in-place "safe".  I am not sure I have used it correctly.  I have 
 * only done light testing (R2010A) but it appears to work, i.e. without 
 * the call to unshare, in-place editing a shared variable causes all 
 * copies to be edited and often crashes Matlab, whereas with the call to 
 * unshare everything seems fine.
 */
#ifndef UNSHARE_HPP
#define UNSHARE_HPP

#include "mex.h"
extern "C" int mxUnshareArray(mxArray *array_ptr, int level);

#endif
