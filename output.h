#ifndef __OUTPUT_H__
#define __OUTPUT_H__

#include "domain.h"

typedef int (*output_selector)(domain* d, int subdomainIdx, vertex* v);

typedef void (*output_processor)(domain* d, int subdomainIdx, int itrNum, output_selector selector);

void file_output_processor(domain* d, int subdomainIdx, int itrNum, output_selector selector);

void noop_processor(domain* d, int subdomainIdx, int itrNum, output_selector selector);

#endif
