#include <stdio.h>
#include <string.h>
#include "output.h"

#define FILE_FORMAT "output/subdomain%d.csv.%d"

void file_output_processor(domain* d, int subdomainIdx, int itrNum, output_selector selector) {
  char fileName[2 * strlen(FILE_FORMAT)];
  sprintf(fileName, FILE_FORMAT, subdomainIdx, itrNum);
  FILE* file = fopen(fileName, "w+");

  subdomain sd = d->subdomains[subdomainIdx];
  int nV = sd.dimX * sd.dimY;
  int i;
  for (i = 0; i < nV; i++) {
    vertex* v = sd.subdomain_vertices[i];
    if (selector(d, subdomainIdx, v)) {
      fprintf(file, "%f, %f, %f\n", v->x, v->y, sd.subdomain_solution.elements[i]);
    }
  }

  fclose(file);
}

void noop_processor(domain* d, int subdomainIdx, int itrNum, output_selector selector) {}
