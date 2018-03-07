#ifndef __ERROR_H__
#define __ERROR_H__

#include <stdio.h>
#include <stdlib.h>

#define error(errorMsg) {                                               \
    fprintf(stderr, "Error on line %d of %s: %s\n", __LINE__, __FILE__, errorMsg); \
    exit(EXIT_FAILURE);                                                 \
  }

#define warn(warnMsg) {                                                 \
    fprintf(stderr, "Warning on line %d of %s: %s\n", __LINE__, __FILE__, warnMsg); \
  }

#define log(logMsg) {                                                   \
    printf("Log Message on line %d of %s: %s\n", __LINE__, __FILE__, logMsg); \
  }

#endif
