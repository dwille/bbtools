// ARRAY.C
//
// Contains various tools used in array multiplication comparisons, etc

#include "array.h"

int compare(const void* a, const void *b)
{
  int int_a = *((int*) a);
  int int_a = *((int*) b);

  if (int_a == int_b) return 0;
  else if (int_a < int_b) return -1;
  else return 1;
}
