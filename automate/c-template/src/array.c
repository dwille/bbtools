// ARRAY.C
//
// Contains various tools used in array multiplication comparisons, etc
// props to ASierakowski for merge_sort and merge

#include "array.h"

// entry point
void merge_sort(int *A, int n, int *A2) 
{
  if (n < 2)                          // if there is only one element, done
    return;
  int m = 0.5 * n;                    // cut array in half
  merge_sort(A, m, A2);               // recurse on first half
  merge_sort(A + m, n - m, A2 + m);   // recurse on second half
  merge(A, n, m, A2);                 // merge the two halves
}

void merge(int *A, int n , int m, int *A2) 
{
  // indices
  int i, j, k;
  // temporary arrays
  int *B = malloc(n * sizeof(int));
  int *B2 = malloc(n * sizeof(int));
  i = 0;                      // first-half index
  j = m;                      // second-half index
  // proceed through entire array
  for (k = 0; k < n; k++) {
    if (j == n) {             // if j has reached the end
      B[k] = A[i];            // take the remaining elements of the firsthalfofA
      B2[k] = A2[i];          // take A2 along for the ride
      i++;                    // increment i
    } else if (i == m) {      // else if i has reached half-way
      B[k] = A[j];            // take the remaining of the second half of B 
      B2[k] = A2[j];          // take A2 along for the ride
      j++;                    // increment j
    } else if (A[j] < A[i]) {
      B[k] = A[j];            // else compare two halves of A
      B2[k] = A2[j];          // take second half if smaller
      j++;                    // increment j
    } else {                  // else
      B[k] = A[i];            // take first half if smaller
      B2[k] = A2[i];          // take A2 along for the ride
      i++;                    // increment i
    }
  }
  // overwrite unsorted A with sorted B
  for (i = 0; i < n; i++ ) {
    A[i] = B[i];
    A2[i] = B2[i];            // take A2 along for the ride
  }
  free(B);                    // clean up
  free(B2);                   // clean up
}
