#include "main.h"

// Define global variables declared in header file

int main(void) 
{
  // Read input file
  main_read_input();

  // Initialize domain and flow arrays
  domain_init(); 

  // Allocate arrays
  //alloc_arrays();

  // Create output directory
  create_output();

  // Radon

  // write to file
  write_reconstruct();
   
  // Free and exit
  free_vars();
  return EXIT_SUCCESS;
}
