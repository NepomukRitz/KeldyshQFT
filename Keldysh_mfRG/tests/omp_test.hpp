#ifndef KELDYSH_MFRG_OMP_TEST_H
#define KELDYSH_MFRG_OMP_TEST_H

#include <omp.h>    // omp functionality
#include <iostream> // text input/output

/**
 * Test function to check omp functionality.
 * Prints one line for each active thread, and total number of threads.
 */
void test_omp() {
    int nthreads, tid; // number of threads, thread ID

    // Fork a team of threads giving them their own copies of variables
#pragma omp parallel private(nthreads, tid)
    {

        // Obtain thread number
        tid = omp_get_thread_num();
        printf("Hello World from thread = %d\n", tid);

        // Only master thread does this
        if (tid == 0)
        {
            nthreads = omp_get_num_threads();
            printf("Number of threads = %d\n", nthreads);
        }

    }  // All threads join master thread and disband
}


#endif //KELDYSH_MFRG_OMP_TEST_H
