#include "util.hpp"

void print_time_stamp() {
    time_t tt = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now()); // get time stamp
    tm loc = *localtime(&tt);                                            // convert time stamp to readable format

    // get elements of time stamp: Y, M, D, h, m, s
    int tms [6] {loc.tm_year + 1900, loc.tm_mon + 1, loc.tm_mday, loc.tm_hour, loc.tm_min, loc.tm_sec};
    // separators for readable time stamp format
    std::string separators [6] {"-", "-", " | ", ":", ":", " | "};

    for (int i=0; i<6; ++i) {
        if (tms[i] < 10) std::cout << "0";     // make sure to use two-digit format: add zero if necessary
        std::cout << tms[i] << separators[i];  // print time organized by separators
    }
}

double get_time() {
    struct timeval tp;
    gettimeofday(&tp, NULL);
    long int ms = tp.tv_sec * 1000 + tp.tv_usec / 1000;
    double t = ms/1000.;
    return t;
}

void get_time(double t0) {
    struct timeval tp;
    gettimeofday(&tp, NULL);
    long int ms = tp.tv_sec * 1000 + tp.tv_usec / 1000;
    double t = ms / 1000.;
    if (MPI_FLAG) {
        if (mpi_world_rank() == 0) {
            std::cout << "time elapsed: ";
            printf("%.3f", t - t0);
            std::cout << "s" << std::endl;
        }
    }
    else {
        std::cout << "time elapsed: ";
        printf("%.3f", t-t0);
        std::cout << "s" << std::endl;
    }
}

void get_time(double t0, std::string prec) {
    struct timeval tp;
    gettimeofday(&tp, NULL);
    if (prec == "us") {
        long int us = tp.tv_sec * 1000000 + tp.tv_usec;
        double t = us / 1000000.;
        if (MPI_FLAG) {
            if (mpi_world_rank() == 0) {
                std::cout << "time elapsed: ";
                printf("%.6f", t - t0);
                std::cout << "s" << std::endl;
            }
        }
        else {
            std::cout << "time elapsed: ";
            printf("%.6f", t - t0);
            std::cout << "s" << std::endl;
        }
    }
    else {
        long int ms = tp.tv_sec * 1000 + tp.tv_usec / 1000;
        double t = ms / 1000.;
        if (MPI_FLAG) {
            if (mpi_world_rank() == 0) {
                std::cout << "time elapsed: ";
                printf("%.3f", t - t0);
                std::cout << "s" << std::endl;
            }
        }
        else {
            std::cout << "time elapsed: ";
            printf("%.3f", t - t0);
            std::cout << "s" << std::endl;
        }
    }
}

void makedir(const std::string& dir_str) {
#ifdef USE_MPI
    if (mpi_world_rank() == 0 or not MPI_FLAG)  // only the process with ID 0 writes into file to avoid collisions
#endif
    {
        const char *dir = dir_str.c_str();
        // Creating Data directory
        if (mkdir(dir, 0777) == -1)
            std::cerr << "Error when creating directory " << dir << " :  " << strerror(errno) << std::endl;

        else
            std::cout << "Directory " << dir << " created \n";
    }
}

void check_input() {
#ifdef ROTATEK2
    assert (nBOS2 == nFER2);
#endif

#if not defined(KELDYSH_FORMALISM)
    static_assert(nFER % 2 == 0, "nFER must be an odd number.");
#endif

    if (BOSONIC_PARAM_FOR_K3) {
        assert(nBOS3 == nFER3); // Frequency grids must be equal in all three dimensions
    }

#if not defined(KELDYSH_FORMALISM) and not defined(ZERO_TEMP)
    assert(nBOS %2 == 1);
    assert(nBOS2%2 == 1);
    assert(nBOS3%2 == 1);
    assert(nFER %2 == 0);
    assert(nFER2%2 == 0);
    assert(nFER3%2 == 0);
#endif
}
