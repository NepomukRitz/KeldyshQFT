#include "util.hpp"
#include <omp.h>

namespace utils {

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
    #ifdef STATIC_FEEDBACK
        assert(MAX_DIAG_CLASS == 1);
    #endif
        if (MAX_DIAG_CLASS<2) assert(N_LOOPS < 2);
    #ifdef ROTATEK2
        assert (nBOS2 == nFER2);
    #endif

    #if GRID==2
        static_assert(INTERPOLATION!=linear, "Linear interpolation not possible for polar coordinates.");
    #endif

    #if not KELDYSH_FORMALISM
        static_assert(nFER % 2 == 0, "nFER must be an even number.");
    #endif

        if (BOSONIC_PARAM_FOR_K3) {
            assert(nBOS3 == nFER3); // Frequency grids must be equal in all three dimensions
        }

    #if not KELDYSH_FORMALISM and not defined(ZERO_TEMP)
        static_assert(nBOS %2 == 1, "Number of frequency points inconsistent for Matsubara T>0");
        static_assert(nBOS2%2 == 1, "Number of frequency points inconsistent for Matsubara T>0");
        static_assert(nBOS3%2 == 1, "Number of frequency points inconsistent for Matsubara T>0");
        static_assert(nFER %2 == 0, "Number of frequency points inconsistent for Matsubara T>0");
        static_assert(nFER2%2 == 0, "Number of frequency points inconsistent for Matsubara T>0");
        static_assert(nFER3%2 == 0, "Number of frequency points inconsistent for Matsubara T>0");
    #endif

    #if SBE_DECOMPOSITION
    static_assert(SWITCH_SUM_N_INTEGRAL, "SBE requires preprocessing of vertex data.");
    #endif
    }

    std::string generate_data_directory(std::string& job) {
    #if KELDYSH_FORMALISM
    #if DEBUG_SYMMETRIES
        std::string data_directory = "../Data_KF_debug/";
    #else
        std::string data_directory = "../Data_KF" + job + "/";
    #endif
    #else
        #if DEBUG_SYMMETRIES
        std::string data_directory = "../Data_MF_debug/";

    #else
        std::string data_directory = "../Data_MF"+ job +"/";
    #endif
    #endif

    #if SBE_DECOMPOSITION
        data_directory = "../Data_KF_debug_SBE/";
    #endif

        makedir(data_directory);

        return data_directory;
    }

    std::string generate_filename(const fRG_config& config) {
        std::string klass = "K" + std::to_string(MAX_DIAG_CLASS) + "_";
        std::string loops = std::to_string(config.nloops) + "LF_";
        std::string n1 = "n1=" + std::to_string(nBOS) + "_";
        std::string n2 = "n2=" + std::to_string(nBOS2) + "_";
        std::string n3 = "n3=" + std::to_string(nBOS3) + "_";
        std::string gamma = "Gamma=" + std::to_string(glb_Gamma) + "_";
        std::string voltage = "V=" + std::to_string(glb_V) + "_";
        std::string temp = "T=" + std::to_string(glb_T) + "_";
        std::string lambda = "L_ini=" + std::to_string((int)Lambda_ini)+"_";
        std::string ode = "nODE=" + std::to_string(config.nODE_);
        std::string extension = ".h5";

        std::string filename = klass + loops + n1;
        if (MAX_DIAG_CLASS >= 2) filename += n2;
    #if defined(STATIC_FEEDBACK)
        filename += "static_";
    #endif
        if (MAX_DIAG_CLASS >= 3) filename += n3;
        filename += gamma;
        if(glb_V != 0.)
            filename += voltage;
        if(glb_T != 0.01)
            filename += temp;
        filename += lambda + ode + extension;

        return filename;
    }

    void print_job_info() {
        if (KELDYSH){
            if (HUBBARD_MODEL) print("Hubbard model in Keldysh formalism: \n");
            else               print("SIAM in Keldysh formalism: \n");
        }
        else{
            if (HUBBARD_MODEL) print("Hubbard model in Matsubara formalism: \n");
            else               print("SIAM in Matsubara formalism: \n");
        }

        if (PARTICLE_HOLE_SYMMETRY) print("Using PARTICLE HOLE Symmetry\n");

        print("U for this run is: ", glb_U, true);
        print("T for this run is: ", glb_T, true);
        print("Lambda flows from ", Lambda_ini);
        print_add(" to ", Lambda_fin, true);
        print("nODE for this run: ", nODE, true);
        if constexpr (MPI_FLAG) print("MPI World Size = " + std::to_string(mpi_world_size()), true);
    #pragma omp parallel default(none)
        {
    #pragma omp master
            print("OMP Threads = " + std::to_string(omp_get_num_threads()), true);
        }
        print("nBOS1 = ", nBOS, true);
        print("nFER1 = ", nFER, true);
        if(MAX_DIAG_CLASS > 1) {
            print("nBOS2 = ", nBOS2, true);
            print("nFER2 = ", nFER2, true);
        }
        if(MAX_DIAG_CLASS > 2) {
            print("nBOS3 = ", nBOS3, true);
            print("nFER3 = ", nFER3, true);
        }
        if (HUBBARD_MODEL) print("n_in = ", n_in, true);
    }

    void hello_world() {
        print("Hello world \n");
    #ifdef __linux__
        print("on linux.\n");
    #elif __APPLE__
        print("on apple.\n");
        print("on apple.\n");
    #endif
    }

}