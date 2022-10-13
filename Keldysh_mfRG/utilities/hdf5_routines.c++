#include "hdf5_routines.hpp"

// Create the memory data type for storing complex numbers in file


H5::H5File create_hdf_file(const std::string & filename) {
    if (mpi_world_rank() == 0) {
        H5::H5File outfile(filename, H5F_ACC_TRUNC);
        return outfile;
    }
    else {
        H5::H5File file;
        return file;
    }
}
H5::H5File open_hdf_file_readOnly(const std::string & filename) {
    // no problem with reading same file from different nodes
    H5::H5File outfile(filename, H5F_ACC_RDONLY);
    return outfile;
}
H5::H5File open_hdf_file_readWrite(const std::string & filename) {
    if (mpi_world_rank() == 0) {
        H5::H5File outfile(filename, H5F_ACC_RDWR);
        return outfile;
    }
    else {
        H5::H5File file;
        return file;
    }
}

void close_hdf_file(H5::H5File & file) {
    if (mpi_world_rank() == 0) {
        file.close();
    }
}


H5::CompType def_mtype_comp() {
    H5::CompType mtype(sizeof(h5_comp));
    mtype.insertMember(RE, HOFFSET(h5_comp, re), H5::PredType::NATIVE_DOUBLE);
    mtype.insertMember(IM, HOFFSET(h5_comp, im), H5::PredType::NATIVE_DOUBLE);
    return mtype;
}
H5::DSetCreatPropList def_proplist_comp() {
    h5_comp fillvalue_vert;
    fillvalue_vert.re = 0;
    fillvalue_vert.im = 0;
    H5::DSetCreatPropList plist_vert;
    H5::CompType mtype_comp = def_mtype_comp();
    plist_vert.setFillValue(mtype_comp, &fillvalue_vert);
    return plist_vert;
}

hsize_t h5_cast(int dim) {
    return static_cast<hsize_t>(dim);
}

namespace hdf5_impl {
    template <typename gridType>
    void write_freqparams_to_hdf_LambdaLayer(H5::Group& group, const gridType& freqgrid, const int Lambda_it, const int numberLambdaLayers, const bool file_exists, const bool verbose) {
        write_to_hdf_LambdaLayer<char>(group, "type", std::vector<char>({freqgrid.get_type()}), Lambda_it, numberLambdaLayers, file_exists);
        write_to_hdf_LambdaLayer<int>(group, "diag_class", std::vector<int>({freqgrid.get_diag_class()}), Lambda_it, numberLambdaLayers, file_exists);
        write_to_hdf_LambdaLayer<int>(group, "purely_positive", std::vector<int>({freqgrid.purely_positive}), Lambda_it, numberLambdaLayers, file_exists);
        write_to_hdf_LambdaLayer<int>(group, "number_of_gridpoints", std::vector<int>({freqgrid.number_of_gridpoints}), Lambda_it, numberLambdaLayers, file_exists);
        write_to_hdf_LambdaLayer<freqType>(group, "w_upper", std::vector<freqType>({freqgrid.w_upper}), Lambda_it, numberLambdaLayers, file_exists);
        write_to_hdf_LambdaLayer<freqType>(group, "w_lower", std::vector<freqType>({freqgrid.w_lower}), Lambda_it, numberLambdaLayers, file_exists);
        if constexpr(std::is_same_v<gridType,FrequencyGrid<eliasGrid>>)
        {
            write_to_hdf_LambdaLayer<double>(group, "Delta_factor", std::vector<double>({freqgrid.Delta_factor}),
                                             Lambda_it, numberLambdaLayers, file_exists);
            write_to_hdf_LambdaLayer<double>(group, "U_factor", std::vector<double>({freqgrid.U_factor}), Lambda_it,
                                             numberLambdaLayers, file_exists);
            write_to_hdf_LambdaLayer<freqType>(group, "W_scale", std::vector<freqType>({freqgrid.W_scale}), Lambda_it,
                                             numberLambdaLayers, file_exists);

        }
        else if constexpr(std::is_same_v<gridType,FrequencyGrid<hybridGrid>>) {
            write_to_hdf_LambdaLayer<double>(group, "pos_section_boundaries0", std::vector<double>({freqgrid.pos_section_boundaries[0]}),
                                             Lambda_it, numberLambdaLayers, file_exists);
            write_to_hdf_LambdaLayer<double>(group, "pos_section_boundaries1", std::vector<double>({freqgrid.pos_section_boundaries[1]}), Lambda_it,
                                             numberLambdaLayers, file_exists);
        }
        else if constexpr(std::is_same_v<gridType,FrequencyGrid<angularGrid>>) {
            write_to_hdf_LambdaLayer<int>(group, "number_of_intervals", std::vector<int>({freqgrid.number_of_intervals}),
                                             Lambda_it, numberLambdaLayers, file_exists);
        }
        else {
            assert(false);
        }
    }

    template<typename gridType>
    void init_freqgrid_from_hdf_LambdaLayer(H5::Group& group, gridType& freqgrid, const int Lambda_it, const double Lambda) {
        std::vector<char> type;
        std::vector<int> diag_class;
        std::vector<int> purely_positive;
        std::vector<int> number_of_gridpoints;
        std::vector<double> w_upper;
        std::vector<double> w_lower;
        read_from_hdf_LambdaLayer<char>(group, "type", type, Lambda_it);
        read_from_hdf_LambdaLayer<int>(group, "diag_class", diag_class, Lambda_it);
        read_from_hdf_LambdaLayer<int>(group, "purely_positive", purely_positive, Lambda_it);
        read_from_hdf_LambdaLayer<int>(group, "number_of_gridpoints", number_of_gridpoints, Lambda_it);
        read_from_hdf_LambdaLayer<double>(group, "w_upper", w_upper, Lambda_it);
        read_from_hdf_LambdaLayer<double>(group, "w_lower", w_lower, Lambda_it);

        gridType freqgrid_new(type[0], diag_class[0], Lambda, fRG_config(), purely_positive[0]);
        freqgrid_new.w_upper = w_upper[0];
        freqgrid_new.w_lower = w_lower[0];

        if constexpr(std::is_same_v<gridType,FrequencyGrid<eliasGrid>>) {
            std::vector<double> Delta_factor;
            std::vector<double> U_factor;
            std::vector<double> W_scale;
            read_from_hdf_LambdaLayer<double>(group, "Delta_factor", Delta_factor, Lambda_it);
            read_from_hdf_LambdaLayer<double>(group, "U_factor", U_factor, Lambda_it);
            read_from_hdf_LambdaLayer<double>(group, "W_scale", W_scale, Lambda_it);
            freqgrid_new.Delta_factor = Delta_factor[0];
            freqgrid_new.U_factor = U_factor[0];
            freqgrid_new.W_scale = W_scale[0];
        }
        else if constexpr(std::is_same_v<gridType,FrequencyGrid<hybridGrid>>){
            std::vector<double> pos_section_boundaries0;
            std::vector<double> pos_section_boundaries1;
            read_from_hdf_LambdaLayer<double>(group, "pos_section_boundaries0", pos_section_boundaries0, Lambda_it);
            read_from_hdf_LambdaLayer<double>(group, "pos_section_boundaries1", pos_section_boundaries1, Lambda_it);
            freqgrid_new.pos_section_boundaries = std::array<double,2>({pos_section_boundaries0[0], pos_section_boundaries1[0]});
        }
        else if constexpr(std::is_same_v<gridType,FrequencyGrid<angularGrid>>){
            std::vector<int> number_of_intervals;
            read_from_hdf_LambdaLayer<int>(group, "number_of_intervals", number_of_intervals, Lambda_it);
            freqgrid_new.number_of_intervals = number_of_intervals[0];
        }
        else {
            assert(false);
        }


        freqgrid_new.initialize_grid();
        freqgrid = freqgrid_new;
    }

}


State<state_datatype,false> read_state_from_hdf(const H5std_string& filename, const int Lambda_it) {
    H5::H5File file_out = open_hdf_file_readOnly(filename);

    std::vector<double> Lambda;
    read_from_hdf_LambdaLayer<double>(file_out, LAMBDA_LIST, Lambda, Lambda_it);
    fRG_config config = read_config_from_hdf(filename);
    State<state_datatype,false> state(Lambda[0], config);

    std::vector<state_datatype> Sigma_H;
    read_from_hdf_LambdaLayer<state_datatype>(file_out, SELF_LIST, state.selfenergy.Sigma.data, Lambda_it);
    read_from_hdf_LambdaLayer<state_datatype>(file_out, HARTREE, Sigma_H, Lambda_it);
    state.selfenergy.asymp_val_R = Sigma_H[0];
    read_from_hdf_LambdaLayer<state_datatype>(file_out, DATASET_irred, state.vertex.irred().bare, Lambda_it);
    read_from_hdf_LambdaLayer<state_datatype>(file_out, DATASET_K1_a, state.vertex.avertex().K1.data, Lambda_it);
    read_from_hdf_LambdaLayer<state_datatype>(file_out, DATASET_K1_p, state.vertex.pvertex().K1.data, Lambda_it);
    read_from_hdf_LambdaLayer<state_datatype>(file_out, DATASET_K1_t, state.vertex.tvertex().K1.data, Lambda_it);
#if MAX_DIAG_CLASS>1
    read_from_hdf_LambdaLayer<state_datatype>(file_out, DATASET_K2_a, state.vertex.avertex().K2.data, Lambda_it);
    read_from_hdf_LambdaLayer<state_datatype>(file_out, DATASET_K2_p, state.vertex.pvertex().K2.data, Lambda_it);
    read_from_hdf_LambdaLayer<state_datatype>(file_out, DATASET_K2_t, state.vertex.tvertex().K2.data, Lambda_it);
#if DEBUG_SYMMETRIES
    read_from_hdf_LambdaLayer<state_datatype>(file_out, DATASET_K2b_a, state.vertex.avertex().K2b.data, Lambda_it);
    read_from_hdf_LambdaLayer<state_datatype>(file_out, DATASET_K2b_p, state.vertex.pvertex().K2b.data, Lambda_it);
    read_from_hdf_LambdaLayer<state_datatype>(file_out, DATASET_K2b_t, state.vertex.tvertex().K2b.data, Lambda_it);
#endif
#endif
#if MAX_DIAG_CLASS>2
    read_from_hdf_LambdaLayer<state_datatype>(file_out, DATASET_K3_a, state.vertex.avertex().K3.data, Lambda_it);
    read_from_hdf_LambdaLayer<state_datatype>(file_out, DATASET_K3_p, state.vertex.pvertex().K3.data, Lambda_it);
    read_from_hdf_LambdaLayer<state_datatype>(file_out, DATASET_K3_t, state.vertex.tvertex().K3.data, Lambda_it);
#endif
    H5::Group group_freqparams(file_out.openGroup(FREQ_PARAMS));
    H5::Group group_freqparams_ffreqs (group_freqparams.openGroup(FFREQS_LIST ));
    H5::Group group_freqparams_bfreqsa (group_freqparams.openGroup(BFREQS_LISTa ));
    H5::Group group_freqparams_bfreqsp (group_freqparams.openGroup(BFREQS_LISTp ));
    H5::Group group_freqparams_bfreqst (group_freqparams.openGroup(BFREQS_LISTt ));
    H5::Group group_freqparams_bfreqs2a(group_freqparams.openGroup(BFREQS2_LISTa));
    H5::Group group_freqparams_bfreqs2p(group_freqparams.openGroup(BFREQS2_LISTp));
    H5::Group group_freqparams_bfreqs2t(group_freqparams.openGroup(BFREQS2_LISTt));
    H5::Group group_freqparams_bfreqs2ba(group_freqparams.openGroup(BFREQS2b_LISTa));
    H5::Group group_freqparams_bfreqs2bp(group_freqparams.openGroup(BFREQS2b_LISTp));
    H5::Group group_freqparams_bfreqs2bt(group_freqparams.openGroup(BFREQS2b_LISTt));
    H5::Group group_freqparams_bfreqs3a(group_freqparams.openGroup(BFREQS3_LISTa));
    H5::Group group_freqparams_bfreqs3p(group_freqparams.openGroup(BFREQS3_LISTp));
    H5::Group group_freqparams_bfreqs3t(group_freqparams.openGroup(BFREQS3_LISTt));
    H5::Group group_freqparams_ffreqs2a(group_freqparams.openGroup(FFREQS2_LISTa));
    H5::Group group_freqparams_ffreqs2p(group_freqparams.openGroup(FFREQS2_LISTp));
    H5::Group group_freqparams_ffreqs2t(group_freqparams.openGroup(FFREQS2_LISTt));
    H5::Group group_freqparams_ffreqs2ba(group_freqparams.openGroup(FFREQS2b_LISTa));
    H5::Group group_freqparams_ffreqs2bp(group_freqparams.openGroup(FFREQS2b_LISTp));
    H5::Group group_freqparams_ffreqs2bt(group_freqparams.openGroup(FFREQS2b_LISTt));
    H5::Group group_freqparams_ffreqs3a(group_freqparams.openGroup(FFREQS3_LISTa));
    H5::Group group_freqparams_ffreqs3p(group_freqparams.openGroup(FFREQS3_LISTp));
    H5::Group group_freqparams_ffreqs3t(group_freqparams.openGroup(FFREQS3_LISTt));
    H5::Group group_freqparams_ffreqs3a2(group_freqparams.openGroup(FFREQS3_LISTa2));
    H5::Group group_freqparams_ffreqs3p2(group_freqparams.openGroup(FFREQS3_LISTp2));
    H5::Group group_freqparams_ffreqs3t2(group_freqparams.openGroup(FFREQS3_LISTt2));

    hdf5_impl::init_freqgrid_from_hdf_LambdaLayer(group_freqparams_ffreqs , state.selfenergy.Sigma.frequencies.  primary_grid, Lambda_it, Lambda[0]);
    hdf5_impl::init_freqgrid_from_hdf_LambdaLayer(group_freqparams_bfreqsa, state.vertex.avertex().K1.frequencies.  primary_grid, Lambda_it, Lambda[0]);
    hdf5_impl::init_freqgrid_from_hdf_LambdaLayer(group_freqparams_bfreqsp, state.vertex.pvertex().K1.frequencies.  primary_grid, Lambda_it, Lambda[0]);
    hdf5_impl::init_freqgrid_from_hdf_LambdaLayer(group_freqparams_bfreqst, state.vertex.tvertex().K1.frequencies.  primary_grid, Lambda_it, Lambda[0]);
#if MAX_DIAG_CLASS>1
    hdf5_impl::init_freqgrid_from_hdf_LambdaLayer(group_freqparams_bfreqs2a,state.vertex.avertex().K2.frequencies.  primary_grid, Lambda_it, Lambda[0]);
    hdf5_impl::init_freqgrid_from_hdf_LambdaLayer(group_freqparams_bfreqs2p,state.vertex.pvertex().K2.frequencies.  primary_grid, Lambda_it, Lambda[0]);
    hdf5_impl::init_freqgrid_from_hdf_LambdaLayer(group_freqparams_bfreqs2t,state.vertex.tvertex().K2.frequencies.  primary_grid, Lambda_it, Lambda[0]);
    hdf5_impl::init_freqgrid_from_hdf_LambdaLayer(group_freqparams_ffreqs2a,state.vertex.avertex().K2.frequencies.secondary_grid, Lambda_it, Lambda[0]);
    hdf5_impl::init_freqgrid_from_hdf_LambdaLayer(group_freqparams_ffreqs2p,state.vertex.pvertex().K2.frequencies.secondary_grid, Lambda_it, Lambda[0]);
    hdf5_impl::init_freqgrid_from_hdf_LambdaLayer(group_freqparams_ffreqs2t,state.vertex.tvertex().K2.frequencies.secondary_grid, Lambda_it, Lambda[0]);
#if DEBUG_SYMMETRIES
    hdf5_impl::init_freqgrid_from_hdf_LambdaLayer(group_freqparams_bfreqs2ba,state.vertex.avertex().K2b.frequencies.  primary_grid, Lambda_it, Lambda[0]);
    hdf5_impl::init_freqgrid_from_hdf_LambdaLayer(group_freqparams_bfreqs2bp,state.vertex.pvertex().K2b.frequencies.  primary_grid, Lambda_it, Lambda[0]);
    hdf5_impl::init_freqgrid_from_hdf_LambdaLayer(group_freqparams_bfreqs2bt,state.vertex.tvertex().K2b.frequencies.  primary_grid, Lambda_it, Lambda[0]);
    hdf5_impl::init_freqgrid_from_hdf_LambdaLayer(group_freqparams_ffreqs2ba,state.vertex.avertex().K2b.frequencies.secondary_grid, Lambda_it, Lambda[0]);
    hdf5_impl::init_freqgrid_from_hdf_LambdaLayer(group_freqparams_ffreqs2bp,state.vertex.pvertex().K2b.frequencies.secondary_grid, Lambda_it, Lambda[0]);
    hdf5_impl::init_freqgrid_from_hdf_LambdaLayer(group_freqparams_ffreqs2bt,state.vertex.tvertex().K2b.frequencies.secondary_grid, Lambda_it, Lambda[0]);
#endif
#endif
#if MAX_DIAG_CLASS>2
    hdf5_impl::init_freqgrid_from_hdf_LambdaLayer(group_freqparams_bfreqs3a, state.vertex.avertex().K3.frequencies.  primary_grid, Lambda_it, Lambda[0]);
    hdf5_impl::init_freqgrid_from_hdf_LambdaLayer(group_freqparams_bfreqs3p, state.vertex.pvertex().K3.frequencies.  primary_grid, Lambda_it, Lambda[0]);
    hdf5_impl::init_freqgrid_from_hdf_LambdaLayer(group_freqparams_bfreqs3t, state.vertex.tvertex().K3.frequencies.  primary_grid, Lambda_it, Lambda[0]);
    hdf5_impl::init_freqgrid_from_hdf_LambdaLayer(group_freqparams_ffreqs3a, state.vertex.avertex().K3.frequencies.secondary_grid, Lambda_it, Lambda[0]);
    hdf5_impl::init_freqgrid_from_hdf_LambdaLayer(group_freqparams_ffreqs3p, state.vertex.pvertex().K3.frequencies.secondary_grid, Lambda_it, Lambda[0]);
    hdf5_impl::init_freqgrid_from_hdf_LambdaLayer(group_freqparams_ffreqs3t, state.vertex.tvertex().K3.frequencies.secondary_grid, Lambda_it, Lambda[0]);
    hdf5_impl::init_freqgrid_from_hdf_LambdaLayer(group_freqparams_ffreqs3a2,state.vertex.avertex().K3.frequencies. tertiary_grid, Lambda_it, Lambda[0]);
    hdf5_impl::init_freqgrid_from_hdf_LambdaLayer(group_freqparams_ffreqs3p2,state.vertex.pvertex().K3.frequencies. tertiary_grid, Lambda_it, Lambda[0]);
    hdf5_impl::init_freqgrid_from_hdf_LambdaLayer(group_freqparams_ffreqs3t2,state.vertex.tvertex().K3.frequencies. tertiary_grid, Lambda_it, Lambda[0]);
#endif

    /// check whether parameters in HDF file agree with parameters/master_parameters.hpp
        /// Write used parameters for documentation purpose
        int REG_loaded, MAX_DIAG_CLASS_loaded, N_LOOPS_loaded, GRID_loaded;
        double glb_Gamma_loaded, glb_T_loaded, glb_mu_loaded, glb_U_loaded, glb_epsilon_loaded, glb_V_loaded;

    H5::Group group_params(file_out.openGroup(PARAM_LIST));
    read_from_hdf(group_params, "REG", REG_loaded);
    read_from_hdf(group_params, "Gamma", glb_Gamma_loaded);
    read_from_hdf(group_params, "MAX_DIAG_CLASS", MAX_DIAG_CLASS_loaded);
    read_from_hdf(group_params, "N_LOOPS", N_LOOPS_loaded);
    read_from_hdf(group_params, "T", glb_T_loaded);
    read_from_hdf(group_params, "mu", glb_mu_loaded);
    read_from_hdf(group_params, "U", glb_U_loaded);
    read_from_hdf(group_params, "epsilon", glb_epsilon_loaded);
    read_from_hdf(group_params, "V", glb_V_loaded);
    //read_from_hdf(group_params, "ODEsolver", ODEsolver_loaded);
    read_from_hdf(group_params, "GRID", GRID_loaded);
    bool are_parameters_identical = REG == REG_loaded and
            MAX_DIAG_CLASS == MAX_DIAG_CLASS_loaded and
            //glb_T == glb_T_loaded and
            glb_mu == glb_mu_loaded and
            //glb_U == glb_U_loaded and
            //glb_epsilon == glb_epsilon_loaded and
            glb_V == glb_V_loaded and
            GRID == GRID_loaded;
    if (!are_parameters_identical) {
        utils::print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
        utils::print("\t Warning!: \t Parameters of executable do not agree with those in HDF file ", filename, "\n");
        const vec<bool> is_identical_parameter = {
                MAX_DIAG_CLASS == MAX_DIAG_CLASS_loaded
                //,glb_T == glb_T_loaded
                ,glb_mu == glb_mu_loaded
                //,glb_U == glb_U_loaded
                //,glb_epsilon == glb_epsilon_loaded
                ,glb_V == glb_V_loaded
                ,GRID == GRID_loaded
        };
        const vec<std::string> parameter_names = {
                "MAX_DIAG_CLASS"
                //,"glb_T"
                ,"glb_mu"
                //,"glb_U"
                //,"glb_epsilon"
                ,"glb_V"
                ,"GRID"
        };
        for (int i = 0; i < is_identical_parameter.size(); i++){
            if (!is_identical_parameter[i]) utils::print("\t parameter ", parameter_names[i], " is different.\n");
        }
        utils::print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    }

    file_out.close();

    return state;
}


bool check_convergence_hdf(const H5std_string& filename, int& Lambda_it) {
    H5::H5File file_out;
    H5::Exception::dontPrint();
    bool is_converged;
    try {
        file_out = H5::H5File(filename, H5F_ACC_RDONLY);

        H5::Group group_params(file_out.openGroup(PARAM_LIST));
        read_from_hdf(group_params, IS_CONVERGED, is_converged);
        read_from_hdf(group_params, "last_Lambda_it", Lambda_it);

    } catch(const H5::FileIException&) {
        is_converged = false;
        Lambda_it = -1;
    }

    file_out.close();

    return is_converged;
}



rvec read_Lambdas_from_hdf(const H5std_string FILE_NAME){



    // Open the file. Access rights: read-only
    H5::H5File file = open_hdf_file_readOnly(FILE_NAME);

    H5::DataSet lambda_dataset = file.openDataSet("lambdas");


    // Create the data spaces for the data sets in file and for buffer objects
    H5::DataSpace dataSpace_Lambda = lambda_dataset.getSpace();
    //
    //Get the dimension size of each dimension in the dataspace and
    //display them.
    //
    hsize_t dims_out[2];
    dataSpace_Lambda.getSimpleExtentDims( dims_out, NULL);
    //std::cout << "rank " << rank << ", dimensions " <<
    //     (unsigned long)(dims_out[0]) << " x " <<
    //     (unsigned long)(dims_out[1]) << std::endl;


    H5::DataSpace dataSpace_Lambda_buffer(1, dims_out);

    /// load data into buffer
    double Lambdas_arr[dims_out[0]];
    lambda_dataset.read(Lambdas_arr, H5::PredType::NATIVE_DOUBLE,
                        dataSpace_Lambda_buffer, dataSpace_Lambda);

    rvec Lambdas (Lambdas_arr, Lambdas_arr + sizeof(Lambdas_arr) / sizeof(double));

    // Terminate
    lambda_dataset.close();
    close_hdf_file(file);

    return Lambdas;


}


bool test_read_write_data_hdf(bool verbose) {
    //if (verbose) utils::print("Testing HDF input and output:", true);

    vec<comp> vector_original(10, glb_i);
    std::array<std::size_t,3> dims_mularr = {1,2,3};
    multidimensional::multiarray<comp,3> mularr_original(dims_mularr, glb_i);

    std::string testfile_name = "test.h5";
    //if (verbose) utils::print("File name: " + testfile_name, true);

    // write to hdf file
    H5::H5File file_out = create_hdf_file(testfile_name.c_str());
    //if (verbose) utils::print("Created file.", true);
    //H5::H5std_string FILE_NAME(data_dir + "test.h5");
    H5::Group group_out(file_out.createGroup("/data"));
    //if (verbose) utils::print("Created group.", true);

    write_to_hdf(group_out, "vector", vector_original, false);
    write_to_hdf(group_out, "mularr", mularr_original, false);
    //if (verbose) utils::print("Written to hdf.", true);

    // read from hdf file
    hid_t file_id = H5Fopen(testfile_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    H5::Group group = H5Gopen(file_id, "data", H5P_DEFAULT);
    //if (verbose) utils::print("Loaded group.", true);
    vec<comp> vector_input;
    read_from_hdf<comp>(group, "vector", vector_input);
    //if (verbose) utils::print("Loaded vector.", true);
    multidimensional::multiarray<comp,3> mularr_input;
    read_from_hdf<comp,3>(group, "mularr", mularr_input);
    //if (verbose) utils::print("Loaded multiarray.", true);

    //if (verbose) utils::print("Length of original vec: ", vector_original.size(), true);
    //if (verbose) utils::print("Length of loaded vec: ", vector_input.size(), true);

    comp deviation_vec = (vector_input - vector_original).max_norm();
    if (std::abs(deviation_vec) < 1e-10) utils::print("Read / write of vector to HDF file successful.", true);
    else utils::print("PROBLEM during read / write of vector to HDF file.", true);

    comp deviation_mularr=(mularr_input -mularr_original).maxabs();
    if (std::abs(deviation_mularr) < 1e-10) utils::print("Read / write of multiarray to HDF file successful.", true);
    else utils::print("PROBLEM during read / write of multiarray to HDF file.", true);



    //H5::H5std_string FILE_NAME(data_dir + "test.h5");
    H5::Group group_out_L(file_out.createGroup("/data_with_LambdaLayers"));
    //if (verbose) putils::print("Created group.", true);

    write_to_hdf_LambdaLayer(group_out_L, "vector", vector_original, 1, 10,  false);
    write_to_hdf_LambdaLayer(group_out_L, "mularr", mularr_original, 1, 10,  false);
    //if (verbose) putils::print("Written to hdf LambdaLayer.", true);

    // read from hdf file
    //hid_t file_id = H5Fopen(testfile_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    H5::Group group_L = H5Gopen(file_id, "data_with_LambdaLayers", H5P_DEFAULT);
    vec<comp> vector_input_L;
    read_from_hdf_LambdaLayer<comp>(group_L, "vector", vector_input_L, 1);
    multidimensional::multiarray<comp,3> mularr_input_L;
    read_from_hdf_LambdaLayer<comp>(group_L, "mularr", mularr_input_L, 1);


    comp deviation_vec_L = (vector_input_L - vector_original).max_norm();
    comp deviation_mularr_L=(mularr_input_L -mularr_original).maxabs();


    bool passed = true;
    if (verbose) {
        if (std::abs(deviation_vec_L) < 1e-10)
            utils::print("Read / write of vector to LambdaLayer of HDF file successful.", true);
        else {
            passed = false;
            utils::print("PROBLEM during read / write of vector to LambdaLayer of HDF file.", true);
        }
        if (std::abs(deviation_mularr_L) < 1e-10)
            utils::print("Read / write of multiarray to LambdaLayer ofHDF file successful.", true);
        else {
            passed = false;
            utils::print("PROBLEM during read / write to LambdaLayer of multiarray to HDF file.", true);
        }
    }

    return passed;
}


bool test_read_write_state_hdf(bool verbose) {
    //if (verbose) utils::print("Testing HDF input and output of states:", true);
    const double Lambda = 1.;
    const int Lambda_it = 1;
    const int numberLambdaLayers = 10;
    fRG_config config;
    State<state_datatype> state_output(Lambda, config);
    state_output.initialize();
    state_output = state_output + 1.;

    write_state_to_hdf("test_state.h5", 0. , numberLambdaLayers, state_output);
    //utils::print("Written state", true);


    State<state_datatype> state_input = read_state_from_hdf("test_state.h5", 0);
    //utils::print("Read state", true);

    State<state_datatype> state_diff = state_output - state_input;

    add_state_to_hdf("test_state.h5", Lambda_it  , state_output*2);
    add_state_to_hdf("test_state.h5", Lambda_it+1, state_input);

    write_state_to_hdf("test_statediff.h5", 0., numberLambdaLayers, state_diff);

    bool passed = true;

    if (verbose) {
        if (state_diff.norm() < 1e-10) utils::print("Read / write of state data to LambdaLayer of HDF file successful.", true);
        else {
            passed = false;
            utils::print("PROBLEM during read / write of state data to LambdaLayer of HDF file. deviation : ",
                  state_diff.norm(), true);
        }
        if ((state_output.selfenergy.Sigma.frequencies.  primary_grid.get_all_frequencies() -
             state_input.selfenergy.Sigma.frequencies.  primary_grid.get_all_frequencies()).max_norm() < 1e-10)
            utils::print("Read / write of frequency grid to LambdaLayer of HDF file successful.", true);
        else {
            passed = false;
            utils::print("PROBLEM during read / write of frequency grid to LambdaLayer to HDF file.", true);
        }
    }

    return passed;

}


fRG_config read_config_from_hdf(const H5std_string FILE_NAME) {
    H5::H5File file_out(FILE_NAME, H5F_ACC_RDONLY);

    fRG_config config;
    int REG_loaded, MAX_DIAG_CLASS_loaded, N_LOOPS_loaded, GRID_loaded;
    double glb_Gamma_loaded, glb_T_loaded, glb_mu_loaded, glb_U_loaded, glb_epsilon_loaded, glb_V_loaded;

    H5::Group group_params(file_out.openGroup(PARAM_LIST));
    read_from_hdf(group_params, "REG", REG_loaded);
    read_from_hdf(group_params, "Gamma", glb_Gamma_loaded);
    read_from_hdf(group_params, "MAX_DIAG_CLASS", MAX_DIAG_CLASS_loaded);
    read_from_hdf(group_params, "N_LOOPS", N_LOOPS_loaded);
    read_from_hdf(group_params, "T", glb_T_loaded);
    read_from_hdf(group_params, "mu", glb_mu_loaded);
    read_from_hdf(group_params, "U", glb_U_loaded);
    read_from_hdf(group_params, "epsilon", glb_epsilon_loaded);
    read_from_hdf(group_params, "V", glb_V_loaded);

    // write values in fRG_config
    config.nloops = N_LOOPS_loaded;
    config.U = glb_U_loaded;
    config.Gamma = glb_Gamma_loaded;
    config.T = glb_T_loaded;
    config.epsilon = glb_epsilon_loaded;

    return config;

}