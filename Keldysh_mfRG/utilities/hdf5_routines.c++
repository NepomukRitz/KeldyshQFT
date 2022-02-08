#include "hdf5_routines.hpp"

// Create the memory data type for storing complex numbers in file


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
    void write_freqparams_to_hdf_LambdaLayer(H5::Group& group, const FrequencyGrid& freqgrid, const int Lambda_it, const int numberLambdaLayers, const bool file_exists, const bool verbose) {
        write_to_hdf_LambdaLayer<char>(group, "type", std::vector<char>({freqgrid.get_type()}), Lambda_it, numberLambdaLayers, file_exists);
        write_to_hdf_LambdaLayer<int>(group, "diag_class", std::vector<int>({freqgrid.get_diag_class()}), Lambda_it, numberLambdaLayers, file_exists);
        write_to_hdf_LambdaLayer<int>(group, "N_w", std::vector<int>({freqgrid.N_w}), Lambda_it, numberLambdaLayers, file_exists);
        write_to_hdf_LambdaLayer<double>(group, "w_upper", std::vector<double>({freqgrid.w_upper}), Lambda_it, numberLambdaLayers, file_exists);
        write_to_hdf_LambdaLayer<double>(group, "w_lower", std::vector<double>({freqgrid.w_lower}), Lambda_it, numberLambdaLayers, file_exists);
        write_to_hdf_LambdaLayer<double>(group, "Delta_factor", std::vector<double>({freqgrid.Delta_factor}), Lambda_it, numberLambdaLayers, file_exists);
        write_to_hdf_LambdaLayer<double>(group, "U_factor", std::vector<double>({freqgrid.U_factor}), Lambda_it, numberLambdaLayers, file_exists);
        write_to_hdf_LambdaLayer<double>(group, "W_scale", std::vector<double>({freqgrid.W_scale}), Lambda_it, numberLambdaLayers, file_exists);
    }

    void init_freqgrid_from_hdf_LambdaLayer(H5::Group& group, FrequencyGrid& freqgrid, const int Lambda_it, const double Lambda) {
        std::vector<char> type;
        std::vector<int> diag_class;
        std::vector<int> N_w;
        std::vector<double> w_upper;
        std::vector<double> w_lower;
        std::vector<double> Delta_factor;
        std::vector<double> U_factor;
        std::vector<double> W_scale;
        read_from_hdf_LambdaLayer<char>(group, "type", type, Lambda_it);
        read_from_hdf_LambdaLayer<int>(group, "diag_class", diag_class, Lambda_it);
        read_from_hdf_LambdaLayer<int>(group, "N_w", N_w, Lambda_it);
        read_from_hdf_LambdaLayer<double>(group, "w_upper", w_upper, Lambda_it);
        read_from_hdf_LambdaLayer<double>(group, "w_lower", w_lower, Lambda_it);
        read_from_hdf_LambdaLayer<double>(group, "Delta_factor", Delta_factor, Lambda_it);
        read_from_hdf_LambdaLayer<double>(group, "U_factor", U_factor, Lambda_it);
        read_from_hdf_LambdaLayer<double>(group, "W_scale", W_scale, Lambda_it);

        FrequencyGrid freqgrid_new(type[0], diag_class[0], Lambda);
        freqgrid_new.w_upper = w_upper[0];
        freqgrid_new.w_lower = w_lower[0];
        freqgrid_new.Delta_factor = Delta_factor[0];
        freqgrid_new.U_factor = U_factor[0];
        freqgrid_new.W_scale = W_scale[0];
        freqgrid_new.initialize_grid();
        freqgrid = freqgrid_new;
    }

}


State<state_datatype> read_state_from_hdf_LambdaLayer(const H5std_string& filename, const int Lambda_it) {
    H5::H5File file_out(filename, H5F_ACC_RDONLY);

    std::vector<double> Lambda;
    read_from_hdf_LambdaLayer<state_datatype>(file_out, LAMBDA_LIST, Lambda, Lambda_it);
    State<state_datatype> state(Lambda[0]);

    std::vector<state_datatype> Sigma_H;
    read_from_hdf_LambdaLayer<state_datatype>(file_out, SELF_LIST, state.selfenergy.Sigma, Lambda_it);
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
#endif
#if MAX_DIAG_CLASS>2
    read_from_hdf_LambdaLayer<state_datatype>(file_out, DATASET_K3_a, state.vertex.avertex().K3.data, Lambda_it);
    read_from_hdf_LambdaLayer<state_datatype>(file_out, DATASET_K3_p, state.vertex.pvertex().K3.data, Lambda_it);
    read_from_hdf_LambdaLayer<state_datatype>(file_out, DATASET_K3_t, state.vertex.tvertex().K3.data, Lambda_it);
#endif
    H5::Group group_freqparams(file_out.openGroup(FREQ_PARAMS));
    H5::Group group_freqparams_bfreqsa (group_freqparams.openGroup(BFREQS_LISTa ));
    H5::Group group_freqparams_bfreqsp (group_freqparams.openGroup(BFREQS_LISTp ));
    H5::Group group_freqparams_bfreqst (group_freqparams.openGroup(BFREQS_LISTt ));
    H5::Group group_freqparams_bfreqs2a(group_freqparams.openGroup(BFREQS2_LISTa));
    H5::Group group_freqparams_bfreqs2p(group_freqparams.openGroup(BFREQS2_LISTp));
    H5::Group group_freqparams_bfreqs2t(group_freqparams.openGroup(BFREQS2_LISTt));
    H5::Group group_freqparams_bfreqs3a(group_freqparams.openGroup(BFREQS3_LISTa));
    H5::Group group_freqparams_bfreqs3p(group_freqparams.openGroup(BFREQS3_LISTp));
    H5::Group group_freqparams_bfreqs3t(group_freqparams.openGroup(BFREQS3_LISTt));
    H5::Group group_freqparams_ffreqs2a(group_freqparams.openGroup(FFREQS2_LISTa));
    H5::Group group_freqparams_ffreqs2p(group_freqparams.openGroup(FFREQS2_LISTp));
    H5::Group group_freqparams_ffreqs2t(group_freqparams.openGroup(FFREQS2_LISTt));
    H5::Group group_freqparams_ffreqs3a(group_freqparams.openGroup(FFREQS3_LISTa));
    H5::Group group_freqparams_ffreqs3p(group_freqparams.openGroup(FFREQS3_LISTp));
    H5::Group group_freqparams_ffreqs3t(group_freqparams.openGroup(FFREQS3_LISTt));

    hdf5_impl::init_freqgrid_from_hdf_LambdaLayer(group_freqparams_bfreqsa, state.vertex.avertex().K1.frequencies_K1.b, Lambda_it, Lambda[0]);
    hdf5_impl::init_freqgrid_from_hdf_LambdaLayer(group_freqparams_bfreqsp, state.vertex.pvertex().K1.frequencies_K1.b, Lambda_it, Lambda[0]);
    hdf5_impl::init_freqgrid_from_hdf_LambdaLayer(group_freqparams_bfreqst, state.vertex.tvertex().K1.frequencies_K1.b, Lambda_it, Lambda[0]);
#if MAX_DIAG_CLASS>1
    hdf5_impl::init_freqgrid_from_hdf_LambdaLayer(group_freqparams_bfreqs2a,state.vertex.avertex().K2.frequencies_K2.b, Lambda_it, Lambda[0]);
    hdf5_impl::init_freqgrid_from_hdf_LambdaLayer(group_freqparams_bfreqs2p,state.vertex.pvertex().K2.frequencies_K2.b, Lambda_it, Lambda[0]);
    hdf5_impl::init_freqgrid_from_hdf_LambdaLayer(group_freqparams_bfreqs2t,state.vertex.tvertex().K2.frequencies_K2.b, Lambda_it, Lambda[0]);
    hdf5_impl::init_freqgrid_from_hdf_LambdaLayer(group_freqparams_ffreqs2a,state.vertex.avertex().K2.frequencies_K2.f, Lambda_it, Lambda[0]);
    hdf5_impl::init_freqgrid_from_hdf_LambdaLayer(group_freqparams_ffreqs2p,state.vertex.pvertex().K2.frequencies_K2.f, Lambda_it, Lambda[0]);
    hdf5_impl::init_freqgrid_from_hdf_LambdaLayer(group_freqparams_ffreqs2t,state.vertex.tvertex().K2.frequencies_K2.f, Lambda_it, Lambda[0]);
#endif
#if MAX_DIAG_CLASS>2
    hdf5_impl::init_freqgrid_from_hdf_LambdaLayer(group_freqparams_bfreqs3a,state.vertex.avertex().K3.frequencies_K3.b, Lambda_it, Lambda[0]);
    hdf5_impl::init_freqgrid_from_hdf_LambdaLayer(group_freqparams_bfreqs3p,state.vertex.pvertex().K3.frequencies_K3.b, Lambda_it, Lambda[0]);
    hdf5_impl::init_freqgrid_from_hdf_LambdaLayer(group_freqparams_bfreqs3t,state.vertex.tvertex().K3.frequencies_K3.b, Lambda_it, Lambda[0]);
    hdf5_impl::init_freqgrid_from_hdf_LambdaLayer(group_freqparams_ffreqs3a,state.vertex.avertex().K3.frequencies_K3.f, Lambda_it, Lambda[0]);
    hdf5_impl::init_freqgrid_from_hdf_LambdaLayer(group_freqparams_ffreqs3p,state.vertex.pvertex().K3.frequencies_K3.f, Lambda_it, Lambda[0]);
    hdf5_impl::init_freqgrid_from_hdf_LambdaLayer(group_freqparams_ffreqs3t,state.vertex.tvertex().K3.frequencies_K3.f, Lambda_it, Lambda[0]);
#endif

    file_out.close();

    return state;
}


/*
void write_to_hdf_group(const std::vector<double>& data, H5::Group& group, const std::string& dataset_name) {
    hsize_t dims[1] = {data.size()};
    H5::DataSpace mydataspace(1, dims);
    H5::DataSet mydataset = group.createDataSet(dataset_name, H5::PredType::NATIVE_DOUBLE, mydataspace);
    mydataset.write(data.data(), H5::PredType::NATIVE_DOUBLE);
    mydataset.close();
    mydataspace.close();
}
void write_to_hdf_group(const std::vector<comp>& data, H5::Group& group, const std::string& dataset_name) {

    // Create the memory data type for storing complex numbers in file
    H5::CompType mtype_comp = def_mtype_comp();                                 /// TODO: Can I make this static?
    h5_comp fillvalue_vert;
    fillvalue_vert.re = 0;
    fillvalue_vert.im = 0;
    H5::DSetCreatPropList plist_vert;
    plist_vert.setFillValue(mtype_comp, &fillvalue_vert);

    hsize_t dims[1] = {data.size()};
    H5::DataSpace mydataspace(1, dims);
    H5::DataSet mydataset = group.createDataSet(dataset_name, H5::PredType::NATIVE_DOUBLE, mydataspace, fillvalue_vert);
    mydataset.write(data.data(), H5::PredType::NATIVE_DOUBLE);
    mydataset.close();
    mydataspace.close();
}

void write_to_hdf_group_LambdaLayer(const std::vector<double>& data, H5::Group& group, const std::string& dataset_name, const hsize_t Lambda_it, const hsize_t nLambda_layers) {
    hsize_t start[2];
    hsize_t stride[2];
    hsize_t count[2];
    hsize_t block[2];

    start[0] = Lambda_it;
    start[1] = 0;
    for (int i = 0; i < 2; i++) {
        stride[i] = 1;
        block[i] = 1;
    }
    count[0] = 1;           // dimension of Lambda layer
    count[1] = data.size(); // flat dimension

    hsize_t RANK = 1;
    hsize_t dims[RANK+1] = {nLambda_layers, data.size()};
    H5::DataSpace file_space(RANK+1, dims);
    H5::DataSpace mem_space(RANK, dims);
    file_space.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
    H5::DataSet mydataset = group.createDataSet(dataset_name, H5::PredType::NATIVE_DOUBLE, file_space);
    mydataset.write(data.data(), H5::PredType::NATIVE_DOUBLE, mem_space, file_space);
    mydataset.close();
    mem_space.close();
    file_space.close();
}
*/


rvec read_Lambdas_from_hdf(const H5std_string FILE_NAME){



    // Open the file. Access rights: read-only
    H5::H5File *file = 0;
    file = new H5::H5File(FILE_NAME, H5F_ACC_RDONLY);

    H5::DataSet lambda_dataset = file->openDataSet("lambdas");

    // Prepare a buffer to which data from file is written. Buffer is copied to result eventually.
    //Buffer buffer;

    // Create the memory data type for storing complex numbers in file
    //H5::CompType mtype_comp = def_mtype_comp();

    // Create the dimension arrays for objects in file and in buffer
    //Dims dims(buffer, Lambda_size);

    // Read the data sets from the file to copy their content into the buffer
    //DataSets dataSets(file);

    // Create the data spaces for the data sets in file and for buffer objects
    H5::DataSpace dataSpace_Lambda = lambda_dataset.getSpace();
    /*
     * Get the number of dimensions in the dataspace.
     */
    int rank = dataSpace_Lambda.getSimpleExtentNdims();
    /*
     * Get the dimension size of each dimension in the dataspace and
     * display them.
     */
    hsize_t dims_out[2];
    int ndims = dataSpace_Lambda.getSimpleExtentDims( dims_out, NULL);
    //std::cout << "rank " << rank << ", dimensions " <<
    //     (unsigned long)(dims_out[0]) << " x " <<
    //     (unsigned long)(dims_out[1]) << std::endl;


    H5::DataSpace dataSpace_Lambda_buffer(1, dims_out);

    /// load data into buffer

    //hsize_t start_1D[1] = {Lambda_it};
    //hsize_t stride_1D[1]= {1};
    //hsize_t count_1D[1] = {1};
    //hsize_t block_1D[1] = {1};
    //dataSpaces_Lambda.selectHyperslab(H5S_SELECT_SET, count_1D, start_1D, stride_1D, block_1D);
    double Lambdas_arr[dims_out[0]];
    lambda_dataset.read(Lambdas_arr, H5::PredType::NATIVE_DOUBLE,
                        dataSpace_Lambda_buffer, dataSpace_Lambda);

    rvec Lambdas (Lambdas_arr, Lambdas_arr + sizeof(Lambdas_arr) / sizeof(double));

    // Terminate
    lambda_dataset.close();
    file->close();
    delete file;

    return Lambdas;


}

State<state_datatype> read_hdf(const H5std_string FILE_NAME, size_t Lambda_it){
    State<state_datatype> result(Lambda_ini);   // Initialize with ANY frequency grid, read grid from HDF file later
    long Lambda_size = read_Lambdas_from_hdf(FILE_NAME).size();
    if (Lambda_it < Lambda_size) {

        // Open the file. Access rights: read-only
        H5::H5File *file = 0;
        file = new H5::H5File(FILE_NAME, H5F_ACC_RDONLY);

        // Prepare a buffer to which data from file is written. Buffer is copied to result eventually.
        Buffer buffer;

        // Create the memory data type for storing complex numbers in file
        H5::CompType mtype_comp = def_mtype_comp();

        // Create the dimension arrays for objects in file and in buffer
        Dims dims(buffer, Lambda_size);

        // Read the data sets from the file to copy their content into the buffer
        DataSets dataSets(file);

        // Create the data spaces for the data sets in file and for buffer objects
        H5::DataSpace dataSpaces_Lambda = dataSets.lambda.getSpace();
        H5::DataSpace dataSpaces_freq_params;
        try {   // storing frequency gri parameters was implemented later --> old files do not have it
            dataSpaces_freq_params = dataSets.freq_params.getSpace();
        }
        catch (H5::DataSetIException error) {
            error.printErrorStack();
        }
        H5::DataSpace dataSpaces_selfenergy = dataSets.self.getSpace();
        H5::DataSpace dataSpaces_irreducible = dataSets.irred.getSpace();

        H5::DataSpace dataSpaces_Lambda_buffer(1, dims.Lambda);
        H5::DataSpace dataSpaces_freq_params_buffer(RANK_freqs-1, dims.freq_params_buffer_dims);
        H5::DataSpace dataSpaces_selfenergy_buffer(RANK_self-1, dims.selfenergy_buffer);
        H5::DataSpace dataSpaces_irreducible_buffer(RANK_irreducible-1, dims.irreducible_buffer);

#if MAX_DIAG_CLASS >= 1
        H5::DataSpace dataSpaces_K1_a = dataSets.K1_a.getSpace();
        H5::DataSpace dataSpaces_K1_p = dataSets.K1_p.getSpace();
        H5::DataSpace dataSpaces_K1_t = dataSets.K1_t.getSpace();

        H5::DataSpace dataSpaces_K1_a_buffer(RANK_K1-1, dims.K1_buffer);
        H5::DataSpace dataSpaces_K1_p_buffer(RANK_K1-1, dims.K1_buffer);
        H5::DataSpace dataSpaces_K1_t_buffer(RANK_K1-1, dims.K1_buffer);
#endif

#if MAX_DIAG_CLASS >= 2
        H5::DataSpace dataSpaces_K2_a = dataSets.K2_a.getSpace();
        H5::DataSpace dataSpaces_K2_p = dataSets.K2_p.getSpace();
        H5::DataSpace dataSpaces_K2_t = dataSets.K2_t.getSpace();

        H5::DataSpace dataSpaces_K2_a_buffer(RANK_K2-1, dims.K2_buffer);
        H5::DataSpace dataSpaces_K2_p_buffer(RANK_K2-1, dims.K2_buffer);
        H5::DataSpace dataSpaces_K2_t_buffer(RANK_K2-1, dims.K2_buffer);
#endif

#if MAX_DIAG_CLASS >= 3
        H5::DataSpace dataSpaces_K3_a = dataSets.K3_a.getSpace();
        H5::DataSpace dataSpaces_K3_p = dataSets.K3_p.getSpace();
        H5::DataSpace dataSpaces_K3_t = dataSets.K3_t.getSpace();

        H5::DataSpace dataSpaces_K3_a_buffer(RANK_K3-1, dims.K3_buffer);
        H5::DataSpace dataSpaces_K3_p_buffer(RANK_K3-1, dims.K3_buffer);
        H5::DataSpace dataSpaces_K3_t_buffer(RANK_K3-1, dims.K3_buffer);
#endif


        //Select hyperslab in the file where the data should be located
        hsize_t start[2];
        hsize_t stride[2];
        hsize_t count[2];
        hsize_t block[2];

        start[0] = Lambda_it;
        start[1] = 0;
        for (int i = 0; i < 2; i++) {
            stride[i] = 1;
            block[i] = 1;
        }
        count[0] = 1;

        /// load data into buffer

        //hsize_t start_1D[1] = {Lambda_it};
        //hsize_t stride_1D[1]= {1};
        //hsize_t count_1D[1] = {1};
        //hsize_t block_1D[1] = {1};
        //dataSpaces_Lambda.selectHyperslab(H5S_SELECT_SET, count_1D, start_1D, stride_1D, block_1D);
        //dataSets.lambda.read(buffer.lambda, H5::PredType::NATIVE_DOUBLE,
        //                   dataSpaces_Lambda_buffer, dataSpaces_Lambda);

        int rank = dataSpaces_Lambda.getSimpleExtentNdims();
        /*
         * Get the dimension size of each dimension in the dataspace and
         * display them.
         */
        hsize_t dims_out[2];
        int ndims = dataSpaces_Lambda.getSimpleExtentDims( dims_out, NULL);
        //std::cout << "rank " << rank << ", dimensions " <<
        //     (unsigned long)(dims_out[0]) << " x " <<
        //     (unsigned long)(dims_out[1]) << std::endl;


        H5::DataSpace dataSpace_Lambda_buffer(1, dims_out);

        /// load data into buffer
        double Lambdas_arr[dims_out[0]];
        dataSets.lambda.read(Lambdas_arr, H5::PredType::NATIVE_DOUBLE,
                             dataSpace_Lambda_buffer, dataSpaces_Lambda);
        *buffer.lambda = Lambdas_arr[Lambda_it];


        count[1] = N_freq_params;
        try {   // storing frequency gri parameters was implemented later --> old files do not have it
            dataSpaces_freq_params.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
            dataSets.freq_params.read(buffer.freq_params, H5::PredType::NATIVE_DOUBLE,
                                      dataSpaces_freq_params_buffer, dataSpaces_freq_params);
        }
        catch (H5::DataSpaceIException error) {
            error.printErrorStack();
        }

        count[1] = buffer.self_dim;
        dataSpaces_selfenergy.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
        dataSets.self.read(buffer.selfenergy, mtype_comp,
                           dataSpaces_selfenergy_buffer, dataSpaces_selfenergy);

        count[1] = buffer.irred_dim;
        dataSpaces_irreducible.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
        dataSets.irred.read(buffer.irreducible_class, mtype_comp,
                            dataSpaces_irreducible_buffer, dataSpaces_irreducible);


#if MAX_DIAG_CLASS >= 1
        count[1] = buffer.K1_dim;
        dataSpaces_K1_a.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
        dataSpaces_K1_p.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
        dataSpaces_K1_t.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);

        dataSets.K1_a.read(buffer.K1_class_a, mtype_comp, dataSpaces_K1_a_buffer, dataSpaces_K1_a);
        dataSets.K1_p.read(buffer.K1_class_p, mtype_comp, dataSpaces_K1_p_buffer, dataSpaces_K1_p);
        dataSets.K1_t.read(buffer.K1_class_t, mtype_comp, dataSpaces_K1_t_buffer, dataSpaces_K1_t);
#endif
#if MAX_DIAG_CLASS >= 2
        count[1] = buffer.K2_dim;
        dataSpaces_K2_a.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
        dataSpaces_K2_p.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
        dataSpaces_K2_t.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);

        dataSets.K2_a.read(buffer.K2_class_a, mtype_comp, dataSpaces_K2_a_buffer, dataSpaces_K2_a);
        dataSets.K2_p.read(buffer.K2_class_p, mtype_comp, dataSpaces_K2_p_buffer, dataSpaces_K2_p);
        dataSets.K2_t.read(buffer.K2_class_t, mtype_comp, dataSpaces_K2_t_buffer, dataSpaces_K2_t);
#endif
#if MAX_DIAG_CLASS >= 3
        count[1] = buffer.K3_dim;
        dataSpaces_K3_a.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
        dataSpaces_K3_p.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
        dataSpaces_K3_t.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);

        dataSets.K3_a.read(buffer.K3_class_a, mtype_comp, dataSpaces_K3_a_buffer, dataSpaces_K3_a);
        dataSets.K3_p.read(buffer.K3_class_p, mtype_comp, dataSpaces_K3_p_buffer, dataSpaces_K3_p);
        dataSets.K3_t.read(buffer.K3_class_t, mtype_comp, dataSpaces_K3_t_buffer, dataSpaces_K3_t);

#endif

        // Initialize the frequency grids of the result State using the parameters stored in buffer
        result_set_frequency_grids(result, buffer);

        // Copy the buffered result into State object
        copy_buffer_to_result(result, buffer);

        // Terminate
        dataSets.close(true);
        file->close();
        delete file;

        return result;

    } else {
        throw std::runtime_error("Cannot read from file " + FILE_NAME + " since Lambda layer out of range");
    }
}



void test_hdf5(H5std_string FILE_NAME, int i, State<state_datatype>& state) {
    // test hdf5: read files and compare to original file
    int cnt = 0;
    State<state_datatype> out = read_hdf(FILE_NAME, i);
    for (int iK=0; iK<2; ++iK) {
        for (int iSE = 0; iSE < nSE; ++iSE) {
            if (state.selfenergy.val(iK, iSE, 0) != out.selfenergy.val(iK, iSE, 0)) {
                std::cout << "Self-energy not equal, " << iK << ", " << iSE << std::endl;
                cnt += 1;
            }
        }
    }

    int spin = 0;
    for (int iK=0; iK<nK_K1; ++iK) {
        for (int i_in=0; i_in<n_in; ++i_in) {
#if MAX_DIAG_CLASS >= 1
            for (int iw1=0; iw1<nBOS; ++iw1) {
                if (state.vertex.avertex().K1.val(iK, spin, iw1, i_in) != out.vertex.avertex().K1.val(iK, spin, iw1, i_in)) {
                    std::cout << "Vertex not equal, " << iK << ", " << iw1 << std::endl;
                    cnt += 1;
                }
                if (state.vertex.pvertex().K1.val(iK, spin, iw1, i_in) != out.vertex.pvertex().K1.val(iK, spin, iw1, i_in)) {
                    std::cout << "Vertex not equal, " << iK << ", " << iw1 << std::endl;
                    cnt += 1;
                }
                if (state.vertex.tvertex().K1.val(iK, spin, iw1, i_in) != out.vertex.tvertex().K1.val(iK, spin, iw1, i_in)) {
                    std::cout << "Vertex not equal, " << iK << ", " << iw1 << std::endl;
                    cnt += 1;
                }
#if MAX_DIAG_CLASS >= 2
                for (int iw2=0; iw2<nFER; ++iw2) {
                    if (state.vertex.avertex().K2.val(iK, spin, iw1, iw2, i_in) != out.vertex.avertex().K2.val(iK, spin, iw1, iw2, i_in)) {
                        std::cout << "Vertex not equal, " << iK << ", " << iw1 << ", " << iw2 << std::endl;
                        cnt += 1;
                    }
                    if (state.vertex.pvertex().K2.val(iK, spin, iw1, iw2, i_in) != out.vertex.pvertex().K2.val(iK, spin, iw1, iw2, i_in)) {
                        std::cout << "Vertex not equal, " << iK << ", " << iw1 << ", " << iw2 << std::endl;
                        cnt += 1;
                    }
                    if (state.vertex.tvertex().K2.val(iK, spin, iw1, iw2, i_in) != out.vertex.tvertex().K2.val(iK, spin, iw1, iw2, i_in)) {
                        std::cout << "Vertex not equal, " << iK << ", " << iw1 << ", " << iw2 << std::endl;
                        cnt += 1;
                    }
#if MAX_DIAG_CLASS == 3
                    for (int iw3=0; iw3<nFER; ++iw3) {
                        if (state.vertex.avertex().K3.val(iK, spin, iw1, iw2, iw3, i_in) != out.vertex.avertex().K3.val(iK, spin, iw1, iw2, iw3, i_in)) {
                            std::cout << "Vertex not equal, " << iK << ", " << iw1 << ", " << iw2 << ", " << iw3 << std::endl;
                            cnt += 1;
                        }
                        if (state.vertex.pvertex().K3.val(iK, spin, iw1, iw2, iw3, i_in) != out.vertex.pvertex().K3.val(iK, spin, iw1, iw2, iw3, i_in)) {
                            std::cout << "Vertex not equal, " << iK << ", " << iw1 << ", " << iw2 << ", " << iw3 << std::endl;
                            cnt += 1;
                        }
                        if (state.vertex.tvertex().K3.val(iK, spin, iw1, iw2, iw3, i_in) != out.vertex.tvertex().K3.val(iK, spin, iw1, iw2, iw3, i_in)) {
                            std::cout << "Vertex not equal, " << iK << ", " << iw1 << ", " << iw2 << ", " << iw3 << std::endl;
                            cnt += 1;
                        }
                    }
#endif
                }
#endif
            }
#endif
        }
    }
    if (cnt == 0) print("HDF5 test successful.", true);
    else print("HDF5 test failed. Number of differences: ", cnt, true);
}


bool test_read_write_data_hdf(bool verbose) {
    //if (verbose) print("Testing HDF input and output:", true);

    vec<comp> vector_original(10, glb_i);
    std::array<std::size_t,3> dims_mularr = {1,2,3};
    multidimensional::multiarray<comp,3> mularr_original(dims_mularr, glb_i);

    std::string testfile_name = "test.h5";
    //if (verbose) print("File name: " + testfile_name, true);

    // write to hdf file
    H5::H5File file_out(testfile_name.c_str(), H5F_ACC_TRUNC);
    //if (verbose) print("Created file.", true);
    //H5::H5std_string FILE_NAME(data_dir + "test.h5");
    H5::Group group_out(file_out.createGroup("/data"));
    //if (verbose) print("Created group.", true);

    write_to_hdf(group_out, "vector", vector_original, false);
    write_to_hdf(group_out, "mularr", mularr_original, false);
    //if (verbose) print("Written to hdf.", true);

    // read from hdf file
    hid_t file_id = H5Fopen(testfile_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    H5::Group group = H5Gopen(file_id, "data", H5P_DEFAULT);
    //if (verbose) print("Loaded group.", true);
    vec<comp> vector_input;
    read_from_hdf<comp>(group, "vector", vector_input);
    //if (verbose) print("Loaded vector.", true);
    multidimensional::multiarray<comp,3> mularr_input;
    read_from_hdf<comp,3>(group, "mularr", mularr_input);
    //if (verbose) print("Loaded multiarray.", true);

    //if (verbose) print("Length of original vec: ", vector_original.size(), true);
    //if (verbose) print("Length of loaded vec: ", vector_input.size(), true);

    comp deviation_vec = (vector_input - vector_original).max_norm();
    if (std::abs(deviation_vec) < 1e-10) print("Read / write of vector to HDF file successful.", true);
    else print("PROBLEM during read / write of vector to HDF file.", true);

    comp deviation_mularr=(mularr_input -mularr_original).maxabs();
    if (std::abs(deviation_mularr) < 1e-10) print("Read / write of multiarray to HDF file successful.", true);
    else print("PROBLEM during read / write of multiarray to HDF file.", true);



    //H5::H5std_string FILE_NAME(data_dir + "test.h5");
    H5::Group group_out_L(file_out.createGroup("/data_with_LambdaLayers"));
    //if (verbose) pprint("Created group.", true);

    write_to_hdf_LambdaLayer(group_out_L, "vector", vector_original, 1, 10,  false);
    write_to_hdf_LambdaLayer(group_out_L, "mularr", mularr_original, 1, 10,  false);
    //if (verbose) pprint("Written to hdf LambdaLayer.", true);

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
            print("Read / write of vector to LambdaLayer of HDF file successful.", true);
        else {
            passed = false;
            print("PROBLEM during read / write of vector to LambdaLayer of HDF file.", true);
        }
        if (std::abs(deviation_mularr_L) < 1e-10)
            print("Read / write of multiarray to LambdaLayer ofHDF file successful.", true);
        else {
            passed = false;
            print("PROBLEM during read / write to LambdaLayer of multiarray to HDF file.", true);
        }
    }

    return passed;
}


bool test_read_write_state_hdf(bool verbose) {
    //if (verbose) print("Testing HDF input and output of states:", true);
    const double Lambda = 1.;
    const int Lambda_it = 1;
    const int numberLambdaLayers = 10;
    State<state_datatype> state_output(Lambda);
    state_output.initialize();
    state_output = state_output + 1.;

    write_state_to_hdf("test_state.h5", state_output, numberLambdaLayers);
    //print("Written state", true);


    State<state_datatype> state_input = read_state_from_hdf_LambdaLayer("test_state.h5", 0);
    //print("Read state", true);

    State<state_datatype> state_diff = state_output - state_input;

    add_state_to_hdf("test_state.h5", state_output*2, Lambda_it);
    add_state_to_hdf("test_state.h5", state_input, Lambda_it+1);

    write_state_to_hdf("test_statediff.h5", state_diff, numberLambdaLayers);

    bool passed = true;

    if (verbose) {
        if (state_diff.norm() < 1e-10) print("Read / write of state data to LambdaLayer of HDF file successful.", true);
        else {
            passed = false;
            print("PROBLEM during read / write of state data to LambdaLayer of HDF file. deviation : ",
                  state_diff.norm(), true);
        }
        if ((state_output.selfenergy.frequencies.get_ws_vec() -
             state_input.selfenergy.frequencies.get_ws_vec()).max_norm() < 1e-10)
            print("Read / write of frequency grid to LambdaLayer of HDF file successful.", true);
        else {
            passed = false;
            print("PROBLEM during read / write of frequency grid to LambdaLayer to HDF file.", true);
        }
    }

    return passed;

}

