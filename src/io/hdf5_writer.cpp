#include <iostream>
#include <tuple>

#include "../library/implementation_interface.hpp"
#include "../library/sonatareport.h"
#include "hdf5_writer.h"

namespace bbp {
namespace sonata {

template void HDF5Writer::write<uint32_t>(const std::string& dataset_name,
                                          const std::vector<uint32_t>& buffer);
template void HDF5Writer::write<uint64_t>(const std::string& dataset_name,
                                          const std::vector<uint64_t>& buffer);
template void HDF5Writer::write<float>(const std::string& dataset_name,
                                       const std::vector<float>& buffer);
template void HDF5Writer::write<double>(const std::string& dataset_name,
                                        const std::vector<double>& buffer);

HDF5Writer::HDF5Writer(const std::string& report_name, hid_t file_handler)
    : report_name_(report_name)
    , file_(file_handler) {
    collective_list_ = Implementation::initialize_colective();
}

HDF5Writer::HDF5Writer(const std::string& report_name)
    : report_name_(report_name) {
    file_ = Implementation::prepare_write(report_name);
    collective_list_ = Implementation::initialize_colective();

    // Create enum type for the ordering of the spikes
    spikes_attr_type_ = H5Tenum_create(H5T_STD_U8LE);
    uint8_t val;
    H5Tenum_insert(spikes_attr_type_, "none", (val = 0, &val));
    H5Tenum_insert(spikes_attr_type_, "by_id", (val = 1, &val));
    H5Tenum_insert(spikes_attr_type_, "by_time", (val = 2, &val));
}

void HDF5Writer::configure_group(const std::string& group_name) {
    // create group if it doestn't exist
    if (!H5Lexists(file_, group_name.c_str(), H5P_DEFAULT)) {
        hid_t group = H5Gcreate2(file_, group_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Gclose(group);
    }
}

void HDF5Writer::configure_attribute(const std::string& dataset_name,
                                     const std::string& attribute_name,
                                     const std::string& attribute_value) {
    logger->trace("Configuring attribute '{}' for group name '{}' with value {}",
                  attribute_name,
                  dataset_name,
                  attribute_value);
    hid_t dataset_id = H5Dopen(file_, dataset_name.data(), H5P_DEFAULT);
    hid_t attr_space = H5Screate(H5S_SCALAR);
    hid_t type = H5Tcopy(H5T_C_S1);
    H5Tset_size(type, H5T_VARIABLE);
    H5Tset_cset(type, H5T_CSET_UTF8);
    const char* value[] = {attribute_value.data()};
    hid_t attr_id =
        H5Acreate2(dataset_id, attribute_name.data(), type, attr_space, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr_id, type, value);

    H5Aclose(attr_id);
    H5Sclose(attr_space);
    H5Dclose(dataset_id);
}

void HDF5Writer::configure_enum_attribute(const std::string& group_name,
                                          const std::string& attribute_name,
                                          const std::string& attribute_value) {
    hid_t group_id = H5Gopen(file_, group_name.data(), H5P_DEFAULT);
    hid_t attr_space = H5Screate(H5S_SCALAR);

    hid_t attr_id = H5Acreate2(
        group_id, attribute_name.c_str(), spikes_attr_type_, attr_space, H5P_DEFAULT, H5P_DEFAULT);
    uint8_t attr_value;
    H5Tenum_valueof(spikes_attr_type_, attribute_value.data(), &attr_value);
    H5Awrite(attr_id, spikes_attr_type_, &attr_value);

    H5Aclose(attr_id);
    H5Sclose(attr_space);
    H5Gclose(group_id);
}

void HDF5Writer::configure_dataset(const std::string& dataset_name,
                                   uint32_t total_steps,
                                   uint32_t total_elements) {
    std::array<hsize_t, 2> dims = {total_steps,
                                   Implementation::get_global_dims(report_name_, total_elements)};
    hid_t data_space = H5Screate_simple(2, dims.data(), nullptr);
    dataset_ = H5Dcreate(file_,
                         dataset_name.c_str(),
                         H5T_IEEE_F32LE,
                         data_space,
                         H5P_DEFAULT,
                         H5P_DEFAULT,
                         H5P_DEFAULT);

    // Caculate the offset of each rank
    offset_[1] = Implementation::get_offset(report_name_, total_elements);
    H5Sclose(data_space);
}

void HDF5Writer::write_2D(const std::vector<float>& buffer,
                          uint32_t steps_to_write,
                          uint32_t total_elements) {
    hsize_t dims = buffer.size();
    hsize_t global_dims = Implementation::get_global_dims(report_name_, dims);

    if (global_dims > 0) {
        std::array<hsize_t, 2> count = {steps_to_write, total_elements};
        hid_t memspace = H5Screate_simple(2, count.data(), nullptr);
        hid_t filespace = H5Dget_space(dataset_);

        H5Sselect_hyperslab(
            filespace, H5S_SELECT_SET, offset_.data(), nullptr, count.data(), nullptr);
        H5Dwrite(dataset_, H5T_NATIVE_FLOAT, memspace, filespace, collective_list_, buffer.data());

        H5Sclose(filespace);
        H5Sclose(memspace);
    }
    offset_[0] += steps_to_write;
}

template <typename T>
void HDF5Writer::write(const std::string& dataset_name, const std::vector<T>& buffer) {
    if (SonataReport::rank_ == 0) {
        logger->debug("WRITING dataset {} in report {}", dataset_name, report_name_);
    }
    hid_t type = h5typemap::get_h5_type(T(0));

    hsize_t dims = buffer.size();
    hsize_t global_dims = Implementation::get_global_dims(report_name_, dims);
    hsize_t offset = Implementation::get_offset(report_name_, dims);

    hid_t data_space = H5Screate_simple(1, &global_dims, nullptr);
    hid_t data_set = H5Dcreate(
        file_, dataset_name.c_str(), type, data_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (global_dims > 0) {
        hid_t memspace = H5Screate_simple(1, &dims, nullptr);
        hid_t filespace = H5Dget_space(data_set);

        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &offset, nullptr, &dims, nullptr);
        H5Dwrite(data_set, type, memspace, filespace, collective_list_, buffer.data());

        H5Sclose(filespace);
        H5Sclose(memspace);
    }
    H5Dclose(data_set);
    H5Sclose(data_space);
}

void HDF5Writer::write_time(const std::string& dataset_name, const std::array<double, 3>& buffer) {
    hsize_t dims = buffer.size();
    hid_t data_space = H5Screate_simple(1, &dims, nullptr);
    hid_t data_set = H5Dcreate(file_,
                               dataset_name.c_str(),
                               H5T_NATIVE_DOUBLE,
                               data_space,
                               H5P_DEFAULT,
                               H5P_DEFAULT,
                               H5P_DEFAULT);

    H5Dwrite(data_set, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer.data());
    H5Dclose(data_set);
    H5Sclose(data_space);
}

void HDF5Writer::close() {
    // We close the dataset "/data", the spike enum type and the hdf5 file
    if (dataset_) {
        H5Dclose(dataset_);
        dataset_ = 0;
    }
    if (spikes_attr_type_) {
        H5Tclose(spikes_attr_type_);
        spikes_attr_type_ = 0;
        if (file_) {
            H5Fclose(file_);
            file_ = 0;
        }
    }
}

}  // namespace sonata
}  // namespace bbp
