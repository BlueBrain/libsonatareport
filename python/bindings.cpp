#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <bbp/sonata/report_reader.h>

#include "generated/docstrings.h"

#include <cstdint>
#include <memory>
#include <string>

namespace py = pybind11;

using namespace pybind11::literals;

using namespace bbp::sonata;


namespace {

// Return a new Numpy array with data owned by another python object
// This avoids copies, and enables correct reference counting for memory keep-alive
template <typename DATA_T, typename DIMS_T, typename OWNER_T>
py::array managedMemoryArray(const DATA_T* data, const DIMS_T& dims, const OWNER_T& owner) {
    const auto& tinfo = py::detail::get_type_info(typeid(OWNER_T));
    const auto& handle = py::detail::get_object_handle(&owner, tinfo);
    return py::array(dims, data, handle);
}

// create a macro to reduce repetition for docstrings
#define DOC_SPIKEREADER_POP(x) DOC(bbp, sonata, SpikeReader, Population, x)
#define DOC_SPIKEREADER(x) DOC(bbp, sonata, SpikeReader, x)
#define DOC_REPORTREADER_POP(x) DOC(bbp, sonata, ReportReader, Population, x)
}

namespace pybind11 {
namespace detail {
template <typename T>
struct type_caster<nonstd::optional<T>>: optional_caster<nonstd::optional<T>> {};

template <>
struct type_caster<nonstd::nullopt_t>: public void_caster<nonstd::nullopt_t> {};
}  // namespace detail
}  // namespace pybind11


template <typename ReportType, typename KeyType>
void bindReportReader(py::module& m, const std::string& prefix) {
    py::class_<DataFrame<KeyType>>(m,
                                   (prefix + "DataFrame").c_str(),
                                   "A container of raw reporting data, compatible with Pandas")
        .def_readonly("ids", &DataFrame<KeyType>::ids)

        // .data and .time members are owned by this c++ object. We can't do std::move.
        // To avoid copies we must declare the owner of the data as the current python
        // object. Numpy will adjust owner reference count according to returned arrays
        // clang-format off
        .def_property_readonly("data", [](const DataFrame<KeyType>& dframe) {
            std::array<ssize_t, 2> dims {0l, ssize_t(dframe.ids.size())};
            if (dims[1] > 0) {
                dims[0] = dframe.data.size() / dims[1];
            }
            return managedMemoryArray(dframe.data.data(), dims, dframe);
        })
        // clang-format on
        .def_property_readonly("times", [](DataFrame<KeyType>& dframe) {
            return managedMemoryArray(dframe.times.data(), dframe.times.size(), dframe);
        });

    py::class_<typename ReportType::Population>(m,
                                                (prefix + "ReportPopulation").c_str(),
                                                "A population inside a ReportReader")
        .def("get",
             &ReportType::Population::get,
             "Return reports with all those node_ids between 'tstart' and 'tstop' with a stride "
             "tstride",
             "node_ids"_a = nonstd::nullopt,
             "tstart"_a = nonstd::nullopt,
             "tstop"_a = nonstd::nullopt,
             "tstride"_a = nonstd::nullopt)
        .def("get_node_ids",
             &ReportType::Population::getNodeIds,
             "Return the list of nodes ids for this population")
        .def_property_readonly("sorted",
                               &ReportType::Population::getSorted,
                               DOC_REPORTREADER_POP(getSorted))
        .def_property_readonly("times",
                               &ReportType::Population::getTimes,
                               DOC_REPORTREADER_POP(getTimes))
        .def_property_readonly("time_units",
                               &ReportType::Population::getTimeUnits,
                               DOC_REPORTREADER_POP(getTimeUnits))
        .def_property_readonly("data_units",
                               &ReportType::Population::getDataUnits,
                               DOC_REPORTREADER_POP(getDataUnits));
    py::class_<ReportType>(m, (prefix + "ReportReader").c_str(), "Used to read somas files")
        .def(py::init([](py::object h5_filepath) { return ReportType(py::str(h5_filepath)); }),
             "h5_filepath"_a)
        .def("get_population_names", &ReportType::getPopulationNames, "Get list of all populations")
        .def("__getitem__", &ReportType::openPopulation);
}


PYBIND11_MODULE(_libsonatareport, m) {
    py::class_<SpikeReader::Population>(m, "SpikePopulation", "A population inside a SpikeReader")
        .def("get",
             &SpikeReader::Population::get,
             "Return spikes with all those node_ids between 'tstart' and 'tstop'",
             "node_ids"_a = nonstd::nullopt,
             "tstart"_a = nonstd::nullopt,
             "tstop"_a = nonstd::nullopt)
        .def_property_readonly(
            "sorting",
            [](const SpikeReader::Population& self) {
                auto s = self.getSorting();
                if (s == SpikeReader::Population::Sorting::by_id)
                    return "by_id";
                if (s == SpikeReader::Population::Sorting::by_time)
                    return "by_time";
                return "none";
            },
            DOC_SPIKEREADER_POP(getSorting));
    py::class_<SpikeReader>(m, "SpikeReader", "Used to read spike files")
        .def(py::init([](py::object h5_filepath) { return SpikeReader(py::str(h5_filepath)); }),
             "h5_filepath"_a)
        .def("get_population_names",
             &SpikeReader::getPopulationNames,
             DOC_SPIKEREADER(getPopulationNames))
        .def("__getitem__", &SpikeReader::openPopulation);

    bindReportReader<SomaReportReader, NodeID>(m, "Soma");
    bindReportReader<ElementReportReader, std::pair<NodeID, uint32_t>>(m, "Element");
}
