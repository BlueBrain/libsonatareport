#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>

#ifdef SONATA_REPORT_HAVE_MPI
#include <mpi.h>
#endif

#include "binary_reader/ReadManager.h"
#include <bbp/sonata/reports.h>
#include <utils/logger.h>

using namespace bbpReader;

static void show_usage(std::string name) {
    std::cerr << "Usage: " << name << " <report_filename> <option(s)>\n"
              << "Options:\n"
              << "\t-h,--help\t\t\tShow this help message\n"
              << "\t<[--soma, --compartment]>\tSelect soma or compartment report\n"
              << "Examples:\n"
              << "\t " << name << " soma.bbp --soma\n"
              << "\t " << name << " compartment.bbp --compartment" << std::endl;
}

bool help(int argc, char* argv[]) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            return true;
        }
    }
    return false;
}

int main(int argc, char* argv[]) {
#ifdef SONATA_REPORT_HAVE_MPI
    MPI_Init(nullptr, nullptr);
#endif
    if (argc != 3 || help(argc, argv)) {
        show_usage(argv[0]);
        return -1;
    }
    const std::string file_name = argv[1];
    const std::string report_type = argv[2];
    logger->info("Trying to convert '{}' binary report...'", file_name);
    std::ifstream f(file_name);
    if (!f.good()) {
        logger->error("File '{}' not found!", file_name);
        return -2;
    }

    FrameParser frame_parser(file_name);
    Header header = frame_parser.getHeader();

    std::string report_name = file_name.substr(file_name.find_last_of("/\\") + 1);

    // Get header information in order to create the report
    double tstart = header.getStartTime();
    double tstop = header.getEndTime();
    double dt = header.getTimeStepSize();
    logger->info("Report info: tstart = '{}', tstop = '{}', dt = '{}'", tstart, tstop, dt);
    sonata_create_report(report_name.data(), tstart, tstop, dt, "mV", "compartment");

    // Get Cell information to create node/element structure
    // TODO: check if it could be done without getting the FrameParser per gid
    ReadManager file(file_name);
    FrameInfo cells = file.retrieveAllCellInfo();
    std::vector<CellID> node_ids(1);
    for (FrameInfo::const_iterator it = cells.begin(); it != cells.end(); ++it) {
        node_ids[0] = it->getCellNum();
        sonata_add_node(report_name.data(), "All", 0, it->getCellNum());

        if (report_type == "--soma") {
            sonata_add_element(report_name.data(), "All", node_ids[0], 0, nullptr);
        } else {  // --compartment
            FrameParser frame_parser_gid(file_name, node_ids);
            int num_element_ids = frame_parser_gid.getBufferSize_elements();
            DataItem* element_ids_buffer = new DataItem[num_element_ids];
            frame_parser_gid.readFrameMapping(element_ids_buffer);

            std::vector<uint32_t> element_ids(element_ids_buffer,
                                              element_ids_buffer + num_element_ids);
            for (auto& element : element_ids) {
                sonata_add_element(report_name.data(), "All", node_ids[0], element, nullptr);
            }
        }
    }
    // Generate the initial structure
    sonata_prepare_datasets();

    uint64_t element_ids_per_frame = frame_parser.getBufferSize_elements();
    DataItem* element_ids_buffer = new DataItem[element_ids_per_frame];
    uint32_t timestep = 0;
    // Write the timestep frames
    while (frame_parser.hasMore()) {
        frame_parser.readFrameData(element_ids_buffer);
        sonata_write_buffered_data(report_name.data(),
                                   element_ids_buffer,
                                   element_ids_per_frame,
                                   1);  // Number of frames to write
        if (timestep % 1000 == 0 && timestep > 0) {
            logger->info("Writing timestep '{}'", timestep);
        }
        timestep++;
    }
    logger->info("'{}' timesteps written!", timestep);
    logger->info("File '{}' successfully converted to '{}.h5'", file_name, report_name);

#ifdef SONATA_REPORT_HAVE_MPI
    MPI_Finalize();
#endif
    return 0;
}