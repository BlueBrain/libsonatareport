#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>

#include "binary_reader/ReadManager.h"
#include <bbp/sonata/reports.h>
#include <utils/logger.h>

using namespace bbpReader;

int get_steps_to_write(const char* var_name) {
    // In general, writing 10 steps per iteration is best performance wise
    // The environment variable STEPS_TO_WRITE could be overwritten for testing
    int steps_to_write = 10;
    char* var_value = getenv(var_name);
    if (var_value != nullptr && atoi(var_value) > 0) {
        steps_to_write = atoi(var_value);
    }
    return steps_to_write;
}

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

    std::string filename = file_name;
    std::string report_name = filename.substr(filename.find_last_of("/\\") + 1);

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

    uint32_t num_steps = header.getNumberOfSteps();
    uint32_t num_frames = get_steps_to_write("STEPS_TO_WRITE");
    uint64_t element_ids_per_frame = frame_parser.getBufferSize_elements();
    uint64_t total_element_ids = element_ids_per_frame * num_frames;
    logger->info("Number of frames to write per iteration: '{}'", num_frames);
    DataItem* all_element_ids_buffer = new DataItem[total_element_ids];
    uint32_t timestep = 0;
    // Write the timestep frames
    while (frame_parser.hasMore()) {
        // Check if number of frames to write is not bigger than the remaining
        if (num_steps - timestep < num_frames) {
            num_frames = num_steps - timestep;
        }
        frame_parser.readMultipleFrameData(all_element_ids_buffer, num_frames);
        sonata_write_buffered_data(report_name.data(),
                                   all_element_ids_buffer,
                                   total_element_ids,
                                   num_frames);
        if (timestep % 1000 == 0 && timestep > 0) {
            logger->info("Writing timestep '{}'", timestep);
        }
        timestep += num_frames;
    }
    logger->info("'{}' timesteps written!", timestep);
    logger->info("File '{}' successfully converted to '{}.h5'", file_name, report_name);
    return 0;
}