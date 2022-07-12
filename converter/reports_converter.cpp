#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <map>

#include <bbp/sonata/reports.h>
#include <utils/logger.h>
#include "binary_reader/ReadManager.h"

using namespace bbpReader;

int main(int argc, char *argv[]) {

    if(argc!=2) {
        logger->error("Wrong number of arguments.");
        logger->info("Try: ./reports_converter example.bbp");
        return 0;
    }
    const char* file_name = argv[1];
    logger->info("Trying to convert '{}' binary report...'", file_name);
    std::ifstream f(file_name);
    if(!f.good()) {
        logger->error("File '{}' not found!", file_name);
        return 0;
    }

    FrameParser frame_parser(file_name);
    Header header = frame_parser.getHeader();

    std::string filename = file_name;
    std::string report_name = filename.substr(0, filename.find("."));

    // Get header information in order to create the report
    double tstart = header.getStartTime();
    double tstop = header.getEndTime();
    double dt = header.getTimeStepSize();
    sonata_create_report(report_name.data(), tstart, tstop, dt, "mV", "compartment");

    // Get Cell information to create node/element structure
    // TODO: check if it could be done without getting the FrameParser per gid
    ReadManager file(file_name);
    FrameInfo cells = file.retrieveAllCellInfo();
    std::vector<CellID> node_ids(1);
    for (FrameInfo::const_iterator it = cells.begin(); it != cells.end(); ++it) {
        node_ids[0] = it->getCellNum();
        sonata_add_node(report_name.data(), "All", 0, it->getCellNum());

        FrameParser frame_parser_gid(file_name, node_ids);
        int num_element_ids = frame_parser_gid.getBufferSize_elements();
        DataItem* element_ids_buffer = new DataItem[num_element_ids];
        frame_parser_gid.readFrameMapping(element_ids_buffer);

        std::vector<uint32_t> element_ids(element_ids_buffer, element_ids_buffer + num_element_ids);
        for (auto& element : element_ids) {
            sonata_add_element(report_name.data(), "All", node_ids[0], element, nullptr);
        }
    }
    // Generate the initial structure
    sonata_prepare_datasets();

    int total_num_element_ids = frame_parser.getBufferSize_elements();
    DataItem* all_element_ids_buffer = new DataItem[total_num_element_ids];
    // Write the timestep frames
    while (frame_parser.hasMore()) {
        frame_parser.readFrameData(all_element_ids_buffer);
        std::vector<float> data_element_ids(all_element_ids_buffer, all_element_ids_buffer + total_num_element_ids);
        sonata_write_buffered_data(report_name.data(), data_element_ids.data(), data_element_ids.size(), 1);
    }
    sonata_clear();
    logger->info("File '{}' successfully converted to '{}.h5'", file_name, report_name);
    return 0;
}