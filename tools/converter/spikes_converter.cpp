#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>

#include <bbp/sonata/reports.h>
#include <utils/logger.h>

static void show_usage(std::string name) {
    std::cerr << "Usage: " << name << " <spike_filename>\n"
              << "Options:\n"
              << "\t-h,--help\tShow this help message\n"
              << "Example:\n"
              << "\t " << name << " out.dat" << std::endl;
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
    if (argc != 2 || help(argc, argv)) {
        show_usage(argv[0]);
        return -1;
    }

    const std::string file_name = argv[1];
    logger->info("Trying to convert '{}' binary report...'", file_name);
    std::ifstream infile(file_name);
    if (!infile) {
        logger->error("File {} not found!", file_name);
        return -2;
    }

    // remove /scatter
    std::string scatter;
    getline(infile, scatter);
    if (scatter != "/scatter") {
        infile.seekg(0, std::ios::beg);
    }

    std::vector<double> spike_timestamps;
    std::vector<int> spike_node_ids;
    double timestamp;
    int node_id;

    while (infile >> timestamp >> node_id) {
        spike_timestamps.push_back(timestamp);
        spike_node_ids.push_back(node_id);
    }

    // Create a spike file
    sonata_create_spikefile(".", file_name.data());

    std::string population_name = "All";
    uint64_t population_offset = 0;
    sonata_add_spikes_population(population_name.data(),
                                 population_offset,
                                 spike_timestamps.data(),
                                 spike_timestamps.size(),
                                 spike_node_ids.data(),
                                 spike_node_ids.size());
    sonata_write_spike_populations();
    // Close the spike file
    sonata_close_spikefile();

    logger->info("File '{}' successfully converted to '{}.h5'", file_name, file_name);

    return 0;
}