#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>

#ifdef SONATA_REPORT_HAVE_MPI
#include <mpi.h>
#endif

#include <bbp/sonata/reports.h>
#include <utils/logger.h>

static void show_usage(const std::string& name) {
    std::cerr << "Usage: " << name << " <spike_filename> [population_name]\n"
              << "Options:\n"
              << "\t-h,--help\t\tShow this help message\n"
              << "\tspike_filename\t\tName of the report to convert\n"
              << "\tpopulation_name\t\tName of the report population (Default 'All')\n"
              << "Example:\n"
              << "\t " << name << " out.dat PopulationA" << std::endl;
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
    if ((argc != 2 && argc != 3) || help(argc, argv)) {
        show_usage(argv[0]);
        return -1;
    }

    const std::string file_name = argv[1];
    std::string population_name = "All";
    if (argc == 3) {
        population_name = argv[2];
    }
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
    const std::string report_name = file_name.substr(file_name.find_last_of("/\\") + 1);
    sonata_create_spikefile(".", report_name.data());

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

    logger->info("File '{}' successfully converted to '{}.h5'", file_name, report_name);

#ifdef SONATA_REPORT_HAVE_MPI
    MPI_Finalize();
#endif
    return 0;
}