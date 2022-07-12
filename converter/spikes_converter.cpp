#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>

#include <bbp/sonata/reports.h>
#include <utils/logger.h>

int main() {

    std::ifstream infile("out.dat");
    if (!infile) {
        logger->error("File out.dat not found!");
        return 0;
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
    sonata_create_spikefile(".", "out");

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

    logger->info("File 'out.dat' successfully converted to 'out.h5'");

    return 0;
}