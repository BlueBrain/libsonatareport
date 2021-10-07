#include <chrono>
#include <cmath>
#include <iostream>
#include <iterator>
#include <thread>
#include <vector>

#ifdef SONATA_REPORT_HAVE_MPI
#include <mpi.h>
#endif

#include <bbp/sonata/reports.h>
#include <utils/logger.h>

struct Cell {
    uint64_t node_id;
    std::string kind;  // soma / element
    std::vector<double> voltages;
};

void generate_spikes(const std::vector<uint64_t>& nodeids,
                     std::vector<double>& spike_timestamps,
                     std::vector<int>& spike_node_ids,
                     double tstart,
                     double tstop,
                     int seed,
                     int max_size) {
    // Generate 10,30,50,70,90 spikes
    size_t num_spikes = (10 + (20 * seed)) % 100;
    spike_timestamps.reserve(num_spikes);
    spike_node_ids.reserve(num_spikes);
    for (size_t i = 0; i < num_spikes; i++) {
        // timestamp between tstart and tstop
        double timestamp = tstart + (0.5 + seed) / (max_size / (tstop - tstart));
        // get an index to the nodeids
        size_t index = seed % nodeids.size();
        uint32_t node_id = nodeids[index];
        spike_timestamps.push_back(timestamp);
        spike_node_ids.push_back(node_id);
    }
}

void generate_elements(Cell& cell, int seed) {
    // 50+-5 elements
    size_t num_elements = 50 + ((seed % 10) - 5);
    if (cell.kind == "soma") {
        num_elements = 1;
    }
    cell.voltages.reserve(num_elements);
    for (size_t j = 0; j < num_elements; j++) {
        cell.voltages.push_back(seed % 10);
    }
}

std::vector<uint64_t> generate_data(std::vector<Cell>& cells,
                                    const std::string& kind,
                                    int seed) {
    std::vector<uint64_t> nodeids;
    // Each nodeid starts with the 1000 + rank*10 (i.e. rank 5 will have nodeids: 1051, 1052,
    // 1053...)
    uint64_t next_nodeid = 1000 + 1 + seed * 10;

    // 5+-5 cells
    uint32_t num_cells = 5 + ((2 + (seed % 10)) - 5);
    nodeids.reserve(num_cells);
    for (uint32_t i = 0; i < num_cells; i++) {
        Cell tmp_cell;
        tmp_cell.kind = kind;

        nodeids.push_back(next_nodeid);
        tmp_cell.node_id = next_nodeid++;

        // element or soma
        generate_elements(tmp_cell, seed);
        cells.push_back(tmp_cell);
    }
    return nodeids;
}

void init(const char* report_name,
          const char* population_name,
          uint64_t population_offset,
          double tstart,
          double tstop,
          double dt,
          std::vector<Cell>& cells,
          const std::string& kind,
          const std::string& units) {
    // logic for registering soma and element reports with reportinglib
    sonata_create_report(report_name, tstart, tstop, dt, units.c_str(), kind.c_str());
    for (auto& cell : cells) {
        sonata_add_node(report_name, population_name, population_offset, cell.node_id);
        int element_id = cell.node_id * 1000;

        for (auto& element : cell.voltages) {
            sonata_add_element(report_name, population_name, cell.node_id, element_id, &element);
            ++element_id;
        }
    }
}

void change_data(std::vector<Cell>& cells, std::vector<float>& buffered_data) {
    // Increment in 1 per timestep every voltage
    for (auto& cell : cells) {
        for (auto& element : cell.voltages) {
            buffered_data.push_back(element);
            element++;
        }
    }
}

void print_data(std::vector<Cell>& cells) {
    for (auto& cell : cells) {
        std::cout << "++NEURON node_id: " << cell.node_id << "\nelements:\n";
        std::copy(cell.voltages.begin(),
                  cell.voltages.end(),
                  std::ostream_iterator<double>(std::cout, ", "));
        std::cout << "\n\n";
    }
}

int main() {
    logger->set_level(spdlog::level::trace);
    int global_rank = 0;
    int global_size = 1;
#ifdef SONATA_REPORT_HAVE_MPI
    MPI_Init(nullptr, nullptr);
    MPI_Comm_rank(MPI_COMM_WORLD, &global_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &global_size);
#endif
    if (global_rank == 0) {
        logger->info("Starting...");
    }

    const double dt = 0.1;
    const double tstart = 0.0;
    const double tstop = 0.3;

    std::vector<Cell> element_cells;
    std::vector<Cell> soma_cells;
    std::vector<uint64_t> element_nodeids;
    std::vector<uint64_t> soma_nodeids;
    std::vector<double> spike_timestamps;
    std::vector<int> spike_node_ids;

    // Each rank will get different number of nodes (some even 0, so will be idle ranks)
    element_nodeids = generate_data(element_cells, "compartment", global_rank);
    soma_nodeids = generate_data(soma_cells, "soma", global_rank);
    generate_spikes(
        soma_nodeids, spike_timestamps, spike_node_ids, tstart, tstop, global_rank, global_size);

    std::vector<int> int_element_nodeids(begin(element_nodeids), end(element_nodeids));
    std::vector<int> int_soma_nodeids(begin(soma_nodeids), end(soma_nodeids));

    if (global_rank == 0) {
        logger->info("Initializing data structures (reports, nodes, elements)");
    }
    const char* element_report = "compartment_report";
    const char* soma_report = "soma_report";
    const char* population_name = "All";
    const char* units = "mV";
    uint64_t population_offset = 0;
    std::vector<float> element_buffered_data = {};
    std::vector<float> soma_buffered_data = {};

    init(element_report,
         population_name,
         population_offset,
         tstart,
         tstop,
         dt,
         element_cells,
         "compartment",
         units);
    init(soma_report,
         population_name,
         population_offset,
         tstart,
         tstop,
         dt,
         soma_cells,
         "soma",
         units);
    sonata_set_max_buffer_size_hint(20);
    sonata_set_atomic_step(dt);

    sonata_setup_communicators();
    sonata_prepare_datasets();
    sonata_time_data();

    if (global_rank == 0) {
        logger->info("Starting the simulation!");
    }
    // Calculate number of steps of the simulation
    double sim_steps = (tstop - tstart) / dt;
    int num_steps = static_cast<int>(std::ceil(sim_steps));
    double t = 0.0;
    for (int i = 0; i < num_steps; i++) {
        if (global_rank == 0) {
            logger->info("Recording data for step = {}", i);
        }
        sonata_record_node_data(i, element_nodeids.size(), &int_element_nodeids[0], element_report);
        sonata_record_node_data(i, soma_nodeids.size(), &int_soma_nodeids[0], soma_report);
        // Also works
        // sonata_rec(i);
        std::this_thread::sleep_for(std::chrono::milliseconds(500));

        // Update timestep on reportinglib
        sonata_check_and_flush(t);
        t += dt;
        // Change data every timestep
        change_data(element_cells, element_buffered_data);
        change_data(soma_cells, soma_buffered_data);
    }
    sonata_flush(t);
    sonata_clear();

    // Write pre-buffered data (GPU usecase)
    if (global_rank == 0) {
        logger->info("GPU usecase: Writting prebuffered steps");
    }
    const char* buffered_soma_report = "buffered_soma_report";
    const char* buffered_population_name = "buffered";

    init(buffered_soma_report,
         buffered_population_name,
         population_offset,
         tstart,
         tstop,
         dt,
         soma_cells,
         "soma",
         units);

    sonata_setup_communicators();
    sonata_prepare_datasets();
    sonata_write_buffered_data(buffered_soma_report, soma_buffered_data.data(), soma_buffered_data.size(), num_steps);
    sonata_clear();

    const std::string output_dir = ".";

    // Create a spike file
    sonata_create_spikefile(output_dir.data());

    std::vector<std::string> population_names{"NodeA", "NodeB"};
    std::vector<uint64_t> population_offsets{0, 1000};
    // Write the spikes for populations "NodeA" and "NodeB"
    for (size_t i = 0; i < population_names.size(); i++) {
        sonata_add_spikes_population(population_names[i].data(),
                                     population_offsets[i],
                                     spike_timestamps.data(),
                                     spike_timestamps.size(),
                                     spike_node_ids.data(),
                                     spike_node_ids.size());
    }

    sonata_write_spike_populations();
    // Close the spike file
    sonata_close_spikefile();

    if (global_rank == 0) {
        logger->info("Finalizing...");
    }

#ifdef SONATA_REPORT_HAVE_MPI
    MPI_Finalize();
#endif
    return 0;
}
