#include <cmath>
#include <iostream>
#include <limits>

#include "../utils/logger.h"
#include "implementation_interface.hpp"
#include "report.h"
#include "sonatareport.h"

namespace bbp {
namespace sonata {

// default buffer size of 4MB
constexpr uint64_t default_max_buffer_size = 4194304;

Report::Report(
    const std::string& report_name, double tstart, double tend, double dt, const std::string& units)
    : populations_(std::make_shared<populations_t>())
    , report_name_(report_name)
    , tstart_(tstart)
    , tend_(tend)
    , dt_(dt)
    , units_(units)
    , max_buffer_size_(default_max_buffer_size)
    , report_is_closed_(false) {
    // Calculate number of reporting steps, rounding the tstart value in case of save-restore
    tstart = round(tstart / dt) * dt;
    num_steps_ = static_cast<int>(std::ceil((tend - tstart) / dt));
}

void Report::add_node(const std::string& population_name,
                      uint64_t population_offset,
                      uint64_t node_id) {
    if (population_exists(population_name)) {
        if (node_exists(population_name, node_id)) {
            throw std::runtime_error("Warning: attempted to add node " + std::to_string(node_id) +
                                     " to the target multiple times on same node. Ignoring.");
        }
        std::shared_ptr<nodes_t> nodes = populations_->at(population_name);
        nodes->emplace(node_id, std::make_shared<Node>(node_id));
    } else {
        // node is new insert it into the map
        std::shared_ptr<nodes_t> nodes = std::make_shared<nodes_t>();
        nodes->emplace(node_id, std::make_shared<Node>(node_id));
        populations_->emplace(population_name, nodes);
        population_offsets_.emplace(population_name, population_offset);
    }
}

bool Report::node_exists(const std::string& population_name, uint64_t node_id) const {
    std::shared_ptr<nodes_t> nodes = populations_->at(population_name);
    return nodes->find(node_id) != nodes->end();
}

bool Report::population_exists(const std::string& population_name) const {
    return populations_->find(population_name) != populations_->end();
}

std::shared_ptr<Node> Report::get_node(const std::string& population_name, uint64_t node_id) const {
    return populations_->at(population_name)->at(node_id);
}
void print_comm_ranks(MPI_Comm comm)
{
   MPI_Group grp, world_grp;

   MPI_Comm_group(MPI_COMM_WORLD, &world_grp);
   MPI_Comm_group(comm, &grp);

   int grp_size;

   MPI_Group_size(grp, &grp_size);

   int *ranks = (int*)malloc(grp_size * sizeof(int));
   int *world_ranks =(int*) malloc(grp_size * sizeof(int));

   for (int i = 0; i < grp_size; i++)
      ranks[i] = i;

   MPI_Group_translate_ranks(grp, grp_size, ranks, world_grp, world_ranks);

   for (int i = 0; i < grp_size; i++)
      printf("comm[%d] has world rank %d\n", i, world_ranks[i]);

   free(ranks); free(world_ranks);

   MPI_Group_free(&grp);
   MPI_Group_free(&world_grp);
}

int Report::prepare_dataset() {
    Implementation::add_communicator(report_name_);
    file_handler_ = Implementation::prepare_write(report_name_);
    
    for (const auto& population : *populations_) {
        const std::string& population_name = population.first;
        std::string comm_name = report_name_ + "_" + population_name;
        Implementation::add_communicator(comm_name);
    }
    if (SonataReport::rank_ == 0) {
        std::shared_ptr<nodes_t> nodes = std::make_shared<nodes_t>();
        populations_->emplace("NodeB", nodes);
    } else {
        std::shared_ptr<nodes_t> nodes = std::make_shared<nodes_t>();
        populations_->emplace("NodeA", nodes);
    }
    for (const auto& population : *populations_) {
        const std::string& population_name = population.first;
        /*std::string comm_name = report_name_ + "_" + population_name;
        Implementation::add_communicator(comm_name);*/
        std::shared_ptr<nodes_t> nodes = population.second;
        sonata_populations_.push_back(
            std::make_unique<SonataData>(report_name_,
                                         population_name,
                                         population_offsets_[population_name],
                                         max_buffer_size_,
                                         num_steps_,
                                         dt_,
                                         tstart_,
                                         tend_,
                                         units_,
                                         nodes,
                                         file_handler_));
        sonata_populations_.back()->prepare_dataset();
        MPI_Barrier(MPI_COMM_WORLD);
    }
    print_comm_ranks(SonataReport::has_nodes_);
    MPI_Barrier(MPI_COMM_WORLD);
    //H5Fclose(file_handler_);
    //MPI_Abort(MPI_COMM_WORLD, 0);
    return 0;
}

void Report::record_data(double step, const std::vector<uint64_t>& node_ids) {
    for (const auto& sonata_data : sonata_populations_) {
        logger->trace("Rank {} - is due to report?", SonataReport::rank_);
        if (sonata_data->is_due_to_report(step)) {
            sonata_data->record_data(step, node_ids);
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }
}

void Report::record_data(double step) {
    for (const auto& sonata_data : sonata_populations_) {
        if (sonata_data->is_due_to_report(step)) {
            sonata_data->record_data(step);
        }
    }
}

void Report::check_and_flush(double timestep) {
    for (const auto& sonata_data : sonata_populations_) {
        sonata_data->check_and_write(timestep);
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

void Report::refresh_pointers(std::function<double*(double*)> refresh_function) {
    for (auto& population : *populations_) {
        std::shared_ptr<nodes_t> nodes = population.second;
        for (auto& kv : *nodes) {
            kv.second->refresh_pointers(refresh_function);
        }
    }
}

void Report::flush(double time) {
    if (SonataReport::rank_ == 0) {
        logger->trace("Flush() called at t={} for report {}", time, report_name_);
    }
    for (const auto& sonata_data : sonata_populations_) {
        // Write if there are any remaining steps to write
        sonata_data->write_data();
        if (time - tend_ + dt_ / 2 > 1e-6) {
            sonata_data->close();
        }
    }
    if (!report_is_closed_) {
        logger->trace("CLOSING FILE (flush)");
        H5Fclose(file_handler_);
        report_is_closed_ = true;
    }
}

void Report::set_max_buffer_size(size_t buffer_size) {
    logger->trace("Setting buffer size to {}", buffer_size);
    max_buffer_size_ = buffer_size;
}

}  // namespace sonata
}  // namespace bbp
