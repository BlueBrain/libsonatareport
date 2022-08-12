#include <algorithm>
#include <iostream>
#include <memory>

#include "../library/implementation_interface.hpp"
#include "../library/sonatareport.h"
#include "sonata_data.h"

namespace bbp {
namespace sonata {

SonataData::SonataData(const std::string& report_name,
                       const std::string& population_name,
                       uint64_t population_offset,
                       size_t max_buffer_size,
                       int num_steps,
                       double dt,
                       double tstart,
                       double tend,
                       const std::string& units,
                       std::shared_ptr<nodes_t> nodes,
                       hid_t file_handler)
    : report_name_(report_name)
    , population_name_(population_name)
    , report_units_(units)
    , population_offset_(population_offset)
    , num_steps_(num_steps)
    , hdf5_writer_(std::make_unique<HDF5Writer>(report_name, file_handler))
    , nodes_(nodes) {
    prepare_buffer(max_buffer_size);
    index_pointers_.resize(nodes->size());

    // Round the tstart value to avoid conflicts in case of save-restore
    time_ = {round(tstart / dt) * dt, tend, dt};
    reporting_period_ = static_cast<int>(dt / SonataReport::atomic_step_);
    last_step_recorded_ = tstart / SonataReport::atomic_step_;
    last_step_ = tend / SonataReport::atomic_step_;
}

SonataData::SonataData(const std::string& report_name)
    : report_name_(report_name)
    , population_name_("All")
    , population_offset_(0)
    , hdf5_writer_(std::make_unique<HDF5Writer>(report_name)) {}

void SonataData::prepare_buffer(size_t max_buffer_size) {
    if (SonataReport::rank_ == 0) {
        logger->debug("PREPARING BUFFER for report {} and population {}",
                      report_name_,
                      population_name_);
    }

    for (auto& kv : *nodes_) {
        total_elements_ += kv.second->get_num_elements();
    }

    // Calculate the timesteps that fit given the buffer size
    {
        // Ranks without elements still need to participate in the writings
        uint32_t max_steps_to_write = (total_elements_ == 0)
                                          ? std::numeric_limits<uint32_t>::max()
                                          : max_buffer_size / (sizeof(float) * total_elements_);
        uint32_t common_max_steps_to_write =
            Implementation::get_max_steps_to_write(report_name_, max_steps_to_write);
        if (common_max_steps_to_write < num_steps_) {  // More steps asked that buffer can contain
            if (common_max_steps_to_write < SonataReport::min_steps_to_record_) {
                steps_to_write_ = SonataReport::min_steps_to_record_;
            } else {
                // Minimum 1 timestep required to write
                uint32_t min_steps_to_write = 1;
                steps_to_write_ = std::max(common_max_steps_to_write, min_steps_to_write);
            }
        } else {  // all the steps asked fit into the given buffer
            // If the buffer size is bigger that all the timesteps needed to record we allocate only
            // the amount of timesteps
            steps_to_write_ = num_steps_;
        }
    }

    remaining_steps_ = num_steps_;

    if (SonataReport::rank_ == 0) {
        logger->debug("\t- Total elements: {}", total_elements_);
        logger->debug("\t- Num steps: {}", num_steps_);
        logger->debug("\t- Steps to write: {}", steps_to_write_);
        logger->debug("\t- Max Buffer size: {}", max_buffer_size);
    }

    size_t buffer_size = total_elements_ * steps_to_write_;
    report_buffer_.resize(buffer_size);

    if (SonataReport::rank_ == 0) {
        logger->debug("\t-Buffer size: {} (count={})",
                      buffer_size * sizeof(report_buffer_[0]),
                      buffer_size);
    }
}

bool SonataData::is_due_to_report(double step) const noexcept {
    // Dont record data if current step < tstart
    if (step < last_step_recorded_) {
        return false;
        // Dont record data if current step > tend
    } else if (step > last_step_) {
        return false;
        // Dont record data if is not a reporting step (step%period)
    } else if (static_cast<int>(step - last_step_recorded_) % reporting_period_ != 0) {
        return false;
    }
    return true;
}

void SonataData::record_data(double step, const std::vector<uint64_t>& node_ids) {
    // Calculate the offset to write into the buffer
    uint32_t offset = static_cast<uint32_t>((step - last_step_recorded_) / reporting_period_);
    uint32_t local_position = last_position_ + total_elements_ * offset;

    if (steps_recorded_ >= num_steps_) {
        if (SonataReport::rank_ == 0) {
            logger->trace("Already recorded {} steps, skipping report recording for {} population",
                          steps_recorded_,
                          population_name_);
        }
        return;
    }

    if (SonataReport::rank_ == 0) {
        logger->trace(
            "Recording data for population {}, step={} last_step_recorded={} steps recorded {} "
            "first node_id={} "
            "buffer_size={} "
            "and offset={}",
            population_name_,
            step,
            last_step_recorded_,
            steps_recorded_,
            node_ids[0],
            report_buffer_.size(),
            local_position);
    }
    for (auto& kv : *nodes_) {
        uint64_t current_node_id = kv.second->get_node_id();
        // Check if node is set to be recorded (found in nodeids)
        if (std::find(node_ids.begin(), node_ids.end(), current_node_id) != node_ids.end()) {
            kv.second->fill_data(report_buffer_.begin() + local_position);
            nodes_recorded_.insert(current_node_id);
        }
        local_position += kv.second->get_num_elements();
    }

    // Increase steps recorded when all nodes from specific rank has been already recorded
    if (nodes_recorded_.size() == nodes_->size()) {
        steps_recorded_++;
    }
}

void SonataData::record_data(double step) {
    uint32_t local_position = last_position_;

    if (steps_recorded_ >= num_steps_) {
        if (SonataReport::rank_ == 0) {
            logger->trace("Already recorded {} steps, skipping report recording for {} population",
                          steps_recorded_,
                          population_name_);
        }
        return;
    }

    if (SonataReport::rank_ == 0) {
        logger->trace(
            "Recording data for step={} last_step_recorded={} buffer_size={} and offset={}",
            step,
            last_step_recorded_,
            report_buffer_.size(),
            local_position);
    }
    for (auto& kv : *nodes_) {
        kv.second->fill_data(report_buffer_.begin() + local_position);
        local_position += kv.second->get_num_elements();
    }
    current_step_++;
    last_position_ += total_elements_;
    last_step_recorded_ += reporting_period_;

    if (current_step_ == steps_to_write_) {
        write_data();
    }
}

void SonataData::check_and_write(double timestep) {
    if (remaining_steps_ <= 0) {
        return;
    }

    if (SonataReport::rank_ == 0) {
        logger->debug("UPDATING timestep t={} for report {} and population {}",
                      timestep,
                      report_name_,
                      population_name_);
    }
    current_step_ += steps_recorded_;
    last_position_ += total_elements_ * steps_recorded_;
    last_step_recorded_ += reporting_period_ * steps_recorded_;
    nodes_recorded_.clear();

    // Write when buffer is full, finish all remaining recordings or when record several steps in a
    // row
    if (current_step_ == steps_to_write_ || current_step_ == remaining_steps_ ||
        steps_recorded_ > 1) {
        if (SonataReport::rank_ == 0) {
            logger->trace(
                "Writing to file {}! population {} steps_to_write={}, current_step={}, "
                "remaining_steps={} "
                "steps_recorded={}",
                report_name_,
                population_name_,
                steps_to_write_,
                current_step_,
                remaining_steps_,
                steps_recorded_);
        }
        write_data();
    }
    steps_recorded_ = 0;
}

void SonataData::prepare_dataset() {
    if (SonataReport::rank_ == 0) {
        logger->debug("PREPARING HEADER for report {} and population {}",
                      report_name_,
                      population_name_);
    }
    // Prepare /report
    for (auto& kv : *nodes_) {
        // /report
        const std::vector<uint32_t> element_ids = kv.second->get_element_ids();
        element_ids_.insert(element_ids_.end(), element_ids.begin(), element_ids.end());
        node_ids_.push_back(kv.second->get_node_id());
    }
    int element_offset = Implementation::get_offset(report_name_, total_elements_);
    logger->trace("\tRank {} - Total elements are: {} and element offset is: {}",
                  SonataReport::rank_,
                  total_elements_,
                  element_offset);

    int last_writer = (total_elements_ == 0) ? -1 : SonataReport::rank_;
    int last_rank = Implementation::get_last_rank(report_name_, last_writer);
    if (SonataReport::rank_ == last_rank) {
        index_pointers_.resize(nodes_->size() + 1);
    }

    // Prepare index pointers
    if (!index_pointers_.empty()) {
        index_pointers_[0] = element_offset;
    }
    for (size_t i = 1; i < index_pointers_.size(); i++) {
        uint64_t previous_node_id = node_ids_[i - 1];
        index_pointers_[i] = index_pointers_[i - 1] +
                             nodes_->at(previous_node_id)->get_num_elements();
    }

    // All ranks need to participate in the write
    write_report_header();
}

void SonataData::convert_gids_to_sonata(std::vector<uint64_t>& node_ids,
                                        uint64_t population_offset) {
    if (getenv("LIBSONATA_ZERO_BASED_GIDS") == nullptr) {
        std::transform(std::begin(node_ids),
                       std::end(node_ids),
                       std::begin(node_ids),
                       [& population_offset = population_offset](int x) {
                           if (x == 0) {
                               throw std::runtime_error(
                                   "Error: node_id is 0 and input data is reported as 1-based");
                           }
                           return x - population_offset - 1;
                       });
    }
}

void SonataData::write_report_header() {
    // TODO: remove configure_group and add it to write_any()
    const std::string reports_population_group = "/report/" + population_name_;
    hdf5_writer_->configure_group("/report");
    hdf5_writer_->configure_group(reports_population_group);
    hdf5_writer_->configure_group(reports_population_group + "/mapping");
    hdf5_writer_->configure_dataset(reports_population_group + "/data",
                                    num_steps_,
                                    total_elements_);
    hdf5_writer_->configure_attribute(reports_population_group + "/data", "units", report_units_);

    std::vector<uint64_t> sonata_node_ids(node_ids_);
    convert_gids_to_sonata(sonata_node_ids, population_offset_);
    hdf5_writer_->write(reports_population_group + "/mapping/node_ids", sonata_node_ids);
    hdf5_writer_->write(reports_population_group + "/mapping/index_pointers", index_pointers_);
    hdf5_writer_->write(reports_population_group + "/mapping/element_ids", element_ids_);
    hdf5_writer_->write_time(reports_population_group + "/mapping/time", time_);
    hdf5_writer_->configure_attribute(reports_population_group + "/mapping/time", "units", "ms");
}

void SonataData::write_spikes_header(Population& population) {
    const std::string order_by = population.get_sorting();
    if (order_by != "by_time" && order_by != "by_id" && order_by != "none") {
        throw std::runtime_error("Order method " + order_by + "does not exists");
    }

    const std::string spikes_population_group = "/spikes/" + population.get_population_name();
    hdf5_writer_->configure_group("/spikes");
    hdf5_writer_->configure_group(spikes_population_group);
    hdf5_writer_->configure_enum_attribute(spikes_population_group, "sorting", order_by);
    hsize_t timestamps_size =
        Implementation::get_global_dims(report_name_, population.get_spike_timestamps().size());
    Implementation::sort_spikes(population.get_spike_timestamps(),
                                population.get_spike_node_ids(),
                                order_by);
    hdf5_writer_->write(spikes_population_group + "/timestamps", population.get_spike_timestamps());
    if (timestamps_size > 0) {
        hdf5_writer_->configure_attribute(spikes_population_group + "/timestamps", "units", "ms");
    }
    std::vector<uint64_t> sonata_spike_node_ids(population.get_spike_node_ids());
    convert_gids_to_sonata(sonata_spike_node_ids, population.get_population_offset());
    hdf5_writer_->write(spikes_population_group + "/node_ids", sonata_spike_node_ids);
}

void SonataData::add_population(std::unique_ptr<Population>&& population) {
    populations_.push_back(std::move(population));
}

void SonataData::write_spike_populations() {
    for (const auto& population : populations_) {
        write_spikes_header(*population);
    }
}

void SonataData::write_data() {
    if (remaining_steps_ <= 0) {  // Nothing left to write
        return;
    }
    if (current_step_ >= remaining_steps_) {  // Avoid writing out of bounds
        current_step_ = remaining_steps_;
    }
    if (SonataReport::rank_ == 0) {
        logger->debug("WRITING timestep data to file {} in population {}",
                      report_name_,
                      population_name_);
    }
    hdf5_writer_->write_2D(report_buffer_, current_step_, total_elements_);
    remaining_steps_ -= current_step_;
    if (SonataReport::rank_ == 0) {
        logger->debug("\t-Steps written: {}", current_step_);
        logger->debug("\t-Remaining steps: {}", remaining_steps_);
    }
    last_position_ = 0;
    current_step_ = 0;
}

void SonataData::close() {
    hdf5_writer_->close();
}

Population::Population(const std::string& population_name,
                       uint64_t population_offset,
                       const std::string& order_by,
                       const std::vector<double>& spike_timestamps,
                       const std::vector<uint64_t>& spike_node_ids)
    : population_name_(population_name.empty() ? "All" : population_name)
    , population_offset_(population_offset)
    , order_by_(order_by.empty() ? "by_time" : order_by)
    , spike_timestamps_(spike_timestamps)
    , spike_node_ids_(spike_node_ids) {}

}  // namespace sonata
}  // namespace bbp
