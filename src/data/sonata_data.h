#pragma once
#include <map>
#include <set>

#include "../io/hdf5_writer.h"
#include "node.h"

namespace bbp {
namespace sonata {

class Population
{
  public:
    Population(const std::string& population_name,
               uint64_t population_offset,
               const std::string& order_by,
               const std::vector<double>& spike_timestamps,
               const std::vector<uint64_t>& spike_node_ids);

    const std::string& get_population_name() const noexcept {
        return population_name_;
    }
    uint64_t get_population_offset() const noexcept {
        return population_offset_;
    }
    const std::string& get_sorting() const noexcept {
        return order_by_;
    }
    void set_sorting(const std::string& order_by) noexcept {
        order_by_ = order_by;
    }
    const std::vector<double>& get_spike_timestamps() const noexcept {
        return spike_timestamps_;
    }
    const std::vector<uint64_t>& get_spike_node_ids() const noexcept {
        return spike_node_ids_;
    }
    std::vector<double>& get_spike_timestamps() {
        return spike_timestamps_;
    }
    std::vector<uint64_t>& get_spike_node_ids() {
        return spike_node_ids_;
    }

  private:
    std::string population_name_;
    uint64_t population_offset_;
    std::string order_by_;
    std::vector<double> spike_timestamps_;
    std::vector<uint64_t> spike_node_ids_;
};

class SonataData
{
  public:
    SonataData(const std::string& report_name,
               const std::string& population_name,
               uint64_t population_offset,
               size_t max_buffer_size,
               int num_steps,
               double dt,
               double tstart,
               double tend,
               const std::string& units,
               std::shared_ptr<nodes_t> nodes);

    SonataData(const std::string& report_name);

    void prepare_dataset();
    void write_report_header();
    void write_spikes_header(Population& population);
    void write_spike_populations();
    void add_population(std::unique_ptr<Population>&& population);

    void write_data(const std::vector<float>& buffered_data, uint32_t steps_to_write);
    void flush();
    void close();

    bool is_due_to_report(double step) const noexcept;
    void record_data(double step, const std::vector<uint64_t>& node_ids);
    void record_data(double step);
    void check_and_write(double timestep);
    void convert_gids_to_sonata(std::vector<uint64_t>& node_ids, uint64_t population_offset);

    const std::vector<float>& get_report_buffer() const noexcept {
        return report_buffer_;
    }

    const std::vector<uint64_t>& get_node_ids() const noexcept {
        return node_ids_;
    }
    const std::vector<uint64_t>& get_index_pointers() const noexcept {
        return index_pointers_;
    }
    const std::vector<uint32_t>& get_element_ids() const noexcept {
        return element_ids_;
    }

  private:
    std::string report_name_;
    std::string population_name_;
    std::string report_units_;
    uint64_t population_offset_;
    std::vector<float> report_buffer_;
    uint32_t total_elements_ = 0;
    uint32_t num_steps_ = 0;
    uint32_t steps_to_write_ = 0;
    uint32_t current_step_ = 0;
    uint32_t steps_recorded_ = 0;
    uint32_t previous_position_ = 0;
    uint32_t remaining_steps_ = 0;
    uint32_t reporting_period_ = 0;
    double previous_step_recorded_ = 0.;
    double final_step_ = 0.;

    std::vector<uint64_t> node_ids_;
    std::vector<uint64_t> index_pointers_;
    std::vector<uint32_t> element_ids_;
    std::array<double, 3> time_;

    std::set<uint64_t> nodes_recorded_;
    const std::unique_ptr<HDF5Writer> hdf5_writer_;
    std::shared_ptr<nodes_t> nodes_;
    std::vector<std::unique_ptr<Population>> populations_;

    void prepare_buffer(size_t max_buffer_size);
};

}  // namespace sonata
}  // namespace bbp
