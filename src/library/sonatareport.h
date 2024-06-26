#pragma once
#include <spdlog/spdlog.h>
#include <string>
#include <unordered_map>
#include <vector>

#ifdef SONATA_REPORT_HAVE_MPI
#include <mpi.h>
#endif

#include "report.h"

namespace bbp {
namespace sonata {

/**
 *  \brief Contains and manages the reports
 */
class SonataReport
{
    using reports_t = std::map<std::string, std::shared_ptr<Report>>;
#ifdef SONATA_REPORT_HAVE_MPI
    using communicators_t = std::unordered_map<std::string, MPI_Comm>;
#endif
  public:
    static double atomic_step_;
    static double min_steps_to_record_;
#ifdef SONATA_REPORT_HAVE_MPI
    static MPI_Comm has_nodes_;
    static communicators_t communicators_;
#endif
    static int rank_;
    static bool first_report;

    /**
     * \brief Destroy all report objects.
     * This should invoke their destructor which will close the report
     * file and clean up
     */
    void clear();
    bool is_empty();

    int get_num_reports() const noexcept {
        return reports_.size();
    }

    std::shared_ptr<Report> create_report(const std::string& name,
                                          const std::string& kind,
                                          double tstart,
                                          double tend,
                                          double dt,
                                          const std::string& units);

    std::shared_ptr<Report> get_report(const std::string& name) const;

    bool report_exists(const std::string& name) const;

    void create_communicators();
    void prepare_datasets();
    void setup_reports();

    void create_spikefile(const std::string& output_dir, const std::string& filename = "out");
    void add_spikes_population(const std::string& population_name,
                               uint64_t population_offset,
                               const std::vector<double>& spike_timestamps,
                               const std::vector<uint64_t>& spike_node_ids,
                               const std::string& order_by = "by_time");
    void write_spike_populations();
    void close_spikefile();

    template <typename Functor, typename T>
    void apply_all(const Functor& functor, T data) {
        std::for_each(reports_.begin(), reports_.end(), [&](reports_t::value_type arg) {
            functor(arg.second, data);
        });
    }

  private:
    reports_t reports_;
    std::unique_ptr<SonataData> spike_data_;
    bool reports_initialized_ = false;
};

}  // namespace sonata
}  // namespace bbp
