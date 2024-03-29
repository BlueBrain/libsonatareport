#include <bbp/sonata/reports.h>
#include <catch2/catch.hpp>
#include <data/sonata_data.h>
#include <iostream>
#include <library/implementation_interface.hpp>
#include <memory>
#ifdef SONATA_REPORT_HAVE_MPI
#include <mpi.h>
#endif

using namespace bbp::sonata;

static const std::string population_name = "All";
static const std::string report_units = "nA";
static const uint64_t population_offset = 100;

SCENARIO("Test SonataData class", "[SonataData][IOWriter]") {
    GIVEN("A node map structure") {
        double dt = 1.0;
        double tstart = 0.0;
        double tend = 3.0;
        sonata_set_atomic_step(dt);
        using nodes_t = std::map<uint64_t, std::shared_ptr<Node>>;
        auto node = std::make_shared<Node>(101);
        double element = 10;
        double element2 = 12;
        node->add_element(&element, 0);
        node->add_element(&element2, 1);
        auto node2 = std::make_shared<Node>(102);
        node2->add_element(&element, 10);
        node2->add_element(&element2, 11);
        node2->add_element(&element2, 12);
        auto node42 = std::make_shared<Node>(142);
        std::vector<double> elements{34.1, 55.21, 3.141592, 44, 2124, 42.42};
        int i = 20;
        for (double& elem : elements) {
            node42->add_element(&elem, i);
            ++i;
        }
        WHEN("We record some data and prepare the dataset for a big enough max buffer size") {
            auto nodes = std::make_shared<nodes_t>(std::initializer_list<nodes_t::value_type>{
                {101, node}, {102, node2}, {142, node42}});

            int num_steps = 3;
            size_t max_buffer_size = 1024;
            std::string report_name = "test_sonatadata";
            hid_t file_handler = Implementation::prepare_write(report_name);
            auto sonata = std::make_unique<SonataData>(report_name,
                                                       population_name,
                                                       population_offset,
                                                       max_buffer_size,
                                                       num_steps,
                                                       dt,
                                                       tstart,
                                                       tend,
                                                       report_units,
                                                       nodes,
                                                       file_handler);
            std::vector<uint64_t> nodeids_1 = {101, 142};
            std::vector<uint64_t> nodeids_2 = {102};

            sonata->prepare_dataset();
            for (int i = 0; i < num_steps; i++) {
                sonata->record_data(i, nodeids_1);
                sonata->record_data(i, nodeids_2);
                sonata->check_and_write(i);
            }
            sonata->close();
            H5Fclose(file_handler);

            THEN(
                "The buffer size is the total number of steps times the total number of elements") {
                // 1024 / (sizeof(float) * 11) = 23.27 > 3 (total number of steps)
                // buffer_size = 11 * 3
                REQUIRE(sonata->get_report_buffer().size() == 33);
            }

            THEN("We check the node ids of the sonata report") {
                const std::vector<uint64_t> node_ids = sonata->get_node_ids();
                std::vector<uint64_t> compare = {101, 102, 142};
                REQUIRE(node_ids == compare);
            }

            THEN("We check the node ids of the sonata report after applying offset") {
                const std::vector<uint64_t> node_ids = sonata->get_node_ids();
                std::vector<uint64_t> sonata_node_ids(node_ids);
                sonata->convert_gids_to_sonata(sonata_node_ids, population_offset);
                std::vector<uint64_t> compare = {0, 1, 41};
                REQUIRE(sonata_node_ids == compare);
                REQUIRE_THROWS(sonata->convert_gids_to_sonata(sonata_node_ids, population_offset));
            }

            THEN("We check the element ids of the sonata report") {
                const std::vector<uint32_t> element_ids = sonata->get_element_ids();
                std::vector<uint32_t> compare = {0, 1, 10, 11, 12, 20, 21, 22, 23, 24, 25};
                REQUIRE(element_ids == compare);
            }

            THEN("We check the index pointers of the sonata report") {
                const std::vector<uint64_t> index_pointers = sonata->get_index_pointers();
                std::vector<uint64_t> compare = {0, 2, 5, 11};
                REQUIRE(index_pointers == compare);
            }
        }
        WHEN("We record some other data and prepare the dataset for a small max buffer size") {
            auto nodes = std::make_shared<nodes_t>(std::initializer_list<nodes_t::value_type>{
                {101, node}, {102, node2}, {142, node42}});
            int num_steps = 3;
            size_t max_buffer_size = 128;
            std::string report_name = "test_sonatadata2";
            hid_t file_handler = Implementation::prepare_write(report_name);
            auto sonata2 = std::make_unique<SonataData>(report_name,
                                                        population_name,
                                                        population_offset,
                                                        max_buffer_size,
                                                        num_steps,
                                                        dt,
                                                        tstart,
                                                        tend,
                                                        report_units,
                                                        nodes,
                                                        file_handler);

            sonata2->prepare_dataset();
            for (int i = 0; i < num_steps; i++) {
                sonata2->record_data(i);
                sonata2->check_and_write(i);
            }
            sonata2->close();
            H5Fclose(file_handler);

            THEN(
                "The buffer size is the number of steps to write that fit on the buffer times the "
                "total elements") {
                // 128 / (sizeof(float) * 11) = 2 < 3 (total number of steps)
                // buffer_size = 11 * 2
                REQUIRE(sonata2->get_report_buffer().size() == 22);
            }
        }
    }
    GIVEN("Spike data") {
        std::vector<double> spike_timestamps{0.3, 0.1, 0.2, 1.3, 0.7};
        std::vector<uint64_t> spike_node_ids{3, 5, 2, 3, 2};
        uint64_t population_offset = 0;
        std::string report_name = "spikes";
        auto sonata_spikes = std::make_unique<SonataData>(report_name);
        auto sonata_population = std::make_unique<Population>(
            population_name, population_offset, "by_time", spike_timestamps, spike_node_ids);
        Population* pop_p = sonata_population.get();
        WHEN("We add the population and write the spikes ordered by time") {
            sonata_spikes->add_population(std::move(sonata_population));
            sonata_spikes->write_spike_populations();
            THEN("We check that the spike nodes ids are ordered according to timestamps") {
                const std::vector<uint64_t> node_ids = pop_p->get_spike_node_ids();
                std::vector<uint64_t> compare = {5, 2, 3, 2, 3};
                REQUIRE(node_ids == compare);
            }
            THEN("We check that the spike timestamps are in order") {
                const std::vector<double> timestamps = pop_p->get_spike_timestamps();
                std::vector<double> compare = {0.1, 0.2, 0.3, 0.7, 1.3};
                REQUIRE(timestamps == compare);
            }
        }
        WHEN("We write the spikes ordered by id") {
            sonata_population->set_sorting("by_id");
            sonata_spikes->write_spikes_header(*sonata_population);
            THEN("We check that the spike node ids are in order") {
                const std::vector<uint64_t> node_ids = pop_p->get_spike_node_ids();
                std::vector<uint64_t> compare = {2, 2, 3, 3, 5};
                REQUIRE(node_ids == compare);
            }
            THEN("We check that the spike timestamps are ordered according to node ids") {
                const std::vector<double> timestamps = pop_p->get_spike_timestamps();
                std::vector<double> compare = {0.2, 0.7, 0.3, 1.3, 0.1};
                REQUIRE(timestamps == compare);
            }
        }
        WHEN("We dont order the spikes before writing") {
            sonata_population->set_sorting("none");
            sonata_spikes->write_spikes_header(*sonata_population);
            THEN("We check that the spike node ids are unordered") {
                const std::vector<uint64_t> node_ids = pop_p->get_spike_node_ids();
                std::vector<uint64_t> compare = {3, 5, 2, 3, 2};
                REQUIRE(node_ids == compare);
            }
            THEN("We check that the spike timestamps are unordered") {
                const std::vector<double> timestamps = pop_p->get_spike_timestamps();
                std::vector<double> compare = {0.3, 0.1, 0.2, 1.3, 0.7};
                REQUIRE(timestamps == compare);
            }
        }
        WHEN("We write the spikes ordered by weird string") {
            THEN("It throws an exception") {
                sonata_population->set_sorting("wrong_order");
                REQUIRE_THROWS(sonata_spikes->write_spikes_header(*sonata_population));
            }
        }
        sonata_spikes->close();
    }
}
