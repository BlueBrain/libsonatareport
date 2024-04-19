#include <catch2/catch.hpp>
#include <data/node.h>
#include <data/soma_node.h>
#include <memory>
#include <spdlog/spdlog.h>

using namespace bbp::sonata;

double* square(double* elem) {
    *elem *= *elem;
    return elem;
}

SCENARIO("Test Node class", "[Node]") {
    GIVEN("An instance of a Node") {
        Node node{1};
        {  // Initial conditions
            REQUIRE(node.get_node_id() == 1);
            REQUIRE(node.get_num_elements() == 0);
            REQUIRE(node.get_element_ids().empty());

            std::vector<float> res;
            node.fill_data(res.begin());
            REQUIRE(res.empty());

            REQUIRE_NOTHROW(node.refresh_pointers(&square));
        }

        WHEN("We add a raw pointer element") {
            std::vector<double> elements = {10, 11, 12, 13, 14};
            size_t i = 0;
            for (auto& element : elements) {
                node.add_element(&element, i);
                ++i;
            }
            THEN("Number of elements is 5") {
                REQUIRE(node.get_num_elements() == 5);
            }
            THEN("The element_ids are") {
                std::vector<uint32_t> compare = {0, 1, 2, 3, 4};
                REQUIRE(node.get_element_ids() == compare);
            }
            THEN("fill_data will return something correct") {
                std::vector<float> result(5, -1.);
                node.fill_data(result.begin());
                REQUIRE(result == std::vector<float>(elements.begin(), elements.end()));
            }
            THEN("refresh_pointers will be call on all elements") {
                node.refresh_pointers(&square);
                std::vector<float> compare{100, 121, 144, 169, 196};
                std::vector<float> result(5, -1);
                node.fill_data(result.begin());
                REQUIRE(result == compare);
            }
            THEN("update_elements with different sizes in elements should fail") {
                // Sizes of element_ids and element_values differ
                std::vector<uint32_t> element_ids = {0, 1, 2, 3, 4};
                std::vector<double> new_elements = {100, 111, 122};
                std::vector<double*> element_values;
                for (auto& element : new_elements) {
                    element_values.push_back(&element);
                }
                REQUIRE_THROWS(node.update_elements(element_ids, element_values));
            }
            THEN("update_elements with wrong size should fail") {
                // Initial node has 5 elements
                std::vector<uint32_t> element_ids = {0, 1, 2};
                std::vector<double> new_elements = {100, 111, 122};
                std::vector<double*> element_values;
                for (auto& element : new_elements) {
                    element_values.push_back(&element);
                }
                REQUIRE_THROWS(node.update_elements(element_ids, element_values));
            }
            THEN("update_elements with different element_ids should fail") {
                // Initial node has {0,1,2,3,4} element_ids
                std::vector<uint32_t> element_ids = {0, 1, 2, 3, 5};
                std::vector<double> new_elements = {100, 111, 122, 133, 144};
                std::vector<double*> element_values;
                for (auto& element : new_elements) {
                    element_values.push_back(&element);
                }
                REQUIRE_THROWS(node.update_elements(element_ids, element_values));
            }
            THEN("update_elements should update the pointers") {
                std::vector<uint32_t> element_ids = {0, 1, 2, 3, 4};
                std::vector<double> new_elements = {100, 111, 122, 133, 144};
                std::vector<double*> element_values;
                for (auto& element : new_elements) {
                    element_values.push_back(&element);
                }
                node.update_elements(element_ids, element_values);
                std::vector<float> result(5, -1.);
                node.fill_data(result.begin());
                REQUIRE(result == std::vector<float>(new_elements.begin(), new_elements.end()));
            }
        }

        WHEN("We add a std::function<double()> element") {
            std::vector<std::function<double()>> elements = {[]() { return 10.0; },
                                                             []() { return 11.0; },
                                                             []() { return 12.0; }};
            size_t i = 0;
            for (auto& element_handle : elements) {
                node.add_element(element_handle, i);
                ++i;
            }
            THEN("Number of elements is 3") {
                REQUIRE(node.get_num_elements() == 3);
            }
            THEN("The element_ids are") {
                std::vector<uint32_t> compare = {0, 1, 2};
                REQUIRE(node.get_element_ids() == compare);
            }
            THEN("fill_data will return something correct") {
                std::vector<float> result(3, -1.0);
                node.fill_data(result.begin());
                std::vector<float> compare = {10.0, 11.0, 12.0};
                REQUIRE(result == compare);
            }
            THEN("refresh_pointers will call the function on all elements") {
                node.refresh_pointers(&square);
                std::vector<float> compare{100, 121, 144};
                std::vector<float> result(3, -1);
                node.fill_data(result.begin());
                REQUIRE(result == compare);
            }
        }
    }

    GIVEN("An instance of a soma node") {
        SomaNode node{1};
        WHEN("We add one element") {
            double elem1 = 1;
            node.add_element(&elem1, 1);

            THEN("Number of elements is 1") {
                REQUIRE(node.get_num_elements() == 1);
            }

            THEN("Adding an other element throw") {
                REQUIRE_THROWS(node.add_element(&elem1, 2));
                REQUIRE(node.get_num_elements() == 1);
            }
        }
    }
}
