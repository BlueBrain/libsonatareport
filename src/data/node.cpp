#include <algorithm>
#include <cassert>
#include <stdexcept>
#include <iostream>

#include "../library/sonatareport.h"
#include "node.h"

namespace bbp {
namespace sonata {

Node::Node(uint64_t node_id)
    : node_id_(node_id) {}

void Node::add_element(double* element_value, uint32_t element_id) {
    if (!element_handles_.empty()) {
        throw std::runtime_error(
            "bbp::sonata::Node::add_element: mixing raw pointers and generic handles is not "
            "supported");
    }
    elements_.push_back(element_value);
    element_ids_.push_back(element_id);
}

void Node::add_element(std::function<double()> element_value, uint32_t element_id) {
    if (!elements_.empty()) {
        throw std::runtime_error(
            "bbp::sonata::Node::add_element: mixing raw pointers and generic handles is not "
            "supported");
    }
    element_handles_.push_back(std::move(element_value));
    element_ids_.push_back(element_id);
}

void Node::fill_data(std::vector<float>& report, size_t position) {
    assert(elements_.empty() || element_handles_.empty());

    {
        size_t num_elements_to_copy = 0;

        if (!elements_.empty()) {
            num_elements_to_copy = elements_.size();
        } else if (!element_handles_.empty()) {
            num_elements_to_copy = element_handles_.size();
        }

        size_t total_size_report = report.size();
        size_t start_position = position;

        if (start_position + num_elements_to_copy > total_size_report) {

            size_t overflow = (start_position + num_elements_to_copy) - total_size_report;

            std::cerr << "ERROR: Overflow detected in Node::fill_data()\n";
            std::cerr << "Rank: " << SonataReport::rank_ << "\n";
            std::cerr << "Report Allocated Size: " << total_size_report << "\n";
            std::cerr << "Elements to Write: " << num_elements_to_copy << "\n";
            std::cerr << "Overflow: " << overflow << "\n";

            abort();
        }
    }

    std::vector<float>::iterator it = report.begin() + position;

    if (!elements_.empty()) {
        std::transform(elements_.begin(), elements_.end(), it, [](auto elem) -> float {
            return *elem;
        });
    } else if (!element_handles_.empty()) {
        std::transform(element_handles_.begin(),
                       element_handles_.end(),
                       it,
                       [](auto const& elem) -> float { return elem(); });
    }
}

void Node::refresh_pointers(std::function<double*(double*)> refresh_function) {
    if (!elements_.empty()) {
        std::transform(elements_.begin(), elements_.end(), elements_.begin(), refresh_function);
    } else if (!element_handles_.empty()) {
        std::transform(element_handles_.begin(),
                       element_handles_.end(),
                       element_handles_.begin(),
                       [&refresh_function](auto const& elem) -> std::function<double()> {
                           return [elem, refresh_function]() -> double {
                               double value = elem();
                               double* refreshed_value = refresh_function(&value);
                               return *refreshed_value;
                           };
                       });
    }
}

}  // namespace sonata
}  // namespace bbp
