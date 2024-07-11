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

void Node::fill_data(std::vector<float>& report, size_t position, double step) {
    assert(elements_.empty() || element_handles_.empty());

    std::set<uint64_t> specific_node_ids = {2482561, 2537489, 2584861, 2644780, 2705687, 2761049, 2807052, 2854979, 2920895, 2974962, 3028685, 3074604, 3098646, 3116703};
    if (step == 1360.0) {
        if (specific_node_ids.find(node_id_) != specific_node_ids.end()) {
            for (const auto& elem : elements_) {
                if (std::abs(*elem) < 1e-3) {
                    std::cerr << "Aborting due to a pointer value being too close to zero." << std::endl;
                    abort();
                }
            }
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
