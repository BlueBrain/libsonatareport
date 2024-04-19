#include <algorithm>
#include <cassert>
#include <stdexcept>

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

void Node::update_elements(std::vector<uint32_t> element_ids, std::vector<double*> element_values) {
    if (element_ids.size() != element_values.size()) {
        throw std::runtime_error(
            "bbp::sonata::Node::update_elements: element_ids and element_values must have the same "
            "size");
    } else if (!elements_.empty()) {
        if (element_ids.size() != element_ids_.size()) {
            throw std::runtime_error(
                "bbp::sonata::Node::update_elements: element_ids and internal element_ids_ must "
                "have the same size");
        }
        for (size_t i = 0; i < element_ids.size(); ++i) {
            if (element_ids[i] != element_ids_[i]) {
                throw std::runtime_error(
                    "bbp::sonata::Node::update_elements: element_ids must be in the same order as "
                    "internal element_ids_");
            }
            elements_[i] = element_values[i];
        }
    }
}


void Node::fill_data(std::vector<float>::iterator it) {
    assert(elements_.empty() || element_handles_.empty());
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
