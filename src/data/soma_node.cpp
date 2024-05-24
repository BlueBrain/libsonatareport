#include <stdexcept>

#include "soma_node.h"

namespace bbp {
namespace sonata {

SomaNode::SomaNode(uint64_t node_id)
    : Node(node_id) {}

void SomaNode::add_element(double* element_value, uint32_t element_id) {
    if (!elements_.empty() || !element_handles_.empty()) {
        throw std::runtime_error("ERROR: Soma report nodes can only have 1 element");
    }
    elements_.push_back(element_value);
    element_ids_.push_back(element_id);
}

void SomaNode::add_element(std::function<double()> element_value, uint32_t element_id) {
    if (!elements_.empty() || !element_handles_.empty()) {
        throw std::runtime_error("ERROR: Soma report nodes can only have 1 element");
    }
    element_handles_.push_back(std::move(element_value));
    element_ids_.push_back(element_id);
}

void SomaNode::update_elements(std::vector<uint32_t> element_ids,
                               std::vector<double*> element_values) {
    if (element_ids.size() != 1 || element_ids_[0] != element_ids[0]) {
        throw std::runtime_error(
            "ERROR: Soma report nodes can only have 1 element and IDs must match");
    }
    elements_ = std::move(element_values);
    element_ids_ = std::move(element_ids);
}
}  // namespace sonata
}  // namespace bbp
