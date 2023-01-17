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

}  // namespace sonata
}  // namespace bbp
