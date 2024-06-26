#pragma once
#include "node.h"

namespace bbp {
namespace sonata {

class SomaNode: public Node
{
  public:
    SomaNode(uint64_t node_id);

    void add_element(double* element_value, uint32_t element_id) override;
    void add_element(std::function<double()> element_value, uint32_t element_id) override;
    void update_elements(std::vector<uint32_t> element_ids,
                         std::vector<double*> element_values) override;
    size_t get_num_elements() const noexcept override {
        return elements_.empty() ? 0 : 1;
    };
};

}  // namespace sonata
}  // namespace bbp
