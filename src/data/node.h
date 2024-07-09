#pragma once
#include <cstdint>
#include <functional>
#include <map>
#include <memory>
#include <vector>

namespace bbp {
namespace sonata {

class Node
{
  public:
    explicit Node(uint64_t node_id);
    virtual ~Node() = default;
    void fill_data(std::vector<float>& report, size_t position, double step);
    void refresh_pointers(std::function<double*(double*)> refresh_function);
    virtual void add_element(double* element_value, uint32_t element_id);
    virtual void add_element(std::function<double()> element_value, uint32_t element_id);
    virtual void update_elements(std::vector<uint32_t> element_ids,
                                 std::vector<double*> element_values);

    uint64_t get_node_id() const noexcept {
        return node_id_;
    }
    virtual size_t get_num_elements() const noexcept {
        return element_ids_.size();
    }
    const std::vector<uint32_t>& get_element_ids() const noexcept {
        return element_ids_;
    }

  private:
    uint64_t node_id_;

  protected:
    std::vector<uint32_t> element_ids_;
    std::vector<double*> elements_;
    std::vector<std::function<double()>> element_handles_;
};

using nodes_t = std::map<uint64_t, std::shared_ptr<Node>>;

}  // namespace sonata
}  // namespace bbp
