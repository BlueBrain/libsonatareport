#include <vector>
#include <fstream>

#include "io_writer.hpp"

class AsciiWriter: public IoWriter {
  private:
    std::ofstream m_file;

  public:
    AsciiWriter(const std::string& report_name);
    ~AsciiWriter() = default;

    void write(double* buffer, int steps_to_write, int total_steps, int total_compartments) override;
    void write(const std::string& name, const std::vector<int>& buffer) override;
    void write(const std::string& name, const std::vector<uint32_t>& buffer) override;
    void write(const std::string& name, const std::vector<uint64_t>& buffer) override;
    void write(const std::string& name, const std::vector<float>& buffer) override;
    void write(const std::string& name, const std::vector<double>& buffer) override;

    void close() override;
};