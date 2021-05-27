#pragma once

#include <fstream>
#ifdef SONATA_REPORT_HAVE_MPI
#include <mpi.h>
#endif

namespace bbp {
namespace sonata {

/**
 * Enum that defines the type of backend filesystem from a path.
 */
typedef enum {
    FSTYPE_DEFAULT = 0x6f510ca1,  // Path on the BFS / local storage
    FSTYPE_IME = 0x13e00000,      // Path with "ime:" prefix
    FSTYPE_UNKNOWN = 0xffffffff,  // Error (e.g., file not found)
} fstype_t;

/**
 * Helper class that defines utilities to make better use of IME.
 */
class IMEUtil
{
  public:
    /**
     * Determines the type of backend filesystem from a given path, and
     * provides the optimal IME path to use with libraries (e.g., MPI-IO).
     * \param path Path to the file or folder.
     * \return Filesystem type and optimal path for IME.
     */
    static std::pair<fstype_t, std::string> get_path_info(std::string path);

#ifdef SONATA_REPORT_HAVE_MPI
    /**
     * Defines the MPI Hints necessary to use IME efficiently with MPI-IO.
     * \param info MPI Info object to be created / updated.
     * \return MPI_SUCCESS if successful, or an MPI-based error.
     */
    static int set_mpi_hints(MPI_Info& info);
#endif
  private:
    IMEUtil() {}
    ~IMEUtil() {}
};

}  // namespace sonata
}  // namespace bbp
