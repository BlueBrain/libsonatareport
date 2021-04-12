#pragma once

#include <fstream>
#include <hdf5.h>
#ifdef H5_HAVE_PARALLEL
#include <mpi.h>
#endif

namespace bbp {
namespace sonata {

/**
 * Enum that defines the type of backend filesystem from a path.
 */
typedef enum {
    FSTYPE_DEFAULT = 0x6f510ca1,   // Path on the BFS / local storage
    FSTYPE_IME = 0x13e00000,       // Path with "ime:" prefix
    FSTYPE_IME_MNTP = 0x13e00001,  // Path where IME is mounted on the BFS
    FSTYPE_IME_FUSE = 0x13e00002,  // https://man7.org/linux/man-pages/man2/fstatfs.2.html
    FSTYPE_UNKNOWN = 0xffffffff,   // Error (e.g., file not found)
    // FSTYPE_GPFS = 0x47504653,       //
    // https://www.ibm.com/support/knowledgecenter/STXKQY_4.2.2/com.ibm.spectrum.scale.v4r22.doc/bl1adm_admxafg.htm
    // FSTYPE_LUSTRE = 0x47504653,     //
    // https://wiki.old.lustre.org/index.php/Architecture_-_Wire_Level_Protocol#Lustre_Messages:_RPC.E2.80.99s
    // FSTYPE_NVME = 0xef53            // https://man7.org/linux/man-pages/man2/fstatfs.2.html
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
    static std::pair<fstype_t, std::string> getPathInfo(std::string path);

    /**
     * Converts a filesystem type to the equivalent string.
     * \param type Filesystem type.
     * \return Equivalent string to the type.
     */
    static std::string getFSTypeString(const fstype_t type);

#ifdef H5_HAVE_PARALLEL
    /**
     * Defines the MPI Hints necessary to use IME efficiently with MPI-IO.
     * \param info MPI Info object to be created / updated.
     * \return MPI_SUCCESS if successful, or an MPI-based error.
     */
    static int setMPIHints(MPI_Info& info);
#endif
  private:
    IMEUtil() {}
    ~IMEUtil() {}
};

}  // namespace sonata
}  // namespace bbp
