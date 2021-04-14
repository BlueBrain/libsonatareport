#include "imeutil.h"
#include <limits.h>
#ifdef SONATA_REPORT_CHECK_IME
#include <sys/vfs.h>
#endif

#define IME_PREFIX "ime://"
#define IME_CONF_ENV "IM_CLIENT_CFG_FILE"
#define IME_CONF_PATH "/etc/ddn/ime/ime.conf"
#define IME_CONF_FUSE_PATH "/etc/ddn/ime/ime-fuse.conf"
#define FUSE_SUPER_MAGIC 0x65735546  // https://man7.org/linux/man-pages/man2/fstatfs.2.html

using namespace bbp::sonata;

#ifdef SONATA_REPORT_CHECK_IME
/**
 * Parses the IME config. files to determine the BFS and FUSE mount points.
 */
std::pair<std::string, std::string> getIMEMountPoints() {
    const char* env = getenv(IME_CONF_ENV);
    const char* conf_path = (env != NULL) ? env : IME_CONF_PATH;

    const auto parseFile =
        [](std::ifstream ifs, const std::string keyword, const char sep) -> std::string {
        if (ifs.is_open()) {
            std::string line;
            while (std::getline(ifs, line)) {
                // Look for mount point setting, and ensure it is uncommented
                size_t offset = line.find(keyword);
                if (offset != std::string::npos && line.find('#') > offset) {
                    offset = line.find("/", offset);
                    return line.substr(offset, line.find(sep, offset) - offset);
                }
            }
        }
        return ":Error:-1:";
    };

    return std::pair<std::string, std::string>(
        parseFile(std::ifstream(conf_path), "mount_point", ';'),
        parseFile(std::ifstream(IME_CONF_FUSE_PATH), "MNTDIR", '\''));
}

/**
 * Verifies that a given path is under an active FUSE mount point.
 */
bool isFUSEMountPoint(const std::string& path) {
    struct statfs st;
    return (statfs(path.c_str(), &st) == 0 && st.f_type == FUSE_SUPER_MAGIC);
}
#endif

std::pair<fstype_t, std::string> IMEUtil::getPathInfo(std::string path) {
#ifdef SONATA_REPORT_CHECK_IME
    // If the path begins with "ime:", assume native access
    if (path.find("ime:") == 0) {
        return std::pair<fstype_t, std::string>(FSTYPE_IME, path);
    }

    // Resolve the full path and return an error if not possible (e.g., file,
    // original parent folder, or both do not exist yet)
    if (path[0] != '/') {
        const auto limit = path.find_last_of('/');
        auto path_orig = (limit != std::string::npos) ? path.substr(0, limit) : ".";
        char full_path[PATH_MAX];
        if (realpath(path_orig.c_str(), full_path) == NULL) {
            return std::pair<fstype_t, std::string>(FSTYPE_UNKNOWN, path);
        }
        path = std::string(full_path) + "/" + path.substr((limit != std::string::npos) ? limit : 0);
    }

    // Check if the path contains the IME keyword
    if (path.find("/ime/") != std::string::npos) {
        // Parse config. files and verify FUSE mount point only once for performance
        static auto mnt_paths = getIMEMountPoints();
        static bool ime_fuse_active = isFUSEMountPoint(mnt_paths.second);

        // Check if the path contains the BFS mount point
        if (path.find(mnt_paths.first) == 0) {
            return std::pair<fstype_t, std::string>(FSTYPE_IME, IME_PREFIX + path);
        }

        // Lastly, evaluate if the path is under a FUSE mount point
        if (ime_fuse_active && path.find(mnt_paths.second) == 0) {
            const off_t offset = mnt_paths.second.size();
            path = IME_PREFIX + mnt_paths.first + path.substr(offset);
            return std::pair<fstype_t, std::string>(FSTYPE_IME, path);
        }
    }
#endif
    // At this point, assume a traditional file system
    return std::pair<fstype_t, std::string>(FSTYPE_DEFAULT, path);
}

#ifdef SONATA_REPORT_HAVE_MPI
int IMEUtil::setMPIHints(MPI_Info& info) {
    int status = MPI_SUCCESS;

    // Create the MPI Info objects, if needed
    if (info == MPI_INFO_NULL) {
        status = MPI_Info_create(&info);
    }

    // Set the hints to disable two-phase I/O and data sieving
    if (status == MPI_SUCCESS) {
        status = MPI_Info_set(info, "romio_cb_write", "disable");
        status |= MPI_Info_set(info, "romio_ds_write", "disable");
    }

    return status;
}
#endif
