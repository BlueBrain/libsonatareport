#include "imeutil.h"
#include <limits.h>
#include <sys/vfs.h>

#define IME_PREFIX "ime://"
#define IME_CONF_ENV "IM_CLIENT_CFG_FILE"
#define IME_CONF_PATH "/etc/ddn/ime/ime.conf"
#define IME_CONF_FUSE_PATH "/etc/ddn/ime/ime-fuse.conf"
#define FUSE_SUPER_MAGIC 0x65735546

using namespace bbp::sonata;

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
        return "";
    };

    return std::pair<std::string, std::string>(
        parseFile(std::ifstream(conf_path), "mount_point", ';'),
        parseFile(std::ifstream(IME_CONF_FUSE_PATH), "MNTDIR", '\''));
}

std::pair<fstype_t, std::string> IMEUtil::getPathInfo(std::string path) {
    // If the path begins with "ime:", assume native access
    if (path.find("ime:") == 0) {
        return std::pair<fstype_t, std::string>(FSTYPE_IME, path);
    }

    // Resolve the full path and return an error if not possible (e.g., file,
    // original parent folder, or both do not exist yet)
    const auto limit = path.find_last_of('/');
    auto path_orig = (limit != std::string::npos) ? path.substr(0, limit) : ".";
    char full_path[PATH_MAX];
    if (realpath(path_orig.c_str(), full_path) == NULL) {
        return std::pair<fstype_t, std::string>(FSTYPE_UNKNOWN, path);
    }
    path = std::string(full_path) + "/" + path.substr((limit != std::string::npos) ? limit : 0);

    // Check if the path contains the IME keyword
    if (path.find("ime") != std::string::npos) {
        static auto mnt_paths = getIMEMountPoints();

        // Check if the path contains the BFS mount point
        if (path.find(mnt_paths.first) == 0) {
            return std::pair<fstype_t, std::string>(FSTYPE_IME_MNTP, IME_PREFIX + path);
        }

        // Lastly, evaluate if the path is under a FUSE mount point
        struct statfs st;
        if (path.find(mnt_paths.second) == 0 && statfs(path_orig.c_str(), &st) == 0 &&
            st.f_type == FUSE_SUPER_MAGIC) {
            const off_t offset = mnt_paths.second.size();
            path = IME_PREFIX + mnt_paths.first + path.substr(offset);
            return std::pair<fstype_t, std::string>(FSTYPE_IME_FUSE, path);
        }
    }

    // At this point, assume a traditional file system
    return std::pair<fstype_t, std::string>(FSTYPE_DEFAULT, path);
}

std::string IMEUtil::getFSTypeString(const fstype_t type) {
    switch (type) {
    case FSTYPE_IME:
        return "IME (Native)";
    case FSTYPE_IME_MNTP:
        return "IME (BFS Mount Point)";
    case FSTYPE_IME_FUSE:
        return "IME (FUSE)";
    case FSTYPE_DEFAULT:
        return "BFS / Local Disk";
    default:
        return "Error";
    }
}

#ifdef H5_HAVE_PARALLEL
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
