#define CATCH_CONFIG_RUNNER
#include <catch2/catch.hpp>
#include <hdf5.h>

#ifdef H5_HAVE_PARALLEL
#include <mpi.h>
#endif

int main(int argc, char* argv[]) {
#ifdef H5_HAVE_PARALLEL
    MPI_Init(nullptr, nullptr);
#endif
    int result = Catch::Session().run(argc, argv);
#ifdef H5_HAVE_PARALLEL
    MPI_Finalize();
#endif
    return result;
}
