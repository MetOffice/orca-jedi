/*
 * (C) British Crown Copyright 2023 Met Office
 */

#pragma once

#include <algorithm>
#include <chrono>
#include <exception>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <thread>

#include "eckit/config/LibEcKit.h"
#include "eckit/config/Resource.h"
#include "eckit/eckit.h"
#include "eckit/mpi/Comm.h"
#include "eckit/runtime/Main.h"
#include "eckit/testing/Test.h"

#include "atlas/library/Library.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"
#include "atlas/runtime/trace/StopWatch.h"
#include "atlas/util/Config.h"

namespace orcamodel {
namespace test {

void setEnv(const std::string& env, bool value) {
  constexpr int DO_NOT_REPLACE_IF_EXISTS = 0;
  ::setenv(env.c_str(), eckit::Translator<bool, std::string>()(value).c_str(),
      DO_NOT_REPLACE_IF_EXISTS);
}

struct OrcaModelTestEnvironment {
    using Config = atlas::util::Config;

    OrcaModelTestEnvironment(int argc, char* argv[]) {
        eckit::Main::initialise(argc, argv);
        eckit::Main::instance().taskID(eckit::mpi::comm("world").rank());


        setEnv("ATLAS_FPE", true);
        setEnv("ATLAS_SIGNAL_HANDLER", true);

        atlas::initialize();
        eckit::mpi::comm().barrier();
    }
    ~OrcaModelTestEnvironment() {
        atlas::finalize();
        atlas::mpi::finalize();
    }
};

static double ATLAS_MPI_BARRIER_TIMEOUT() {
    static double v = eckit::Resource<double>("$ATLAS_MPI_BARRIER_TIMEOUT", 3.);
    return v;
}

static int barrier_timeout(double seconds) {
    auto barrier = eckit::mpi::comm().iBarrier();
    atlas::runtime::trace::StopWatch watch;
    while (!barrier.test()) {
      watch.start();
      std::this_thread::sleep_for(std::chrono::milliseconds(100));
      watch.stop();
      if (watch.elapsed() > seconds) {
        return 1;
      }
    }
    return 0;
}

//------------------------------------------------------------------------------

template <typename Environment>
int run(int argc, char* argv[]) {
    Environment env(argc, argv);
    int errors = eckit::testing::run_tests(argc, argv, false);
    if (eckit::mpi::comm().size() > 1) {
        if (barrier_timeout(ATLAS_MPI_BARRIER_TIMEOUT())) {
            eckit::Log::warning() << "\nWARNING: Tests failed with MPI deadlock"
                                  << " (${ATLAS_MPI_BARRIER_TIMEOUT}="
                                  << ATLAS_MPI_BARRIER_TIMEOUT()
                                  << ").\nCalling MPI_Abort..." << std::endl;
            eckit::mpi::comm().abort();
        }
        eckit::mpi::comm().allReduceInPlace(errors, eckit::mpi::max());
    }
    return errors;
}

int run(int argc, char* argv[]) {
    return run<orcamodel::test::OrcaModelTestEnvironment>(argc, argv);
}

//------------------------------------------------------------------------------

}  // namespace test
}  // namespace orcamodel

