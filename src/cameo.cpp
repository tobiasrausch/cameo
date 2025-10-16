#define _SECURE_SCL 0
#define _SCL_SECURE_NO_WARNINGS
#include <iostream>
#include <fstream>

#define BOOST_DISABLE_ASSERTS
#define BOOST_UUID_RANDOM_PROVIDER_FORCE_POSIX

#ifdef PROFILE
#include "gperftools/profiler.h"
#endif

#include <boost/filesystem.hpp>
#include <htslib/sam.h>

#include "version.h"
#include "threadpool.h"
#include "util.h"
#include "cameo.h"

using namespace cameo;


inline void
displayUsage() {
  std::cerr << "Usage: cameo <command> <arguments>" << std::endl;
  std::cerr << std::endl;
  std::cerr << "    scan        find demethylated promoter" << std::endl;
  std::cerr << std::endl;
}

int main(int argc, char **argv) {
    if (argc < 2) { 
      printTitle("cameo");
      displayUsage();
      return 0;
    }

    if ((std::string(argv[1]) == "version") || (std::string(argv[1]) == "--version") || (std::string(argv[1]) == "--version-only") || (std::string(argv[1]) == "-v")) {
      std::cerr << "cameo version: v" << cameoVersionNumber << std::endl;
      std::cerr << " using Boost: v" << BOOST_VERSION / 100000 << "." << BOOST_VERSION / 100 % 1000 << "." << BOOST_VERSION % 100 << std::endl;
      std::cerr << " using HTSlib: v" << hts_version() << std::endl;
      return 0;
    }
    else if ((std::string(argv[1]) == "help") || (std::string(argv[1]) == "--help") || (std::string(argv[1]) == "-h") || (std::string(argv[1]) == "-?")) {
      printTitle("cameo");
      displayUsage();
      return 0;
    }
    else if ((std::string(argv[1]) == "warranty") || (std::string(argv[1]) == "--warranty") || (std::string(argv[1]) == "-w")) {
      displayWarranty();
      return 0;
    }
    else if ((std::string(argv[1]) == "license") || (std::string(argv[1]) == "--license") || (std::string(argv[1]) == "-l")) {
      bsd();
      return 0;
    }
    else if ((std::string(argv[1]) == "scan")) {
      return scan(argc-1,argv+1);
    }
    std::cerr << "Unrecognized command " << std::string(argv[1]) << std::endl;
    return 1;
}
