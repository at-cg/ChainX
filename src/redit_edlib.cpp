#include <iostream>
#include <cassert>
#include <vector>
#include <tuple>
#include <algorithm>
#include <zlib.h>  
#include <string>
#include <chrono>

//third-party lib
#include "mummer/sparseSA.hpp"
#include "kseq/kseq.h"
#include "prettyprint/prettyprint.hpp"
#include "edlib/edlib.h"

//own includes
#include "parseCmdArgs.hpp"
#include "algo.hpp"

#undef VERBOSE
#define VERBOSE 0

int main(int argc, char **argv) 
{
  redit::Parameters parameters;
  redit::parseandSave(argc, argv, parameters);

  std::vector<std::string> queries; //one or multiple sequences
  std::vector<std::string> target; //single sequence

  redit::readSequences(parameters.qfile, queries);
  redit::readSequences(parameters.tfile, target);

  int queryLenSum = 0;
  for (auto &q: queries) queryLenSum += q.length();
  std::cerr << "INFO, redit::main, read " << queries.size() << " queries, " << queryLenSum << " residues\n";
  std::cerr << "INFO, redit::main, read target, " << target[0].length() << " residues\n";

  for (int i = 0; i < queries.size(); i++) 
  {
    //Start timer
    std::cerr << "\nINFO, redit::main, timer reset\n";
    auto tStart = std::chrono::system_clock::now();

    //compute anchor-restricted edit distance
    std::cout << "INFO, redit::main, query #" << i << " (" << queries[i].length() << " residues), ";
    if (parameters.mode == "g")
    {
      EdlibAlignResult result = edlibAlign(queries[i].data(), queries[i].length(), target[0].data(), target[0].length(), edlibDefaultAlignConfig());
      if (result.status == EDLIB_STATUS_OK) {
        std::cout << "distance = " << result.editDistance << "\n";
      }
      edlibFreeAlignResult(result);
    }
    else if (parameters.mode == "sg")
    {
      EdlibAlignResult result = edlibAlign(queries[i].data(), queries[i].length(), target[0].data(), target[0].length(), 
          edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0));
      if (result.status == EDLIB_STATUS_OK) {
        std::cout << "distance = " << result.editDistance << "\n";
      }
    }
    else 
      std::cerr << "ERROR, redit::main, incorrect mode specified" << "\n";


    std::chrono::duration<double> wctduration = (std::chrono::system_clock::now() - tStart);
    std::cerr << "INFO, redit::main, distance computation finished (" << wctduration.count() << " seconds elapsed)\n";
  }

  return 0;
}
