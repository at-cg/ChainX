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

//own includes
#include "parseCmdArgs.hpp"
#include "algo.hpp"

#undef VERBOSE
#define VERBOSE 0

int main(int argc, char **argv) 
{
  chainx::Parameters parameters;
  chainx::parseandSave_printanchors(argc, argv, parameters);

  std::vector<std::string> queries; //one or multiple sequences
  std::vector<std::string> query_ids;
  std::vector<std::string> target; //single sequence
  std::vector<std::string> target_ids;

  chainx::readSequences(parameters.qfile, queries, query_ids);
  chainx::readSequences(parameters.tfile, target, target_ids);

  int queryLenSum = 0;
  for (auto &q: queries) queryLenSum += q.length();
  std::cerr << "INFO, printanchors::main, read " << queries.size() << " queries, " << queryLenSum << " residues\n";
  if (!parameters.all2all) std::cerr << "INFO, printanchors::main, read target, " << target[0].length() << " residues\n";

  //Start timer
  auto tStart = std::chrono::system_clock::now();
  std::cerr << "\nINFO, printanchors::main, timer set\n";

  std::vector<std::tuple<int, int, int>> fwd_matches;
  //lambda function
  auto append_matches = [&](const mummer::mummer::match_t& m) { fwd_matches.emplace_back(m.ref, m.query, m.len); }; //0-based coordinates

  //Compute anchors
  mummer::mummer::sparseSA sa (mummer::mummer::sparseSA::create_auto(target[0].data(), target[0].length(), parameters.minLen, true));

  std::chrono::duration<double> wctduration = (std::chrono::system_clock::now() - tStart);
  std::cerr << "INFO, printanchors::main, suffix array computed in " << wctduration.count() << " seconds\n";

  for (int i = 0; i < queries.size(); i++)
  {
    std::cerr << "\nINFO, printanchors::main, timer reset\n";
    tStart = std::chrono::system_clock::now();
    fwd_matches.clear();
    if (parameters.matchType == "MEM")
      sa.findMEM_each(queries[i].data(), queries[i].length(), parameters.minLen, false, append_matches);
    else if (parameters.matchType == "MUM")
      sa.findMUM_each(queries[i].data(), queries[i].length(), parameters.minLen, false, append_matches);
    else
      std::cerr << "ERROR, printanchors::main, incorrect anchor type specified" << "\n";

    wctduration = (std::chrono::system_clock::now() - tStart);
    std::cerr << "INFO, printanchors::main, anchor computation finished (" << wctduration.count() << " seconds elapsed)\n";
    std::cerr << "INFO, printanchors::main, printing anchor coorindates (1-based) to stdout\n";
    std::cerr << "INFO, printanchors::main, format: <qry_id> <qry_st> <qry_end> <target_st> <target_end> <length>\n";

    //print 1-based coordinates
    for (auto &e:fwd_matches)
    {
      std::cout << i << "\t" \
        << std::get<1>(e)+1 << "\t" \
        << std::get<1>(e) + std::get<2>(e) << "\t" \
        << std::get<0>(e)+1 << "\t" \
        << std::get<0>(e) + std::get<2>(e) << "\t" \
        << std::get<2>(e) << "\n";
    }
  }

  return 0;
}
