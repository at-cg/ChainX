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

  //Start timer
  auto tStart = std::chrono::system_clock::now();
  std::cerr << "\nINFO, redit::main, timer set\n";

  //Compute anchors
  mummer::mummer::sparseSA sa = mummer::mummer::sparseSA::create_auto(target[0].data(), target[0].length(), parameters.minLen, true);

  std::chrono::duration<double> wctduration = (std::chrono::system_clock::now() - tStart);

  if (parameters.load)
    std::cerr << "INFO, redit::main, suffix array loaded  in " << wctduration.count() << " seconds\n";
  else
    std::cerr << "INFO, redit::main, suffix array computed in " << wctduration.count() << " seconds\n";

  for (int i = 0; i < queries.size(); i++) 
  {
    std::cerr << "\nINFO, redit::main, timer reset\n";
    tStart = std::chrono::system_clock::now();
    std::vector<std::tuple<int, int, int>> fwd_matches;
    auto append_matches = [&](const mummer::mummer::match_t& m) { fwd_matches.emplace_back(m.ref, m.query, m.len); }; //0-based coordinates
    sa.findMEM_each(queries[i].data(), queries[i].length(), parameters.minLen, false, append_matches);

    wctduration = (std::chrono::system_clock::now() - tStart);
    std::cerr << "INFO, redit::main, MEMs identified (" << wctduration.count() << " seconds elapsed)\n";

    //place dummy MEMs and then sort
    fwd_matches.emplace_back(-1,-1,1);
    fwd_matches.emplace_back(target[0].length(), queries[i].length(), 1);
    std::sort (fwd_matches.begin(), fwd_matches.end(),
        [](const std::tuple<int,int,int>& a,
          const std::tuple<int,int,int>& b) -> bool
        {
        return std::get<0>(a) < std::get<0>(b);
        });

    std::size_t sum_anchor_len = 0;
    for (auto &e: fwd_matches) sum_anchor_len += std::get<2>(e);
    std::cerr << "INFO, redit::main, count of anchors (including dummy) = " << fwd_matches.size() << ", average length = " << sum_anchor_len * 1.0 / fwd_matches.size() << "\n";

    if (VERBOSE)
      std::cerr << "List of sorted anchors = " << fwd_matches << "\n";


    //compute anchor-restricted edit distance
    std::cout << "INFO, redit::main, query #" << i << " (" << queries[i].length() << " residues), ";
    if (parameters.mode == "g")
      std::cout << "distance = " << redit::compute_global(fwd_matches) << "\n";
    else if (parameters.mode == "sg")
      std::cout << "distance = " << redit::compute_semiglobal(fwd_matches) << "\n";
    else if (parameters.mode == "dp-g")
      std::cout << "distance = " << redit::DP_global(fwd_matches) << "\n";
    else if (parameters.mode == "dp-sg")
      std::cout << "distance = " << redit::DP_semiglobal(fwd_matches) << "\n";
    else 
      std::cerr << "ERROR, redit::main, incorrect mode specified" << "\n";


    wctduration = (std::chrono::system_clock::now() - tStart);
    std::cerr << "INFO, redit::main, distance computation finished (" << wctduration.count() << " seconds elapsed)\n";
  }

	return 0;
}
