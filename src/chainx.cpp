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
  chainx::parseandSave_chainx(argc, argv, parameters);

  std::vector<std::string> queries; //one or multiple sequences
  std::vector<std::string> query_ids;
  std::vector<std::string> target; //single sequence
  std::vector<std::string> target_ids;

  chainx::readSequences(parameters.qfile, queries, query_ids);
  chainx::readSequences(parameters.tfile, target, target_ids);

  int queryLenSum = 0;
  for (auto &q: queries) queryLenSum += q.length();
  std::cerr << "INFO, chainx::main, read " << queries.size() << " queries, " << queryLenSum << " residues\n";
  if (!parameters.all2all) std::cerr << "INFO, chainx::main, read target, " << target[0].length() << " residues\n";

  //Start timer
  auto tStart = std::chrono::system_clock::now();
  std::cerr << "\nINFO, chainx::main, timer set\n";

  std::vector<std::tuple<int, int, int>> fwd_matches;
  //lambda function
  auto append_matches = [&](const mummer::mummer::match_t& m) { fwd_matches.emplace_back(m.ref, m.query, m.len); }; //0-based coordinates

  if (!parameters.all2all)
  {
    //Compute anchors
    mummer::mummer::sparseSA sa (mummer::mummer::sparseSA::create_auto(target[0].data(), target[0].length(), parameters.minLen, true));

    std::chrono::duration<double> wctduration = (std::chrono::system_clock::now() - tStart);
    std::cerr << "INFO, chainx::main, suffix array computed in " << wctduration.count() << " seconds\n";

    for (int i = 0; i < queries.size(); i++)
    {
      std::cerr << "\nINFO, chainx::main, timer reset\n";
      tStart = std::chrono::system_clock::now();
      fwd_matches.clear();
      if (parameters.matchType == "MEM")
        sa.findMEM_each(queries[i].data(), queries[i].length(), parameters.minLen, false, append_matches);
      else if (parameters.matchType == "MUM")
        sa.findMUM_each(queries[i].data(), queries[i].length(), parameters.minLen, false, append_matches);
      else
        std::cerr << "ERROR, chainx::main, incorrect anchor type specified" << "\n";


      wctduration = (std::chrono::system_clock::now() - tStart);
      if (VERBOSE && parameters.matchType == "MEM") std::cerr << "INFO, chainx::main, MEMs identified (" << wctduration.count() << " seconds elapsed)\n";
      if (VERBOSE && parameters.matchType == "MUM") std::cerr << "INFO, chainx::main, MUMs identified (" << wctduration.count() << " seconds elapsed)\n";

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
      std::cerr << "INFO, chainx::main, count of anchors (including dummy) = " << fwd_matches.size() << ", average length = " << sum_anchor_len * 1.0 / fwd_matches.size() << "\n";

      if (VERBOSE)
        std::cerr << "List of sorted anchors = " << fwd_matches << "\n";

      //compute anchor-restricted edit distance
      std::cerr << "INFO, chainx::main, query #" << i << " (" << queries[i].length() << " residues), ";
      if (parameters.mode == "g")
      {
        if (parameters.naive)
          std::cout << "distance = " << chainx::DP_global(fwd_matches) << "\n";
        else
          std::cout << "distance = " << chainx::compute_global(fwd_matches) << "\n";
      }
      else if (parameters.mode == "sg")
      {
        if (parameters.naive)
          std::cout << "distance = " << chainx::DP_semiglobal(fwd_matches) << "\n";
        else
          std::cout << "distance = " << chainx::compute_semiglobal(fwd_matches) << "\n";
      }
      else
        std::cerr << "ERROR, chainx::main, incorrect mode specified" << "\n";


      wctduration = (std::chrono::system_clock::now() - tStart);
      std::cerr << "INFO, chainx::main, distance computation finished (" << wctduration.count() << " seconds elapsed)\n";
    }
  }
  else
  {
    std::vector<std::vector<int>> costs (queries.size());
    for(std::size_t i = 0; i < queries.size(); i++) costs[i] = std::vector<int>(queries.size(), -1);

    for (std::size_t i = 0; i < queries.size(); i++)
    {
      //build SA of queries[i]
      mummer::mummer::sparseSA sa = mummer::mummer::sparseSA::create_auto(queries[i].data(), queries[i].length(), parameters.minLen, true);

      for (std::size_t j = 0; j < i; j++)
      {
        //compute costs[i][j] && costs[j][i]

        fwd_matches.clear();
        if (parameters.matchType == "MEM")
          sa.findMEM_each(queries[j].data(), queries[j].length(), parameters.minLen, false, append_matches);
        else if (parameters.matchType == "MUM")
          sa.findMUM_each(queries[j].data(), queries[j].length(), parameters.minLen, false, append_matches);
        else
        {
          std::cerr << "ERROR, chainx::main, incorrect anchor type specified" << "\n";
          exit(1);
        }

        //place dummy MEMs and then sort
        fwd_matches.emplace_back(-1,-1,1);
        fwd_matches.emplace_back(queries[i].length(), queries[j].length(), 1);
        std::sort (fwd_matches.begin(), fwd_matches.end(),
            [](const std::tuple<int,int,int>& a,
              const std::tuple<int,int,int>& b) -> bool
            {
            return std::get<0>(a) < std::get<0>(b);
            });

        if (parameters.naive)
          costs[j][i] = costs[i][j] = chainx::DP_global(fwd_matches);
        else
          costs[j][i] = costs[i][j] = chainx::compute_global(fwd_matches);

      }

      costs[i][i] = 0;
    }

    std::cerr << "\nINFO, chainx::main, printing distance matrix to stdout\n";

    //phylip-formatted output
    {
      std::cout << queries.size() << "\n";
      for (std::size_t i = 0; i < queries.size(); i++)
      {
        std::cout << query_ids[i];
        for (std::size_t j = 0; j < queries.size(); j++)
        {
          std::cout << "  " << costs[i][j];
        }
        std::cout << "\n";
      }
    }

    std::chrono::duration<double> wctduration = (std::chrono::system_clock::now() - tStart);
    std::cerr << "INFO, chainx::main, all-to-all distance computation took " << wctduration.count() << " seconds\n";
  }

  return 0;
}
