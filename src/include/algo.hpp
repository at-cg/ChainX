#ifndef REDIT_ALGO_H
#define REDIT_ALGO_H

#include <iostream>
#include <cassert>
#include <vector>
#include <tuple>
#include <algorithm>
#include <zlib.h>  
#include <string>
#include <chrono>

#undef VERBOSE
#define VERBOSE 0

namespace redit
{
  /**
   * @brief   compute anchor-restricted edit distance using strong precedence criteria
   * 			optimized to run faster using engineering trick(s)
   **/
  int compute_global(const std::vector<std::tuple<int, int, int>> &anchors)
  {
    int n = anchors.size();
    std::vector<int> costs(n, 0);

    int bound_redit = 100; //distance assumed to be <= 100
    int revisions = 0;
    //with this assumption on upper bound of distance, a gap of >bound_redit will not be allowed between adjacent anchors

    while (true) 
    {
      int inner_loop_start = 0;

      for(int j=1; j<n; j++)
      {
        //compute cost[i] here
        int find_min_cost = std::numeric_limits<int>::max();

        int j_a = std::get<0>(anchors[j]);
        int j_b = std::get<0>(anchors[j]) + std::get<2>(anchors[j]) - 1;
        int j_c = std::get<1>(anchors[j]);
        int j_d = std::get<1>(anchors[j]) + std::get<2>(anchors[j]) - 1;

        // anchor i < anchor j 

        while (j_a - std::get<0>(anchors[inner_loop_start]) - 1 > bound_redit)
          inner_loop_start++;

        for(int i=j-1; i>=inner_loop_start; i--)
        {
          int i_a = std::get<0>(anchors[i]);
          int i_b = std::get<0>(anchors[i]) + std::get<2>(anchors[i]) - 1;
          int i_c = std::get<1>(anchors[i]);
          int i_d = std::get<1>(anchors[i]) + std::get<2>(anchors[i]) - 1;

          if (costs[i] < std::numeric_limits<int>::max() && i_a < j_a && i_b < j_b && i_c < j_c && i_d < j_d)
          {
            int gap1 = std::max(0, j_a - i_b - 1);
            int gap2 = std::max(0, j_c - i_d - 1);
            int g = std::max(gap1,gap2);

            int overlap1 = std::max(0, i_b - j_a + 1);
            int overlap2 = std::max(0, i_d - j_c + 1);
            int o = std::abs(overlap1 - overlap2);

            find_min_cost = std::min(find_min_cost, costs[i] + g + o);
          }
        }
        //save optimal cost at offset j
        costs[j] = find_min_cost;
      }

      if (costs[n-1] > bound_redit)
      {
        bound_redit = bound_redit * 4;
        revisions++;
      }
      else
        break;
    }

    if (VERBOSE)
      std::cerr << "Cost array = " << costs << "\n";

    if (VERBOSE)
      std::cerr << "Chaining cost computed " << revisions + 1 << " times" << "\n";
    return costs[n-1];
  }

  /**
   * @brief   compute anchor-restricted (semi-global) edit distance using strong precedence criteria
   * 			optimized to run faster using engineering trick(s)
   **/
  int compute_semiglobal(const std::vector<std::tuple<int, int, int>> &anchors)
  {
    int n = anchors.size();
    std::vector<int> costs(n, 0);

    int bound_redit = 100; //distance assumed to be <= 100
    int revisions = 0;
    //with this assumption on upper bound of distance, a gap of >bound_redit will not be allowed between adjacent anchors

    while (true) 
    {
      int inner_loop_start = 0;

      for(int j=1; j<n; j++)
      {
        //compute cost[i] here
        int find_min_cost = std::numeric_limits<int>::max();

        int j_a = std::get<0>(anchors[j]);
        int j_b = std::get<0>(anchors[j]) + std::get<2>(anchors[j]) - 1;
        int j_c = std::get<1>(anchors[j]);
        int j_d = std::get<1>(anchors[j]) + std::get<2>(anchors[j]) - 1;

        // anchor i < anchor j 

        while (j_a - std::get<0>(anchors[inner_loop_start]) - 1 > bound_redit)
          inner_loop_start++;

        {
          //always consider the first dummy anchor 
          //connection to first dummy anchor is done with modified cost to allow free gaps
          int i_d = std::get<1>(anchors[0]) + std::get<2>(anchors[0]) - 1;
          int qry_gap = j_c - i_d - 1;
          find_min_cost = std::min(find_min_cost, costs[0] + qry_gap);
        }

        //process all anchors in array for the final last dummy anchor
        if (j == n-1) inner_loop_start=0;

        for(int i=j-1; i>=inner_loop_start; i--)
        {
          int i_a = std::get<0>(anchors[i]);
          int i_b = std::get<0>(anchors[i]) + std::get<2>(anchors[i]) - 1;
          int i_c = std::get<1>(anchors[i]);
          int i_d = std::get<1>(anchors[i]) + std::get<2>(anchors[i]) - 1;

          if (costs[i] < std::numeric_limits<int>::max() && i_a < j_a && i_b < j_b && i_c < j_c && i_d < j_d)
          {
            int gap1 = std::max(0, j_a - i_b - 1);
            int gap2 = std::max(0, j_c - i_d - 1);

            if (j == n-1) gap1=0; //modified cost for the last dummy anchor to allow free gaps
            int g = std::max(gap1,gap2);

            int overlap1 = std::max(0, i_b - j_a + 1);
            int overlap2 = std::max(0, i_d - j_c + 1);
            int o = std::abs(overlap1 - overlap2);

            find_min_cost = std::min(find_min_cost, costs[i] + g + o);
          }
        }

        //save optimal cost at offset j
        costs[j] = find_min_cost;
      }

      if (costs[n-1] > bound_redit)
      {
        bound_redit = bound_redit * 4;
        revisions++;
      }
      else
        break;
    }

    if (VERBOSE)
      std::cerr << "Cost array = " << costs << "\n";

    if (VERBOSE)
      std::cerr << "Chaining cost computed " << revisions + 1 << " times" << "\n";
    return costs[n-1];
  }

  /**
   * @brief   compute anchor-restricted edit distance using standard edit-distance like dynamic programming 
   **/
  int DP_global(const std::vector<std::tuple<int, int, int>> &anchors)
  {
    int n = anchors.size();

    //get sequence lengths from end dummy anchor
    //assuming anchors are already sorted
    int len_ref = std::get<0>(anchors[n-1]);
    int len_qry = std::get<1>(anchors[n-1]);

    //initialize a boolean matrix (len_ref+1 x len_qry+1)
    //we will offset by 1 to be consistent with DP matrix
    std::vector<std::vector<bool> > matchAllowed(len_ref+1);
    for(int i=0; i<len_ref+1; i++) matchAllowed[i] = std::vector<bool>(len_qry+1, false);
    for(int i = 0; i<n-1; i++) //use all (except end dummy) anchors
    {
      int e_a = std::get<0>(anchors[i]);
      int e_c = std::get<1>(anchors[i]);
      int e_len = std::get<2>(anchors[i]);

      for(int j=0; j<e_len; j++)
      {
        matchAllowed[e_a+1  + j][e_c+1  + j] = true;
      }
    }


    //initialize dp_matrix (len_ref+1 x len_qry+1) 
    std::vector<std::vector<int> > dp_matrix(len_ref+1);
    for(int i=0; i<=len_ref; i++) dp_matrix[i] = std::vector<int>(len_qry+1);
    for(int i=0; i<=len_ref; i++) dp_matrix[i][0] = i;
    for(int j=0; j<=len_qry; j++) dp_matrix[0][j] = j;

    for(int i=1; i<=len_ref; i++)
    {
      for(int j=1; j<=len_qry; j++)
      {
        if (matchAllowed[i][j])
          dp_matrix[i][j] = std::min({dp_matrix[i - 1][j - 1], dp_matrix[i - 1][j] + 1, dp_matrix[i][j - 1] + 1});
        else
          dp_matrix[i][j] = std::min({dp_matrix[i - 1][j - 1] + 1, dp_matrix[i - 1][j] + 1, dp_matrix[i][j - 1] + 1});
      }
    }

    return dp_matrix[len_ref][len_qry];
  }

  /**
   * @brief   compute anchor-restricted (semi-global) edit distance using standard edit-distance like dynamic programming 
   **/
  int DP_semiglobal(const std::vector<std::tuple<int, int, int>> &anchors)
  {
    int n = anchors.size();

    //get sequence lengths from end dummy anchor
    //assuming anchors are already sorted
    int len_ref = std::get<0>(anchors[n-1]);
    int len_qry = std::get<1>(anchors[n-1]);

    //initialize a boolean matrix (len_ref+1 x len_qry+1)
    //we will offset by 1 to be consistent with DP matrix
    std::vector<std::vector<bool> > matchAllowed(len_ref+1);
    for(int i=0; i<len_ref+1; i++) matchAllowed[i] = std::vector<bool>(len_qry+1, false);
    for(int i = 0; i<n-1; i++) //use all (except end dummy) anchors
    {
      int e_a = std::get<0>(anchors[i]);
      int e_c = std::get<1>(anchors[i]);
      int e_len = std::get<2>(anchors[i]);

      for(int j=0; j<e_len; j++)
      {
        matchAllowed[e_a+1  + j][e_c+1  + j] = true;
      }
    }


    //initialize dp_matrix (len_ref+1 x len_qry+1) 
    std::vector<std::vector<int> > dp_matrix(len_ref+1);
    for(int i=0; i<=len_ref; i++) dp_matrix[i] = std::vector<int>(len_qry+1);
    for(int i=0; i<=len_ref; i++) dp_matrix[i][0] = 0; //changed from i to 0 for free gaps
    for(int j=0; j<=len_qry; j++) dp_matrix[0][j] = j;

    for(int i=1; i<=len_ref; i++)
    {
      for(int j=1; j<=len_qry; j++)
      {
        if (matchAllowed[i][j])
          dp_matrix[i][j] = std::min({dp_matrix[i - 1][j - 1], dp_matrix[i - 1][j] + 1, dp_matrix[i][j - 1] + 1});
        else
          dp_matrix[i][j] = std::min({dp_matrix[i - 1][j - 1] + 1, dp_matrix[i - 1][j] + 1, dp_matrix[i][j - 1] + 1});
      }
    }

    int final_distance = std::numeric_limits<int>::max();
    for(int i=0; i<=len_ref; i++) final_distance = std::min (final_distance, dp_matrix[i][len_qry]);
    return final_distance;
  }
}

#endif
