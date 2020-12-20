#include <iostream>
#include <cassert>
#include <vector>
#include <tuple>
#include <algorithm>
#include <zlib.h>  
#include <string>
#include "mummer/sparseSA.hpp"
#include "kseq/kseq.h"
#include "edlib/edlib.h"
#include "prettyprint/prettyprint.hpp"

KSEQ_INIT(gzFile, gzread)
#define VERBOSE 0

/**
 * @brief   reads single sequence from input fasta / fastq file
 **/
void readSeq(const char* filename, std::string &str)
{
	gzFile fp = gzopen(filename, "r");
	kseq_t *seq = kseq_init(fp);
	int len;

	if ((len = kseq_read(seq)) >= 0) 
	{
	  str = seq->seq.s;
  }

	std::transform(str.begin(), str.end(), str.begin(), ::toupper);

	kseq_destroy(seq);  
	gzclose(fp); 
}

/**
 * @brief   compute anchor-restricted edit distance using strong precedence criteria
 **/
int redit_strong_prec(const std::vector<std::tuple<int, int, int>> &anchors)
{
	int n = anchors.size();
	std::vector<int> costs(n, 0);
	for(int j=1; j<n; j++)
	{
	  //compute cost[i] here
	  int find_min_cost = std::numeric_limits<int>::max();

	  int j_a = std::get<0>(anchors[j]);
	  int j_b = std::get<0>(anchors[j]) + std::get<2>(anchors[j]) - 1;
	  int j_c = std::get<1>(anchors[j]);
	  int j_d = std::get<1>(anchors[j]) + std::get<2>(anchors[j]) - 1;

	  // anchor i < anchor j 
	  for(int i=j-1; i>=0; i--)
	  {
		int i_a = std::get<0>(anchors[i]);
		int i_b = std::get<0>(anchors[i]) + std::get<2>(anchors[i]) - 1;
		int i_c = std::get<1>(anchors[i]);
		int i_d = std::get<1>(anchors[i]) + std::get<2>(anchors[i]) - 1;

		if (i_a < j_a && i_b < j_b && i_c < j_c && i_d < j_d)
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

	if (VERBOSE)
		std::cerr << "Cost array = " << costs << "\n";

	return costs[n-1];
}

/**
 * @brief   compute anchor-restricted edit distance using strong precedence criteria
 **/
int redit_strong_prec_optimized(const std::vector<std::tuple<int, int, int>> &anchors)
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

	std::cerr << "Chaining cost computed " << revisions + 1 << " times" << "\n";
	return costs[n-1];
}

/**
 * @brief   compute anchor-restricted edit distance using strong precedence criteria
 **/
int redit_strong_prec_infix_optimized(const std::vector<std::tuple<int, int, int>> &anchors)
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

					if (j==1 || j == n-1) gap1=0; //free terminal gap
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

	std::cerr << "Chaining cost computed " << revisions + 1 << " times" << "\n";
	return costs[n-1];
}

/**
 * @brief   compute anchor-restricted edit distance using weak precedence criteria (begin points only)
 **/
int redit_weak_prec(const std::vector<std::tuple<int, int, int>> &anchors)
{
	int n = anchors.size();
	std::vector<int> costs(n, 0);
	for(int j=1; j<n; j++)
	{
	  //compute cost[i] here
	  int find_min_cost = std::numeric_limits<int>::max();

	  int j_a = std::get<0>(anchors[j]);
	  int j_b = std::get<0>(anchors[j]) + std::get<2>(anchors[j]) - 1;
	  int j_c = std::get<1>(anchors[j]);
	  int j_d = std::get<1>(anchors[j]) + std::get<2>(anchors[j]) - 1;

	  // anchor i < anchor j 
	  for(int i=j-1; i>=0; i--)
	  {
		int i_a = std::get<0>(anchors[i]);
		int i_b = std::get<0>(anchors[i]) + std::get<2>(anchors[i]) - 1;
		int i_c = std::get<1>(anchors[i]);
		int i_d = std::get<1>(anchors[i]) + std::get<2>(anchors[i]) - 1;

		if (i_a < j_a && i_c < j_c)
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

	if (VERBOSE)
		std::cerr << "Cost array = " << costs << "\n";

	return costs[n-1];
}

/**
 * @brief   compute anchor-restricted edit distance using weak precedence criteria (end points only)
 **/
int redit_weak_rev_prec(const std::vector<std::tuple<int, int, int>> &anchors_)
{
	std::vector<std::tuple<int, int, int>> anchors = anchors_; //local copy

	std::sort (anchors.begin(), anchors.end(),
			[](const std::tuple<int,int,int>& a,
				const std::tuple<int,int,int>& b) -> bool
			{
			//sort by end point of each anchor rather than start
			return std::get<0>(a) + std::get<1>(a) < std::get<0>(b) + std::get<1>(b);
			});

	int n = anchors.size();
	std::vector<int> costs(n, 0);
	for(int j=1; j<n; j++)
	{
	  //compute cost[i] here
	  int find_min_cost = std::numeric_limits<int>::max();

	  int j_a = std::get<0>(anchors[j]);
	  int j_b = std::get<0>(anchors[j]) + std::get<2>(anchors[j]) - 1;
	  int j_c = std::get<1>(anchors[j]);
	  int j_d = std::get<1>(anchors[j]) + std::get<2>(anchors[j]) - 1;

	  // anchor i < anchor j 
	  for(int i=j-1; i>=0; i--)
	  {
		int i_a = std::get<0>(anchors[i]);
		int i_b = std::get<0>(anchors[i]) + std::get<2>(anchors[i]) - 1;
		int i_c = std::get<1>(anchors[i]);
		int i_d = std::get<1>(anchors[i]) + std::get<2>(anchors[i]) - 1;

		if (i_b < j_b && i_d < j_d)
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

	if (VERBOSE)
		std::cerr << "Cost array = " << costs << "\n";

	return costs[n-1];
}

/**
 * @brief   compute anchor-restricted edit distance using standard edit-distance like dynamic programming 
 **/
int redit_dp(const std::vector<std::tuple<int, int, int>> &anchors)
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

int main(int argc, char **argv) 
{
	assert (argc == 5);

	std::string seq_ref, seq_qry;
	readSeq(argv[1], seq_ref);
	readSeq(argv[2], seq_qry);
	int k = std::stoi(argv[3]);
	int method = std::stoi(argv[4]); //0=edit, 1=redit_s, 2=redit_w1; 3=redit_w2, 4=redit_dp, 5 = redit_s_optimized, 6 = redit_s_infix_optimized

	assert (seq_ref.length() > 0 && seq_qry.length() > 0);

	std::cerr << "Length of ref sequence = " << seq_ref.length() << "\n";
	std::cerr << "Length of qry sequence = " << seq_qry.length() << "\n";

	//compute exact edit distance using external library
	if (method == 0)
	{
	  EdlibAlignResult result = edlibAlign(seq_ref.c_str(), seq_ref.length(), seq_qry.c_str(), seq_qry.length(), edlibDefaultAlignConfig());
	  if (result.status == EDLIB_STATUS_OK) {
		  std::cout << result.editDistance << "\n";
	  }
	  edlibFreeAlignResult(result);
  }
	else if (method > 0)
	{
	  std::vector<std::tuple<int, int, int>> fwd_matches;
	  mummer::mummer::sparseSA sa = mummer::mummer::sparseSA::create_auto(seq_ref.c_str(), seq_ref.length(), k, true);
	  auto append_matches = [&](const mummer::mummer::match_t& m) { fwd_matches.emplace_back(m.ref, m.query, m.len); }; //0-based coordinates
	  sa.findMEM_each(seq_qry.c_str(), seq_qry.length(), k, false, append_matches);


	  //place dummy MEMs and then sort
	  fwd_matches.emplace_back(-1,-1,1);
	  fwd_matches.emplace_back(seq_ref.length(), seq_qry.length(), 1);
	  std::sort (fwd_matches.begin(), fwd_matches.end(),
			  [](const std::tuple<int,int,int>& a,
				  const std::tuple<int,int,int>& b) -> bool
			  {
			  return std::get<0>(a) < std::get<0>(b);
			  });

	  std::cerr << "Count of anchors including dummy = " << fwd_matches.size() << "\n";

	  if (VERBOSE)
		  std::cerr << "List of sorted anchors = " << fwd_matches << "\n";

	  //compute anchor-restricted edit distance
	  if (method == 1)
		  std::cout << redit_strong_prec(fwd_matches) << "\n";
	  else if (method == 2)
		  std::cout << redit_weak_prec(fwd_matches) << "\n";
	  else if (method == 3)
		  std::cout << redit_weak_rev_prec(fwd_matches) << "\n";
	  else if (method == 4)
		  std::cout << redit_dp(fwd_matches) << "\n";
	  else if (method == 5)
		  std::cout << redit_strong_prec_optimized(fwd_matches) << "\n";
	  else if (method == 6)
		  std::cout << redit_strong_prec_infix_optimized(fwd_matches) << "\n";
  }

	return 0;
}
