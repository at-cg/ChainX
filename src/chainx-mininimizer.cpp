#include <iostream>
#include <cassert>
#include <vector>
#include <tuple>
#include <algorithm>
#include <zlib.h>  
#include <string>
#include <chrono>
#include <stdlib.h>
#include <assert.h>

//third-party lib
#include "kseq/kseq.h"
#include "prettyprint/prettyprint.hpp"
#include "minimap2-2.24/minimap.h"

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

  /*********** MINIMAP2 API ****************/
  mm_idxopt_t iopt;
  mm_mapopt_t mopt;
  int n_threads = 1;

  mm_verbose = 2; // disable message output to stderr
  mm_set_opt(0, &iopt, &mopt);
  mopt.flag |= MM_F_CIGAR; // perform alignment
  mopt.flag |= MM_F_FOR_ONLY; //forward strand only

  // open query file for reading; you may use your favorite FASTA/Q parser
  gzFile f = gzopen(parameters.qfile.c_str(), "r");
  assert(f);
  kseq_t *ks = kseq_init(f);

  mm_idx_reader_t *r = mm_idx_reader_open(parameters.tfile.c_str(), &iopt, 0); //index target
  mm_idx_t *mi;
  while ((mi = mm_idx_reader_read(r, n_threads)) != 0) {
    mm_mapopt_update(&mopt, mi);
    mopt.mid_occ = 0; //do not ignore high frequency minimizers
    mm_tbuf_t *tbuf = mm_tbuf_init(); // thread buffer; for multi-threading, allocate one tbuf for each thread
    gzrewind(f);
    kseq_rewind(ks);
    while (kseq_read(ks) >= 0) { // each kseq_read() call reads one query sequence
      int j, i, n_reg;
      mm_reg1_t *reg;
      reg = mm_map(mi, ks->seq.l, ks->seq.s, &n_reg, tbuf, &mopt, 0); // get all hits for the query
      free(reg);
    }
    mm_tbuf_destroy(tbuf);
    mm_idx_destroy(mi);
  }
  mm_idx_reader_close(r); // close the index reader
  kseq_destroy(ks); // close the query file
  gzclose(f);

  return 0;
}
