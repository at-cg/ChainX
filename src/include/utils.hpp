#ifndef COMMON_UTILS_HPP
#define COMMON_UTILS_HPP

#include "kseq/kseq.h"
#include <fstream>
KSEQ_INIT(gzFile, gzread)

namespace chainx
{

  /**
   * @brief   reads sequences from input fasta / fastq file
   **/
  void readSequences(const std::string &path, std::vector<std::string> &seqs, std::vector<std::string> &ids)
  {
    gzFile fp = gzopen(path.data(), "r");

    if (fp == NULL) {
      fprintf(stderr, "gzopen failed to read input file\n");
      exit(1);
    }

    assert (seqs.size() == 0);

    kseq_t *seq = kseq_init(fp);
    int len;

    while ((len = kseq_read(seq)) >= 0) 
    {
      std::string str (seq->seq.s);
      std::string name (seq->name.s);

      //convert to upper case
      std::transform(str.begin(), str.end(), str.begin(), ::toupper);

      seqs.push_back(str);
      ids.push_back(name);
    }

    assert (seqs.size() > 0);

    kseq_destroy(seq);  
    gzclose(fp); 
  }

  inline bool exists (const std::string& filename) {
  std::ifstream f(filename.c_str());
  return f.good();
  }
}

#endif
