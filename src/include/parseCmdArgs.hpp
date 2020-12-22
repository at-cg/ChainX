#ifndef PARSE_CMD_HPP 
#define PARSE_CMD_HPP

#include "clipp/clipp.h"

//own includes
#include "utils.hpp"

namespace redit
{
  struct Parameters
  {
    std::string tfile;        //target sequence file (fasta/q)
    std::string qfile;        //file specifying query sequences
    int minLen;               //minimum MEM to consider
    std::string mode;         //"G" -> global, "SG" -> semi-global
                              //"DP_G" -> global using slow DP, "DP_SG" -> semi-global using slow DP
    std::string matchType;    //all MEMs or just consider MUMs (i.e., single occurence in query and ref)
    std::string indexpathl;   //file name from where index is loaded
    std::string indexpaths;   //file name where index is saved
    bool save = false;        //user opted to save index
    bool load = false;        //user opted to load index
    bool naive = false;       //use naive 2d dynamic programming algorithm similar to edit distance
  };

  void parseandSave(int argc, char** argv, Parameters &param)
  {
    //define all arguments
    auto cli =
      (
       clipp::required("-q") & clipp::value("path", param.qfile).doc("query sequences in fasta or fastq format"),
       clipp::required("-t") & clipp::value("path", param.tfile).doc("target sequence in fasta format"),
       clipp::required("-l") & clipp::value("length", param.minLen).doc("minimum match length (e.g., 20)"),
       clipp::required("-m") & (clipp::required("g").set(param.mode) | clipp::required("sg").set(param.mode)).doc("distance function (e.g., global or semi-global)"),
       clipp::required("-a") & (clipp::required("MEM").set(param.matchType) | clipp::required("MUM").set(param.matchType)).doc("filter maximal exact matches, 'MEM' consider all, 'MUM' must be single-copy"),
       clipp::option("--naive").set(param.naive).doc("use slow 2d dynamic programming algorithm for correctness check"),
       clipp::option("--load") & clipp::value("prefix1", param.indexpathl).doc("load suffix array from files starting with specified prefix"),
       clipp::option("--save") & clipp::value("prefix2", param.indexpaths).doc("save suffix array to files starting with specified prefix")
      );

    if(!clipp::parse(argc, argv, cli))
    {
      //print help page
      clipp::operator<<(std::cout, clipp::make_man_page(cli, argv[0])) << std::endl;
      exit(1);
    }

    param.load = param.indexpathl.size() > 0;
    param.save = param.indexpaths.size() > 0;

    if (param.load && param.save)
    {
      std::cerr << "ERROR, redit::parseandSave, loading and saving of index simultaneously is not allowed" << std::endl;
      exit(1);
    }

    //print all input parameters
    std::cerr << "INFO, redit::parseandSave, target sequence file = " << param.tfile << std::endl;
    std::cerr << "INFO, redit::parseandSave, query sequences file = " << param.qfile << std::endl;
    if (!param.naive) std::cerr << "INFO, redit::parseandSave, mode = " << param.mode << std::endl;
    if (param.naive) std::cerr << "INFO, redit::parseandSave, mode = " << param.mode << " (naive 2d DP)" << std::endl;
    std::cerr << "INFO, redit::parseandSave, anchor: min. length = " << param.minLen << ", type = " << param.matchType << std::endl;
    if (param.load) std::cerr << "INFO, redit::parseandSave, index will be loaded using prefix = " << param.indexpathl << std::endl;
    if (param.save) std::cerr << "INFO, redit::parseandSave, index will be saved using prefix = " << param.indexpaths << std::endl;

    if (! exists(param.tfile))
    {
      std::cerr << "ERROR, redit::parseandSave, target sequence file could not be opened" << std::endl;
      exit(1);
    }

    if (! exists(param.qfile))
    {
      std::cerr << "ERROR, redit::parseandSave, query sequence file could not be opened" << std::endl;
      exit(1);
    }
  }
}

#endif
