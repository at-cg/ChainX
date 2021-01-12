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
    std::string mode;         //"g" -> global, "sg" -> semi-global
    std::string matchType;    //all MEMs or just consider MUMs (i.e., single occurence in query and ref)
    bool naive = false;       //use naive 2d dynamic programming algorithm similar to edit distance
    bool all2all = false;     //compute all to all global distance among query sequences
  };

  void parseandSave_redit(int argc, char** argv, Parameters &param)
  {
    //define all arguments
    auto cli =
      (
       clipp::required("-q") & clipp::value("path", param.qfile).doc("query sequences in fasta or fastq format"),
       clipp::required("-t") & clipp::value("path", param.tfile).doc("target sequence in fasta format"),
       clipp::required("-l") & clipp::value("length", param.minLen).doc("minimum match length (e.g., 20)"),
       clipp::required("-m") & (clipp::required("g").set(param.mode) | clipp::required("sg").set(param.mode)).doc("distance function (e.g., global or semi-global)"),
       clipp::required("-a") & (clipp::required("MEM").set(param.matchType) | clipp::required("MUM").set(param.matchType)).doc("filter maximal exact matches, 'MEM' consider all, 'MUM' must be single-copy"),
       clipp::option("--all2all").set(param.all2all).doc("output all to all global distances among query sequences in phylip format"),
       clipp::option("--naive").set(param.naive).doc("use slow 2d dynamic programming algorithm for correctness check")
      );

    if(!clipp::parse(argc, argv, cli))
    {
      //print help page
      clipp::operator<<(std::cout, clipp::make_man_page(cli, argv[0])) << std::endl;
      exit(1);
    }

    //print all input parameters
    std::cerr << "INFO, redit::parseandSave, target sequence file = " << param.tfile << std::endl;
    std::cerr << "INFO, redit::parseandSave, query sequences file = " << param.qfile << std::endl;
    if (!param.naive) std::cerr << "INFO, redit::parseandSave, mode = " << param.mode << std::endl;
    if (param.naive) std::cerr << "INFO, redit::parseandSave, mode = " << param.mode << " (naive 2d DP)" << std::endl;
    if (param.all2all) std::cerr << "INFO, redit::parseandSave, computing all-to-all distances" << std::endl;
    std::cerr << "INFO, redit::parseandSave, anchor : minimim length = " << param.minLen << ", type = " << param.matchType << std::endl;
    //if (param.load) std::cerr << "INFO, redit::parseandSave, index will be loaded using prefix = " << param.indexpathl << std::endl;
    //if (param.save) std::cerr << "INFO, redit::parseandSave, index will be saved using prefix = " << param.indexpaths << std::endl;

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

    if (param.all2all)
    {
      if (param.qfile != param.tfile)
      {
        std::cerr << "ERROR, redit::parseandSave, query and target sequence file paths must be same in all-to-all mode" << std::endl;
        exit(1);
      }

      if (param.mode != "g")
      {
        std::cerr << "ERROR, redit::parseandSave, only global distance function [ -m g ] can be used in all-to-all mode" << std::endl;
        exit(1);
      }
    }
  }

  void parseandSave_edlib(int argc, char** argv, Parameters &param)
  {
    //define all arguments
    auto cli =
      (
       clipp::required("-q") & clipp::value("path", param.qfile).doc("query sequences in fasta or fastq format"),
       clipp::required("-t") & clipp::value("path", param.tfile).doc("target sequence in fasta format"),
       clipp::required("-m") & (clipp::required("g").set(param.mode) | clipp::required("sg").set(param.mode)).doc("distance function (e.g., global or semi-global)"),
       clipp::option("--all2all").set(param.all2all).doc("output all to all global distances among query sequences in phylip format")
      );

    if(!clipp::parse(argc, argv, cli))
    {
      //print help page
      clipp::operator<<(std::cout, clipp::make_man_page(cli, argv[0])) << std::endl;
      exit(1);
    }

    //print all input parameters
    std::cerr << "INFO, redit::parseandSave, target sequence file = " << param.tfile << std::endl;
    std::cerr << "INFO, redit::parseandSave, query sequences file = " << param.qfile << std::endl;
    std::cerr << "INFO, redit::parseandSave, mode = " << param.mode << std::endl;
    if (param.all2all) std::cerr << "INFO, redit::parseandSave, computing all-to-all distances" << std::endl;

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

    if (param.all2all)
    {
      if (param.qfile != param.tfile)
      {
        std::cerr << "ERROR, redit::parseandSave, query and target sequence file paths must be same in all-to-all mode" << std::endl;
        exit(1);
      }

      if (param.mode != "g")
      {
        std::cerr << "ERROR, redit::parseandSave, only global distance function [ -m g ] can be used in all-to-all mode" << std::endl;
        exit(1);
      }
    }
  }
}

#endif
