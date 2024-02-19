#include "arguments.hpp"

#include <seqan/arg_parse.h>

using namespace seqan;

bool getArguments(Arguments &a, int argc, char const ** argv) {

  ArgumentParser parser("countT2C");
  addArgument(parser, seqan::ArgParseArgument(
		  seqan::ArgParseArgument::STRING, "bam file name"));
  addArgument(parser, seqan::ArgParseArgument(
                  seqan::ArgParseArgument::STRING, "reference fasta file name"));
  addOption(parser, seqan::ArgParseOption(
		  "e", "exclude-file", "A file with positions to exclude. Tab-delimited,"
		    "first column contig name, second column 1-based positon).",
		    seqan::ArgParseArgument::STRING, "STRING"));
  addOption(parser, seqan::ArgParseOption(
		  "l", "readlength", "Read length (or maximal number of Ts per read).",
		  seqan::ArgParseArgument::INTEGER, "INT"));
  addOption(parser, seqan::ArgParseOption(
		  "o", "output", "Output file for the gene level analysis",
		  seqan::ArgParseArgument::STRING, "STRING"));
  addOption(parser, seqan::ArgParseOption(
		  "m", "mutation-position", "Output file to report position-specific mutation rates",
		  seqan::ArgParseArgument::STRING, "STRING"));
  addOption(parser, seqan::ArgParseOption(
		  "d", "mutation-distribution", "Output file to report mutation distribution",
		  seqan::ArgParseArgument::STRING, "STRING"));
  addOption(parser, seqan::ArgParseOption(
		  "g", "gff", "Annotation file (gff or gtf)",
		  seqan::ArgParseArgument::STRING, "STRING"));
  addOption(parser, seqan::ArgParseOption(
		  "c", "clip", "Clip c nucleotides from both side of each alignment in bam file", 
		  seqan::ArgParseArgument::INTEGER, "INT"));
  addOption(parser, seqan::ArgParseOption(
					  "c15", "c_r1_5p", "Clip nucleotides from read 1, 5 prime", 
					  seqan::ArgParseArgument::INTEGER, "INT"));
  addOption(parser, seqan::ArgParseOption(
					  "c13", "c_r1_3p", "Clip nucleotides from read 1, 3 prime", 
					  seqan::ArgParseArgument::INTEGER, "INT"));
  addOption(parser, seqan::ArgParseOption(
					  "c25", "c_r2_5p", "Clip nucleotides from read 2, 5 prime", 
					  seqan::ArgParseArgument::INTEGER, "INT"));
  addOption(parser, seqan::ArgParseOption(
					  "c23", "c_r2_3p", "Clip nucleotides from read 2, 3 prime", 
					  seqan::ArgParseArgument::INTEGER, "INT"));
  addOption(parser, seqan::ArgParseOption(
					  "c1", "clip_set_1", "Positions to be skipped in read 1, comma separated",
					  seqan::ArgParseArgument::STRING, "STRING"));
  addOption(parser, seqan::ArgParseOption(
					  "c2", "clip_set_2", "Positions to be skipped in read 2, comma separated",
					  seqan::ArgParseArgument::STRING, "STRING"));
  addOption(parser, seqan::ArgParseOption(
		  "q", "quality", "Minimal phread quality", 
		    seqan::ArgParseArgument::INTEGER, "INT"));  
  addOption(parser, seqan::ArgParseOption("p", "paired-end", "Is paired end sequencing."));
  addOption(parser, seqan::ArgParseOption("r", "reverse-complement", "primary read binds reverse strand."));
  addOption(parser, seqan::ArgParseOption("x", "output-exons", "Also output T2Cs for exon annotations."));
  addOption(parser, seqan::ArgParseOption("a", "include-ambiguous", "Include counts for ambiguous gene alignments."));


  setDefaultValue(parser, "readlength", "75");
  setDefaultValue(parser, "clip", "0");
  setDefaultValue(parser, "quality", "0");
  setDefaultValue(parser, "output", "gene-level-T2C.tsv");
  setDefaultValue(parser, "mutation-position","mutation-position.tsv");

  ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

  if (res != seqan::ArgumentParser::PARSE_OK) {
    return false;
  }
  
  getArgumentValue(a.bam_filename, parser, 0);
  
  getArgumentValue(a.fasta_filename, parser, 1);

  getOptionValue(a.gff_filename, parser, "gff");

  getOptionValue(a.quality_threshold, parser, "quality");

  getOptionValue(a.exclude_file, parser, "exclude-file");

  getOptionValue(a.output_file, parser, "output");

  getOptionValue(a.mutations_frequency_distribution_file, parser, "mutation-distribution");

  getOptionValue(a.position_specific_mutation_rate_file, parser, "mutation-position");

  getOptionValue(a.read_length, parser, "readlength");

  unsigned clip;
  getOptionValue(clip, parser, "clip");
  a.clip_read1_5p=clip;
  a.clip_read2_5p=clip;
  a.clip_read1_3p=clip;
  a.clip_read2_5p=clip;
  if (isSet(parser,"c_r1_5p"))  getOptionValue(a.clip_read1_5p, parser, "c_r1_5p");
  if (isSet(parser,"c_r2_5p"))  getOptionValue(a.clip_read2_5p, parser, "c_r2_5p");
  if (isSet(parser,"c_r1_3p"))  getOptionValue(a.clip_read1_3p, parser, "c_r1_3p");
  if (isSet(parser,"c_r2_3p"))  getOptionValue(a.clip_read2_3p, parser, "c_r2_3p");
  if (isSet(parser,"clip_set_1")) {
    string clip_set_tmp;
    getOptionValue(clip_set_tmp, parser, "clip_set_1");
    std::istringstream iss(clip_set_tmp);
    std::string token;
    while (std::getline(iss, token, ',')) 
      a.clip_set_read1.insert(stoi( token ) );
  }
  if (isSet(parser,"clip_set_2")) {
    string clip_set_tmp;
    getOptionValue(clip_set_tmp, parser, "clip_set_2");
    std::istringstream iss(clip_set_tmp);
    std::string token;
    while (std::getline(iss, token, ',')) 
      a.clip_set_read2.insert(stoi( token ) );
  }
  
  a.reverse_strand = isSet(parser, "reverse-complement");
  a.output_exons = isSet(parser, "output-exons");

  a.paired_end = isSet(parser, "paired-end");
  a.include_ambiguous = isSet(parser, "include-ambiguous");
  for (auto iter=a.clip_set_read1.begin(); iter!=a.clip_set_read1.end(); ++iter) {
    std::cerr << *iter << " " << std::endl;
  }
return true;
}
		   
