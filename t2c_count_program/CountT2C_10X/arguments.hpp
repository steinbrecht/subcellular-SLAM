#ifndef ARGUMENTS_HPP
#define ARGUMENTS_HPP

#include <string>
#include <set>


using namespace std;

struct Arguments {
  string bam_filename;
  string fasta_filename;
  string exclude_file;
  string output_file;
  string position_specific_mutation_rate_file;
  string mutations_frequency_distribution_file;
  string gff_filename;

  unsigned read_length;  
  unsigned clip_read1_5p;
  unsigned clip_read2_5p;
  unsigned clip_read1_3p;
  unsigned clip_read2_3p;
  int quality_threshold;
  bool reverse_strand;
  bool paired_end;
  bool output_exons;
  bool include_ambiguous;
  std::set<unsigned> clip_set_read1;
  std::set<unsigned> clip_set_read2;
};

bool getArguments(Arguments &a, int argc, char const ** argv);

#endif // ARGUMENTS_HPP
