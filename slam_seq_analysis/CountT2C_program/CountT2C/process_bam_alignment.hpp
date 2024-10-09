#ifndef PROCESS_BAM_ALIGNMENT_HPP
#define PROCESS_BAM_ALIGNMENT_HPP

#include <vector>
#include <seqan/bam_io.h>
#include <seqan/align.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include "annotation_file.hpp"
#include "arguments.hpp"


using namespace seqan;
using namespace std;


typedef FormattedFileContext<BamFileIn, void>::Type TBamContext;


struct counter {
  unsigned pos; // position in annotation 
  unsigned Ts; // Number Ts;
  unsigned T2Cs; // Number T2Cs;
};

bool process_bam_alignment(const BamAlignmentRecord &record,
			   const StringSet<CharString> &contig_names,
			   StringSet<Dna5String> &seqs,
			   int &current_contig,
			   const FaiIndex &faiIndex,
			   const seqan::String<seqan::String<GeneAnnotation> > & gene_annotation,
			   const seqan::String<TIntervalTree> & intervalTrees,
			   const vector<vector<int> > &positions_to_exclude,
			   const Arguments &a,
			   vector<unsigned> &mutation_rate_read1,
			   vector<unsigned> &mutation_rate_read2,
			   std::vector<counter>  &gene_counter,
			   std::vector<counter> &exon_counter,
			   bool &is_intronic,
			   CharString &cellID
			   );


#endif
