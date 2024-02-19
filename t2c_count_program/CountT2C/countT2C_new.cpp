#include <string>
#include <iostream>
#include <vector>
#include <cassert>
#include "arguments.hpp"
#include "annotation_file.hpp"
#include "process_bam_alignment.hpp"

using namespace std;
using namespace seqan;


void unify_gene_counter(std::vector<counter> &gene_counter1, const std::vector<counter> gene_counter2) {
  for (auto iter1=gene_counter1.begin(); iter1!=gene_counter1.end(); /* empty */ ) {
    bool found=false;
    for (auto iter2=gene_counter2.begin(); iter2!=gene_counter2.end(); iter2++) {
      if (iter1->pos == iter2->pos) { 
	iter1->Ts+=iter2->Ts; 
	iter1->T2Cs+=iter2->T2Cs; 
	found=true;
      }
    }
    if (!found) {
      gene_counter1.erase(iter1);
    } else {
      iter1++;
    }
  }
}



int main(int argc, char const ** argv) {

  Arguments a;
  if (!getArguments(a,argc,argv)) return 0;

  // --- BAM File Stuff
  unsigned int total_read_length=a.read_length;
  if (a.paired_end) total_read_length=total_read_length*2;
  vector<int> matrix(total_read_length*total_read_length,0);
  vector<int> intronic_matrix(total_read_length*total_read_length,0);
  String<String<GeneAnnotation> > gene_annotation;
  vector<unsigned> mutation_rate_read1(a.read_length+1,0);
  vector<unsigned> mutation_rate_read2(a.read_length+1,0);

  unsigned ambigous_alignments=0;
  unsigned not_aligned_to_genes=0;
  unsigned uniquely_aligned_to_genes=0;
  unsigned uniquely_aligned_to_introns=0;

  try   {
    cerr << "loading bam file: " << a.bam_filename << endl;
    
    BamFileIn bamFileIn(a.bam_filename.c_str());
    BamHeader header;
    readHeader(header, bamFileIn);
    TBamContext const & bamContext = context(bamFileIn);
    StringSet<CharString> contig_names = contigNames(bamContext);
    TNameStoreCache contig_names_cache = contigNamesCache(bamContext);
    cerr << length(contig_names) << endl;
    int numPrinted = 0, num=10000;
    
    // --- Loaad Reference Fasta index
    cerr << "loading fasta index: " << a.fasta_filename << endl;
    FaiIndex faiIndex;
    if (!open(faiIndex, a.fasta_filename.c_str())) {
      cerr << "ERROR: Could not load FAI index " << a.fasta_filename << ".fai\n";
      return -1;
    }

    
    vector<vector<int> > positions_to_exclude(length(contig_names));
    if (a.exclude_file!="") {
      cerr << "loading exclude_file file: " << a.exclude_file << endl;
      int pos;
      string contig;
      unsigned idx = 0;
      ifstream file(a.exclude_file.c_str(),ios_base::in);
      while (file >> contig >> pos) {
	if (getIdByName(idx, contig_names_cache, contig))
	  positions_to_exclude[idx].push_back(pos-1);
	else
	  cerr << "problem, cannot find contig " << contig << " from exclude file in bam file header!" << endl;
      }
      int total=0;
      for (auto iter=positions_to_exclude.begin(); iter!=positions_to_exclude.end(); iter++) {
	sort(iter->begin(), iter->end());
	total += iter->size();
      }
      std::cerr << "using " << a.exclude_file << " file with " << total <<  " positions of mutations that should be ignored." << endl;   
    }
    StringSet<Dna5String> seqs;
    String<TIntervalTree> intervalTrees;
    if (a.gff_filename!="") {
      loadAnnotationFile(intervalTrees, gene_annotation, contig_names_cache, a.gff_filename);
    }
    // prepare data
    int current_contig = -1; 
    // Read header.
    
    BamAlignmentRecord record;

    bool is_intronic, ambigous, found_gene, is_intronic_second=false;
    std::vector<counter> gene_counter;
    std::vector<counter> exon_counter;

    CharString past_name = "";
    
    while (!atEnd(bamFileIn) ) {
      readRecord(record, bamFileIn);
      numPrinted++;
      
      // Quality filtering: skip secondary, unmapped or not properly mapped paired sequences
      
      if (hasFlagSecondary(record) | hasFlagUnmapped(record) ) continue;
      if (a.paired_end & (!hasFlagAllProper(record))) continue;
      
      if (record.rID != BamAlignmentRecord::INVALID_REFID && record.rID<(int)length(contig_names) ) {
	CharString cellID="";
	CharString current_name = record.qName;
	if (current_name == past_name) {
	  std::vector<counter> gene_counter2;
	  process_bam_alignment(record, // input variables
				contig_names,
				seqs,
				current_contig,
				faiIndex,
				gene_annotation,
				intervalTrees,
				positions_to_exclude,
				a,
				mutation_rate_read1, // output variables
				mutation_rate_read2,
				gene_counter2,
				exon_counter,
				is_intronic_second,
				cellID
				);
	  unify_gene_counter(gene_counter,gene_counter2);
	  is_intronic=is_intronic||is_intronic_second;
	  ambigous=gene_counter.size()>1;
	  found_gene=gene_counter.size()>0;
	  
	} else {
	  exon_counter.clear();
          gene_counter.clear();
	  process_bam_alignment(record, // input variables
				contig_names,
				seqs,
				current_contig,
				faiIndex,
				gene_annotation,
				intervalTrees,
				positions_to_exclude,
				a,
				mutation_rate_read1, // output variables
				mutation_rate_read2,
				gene_counter,
				exon_counter,
				is_intronic,
				cellID
				);
	  ambigous=gene_counter.size()>1;
	  found_gene=gene_counter.size()>0;
	}

	if ((current_name == past_name)||(!a.paired_end)) {


	  for (auto iter=gene_counter.begin(); iter!=gene_counter.end(); ++iter) {
	    // sanity check: More Ts then read_length: read_length parameter incorrect!
	    if (iter->Ts>total_read_length) {
	      cerr << "More Ts than readlength, adjust readlength parameter" << endl;
	      exit(-1);
	    }
	  }
	  
	  // when gene annotation was given, add to the counters
	  if (length(intervalTrees)>0) {
	    if (found_gene) {
	      assert(gene_counter.size()>0);
	      if (!ambigous) {
		uniquely_aligned_to_genes++;
	      }	else {
		ambigous_alignments++;
	      }
	      if ((!ambigous) || (a.include_ambiguous)) {

		if (is_intronic) {
		  intronic_matrix[gene_counter[0].Ts*total_read_length+gene_counter[0].T2Cs]++;			
		  uniquely_aligned_to_introns++;
		  for (auto iter=gene_counter.begin(); iter!=gene_counter.end(); ++iter) {
		    gene_annotation[record.rID][iter->pos].counter[cellID].intronicTs+=iter->Ts;
		    gene_annotation[record.rID][iter->pos].counter[cellID].intronicT2Cs+=iter->T2Cs;		
		  }
		} else {
		  for (auto iter=gene_counter.begin(); iter!=gene_counter.end(); ++iter) {
		    gene_annotation[record.rID][iter->pos].counter[cellID].Ts+=iter->Ts;
		    gene_annotation[record.rID][iter->pos].counter[cellID].T2Cs+=iter->T2Cs;
		  }
		  for (auto iter=exon_counter.begin(); iter<exon_counter.end(); iter++) {
		    gene_annotation[record.rID][iter->pos].counter[cellID].Ts+=iter->Ts;
		    gene_annotation[record.rID][iter->pos].counter[cellID].T2Cs+=iter->T2Cs;
		  }

		}
	      }
	    } else {
	      not_aligned_to_genes++;  
	    } 
	  }
	  // add to the T/T2C count matrix
	  if (gene_counter.size()) {
	    matrix[gene_counter[0].Ts*total_read_length+gene_counter[0].T2Cs]++;
	  }
	}
	past_name=current_name;
      }
      if ((numPrinted % num) == 0) { cerr << numPrinted/1000 << "k sequences processed ... " << endl; }
    }
  }
  catch (Exception const & e)
    {
      cerr << "ERROR: " << e.what() << endl;
      return 1;
    }
  
  if (a.mutations_frequency_distribution_file!="") {
    ofstream os(a.mutations_frequency_distribution_file);
    os << "num_Ts\tnum_T2Cs\tfreq_exons\tfreq_introns" << endl;
    for (unsigned i=0; i<total_read_length; ++i) {
      for (unsigned j=0; j<total_read_length; ++j) {
	if ((matrix[i*total_read_length+j]>0)||(intronic_matrix[i*total_read_length+j]>0)) {
	  os << i << "\t" << j << "\t" << matrix[i*total_read_length+j]  << "\t" << intronic_matrix[i*total_read_length+j] << endl;
	}
      }
    }
  }
 
  if ((length(gene_annotation)>0)&(a.output_file!="")) {
    ofstream os(a.output_file);
    os << "gene_id\tcell_barcode\texon_Ts\texon_T2Cs\tintron_Ts\tintron_T2Cs" << endl;
    for (unsigned i=0; i<length(gene_annotation); i++) {
      for (unsigned j=0; j<length(gene_annotation[i]); j++) {
	if ((gene_annotation[i][j].type=="gene") ||(a.output_exons && (gene_annotation[i][j].type=="exon"))  ) {
	  for (auto iter=gene_annotation[i][j].counter.begin(); 
	       iter!=gene_annotation[i][j].counter.end(); iter++) {
	    os << gene_annotation[i][j].name << "\t" << iter->first << "\t" << iter->second.Ts << "\t" << iter->second.T2Cs << "\t" << iter->second.intronicTs << "\t" << iter->second.intronicT2Cs<< endl;
	  }
	}
      }
    }  
  }

  if (a.position_specific_mutation_rate_file!="") {
    ofstream os(a.position_specific_mutation_rate_file);
    os << "position\tmut_read1\tmut_read2" << endl;
    for (unsigned i=0; i<mutation_rate_read1.size(); i++) {
      os << i << "\t" << mutation_rate_read1[i] << "\t" << mutation_rate_read2[i] << endl;
    }  
  }

  
  cerr << "uniquely_aligned_to_genes: " << uniquely_aligned_to_genes << endl;
  cerr << "uniquely_aligned_to_introns: " << uniquely_aligned_to_introns << endl;
      
  cerr << "ambigous_alignments: " << ambigous_alignments << endl;
  cerr << "not_aligned_to_genes: " << not_aligned_to_genes << endl;

  return 0; 
}
