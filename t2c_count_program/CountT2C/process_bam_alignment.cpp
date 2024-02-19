#include "process_bam_alignment.hpp"

#include "annotation_file.hpp"
#include <sstream> 


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
			   std::vector<counter> &gene_counter,
			   std::vector<counter> &exon_counter,
			   bool &is_intronic,
			   CharString &cellID
			   ) {
  cellID="";
  //   Extract Cell ID
  BamTagsDict tagsDict(record.tags);
  unsigned tagIdx = 0;
  if (findTagKey(tagIdx, tagsDict, "CB")) {
    if (!extractTagValue(cellID, tagsDict, tagIdx)) {
      std::cerr << "ERROR: Cannot extract key!\n";
    }
  }
 
  
  // ----
  // load reference sequence for contig if not allready done
  if (record.rID!=current_contig) {
    current_contig=record.rID;
    cerr << "read reference for " << contig_names[current_contig] << " ... ";
    unsigned idx = 0;
    if (!getIdByName(idx, faiIndex, contig_names[current_contig])) {
      cerr << "ERROR: FAI index has no entry for contig " << contig_names[current_contig] << endl; return(-1);
    }
    Dna5String seq;
    readSequence(seq, faiIndex, idx);
    clear(seqs);
    appendValue(seqs, seq);
    cerr << "done" << endl;
  }

  // ----
  // process current bam alignment
  Align<Dna5String> align;
  bamRecordToAlignment(align, seqs[0], record);
  Dna5 T='T';
  Dna5 C='C';
  char strand = '+';
  if (!a.reverse_strand) {
    if ( (hasFlagRC(record) & hasFlagFirst(record)) |  
	 (!hasFlagRC(record) & hasFlagLast(record)) |
	 (!a.paired_end & hasFlagRC(record)) ) { T='A'; C='G'; strand = '-';  }
  } else {
    if ( (!hasFlagRC(record) & hasFlagFirst(record)) |  
	 (hasFlagRC(record) & hasFlagLast(record)) |
	 (!a.paired_end & !hasFlagRC(record)) ) { T='A'; C='G'; strand = '-'; }
  }
 
  
  auto & row1 = row(align, 0);	  // reference sequence
  auto & row2 = row(align, 1);	  // read sequence
  
  // now find the gene (if annotation is given)
  is_intronic=false;
  String<unsigned> annotation_index;
  if (length(intervalTrees)>(unsigned int)record.rID) {
    findIntervals(annotation_index, intervalTrees[record.rID], toSourcePosition(row1, 0), toSourcePosition(row1, length(row1)-1));
    // increase read counter for each overlapping annotation given the id in the interval tree
    if (length(annotation_index)>0) {
      for (unsigned j = 0; j < length(annotation_index); ++j) {
	if (gene_annotation[record.rID][annotation_index[j]].type=="gene") {
	  if ((gene_annotation[record.rID][annotation_index[j]].beginPos<=toSourcePosition(row1, 0)) &&
	      (gene_annotation[record.rID][annotation_index[j]].endPos>=toSourcePosition(row1, length(row1)-1)) &&
	      (gene_annotation[record.rID][annotation_index[j]].strand==strand)) {
	    bool gene_known=false;
	    for (auto iter=gene_counter.begin(); iter!=gene_counter.end(); ++iter) {
	      if (iter->pos==annotation_index[j]) { gene_known=true; break; }
	    }
	    if (!gene_known) {
	      counter tmp;
	      tmp.pos=annotation_index[j];
	      tmp.Ts=0;
	      tmp.T2Cs=0;
	      gene_counter.push_back(tmp);
	    }
	  }
	}
      }
    }
  } 

  for (unsigned i = 0; i < length(row1); ++i) {
    if (row1[i]==T) {
      if (row2[i]!='-') {          // is aligned (read sequence is not "-")
	// check if the reference position is not in the clipped region
	auto read_position=toSourcePosition(row2, i);
	if (hasFlagRC(record)) { read_position=a.read_length-toSourcePosition(row2,i); }
	
	SEQAN_ASSERT_LEQ(read_position,a.read_length);
	SEQAN_ASSERT_GEQ(read_position,(unsigned)0);

	if (((hasFlagFirst(record)) & (read_position>=a.clip_read1_5p) & (read_position<a.read_length-a.clip_read1_3p) & (a.clip_set_read1.find(read_position)==a.clip_set_read1.end())) |
	    ((hasFlagLast(record)) & (read_position>=a.clip_read2_5p) & (read_position<a.read_length-a.clip_read2_3p) & (a.clip_set_read2.find(read_position)==a.clip_set_read2.end())) |
	    (!a.paired_end & (read_position>=a.clip_read1_5p) & (read_position<a.read_length-a.clip_read1_3p) & (a.clip_set_read1.find(read_position)==a.clip_set_read1.end()))) {

	  // check if the reference position is not in the exclude list.
	  auto ref_position=toSourcePosition(row1, i);
	  if (!binary_search(positions_to_exclude[record.rID].begin(),
			     positions_to_exclude[record.rID].end(),
			     ref_position )) {
	    // is the phred score OK? 

	    if (((char)record.qual[toSourcePosition(row2, i)])-33>=a.quality_threshold) {

	      for (auto iter=gene_counter.begin(); iter!=gene_counter.end(); iter++)
		iter->Ts++;
	      if (row2[i]==C) {
                for (auto iter=gene_counter.begin(); iter!=gene_counter.end(); iter++)
		  iter->T2Cs++;
		if (hasFlagFirst(record)) {
		  mutation_rate_read1[read_position]++;
		} else {
		  mutation_rate_read2[read_position]++;
		}
	      }
	      if (length(annotation_index)>0) {
		// check if the position maps to any exon...
		bool maps_to_exon=false;
		for (unsigned j=0; j<length(annotation_index); ++j) {
		  if ((gene_annotation[record.rID][annotation_index[j]].beginPos<= ref_position) &&
		      (gene_annotation[record.rID][annotation_index[j]].endPos>= ref_position) &&
		      (gene_annotation[record.rID][annotation_index[j]].strand==strand) &&
		      (gene_annotation[record.rID][annotation_index[j]].type=="exon") ){
		    maps_to_exon=true;
                    bool found_exon=false;
		    for (vector<counter>::iterator iter=exon_counter.begin(); iter<exon_counter.end(); iter++) {
		      if (iter->pos==annotation_index[j]) {
                        found_exon=true;
			iter->Ts++;
			if (row2[i]==C) { iter->T2Cs++; }
		      }
		    }
                    if (!found_exon) {
		      counter tmp;
		      tmp.pos=annotation_index[j];
		      tmp.Ts=1;
		      tmp.T2Cs=0;
		      if (row2[i]==C) { tmp.T2Cs++; }
		      exon_counter.push_back(tmp);
                    }
		  }
		}
		if (!maps_to_exon) {
		  is_intronic=true;
		}
	      }
	    }
	  }
	}
      } 
    }
  }
  
  return true;
}
