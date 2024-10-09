#include "annotation_file.hpp"

#include <seqan/gff_io.h>
#include <seqan/seq_io.h>
#include <seqan/misc/interval_tree.h>

using namespace std;
using namespace seqan;

typedef IntervalAndCargo<unsigned, unsigned> TInterval;

bool loadAnnotationFile(String<TIntervalTree> & intervalTrees,
			String<String<GeneAnnotation> > & gene_annotation,
			TNameStoreCache &nameCache,
			std::string const & annotationFileName)
{
  cerr << "load annotation file: " << annotationFileName << endl;
  String<String<TInterval> >intervals;
  clear(gene_annotation);
  resize(intervals,length(nameCache));
  resize(gene_annotation,length(nameCache));
  GffFileIn gffIn(annotationFileName.c_str());
  GffRecord record;
  while (!atEnd(gffIn))
    {
      readRecord(record, gffIn);
      if ((record.type=="gene") | (record.type=="exon")) {
        int32_t beginPos = record.beginPos;
        int32_t endPos = record.endPos;
	if (beginPos > endPos)
	  std::swap(beginPos, endPos);
        unsigned contigId = 0;
	GeneAnnotation anno;
	anno.strand=record.strand;
	anno.type=record.type;
	anno.beginPos=beginPos;
	anno.endPos=endPos;
	for (unsigned i=0; i<length(record.tagValues); i++) {
	  if ((record.type=="gene") && (record.tagNames[i] == "gene_id")) {
	    anno.name=record.tagValues[i];
	  }
	  if ((record.type=="exon") && (record.tagNames[i] == "exon_id")) {
	    anno.name=record.tagValues[i];
	  }
	}
	getIdByName(contigId, nameCache, record.ref);
	if (contigId>=length(intervals)) { resize(intervals,contigId+1); }
	if (contigId>=length(gene_annotation)) { resize(gene_annotation,contigId+1); }
        // insert forward-strand interval of the gene and its annotation id
        appendValue(intervals[contigId], TInterval(beginPos, endPos,length(gene_annotation[contigId])));
	appendValue(gene_annotation[contigId], anno);
      }
    }

  int numContigs = length(intervals);
  resize(intervalTrees, numContigs);
  std::cerr << length(intervals) << " " << length(intervalTrees) << " " << length(nameCache) << std::endl;
  
  for (int i = 0; i < numContigs; ++i)
    createIntervalTree(intervalTrees[i], intervals[i]);

  return true;
}
