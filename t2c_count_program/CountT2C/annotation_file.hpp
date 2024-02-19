
#ifndef ANNOTATION_FILE_HPP
#define ANNOTATION_FILE_HPP
#include <string>
#include <seqan/sequence.h>
#include <seqan/misc/interval_tree.h>
#include <seqan/misc/name_store_cache.h>

typedef seqan::IntervalTree<unsigned, unsigned> TIntervalTree;
typedef seqan::NameStoreCache<seqan::StringSet<seqan::CharString> >  TNameStoreCache;

struct Counter {
  Counter(): Ts(0), T2Cs(0), intronicTs(0), intronicT2Cs(0) {}
  unsigned Ts;
  unsigned T2Cs;
  unsigned intronicTs;
  unsigned intronicT2Cs;
};

struct GeneAnnotation {
  GeneAnnotation() : type(), strand('+'), name(), beginPos(0), endPos(0) {}
  // Feature
  seqan::CharString type;
  char strand;
  seqan::CharString name;
  unsigned beginPos;
  unsigned endPos;

  // Counters
  std::map<seqan::CharString,Counter> counter;
};


bool loadAnnotationFile(seqan::String<TIntervalTree> & intervalTrees,
			seqan::String<seqan::String<GeneAnnotation> > & gene_annotation,
			TNameStoreCache &nameCache,
			std::string const & annotationFileName);

#endif
