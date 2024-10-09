OPTS=-g -Iseqan/include/ -O3 -DSEQAN_HAS_ZLIB=1 -Wall -std=c++14

OBJ=countT2C_new.o arguments.o annotation_file.o process_bam_alignment.o


default: all

clean: 
	rm *.o  countT2C

all: countT2C

countT2C: $(OBJ)
	c++ -o countT2C $(OBJ) -lz -lpthread

countT2C_new.o: countT2C_new.cpp arguments.hpp annotation_file.hpp process_bam_alignment.hpp
	c++ $(OPTS) -c countT2C_new.cpp

arguments.o: arguments.cpp arguments.hpp 
	c++ $(OPTS) -c arguments.cpp

annotation_file.o: annotation_file.cpp annotation_file.hpp 
	c++ $(OPTS) -c annotation_file.cpp

process_bam_alignment.o: process_bam_alignment.cpp process_bam_alignment.hpp annotation_file.hpp
	c++ $(OPTS) -c process_bam_alignment.cpp
