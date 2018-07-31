OBJS = fastq.o options.o file_util.o trim.o plot.o seq_overlap.o

CC = g++
cc = gcc

PROFILE = #-pg
OPENMP = -fopenmp

# Add -std=c++0x to enable <unordered_map>
#FLAGS = $(PROFILE) -O3 -Wall $(OPENMP) -std=c++0x -msse4.1 -DATA32
FLAGS = $(PROFILE) -O3 -Wall $(OPENMP) -std=c++0x -msse2 -DATA16

INC = -I.
LIBS = -lm -lz

.SUFFIXES : .o .cpp .c
.cpp.o:
	$(CC) $(FLAGS) $(PROFILE) $(INC) -c $<
.c.o:
	$(cc) $(FLAGS) $(PROFILE) $(INC) -c $<

all: FaQCs VTrim

FaQCs : $(OBJS) FaQCs.o
	$(CC) $(PROFILE) -o FaQCs $(OBJS) FaQCs.o $(LIBS) $(OPENMP)

VTrim : VTrim.o fastq.o
	$(CC) $(PROFILE) -o VTrim fastq.o VTrim.o $(LIBS) $(OPENMP)
	
clean:
	-rm -f *.o


