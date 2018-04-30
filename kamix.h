#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

#include "bgzf.h"

// we only want to store the offset for every 10000th
// line. otherwise, were we to store the position of every
// line in the file, the index could become very large for
// files with many records.
#define CHUNK_SIZE 4096
#define KMER_LENGTH 32

typedef struct {
  uint64_t offset;
  uint64_t kmer;
} chunk_t;

struct index_info {
    uint64_t header_end;
    vector<chunk_t> chunk_offsets;
    uint64_t num_lines;
};

bool bgzf_getline_counting(BGZF * stream);

/*
Create a gbi index for the file to facilitate
random access via the BGZF seek utility
*/
int create_kamix_index(string bgzf_file, int argc, char **argv);

/*
Load an existing gbi index for the file to facilitate
random access via the BGZF seek utility
*/
void load_index(string bgzf_file, index_info &index);


int kamix_query(string bgzf_file, int argc, char **argv);

int get_kmer(BGZF *bgzf_fp, index_info index, uint64_t kmer_query, bool print_header = 0);

int get_kmer_length(string bgzf_file);

/*
Extract K random lines from file using reservoir sampling
*/
int random(string bgzf_file, uint64_t K);

/*
Return the total number of records in the index.
*/
int size(string bgzf_file);
