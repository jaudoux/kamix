#include <cstdlib>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <vector>
#include <unistd.h>
#include <string.h>
using namespace std;

#include "bgzf.h"
#include "kamix.h"

/*enum DNA_MAP {A, C, G, T};  // A=1, C=0, T=2, G=3*/
uint64_t str_to_int(char* str, size_t l)
{
  uint64_t strint = 0;
  for (size_t i = 0; i < l; i++) {
    uint8_t curr = 0;
    switch (str[i]) {
      case 'A': { curr = 0; break; }
      case 'T': { curr = 3; break; }
      case 'C': { curr = 1; break; }
      case 'G': { curr = 2; break; }
    }
    strint = strint << 2;
    strint = strint | curr;
  }
  return strint;
}

static inline uint64_t int_revcomp(uint64_t factor, uint32_t length) {
  uint64_t mask;
  if (length == 32)
    mask = ~0;
  else
    mask =  ((uint64_t) 1 << (2*length)) - 1;

  factor ^= mask;

  uint64_t mask_lsb;
  // Corresponds to the rightmost nucleotide
  mask_lsb = 3;
  uint64_t shift = 0;
  uint64_t result = 0;
  for(size_t j = 0; j < length; j++){
    result <<= 2;
    // get the leftmost nucleotide and put it at the end
    result |= (factor & mask_lsb) >> shift;
    mask_lsb <<= 2;
    shift += 2;
  }

  return result;
}

uint64_t canonical_kmer(uint64_t kmer, uint32_t length) {
  uint64_t rc = int_revcomp(kmer, length);
  if(rc < kmer) {
    return rc;
  }
  return kmer;
}

bool bgzf_getline_counting(BGZF * stream)
{
  int c = -1;
  while (true)
  {
    c = bgzf_getc (stream);
    if (c == -1)
    return true;
    else if (c == 10) // \n
    return false;
  }
}

/*
Create a gbi index for the file to facilitate
random access via the BGZF seek utility
*/
int create_kamix_index(string bgzf_file, int argc, char **argv)
{
  int c, help = 0;
  size_t chunk_size = CHUNK_SIZE;

  while ((c = getopt(argc, argv, "hc:")) >= 0) {
    switch (c) {
      case 'c': chunk_size = atoi(optarg); break;
      case 'h': help = 1; break;
    }
  }

  if (argc == 0 || help) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   kamix index [options] file:path\n\n");
    fprintf(stderr, "Options: -c INT      Chunk size (default %d)\n", CHUNK_SIZE);
    fprintf(stderr, "         -h          print this help message\n");
		return 1;
  }

  if (!bgzf_is_bgzf(bgzf_file.c_str()))
  {
    cerr << "[kamix] " << bgzf_file << " doesn't exist or wasn't compressed with bgzip" << endl;
    exit (1);
  }

  BGZF *bgzf_fp = bgzf_open(bgzf_file.c_str(), "r");
  if (bgzf_fp == NULL)
  {
    cerr << "[kamix] could not open file: " << bgzf_file << endl;
    exit (1);
  }

  // create an index for writing
  string index_file_name = bgzf_file + ".gbi.tmp";
  ofstream index_file(index_file_name.c_str(), ios::out);

  // add the offset for the end of the header to the index

  int status;
  kstring_t *line = new kstring_t;
  char *kmer;
  line->s = (char*)malloc(sizeof(char));
  line->s[0] = '\0';
  line->l = 0;
  line->m = 0;

  uint64_t prev_offset = 0;
  uint64_t offset = 0;
  while ((status = bgzf_getline(bgzf_fp, '\n', line)) >= 0)
  {
    offset = bgzf_tell(bgzf_fp);
    if (line->s[0] != '#')
    break;
    prev_offset = offset;
  }
  index_file << prev_offset << endl;

  // add the offsets for each CHUNK_SIZE
  // set of records to the index
  size_t chunk_count = 0;
  uint64_t total_lines = 1;
  chunk_t chunk;
  //vector<uint64_t> chunk_positions;
  vector<chunk_t> chunk_positions;

  chunk.offset = prev_offset;
  chunk.kmer = 0;
  chunk_positions.push_back(chunk);
  int eof = 1;
  while (true)
  {
    // grab the next line and store the offset
    eof = bgzf_getline(bgzf_fp, '\n', line);
    offset = bgzf_tell(bgzf_fp);
    chunk_count++;
    // stop if we have encountered an empty line
    if (eof < 0 || offset == prev_offset)
    {
      if (bgzf_check_EOF(bgzf_fp) == 1) {
        if (offset > prev_offset) {
          total_lines++;
          prev_offset = offset;
        }
        break;
      }
      break;
    }
    // store the offset of this chunk start
    else if (chunk_count == chunk_size)
    {
      kmer = strtok(line->s,"\t");
      chunk.offset = prev_offset;
      chunk.kmer = str_to_int(kmer,strlen(kmer));
      chunk_positions.push_back(chunk);
      chunk_count = 0;
    }
    total_lines++;
    prev_offset = offset;
  }
  chunk.offset = prev_offset;
  kmer = strtok(line->s,"\t");
  chunk.kmer = str_to_int(kmer,strlen(kmer));
  chunk_positions.push_back(chunk);
  bgzf_close(bgzf_fp);

  index_file << total_lines << endl;
  for (size_t i = 0; i < chunk_positions.size(); ++i)
  {
    index_file << chunk_positions[i].offset << "\t" << chunk_positions[i].kmer << endl;
    //index_file << chunk_positions[i].offset << endl;
  }
  index_file.close();

  return std::rename((bgzf_file + ".gbi.tmp").c_str(), (bgzf_file + ".gbi").c_str());
}

/*
Load an existing gbi index for the file to facilitate
random access via the BGZF seek utility
*/
void load_index(string bgzf_file, index_info &index)
{
  string index_file_name = bgzf_file + ".gbi";
  chunk_t chunk;
  // open the index file for reading
  ifstream index_file(index_file_name.c_str(), ios::in);

  if ( !index_file ) {
    cerr << "[kamix] could not find index file: " << index_file_name << ". Exiting!" << endl;
    exit (1);
  }
  else {
    string line;
    getline (index_file, line);
    index.header_end = atol(line.c_str());

    getline (index_file, line);
    index.num_lines = atol(line.c_str());

    while (index_file >> line)
    {
      chunk.offset = atol(line.c_str());
      index_file >> line;
      chunk.kmer = atol(line.c_str());
      index.chunk_offsets.push_back(chunk);
    }
  }
  index_file.close();
}

int get_kmer_length(string bgzf_file) {

  BGZF *bgzf_fp = bgzf_open(bgzf_file.c_str(), "r");

  bgzf_seek (bgzf_fp, 0, SEEK_SET);

  int status;
  kstring_t *line = (kstring_t*)calloc(1, sizeof(kstring_t));

  // Skip header lines starting with a # character
  while ((status = bgzf_getline(bgzf_fp, '\n', line)) > 0)
  {
    if (line->s[0] == '#')
      printf("%s\n", line->s);
    else break;
  }

  // Skip header line
  status = bgzf_getline(bgzf_fp, '\n', line);

  char *kmer = strtok(line->s,"\t");
  int kmer_length = strlen(kmer);

  free(line->s);
  free(line);

  return kmer_length;
}

int kamix_query(string bgzf_file, int argc, char **argv) {

  int c, use_canonical_kmer = 0;

  while ((c = getopt(argc, argv, "C")) >= 0) {
    switch (c) {
      case 'C': use_canonical_kmer = 1; break;
    }
  }

  if (argc == 0) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   kamix query [options] file:path mers:string+\n\n");
    fprintf(stderr, "Options: -C      use canonical kmers\n");
		return 1;
  }

  // load the index
  index_info index;
  load_index(bgzf_file, index);

  // load the BGZF file
  BGZF *bgzf_fp = bgzf_open(bgzf_file.c_str(), "r");
  if (bgzf_fp == NULL)
  {
    cerr << "[kamix] could not open file: " << bgzf_file << endl;
    exit (1);
  }

  int kmer_length = get_kmer_length(bgzf_file);

  for(int i = optind; i < argc; i++) {
    int query_length = strlen(argv[i]);
    if(query_length < kmer_length) {
      cerr << "Query is smaller that the k-mer size" << endl;
    } else {
      for(size_t j = 0; j <= strlen(argv[i]) - kmer_length; j++) {
        uint64_t kmer_query = str_to_int(&argv[i][j], kmer_length);
        if(use_canonical_kmer) {
          kmer_query = canonical_kmer(kmer_query, kmer_length);
        }

        if(i == optind && j == 0)
          get_kmer(bgzf_fp, index, kmer_query, 1);
        else
          get_kmer(bgzf_fp, index, kmer_query, 0);

      }
    }
  }

  return 1;
}

int get_kmer(BGZF *bgzf_fp, index_info index, uint64_t kmer_query, bool print_header)
{
  // FIXME check that k-mer is <= 32
  //uint64_t kmer_query = str_to_int(kmer,strlen(kmer));

  // dump the header if there is one
  int status;
  kstring_t *line = (kstring_t*)calloc(1, sizeof(kstring_t));
  char *kmer_tmp;
  uint64_t kmer_int = 0;

  bgzf_seek (bgzf_fp, 0, SEEK_SET);

  // Skip header lines starting with a # character
  while ((status = bgzf_getline(bgzf_fp, '\n', line)) > 0)
  {
    if (line->s[0] == '#')
      printf("%s\n", line->s);
    else break;
  }

  // Print header line
  if(print_header)
    printf("%s\n", line->s);

  status = bgzf_getline(bgzf_fp, '\n', line);

  // Find chunk dichotomy-style
  int i = index.chunk_offsets.size() / 2; // Start on the middle
  int min = 0;
  int max = index.chunk_offsets.size() - 1; // Max

  //fprintf(stderr, "dichotomy search\n");
  while(max > min) {
    //cerr << "min: " << min << " max: " << max << " i: " << i << endl;
    // We go up
    if(index.chunk_offsets[i].kmer > kmer_query) {
      max = i - 1;
    // We go down
    } else if(index.chunk_offsets[i+1].kmer <= kmer_query) {
      min = i + 1;
    // We found the right interval
    } else {
      break;
    }
    i = min + (max - min) / 2; // minddle
  }

  // We found the right interval
  char* tempstr = (char*)malloc(sizeof(char) * 33);
  if(max >= min) {
    bgzf_seek (bgzf_fp, index.chunk_offsets[i].offset, SEEK_SET);
    // Go to next line while the current k-mer is lower in lexicographic order
    while (kmer_int < kmer_query && status > 0)
    {
      status = bgzf_getline(bgzf_fp, '\n', line);
      strncpy(tempstr, line->s, 32);
      kmer_tmp = strtok(tempstr,"\t");
      kmer_int = str_to_int(kmer_tmp,strlen(kmer_tmp));
    }
    // If we have found the k-mer we print it, otherwise
    // the queried k-mers is not present in the file
    if(kmer_int == kmer_query) {
      printf("%s\n", line->s);
    } else {
      // FIXME We should print the kmer here
      // printf("%s\tnot_found\n", kmer);

      printf("kmer\tnot_found\n");
    }
  }
  return 0;
}

/*
Extract K random lines from file using reservoir sampling
*/
int random(string bgzf_file, uint64_t K)
{
  // load index into vector of offsets
  index_info index;
  load_index(bgzf_file, index);

  if ((uint64_t) K > index.num_lines)
  {
    cerr << "[kamix] warning: requested more lines than in the file." << endl;
    exit(1);
  }
  else {
    // load the BGZF file
    BGZF *bgzf_fp = bgzf_open(bgzf_file.c_str(), "r");
    if (bgzf_fp == NULL)
    {
      cerr << "[kamix] could not open file: " << bgzf_file << endl;
      exit (1);
    }

    // seed our random number generator
    size_t seed = (unsigned)time(0)+(unsigned)getpid();
    srand(seed);

    // reservoir sample
    uint64_t s = 0;
    uint64_t N = 0;
    uint64_t result_size = 0;
    vector<string> sample;
    int status;
    kstring_t *line = new kstring_t;
    line->s = (char*)malloc(sizeof(char));
    line->s[0] = '\0';
    line->l = 0;
    line->m = 0;

    while ((status = bgzf_getline(bgzf_fp, '\n', line)) != 0)
    {
      if (line->s[0] == '#')
      printf("%s\n", line->s);
      else break;
    }

    while ((status = bgzf_getline(bgzf_fp, '\n', line)) != 0)
    {
      N++;

      if (status < 0)
      break;

      if (result_size < K)
      {
        sample.push_back(line->s);
        result_size++;
      }
      else
      {
        s = (int) ((double)rand()/(double)RAND_MAX * N);
        if (s < K)
        sample[s] = line->s;
      }
    }
    bgzf_close(bgzf_fp);

    // report the sample
    for (size_t i = 0; i < sample.size(); ++i)
    printf("%s\n", sample[i].c_str());

  }
  return EXIT_SUCCESS;
}

/*
Return the total number of records in the index.
*/
int size(string bgzf_file)
{
  index_info index;
  load_index(bgzf_file, index);

  return index.num_lines;
}
