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
int create_kamix_index(string bgzf_file)
{

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
    else if (chunk_count == CHUNK_SIZE)
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

int kamix_query(string bgzf_file, int argc, char **argv) {

  if (argc == 0) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   kamix query [options] file:path mers:string+\n\n");
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

  for(int i = 0; i < argc; i++) {
    //cerr << "counting : " << argv[i] << endl;
    if(i == 0)
      get_kmer(bgzf_fp, index, argv[i], 1);
    else
      get_kmer(bgzf_fp, index, argv[i], 0);
  }

  return 1;
}

int get_kmer(BGZF *bgzf_fp, index_info index, char *kmer, bool print_header)
{
  // FIXME check that k-mer is <= 32
  uint64_t kmer_query = str_to_int(kmer,strlen(kmer));

  // dump the header if there is one
  int status;
  kstring_t *line = (kstring_t*)calloc(1, sizeof(kstring_t));
  char *kmer_tmp;
  uint64_t kmer_int = 0;

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
    while (kmer_int < kmer_query)
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
    } else if(kmer_int > kmer_query) {
      printf("kmer %s not found\n", kmer);
    }
  }
  return 0;
}

/*
Extract lines [FROM, TO] from file.
*/
int grab(string bgzf_file, uint64_t from_line, uint64_t to_line)
{
  // load index into vector of offsets
  index_info index;
  load_index(bgzf_file, index);

  if ((from_line > index.num_lines)
  ||
  (to_line > index.num_lines))
  {
    cerr << "[kamix] requested lines exceed the number of lines in the file." << endl;
    exit(1);
  }
  else if (from_line < 0)
  {
    cerr << "[kamix] indexes must be positive numbers." << endl;
    exit(1);
  }
  else if (from_line > to_line)
  {
    cerr << "[kamix] requested end line is less than the requested begin line." << endl;
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

    // dump the header if there is one
    int status;
    kstring_t *line = new kstring_t;
    line->s = (char*)malloc(sizeof(char));
    line->s[0] = '\0';
    line->l = 0;
    line->m = 0;

    while ((status = bgzf_getline(bgzf_fp, '\n', line)) > 0)
    {
      if (line->s[0] == '#')
      printf("%s\n", line->s);
      else break;
    }

    // easier to work in 0-based space
    uint64_t from_line_0  = from_line - 1;
    // get the chunk index for the requested line
    uint64_t requested_chunk = from_line_0 / CHUNK_SIZE;
    // derive the first line in that chunk
    uint64_t chunk_line_start = (requested_chunk * CHUNK_SIZE);

    // jump to the correct offset for the relevant chunk
    // and fast forward until we find the requested line
    bgzf_seek (bgzf_fp, index.chunk_offsets[requested_chunk].offset, SEEK_SET);
    while (chunk_line_start <= from_line_0)
    {
      status = bgzf_getline(bgzf_fp, '\n', line);
      chunk_line_start++;
    }
    // now, print each line until we reach the end of the requested block
    do
    {
      printf("%s\n", line->s);
      status = bgzf_getline(bgzf_fp, '\n', line);
      chunk_line_start++;
    } while (chunk_line_start <= to_line);
  }
  return EXIT_SUCCESS;
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
