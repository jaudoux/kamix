#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

#include "bgzf.h"
#include "kamix.h"

#define VERSION "0.0.2"

static int usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:   kamix <command> file:path <arguments>\n");
	fprintf(stderr, "Version: %s\n\n", VERSION);
	fprintf(stderr, "Command: index      Create a kamix index (kmers.gz.gbi)\n");
	fprintf(stderr, "         query      Query a k-mer\n");
	fprintf(stderr, "         random     Extract the 100 random k-mers\n");
	fprintf(stderr, "         check      Is the file bgzipped?\n");
	fprintf(stderr, "         size       Get total number of k-mers in the file\n");
	fprintf(stderr, "\n");
	return 1;
}

int main (int argc, char **argv)
{
    if (argc == 1) {
      return usage();
    } else if (argc == 2) {
      fprintf(stderr, "Missing kmers matrice\n");
      return usage();
    } else if (argc >= 3) {
        // create input file for the purpose of the example
        string bgzf_file = argv[2];

        string sub_command = argv[1];

        if (sub_command == "index")
        {
            create_kamix_index(bgzf_file);
        }
        else if (sub_command == "grab")
        {
            int64_t from_line = atol(argv[3]);
            int64_t to_line = from_line;
            if (argc == 5)
                to_line = atol(argv[4]);

            grab(bgzf_file, from_line, to_line);
        }
        else if (sub_command == "query")
        {
            kamix_query(bgzf_file, argc-3, argv+3);
        }
        else if (sub_command == "random")
        {
            size_t N = atoi(argv[3]);
            random(bgzf_file, N);
        }
        else if (sub_command == "check")
        {
          cout << ((bgzf_is_bgzf(bgzf_file.c_str()) == 1) ? "yes" : "no") << "\n";
        }
        else if (sub_command == "size")
        {
          cout << size(bgzf_file) << "\n";
        } else {
          cout << "unknown command: " << sub_command << endl;
          cout << "available commands are: index, grab, random, check, size" << endl;
        }
    }

    return EXIT_SUCCESS;
}
