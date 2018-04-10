all:
	gcc -Wall -O2 -o kamix kamix_main.cpp kamix.cpp bgzf.c -lstdc++ -lz
