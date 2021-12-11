#include <simplicial_arrangement/simplicial_arrangement.h>
#include "look_up_table.h"

#include <stdlib.h>
#include <iostream>

using namespace simplicial_arrangement;

int main(int argc, const char* argv[]) {
    std::cout << "going to load table ..." << std::endl;
    bool loaded = load_lookup_table();

	std::cout << "loaded : " << loaded << std::endl;

	std::cout << "1 func table size = " << one_func_lookup_table->size() << std::endl;
    std::cout << "2 func table size = " << two_func_lookup_table->size() << std::endl;

	// test



	return 0;
}