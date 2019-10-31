/**
 * myzcat.cpp
 * A zcat replacement, for educational purpose.
 * Uses boost::iostream and zlib.
 * 
 * Compile like this:
 *   clang++ -o myzcat myzcat.cpp -lz -lboost_iostreams
 * 
 * This code is published as public domain.
 */
#include <fstream>
#include <iostream>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

int main(int argc, char** argv) {
    if(argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <gzipped input file>" << std::endl;
    }
    //Read from the first command line argument, assume it's gzipped
    std::ifstream file(argv[1], std::ios_base::in | std::ios_base::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
    inbuf.push(boost::iostreams::gzip_decompressor());
    inbuf.push(file);
    //Convert streambuf to istream
    std::istream instream(&inbuf);
    //Iterate lines
    std::string line;
    while(std::getline(instream, line)) {
        std::cout << line << std::endl;
    }
    //Cleanup
    file.close();
}
