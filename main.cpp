#include <iostream>
#include "alignmentIO.hpp"
#include "localAligner.hpp"

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        std::cerr << "Use: " << argv[0] << " <path_txt>\n";
        return 1;
    }
    std::string filename = argv[1];
    try
    {
        auto [s1, s2] = AlignmentIO::loadSequences(filename);
        LocalAligner aligner(1, -1, -2);
        auto result = aligner.align(s1, s2);
        AlignmentIO::writeToFile("alignment.txt", result);
        std::cout << "Completed alignment. Results in 'alignment.txt'\n";
    }
    catch (const std::exception &ex)
    {
        std::cerr << "Error: " << ex.what() << "\n";
        return 1;
    }
    return 0;
}
