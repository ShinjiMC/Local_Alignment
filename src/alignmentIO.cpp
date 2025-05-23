#include "alignmentIO.hpp"
#include <fstream>
#include <iostream>
#include <stdexcept>

std::pair<std::string, std::string> AlignmentIO::loadSequences(const std::string &filename)
{
    std::ifstream file(filename);
    if (!file)
        throw std::runtime_error("Could not open file: " + filename);
    std::string line1, line2;
    if (!std::getline(file, line1) || !std::getline(file, line2))
        throw std::runtime_error("Error reading lines from file: " + filename);
    return {line1, line2};
}

void AlignmentIO::writeToFile(const std::string &filename,
                              const LocalAlignmentResult &result)
{
    std::ofstream out(filename);
    if (!out)
    {
        std::cerr << "Could not open file: " << filename << "\n";
        return;
    }
    out << "Maximum Score: " << result.maxScore << "\n\n";
    out << "Best Local Alignment:\n";
    out << result.bestAlignment << "\n\n";
    out << "Position in Sequence A: [" << result.startPosA.first
        << ", " << result.startPosA.second << "]\n";
    out << "Position in Sequence B: [" << result.startPosB.first
        << ", " << result.startPosB.second << "]\n\n";
    out << "Matrix:\n";
    for (const auto &row : result.scoreMatrix)
    {
        for (int val : row)
            out << val << "\t";
        out << "\n";
    }
    out << "\nAlternative subsequence with same max score exists: "
        << (result.hasAlternative ? "Yes" : "No") << "\n\n";
    if (result.hasAlternative)
    {
        int idx = 1;
        for (const auto &a1 : result.alternativeAlignments)
        {
            out << "Alternative #" << idx++ << ":\n";
            out << a1 << "\n\n";
        }
    }
}
