#pragma once
#include "structs.hpp"
#include <string>
#include <vector>
#include <utility>
#include <string>

class AlignmentIO
{
public:
    static std::pair<std::string, std::string> loadSequences(const std::string &filename);
    static void writeToFile(const std::string &filename,
                            const LocalAlignmentResult &result);
};
