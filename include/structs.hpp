#pragma once
#include <string>
#include <vector>
#include <utility>
#include <string>

struct LocalAlignmentResult
{
    int maxScore;
    std::vector<std::vector<int>> scoreMatrix;
    std::string bestAlignment;
    std::pair<int, int> startPosA;
    std::pair<int, int> startPosB;
    bool hasAlternative;
    std::vector<std::string> alternativeAlignments;
};