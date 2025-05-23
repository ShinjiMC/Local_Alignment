#pragma once
#include "structs.hpp"
#include <string>
#include <vector>
#include <algorithm>

class LocalAligner
{
private:
    int matchScore, mismatchPenalty, gapPenalty;
    std::pair<int, int> fillScoreMatrix(const std::string &a, const std::string &b,
                                        std::vector<std::vector<int>> &dp, int &maxScore);
    std::pair<std::string, std::string> backtrackAlignment(const std::string &a, const std::string &b,
                                                           const std::vector<std::vector<int>> &dp,
                                                           std::pair<int, int> maxPos);
    bool checkAlternativeMax(const std::vector<std::vector<int>> &dp, std::pair<int, int> maxPos, int maxScore);
    std::vector<std::pair<int, int>> findAlternativeMaxPositions(
        const std::vector<std::vector<int>> &dp,
        std::pair<int, int> maxPos, int maxScore);

public:
    LocalAligner(int match = 2, int mismatch = -1, int gap = -2)
        : matchScore(match), mismatchPenalty(mismatch), gapPenalty(gap) {}
    LocalAlignmentResult align(const std::string &a, const std::string &b);
};
