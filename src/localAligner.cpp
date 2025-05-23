#include "localAligner.hpp"

std::pair<int, int> LocalAligner::fillScoreMatrix(const std::string &a, const std::string &b,
                                                  std::vector<std::vector<int>> &dp, int &maxScore)
{
    int m = a.size();
    int n = b.size();
    std::pair<int, int> maxPos{0, 0};
    maxScore = 0;
    for (int i = 1; i <= m; ++i)
        for (int j = 1; j <= n; ++j)
        {
            int match = (a[i - 1] == b[j - 1]) ? matchScore : mismatchPenalty;
            int diag = dp[i - 1][j - 1] + match;
            int up = dp[i - 1][j] + gapPenalty;
            int left = dp[i][j - 1] + gapPenalty;
            dp[i][j] = std::max(0, std::max({diag, up, left}));
            if (dp[i][j] > maxScore)
            {
                maxScore = dp[i][j];
                maxPos = {i, j};
            }
        }
    return maxPos;
}

std::pair<std::string, std::string> LocalAligner::backtrackAlignment(const std::string &a, const std::string &b,
                                                                     const std::vector<std::vector<int>> &dp,
                                                                     std::pair<int, int> maxPos)
{
    std::string alignedA, alignedB;
    int i = maxPos.first;
    int j = maxPos.second;

    while (i > 0 && j > 0 && dp[i][j] != 0)
    {
        int score = dp[i][j];
        int scoreDiag = dp[i - 1][j - 1];
        int scoreUp = dp[i - 1][j];
        int scoreLeft = dp[i][j - 1];
        int match = (a[i - 1] == b[j - 1]) ? matchScore : mismatchPenalty;
        if (score == scoreDiag + match)
        {
            alignedA = a[i - 1] + alignedA;
            alignedB = b[j - 1] + alignedB;
            --i;
            --j;
        }
        else if (score == scoreUp + gapPenalty)
        {
            alignedA = a[i - 1] + alignedA;
            alignedB = '-' + alignedB;
            --i;
        }
        else
        {
            alignedA = '-' + alignedA;
            alignedB = b[j - 1] + alignedB;
            --j;
        }
    }
    return {alignedA, alignedB};
}

bool LocalAligner::checkAlternativeMax(const std::vector<std::vector<int>> &dp, std::pair<int, int> maxPos, int maxScore)
{
    int m = dp.size();
    int n = dp[0].size();
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < n; ++j)
            if ((i != maxPos.first || j != maxPos.second) && dp[i][j] == maxScore)
                return true;
    return false;
}

std::vector<std::pair<int, int>> LocalAligner::findAlternativeMaxPositions(
    const std::vector<std::vector<int>> &dp,
    std::pair<int, int> maxPos, int maxScore)
{
    int m = dp.size();
    int n = dp[0].size();
    std::vector<std::pair<int, int>> altPositions;
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < n; ++j)
            if ((i != maxPos.first || j != maxPos.second) && dp[i][j] == maxScore)
                altPositions.push_back({i, j});
    return altPositions;
}

LocalAlignmentResult LocalAligner::align(const std::string &a, const std::string &b)
{
    int m = a.size();
    int n = b.size();
    std::vector<std::vector<int>> dp(m + 1, std::vector<int>(n + 1, 0));
    int maxScore;
    std::pair<int, int> maxPos = fillScoreMatrix(a, b, dp, maxScore);
    auto bestAlignment = backtrackAlignment(a, b, dp, maxPos);
    int startI = maxPos.first - bestAlignment.first.size();
    int startJ = maxPos.second - bestAlignment.second.size();
    std::pair<int, int> startPosA = {startI, maxPos.first - 1};
    std::pair<int, int> startPosB = {startJ, maxPos.second - 1};
    std::string bestA = bestAlignment.first;
    std::vector<std::string> altAlignments;
    bool hasAlternative = false;
    if (maxScore > 0)
    {
        std::vector<std::pair<int, int>> altPositions = findAlternativeMaxPositions(dp, maxPos, maxScore);
        for (auto &pos : altPositions)
        {
            auto alignment = backtrackAlignment(a, b, dp, pos);
            altAlignments.push_back(alignment.first);
        }
        hasAlternative = !altAlignments.empty();
    }
    return LocalAlignmentResult{
        maxScore,
        dp,
        bestA,
        startPosA,
        startPosB,
        hasAlternative,
        altAlignments};
}
