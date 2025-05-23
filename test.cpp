#include <gtest/gtest.h>
#include "localAligner.hpp"

LocalAligner aligner(1, -1, -1);

TEST(LocalAlignerTest, BasicMatchTest)
{
    auto result = aligner.align("GATACAA", "CGATACAT");
    EXPECT_EQ(result.maxScore, 6);
    EXPECT_EQ(result.bestAlignment, "GATACA");
    EXPECT_EQ(result.startPosA, std::make_pair(0, 5));
    EXPECT_EQ(result.startPosB, std::make_pair(1, 6));
    ASSERT_EQ(result.scoreMatrix[3][4], 3);
}

TEST(LocalAlignerTest, AlternativeAlignmentsTest)
{
    auto result = aligner.align("GACCTACA", "GATCA");
    EXPECT_TRUE(result.hasAlternative);
    ASSERT_FALSE(result.alternativeAlignments.empty());
    EXPECT_EQ(result.alternativeAlignments[0], "GA-C");
}

TEST(LocalAlignerTest, MatrixDimensions)
{
    auto result = aligner.align("GATTACA", "CGATACAT");
    EXPECT_EQ(result.scoreMatrix.size(), 8);
    EXPECT_EQ(result.scoreMatrix[0].size(), 9);
}

TEST(LocalAlignerTest, ZeroScoreForNonMatchingSequences)
{
    auto result = aligner.align("AAAA", "TTTT");
    EXPECT_EQ(result.maxScore, 0);
    EXPECT_TRUE(result.bestAlignment.empty());
    EXPECT_EQ(result.startPosA.first, 0);
    EXPECT_EQ(result.startPosB.first, 0);
    EXPECT_FALSE(result.hasAlternative);
    EXPECT_TRUE(result.alternativeAlignments.empty());
    EXPECT_EQ(result.alternativeAlignments.size(), 0);
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}