#include <gtest/gtest.h>
#include "intersection.hpp"

using namespace triangle_intersection;

TEST(Intersection, Intersection) {
    double t1[] = { -78, 99, 40, -21, -72, 63, -19, -78, -83 };
    double t2[] = { 9, 5, -21, 96, 77, -51, -95, -1, -16 };

    ASSERT_TRUE(have_intersection(t1, t2));
}

TEST(Intersection, NoIntersection) {
    double t1[] = { -1, 4, 3, 3, 5, -2, -4, -7, 1 };
    double t2[] = { 3, -5, -4, -2.5, -5.7, 0, 6, 1.6, 2 };

    ASSERT_FALSE(have_intersection(t1, t2));
}

TEST(Intersection, NoIntersectionParallel) {
    double t1[] = { -1, 4, 3, 3, 5, -2, -4, -7, 1 };
    double t2[] = { -1, 4, 5, 3, 5, 0, -4, -7, 3 };

    ASSERT_FALSE(have_intersection(t1, t2));
}

TEST(Intersection, IntersectionCoplanar) {
    double t1[] = { 3, 3, 4, 3, -1, 4, 1.5, 1.5, 4 };
    double t2[] = { 10, 5, 4, -1, 0, 4, 7, -3, 4 };

    ASSERT_TRUE(have_intersection(t1, t2));
}

TEST(Intersection, IntersectionCoplanarContainement) {
    double t1[] = { 3, 3, 4, 3, -1, 4, 1.5, 1.5, 4 };
    double t2[] = { 5, 10, 4, -1, 0, 4, 7, -3, 4 };

    ASSERT_TRUE(have_intersection(t1, t2));
}

TEST(Intersection, SharedVertexCoplanar) {
    double t1[] = { 3, 3, 4, 3, -1, 4, 1.5, 1.5, 4 };
    double t2[] = { 5, 10, 4, -1, 0, 4, 1.5, 1.5, 4 };

    ASSERT_TRUE(have_intersection(t1, t2));
}

TEST(Intersection, SharedEdgeCoplanar) {
    double t1[] = { 3, 3, 4, 3, -1, 4, 1.5, 1.5, 4 };
    double t2[] = { 3, 3, 4, -1, 0, 4, 1.5, 1.5, 4 };

    ASSERT_TRUE(have_intersection(t1, t2));
}

TEST(Intersection, NoIntersectionCoplanar) {
    double t1[] = { 3, 3, 4, 3, -1, 4, 1.5, 1.5, 4 };
    double t2[] = { 5, 10, 4, 4, 0, 4, 7, -3, 4 };

    ASSERT_FALSE(have_intersection(t1, t2));
}

TEST(Intersection, SharedVertex) {
    double t1[] = { -1, 4, 3, 3, 5, -2, -4, -7, 1 };
    double t2[] = { 3, -5, -4, 3, 5, -2, 6, 1.6, 2 };

    ASSERT_TRUE(have_intersection(t1, t2));
}

TEST(Intersection, VertexOnEdge) {
    double t1[] = { 10, -25.5, 3, 6.3, -25.5, 1, 5, 0, 0 };
    double t2[] = { 5, -5, -4, 5, 5, 4, 6, 1.6, 2 };

    ASSERT_TRUE(have_intersection(t1, t2));
}

TEST(Intersection, SharedEdge) {
    double t1[] = { 10, -25.5, 3, 12.3, 9.3, -4, 2.5, -2.7, 2 };
    double t2[] = { 5, -5, -4, 12.3, 9.3, -4, 2.5, -2.7, 2 };

    ASSERT_TRUE(have_intersection(t1, t2));
}

TEST(Intersection_Degenerate, IntersectionTrianglePoint) {
    double t1[] = { 3, 3, 4, -3, -3, -4, 3, -1, 4 };
    double t2[9] = { 0, -1, 0, 0, -1, 0, 0, -1, 0 };

    ASSERT_TRUE(have_intersection(t1, t2));
}

TEST(Intersection_Degenerate, NoIntersectionTrianglePoint) {
    double t1[] = { 5.45, -1.77, -0.68, 5.45, -1.77, -0.68, 5.45, -1.77, -0.68 };
    double t2[] = { 5, -5, -4, 12.3, 9.3, -4, 2.5, -2.7, 2 };

    ASSERT_FALSE(have_intersection(t1, t2));
}

TEST(Intersection_Degenerate, IntersectionTriangleSegment) {
    double t1[] = { 3, 3, 4, -3, -3, -4, 3, -1, 4 };
    double t2[] = { -5, 2, -2, -5, 2, -2, 5, -4, 3 };

    ASSERT_TRUE(have_intersection(t1, t2));
}

TEST(Intersection_Degenerate, NoIntersectionTriangleSegment) {
    double t1[] = { 3, 3, 4, -3, -3, -4, 3, -1, 4 };
    double t2[] = { -1, 2, -2, 5, -4, 3, -1, 2, -2 };

    ASSERT_FALSE(have_intersection(t1, t2));
}

TEST(Intersection_Degenerate, NoIntersectionSegmentSegment) {
    double t1[] = { 3, 3, 4, -3, -3, -4, -3, -3, -4 };
    double t2[] = { -1, 2, -2, 5, -4, 3, -1, 2, -2 };

    ASSERT_FALSE(have_intersection(t1, t2));
}

TEST(Intersection_Degenerate, IntersectionSegmentSegment) {
    double t1[] = { 3, 3, 4, -3, -3, -4, -3, -3, -4 };
    double t2[] = { 3, 3, -4, -3, -3, 44, -3, -3, 4 };

    ASSERT_TRUE(have_intersection(t1, t2));
}

TEST(Intersection_Degenerate, NoIntersectionSegmentSegmentCoplanar) {
    double t1[] = { 3, 3, -4, -3, -3, 4, -3, -3, 4 };
    double t2[] = { 4, 4, -4, -2, -2, 4, -2, -2, 4 };

    ASSERT_FALSE(have_intersection(t1, t2));
}

TEST(Intersection_Degenerate, NoIntersectionSegmentSegmentParallel) {
    double t1[] = { 3, 3, -4, -3, -3, 4, -3, -3, 4 };
    double t2[] = { 4, 3, -4, -2, -1, 4, -2, -1, 4 };

    ASSERT_FALSE(have_intersection(t1, t2));
}

TEST(Intersection_Degenerate, IntersectionTriangleSegmentCoplanar) {
    double t1[] = { 10, -10, 3, 12, 15, 3, -5, -3, 3 };
    double t2[] = { 15, 9, 3, -10, 5, 3, 2.5, 7, 3 };

    ASSERT_TRUE(have_intersection(t1, t2));
}

TEST(Intersection_Degenerate, IntersectionSegmentPoint) {
    double t1[] = { 2.5, 7, 3, 2.5, 7, 3, 2.5, 7, 3 };
    double t2[] = { 15, 9, 3, -10, 5, 3, 2.5, 7, 3 };

    ASSERT_TRUE(have_intersection(t1, t2));
}

TEST(Intersection_Degenerate, NoIntersectionSegmentPoint) {
    double t1[] = { 5, 7, 3, 5, 7, 3, 5, 7, 3 };
    double t2[] = { 15, 9, 3, -10, 5, 3, 2.5, 7, 3 };

    ASSERT_FALSE(have_intersection(t1, t2));
}

TEST(Intersection_Degenerate, IntersectionPointPoint) {
    double t1[] = { 5, 7, 3, 5, 7, 3, 5, 7, 3 };
    double t2[] = { 5, 7, 3, 5, 7, 3, 5, 7, 3 };

    ASSERT_TRUE(have_intersection(t1, t2));
}

TEST(Intersection_Degenerate, NoIntersectionPointPoint) {
    double t1[] = { 10, -25.5, 3, 10, -25.5, 3, 10, -25.5, 3 };
    double t2[] = { 5.45, -1.77, -0.68, 5.45, -1.77, -0.68, 5.45, -1.77, -0.68 };

    ASSERT_FALSE(have_intersection(t1, t2));
}