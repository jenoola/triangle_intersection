#include <utility>
#include "intersection.hpp"

namespace triangle_intersection {
    constexpr double EPS = 1e-12;

    bool have_intersection_t_p(double t[9], double p[3]);
    bool have_intersection_t_s(double t[9], double s[6]);
    bool have_intersection_s_s(double s1[6], double s2[6]);
    bool have_intersection_s_p(double s[6], double p[3]);
    bool have_intersection_p_p(double p1[3], double p2[3]);
    bool have_intersection_coplanar_t_t(double t1[9], double t2[9]);
    
    double get_determinant_3d(double p1[3], double p2[3], double p3[3], double p4[3]);
    void reorder_points(double t1[9], double t2[9], double d1[3], double d2[3]);
    bool is_point(double t[9]);
    bool is_segment(double t[9]);
    void cross_product (double v1[3], double v2[3], double cp[3]);
    double dot_product(double v1[3], double v2[3]);
    double get_length(double p1[3], double p2[3]);

    void make_couterclockwise_2d(double t[6]);
    int get_dominant_axis(double t[9]);
    void project_t_2d(double t1[9], double t2[6], int drop);
    void project_p_2d(double p1[3], double p2[2], int drop);
    bool have_intersection_s_s_2d(double s1[4], double s2[4]);
    double get_determinant_2d(double p1[2], double p2[2], double p3[2]);
    bool have_intersection_t_p_2d(double t[6], double p[2]);
    bool is_same_side(double sp1[2], double sp2[2], double p1[2], double p2[2]);
}