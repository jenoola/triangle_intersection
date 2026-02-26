#include "intersection_private.hpp"

namespace triangle_intersection {
    double get_determinant_3d(double p1[3], double p2[3], double p3[3], double p4[3]) noexcept {
        /*
        a b c
        d e f = aei + bfg + cdh - ceg - bdi - afh
        g h i
        */
        double m[3][3] = {
            { p1[0] - p4[0], p1[1] - p4[1], p1[2] - p4[2] },
            { p2[0] - p4[0], p2[1] - p4[1], p2[2] - p4[2] },
            { p3[0] - p4[0], p3[1] - p4[1], p3[2] - p4[2] }
        };

        return (m[0][0] * m[1][1] * m[2][2]) + (m[0][1] * m[1][2] * m[2][0]) + (m[0][2] * m[1][0] * m[2][1]) 
            - (m[0][2] * m[1][1] * m[2][0]) - (m[0][1] * m[1][0] * m[2][2]) - (m[0][0] * m[1][2] * m[2][1]);
    }

    double get_determinant_2d(double p1[2], double p2[2], double p3[2]) noexcept {
        double m[2][2] = {
            { p1[0] - p3[0], p1[1] - p3[1] },
            { p2[0] - p3[0], p2[1] - p3[1] }
        };

        return m[0][0] * m[1][1] - m[0][1] * m[1][0];
    }

    void reorder_points(double t1[9], double t2[9], double d1[3], double d2[3]) noexcept {

        auto circular_permutation = [] (double t[9], double d[3]) {
            if ((d[1] == 0 && ((d[0] > 0 && d[2] > 0) || (d[0] < 0 && d[2] < 0)))
                || (d[0] == 0 && d[2] == 0)
                || (d[1] != 0 && ((d[0] > 0 && d[2] > 0) || (d[0] < 0 && d[2] < 0)))) {
                std::swap(t[3], t[0]);
                std::swap(t[4], t[1]);
                std::swap(t[5], t[2]);

                std::swap(t[6], t[3]);
                std::swap(t[7], t[4]);
                std::swap(t[8], t[5]);
                return;
            }
            
            if ((d[2] == 0 && ((d[0] > 0 && d[1] > 0) || (d[0] < 0 && d[1] < 0)))
                || (d[0] == 0 && d[1] == 0)
                || (d[2] != 0 && ((d[0] > 0 && d[1] > 0) || (d[0] < 0 && d[1] < 0)))) {
                std::swap(t[3], t[0]);
                std::swap(t[4], t[1]);
                std::swap(t[5], t[2]);

                std::swap(t[6], t[0]);
                std::swap(t[7], t[1]);
                std::swap(t[8], t[2]);
                return;
            }
        };

        circular_permutation(t1, d1);
        circular_permutation(t2, d2);
        
        if (get_determinant_3d(t1 + 6, t1 + 3, t1, t2) < 0) {
            std::swap(t1[3], t1[6]);
            std::swap(t1[4], t1[7]);
            std::swap(t1[5], t1[8]);
        }
        
        if (get_determinant_3d(t2 + 6, t2 + 3, t2, t1) < 0) {
            std::swap(t2[3], t2[6]);
            std::swap(t2[4], t2[7]);
            std::swap(t2[5], t2[8]);
        }
    }

    bool is_point(double t[9]) noexcept {
        return t[0] == t[3] && t[0] == t[6] && t[1] == t[4] && t[1] == t[7] && t[2] == t[5] && t[2] == t[8];
    }

    bool is_segment(double t[9]) noexcept {
        double l1 = get_length(t, t + 3);
        double l2 = get_length(t, t + 6);
        double l3 = get_length(t + 3, t + 6);

        if (l2 == l1 + l3) {
            std::swap(t[3], t[6]);
            std::swap(t[4], t[7]);
            std::swap(t[5], t[8]);

            return true;
        }

        if (l3 == l1 + l2) {
            std::swap(t[0], t[6]);
            std::swap(t[1], t[7]);
            std::swap(t[2], t[8]);

            return true;
        }

        return l1 == l2 + l3;
    }

    double get_length(double p1[3], double p2[3]) noexcept {
        // L == sqrt((x2 - x1) ^ 2 + (y2 - y1) ^ 2 + (z2 - z1) ^ 2)
        return std::sqrt(pow(p2[0] - p1[0], 2.0) + pow(p2[1] - p1[1], 2.0) + pow(p2[2] - p1[2], 2.0));
    };

    void cross_product(double v1[3], double v2[3], double cp[3]) noexcept {
        cp[0] = v1[1] * v2[2] - v1[2] * v2[1];
        cp[1] = v1[2] * v2[0] - v1[0] * v2[2];
        cp[2] = v1[0] * v2[1] - v1[1] * v2[0];
    };

    double dot_product(double v1[3], double v2[3]) noexcept {
        return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
    };

    void make_couterclockwise_2d(double t[6]) noexcept {
        if (get_determinant_2d(t, t + 2, t + 4) < 0) {
            std::swap(t[2], t[4]);
            std::swap(t[3], t[5]);
        }
    }

    int get_dominant_axis(double t[9]) noexcept {
        double n[3] = {
            (t[4] - t[1]) * (t[8] - t[2]) - (t[5] - t[2]) * (t[7] - t[1]),
            (t[5] - t[2]) * (t[6] - t[0]) - (t[3] - t[0]) * (t[8] - t[2]),
            (t[3] - t[0]) * (t[7] - t[1]) - (t[4] - t[1]) * (t[6] - t[0])
        };

        if (abs(n[0]) > abs(n[1])) {
            return abs(n[0]) > abs(n[2]) ? 0 : 2;
        } 

        return abs(n[1]) > abs(n[2]) ? 1 : 2;
    };

    void project_t_2d(double t1[9], double t2[6], int drop) noexcept {
        int j = 0;
        for (int i = 0; i < 9; i++) {
            if (i % 3 == drop) {
                continue;
            }
            t2[j] = t1[i];
            j++;
        }
    }
    
    void project_p_2d(double p1[3], double p2[2], int drop) noexcept {
        if (drop == 0) {
            p2[0] = p1[1];
            p2[1] = p1[2];
        } else  {
            p2[0] = p1[0];
            p2[1] = drop == 1 ? p1[2] : p2[1];
        }
    }
    
    bool is_same_side(double sp1[2], double sp2[2], double p1[2], double p2[2]) noexcept {
        double d1 = get_determinant_2d(sp1, sp2, p1);
        double d2 = get_determinant_2d(sp1, sp2, p2);

        return (d1 < 0 && d2 < 0) || (d1 > 0 && d2 > 0);
    }
}