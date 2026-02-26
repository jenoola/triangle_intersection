#include "intersection_private.hpp"

namespace triangle_intersection {
    bool have_intersection(double t1[9], double t2[9]) noexcept {
        if (is_point(t1)) {
            if (is_point(t2)) {
                return have_intersection_p_p(t1, t2);
            }
            if (is_segment(t2)) {
                return have_intersection_s_p(t2, t1);
            }
            return have_intersection_t_p(t2, t1);
        }

        if (is_point(t2)) {
            if (is_segment(t1)) {
                return have_intersection_s_p(t1, t2);
            }
            return have_intersection_t_p(t1, t2);
        }

        if (is_segment(t1)) {
            if (is_segment(t2)) {
                return have_intersection_s_s(t1, t2);
            }
            return have_intersection_t_s(t2, t1);
        }

        if (is_segment(t2)) {
            return have_intersection_t_s(t1, t2);
        }

        double d1[3] = { 
            get_determinant_3d(t2, t2 + 3, t2 + 6, t1), 
            get_determinant_3d(t2, t2 + 3, t2 + 6, t1 + 3), 
            get_determinant_3d(t2, t2 + 3, t2 + 6, t1 + 6)
        };

        if (d1[0] == 0 && d1[1] == 0 && d1[2] == 0) {
            return have_intersection_coplanar_t_t(t1, t2);
        }

        if ((d1[0] < 0 && d1[1] < 0 && d1[2] < 0) || (d1[0] > 0 && d1[1] > 0 && d1[2] > 0)) {
            return false;
        }

        double d2[3] = {
            get_determinant_3d(t1, t1 + 3, t1 + 6, t2), 
            get_determinant_3d(t1, t1 + 3, t1 + 6, t2 + 3),
            get_determinant_3d(t1, t1 + 3, t1 + 6, t2 + 6)
        };

        if ((d2[0] < 0 && d2[1] < 0 && d2[2] < 0) || (d2[0] > 0 && d2[1] > 0 && d2[2] > 0)) {
            return false;
        }

        reorder_points(t1, t2, d1, d2);

        double d[2] = {
            get_determinant_3d(t1, t1 + 3, t2, t2 + 3),
            get_determinant_3d(t1, t1 + 6, t2 + 6, t2)
        };

        return d[0] >= 0 && d[1] >= 0;
    }

    bool have_intersection_t_s(double t[9], double s[6]) noexcept {
        double edge1[3] = { t[3] - t[0], t[4] - t[1], t[5] - t[2] };
        double edge2[3] = { t[6] - t[0], t[7] - t[1], t[8] - t[2] };

        double dir[3] = {
            s[3] - s[0],
            s[4] - s[1],
            s[5] - s[2]
        };

        double v1[3];
        cross_product(dir, edge2, v1);

        double det = dot_product(edge1, v1);
        if (abs(det) < EPS) {
            double cp[3];
            cross_product(edge1, edge2, cp);
            
            if (abs(dot_product(dir, cp)) >= EPS) {
                return false;
            }

            double st3[6] = { t[0], t[1], t[2], t[6], t[7], t[8] };

            return have_intersection_s_s(t, s) || have_intersection_s_s(t + 3, s) || have_intersection_s_s(st3, s);
        }

        double invDet = 1.0 / det;

        double v2[3] = {
            s[0] - t[0],
            s[1] - t[1],
            s[2] - t[2]
        };

        double u = dot_product(v2, v1) * invDet;
        if (u < -EPS || u > 1 + EPS) {
            return false;
        }

        double v3[3];
        cross_product(v2, edge1, v3);

        double v = dot_product(dir, v3) * invDet;
        if (v < -EPS || u + v > 1 + EPS) {
            return false;
        }

        double r = dot_product(edge2, v3) * invDet;

        if (r < -EPS || r > 1 + EPS) {
            return false;
        }

        double cross_point[3] = { s[0] + dir[0] * r, s[1] + dir[1] * r, s[2] + dir[2] * r };
        return have_intersection_s_p(s, cross_point);
    }
    
    bool have_intersection_s_s(double s1[6], double s2[6]) noexcept {
        double dir1[3] = { s1[3] - s1[0], s1[4] - s1[1], s1[5] - s1[2] };
        double dir2[3] = { s2[3] - s2[0], s2[4] - s2[1], s2[5] - s2[2] };
        double v[3] = { s1[0] - s2[0], s1[1] - s2[1], s1[2] - s2[2] };

        double a = dot_product(dir1, dir1); // always >= 0
        double b = dot_product(dir1, dir2);
        double c = dot_product(dir2, dir2); // always >= 0
        double d = dot_product(dir1, v);
        double e = dot_product(dir2, v);
        double det = a * c - b * b;

        double coef1, coef2;
        
        if (det < EPS) {
            // lines are almost parallel
            coef1 = 0.0;
            coef2 = (b > c ? d / b : e / c);
        } else {
            coef1 = (b * e - c * d) / det;
            coef2 = (a * e - b * d) / det;
        }

        if (coef1 < 0 || coef1 > 1 || coef2 < 0 || coef2 > 1) {
            return false;
        }

        double dist[3] = {
            v[0] + coef1 * dir1[0] - coef2 * dir2[0],
            v[1] + coef1 * dir1[1] - coef2 * dir2[1],
            v[2] + coef1 * dir1[2] - coef2 * dir2[2]
        };

        return dot_product(dist, dist) <= EPS * EPS;
    }

    bool have_intersection_t_p(double t[9], double p[3]) noexcept {
        double v1[3] = { t[3] - t[0], t[4] - t[1], t[5] - t[2] };
        double v2[3] = { t[6] - t[0], t[7] - t[1], t[8] - t[2] };
        double u[3] = { p[0] - t[0], p[1] - t[1], p[2] - t[2] };

        double cp[3];
        cross_product(v1, v2, cp);
        
        if (abs(dot_product(u, cp)) >= EPS) {
            return false;
        }

        int k = get_dominant_axis(t);
        double t1[6];
        project_t_2d(t, t1, k);
        double p1[2];
        project_p_2d(p, p1, k);

        return have_intersection_t_p_2d(t1, p1);
    }

    bool have_intersection_s_p(double s[6], double p[3]) noexcept {
        return get_length(s, s + 3) == get_length(s, p) + get_length(p, s + 3); 
    }

    bool have_intersection_p_p(double p1[3], double p2[3]) noexcept {
        return p1[0] == p2[0] && p1[1] == p2[1] && p1[2] == p2[2];
    }
    
    bool have_intersection_coplanar_t_t(double t1[9], double t2[9]) noexcept {
        int k = get_dominant_axis(t1);
        double t11[6], t21[6];
        project_t_2d(t1, t11, k);
        project_t_2d(t2, t21, k);

        make_couterclockwise_2d(t11);
        make_couterclockwise_2d(t21);

        double p1[2] = { t11[0], t11[1] };
        double q1[2] = { t11[2], t11[3] };
        double r1[2] = { t11[4], t11[5] };
        double p2[2] = { t21[0], t21[1] };
        double q2[2] = { t21[2], t21[3] };
        double r2[2] = { t21[4], t21[5] };

        double d[3] = {
            get_determinant_2d(p2, q2, p1),
            get_determinant_2d(q2, r2, p1),
            get_determinant_2d(r2, p2, p1)
        };

        int count_zero = 0;
        int count_pos = 0;
        for (auto d : d) {
            if (d == 0) {
                count_zero++;
            } else if (d > 0) {
                count_pos++;
            }
        };

        if (count_pos == 3) {
            return true;
        }

        if (count_zero == 2) {
            return true;
        }

        if (count_pos == 2 && count_zero == 1) {
            return true;
        }

        while (!((d[0] > 0 && d[1] >= 0 && d[2] <= 0) || (d[0] > 0 && d[1] <= 0 && d[2] <= 0))) {
            std::swap(q2[0], p2[0]);
            std::swap(q2[1], p2[1]);
            std::swap(d[1], d[0]);

            std::swap(r2[0], q2[0]);
            std::swap(r2[1], q2[1]);
            std::swap(d[2], d[1]);
        }

        if (d[0] > 0 && d[1] >= 0 && d[2] <= 0) {
            if (get_determinant_2d(r2, p2, q1) >= 0) {
                if (get_determinant_2d(r2, p1, q1) < 0
                    || get_determinant_2d(p1, p2, q1) >= 0
                    || get_determinant_2d(p1, p2, r1) < 0) {
                    return true;
                }

                return get_determinant_2d(q1, r1, p2) >= 0;
            } else {
                if (get_determinant_2d(r2, p2, r1) < 0
                    || get_determinant_2d(q1, r1, r2) < 0) {
                    return false;
                }

                return get_determinant_2d(p1, p1, r1) <= 0;
            }
        }
        
        if (get_determinant_2d(r2, p2, q1) >= 0) {
            if (get_determinant_2d(q2, r2, q1) >= 0) {
                if (get_determinant_2d(p1, p2, q1) >= 0) {
                    return get_determinant_2d(p1, q2, q1) <= 0;
                } else {
                    if (get_determinant_2d(p1, p2, r1) < 0) {
                        return false;
                    }
                    return get_determinant_2d(r2, p2, r1) >= 0;
                }
            } else {
                if (get_determinant_2d(p1, q2, q1) > 0
                    || get_determinant_2d(q2, r2, r1) < 0) {
                    return false;
                }

                return get_determinant_2d(q1, r1, q2) >= 0;
            }
        }

        if (get_determinant_2d(r2, p2, r1) < 0) {
            return false;
        }

        if (get_determinant_2d(q1, r1, r2) >= 0) {
            return get_determinant_2d(r1, p1, p2) >= 0;
        }

        if (get_determinant_2d(q1, r1, q2) >= 0) {
            return get_determinant_2d(q2, r2, r1) >= 0;
        }

        return false;
    }

    bool have_intersection_s_s_2d(double s1[4], double s2[4]) noexcept {
        double d1[2] = {
            get_determinant_2d(s1, s1 + 2, s2),
            get_determinant_2d(s1, s1 + 2, s2 + 2)
        };
        double d2[2] = {
            get_determinant_2d(s2, s2 + 2, s1),
            get_determinant_2d(s2, s2 + 2, s1 + 2)
        };

        return ((d1[0] > 0 && d1[1] < 0) || (d1[0] < 0 && d1[1] > 0)
                && (d2[0] > 0 && d2[1] < 0) || (d2[0] < 0 && d2[1] > 0));
    }
    
    bool have_intersection_t_p_2d(double t[6], double p[2]) noexcept {
        return is_same_side(t, t + 2, t + 4, p) 
            && is_same_side(t + 2, t + 4, t, p) 
            && is_same_side(t + 4, t, t + 2, p);
    }
}