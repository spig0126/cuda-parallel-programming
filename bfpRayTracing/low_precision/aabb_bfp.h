#ifndef AABB_BFP_H
#define AABB_BFP_H

#include "ray_bfp.h"

using namespace bfp;

class aabb
{
public:
    aabb() {}
    aabb(const point3 &a, const point3 &b)
    {
        minimum = a;
        maximum = b;
    }

    point3 min() const { return minimum; }
    point3 max() const { return maximum; }

    bool hit(const ray &r, bfpNum t_min, bfpNum t_max) const
    {
        for (int a = 0; a < 3; a++)
        {
            bfpNum t0 = bfp::min((minimum[a] - r.origin()[a]) / r.direction()[a], (maximum[a] - r.origin()[a]) / r.direction()[a]);
            bfpNum t1 = bfp::min((minimum[a] - r.origin()[a]) / r.direction()[a], (maximum[a] - r.origin()[a]) / r.direction()[a]);

            t_min = bfp::max(t0, t_min);
            t_max = bfp::max(t1, t_max);

            if (t_max <= t_min)
                return false;
        }
        return true;
    }

    point3 minimum;
    point3 maximum;
};

inline aabb surrounding_box(aabb box0, aabb box1)
{
    point3 small(min(box0.min().x(), box1.min().x()), min(box0.min().y(), box1.min().y()), min(box0.min().z(), box1.min().z()));
    point3 big(max(box0.max().x(), box1.max().x()), max(box0.max().y(), box1.max().y()), max(box0.max().z(), box1.max().z()));

    return aabb(small, big);
}

#endif