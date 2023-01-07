/*
 * ===================================================
 *
 *       Filename:  bvh.h
 *    Description:  Ray Tracing: The Next Week (RTTNW): ~BVH
 *        Created:  2022/07/13
 *
 * ===================================================
 */

// Preprocessors

#ifndef BVH_H
#define BVH_H

#include <algorithm>
#include "utility.h"
#include "hittable.h"
#include "hittable_list.h"

// Function Prototypes

inline bool box_x_compare(const shared_ptr<hittable> a, const shared_ptr<hittable> b);
inline bool box_y_compare(const shared_ptr<hittable> a, const shared_ptr<hittable> b);
inline bool box_z_compare(const shared_ptr<hittable> a, const shared_ptr<hittable> b);

// Classes

class bvh_node : public hittable
{
public:
    bvh_node();

    bvh_node(const hittable_list &list, float time0, float time1)
        : bvh_node(1, list.objects, 0, list.objects.size(), time0, time1) // 1: root node의 nodeIdx값
    {
    }

    bvh_node(
        int idx,
        const std::vector<shared_ptr<hittable>> &src_objects,
        int start, int end, float time0, float time1);

    virtual bool hit(
        const ray &r, float t_min, float t_max, hit_record &rec) const override;

    virtual bool bounding_box(float time0, float time1, aabb &output_box) const override;

public:
    shared_ptr<hittable> left;
    shared_ptr<hittable> right;
    aabb box;

    int nodeIdx;
};

bool bvh_node::bounding_box(float time0, float time1, aabb &output_box) const
{
    output_box = box;
    return true;
}

bool bvh_node::hit(const ray &r, float t_min, float t_max, hit_record &rec) const
{
    if (!box.hit(r, t_min, t_max))
        return false;

    // Node Information
    point3 pm = box.minimum, pM = box.maximum;
    // printf("------------------------------ NODE %d HIT! ------------------------------\n", nodeIdx);
    // printf("  -%-8s: (%.1lf,%.1lf,%.1lf) ~ (%.1lf,%.1lf,%.1lf)\n", "AABB", pm.e[0], pm.e[1], pm.e[2], pM.e[0], pM.e[1], pM.e[2]);
    // printf("--------------------------------------------------------------------------\n");

    //	printf("BVH is executed\n");
    bool hit_left = left->hit(r, t_min, t_max, rec);
    bool hit_right = right->hit(r, t_min, hit_left ? rec.t : t_max, rec);

    return hit_left || hit_right;
}

bvh_node::bvh_node(
    int idx,
    const std::vector<shared_ptr<hittable>> &src_objects,
    int start, int end, float time0, float time1)
{

    //	printf("Creating BVH Nodes\n");

    // 노드마다 고유번호 부여
    nodeIdx = idx;

    auto objects = src_objects; // Create a modifiable array of the source scene objects

    int axis = 0; // 일단 x축으로 고정!
                  //    	int axis = random_int(0,2);
    auto comparator = (axis == 0)   ? box_x_compare
                      : (axis == 1) ? box_y_compare
                                    : box_z_compare;

    int object_span = end - start;

    if (object_span == 1)
    {
        left = right = objects[start];
    }
    else if (object_span == 2)
    {
        if (comparator(objects[start], objects[start + 1]))
        {
            left = objects[start];
            right = objects[start + 1];
        }
        else
        {
            left = objects[start + 1];
            right = objects[start];
        }
    }
    else
    {
        std::sort(objects.begin() + start, objects.begin() + end, comparator);

        auto mid = start + object_span / 2;
        left = make_shared<bvh_node>(idx + 1, objects, start, mid, time0, time1);
        right = make_shared<bvh_node>(idx + 2, objects, mid, end, time0, time1);
    }

    aabb box_left, box_right;

    if (!left->bounding_box(time0, time1, box_left) || !right->bounding_box(time0, time1, box_right))
        std::cerr << "No bounding box in bvh_node constructor.\n";
    box = surrounding_box(box_left, box_right);

    // Node Information
    point3 pm = box.minimum, pM = box.maximum;
    // printf("------------------------------ NODE %d CREATED (obj# = %d) ------------------------------\n", nodeIdx, object_span);
    // printf("  -%-8s: (%d ~ %d)\n", "Objects", (int)start, (int)end - 1);
    // printf("  -%-8s: (%.1lf,%.1lf,%.1lf) ~ (%.1lf,%.1lf,%.1lf)\n", "AABB", pm.e[0], pm.e[1], pm.e[2], pM.e[0], pM.e[1], pM.e[2]);
}

inline bool box_compare(const shared_ptr<hittable> a, const shared_ptr<hittable> b, int axis)
{
    aabb box_a;
    aabb box_b;

    if (!a->bounding_box(0, 0, box_a) || !b->bounding_box(0, 0, box_b))
        std::cerr << "No bounding box in bvh_node constructor.\n";

    return box_a.min().e[axis] < box_b.min().e[axis];
}

inline bool box_x_compare(const shared_ptr<hittable> a, const shared_ptr<hittable> b)
{
    return box_compare(a, b, 0);
}

inline bool box_y_compare(const shared_ptr<hittable> a, const shared_ptr<hittable> b)
{
    return box_compare(a, b, 1);
}

inline bool box_z_compare(const shared_ptr<hittable> a, const shared_ptr<hittable> b)
{
    return box_compare(a, b, 2);
}

#endif