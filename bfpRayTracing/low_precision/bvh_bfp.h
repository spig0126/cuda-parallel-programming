#ifndef BVH_BFP_H
#define BVH_BFP_H

#include <algorithm>

#include "hittable_list_bfp.h"

inline bool box_x_compare(const shared_ptr<hittable> a, const shared_ptr<hittable> b);
inline bool box_y_compare(const shared_ptr<hittable> a, const shared_ptr<hittable> b);
inline bool box_z_compare(const shared_ptr<hittable> a, const shared_ptr<hittable> b);

class bvh_node : public hittable
{
public:
    bvh_node();
    bvh_node(const hittable_list &list, bfpNum time0, bfpNum time1) : bvh_node(1, list.objects, 0, list.objects.size(), time0, time1) {} // 1: root node의 nodeIdx값
    bvh_node(int idx, const std::vector<shared_ptr<hittable>> &src_objects, int start, int end, bfpNum time0, bfpNum time1);

    virtual bool hit(const ray &r, bfpNum t_min, bfpNum t_max, hit_record &rec) const override;
    virtual bool bounding_box(bfpNum time0, bfpNum time1, aabb &output_box) const override;

public:
    shared_ptr<hittable> left;
    shared_ptr<hittable> right;
    aabb box;

    int nodeIdx;
};

bool bvh_node::bounding_box(bfpNum time0, bfpNum time1, aabb &output_box) const
{
    output_box = box;
    return true;
}

bool bvh_node::hit(const ray &r, bfpNum t_min, bfpNum t_max, hit_record &rec) const
{
    if (!box.hit(r, t_min, t_max))
        return false;

    // Node Information
    // point3 pm = box.minimum, pM = box.maximum;

    bool hit_left = left->hit(r, t_min, t_max, rec);
    bool hit_right = right->hit(r, t_min, hit_left ? rec.t : t_max, rec);

    return hit_left || hit_right;
}

bvh_node::bvh_node(int idx,
                   const std::vector<shared_ptr<hittable>> &src_objects,
                   int start, int end, bfpNum time0, bfpNum time1)
{
    // 노드마다 고유번호 부여
    nodeIdx = idx;
    auto objects = src_objects; // Create a modifiable array of the source scene objects
    int axis = 0;               // 일단 x축으로 고정!
    int object_span = end - start;
    auto comparator = (axis == 0)   ? box_x_compare
                      : (axis == 1) ? box_y_compare
                                    : box_z_compare;
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

    // bfp 적용
    aabb box_left, box_right;

    if (!left->bounding_box(time0, time1, box_left) // 자식 노드의 aabb를 box_left에 저장하는 부분
        || !right->bounding_box(time0, time1, box_right))
        std::cerr << "No bounding box in bvh_node constructor.\n";
    box = surrounding_box(box_left, box_right); // 현재 노드의 aabb 생성 fuction
 
    // Node Information
    // point3 pm = box.minimum, pM = box.maximum;
}

inline bool box_compare(const shared_ptr<hittable> a,
                        const shared_ptr<hittable> b, int axis)
{
    aabb box_a;
    aabb box_b;

    if (!a->bounding_box(b_0, b_0, box_a) || !b->bounding_box(b_0, b_0, box_b))
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