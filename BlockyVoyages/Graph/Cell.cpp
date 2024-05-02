#include "stdafx.h"

#include "Cell.h"
#include "Edge.h"
#include "Node.h"

namespace BlockyVoyages {
namespace Graph {

Cell::Cell()
    : m_center(0.0f, 0.0f),
      m_on_boundary(false)
{}

Cell::Cell(const Cell* other)
    : m_edges(other->m_edges),
      m_center(other->m_center),
      m_on_boundary(other->m_on_boundary),
      m_min_bounds(other->m_min_bounds),
      m_max_bounds(other->m_max_bounds)
{}

void Cell::UpdatePointers(const std::map<Edge*, Edge*>& edge_map) {
    for (uint32 ind = 0; ind < m_edges.size(); ++ind) {
        m_edges[ind] = edge_map.find(m_edges[ind])->second;
    }
}

void Cell::CalculateCenterStats() {
    m_center.Zero();
    for (auto* edge : m_edges) {
        m_center += edge->GetStart(this)->GetPosition();
    }
    m_center /= static_cast<float32>(m_edges.size());
}

bool Cell::PointInside(const Vector2f& point) {
    // find how many of the edge segments intersect a line to the positive x
    // direction from the point.
    int32 n_intersects = 0;

    // first perform a bounding box test on the point
    if (point.x < m_min_bounds.x || point.x > m_max_bounds.x) {
        return false;
    }

    if (point.y < m_min_bounds.y || point.y > m_max_bounds.y) {
        return false;
    }

    for (Edge* edge : m_edges) {
        const Vector2f& edge_start = edge->GetStart(this)->GetPosition();
        Vector2f edge_norm = edge->GetSegment(this).GetNormal();

        Vector2f dir_to_edge;
        dir_to_edge.Sub(point, edge_start);

        float dp = dir_to_edge.Dot(edge_norm);

        if (dp <= 0.0f)  {
            return false;
        }
    }
    return true;

    /*
    for (Edge* edge : m_edges) {
        if (start.x > end.x) {
            std::swap(start, end);
        }
        if (point.x > end.x && point.x > start.x) {
            continue;
        }
        if (point.y < end.y && point.y < start.y) {
            continue;
        }
        if (point.y > end.y && point.y > start.y) {
            continue;
        }

        Vector2f edge_dir(end.x - start.x, end.y - start.y);

        if (edge_dir.y == 0.0f) {
            // if this line is horizontal, the only way for it to intersect the
            // ray is for it to have the same y coordinate and be to the right.
            if (start.y == point.y && (start.x > point.x || end.x > point.x)) {
                ++n_intersects;
            }
        } else {
            float32 inv_slope = edge_dir.x / edge_dir.y;
            float32 x_inter = inv_slope * (point.y - start.y) + start.x;

            if (x_inter >= point.x && start.x <= x_inter && x_inter <= end.x) {
                ++n_intersects;
            }
        }
    }

    // if the line intersects an odd number of times, the point is inside.
    if (n_intersects > 2) {
        LOG(Debugability::DBG_LVL_MESSAGE, "N_INTersects: %d", n_intersects);
    }
    return n_intersects % 2;*/
}

Cell* Cell::GetNextCell(const Geometry::Segment& segment) {
    for (Edge* edge : m_edges) {
        const Vector2f& start = edge->GetStart(this)->GetPosition();
        const Vector2f& end = edge->GetEnd(this)->GetPosition();
        Vector2f normal(end.y - start.y, start.x - end.x);
        Geometry::Segment edge_seg = edge->GetSegment(this);
        Vector2f segment_dir;
        segment_dir.Sub(segment.GetEnd(), segment.GetStart());
        if (normal.Dot(segment_dir) < 0) {
            continue;
        }
        Vector2f int_point;
        //return edge->GetOtherCell(this);
        if (segment.Intersects(edge_seg, int_point)) {
            return edge->GetOtherCell(this);
        }
    }
    return nullptr;
}

}
}