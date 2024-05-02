#include "stdafx.h"

#include "Arc.h"

#include <assert.h>

#include "VoronoiEdge.h"
#include "../Node.h"
#include "VoronoiCell.h"
#include "../Graph.h"

namespace BlockyVoyages {
namespace Graph {

Arc::Arc(const Vector2d& site, VoronoiCell* cell)
    : m_site(site),
      m_edge(nullptr),
      m_cell(cell),
      m_parent(nullptr),
      m_right(nullptr),
      m_left(nullptr),
      m_next(nullptr),
      m_prev(nullptr),
      m_event(nullptr)
{}

Arc::Arc(const Vector2d& site, VoronoiEdge* edge, Arc* left, Arc* right)
    : m_site(site),
      m_edge(edge),
      m_parent(nullptr),
      m_left(left),
      m_right(right),
      m_next(nullptr),
      m_prev(nullptr),
      m_event(nullptr)
{
    m_left->m_parent = this;
    m_right->m_parent = this;
}

void Arc::SetLeft(Arc* new_left) {
    m_left = new_left;
    m_left->m_parent = this;
}

void Arc::SetRight(Arc* new_right) {
    m_right = new_right;
    m_right->m_parent = this;
}

Arc* Arc::GetLeftBreakpoint() {
    // by the definition of a tree, a node to the left occurs before this node.
    // if there is a left node, then we can just go left without worrying
    // about finding an ancestor.
    if (nullptr != m_left) {
        return this;
    }

    Arc* ancestor = m_parent;
    Arc* prev = this;
    while (nullptr != ancestor && ancestor->m_left == prev) {
        prev = ancestor;
        ancestor = ancestor->m_parent;
    }

    return ancestor;
}

Arc* Arc::GetRightBreakpoint() {
    // by the definition of a tree, a node to the right occurs before this node.
    // if there is a right node, then we can just go right without worrying
    // about finding an ancestor.
    if (nullptr != m_right) {
        return this;
    }

    Arc* ancestor = m_parent;
    Arc* prev = this;
    while (nullptr != ancestor && ancestor->m_right == prev) {
        prev = ancestor;
        ancestor = ancestor->m_parent;
    }

    return ancestor;
}

double Arc::GetYAtX(double x, double sweep_y) const {
    double site_x = static_cast<double>(m_site.x);
    double site_y = static_cast<double>(m_site.y);

    double x_diff = site_x - static_cast<double>(x);
    double out_y = 0.5 * (site_y * site_y + x_diff * x_diff - sweep_y * sweep_y) / (site_y - sweep_y);
    return out_y;
}

void Arc::RemoveFromTree() {
    assert (IsLeaf());
    // patch up the arc neighbor pointers
    if (m_prev) {
        m_prev->m_next = m_next;
    }
    if (m_next) {
        m_next->m_prev = m_prev;
    }

    Arc* gParent = m_parent->m_parent;
    if (m_parent == gParent->m_left) {
        if (this == m_parent->m_left) {
            gParent->SetLeft(m_parent->m_right);
        } else {
            gParent->SetLeft(m_parent->m_left);
        }
    } else {
        if (this == m_parent->m_left) {
            gParent->SetRight(m_parent->m_right);
        } else {
            gParent->SetRight(m_parent->m_left);
        }
    }
    delete m_parent;
    delete this;
}

Arc* Arc::GetLeftLeaf() {
    Arc* left = this;
	while (!left->IsLeaf()) {
        left = left->m_left;
    }
	return left;
}

Arc* Arc::GetRightLeaf() {
    Arc* right = this;
	while (!right->IsLeaf()) {
        right = right->m_right;
    }
	return right;
}

double Arc::FindBreakpointX(double y) const {
    assert (!IsLeaf());

    Arc* left_arc = m_left->GetRightLeaf();
    Arc* right_arc = m_right->GetLeftLeaf();

    const Vector2d& left_site = left_arc->GetSite();
    const Vector2d& right_site = right_arc->GetSite();

    double result;
    if (left_site.y == right_site.y) {
        result = (left_site.x + right_site.x) / 2;
    } else if (right_site.y == y) {
        result = right_site.x;
    } else if (left_site.y == y) {
        result = left_site.x;
    } else {
        double left_y = left_site.y - y;
        double right_y = right_site.y - y;
        double right_x = right_site.x - left_site.x;

        // Use the quadratic formula.
        double z0 = 2.0f * left_y;
        double z1 = 2.0f * right_y;

        double a = 1.0f / z0 - 1.0f / z1;
        double b = -2.0f * (0.0f / z0 - right_x / z1);
        double c = (left_y * left_y) / z0 -
                    (right_x * right_x + right_y * right_y) / z1;

        double det = b * b - 4.0f * a * c;
        double sqrt_det = sqrt(det);

        double x1 = (-b - sqrt_det) / (2.0f * a);
        double x2 = (-b + sqrt_det) / (2.0f * a);
        if (left_site.y < right_site.y) {
            result = max(x1, x2);
        } else {
            result = min(x1, x2);
        }
        result += left_site.x;
    }
    // Plug back into one of the parabola equations.
    return result;
}

}
}