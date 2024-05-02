#if !defined(__ARC_H__)
#define __ARC_H__

#include "../../Math/MathMain.h"

namespace BlockyVoyages {
namespace Graph {

class VoronoiEdge;
class Node;
class Event;
class VoronoiCell;

class Arc {
public:
    Arc(const Vector2d& site, VoronoiCell* cell);
    Arc(const Vector2d& site, VoronoiEdge* edge, Arc* left, Arc* right);

    bool IsLeaf() const { return (nullptr == m_right && nullptr == m_left); }

    void SetLeft(Arc* new_left);
    void SetRight(Arc* new_right);

    Arc* GetParent() const { return m_parent; }
    Arc* GetLeft() const { return m_left; }
    Arc* GetRight() const { return m_right; }

    const Vector2d& GetSite() const { return m_site; }
    Event* GetEvent() const { return m_event; }
    VoronoiEdge* GetEdge() const { return m_edge; }
    VoronoiCell* GetCell() const { return m_cell; }

    void SetEdge(VoronoiEdge* edge) { m_edge = edge; }
    void SetEvent(Event* event) { m_event = event; }
    void SetCell(VoronoiCell* cell) { m_cell = cell; }

    ///////////////////////////////////////////////////////////////////////////
    // The following methods are valid only for leaves.
    ///////////////////////////////////////////////////////////////////////////
    double GetYAtX(double x, double sweep_y) const;

    void RemoveFromTree();

    Arc* GetLeftBreakpoint();
    Arc* GetRightBreakpoint();

    Arc* GetPrevious() const { return m_prev; }
    Arc* GetNext() const { return m_next; }

    void SetNext(Arc* new_next) { m_next = new_next; }
    void SetPrev(Arc* new_prev) { m_prev = new_prev; }

    ///////////////////////////////////////////////////////////////////////////
    // The following methods are only valid for nodes that are not leaves
    ///////////////////////////////////////////////////////////////////////////

    // Get the leftmost leaf of the subtree of the tree at this node
    Arc* GetLeftLeaf();
    // Get the rightmost leaf of the subtree of the tree at this node
    Arc* GetRightLeaf();

    // Find the intersection point of this arc. This is only valid for an
    // internal node of the tree. Leaf arcs only represent a single
    // parabola, so they cannot have a breakpoint
    double FindBreakpointX(double y) const;

private:
    Vector2d m_site;
    
    // the edge this arc is building up (for nodes).
    VoronoiEdge* m_edge;

    // The cell this arc is covering (for leaves).
    VoronoiCell* m_cell;

    // Tree pointers
    Arc* m_parent;
    Arc* m_left;
    Arc* m_right;

    Arc* m_next;
    Arc* m_prev;

    Event* m_event;
};

}
}

#endif