#if !defined(__VORONOI_GRAPH_H__)
#define __VORONOI_GRAPH_H__

#include <queue>
#include <vector>
#include <set>

#include "../../Math/MathMain.h"

#include "../Graph.h"
#include "VoronoiCell.h"
#include "VoronoiEdge.h"

namespace BlockyVoyages {

namespace Graphics {
class Shader;
}

namespace Graph {

class Arc;
class Node;

class Event {
public:
    class CmpEventsFunctor {
    public:
        bool operator()(Event* e1, Event* e2) {
            return ((e1->m_y < e2->m_y) ||
                    (e1->m_y == e2->m_y && e1->m_site.x > e2->m_site.x));
        }
    };

    Event(const Vector2d& site)
        : m_site(site),
          m_arc(nullptr),
          m_valid(true),
          m_y(m_site.y)
    {}

    Event(const Vector2d& site, Arc* arc, double y)
        : m_site(site),
          m_arc(arc),
          m_valid(true),
          m_y(y)
    {}

    const Vector2d& GetSite() const { return m_site; }
    const Arc* GetArc() const { return m_arc; }
    bool IsCircleEvent() const { return nullptr != m_arc; }

    Vector2d m_site;
    double m_y;
    Arc* m_arc;
    bool m_valid;
};

class Voronoi : public Graph {
public:
    Voronoi(const std::vector<Vector2f>& points,
            const Vector2f& min_bounds,
            const Vector2f& max_bounds);
    ~Voronoi();

    // A voronoi diagram can't be copied into
    void CopyFrom(Graph* other);

private:
    enum RegionCodes {
        INSIDE = 0,
        LEFT = 1,
        RIGHT = 2,
        BELOW = 4,
        ABOVE = 8
    };

    bool ProcessNextEvent();

    // The handlers for the two possible types of events
    void HandleSiteEvent(const Vector2d& site);
    void HandleCircleEvent(Event* event);

    Arc* FindArcUnderSite(const Vector2d& site);
    bool CheckForCircleEvent(Arc* arc, double y);
    bool FindCircleSite(const Vector2d& p0, const Vector2d& p1, const Vector2d& p2, double& max_y, Vector2d& circle_site);

    int32 ComputeRegionCode(const Vector2f& pos) const;

    void FinishArc(Arc* arc = nullptr);
    // Cleans up the nodes which only have two edges.
    void RemoveMidNodes();
    // Clips the edges to the bounding area. This does not add the outer edges.
    void ClipEdges();
    // Find the cells in the graph. This will also close any open edges on the boundary.
    void FindCells();

    std::priority_queue<Event*, std::vector<Event*>, Event::CmpEventsFunctor> m_events;

    Arc* m_root;

    double m_sweep_y;
    Vector2d m_last_site;
};

}
}

#endif