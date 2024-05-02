#if !defined(__CELL_H__)
#define __CELL_H__

#include <vector>
#include <map>

#include "../Geometry/Segment.h"
#include "../Math/MathMain.h"

namespace BlockyVoyages {
namespace Graph {

class VoronoiEdge;
class Edge;
class Node;

class Cell {
public:
    Cell();
    Cell(const Cell* other);
    virtual ~Cell() {}

    void AddEdge(Edge* edge) { m_edges.push_back(edge); }
    const std::vector<Edge*>& GetEdges() const { return m_edges; }
    const Vector2f& GetCenter() const { return m_center; }
    void GetBounds(Vector2f& min_bounds, Vector2f& max_bounds) {
        min_bounds = m_min_bounds;
        max_bounds = m_max_bounds;
    }

    void UpdatePointers(const std::map<Edge*, Edge*>& edge_map);
    bool AtBoundary() const { return m_on_boundary; }

    virtual void CalculateCenterStats();

    bool PointInside(const Vector2f& point);
    // gets the edge for which the given segment is leaving this cell. If no
    // such edge exists (the cell is on the end point), returns null.
    Cell* GetNextCell(const Geometry::Segment& segment);

protected:
    std::vector<Edge*> m_edges;
    Vector2f m_center;

    // optimization structures. Can be calculated on the fly if needed
    bool m_on_boundary;
    Vector2f m_min_bounds;
    Vector2f m_max_bounds;
};

}
}

#endif