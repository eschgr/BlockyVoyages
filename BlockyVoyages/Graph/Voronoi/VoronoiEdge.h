#if !defined(__VORONOI_EDGE_H__)
#define __VORONOI_EDGE_H__

#include "../Edge.h"

#include "../../Math/MathMain.h"

namespace BlockyVoyages {
namespace Graph {

class Node;
class Cell;

// An edge for a double connected graph
class VoronoiEdge : public Edge {
public:
    VoronoiEdge(Node* start_node, const Vector2f& dir, Cell* left, Cell* right);
    
    const Vector2f& GetDirection() const { return m_dir; }
    bool IsValid() const { return m_valid; }

    // Get the edge direction in relation to the given cell. The cell must be
    // either the left or right cell.
    Vector2f GetDirection(Cell* cell);

    void Invalidate() { m_valid = false; }

    bool CellOnRight(Cell* cell) { return cell == m_right_cell; }

protected:
    // the direction of the edge.
    Vector2f m_dir;
    bool m_valid;
};

}
}

#endif // __EDGE_H__