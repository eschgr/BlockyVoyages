#if !defined(__EDGE_H__)
#define __EDGE_H__

#include "../Geometry/Segment.h"
#include "../Math/MathMain.h"

#include <map>

namespace BlockyVoyages {
namespace Graph {

class Node;
class Cell;

using Geometry::Segment;

// An edge for a double connected graph
class Edge {
public:
    Edge(Node* start_node, Node* end_node, Cell* left, Cell* right);
    Edge(const Edge* other)
        : m_start(other->m_start),
          m_end(other->m_end),
          m_left_cell(other->m_left_cell),
          m_right_cell(other->m_right_cell),
          m_segment(other->m_segment),
          m_reverse_segment(other->m_reverse_segment)
    {}

    virtual ~Edge() {}
    
    void SetStart(Node* node);
    void SetEnd(Node* node);

    // Get the node marked "start" without considering the cell.
    Node* GetStart() const { return m_start; }
    // Get the node marked "end" without considering the cell.
    Node* GetEnd() const { return m_end; }

    // get the start node when the given cell is on the left. The cell must be
    // either the left or right cell.
    Node* GetStart(Cell* cell);

    // get the end node when the given cell is on the left. The cell must be
    // either the left or right cell.
    Node* GetEnd(Cell* cell);

    Cell* GetOtherCell(Cell* cell);
    Node* GetOtherNode(Node* node);

    const Segment& GetSegment() const { return m_segment; }
    const Segment& GetSegment(Cell* cell) const;
    const Segment& GetSegment(Node* node) const;

    float32 Length() { return m_segment.Length(); }
    float32 LengthSquared() { return m_segment.LengthSquared(); }

    void UpdatePointers(const std::map<Node*, Node*>& node_map, const std::map<Cell*, Cell*>& cell_map);
    void ReplaceNode(Node* old_node, Node* new_node);

protected:
    Node* m_start;
    Node* m_end;

    Cell* m_left_cell;
    Cell* m_right_cell;

    Geometry::Segment m_segment;
    Geometry::Segment m_reverse_segment;

    void BuildSegment();
};

}
}

#endif // __EDGE_H__