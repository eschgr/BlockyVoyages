#include "stdafx.h"

#include "VoronoiEdge.h"
#include "../Node.h"
#include "../Cell.h"

#include <assert.h>

namespace BlockyVoyages {
namespace Graph {

VoronoiEdge::VoronoiEdge(Node* start_node, const Vector2f& dir, Cell* left, Cell* right)
    : Edge(start_node, nullptr, left, right),
      m_dir(dir),
      m_valid(true)
{
    m_dir.Normalize();
}

// Get the edge direction with the given cell on the left side. The cell must be
// either the left or right cell.
Vector2f VoronoiEdge::GetDirection(Cell* cell) {
    assert(cell == m_right_cell || cell == m_left_cell);

    if (cell == m_left_cell) {
        return m_dir;
    }
    return -m_dir;
}

}
}