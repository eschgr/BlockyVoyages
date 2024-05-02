#include "stdafx.h"

#include "Edge.h"
#include "Node.h"
#include "Cell.h"

#include <assert.h>

namespace BlockyVoyages {
namespace Graph {
Edge::Edge(Node* start_node, Node* end_node, Cell* left, Cell* right)
    : m_start(start_node),
      m_end(end_node),
      m_left_cell(left),
      m_right_cell(right),
      m_segment(start_node->GetPosition(), end_node ? end_node->GetPosition() : 
                                            Vector2f(std::numeric_limits<float>::infinity(),
                                                     std::numeric_limits<float>::infinity())),
      m_reverse_segment(end_node ? end_node->GetPosition() : Vector2f(std::numeric_limits<float>::infinity(),
                                                                      std::numeric_limits<float>::infinity()),
                        start_node->GetPosition())
{
    start_node->AddEdge(this);
    if (m_end) {
        end_node->AddEdge(this);
    }
}

void Edge::SetStart(Node* node) {
    if (m_start) {
        m_start->RemoveEdge(this);
    }
    m_start = node;
    node->AddEdge(this);
    if (m_start && m_end) {
        BuildSegment();
    }
}

void Edge::SetEnd(Node* node) {
    if (m_end) {
        m_end->RemoveEdge(this);
    }
    m_end = node;
    node->AddEdge(this);
    if (m_start && m_end) {
        BuildSegment();
    }
}

// get the start node when the given cell is on the left. The cell must be
// either the left or right cell.
Node* Edge::GetStart(Cell* cell) {
    assert(cell == m_right_cell || cell == m_left_cell);

    if (cell == m_left_cell) {
        return m_start;
    }
    return m_end;
}

// get the end node when the given cell is on the left. The cell must be
// either the left or right cell.
Node* Edge::GetEnd(Cell* cell) {
    assert(cell == m_right_cell || cell == m_left_cell);

    if (cell == m_left_cell) {
        return m_end;
    }
    return m_start;
}

Cell* Edge::GetOtherCell(Cell* cell) {
    assert(cell == m_right_cell || cell == m_left_cell);

    if (cell == m_left_cell) {
        return m_right_cell;
    }
    return m_left_cell;
}

Node* Edge::GetOtherNode(Node* node) {
    assert (node == m_start || node == m_end);

    if (node == m_start) {
        return m_end;
    }
    return m_start;
}

void Edge::UpdatePointers(const std::map<Node*, Node*>& node_map, const std::map<Cell*, Cell*>& cell_map) {
    if (m_start) {
        m_start = node_map.find(m_start)->second;
    }
    if (m_end) {
        m_end = node_map.find(m_end)->second;
    }
    if (m_left_cell) {
        m_left_cell = cell_map.find(m_left_cell)->second;
    }
    if (m_right_cell) {
        m_right_cell = cell_map.find(m_right_cell)->second;
    }
    if (m_start && m_end) {
        BuildSegment();
    }
}

void Edge::ReplaceNode(Node* old_node, Node* new_node) {
    assert (old_node == m_start || old_node == m_end);
    if (old_node == m_start) {
        SetStart(new_node);
    } else {
        SetEnd(new_node);
    }
    if (m_start && m_end) {
        BuildSegment();
    }
}

const Segment& Edge::GetSegment(Cell* cell) const {
    assert(cell == m_right_cell || cell == m_left_cell);

    if (cell == m_left_cell) {
        return m_segment;
    }
    return m_reverse_segment;
}

const Segment& Edge::GetSegment(Node* node) const {
    assert(node == m_start || node == m_end);

    if (node == m_start) {
        return m_segment;
    }
    return m_reverse_segment;
}

void Edge::BuildSegment() {
    m_segment.SetEnds(m_start->GetPosition(), m_end->GetPosition());
    m_reverse_segment.MirrorOf(m_segment);
}

}
}