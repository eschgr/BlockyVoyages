#include "stdafx.h"

#include "VoronoiCell.h"

#include <algorithm>
#include <queue>
#include <list>
#include <assert.h>

#include "../Edge.h"
#include "VoronoiEdge.h"
#include "../Node.h"

namespace BlockyVoyages {
namespace Graph {

class MatchBadEdge {
public:
    bool operator()(Edge* edge) {
        if (!edge) {
            return true;
        }

        VoronoiEdge* vedge = dynamic_cast<VoronoiEdge*>(edge);
        return !vedge->IsValid();
    }
};

class CompareLeftEdges {
public:
    CompareLeftEdges(VoronoiCell* cell) : m_cell(cell) {}

    bool operator()(Edge* edge1, Edge* edge2) {
        // find the top point of the two edges
        const Vector2f& start1 = edge1->GetStart(m_cell)->GetPosition();
        const Vector2f& start2 = edge2->GetStart(m_cell)->GetPosition();

        // the reply is from the viewpoint of the left side. The right side
        // must be negated.
        if (start1.y != start2.y) {
            return start1.y < start2.y;
        } else {
            return start1.x > start2.x;
        }
    }

private:
    VoronoiCell* m_cell;
};

class CompareRightEdges {
public:
    CompareRightEdges(VoronoiCell* cell) : m_cell(cell) {}

    bool operator()(Edge* edge1, Edge* edge2) {
        // find the top point of the two edges
        const Vector2f& start1 = edge1->GetStart(m_cell)->GetPosition();
        const Vector2f& start2 = edge2->GetStart(m_cell)->GetPosition();

        // the reply is from the viewpoint of the left side. The right side
        // must be negated.
        if (start1.y != start2.y) {
            return start1.y > start2.y;
        } else {
            return start1.x < start2.x;
        }
    }

private:
    VoronoiCell* m_cell;
};

void VoronoiCell::RemoveBadEdges() { 
    m_edges.erase(std::remove_if(m_edges.begin(), m_edges.end(), MatchBadEdge()), m_edges.end());
}

void VoronoiCell::Finish(const Vector2f& min_coord, const Vector2f& max_coord,
                  std::vector<Edge*>& all_edges, std::vector<Node*>& all_nodes) {
    std::priority_queue<VoronoiEdge*, std::vector<VoronoiEdge*>, CompareLeftEdges> left_edges(CompareLeftEdges(this));
    std::priority_queue<VoronoiEdge*, std::vector<VoronoiEdge*>, CompareRightEdges> right_edges(CompareRightEdges(this));
    
    for (Edge* edge : m_edges) {
        VoronoiEdge* vedge = dynamic_cast<VoronoiEdge*>(edge);
        Vector2f dir = vedge->GetDirection(this);
        dir.Sub(edge->GetEnd(this)->GetPosition(), edge->GetStart(this)->GetPosition());
        /*LOG(Debugability::DBG_LVL_MESSAGE, "Edge %p: Start: %p (%.10f, %.10f) End: %p (%.10f, %.10f) Dir: (%f, %f) Sub (%f, %f)",
            edge,
            edge->GetStart(this), edge->GetStart(this)->GetPosition().x, edge->GetStart(this)->GetPosition().y,
            edge->GetEnd(this), edge->GetEnd(this)->GetPosition().x, edge->GetEnd(this)->GetPosition().y,
            vedge->GetDirection().x, vedge->GetDirection().y,
            dir.x, dir.y);*/
        if (dir.y < 0.0f) {
            left_edges.push(vedge);
            //LOG(Debugability::DBG_LVL_MESSAGE, "Go Left.");
        } else if (dir.y > 0.0f) {
            right_edges.push(vedge);
            //LOG(Debugability::DBG_LVL_MESSAGE, "Go Right.");
        } else if (dir.x > 0.0f) {
            left_edges.push(vedge);
            //LOG(Debugability::DBG_LVL_MESSAGE, "Go Left.");
        } else {
            right_edges.push(vedge);
            //LOG(Debugability::DBG_LVL_MESSAGE, "Go Right.");
        }
    }

    // clear out the edge list since all the edges are now stored in the queues
    m_edges.clear();

    // add the edges to the left to the back
    while (!left_edges.empty() || !right_edges.empty()) {
        Edge* edge;
        if (left_edges.empty()) {
            edge = right_edges.top();
            right_edges.pop();
        } else {
            edge = left_edges.top();
            left_edges.pop();
        }

        //LOG(Debugability::DBG_LVL_MESSAGE, "Edge %p", edge);

        Node* edge_start = edge->GetStart(this);
        if (!m_edges.empty()) {
            CheckNodesConnected(m_edges.back()->GetEnd(this), edge_start,
                                min_coord, max_coord, all_edges, all_nodes);
        }
        m_edges.push_back(edge);
    }

    // check if the beginning and end are connected
    assert (!m_edges.empty());
    CheckNodesConnected(m_edges.back()->GetEnd(this), m_edges.front()->GetStart(this),
                        min_coord, max_coord, all_edges, all_nodes);

    m_center = (*m_edges.begin())->GetStart(this)->GetPosition();
    m_min_bounds = m_center;
    m_max_bounds = m_center;
    for (auto edge = ++m_edges.begin(); edge != m_edges.end(); ++edge) {
        const Vector2f& start = (*edge)->GetStart(this)->GetPosition();
        m_center += start;
        if (start.x < m_min_bounds.x) {
            m_min_bounds.x = start.x;
        } else if (start.x > m_max_bounds.x) {
            m_max_bounds.x = start.x;
        }
        if (start.y < m_min_bounds.y) {
            m_min_bounds.y = start.y;
        } else if (start.y > m_max_bounds.y) {
            m_max_bounds.y = start.y;
        }
    }
    m_center /= static_cast<float32>(m_edges.size());
}

void VoronoiCell::CheckNodesConnected(Node* start_node, Node* end_node,
                              const Vector2f& min_coord, const Vector2f& max_coord,
                              std::vector<Edge*>& all_edges, std::vector<Node*>& all_nodes) {
    if (start_node == end_node) {
        return;
    }

    // start_pos is not constant because it may change as edges are added.
    Vector2f start_pos = start_node->GetPosition();
    // end_pos should never change as it is the target node.
    const Vector2f& end_pos = end_node->GetPosition();
    assert(start_pos.x == min_coord.x || start_pos.x == max_coord.x ||
           start_pos.y == min_coord.y || start_pos.y == max_coord.y);
    assert(end_pos.x == min_coord.x || end_pos.x == max_coord.x ||
           end_pos.y == min_coord.y || end_pos.y == max_coord.y);

    m_on_boundary = true;

    // check to see if the nodes are on the same edge
    if ((start_pos.x == end_pos.x && (start_pos.x == min_coord.x || start_pos.x == max_coord.x)) ||
        (start_pos.y == end_pos.y && (start_pos.y == min_coord.y || start_pos.y == max_coord.y))) {
        VoronoiEdge* new_edge = new VoronoiEdge(start_node, end_pos - start_pos, this, nullptr);
        new_edge->SetEnd(end_node);
        all_edges.push_back(new_edge);
        m_edges.push_back(new_edge);
    } else {
        // categorize the starting point
        int32 start_edge, end_edge;
        if (start_pos.x == min_coord.x) {
            start_edge = 0;
        } else if (start_pos.y == min_coord.y) {
            start_edge = 1;
        } else if (start_pos.x == max_coord.x) {
            start_edge = 2;
        } else { //if (start_pos.y == max_coord.y) {
            start_edge = 3;
        }

        if (end_pos.x == min_coord.x) {
            end_edge = 0;
        } else if (end_pos.y == min_coord.y) {
            end_edge = 1;
        } else if (end_pos.x == max_coord.x) {
            end_edge = 2;
        } else {//if (end_pos.y == max_coord.y) {
            end_edge = 3;
        }

        while (start_edge != end_edge) {
            Vector2f corner_pos;
            switch (start_edge) {
            case 0:
                corner_pos = min_coord;
                break;
            case 1:
                corner_pos.Set(max_coord.x, min_coord.y);
                break;
            case 2:
                corner_pos = max_coord;
                break;
            case 3:
                corner_pos.Set(min_coord.x, max_coord.y);
                break;
            default:
                assert(0);
            }
            Vector2f dir;
            dir.Sub(corner_pos, start_pos);
            Node* new_node = new Node(corner_pos);
            VoronoiEdge* new_edge = new VoronoiEdge(start_node, dir, this, nullptr);
            new_edge->SetEnd(new_node);

            // keep track of the new edge
            all_edges.push_back(new_edge);
            all_nodes.push_back(new_node);
            m_edges.push_back(new_edge);

            // advance to the next node
            start_node = new_node;
            start_pos = corner_pos;
            start_edge = (start_edge + 1) % 4;
        }

        VoronoiEdge* new_edge = new VoronoiEdge(start_node, end_pos - start_pos, this, nullptr);
        new_edge->SetEnd(end_node);
        all_edges.push_back(new_edge);
        m_edges.push_back(new_edge);
    }
}

}
}