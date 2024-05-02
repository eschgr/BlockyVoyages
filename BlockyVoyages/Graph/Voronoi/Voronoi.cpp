#include "stdafx.h"

#include "Voronoi.h"

#include "Arc.h"
#include "VoronoiCell.h"
#include "VoronoiEdge.h"
#include "../Node.h"

#include "../../Graphics/Shader.h"
#include "../../Graphics/VertexArrayObject.h"

#include <algorithm>

#include <assert.h>

namespace BlockyVoyages {
namespace Graph {

class MatchBadGraphElement {
public:
    bool operator()(Node* value) {
        if (!value) {
            return true;
        }

        if (value->Degree() == 0) {
            delete value;
            return true;
        }
        return false;
    }

    bool operator()(Edge* edge) {
        if (!edge) {
            return true;
        }
        VoronoiEdge* vedge = dynamic_cast<VoronoiEdge*>(edge);
        if (!vedge->IsValid()) {
            delete edge;
            return true;
        }
        return false;
    }
};

Voronoi::Voronoi(const std::vector<Vector2f>& points,
                 const Vector2f& min_bounds,
                 const Vector2f& max_bounds)
    : Graph(min_bounds, max_bounds),
      m_root(nullptr),
      m_sweep_y(max_bounds.y)
{
    assert(points.size() >= 2);
    
    // Filling the event queue with the event points
    for (auto& point : points) {
        m_events.push(new Event(Vector2d(point.x, point.y)));
    }
    
    // Process all the events.
    while (ProcessNextEvent());

    // Finish up the points
    FinishArc(m_root);
    RemoveMidNodes();
    ClipEdges();

    // clean up a bit.
    for (auto* cell : m_cells) {
        VoronoiCell* vcell = dynamic_cast<VoronoiCell*>(cell);
        vcell->RemoveBadEdges();
        vcell->Finish(m_min_bounds, m_max_bounds, m_edges, m_nodes);
    }
    m_edges.erase(std::remove_if(m_edges.begin(), m_edges.end(), MatchBadGraphElement()), m_edges.end());
    m_nodes.erase(std::remove_if(m_nodes.begin(), m_nodes.end(), MatchBadGraphElement()), m_nodes.end());
}

Voronoi::~Voronoi() {
}

// A voronoi diagram can't be copied into
void Voronoi::CopyFrom(Graph* other) { 
    assert(0);
}

bool Voronoi::ProcessNextEvent() {
    if (m_events.empty()) {
        return false;
    }

    Event* event = m_events.top();
    m_events.pop();

    m_sweep_y = event->m_y;

    if (event->m_valid) {
        if (event->IsCircleEvent()) {
            HandleCircleEvent(event);
        } else if (event->GetSite() != m_last_site) {
            HandleSiteEvent(event->GetSite());
            m_last_site = event->m_site;
        }
    }
    delete event;
    return true;
}

void Voronoi::FinishArc(Arc* arc) {
    if (nullptr == m_root) {
        return;
    }
    if (nullptr == arc) {
        FinishArc(m_root);
        m_root = nullptr;
        return;
    }
    if(arc->IsLeaf()) {
        delete arc;
        return;
    }

    VoronoiEdge* arc_edge = dynamic_cast<VoronoiEdge*>(arc->GetEdge());
    const Vector2f& start_pos = arc_edge->GetStart()->GetPosition();
    // Skip any edge that is outside the boundary. Edges that start outside the
    // boundary will only be set to the boundary and will not enter the region.
    // Still set the end to the start, though. This will be removed at a later
    // step.
    int region_code = ComputeRegionCode(start_pos);
    const Vector2f& arc_dir = arc_edge->GetDirection();
    if (((region_code & LEFT) && arc_dir.x < 0) ||
        ((region_code & RIGHT) && arc_dir.x > 0) ||
        ((region_code & BELOW) && arc_dir.y < 0) ||
        ((region_code & ABOVE) && arc_dir.y > 0)) {
        // this edge should just be tossed out as it will never affect the
        // graph after clipping and it may leave a hanging edge at the edge
        // of the boundary.
        // Let the clipping handle it, though, as removing the edge will cause
        // the node to look like a mid-point node with degree 2.s
        Node* new_node = new Node(arc_edge->GetStart()->GetPosition());
        arc_edge->SetEnd(new_node);
        m_nodes.push_back(new_node);
    } else {
        Vector2f end_pos;
        if (fabs(arc_dir.x) < fabs(arc_dir.y)) {
            if (arc_dir.y > 0.0f) {
                end_pos.y = m_max_bounds.y;
            } else {
                end_pos.y = m_min_bounds.y;
            }
            end_pos.x = arc_dir.x * (end_pos.y - start_pos.y) / arc_dir.y + start_pos.x;
        } else if (arc_dir.x != 0.0f) {
            if(arc_dir.x > 0.0f) {
                end_pos.x = m_max_bounds.x;
            } else {
                end_pos.x = m_min_bounds.x;
            }
            end_pos.y = arc_dir.y * (end_pos.x - start_pos.x) / arc_dir.x + start_pos.y;
        } else {
            end_pos = start_pos;
        }
        assert(nullptr != arc_edge);
        assert(nullptr == arc_edge->GetEnd());
        Node* new_node = new Node(end_pos);
        arc_edge->SetEnd(new_node);
        m_nodes.push_back(new_node);
    }

    FinishArc(arc->GetLeft());
    FinishArc(arc->GetRight());
    delete arc;
}

void Voronoi::RemoveMidNodes() {
    // remove those nodes with degree two. Since the Voronoi diagram can only
    // have straight edges except where many edges meet at a single node, edges
    // with degree 2 fall on a straight edge
    for (uint32 node_ind = 0; node_ind < m_nodes.size(); ++node_ind) {
        Node* node = m_nodes[node_ind];
        if (node->Degree() == 2) {
            VoronoiEdge* edge0 = dynamic_cast<VoronoiEdge*>(node->GetEdges()[0]);
            VoronoiEdge* edge1 = dynamic_cast<VoronoiEdge*>(node->GetEdges()[1]);

            assert (edge0->GetStart() == edge1->GetStart());

            edge0->SetStart(edge1->GetEnd());
            edge0->GetStart()->RemoveEdge(edge1);
            edge1->Invalidate();
            delete node;
            m_nodes[node_ind] = nullptr;
        }
    }
}

void Voronoi::ClipEdges() {
    std::vector<VoronoiEdge*> new_edges;
    for (auto* edge : m_edges) {
        VoronoiEdge* vedge = dynamic_cast<VoronoiEdge*>(edge);
        if (!vedge->IsValid()) {
            continue;
        }
        
        Vector2f start_pos = vedge->GetStart()->GetPosition();
        Vector2f end_pos = vedge->GetEnd()->GetPosition();

        if (start_pos.DistanceSquared(end_pos) <= 1e-7f) {
            // if the start and end positions are close togeth, collapse this edge.
            Node* end = vedge->GetEnd();
            Node* start = vedge->GetStart();
            end->RemoveEdge(vedge);
            start->RemoveEdge(vedge);
            // move all the edges from the end node to the start node
            std::vector<Edge*> edges = end->GetEdges();
            for (auto* edge : edges) {
                // need to set remove the end from the edge and set it as start
                edge->ReplaceNode(end, start);
            }
            vedge->Invalidate();
            continue;
        }

        // use the Cohen-Sutherland method to cut the line to the bounding region
        int32 start_code = ComputeRegionCode(start_pos);
        int32 end_code = ComputeRegionCode(end_pos);
        bool reject = false;

        while (0 != (start_code | end_code) && !reject) {
            // this line is completely outside the region.
            if (0 != (start_code & end_code)) {
                reject = true;
            } else {
                Vector2f* pos0;
                Vector2f* pos1;
                Vector2f dir;
                int32 out_code;
                if (end_code) {
                    out_code = end_code;
                    pos0 = &start_pos;
                    pos1 = &end_pos;
                } else {
                    out_code = start_code;
                    pos0 = &end_pos;
                    pos1 = &start_pos;
                }
                dir.Sub(*pos1, *pos0);
                if (out_code & (LEFT | RIGHT)) {
                    if (out_code & LEFT) {
                        pos1->x = m_min_bounds.x;
                    } else {
                        pos1->x = m_max_bounds.x;
                    }
                    pos1->y = dir.y * (pos1->x - pos0->x) / dir.x + pos0->y;
                } else if (out_code & (BELOW | ABOVE)) {
                    if (out_code & BELOW) {
                        pos1->y = m_min_bounds.y;
                    } else {
                        pos1->y = m_max_bounds.y;
                    }
                    pos1->x = dir.x * (pos1->y - pos0->y) / dir.y + pos0->x;
                }
            }
            start_code = ComputeRegionCode(start_pos);
            end_code = ComputeRegionCode(end_pos);
        }

        if (reject) {
            vedge->GetStart()->RemoveEdge(vedge);
            vedge->GetEnd()->RemoveEdge(vedge);
            vedge->Invalidate();
        } else {
            // if the start position or end positions moved, change them in the
            // edge.
            if (start_pos != vedge->GetStart()->GetPosition()) {
                Node* new_node = new Node(start_pos);
                vedge->SetStart(new_node);
                m_nodes.push_back(new_node);
            }
            if (end_pos != vedge->GetEnd()->GetPosition()) {
                Node* new_node = new Node(end_pos);
                vedge->SetEnd(new_node);
                m_nodes.push_back(new_node);
            }
        }
    }
}

void Voronoi::HandleSiteEvent(const Vector2d& site) {
    if (nullptr == m_root) {
        VoronoiCell* cell = new VoronoiCell();
        m_cells.push_back(cell);
        m_root = new Arc(site, cell);
        return;
    }

    // find the arc that this is over.
    Arc* arc = m_root;
    while (!arc->IsLeaf()) {
        double intersection_x = arc->FindBreakpointX(site.y);
        if (site.x <= intersection_x) {
            arc = arc->GetLeft();
        } else {
            arc = arc->GetRight();
        }
    }
    assert(nullptr != arc);

    // If the y's match, then it has to be on the same y as the root. All arcs
    // at this point will be degenerate. If the root and incoming site's y
    // don't match, at least the first arc will not be degenerate.
    if (arc->GetSite().y == site.y) {
        // create the edge for the root
        Vector2f edge_start(static_cast<float32>(0.5 * (site.x + arc->GetSite().x)),
                            static_cast<float32>(m_max_bounds.y));
        Node* new_node = new Node(edge_start);
        m_nodes.push_back(new_node);
        // create a new cell for the new site.
        VoronoiCell* cell = new VoronoiCell();
        m_cells.push_back(cell);
        VoronoiEdge* edge1 = new VoronoiEdge(new_node, Vector2f(0.0f, -1.0f), arc->GetCell(), cell);
        m_edges.push_back(edge1);
        arc->SetEdge(edge1);
        cell->AddEdge(edge1);
        arc->GetCell()->AddEdge(edge1);

        // set up the new nodes to be added
        Arc* new_left = new Arc(arc->GetSite(), arc->GetCell());
        Arc* new_right = new Arc(site, cell);

        if (arc->GetPrevious()) {
            arc->GetPrevious()->SetNext(new_left);
        }
        new_left->SetPrev(arc->GetPrevious());
        new_left->SetNext(new_right);
        new_right->SetPrev(new_left);
        new_right->SetNext(arc->GetNext());
        if (arc->GetNext()) {
            arc->GetNext()->SetPrev(new_right);
        }

        arc->SetLeft(new_left);
        arc->SetRight(new_right);
        return;
    }

    // if the arc has an event, then the event is not a real event and should
    // be ignored.
    Event* arc_event = arc->GetEvent();
    if (nullptr != arc_event) {
        arc_event->m_valid = false;
        arc->SetEvent(nullptr);
    }

    VoronoiCell* cell = new VoronoiCell();
    m_cells.push_back(cell);

    // create the edge which will be traced out by the new breakpoint.
    Vector2f edge_start;
    edge_start.Set(static_cast<float32>(site.x),
                   static_cast<float32>(arc->GetYAtX(site.x, site.y)));
    Node* new_node = new Node(edge_start);
    m_nodes.push_back(new_node);
    Vector2f diff(static_cast<float32>(arc->GetSite().x - site.x),
                  static_cast<float32>(arc->GetSite().y - site.y));
    Vector2f dir;
    dir.PerpendicularOf(diff);
    VoronoiEdge* left_edge = new VoronoiEdge(new_node, dir, cell, arc->GetCell());
    dir *= -1.0f;
    VoronoiEdge* right_edge = new VoronoiEdge(new_node, dir, cell, arc->GetCell());
    m_edges.push_back(left_edge);
    m_edges.push_back(right_edge);

    // Create the cell which will be created by the site.
    arc->GetCell()->AddEdge(left_edge);
    arc->GetCell()->AddEdge(right_edge);
    cell->AddEdge(left_edge);
    cell->AddEdge(right_edge);

    // set up the new nodes to be added
    Arc* new_left_left = new Arc(arc->GetSite(), arc->GetCell());
    Arc* new_left_right = new Arc(site, cell);
    Arc* new_right = new Arc(arc->GetSite(), arc->GetCell());
    Arc* new_left = new Arc(site, left_edge, new_left_left, new_left_right);

    arc->SetEdge(right_edge);
    arc->SetLeft(new_left);
    arc->SetRight(new_right);

    // Set up the arc neighbor pointers
    if (arc->GetPrevious()) {
        arc->GetPrevious()->SetNext(new_left_left);
    }
    new_left_left->SetPrev(arc->GetPrevious());
    new_left_left->SetNext(new_left_right);
    new_left_right->SetPrev(new_left_left);
    new_left_right->SetNext(new_right);
    new_right->SetPrev(new_left_right);
    new_right->SetNext(arc->GetNext());
    if (arc->GetNext()) {
        arc->GetNext()->SetPrev(new_right);
    }

    CheckForCircleEvent(new_left_left, site.y);
    CheckForCircleEvent(new_right, site.y);
}

void Voronoi::HandleCircleEvent(Event* event) {
    Arc* center = event->m_arc;

    Arc* prev_breakpoint = center->GetLeftBreakpoint();
    Arc* next_breakpoint = center->GetRightBreakpoint();

    Arc* prev = center->GetPrevious();
    Arc* next = center->GetNext();

    // these states should never occur if everything prior is fine. A circle
    // event can only occur when there are at least three leaves in the tree
    // and can never occur on the left- or right-most leaves.
    assert(nullptr != next);
    assert(nullptr != prev);
    assert(next != prev);

    // delete all circle events which contain the arc from Q
    Event* arc_event = prev->GetEvent();
    if (nullptr != arc_event) {
        arc_event->m_valid = false;
        prev->SetEvent(nullptr);
    }
    arc_event = next->GetEvent();
    if (nullptr != arc_event) {
        arc_event->m_valid = false;
        next->SetEvent(nullptr);
    }

    // add the center of the circle causing the event as a vertex record into the graph
    Node* new_node = new Node(Vector2f(static_cast<float32>(event->m_site.x),
                                       static_cast<float32>(event->m_site.y)));
    m_nodes.push_back(new_node);
    next_breakpoint->GetEdge()->SetEnd(new_node);
    prev_breakpoint->GetEdge()->SetEnd(new_node);
    
    /// attach the three new records to the half edge records that end at the vertex.
    Arc* top_breakpoint;
    Arc* cur_arc = center;
    while (cur_arc) {
        if (cur_arc == prev_breakpoint) {
            top_breakpoint = next_breakpoint;
            break;
        } else if (cur_arc == next_breakpoint) {
            top_breakpoint = prev_breakpoint;
            break;
        }
        cur_arc = cur_arc->GetParent();
    }
    // create two half-edge records that correspond to the new breakpoint
    Vector2f diff(static_cast<float32>(prev->GetSite().x - next->GetSite().x),
                  static_cast<float32>(prev->GetSite().y - next->GetSite().y));
    Vector2f dir;
    dir.PerpendicularOf(diff);
    VoronoiEdge* new_edge = new VoronoiEdge(new_node, dir, next->GetCell(), prev->GetCell());
    m_edges.push_back(new_edge);
    top_breakpoint->SetEdge(new_edge);
    prev->GetCell()->AddEdge(new_edge);
    next->GetCell()->AddEdge(new_edge);

    // Remove the event's arc from the tree. After this call, the center is no
    // longer valid
    center->RemoveFromTree();
    center = nullptr;

    // check for circle events at the previous and next nodes.
    CheckForCircleEvent(prev, event->m_y);
    CheckForCircleEvent(next, event->m_y);
}

Arc* Voronoi::FindArcUnderSite(const Vector2d& site) {
    Arc* arc = m_root;
    while (!arc->IsLeaf()) {
        double x = arc->FindBreakpointX(site.y);
        if (x > site.x) {
            arc = arc->GetLeft();
        } else {
            arc = arc->GetRight();
        }
    }
    return arc;
}

// Look for a new circle event for the given arc.
bool Voronoi::CheckForCircleEvent(Arc* arc, double y) {
    // Invalidate any old event.
    Event* arc_event = arc->GetEvent();
    if (nullptr != arc_event && arc_event->m_site.y != y) {
        arc_event->m_valid = false;
        arc->SetEvent(nullptr);
    }

    Arc* prev = arc->GetPrevious();
    Arc* next = arc->GetNext();
    // if prev or next is null, this is an end arc and there is no way for this
    // to produc an arc event.
    if (nullptr == prev || nullptr == next) {
        return false;
    }

    double min_y;
    Vector2d circle_site;

    if (!FindCircleSite(prev->GetSite(), arc->GetSite(), next->GetSite(), min_y, circle_site)) {
        return false;
    }

    Event* new_event = new Event(circle_site, arc, min_y);
    arc->SetEvent(new_event);
    m_events.push(new_event);

    return true;
}

// Find the center and bottommost point on the circle through p0, p1, and p2.
bool Voronoi::FindCircleSite(const Vector2d& p0, const Vector2d& p1, const Vector2d& p2, double& min_y, Vector2d& circle_site) {
    Vector2d l10, l20;
    l10.Sub(p1, p0);
    l20.Sub(p2, p0);

    double area = (l10.x * l20.y - l10.y * l20.x);

    // ignore all points which makes a left turn or where the determinate is 0.
    if (area >= 0) {
      return false;
    }

    circle_site.Set((l20.y * l10.LengthSquared() - l10.y * l20.LengthSquared()) * 0.5f / area,
                    (l10.x * l20.LengthSquared() - l20.x * l10.LengthSquared()) * 0.5f / area);
    double radius = circle_site.Length();
    circle_site += p0;

    // circle_site.x minus radius equals min y coordinate.
    min_y = circle_site.y - radius;
    return true;
}

int32 Voronoi::ComputeRegionCode(const Vector2f& pos) const {
    int32 code = 0;
    if (pos.x < m_min_bounds.x) {
        code |= LEFT;
    } else if (pos.x > m_max_bounds.x){
        code |= RIGHT;
    }

    if (pos.y < m_min_bounds.y) {
        code |= BELOW;
    } else if (pos.y > m_max_bounds.y) {
        code |= ABOVE;
    }
    return code;
}

}
}