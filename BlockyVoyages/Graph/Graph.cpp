#include "stdafx.h"
#include "Graph.h"

#include "Cell.h"

#include "../Graphics/Shader.h"
#include "../Graphics/VertexArrayObject.h"

#include <algorithm>
#include <map>

#include <assert.h>

namespace BlockyVoyages {
namespace Graph {

Graph::Graph()
{}

Graph::Graph(const Vector2f& min_bounds, const Vector2f& max_bounds)
    : m_min_bounds(min_bounds),
      m_max_bounds(max_bounds)
{}

Graph::~Graph() {
    for (auto* edge : m_edges) {
        delete edge;
    }

    for (auto* node : m_nodes) {
        delete node;
    }

    for (auto* cell : m_cells) {
        delete cell;
    }
}

void Graph::Draw(Graphics::Shader* shader) {
    std::vector<Graphics::BufferComponent> comps(1);
    comps[0].attrib_name = "position";
    comps[0].element_count = 2;

    int32 color = 0;
    /*for (auto* cell : m_cells)
    {
        //Cell* cell = m_cells[1];//m_cells.size() >> 1];
        const std::vector<Edge*>& edges = cell->GetEdges();
        assert(!edges.empty());

        Graphics::VertexArrayObject cell_obj(shader);
        std::vector<Vector2f> cell_verts;

        cell_verts.reserve(edges.size());
        for (Edge* edge : edges) {
            cell_verts.push_back(edge->GetStart(cell)->GetPosition() * 0.9f + (m_max_coord - m_min_coord) * 0.05f);
        }
        int32 ind = cell_obj.AttachAttributeArray(&cell_verts[0], sizeof(cell_verts[0]), cell_verts.size(), false, comps, false);
        shader->setUniform("color", Vector3f(1.0, 0.0f, static_cast<float32>(color++) / m_cells.size()));
        cell_obj.Draw(GL_TRIANGLE_FAN, cell_verts.size());

        shader->setUniform("color", Vector3f(1.0, 1.0f, 0.0f));
        cell_obj.Draw(GL_POINTS, 1);
    }
    */

    Graphics::VertexArrayObject edge_obj(shader);
    std::vector<Vector2f> edge_verts;

    int32 ind = edge_obj.AttachAttributeArray(nullptr, sizeof(Vector2f), 0, false, comps, false);

    edge_verts.clear();
    edge_verts.reserve(m_edges.size() * 2);
    for(auto edge : m_edges) {
        if (edge->GetStart() && edge->GetEnd()) {
            const Vector2f& start_pos = edge->GetStart()->GetPosition();
            const Vector2f& end_pos = edge->GetEnd()->GetPosition();
            edge_verts.push_back(start_pos);// * 0.9f + (m_max_coord - m_min_coord) * 0.05f);
            edge_verts.push_back(end_pos);// * 0.9f + (m_max_coord - m_min_coord) * 0.05f);
        }
    }

    if (!edge_verts.empty()) {
        shader->setUniform("color", Vector3f(0.0, 0.0f, 1.0f));
        edge_obj.ResetBufferData(ind, &edge_verts[0], edge_verts.size());
        edge_obj.Draw(GL_LINES, edge_verts.size());
    }

    /*
    edge_verts.clear();
    edge_verts.reserve(m_cells.size());
    for(auto cell : m_cells) {
        edge_verts.push_back(cell->GetCenter());// * 0.9f + (m_max_coord - m_min_coord) * 0.05f);
    }

    if (!edge_verts.empty()) {
        shader->setUniform("color", Vector3f(0.0, 0.0f, 1.0f));
        edge_obj.ResetBufferData(ind, &edge_verts[0], edge_verts.size());
        edge_obj.Draw(GL_POINTS, edge_verts.size());
    }*/
}

void Graph::CopyFrom(const Graph* other) {
    m_min_bounds = other->m_min_bounds;
    m_max_bounds = other->m_max_bounds;
    // clear out the old data
    for (auto* node : m_nodes) {
        delete node;
    }
    for (auto* edge : m_edges) {
        delete edge;
    }
    for (auto* cell : m_cells) {
        delete cell;
    }
    m_nodes.clear();
    m_edges.clear();
    m_cells.clear();

    m_nodes.reserve(other->m_nodes.size());
    m_edges.reserve(other->m_edges.size());
    m_cells.reserve(other->m_cells.size());

    std::map<Node*, Node*> node_map;
    std::map<Edge*, Edge*> edge_map;
    std::map<Cell*, Cell*> cell_map;

    // copy all the nodes without mapping the data
    for (auto* node : other->m_nodes) {
        Node* new_node = CopyNode(node);
        m_nodes.push_back(new_node);
        node_map[node] = new_node;
    }

    // copy all the edges without mapping the data
    for (auto* edge : other->m_edges) {
        Edge* new_edge = CopyEdge(edge);
        m_edges.push_back(new_edge);
        edge_map[edge] = new_edge;
    }

    // copy all the cells without mapping the data
    for (auto* cell : other->m_cells) {
        Cell* new_cell = CopyCell(cell);
        m_cells.push_back(new_cell);
        cell_map[cell] = new_cell;
    }

    // map the information in the nodes
    for (auto* node : m_nodes) {
        node->UpdatePointers(edge_map);
    }

    for (auto* edge : m_edges) {
        edge->UpdatePointers(node_map, cell_map);
    }

    for (auto* cell : m_cells) {
        cell->UpdatePointers(edge_map);
    }

    // map the information in the edges
    // map the information in the cells
}

Cell* Graph::GetCellAt(const Vector2f& pnt) const {
    // if the point is outside the bounds, then return nullptr
    if (pnt.x < m_min_bounds.x || pnt.x > m_max_bounds.x) {
        return nullptr;
    }
    if (pnt.y < m_min_bounds.y || pnt.y > m_max_bounds.y) {
        return nullptr;
    }
    // TODO(gesch): improve the efficiency of this by using a quad tree or
    // similar structure. Bug #1 filed for this.
    for (auto* cell : m_cells) {
        if (cell->PointInside(pnt)) {
            return cell;
        }
    }
    // We should never be able to reach here since cells should cover the whole
    // region.
    assert(0);
    return nullptr;
}

Node* Graph::CopyNode(Node* other) {
    return new Node(other);
}

Edge* Graph::CopyEdge(Edge* other) {
    return new Edge(other);
}

Cell* Graph::CopyCell(Cell* other) {
    return new Cell(other);
}

}
}