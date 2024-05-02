#if !defined(__GRAPH_H__)
#define __GRAPH_H__

#include <queue>
#include <vector>

#include "../Math/MathMain.h"

#include "Edge.h"
#include "Node.h"
#include "Cell.h"

namespace BlockyVoyages {

namespace Graphics {
class Shader;
}

namespace Graph {

class Graph {
public:
    Graph();
    Graph(const Vector2f& min_bounds, const Vector2f& max_bounds);
    virtual ~Graph() = 0;

    virtual void Draw(Graphics::Shader* shader);

    virtual void CopyFrom(const Graph* other);

    const std::vector<Edge*>& GetEdges() const { return m_edges; }
    const std::vector<Node*>& GetNodes() const { return m_nodes; }
    const std::vector<Cell*>& GetCells() const { return m_cells; }

    virtual Cell* GetCellAt(const Vector2f& pnt) const;

protected:
    std::vector<Edge*> m_edges;
    std::vector<Node*> m_nodes;
    std::vector<Cell*> m_cells;

    virtual Node* CopyNode(Node* other);
    virtual Edge* CopyEdge(Edge* other);
    virtual Cell* CopyCell(Cell* other);

    Vector2f m_min_bounds;
    Vector2f m_max_bounds;
};

}
}

#endif