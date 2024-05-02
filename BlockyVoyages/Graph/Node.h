#if !defined(__NODE_H__)
#define __NODE_H__

#include "../Math/MathMain.h"

#include <map>

namespace BlockyVoyages {
namespace Graph {

class Edge;

class Node {
public:
    Node(const Vector2f& pos)
        : m_pos(pos)
    {}

    Node(const Node* other)
        : m_pos(other->m_pos),
          m_edges(other->m_edges)
    {}

    virtual ~Node() {}

    const Vector2f& GetPosition() const { return m_pos; }
    const std::vector<Edge*>& GetEdges() const { return m_edges; }
    int32 Degree() const { return m_edges.size(); }

    void SetPosition(const Vector2f& pos) { m_pos = pos; }

    void AddEdge(Edge* edge) { m_edges.push_back(edge); }
    void RemoveEdge(Edge* edge) { m_edges.erase(std::find(m_edges.begin(), m_edges.end(), edge)); }

    void UpdatePointers(const std::map<Edge*, Edge*>& edge_map) {
        for (uint32 ind = 0; ind < m_edges.size(); ++ind) {
            m_edges[ind] = edge_map.find(m_edges[ind])->second;
        }
    }

protected:
    std::vector<Edge*> m_edges;
    Vector2f m_pos;
};

}
}

#endif