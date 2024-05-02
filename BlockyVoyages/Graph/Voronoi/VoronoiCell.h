#if !defined(__VORONOI_CELL_H__)
#define __VORONOI_CELL_H__

#include <vector>

#include "../Cell.h"

#include "../../Math/MathMain.h"

namespace BlockyVoyages {
namespace Graph {

class Edge;
class Node;

class VoronoiCell : public Cell {
public:
    void RemoveBadEdges();
    void Finish(const Vector2f& min_coord, const Vector2f& max_coord,
                std::vector<Edge*>& all_edges, std::vector<Node*>& all_nodes);

private:
    void CheckNodesConnected(Node* start_node, Node* end_node,
                             const Vector2f& min_coord, const Vector2f& max_coord,
                             std::vector<Edge*>& all_edges, std::vector<Node*>& all_nodes);
};

}
}

#endif