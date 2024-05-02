#if !defined(__BIOME_MAP_H__)
#define __BIOME_MAP_H__

#include <random>

#include "../../Graph/Graph.h"

#include "../../Geometry/Segment.h"

namespace BlockyVoyages {
namespace Graph {
    class Node;
    class Edge;
    class Cell;
}

namespace World {
namespace Biome {

using ::BlockyVoyages::Graph::Node;
using ::BlockyVoyages::Graph::Edge;
using ::BlockyVoyages::Graph::Cell;

class BiomeFactory;

class BiomeMap : public Graph::Graph {
public:
    BiomeMap(int32 seed, float32 max_height);

    virtual void CopyFrom(const Graph* other);
    virtual Cell* GetCellAt(const Vector2f& pnt) const;

    void AssignWater();
    void AssignHeight();
    void AssignTemperature();
    void AssignMoisture();
    void AssignBiomes(const BiomeFactory& biome_factory);

    void Draw(Graphics::Shader* shader);

protected:
    virtual Node* CopyNode(Node* other);
    virtual Edge* CopyEdge(Edge* other);
    virtual Cell* CopyCell(Cell* other);

    static const int32 kCellSize = 10;

    std::mt19937 m_rng;
    // offset and scaling for figuring out the map features.
    // May be better as locals in the generation methods
    Vector2f m_simplex_offset;
    float32 m_simplex_scale;
    float32 m_max_height;
    float32 m_max_distance;
    bool m_is_south_hemisphere;
    Vector2f m_wind_dir;

    std::vector<Geometry::Segment> m_mountains;

    std::vector<std::vector<Cell*>> m_region_grid;
    int32 m_region_grid_width;
    int32 m_region_grid_height;

    void CalculateDistanceToShore();
    void BuildRegionGrid();
};

}
}
}

#endif