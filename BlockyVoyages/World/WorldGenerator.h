#pragma once

#include "BlockTypes.h"
#include "../Math/MathMain.h"

#include <memory>

namespace BlockyVoyages {
namespace Graph {
class Graph;
}

namespace World {
namespace Biome {
class BiomeMap;
class BiomeRegion;
}

class WorldGenerator {
public:
    WorldGenerator(float32 feature_size);
    ~WorldGenerator();

    int32 GetSeed() { return m_seed; }
    Vector2f GetWorldDimensions();

    void GenerateBiomeMap(int32 seed, float32 max_height, float32 area);

    Biome::BiomeMap* GetBiomeMap() { return m_biome_map.get(); }

    // This starts the process of generating a large number of base blocks. To
    // generate a region of blocks, start by calling InitializeRegion. Next, 
    // call GetNextBlock in a loop until false is returned. Once done, call
    // ReleaseRegion to clean up the working space.
    bool InitializeRegion(const Vector2f& minCoord, const Vector2f& maxCoord);
    bool GetNextBlock(Vector3f& center, BlockType& type);
    void ReleaseRegion();

    // GetBlockType is meant to be used to determine the type of a leaf node.
    // It is not meant to be used to generate a large area of terrain.
    BlockType GetBlockType(const Vector3f& blockCenter);
    // GetBlockLODType returns the type of a large block. It can find the type
    // of a block that covers a large area. Not meant for leaf nodes or near
    // leaf nodes. For near leaf nodes, it is more efficient and just as good
    // to take the type that occurs the most of the child blocks.
    BlockType GetLODBlockType(const Vector3f& blockCenter, float32 sideSize);

private:
#if !defined(_DEBUG)
    static const int32 kNumRegions = 25000;
#else
    static const int32 kNumRegions = 50;
#endif
    int32 m_seed;

    float32 m_feature_size;
    float32 m_half_feature_size;

    float32 m_max_height;
    float32 m_area;

    Vector3f m_min_bounds;
    Vector3f m_max_bounds;
    Vector2i m_int_dims;
    Vector2i m_buffer_dims;
    Vector2i m_cur_pos;
    float32 m_cur_height;
    Biome::BiomeRegion* m_cur_region;

    std::vector<float32> m_height_cache;

    std::unique_ptr<Biome::BiomeMap> m_biome_map;

    // Generates the bounds for the world. Uses the seed from the constructor
    // to find a random apsect ratio and then generates the width and height.
    void GenerateBounds();
    // Generates the base graph to be used for the biome map. This currently
    // uses a Voronoi diagram, but could possibly return other types as well.
    std::unique_ptr<Graph::Graph> GenerateBaseGraph();

    float32 GetSurfaceHeight(const Vector2f& coord);
    float32 GetCachedHeight(int32 x, int32 y);
    void CulculateColumnStats();
};

}
}