#pragma once

#include "WorldGenerator.h"
#include "../Graphics/VertexArrayObject.h"
#include "../Geometry/AABB.h"
#include "../Math/MathMain.h"

namespace BlockyVoyages {
namespace Graphics {
class Texture;
class TextureArray;
}

namespace World {

class Camera;

class MapNode {
public:
    int8 m_level;
    // the center of the node 
    Vector3f m_center;
    BlockType m_type;
    Geometry::AABB m_aabb;
    float32 m_width;

    MapNode* m_childNodes;

    MapNode()
        : m_level(0),
          m_center(0.0f, 0.0f, 0.0f),
          m_type(kAir),
          m_width(0.0f),
          m_childNodes(nullptr)
    {
        memset(&m_aabb, 0, sizeof(m_aabb));
    }

    MapNode(int8 level, const Vector3f& center, float32 featureSize)
        : m_level(level),
          m_center(center),
          m_type(kAir),
          m_childNodes(nullptr)
    {
        m_width = featureSize * static_cast<float32>(pow(2.0f, m_level));
        CalculateAABB();
    }

    ~MapNode() {
        delete[] m_childNodes;
    }

    void CalculateAABB() { m_aabb.fromCenterSize(m_center, m_width); }

    float32 GetChildWidth() { return m_width * 0.5f; }

    float32 GetWidth() { return m_width; }

    // calculates the type of the node based upon the child types.
    void calculateType() {
        // find the type of the node. Since this is not the leaf node, the type is
        // determined by the leaf nodes.
        int32 typeScores[kBlockTypeCount];
        int32 maxScore = 0;
        m_type = kAir;
        memset(typeScores, 0, sizeof(typeScores));
        if (m_childNodes) {
            for(int32 i = 0; i < 8; ++i) {
                BlockType nodeType = m_childNodes[i].m_type;
                if(nodeType != kAir && ++typeScores[nodeType] > maxScore) {
                    maxScore = typeScores[nodeType];
                    m_type = nodeType;
                }
            }
        }
    }
};

struct RenderInfo {
    Vector3f offset;
    float32 scale;
    int32 type;
};

/**
 * The map is the main engine for the game. It holds the terrain information
 * for the world.
 */
class Map {
public:
    Map(float32 featureSize, float32 maxHeight, float32 area, int32 seed);
    ~Map();

    // Adds a region to the map which covers the entire specified region,
    // expanded if the region is not aligned to feature size. 2D coordinates
    // are map space coordinates.
    void addRegion(const Vector2f& minCoords, const Vector2f& maxCoords);

    float32 getFeatureSize() { return m_featureSize; }
    Vector2f GetWorldDimensions() { return m_worldGen.GetWorldDimensions(); }

    void Draw(const Camera* camera);
    void DrawCells(const Camera* camera);

    void AddBlock(const Vector3f& center, BlockType type);

    void getIntersectingLeaves(const Geometry::AABB& cube, std::vector<MapNode*>& outNodes) const;

private:
    // The maximum height of the tallest mountain in the world. In meters.
    float32 m_maxHeight;
    // The size of the terrain features (block size). In meters.
    float32 m_featureSize;
    float32 m_halfFeatureSize;
    // The total area of the map
    float32 m_area;
    
    MapNode* m_root;
    BlockyVoyages::Graphics::VertexArrayObject* m_cubeBuffer;
    BlockyVoyages::Graphics::TextureArray* m_textureAtlas;
    BlockyVoyages::Graphics::Shader* m_mapShader;
    std::unique_ptr<BlockyVoyages::Graphics::Shader> m_outlineShader;

    int32 m_block_info_handle;
    std::vector<RenderInfo> m_render_info;

    WorldGenerator m_worldGen;

    void FreeOctree(MapNode* node);
    // expands the octree so that the given point is inside the map
    void CreateChildNodes(MapNode* node);
    void expandOctree(const Vector3f& point);
    void addBlock(const Vector3f& center, BlockType type, MapNode* node);
    void addBlockNoLOD(const Vector3f& center, BlockType type);
    void calculateRegionLOD(const Vector2f& minCoords, const Vector2f& maxCoords, MapNode* node);
    void GetNodesToDraw(MapNode* node,
                        const std::vector<Vector4f>& frustumPlanes,
                        const Vector3f& cameraPos,
                        int8* traverseOrder);
    void DrawOctreeNodeOutline(MapNode* node);
    void getIntersectingLeaves(const Geometry::AABB& cube, MapNode* node, std::vector<MapNode*>& outNodes) const;
};

}
}