#include "stdafx.h"

#include "Map.h"

#include <assert.h>

#include <algorithm>
#include <stack>
#include <vector>

#include "../Debugability/Logger.h"
#include "../Graphics/OpenGL.h"
#include "../Graphics/Shader.h"
#include "../Graphics/TextureArray.h"
#include "../Graphics/VertexBuffer.h"
#include "../Math/MathMain.h"
#include "../Math/SimplexNoise.h"
#include "Camera.h"
#include "WorldGenerator.h"


namespace BlockyVoyages {
namespace World {

using Geometry::AABB;

Graphics::Vertex cube_vertices[] = {
    // front
    {Vector3f(0.5f, 0.5f, 0.5f), Vector3f(0.8f, 0.8f, 0.8f),
     Vector2f(1.0f, 1.0f)},  // rtn
    {Vector3f(-0.5f, 0.5f, 0.5f), Vector3f(0.8f, 0.8f, 0.8f),
     Vector2f(0.0f, 1.0f)},  // ltn
    {Vector3f(-0.5f, -0.5f, 0.5f), Vector3f(0.8f, 0.8f, 0.8f),
     Vector2f(0.0f, 0.0f)},  // lbn
    {Vector3f(0.5f, -0.5f, 0.5f), Vector3f(0.8f, 0.8f, 0.8f),
     Vector2f(1.0f, 0.0f)},  // rbn
    // back
    {Vector3f(-0.5f, 0.5f, -0.5f), Vector3f(0.8f, 0.8f, 0.8f),
     Vector2f(1.0f, 1.0f)},  // ltf
    {Vector3f(0.5f, 0.5f, -0.5f), Vector3f(0.8f, 0.8f, 0.8f),
     Vector2f(0.0f, 1.0f)},  // rtf
    {Vector3f(0.5f, -0.5f, -0.5f), Vector3f(0.8f, 0.8f, 0.8f),
     Vector2f(0.0f, 0.0f)},  // rbf
    {Vector3f(-0.5f, -0.5f, -0.5f), Vector3f(0.8f, 0.8f, 0.8f),
     Vector2f(1.0f, 0.0f)},  // lbf
    // right
    {Vector3f(0.5f, 0.5f, 0.5f), Vector3f(0.8f, 0.8f, 0.8f),
     Vector2f(1.0f, 1.0f)},  // rtn
    {Vector3f(0.5f, -0.5f, 0.5f), Vector3f(0.8f, 0.8f, 0.8f),
     Vector2f(1.0f, 0.0f)},  // rbn
    {Vector3f(0.5f, -0.5f, -0.5f), Vector3f(0.8f, 0.8f, 0.8f),
     Vector2f(0.0f, 0.0f)},  // rbf
    {Vector3f(0.5f, 0.5f, -0.5f), Vector3f(0.8f, 0.8f, 0.8f),
     Vector2f(0.0f, 1.0f)},  // rtf
    // left
    {Vector3f(-0.5f, 0.5f, -0.5f), Vector3f(0.8f, 0.8f, 0.8f),
     Vector2f(1.0f, 1.0f)},  // ltf
    {Vector3f(-0.5f, -0.5f, -0.5f), Vector3f(0.8f, 0.8f, 0.8f),
     Vector2f(1.0f, 0.0f)},  // lbf
    {Vector3f(-0.5f, -0.5f, 0.5f), Vector3f(0.8f, 0.8f, 0.8f),
     Vector2f(0.0f, 0.0f)},  // lbn
    {Vector3f(-0.5f, 0.5f, 0.5f), Vector3f(0.8f, 0.8f, 0.8f),
     Vector2f(0.0f, 1.0f)},  // ltn
    // top
    {Vector3f(-0.5f, 0.5f, 0.5f), Vector3f(1.0f, 1.0f, 1.0f),
     Vector2f(1.0f, 1.0f)},  // ltn
    {Vector3f(0.5f, 0.5f, 0.5f), Vector3f(1.0f, 1.0f, 1.0f),
     Vector2f(0.0f, 1.0f)},  // rtn
    {Vector3f(0.5f, 0.5f, -0.5f), Vector3f(1.0f, 1.0f, 1.0f),
     Vector2f(0.0f, 0.0f)},  // rtf
    {Vector3f(-0.5f, 0.5f, -0.5f), Vector3f(1.0f, 1.0f, 1.0f),
     Vector2f(1.0f, 0.0f)},  // ltf
    // bottom
    {Vector3f(0.5f, -0.5f, 0.5f), Vector3f(0.5f, 0.5f, 0.5f),
     Vector2f(1.0f, 1.0f)},  // rbn
    {Vector3f(-0.5f, -0.5f, 0.5f), Vector3f(0.5f, 0.5f, 0.5f),
     Vector2f(0.0f, 1.0f)},  // lbn
    {Vector3f(-0.5f, -0.5f, -0.5f), Vector3f(0.5f, 0.5f, 0.5f),
     Vector2f(0.0f, 0.0f)},  // lbf
    {Vector3f(0.5f, -0.5f, -0.5f), Vector3f(0.5f, 0.5f, 0.5f),
     Vector2f(1.0f, 0.0f)}  // rbf
};

int32 inds[] = {0,  1,  2,  0,  2,  3,  4,  5,  6,  4,  6,  7,
                8,  9,  10, 8,  10, 11, 12, 13, 14, 12, 14, 15,
                16, 17, 18, 16, 18, 19, 20, 21, 22, 20, 22, 23};

float32 offsets[] = {10.0f, 1.0f, 10.0f, 0.25f, -1.0f, 1.0f, 10.0f, 0.5f,
                     15.0f, 1.0f, 10.0f, 5.0f,  41.0f, 1.0f, 10.0f, 0.5f,
                     14.0f, 1.0f, 10.0f, 0.75f, 11.0f, 1.0f, 10.0f, 0.25f,
                     16.0f, 1.0f, 10.0f, 0.1f,  18.0f, 1.0f, 10.0f, 1.0f,
                     12.0f, 1.0f, 10.0f, 2.5f,  19.0f, 1.0f, 10.0f, 0.5f};

Map::Map(float32 featureSize, float32 maxHeight, float32 area, int32 seed)
    : m_root(nullptr),
      m_cubeBuffer(nullptr),
      m_textureAtlas(nullptr),
      m_featureSize(featureSize),
      m_halfFeatureSize(featureSize * 0.5f),
      m_maxHeight(maxHeight),
      m_worldGen(featureSize) {
  // load the VA for the cube. It's not a VBO since all vertices are unique and
  // indexing would only waste space.
  std::vector<uint32> cube_inds(sizeof(inds) / sizeof(inds[0]));
  for (size_t ind = 0; ind < cube_inds.size(); ++ind) {
    cube_inds[ind] = inds[ind];
  }

  int32 glError;

  std::vector<std::string> block_textures(kBlockTypeCount);
  block_textures[0] = "../Data/dirt.png";
  block_textures[1] = "../Data/grass.png";
  block_textures[2] = "../Data/dead_grass.png";
  block_textures[3] = "../Data/snowy_grass.png";
  block_textures[4] = "../Data/temparate.png";
  block_textures[5] = "../Data/swamp.png";
  block_textures[6] = "../Data/snow.png";
  block_textures[7] = "../Data/ice.png";
  block_textures[8] = "../Data/sand.png";
  block_textures[9] = "../Data/gravel.png";
  m_textureAtlas =
      new BlockyVoyages::Graphics::TextureArray(block_textures, true);

  m_mapShader = new BlockyVoyages::Graphics::Shader(
      "../Data/Shaders/terrain.vs", "../Data/Shaders/terrain.fs");
  m_cubeBuffer = new BlockyVoyages::Graphics::VertexArrayObject(m_mapShader);
  std::vector<Graphics::BufferComponent> vertex_components(3);
  vertex_components[0].attrib_name = "position";
  vertex_components[0].element_count = 3;
  vertex_components[0].offset =
      Graphics::BUFFER_OFFSET(offsetof(Graphics::Vertex, pos));
  vertex_components[1].attrib_name = "normal";
  vertex_components[1].element_count = 3;
  vertex_components[1].offset =
      Graphics::BUFFER_OFFSET(offsetof(Graphics::Vertex, normal));
  vertex_components[2].attrib_name = "tex_coord";
  vertex_components[2].element_count = 2;
  vertex_components[2].offset =
      Graphics::BUFFER_OFFSET(offsetof(Graphics::Vertex, texCoords));

  m_cubeBuffer->AttachAttributeArray(
      cube_vertices, sizeof(Graphics::Vertex),
      sizeof(cube_vertices) / sizeof(cube_vertices[0]), false,
      vertex_components, false);
  if ((glError = glGetError()) != 0) {
    LOG(Debugability::DBG_LVL_HIGH, "Got OpenGL error %d.", glError);
  }
  m_cubeBuffer->AttachIndexArray(cube_inds);
  std::vector<Graphics::BufferComponent> pos_components(3);
  pos_components[0].attrib_name = "offset";
  pos_components[0].element_count = 3;
  pos_components[0].offset =
      Graphics::BUFFER_OFFSET(offsetof(RenderInfo, RenderInfo::offset));

  pos_components[1].attrib_name = "scale";
  pos_components[1].element_count = 1;
  pos_components[1].offset =
      Graphics::BUFFER_OFFSET(offsetof(RenderInfo, RenderInfo::scale));

  pos_components[2].attrib_name = "type";
  pos_components[2].element_count = 1;
  pos_components[2].data_type = GL_INT;
  pos_components[2].offset =
      Graphics::BUFFER_OFFSET(offsetof(RenderInfo, RenderInfo::type));
  m_block_info_handle = m_cubeBuffer->AttachAttributeArray(
      nullptr, sizeof(RenderInfo), 0, true, pos_components, true);

  // TODO(gesch): Move this to someplace else when menu enabled map creation is
  // completed.
  m_worldGen.GenerateBiomeMap(seed, m_maxHeight, area);
}

Map::~Map() {
  delete m_cubeBuffer;
  delete m_textureAtlas;
  delete m_root;
}

void Map::addRegion(const Vector2f& minCoord, const Vector2f& maxCoord) {
  // calculate the true bounds of the world.
  float32 minX = std::floor(minCoord.x / m_featureSize) * m_featureSize;
  float32 minY = std::floor(minCoord.y / m_featureSize) * m_featureSize;

  float32 maxX = std::ceil(maxCoord.x / m_featureSize) * m_featureSize;
  float32 maxY = std::ceil(maxCoord.y / m_featureSize) * m_featureSize;

  float32 maxHeight = std::ceil(m_maxHeight / m_featureSize) * m_featureSize;

  if (maxX < minX) {
    std::swap(minX, maxX);
  }
  if (maxY < minY) {
    std::swap(minY, maxY);
  }

  if (nullptr == m_root) {
    Vector3f center(maxX + minX, maxHeight, maxY + minY);
    center *= 0.5f;
    float32 maxDim = maxX - minX;
    if (maxY - minY > maxDim) {
      maxDim = maxY - minY;
    }
    if (maxHeight > maxDim) {
      maxDim = maxHeight;
    }
    maxDim /= m_featureSize;
    int32 level = 0;
    while (maxDim > 1) {
      maxDim *= 0.5;
      ++level;
    }
    m_root = new MapNode(level, center, m_featureSize);
  } else {
    // find the center of the furthest away point to the root's center
    Vector3f farPoint(0.0f, maxHeight * 0.5f, 0.0f);
    Vector2f center(maxX + minX, maxY + minY);
    center *= 0.5f;
    Vector2f diff(center.x - m_root->m_center.x, center.y - m_root->m_center.y);
    farPoint.x = diff.x > 0.0f ? maxX : minX;
    farPoint.z = diff.y > 0.0f ? maxY : minY;
    expandOctree(farPoint);
  }

  // ensure the coordinates align to the centers of the new nodes
  m_worldGen.InitializeRegion(minCoord, maxCoord);
  Vector3f blockCenter;
  BlockType blockType;
  while (m_worldGen.GetNextBlock(blockCenter, blockType)) {
    addBlockNoLOD(blockCenter, blockType);
  }

  calculateRegionLOD(minCoord, maxCoord, m_root);
  m_worldGen.ReleaseRegion();
}

struct Ranking {
  int8 origPos;
  float32 ranking;
};
bool compareRankings(Ranking l, Ranking r) { return l.ranking < r.ranking; }

void Map::Draw(const Camera* camera) {
  if (!m_root) {
    return;
  }

  std::vector<Vector4f> frustumPlanes(camera->getViewFrustumPlanes());
  // determine the order the voxels should be traversed to achieve a front
  // to back ordering.
  Vector4f normals[] = {Vector4f(-1.0, -1.0, -1.0), Vector4f(1.0, -1.0, -1.0),
                        Vector4f(-1.0, -1.0, 1.0),  Vector4f(1.0, -1.0, 1.0),
                        Vector4f(-1.0, 1.0, -1.0),  Vector4f(1.0, 1.0, -1.0),
                        Vector4f(-1.0, 1.0, 1.0),   Vector4f(1.0, 1.0, 1.0)};

  std::vector<Ranking> rankings(8);
  for (int i = 0; i < 8; ++i) {
    rankings[i].origPos = i;
    rankings[i].ranking = normals[i].Dot(frustumPlanes[0]);
  }
  std::sort(rankings.begin(), rankings.end(), compareRankings);
  int8 traverseOrder[8];
  for (int i = 0; i < 8; ++i) {
    traverseOrder[i] = rankings[i].origPos;
  }

  m_render_info.clear();
  GetNodesToDraw(m_root, frustumPlanes, camera->getPosition(), traverseOrder);
  // If there is nothing found to draw, just return
  if (m_render_info.size() == 0) {
    return;
  }

  // set up the terrain shader
  Graphics::g_opengl->SetActiveShader(m_mapShader);
  m_mapShader->setUniform("view_proj_mat", camera->getViewProjectionMatrix());
  Matrix44f id_mat;
  id_mat.Identity();

  m_textureAtlas->Bind();
  m_cubeBuffer->Bind();

  m_cubeBuffer->ResetBufferData(m_block_info_handle, &m_render_info[0],
                                m_render_info.size());
  m_cubeBuffer->DrawInstanced(GL_TRIANGLES, sizeof(inds) / sizeof(inds[0]),
                              m_render_info.size());
}

void Map::DrawCells(const Camera* camera) {
  if (nullptr != m_root) {
    if (m_outlineShader == nullptr) {
      m_outlineShader.reset(new BlockyVoyages::Graphics::Shader(
          "../Data/Shaders/cell_octree.vs", "../Data/Shaders/cell_octree.fs"));
    }
    Graphics::g_opengl->SetActiveShader(m_outlineShader.get());
    m_outlineShader->setUniform("view_proj_mat",
                                camera->getViewProjectionMatrix());

    DrawOctreeNodeOutline(m_root);
  }
}

void Map::AddBlock(const Vector3f& center, BlockType type) {
  expandOctree(center);
  assert(nullptr != m_root);

  // go down until we get to the exact node that needs to be added
  addBlock(center, type, m_root);
}

void Map::getIntersectingLeaves(const AABB& aabb,
                                std::vector<MapNode*>& outNodes) const {
  if (nullptr != m_root) {
    getIntersectingLeaves(aabb, m_root, outNodes);
  }
}

void Map::FreeOctree(MapNode* node) {
  if (nullptr != node && node->m_level > 0) {
    // go through each of the subnodes
    for (int i = 0; i < 8; ++i) {
      FreeOctree(&node->m_childNodes[i]);
    }
    delete node;
  }
}

void Map::CreateChildNodes(MapNode* node) {
  node->m_childNodes = new MapNode[8];
  float32 childWidth = node->GetChildWidth();
  float32 childHalfWidth = childWidth * 0.5f;
  for (int32 i = 0; i < 8; ++i) {
    auto& child = node->m_childNodes[i];
    child.m_level = node->m_level - 1;
    child.m_width = childWidth;
    child.m_center.x =
        node->m_center.x + (i % 2 ? childHalfWidth : -childHalfWidth);
    child.m_center.y =
        node->m_center.y + (i >> 2 ? childHalfWidth : -childHalfWidth);
    child.m_center.z =
        node->m_center.z + ((i >> 1) % 2 ? childHalfWidth : -childHalfWidth);
    child.CalculateAABB();
  }
}

void Map::expandOctree(const Vector3f& point) {
  if (nullptr != m_root) {
    // expand the world until the new block fits.
    while (!m_root->m_aabb.isInside(point)) {
      // The point is outside the current world, so we'll need to expand by one
      // in the direction of where the new block needs to go
      int32 offset = 0;
      Vector3f newCenter;
      if (point.x >= m_root->m_center.x) {
        newCenter.x = m_root->m_aabb.getMaximum().x;
        offset += 1;
      } else {
        newCenter.x = m_root->m_aabb.getMinimum().x;
      }

      if (point.y >= m_root->m_aabb.getMinimum().y) {
        newCenter.y = m_root->m_aabb.getMaximum().y;
        offset += 4;
      } else {
        newCenter.y = m_root->m_aabb.getMinimum().y;
      }

      if (point.z >= m_root->m_center.z) {
        newCenter.z = m_root->m_aabb.getMaximum().z;
        offset += 2;
      } else {
        newCenter.z = m_root->m_aabb.getMinimum().z;
      }
      MapNode* newNode =
          new MapNode(m_root->m_level + 1, newCenter, m_featureSize);
      CreateChildNodes(newNode);
      for (int32 i = 0; i < 8; ++i) {
        if (newNode->m_childNodes[i].m_center.y < point.y) {
          newNode->m_childNodes[i].m_type = kDirt;
        }
      }
      newNode->m_childNodes[offset] = *m_root;
      m_root = newNode;
    }
  } else {
    // if the root is null, this is a new map. Use the point as the root.
    m_root = new MapNode(0, point, m_featureSize);
  }
}

void Map::addBlock(const Vector3f& center, BlockType type, MapNode* node) {
  // go down until we get to the exact node that needs to be added
  std::stack<MapNode*> needsUpdate;
  if (node->m_level == 0) {
    node->m_type = type;
    return;
  }

  if (nullptr == node->m_childNodes) {
    CreateChildNodes(node);
    for (int32 i = 0; i < 8; ++i) {
      if (node->m_childNodes[i].m_center.y < center.y) {
        node->m_childNodes[i].m_type = type;
      }
    }
  }

  // use center information to get which sub-node this block will be in
  int offset = 0;
  if (center.x >= node->m_center.x) {
    offset += 1;
  }
  if (center.y >= node->m_center.y) {
    offset += 4;
  }
  if (center.z >= node->m_center.z) {
    offset += 2;
  }

  addBlock(center, type, &node->m_childNodes[offset]);

  node->calculateType();
}

void Map::addBlockNoLOD(const Vector3f& center, BlockType type) {
  // go down until we get to the exact node that needs to be added
  MapNode* node = m_root;

  while (node->m_level > 0) {
    if (nullptr == node->m_childNodes) {
      CreateChildNodes(node);
    }

    // use center information to get which sub-node this block will be in
    int offset = 0;
    if (center.x >= node->m_center.x) {
      offset += 1;
    }
    if (center.y >= node->m_center.y) {
      offset += 4;
    }
    if (center.z >= node->m_center.z) {
      offset += 2;
    }
    node = &node->m_childNodes[offset];
  }
  node->m_type = type;
}

void Map::calculateRegionLOD(const Vector2f& minCoord, const Vector2f& maxCoord,
                             MapNode* node) {
  if (node->m_level > 0) {
    if (node->m_childNodes) {
      for (int32 i = 0; i < 8; ++i) {
        if (maxCoord.x >= node->m_aabb.getMinimum().x &&
            node->m_aabb.getMaximum().x >= minCoord.x &&
            maxCoord.y >= node->m_aabb.getMinimum().z &&
            node->m_aabb.getMaximum().z >= minCoord.y) {
          calculateRegionLOD(minCoord, maxCoord, &node->m_childNodes[i]);
        }
      }
    }

    // calculate the LOD of this node.
    node->calculateType();
  }
}

void Map::GetNodesToDraw(MapNode* node,
                         const std::vector<Vector4f>& frustumPlanes,
                         const Vector3f& cameraPos, int8* traverseOrder) {
  // determine if the node is visible
  // skip the checks on the lower levels because at worst, they'll be just
  // offscreen
  float32 halfWidth = node->GetChildWidth();
  float32 distance = cameraPos.Distance(node->m_center) - halfWidth;
  int32 targetLevel = static_cast<int32>(log((cameraPos.Distance(node->m_center) - halfWidth) / 65.0f + 1) / log(1.5f));

    if (targetLevel >= node->m_level) {
      if (node->m_type != kAir) {
        RenderInfo info;
        info.offset = node->m_center;
        info.scale = node->GetWidth();
        info.type = node->m_type;
        m_render_info.push_back(info);
      }
    } else if (node->m_childNodes) {
      if (node->m_level > 3) {
        for (int32 i = 0; i < 8; ++i) {
          MapNode* nextNode = &node->m_childNodes[traverseOrder[i]];
          bool visible = true;
          for (int i = 0; i < 6 && visible; ++i) {
            visible = visible && nextNode->m_aabb.onFrontHalf(frustumPlanes[i]);
          }
          if (visible) {
            GetNodesToDraw(nextNode, frustumPlanes, cameraPos, traverseOrder);
          }
        }
      } else {
        for (int i = 0; i < 8; ++i) {
          MapNode* nextNode = &node->m_childNodes[traverseOrder[i]];
          GetNodesToDraw(nextNode, frustumPlanes, cameraPos, traverseOrder);
        }
      }
    }
}

void Map::DrawOctreeNodeOutline(MapNode* node) {
    assert(nullptr != node);

    {
      BlockyVoyages::Graphics::VertexArrayObject outline_buffer(m_mapShader);
      std::vector<Graphics::BufferComponent> vertex_components(1);
      vertex_components[0].attrib_name = "position";
      vertex_components[0].element_count = 3;

      const Vector3f& min_bound = node->m_aabb.getMinimum();
      const Vector3f& max_bound = node->m_aabb.getMaximum();

      Vector3f vertices[] = {
          Vector3f(min_bound.x, min_bound.y, min_bound.z),
          Vector3f(max_bound.x, min_bound.y, min_bound.z),
          Vector3f(max_bound.x, min_bound.y, min_bound.z),
          Vector3f(max_bound.x, min_bound.y, max_bound.z),
          Vector3f(max_bound.x, min_bound.y, max_bound.z),
          Vector3f(min_bound.x, min_bound.y, max_bound.z),
          Vector3f(min_bound.x, min_bound.y, max_bound.z),
          Vector3f(min_bound.x, min_bound.y, min_bound.z),
          Vector3f(max_bound.x, min_bound.y, min_bound.z),
          Vector3f(max_bound.x, max_bound.y, min_bound.z),
          Vector3f(max_bound.x, max_bound.y, min_bound.z),
          Vector3f(max_bound.x, max_bound.y, max_bound.z),
          Vector3f(max_bound.x, max_bound.y, max_bound.z),
          Vector3f(max_bound.x, min_bound.y, max_bound.z),
          Vector3f(max_bound.x, max_bound.y, max_bound.z),
          Vector3f(min_bound.x, max_bound.y, max_bound.z),
          Vector3f(min_bound.x, max_bound.y, max_bound.z),
          Vector3f(min_bound.x, max_bound.y, min_bound.z),
          Vector3f(min_bound.x, max_bound.y, min_bound.z),
          Vector3f(max_bound.x, max_bound.y, min_bound.z),
          Vector3f(min_bound.x, max_bound.y, max_bound.z),
          Vector3f(min_bound.x, min_bound.y, max_bound.z),
          Vector3f(min_bound.x, min_bound.y, min_bound.z),
          Vector3f(min_bound.x, max_bound.y, min_bound.z),
      };

      outline_buffer.AttachAttributeArray(
          vertices, sizeof(Vector3f), sizeof(vertices) / sizeof(vertices[0]),
          false, vertex_components, false);
      int32 glError;
      if ((glError = glGetError()) != 0) {
        LOG(Debugability::DBG_LVL_HIGH, "Got OpenGL error %d.", glError);
      }

      outline_buffer.Draw(GL_LINES, sizeof(vertices) / sizeof(vertices[0]));
    }

    /*if (node->m_childNodes && node->m_level > 4) {
        for (int i = 0; i < 8; ++i) {
            DrawOctreeNodeOutline(&node->m_childNodes[i]);
        }
    }*/
}

void Map::getIntersectingLeaves(const AABB& cube, MapNode* node,
                                std::vector<MapNode*>& outNodes) const {
    // if the cube doesn't intersect the node at any point, it never will.
    // Return empty
    if (!cube.intersects(node->m_aabb)) {
      return;
    }

    // We know this intersects, so just return it.
    if (node->m_level == 0) {
      if (node->m_type != kAir) {
        outNodes.push_back(node);
      }
      return;
    }

    if (node->m_childNodes) {
      for (int32 i = 0; i < 8; ++i) {
        getIntersectingLeaves(cube, &node->m_childNodes[i], outNodes);
      }
    }
}

}  // namespace World
}  // namespace BlockyVoyages