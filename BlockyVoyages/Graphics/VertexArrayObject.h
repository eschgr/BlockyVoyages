#pragma once

#include <string>
#include <vector>

#include "../types.h"

namespace BlockyVoyages {
namespace Graphics {

class Shader;
class VertexBuffer;

struct BufferComponent {
    std::string attrib_name;
    GLenum data_type;
    int32 element_count;
    bool normalized;
    void* offset;

    BufferComponent()
        : data_type(GL_FLOAT),
          element_count(3),
          offset(nullptr),
          normalized(false)
    {}
};

class VertexArrayObject {
public:
    VertexArrayObject(Shader* shader);
    ~VertexArrayObject();

    int32 AttachAttributeArray(void* buffer,
                               int32 stride,
                               int32 count,
                               bool dynamic,
                               const std::vector<BufferComponent>& components,
                               bool instanced);
    int32 AttachIndexArray(std::vector<uint32> inds);

    void ResetBufferData(int32 vbo,
                         void* buffer,
                         int32 count);
    
    void* GetBufferPointer(int32 vbo);
    void FinishBufferWrite(int32 vbo);
    
    void Bind();
    void Unbind();

    void Draw(GLenum primType, int32 primCount);
    void DrawInstanced(GLenum primType, int32 primCount, int32 instances);

private:
    uint32 m_vao_handle;
    int32 m_indices_ind;
    Shader* m_shader;
    std::vector<VertexBuffer*> m_vbos;
};

}
}