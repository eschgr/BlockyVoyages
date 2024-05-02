#include "stdafx.h"

#include "gl\glExt.h"

#include "OpenGL.h"
#include "VertexBuffer.h"

namespace BlockyVoyages {
namespace Graphics {

VertexBuffer::VertexBuffer(GLenum type,
                           void* buffer,
                           int32 count,
                           int32 stride,
                           bool dynamic)
    : m_type(type),
      m_stride(stride),
      m_dynamic(dynamic),
      m_instanced(false),
      m_buffer_ptr(nullptr)
{
    glGenBuffers(1, &m_handle);
    ResetBuffer(buffer, count);
}

VertexBuffer::~VertexBuffer() {
    glDeleteBuffers(1, &m_handle);
}

void VertexBuffer::ResetBuffer(void* buffer, int32 count) {
    glBindBuffer(m_type, m_handle);
    glBufferData(m_type, m_stride * count, buffer, m_dynamic ? GL_DYNAMIC_DRAW : GL_STATIC_DRAW);
}

void* VertexBuffer::GetBufferPointer() {
    if (nullptr == m_buffer_ptr) {
        glBindBuffer(m_type, m_handle);
        m_buffer_ptr = glMapBuffer(GL_ARRAY_BUFFER, GL_READ_ONLY);
    }
    return m_buffer_ptr;
}

void VertexBuffer::FinishBufferWrite() {
    if (nullptr != m_buffer_ptr) {
        glUnmapBuffer(GL_ARRAY_BUFFER);
        m_buffer_ptr = nullptr;
    }
}

void VertexBuffer::SetInstanced() {
    for (auto loc : m_attrib_locs) {
        glVertexAttribDivisor(loc, 1);
    }
    m_instanced = true;
}

}
}