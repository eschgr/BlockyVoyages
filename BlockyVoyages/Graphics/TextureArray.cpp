#include "stdafx.h"

#include "TextureArray.h"

#include "PNGReader.h"
#include "Shader.h"
#include "Texture.h"

namespace BlockyVoyages {
namespace Graphics {

TextureArray::TextureArray(const std::vector<std::string>& filenames, bool make_mipsmaps)
    : m_has_mipmaps(make_mipsmaps),
      m_width(0),
      m_height(0),
      m_size(filenames.size()),
      m_tex_handle(0)
{
    if (0 == filenames.size()) {
        LOG(Debugability::DBG_LVL_HIGH, "List of filenames cannot be 0 length for a texture array.");
        return;
    }
    if (1 == filenames.size()) {
        LOG(Debugability::DBG_LVL_WARNING,
                "List of filenames only contains a single entry. "
                "Consider using a Texture instead.");
    }
    // Read in the first texture to get information about the textures. After
    // that, just pass everything to the set texture by index and filename
    // method.
    std::vector<ubyte> buffer;
    PNGReader reader(filenames[0]);
    reader.ReadImage(buffer);
    m_width = reader.GetWidth();
    m_height = reader.GetHeight();
    m_comp_count = Texture::getElementCount(reader.GetComponentType());
    SetupTexture();
    LoadTexture(&buffer[0], 0, reader.GetComponentType(), reader.GetDataType());

    for (uint32 tex = 1; tex < filenames.size(); ++tex) {
        SetTexture(filenames[tex], tex);
    }
}

TextureArray::~TextureArray() {
    glDeleteTextures(1, &m_tex_handle);
}

void TextureArray::Bind() {
    glBindTexture(GL_TEXTURE_2D_ARRAY, m_tex_handle);
}

void TextureArray::LoadTexture(const void* data, int32 index, GLenum format, GLenum dataType) {
    if (m_has_mipmaps) {
        std::vector<void*> out_data;
        Texture::makeMipMaps(data, m_width, m_height, format, dataType, out_data);

        int32 curWidth = m_width;
        int32 curHeight = m_height;
        for (uint32 level = 0; level < out_data.size(); ++level) {
            glTexSubImage3D(GL_TEXTURE_2D_ARRAY, level,
                            0, 0, index,
                            curWidth, curHeight, 1,
                            format, dataType,
                            out_data[level]);
            curWidth >>= 1;
            curHeight >>= 1;
            delete[] out_data[level];
        }
    } else {
        glTexSubImage3D(GL_TEXTURE_2D_ARRAY, 0,
                        0, 0, index,
                        m_width, m_height, 1,
                        format, dataType,
                        data);
    }

    int32 error = glGetError();
    if (error) {
        LOG(Debugability::DBG_LVL_HIGH, "GL Error %d generated during image load.", error);
    }
}

void TextureArray::SetupTexture() {
    glGenTextures(1, &m_tex_handle);

    glBindTexture(GL_TEXTURE_2D_ARRAY, m_tex_handle);

    glPixelStorei(GL_UNPACK_ROW_LENGTH, 0);
    glPixelStorei(GL_UNPACK_SKIP_PIXELS, 0);
    glPixelStorei(GL_UNPACK_SKIP_ROWS, 0);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

    glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    GLenum internal_size_comp;
    int32 mip_levels;
    switch (m_comp_count) {
    case 1:
        internal_size_comp = GL_R8;
        break;
    case 2:
        internal_size_comp = GL_RG8;
        break;
    case 3:
        internal_size_comp = GL_RGB8;
        break;
    case 4:
        internal_size_comp = GL_RGBA8;
    }

    if (m_has_mipmaps) {
        mip_levels = Texture::GetMipMapLevelCount(m_width, m_height);
        glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    } else {
        mip_levels = 1;
        glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    }
    glTexStorage3D(GL_TEXTURE_2D_ARRAY, mip_levels, internal_size_comp, m_width, m_height, m_size);
}

bool TextureArray::SetTexture(const std::string& filename, int32 index) {
    if (index >= m_size) {
        LOG(Debugability::DBG_LVL_WARNING, "Invalid index %d (%d >= %d).", index, index, m_size);
        return false;
    }
    PNGReader reader(filename);
    // this function is a modification of the of the file "example.c" from the libpng library
    // load the png image from file
    std::vector<ubyte> img_buffer;
    reader.ReadImage(img_buffer);

    if (reader.GetWidth() != m_width) {
        LOG(Debugability::DBG_LVL_WARNING, "PNG image '%s' width %d does not match texture array width %d.",
                filename.c_str(), reader.GetWidth(), m_width);
        return false;
    }

    if (reader.GetHeight() != m_height) {
        LOG(Debugability::DBG_LVL_WARNING, "PNG image '%s' width %d does not match texture array width %d.",
                filename.c_str(), reader.GetHeight(), m_height);
        return false;
    }

    LoadTexture(&img_buffer[0],
                index,
                reader.GetComponentType(),
                reader.GetDataType());
    return false;
}

}
}