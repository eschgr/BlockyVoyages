#include "stdafx.h"

#include "Texture.h"

#include "PNGReader.h"
#include "../Debugability/Logger.h"
#include "../Math/MathMain.h"

#include <assert.h>

namespace BlockyVoyages {
namespace Graphics {

Texture::Texture(const std::string& filename, bool make_mips) {
    glGenTextures(1, &m_texHandle);

    PNGReader reader(filename);
    // this function is a modification of the of the file "example.c" from the libpng library
    // load the png image from file
    std::vector<ubyte> img_buffer;
    reader.ReadImage(img_buffer);
    loadTexture(&img_buffer[0],
                reader.GetWidth(),
                reader.GetHeight(),
                reader.GetComponentType(),
                reader.GetDataType(),
                make_mips);
}

Texture::Texture(const void* data, int32 width, int32 height, GLenum format, GLenum dataType, bool makeMips)
    : m_texHandle(0)
{
    glGenTextures(1, &m_texHandle);

    loadTexture(data, width, height, format, dataType, makeMips);
}

Texture::~Texture(void) {
    glDeleteTextures(1, &m_texHandle);
}

void Texture::bind(void) const {
    glBindTexture(GL_TEXTURE_2D, m_texHandle);
}

void Texture::setParameter(GLenum param, uint32 value) {
    glTexParameteri(GL_TEXTURE_2D, param, value);
}

size_t Texture::getDataSize(GLenum type) {
    switch(type) {
        case GL_UNSIGNED_BYTE:
            return 1;
        case GL_INT:
        case GL_FLOAT:
            LOG(Debugability::DBG_LVL_WARNING,
                    "Integer and float types not supported in textures. See DataBuffers instead.");
            return 0;
    }
    return 0;
}

size_t Texture::getElementCount(GLenum format) {
    switch(format) {
    case GL_RED:
        return 1;
    case GL_RG:
        return 2;
    case GL_RGB:
        return 3;
    case GL_RGBA:
        return 4;
    }

    return 0;
}

GLenum Texture::GetSizedFormat(GLenum dataType, GLenum format) {
    switch (format) {
    case GL_RED:
        switch (dataType) {
        case GL_UNSIGNED_BYTE:
            return GL_R8;
        case GL_INT:
            return GL_R32I;
        case GL_FLOAT:
            return GL_R32F;
        }
        break;
    case GL_RG:
        switch (dataType) {
        case GL_UNSIGNED_BYTE:
            return GL_RG8;
        case GL_INT:
            return GL_RG32I;
        case GL_FLOAT:
            return GL_RG32F;
        }
        break;
    case GL_RGB:
        switch (dataType) {
        case GL_UNSIGNED_BYTE:
            return GL_RGB8;
        case GL_INT:
            return GL_RGB32I;
        case GL_FLOAT:
            return GL_RGB32F;
        }
        break;
    case GL_RGBA:
        switch (dataType) {
        case GL_UNSIGNED_BYTE:
            return GL_RGBA8;
        case GL_INT:
            return GL_RGBA32I;
        case GL_FLOAT:
            return GL_RGBA32F;
        }
        break;
    }
    return 0;
}

int32 Texture::GetMipMapLevelCount(int32 width, int32 height) {
    int32 wLog2 = intLog2(width);
    int32 hLog2 = intLog2(height);
    return min(wLog2, hLog2) + 1;
}

void Texture::makeMipMaps(const void* data, int32 width, int32 height,
                          GLenum format, GLenum dataType,
                          std::vector<void*>& out_data) {
    // the following later will need to allow for FBOs

    // if the picture does not need to be adjusted, do something different to create the mipmaps
    // since we are creating mip maps, we hope that we don't have to do anything special
    // since we don't need to adjust for anything, just use the very simple algorithm
    int32 levels = GetMipMapLevelCount(width, height);

    int32 nComps = getElementCount(format);
    int32 dataSize = getDataSize(dataType);

    out_data.resize(levels);
    out_data[0] = new ubyte[dataSize * nComps * width * height];
    memcpy(out_data[0], data, dataSize * nComps * width * height);

    for(int32 level = 1; level < levels; ++level) {
        int32 curWidth = width >> level;
        int32 curHeight = height >> level;

        size_t outPix = 0;
        size_t inPix = 0;
        size_t inElemSpan = nComps * (width >> (level - 1));

        out_data[level] = new ubyte[dataSize * nComps * curWidth * curHeight];
        
        ubyte* curData = reinterpret_cast<ubyte*>(out_data[level]);
        ubyte* prevData0 = &reinterpret_cast<ubyte*>(out_data[level - 1])[0];
        ubyte* prevData1 = &reinterpret_cast<ubyte*>(out_data[level - 1])[nComps * dataSize];
        ubyte* prevData2 = &reinterpret_cast<ubyte*>(out_data[level - 1])[inElemSpan * dataSize];
        ubyte* prevData3 = &reinterpret_cast<ubyte*>(out_data[level - 1])[(inElemSpan + nComps) * dataSize];

        for(int32 row = 0; row < curHeight; ++row, inPix += inElemSpan) {
            for(int32 col = 0; col < curWidth; ++col, inPix += nComps) {
                for(int32 comp = 0; comp < nComps; ++comp) {
                    switch(dataType) {
                        case GL_UNSIGNED_BYTE:
                            *curData = (*prevData0 >> 2) + (*prevData1 >> 2) +
                                       (*prevData2 >> 2) + (*prevData3 >> 2);
                            break;
                        case GL_INT:
                            *reinterpret_cast<int32*>(curData) =
                                (*reinterpret_cast<int32*>(prevData0) >> 2) + 
                                (*reinterpret_cast<int32*>(prevData1) >> 2) +
                                (*reinterpret_cast<int32*>(prevData2) >> 2) +
                                (*reinterpret_cast<int32*>(prevData3) >> 2);
                            break;
                        case GL_FLOAT:
                            *reinterpret_cast<float32*>(curData) = 
                                (*reinterpret_cast<float32*>(prevData0) + 
                                 *reinterpret_cast<float32*>(prevData1) +
                                 *reinterpret_cast<float32*>(prevData2) +
                                 *reinterpret_cast<float32*>(prevData3)) * 0.25f;
                            break;
                    }
                    curData += dataSize;
                    prevData0 += dataSize;
                    prevData1 += dataSize;
                    prevData2 += dataSize;
                    prevData3 += dataSize;
                }
                // skip a pixel of the input data
                prevData0 += nComps * dataSize;
                prevData1 += nComps * dataSize;
                prevData2 += nComps * dataSize;
                prevData3 += nComps * dataSize;
            }
            // skip a row of input data
            prevData0 += inElemSpan * dataSize;
            prevData1 += inElemSpan * dataSize;
            prevData2 += inElemSpan * dataSize;
            prevData3 += inElemSpan * dataSize;
        }
    }
}

void Texture::loadTexture(const void* data, int32 width, int32 height, GLenum format, GLenum dataType, bool makeMips) {
    assert (width > 0 && height > 0);

    GLenum sized_format = GetSizedFormat(dataType, format);

    glBindTexture(GL_TEXTURE_2D, m_texHandle);

    glPixelStorei(GL_UNPACK_ROW_LENGTH, 0);
    glPixelStorei(GL_UNPACK_SKIP_PIXELS, 0);
    glPixelStorei(GL_UNPACK_SKIP_ROWS, 0);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    if (makeMips) {
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
        std::vector<void*> out_data;
        makeMipMaps(data, width, height, format, dataType, out_data);

        int32 curWidth = width;
        int32 curHeight = height;
        glTexStorage2D(GL_TEXTURE_2D, 1, sized_format, width, height);
        for (uint32 level = 0; level < out_data.size(); ++level) {
            glTexSubImage2D(GL_TEXTURE_2D, level,
                        0, 0,
                        width, height,
                        format, dataType,
                        out_data[level]);
            curWidth >>= 1;
            curHeight >>= 1;
            delete[] out_data[level];
        }
    } else {
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        int32 error = glGetError();
        if (error) {
            LOG(Debugability::DBG_LVL_HIGH, "GL Error %d generated during image load.", error);
        }
        glTexStorage2D(GL_TEXTURE_2D, 1, sized_format, width, height);
        error = glGetError();
        if (error) {
            LOG(Debugability::DBG_LVL_HIGH, "GL Error %d generated during image load.", error);
        }
        glTexSubImage2D(GL_TEXTURE_2D, 0,
                        0, 0,
                        width, height,
                        format, dataType,
                        data);
        error = glGetError();
        if (error) {
            LOG(Debugability::DBG_LVL_HIGH, "GL Error %d generated during image load.", error);
        }
    }

    int32 error = glGetError();
    if (error) {
        LOG(Debugability::DBG_LVL_HIGH, "GL Error %d generated during image load.", error);
    }
}

}
}