#include "stdafx.h"

#include "Font.h"

#include <Windows.h>

#include <gl\gl.h>
#include <gl\glu.h>

#include <cstdlib>
#include <fstream>
#include <vector>

#include "OpenGL.h"
#include "Shader.h"
#include "Texture.h"
#include "VertexArrayObject.h"
#include "../MainWindow.h"
#include "../Math/MathMain.h"

namespace BlockyVoyages {
namespace Graphics {

Font::Font(const FT_Library& ftLib, const std::string& name, int32 pnt_size, Shader* font_shader)
    : m_face(nullptr),
      m_name(name),
      m_pntSize(pnt_size),
      m_texture(nullptr),
      m_nRefs(0),
      m_fontShader(font_shader)
{
    // create the font file name.
    std::string val = getenv("windir");
    std::string fullName = val + "\\Fonts\\" + name;

    // load the face
    int32 err = FT_New_Face(ftLib, fullName.c_str(), 0, &m_face);
    if (err) {
        // find out why there was an error
        if (err == FT_Err_Unknown_File_Format) {
            //LOG(DBG_LVL_CRITICAL, "Font file format unknown.");
            return;
        }
        else {
            std::ifstream strm(fullName.c_str());
            if (!strm.is_open()) {
                //LOG(DBG_LVL_CRITICAL, "Unable to find font file %s.", fullName.c_str());
            }
            else {
                //LOG(DBG_LVL_CRITICAL, "Could not load font for unknown reason. ERROR CODE: %d", err);
                strm.close();
            }
            return;
        }
    }

    // one possible glyph for each character
    m_glyphs.resize(256);
    memset(&m_glyphs[0], 0, sizeof(Glyph) * 256);

    // Use the size of the font to mean in pixels
    FT_Set_Pixel_Sizes(m_face, 0, pnt_size);
    m_pntSize = pnt_size;

    // create the font texture
    CreateFontTexture();
}

Font::~Font() {
    delete m_texture;

    FT_Done_Face(m_face);
}

bool Font::changeFontSize(int32 newpnt_size) {
    // if the new and the old are the same size, then we don't have to do anything
    if (newpnt_size == m_pntSize) {
        return true;
    }

    delete m_texture;
    m_pntSize = newpnt_size;

    // clear out the glyphs
    memset(&m_glyphs[0], 0, sizeof(Glyph) * 256);

    FT_Set_Pixel_Sizes(m_face, 0, m_pntSize);

    // create the font texture
    return CreateFontTexture();
}

struct FontVertex {
    Vector2f pos;
    Vector2f texCoord;
};

void Font::drawText(int32 x, int32 y, const Vector3f& color, const std::string& str) const {
    bool useKerning = FT_HAS_KERNING(m_face) != 0;
    int32 prevInd = 0;
    Vector2i pnt(x, y);

    std::vector<FontVertex> verts(str.size() * 6);
    int32 curVert = 0;
    int32 curInd = 0;
    for(size_t cnt = 0; cnt < str.length(); ++cnt) {
        if (m_glyphs[str[cnt]].glyphInd) {
            if (useKerning && prevInd) {
                FT_Vector delta;
                FT_Get_Kerning(m_face, prevInd, m_glyphs[str[cnt]].glyphInd, 0, &delta);

                pnt.x += delta.x >> 6;
            }

            Vector2i tmp;
            tmp.Add(pnt, m_glyphs[str[cnt]].drawStart);

            float32 minX = static_cast<float32>(tmp.x);
            float32 maxX = static_cast<float32>(tmp.x + m_glyphs[str[cnt]].dims.x);
            float32 minY = static_cast<float32>(tmp.y);
            float32 maxY = static_cast<float32>(tmp.y + m_glyphs[str[cnt]].dims.y);

            verts[curVert].texCoord.Set(m_glyphs[str[cnt]].minCoords.x, 
                                         m_glyphs[str[cnt]].minCoords.y);
            verts[curVert++].pos.Set(minX, minY);

            verts[curVert].texCoord.Set(m_glyphs[str[cnt]].maxCoords.x,
                                         m_glyphs[str[cnt]].minCoords.y);
            verts[curVert++].pos.Set(maxX, minY);

            verts[curVert].texCoord.Set(m_glyphs[str[cnt]].maxCoords.x, 
                                         m_glyphs[str[cnt]].maxCoords.y);
            verts[curVert++].pos.Set(maxX, maxY);

            verts[curVert].texCoord.Set(m_glyphs[str[cnt]].minCoords.x, 
                                         m_glyphs[str[cnt]].minCoords.y);
            verts[curVert++].pos.Set(minX, minY);

            verts[curVert].texCoord.Set(m_glyphs[str[cnt]].maxCoords.x, 
                                         m_glyphs[str[cnt]].maxCoords.y);
            verts[curVert++].pos.Set(maxX, maxY);

            verts[curVert].texCoord.Set(m_glyphs[str[cnt]].minCoords.x,
                                         m_glyphs[str[cnt]].maxCoords.y);
            verts[curVert++].pos.Set(minX, maxY);
        }

        prevInd = m_glyphs[str[cnt]].glyphInd;

        pnt += m_glyphs[str[cnt]].advance;
    }
    VertexArrayObject object(m_fontShader);
    std::vector<BufferComponent> components(2);
    components[0].attrib_name = "position";
    components[0].element_count = 2;
    components[0].offset = BUFFER_OFFSET(offsetof(FontVertex, pos));
    components[1].attrib_name = "tex_coord";
    components[1].element_count = 2;
    components[1].offset = BUFFER_OFFSET(offsetof(FontVertex, texCoord));
    object.AttachAttributeArray(&verts[0], sizeof(FontVertex), verts.size(), false, components, false);

    // draw the text
    g_opengl->SetActiveShader(m_fontShader);
    m_fontShader->setUniform("screen_width", static_cast<float32>(g_mainWindow->GetWidth()));
    m_fontShader->setUniform("screen_height", static_cast<float32>(g_mainWindow->GetHeight()));
    m_fontShader->setUniform("color", color);

    m_texture->bind();
    object.Bind();
    object.Draw(GL_TRIANGLES, verts.size());
    object.Unbind();
}

Vector2i Font::getStringDrawSize(const std::string& str) const {
    bool useKerning = FT_HAS_KERNING(m_face) != 0;
    int32 prevInd = 0;

    Vector2i pnt(0, 0);

    Vector2i minDims(0, 0);
    Vector2i maxDims(0, 0);

    for(size_t cnt = 0; cnt <str.length(); ++cnt) {
        if (m_glyphs[str[cnt]].glyphInd) {
            if (useKerning && prevInd) {
                FT_Vector delta;
                FT_Get_Kerning(m_face, prevInd, m_glyphs[str[cnt]].glyphInd, 0, &delta);

                pnt.x += delta.x >> 6;
            }

            Vector2i tmp;
            tmp.Add(pnt, m_glyphs[str[cnt]].drawStart);
            if (tmp.x < minDims.x) {
                minDims.x = tmp.x;
            }

            if (tmp.y < minDims.y) {
                minDims.y = tmp.y;
            }

            int32 dim = tmp.x + m_glyphs[str[cnt]].dims.x;
            if (dim> maxDims.x) {
                maxDims.x = dim;
            }

            dim = tmp.y + m_glyphs[str[cnt]].dims.y;
            if (dim> maxDims.y) {
                maxDims.y = dim;
            }
        }

        prevInd = m_glyphs[str[cnt]].glyphInd;
        pnt += m_glyphs[str[cnt]].advance;
    }

    return maxDims - minDims;
}

Vector2i Font::getCharDrawSize(ubyte letter) const {
    if (!m_glyphs.empty()) {
        return m_glyphs[letter].dims;
    }

    return Vector2i(0, 0);
}

bool Font::CreateFontTexture(void) {
    int32 nonZero = 0;

    FT_Matrix mat;
    mat.xx = (FT_Fixed)(0x10000L);
    mat.xy = mat.yx = 0;
    mat.yy = (FT_Fixed)(-0x10000L);

    FT_Set_Transform(m_face, &mat, nullptr);
        
    for(int32 letter = 0; letter < 256; ++letter) {
        // store the glyph index for later use (mostly kerning)
        m_glyphs[letter].glyphInd = FT_Get_Char_Index(m_face, letter);

        if (m_glyphs[letter].glyphInd != 0) {
            ++nonZero;
        }
    }

    FT_Glyph glyphs[256];
    int32 letsPerRow = static_cast<int32>(sqrt(static_cast<float32>(nonZero))) + 1;
    int32 textureHeight = 0;
    int32 textureWidth = 0;
    int32 rowWidth = 0;
    int32 rowHeight = 0;
    int32 row = 0;
    int32 col = 0;

    // preload every letter
    for(int32 letter = 0; letter < 256; ++letter) {
        // find the dimensions of each row
        if (FT_Load_Char(m_face, letter, FT_LOAD_RENDER)) {
            // LOG(DBG_LVL_WARNING, "Can't load glyph for character %c.", letter);
            continue;
        }

        FT_GlyphSlot slot = m_face->glyph;

        // store the advance, even if there is no glyph for this one
        m_glyphs[letter].advance.x = slot->advance.x >> 6;
        m_glyphs[letter].advance.y = slot->advance.y >> 6;

        // if the letter has an index available to it
        if (m_glyphs[letter].glyphInd > 0) {
            m_glyphs[letter].drawStart.x = slot->bitmap_left;
            m_glyphs[letter].drawStart.y = -slot->bitmap_top;

            FT_BBox glyphBBox;
            if (FT_Get_Glyph(slot, &glyphs[letter])) {
                continue;
            }

            FT_Glyph_Get_CBox(glyphs[letter], FT_GLYPH_BBOX_TRUNCATE, &glyphBBox);
            m_glyphs[letter].maxDims.Set(glyphBBox.xMax, glyphBBox.yMax);
            m_glyphs[letter].minDims.Set(glyphBBox.xMin, glyphBBox.yMin);

            m_glyphs[letter].dims.Sub(m_glyphs[letter].maxDims, m_glyphs[letter].minDims);

            rowHeight = max(rowHeight, m_glyphs[letter].dims.y);

            rowWidth += m_glyphs[letter].dims.x;

            ++col;
            if (col == letsPerRow || letter == 255) {
                ++row;
                col = 0;
                if (rowWidth> textureWidth) {
                    textureWidth = rowWidth;
                }
                textureHeight += rowHeight;
                rowWidth = 0;
                rowHeight = 0;
            }
        }
    }

    textureWidth += 3 * (letsPerRow - 1) + 2;
    textureHeight += 3 * (letsPerRow - 1) + 2;
    int32 drawX = 1;
    int32 drawY = 1;

    rowHeight = 0;
    col = 0;

    std::vector<uint8> buff(textureWidth * textureHeight, 0);

    for(int32 letter = 0; letter < 256; ++letter) {
        if (m_glyphs[letter].glyphInd) {
            // load the text texture into video memory
            FT_BitmapGlyph slot = (FT_BitmapGlyph)glyphs[letter];
                
            m_glyphs[letter].minCoords.Set(static_cast<float32>(drawX), static_cast<float32>(drawY));
            m_glyphs[letter].maxCoords.x = m_glyphs[letter].minCoords.x + 
                                                static_cast<float32>(m_glyphs[letter].dims.x);
            m_glyphs[letter].maxCoords.y = m_glyphs[letter].minCoords.y + 
                                                static_cast<float32>(m_glyphs[letter].dims.y) ;

            for(int32 pixRow = 0; pixRow < slot->bitmap.rows; ++pixRow) {
                int32 dstPix = drawX + (drawY + pixRow) * textureWidth;
                int32 srcPix = pixRow * slot->bitmap.width;
                memcpy(&buff[dstPix], &slot->bitmap.buffer[srcPix], slot->bitmap.width);
            }

            drawX += (m_glyphs[letter].dims.x) + 3;
            rowHeight = max(rowHeight, m_glyphs[letter].dims.y);
            ++col;

            if (col == letsPerRow) {
                drawX = 1;
                drawY += rowHeight + 3;
                rowHeight = 0;
                col = 0;
            }

            FT_Done_Glyph(glyphs[letter]);
        }
    }

    m_texture = new Texture(&buff[0], textureWidth, textureHeight, GL_RED, GL_UNSIGNED_BYTE, false);
    m_texture->setParameter(GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    m_texture->setParameter(GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    m_texture->setParameter(GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    m_texture->setParameter(GL_TEXTURE_MIN_FILTER, GL_LINEAR);

    for(int32 letter = 0; letter < 256; ++letter) {
        if (m_glyphs[letter].glyphInd != 0) {
            m_glyphs[letter].minCoords.x /= textureWidth;
            m_glyphs[letter].minCoords.y /= textureHeight;
            m_glyphs[letter].maxCoords.x /= textureWidth;
            m_glyphs[letter].maxCoords.y /= textureHeight;
        }
    }
        
    return true;
}

}
}