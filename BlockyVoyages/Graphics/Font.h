#pragma once

#include <string>
#include <vector>

#include "../../freetype-2.5.0.1/include/ft2build.h"
#include FT_FREETYPE_H
#include FT_GLYPH_H

#include "../Math/MathMain.h"
//#include "Texture.h"
// #include "FontManager.h"

namespace BlockyVoyages {
namespace Graphics {

class Texture;
class Shader;

struct Glyph {
    int32 glyphInd;
    Vector2i maxDims;
    Vector2i minDims;
    Vector2i advance;
    Vector2i drawStart;

    Vector2f minCoords;
    Vector2f maxCoords;

    Vector2i dims;
};

class Font {
public:
    bool changeFontSize(int32 newPntSize);

    void drawText(int32 x, int32 y, const Vector3f& color, const std::string& str) const;
    Vector2i getStringDrawSize(const std::string& str) const;
    Vector2i getCharDrawSize(ubyte letter) const;

    const std::string& getName() const { return m_name; }
    int32 getPntSize() const { return m_pntSize; }

    bool isInitialized() const { return m_texture != nullptr; }

private:
    FT_Face m_face;

    std::string m_name;

    std::vector<Glyph> m_glyphs;
    int32 m_pntSize;

    Texture* m_texture;
    Shader* m_fontShader;

    bool CreateFontTexture();

    int32 m_nRefs;

    Font(const FT_Library& ftLib, const std::string& name, int32 pntSize, Shader* fontShader);
    ~Font();

    friend class FontManager;
};
}
}