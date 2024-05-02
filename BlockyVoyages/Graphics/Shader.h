#pragma once

#include "../Math/MathMain.h"

#include <map>
#include <string>

namespace BlockyVoyages {
namespace Graphics {

// These shaders don't support tesselation or geometry shaders. That may come
// in a later update.

class Shader {
public:
    Shader();
    Shader(const std::string& vs_file, const std::string& fs_file);
    virtual ~Shader();

    // Loads a shader from file and compiles them. Returns true on success.
    bool loadFromFile(const std::string& vs_file, const std::string& fs_file);

    // Loads a shader from string and compiles them. Returns true on success.
    bool loadFromString(const std::string& vs_str, const std::string& fs_str);

    int32 getAttribLocation(const std::string& attrib_name) const;

    void bind();

    // for not often set uniforms. For those that are set extremely often, it
    // may be faster to implement a way to get the variable location and then
    // pass that in often.
    void setUniform(const std::string& var_name, float32 value);
    void setUniform(const std::string& var_name, Vector2f value);
    void setUniform(const std::string& var_name, Vector3f value);
    void setUniform(const std::string& var_name, Vector4f value);

    void setUniform(const std::string& var_name, int32 value);
    void setUniform(const std::string& var_name, Vector2i value);
    void setUniform(const std::string& var_name, Vector3i value);
    void setUniform(const std::string& var_name, Vector4i value);

    void setUniform(const std::string& var_name, const Matrix22f& matrix);
    void setUniform(const std::string& var_name, const Matrix33f& matrix);
    void setUniform(const std::string& var_name, const Matrix44f& matrix);

private:
    uint32 m_vs_handle;
    uint32 m_fs_handle;
    uint32 m_program;

    std::map<std::string, uint32> m_uniformHandles;

    // A helper function to load a whole text file at once.
    bool loadShaderFile(const std::string& filename, std::string& out_str);
    bool loadShader(const std::string& shader, GLenum shader_type, uint32& shader_handle);

    void getUniformLocations();

    void reportShaderError(uint32 object_id);
    void reportLinkerError();

    // prevent copy
    Shader(const Shader& other);
    const Shader& operator=(const Shader& other);
};

}
}