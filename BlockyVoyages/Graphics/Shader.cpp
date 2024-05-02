#include "stdafx.h"

#include "Shader.h"

#include "OpenGL.h"
#include "../Debugability/Logger.h"

#include <fstream>

namespace BlockyVoyages {
namespace Graphics {

Shader::Shader() 
    : m_vs_handle(0),
      m_fs_handle(0),
      m_program(0)
{}

Shader::Shader(const std::string& vs_file, const std::string& fs_file)
    : m_vs_handle(0),
      m_fs_handle(0),
      m_program(0)
{
    loadFromFile(vs_file, fs_file);
}

Shader::~Shader() {
    if (0 != m_program) {
        glDetachShader(m_program, m_vs_handle);
        glDetachShader(m_program, m_fs_handle);
        glDeleteShader(m_vs_handle);
        glDeleteShader(m_fs_handle);
        glDeleteProgram(m_program);
    }
}

// Loads a shader from file and compiles them. Returns true on success.
bool Shader::loadFromFile(const std::string& vs_file, const std::string& fs_file) {
    std::string vs_text, fs_text;
    if (!loadShaderFile(vs_file, vs_text)) {
        return false;
    }
    if (!loadShaderFile(fs_file, fs_text)) {
        return false;
    }
    bool result = loadFromString(vs_text, fs_text);
    if (!result) {
        LOG(Debugability::DBG_LVL_MEDIUM, "Shaders failed to load");
        LOG(Debugability::DBG_LVL_MEDIUM, "   Vertex shader: '%s'", vs_file.c_str());
        LOG(Debugability::DBG_LVL_MEDIUM, "   Fragment shader: '%s'", fs_file.c_str());
    }
    return result;
}

// Loads a shader from string and compiles them. Returns true on success.
bool Shader::loadFromString(const std::string& vs_str, const std::string& fs_str) {
    // load the individual shaders
    if (!loadShader(vs_str, GL_VERTEX_SHADER, m_vs_handle)) {
        LOG(Debugability::DBG_LVL_MEDIUM, "Failed to load vertex shader.");
        reportShaderError(m_vs_handle);
        return false;
    }
    if (!loadShader(fs_str, GL_FRAGMENT_SHADER, m_fs_handle)) {
        LOG(Debugability::DBG_LVL_MEDIUM, "Failed to load fragment shader.");
        reportShaderError(m_fs_handle);
        return false;
    }
    // create and link the program
    m_program = glCreateProgram();
    glAttachShader(m_program, m_vs_handle);
    glAttachShader(m_program, m_fs_handle);
    
    int32 glError;
    if ((glError = glGetError()) != 0) {
        LOG(Debugability::DBG_LVL_HIGH, "Got OpenGL error %d.", glError);
    }

    glLinkProgram(m_program);

    int32 status = 0;
    glGetProgramiv(m_program, GL_LINK_STATUS, &status);

    if (0 == status) {
        LOG(Debugability::DBG_LVL_MEDIUM, "Program link failed.");
        reportLinkerError();
        return false;
    }
    getUniformLocations();
    return true;
}

int32 Shader::getAttribLocation(const std::string& attrib_name) const {
    return glGetAttribLocation(m_program, attrib_name.c_str());
}

void Shader::bind() {
    glUseProgram(m_program);
}

void Shader::setUniform(const std::string& var_name, float32 value) {
    auto uniform = m_uniformHandles.find(var_name);
    if (uniform != m_uniformHandles.end()) {
        glUniform1f(uniform->second, value);
    }
}

void Shader::setUniform(const std::string& var_name, Vector2f value) {
    auto uniform = m_uniformHandles.find(var_name);
    if (uniform != m_uniformHandles.end()) {
        glUniform2f(uniform->second, value.x, value.y);
    }
}

void Shader::setUniform(const std::string& var_name, Vector3f value) {
    auto uniform = m_uniformHandles.find(var_name);
    if (uniform != m_uniformHandles.end()) {
        glUniform3f(uniform->second, value.x, value.y, value.z);
    }
}

void Shader::setUniform(const std::string& var_name, Vector4f value) {
    auto uniform = m_uniformHandles.find(var_name);
    if (uniform != m_uniformHandles.end()) {
        glUniform4f(uniform->second, value.x, value.y, value.z, value.w);
    }
}

void Shader::setUniform(const std::string& var_name, int32 value) {
    auto uniform = m_uniformHandles.find(var_name);
    if (uniform != m_uniformHandles.end()) {
        glUniform1i(uniform->second, value);
    }
}

void Shader::setUniform(const std::string& var_name, Vector2i value) {
    auto uniform = m_uniformHandles.find(var_name);
    if (uniform != m_uniformHandles.end()) {
        glUniform2i(uniform->second, value.x, value.y);
    }
}

void Shader::setUniform(const std::string& var_name, Vector3i value) {
    auto uniform = m_uniformHandles.find(var_name);
    if (uniform != m_uniformHandles.end()) {
        glUniform3i(uniform->second, value.x, value.y, value.z);
    }
}

void Shader::setUniform(const std::string& var_name, Vector4i value) {
    auto uniform = m_uniformHandles.find(var_name);
    if (uniform != m_uniformHandles.end()) {
        glUniform4i(uniform->second, value.x, value.y, value.z, value.w);
    }
}

void Shader::setUniform(const std::string& var_name, const Matrix22f& matrix) {
    auto uniform = m_uniformHandles.find(var_name);
    if (uniform != m_uniformHandles.end()) {
        glUniformMatrix2fv(uniform->second, 1, false, matrix.pntr);
    }
}

void Shader::setUniform(const std::string& var_name, const Matrix33f& matrix) {
    auto uniform = m_uniformHandles.find(var_name);
    if (uniform != m_uniformHandles.end()) {
        glUniformMatrix3fv(uniform->second, 1, false, matrix.pntr);
    }
}

void Shader::setUniform(const std::string& var_name, const Matrix44f& matrix) {
    auto uniform = m_uniformHandles.find(var_name);
    if (uniform != m_uniformHandles.end()) {
        glUniformMatrix4fv(uniform->second, 1, false, matrix.pntr);
    }
}

bool Shader::loadShaderFile(const std::string& filename, std::string& out_str) {
    std::ifstream in_file(filename);
    if (!in_file) {
        return false;
    }
    out_str.assign((std::istreambuf_iterator<char>(in_file)),
                   (std::istreambuf_iterator<char>()));
    return true;
}

bool Shader::loadShader(const std::string& shader, GLenum shader_type, uint32& shader_handle) {
    shader_handle = glCreateShader(shader_type);
    const char* shader_text = shader.c_str();
    glShaderSource(shader_handle, 1, &shader_text, nullptr);
    glCompileShader(shader_handle);

    int32 status = 0;
    glGetShaderiv(shader_handle, GL_COMPILE_STATUS, &status);
    return 0 != status;
}

void Shader::getUniformLocations() {
    m_uniformHandles.clear();
    int uniform_count = 0;
    glGetProgramiv(m_program, GL_ACTIVE_UNIFORMS, &uniform_count);
    for (uint32 i = 0; i < static_cast<uint32>(uniform_count); ++i)  {
        int32 name_len;
        int32 size;
        GLenum type;
        char name[100];
        glGetActiveUniform(m_program, i, sizeof(name) - 1, &name_len, &size, &type, name);
        name[name_len] = 0;
        m_uniformHandles[name] = glGetUniformLocation(m_program, name);
    }
}

void Shader::reportShaderError(uint32 object_id) {
    int info_log_length = 0;
    int chars_written  = 0;

    glGetShaderiv(object_id, GL_INFO_LOG_LENGTH, &info_log_length);
    if (info_log_length > 0) {
        auto* info_log = new char[info_log_length];
        glGetShaderInfoLog(object_id, info_log_length, &chars_written, info_log);
	    LOG(Debugability::DBG_LVL_NORMAL, "%s", info_log);
        delete[] info_log;
    }
}

void Shader::reportLinkerError() {
    int info_log_length = 0;
    int chars_written  = 0;

    glGetProgramiv(m_program, GL_INFO_LOG_LENGTH, &info_log_length);
    if (info_log_length > 0) {
        auto* info_log = new char[info_log_length];
        glGetProgramInfoLog(m_program, info_log_length, &chars_written, info_log);
	    LOG(Debugability::DBG_LVL_NORMAL, "%s", info_log);
        delete[] info_log;
    }
}

}
}