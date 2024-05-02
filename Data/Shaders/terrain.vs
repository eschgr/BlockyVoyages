#version 330
#extension GL_ARB_draw_instanced : enable

uniform mat4 view_proj_mat;

in vec3 position;
in vec3 normal;
in vec2 tex_coord;

in vec3 offset;
in float scale;
in int type;

out vec2 tex_coord_;
out vec4 color_;
flat out int type_;

void main() {
    tex_coord_ = tex_coord;
    color_ = vec4(normal, 1.0f);
    gl_Position = view_proj_mat * vec4(position * scale + offset, 1.0f);
    type_ = type;
}