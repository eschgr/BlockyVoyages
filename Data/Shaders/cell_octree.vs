#version 330
#extension GL_ARB_draw_instanced : enable

uniform mat4 view_proj_mat;
uniform vec4 color;

in vec3 position;

out vec4 color_;

void main() {
    color_ = color;
    gl_Position = view_proj_mat * vec4(position, 1.0f);
}