#version 330

uniform mat4 model_mat;
uniform mat4 view_proj_mat;

in vec4 position;
in vec3 normal;
in vec2 tex_coord;
in vec3 color;

out vec2 tex_coord_;
out vec4 color_;

void main() {
    tex_coord_ = tex_coord;
    color_ = vec4(color, 1.0);
	gl_Position = view_proj_mat * model_mat * position;
}