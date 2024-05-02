#version 330

uniform float screen_width;
uniform float screen_height;
uniform vec3 color;

in vec2 position;
in vec2 tex_coord;

out vec2 tex_coord_;
out vec4 color_;

void main() {
    tex_coord_ = tex_coord;
    color_ = vec4(color, 1.0);
	gl_Position = vec4(position.x * 2.0f / screen_width - 1.0f,
                       position.y * 2.0f / screen_height - 1.0f,
                       0.0f,
                       1.0f);
}