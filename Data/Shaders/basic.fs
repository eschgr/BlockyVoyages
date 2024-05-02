#version 330

uniform sampler2D tex;
in vec2 tex_coord_;
in vec4 color_;

out vec4 frag_color;

void main() {
    vec4 texture_color = texture(tex, tex_coord_.st);

    frag_color = color_ * texture_color;
}