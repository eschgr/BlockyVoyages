#version 330

uniform sampler2DArray tex;
in vec2 tex_coord_;
in vec4 color_;
flat in int type_;

out vec4 frag_color;

void main() {
    vec4 texture_color = texture(tex, vec3(tex_coord_.st, type_));

    frag_color = color_ * texture_color;
}