#version 330

out vec4 frag_color;

uniform vec3 biome_color;

void main() {
    frag_color = vec4(biome_color, 1.0f);
}