/////////////////////////////////////////////////////////////////////////////
/// Author:      Ahmed Hamdi Boujelben <ahmed.hamdi.boujelben@gmail.com>
/// Created:     2016
/// Copyright:   (c) 2020 Ahmed Hamdi Boujelben
/// Licence:     Attribution-NonCommercial 4.0 International
////////////////////////////////////////////////////////////////////////////


#version 330 core
layout(location = 0) in vec3 cubePos;
layout(location = 1) in vec3 cubeSize;
layout(location = 2) in vec3  cubeColor;

out VertexData
{
    vec3 cubePos;
    vec3 cubeSize;
    vec3 cubeColor;
} outData;

void main()
{
    outData.cubePos = cubePos;
    outData.cubeSize = cubeSize;
    outData.cubeColor = cubeColor;
}
