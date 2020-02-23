/////////////////////////////////////////////////////////////////////////////
/// Author:      Ahmed Hamdi Boujelben <ahmed.hamdi.boujelben@gmail.com>
/// Created:     2016
/// Copyright:   (c) 2020 Ahmed Hamdi Boujelben
/// Licence:     Attribution-NonCommercial 4.0 International
////////////////////////////////////////////////////////////////////////////

#version 330 core
out vec4 FragColor;

in FragData
{
    flat vec3 cameraSpherePos;
    flat float sphereRadius;
    flat vec3 sphereColor;
    smooth vec2 mapping;
};

uniform vec3 lightPos;
uniform vec3 viewPos;
uniform vec3 lightColor;
uniform mat4 projection;
uniform vec3 phase1Color;
uniform vec3 phase2Color;
uniform vec3 phase3Color;
uniform vec3 phase4Color;

void Impostor(out vec3 cameraPos, out vec3 cameraNormal)
{
    float lensqr = dot(mapping, mapping);
    if(lensqr > 1.0)
        discard;

    cameraNormal = vec3(mapping, sqrt(1.0 - lensqr));
    cameraPos = (cameraNormal * sphereRadius) + cameraSpherePos;
}

void main()
{
    vec3 cameraPos;
    vec3 cameraNormal;
    float alpha;

    vec3 objectColor;
    if(sphereColor.x==1)
        objectColor=phase1Color+sphereColor.y*(phase4Color-phase1Color);
    if(sphereColor.x==2)
        objectColor=phase2Color+sphereColor.y*(phase4Color-phase2Color);
    if(sphereColor.x==3)
        objectColor=phase3Color+sphereColor.y*(phase4Color-phase3Color);


    Impostor(cameraPos, cameraNormal);

    //Set the depth based on the new cameraPos.
    vec4 clipPos = projection * vec4(cameraPos, 1.0);
    float ndcDepth = clipPos.z / clipPos.w;
    gl_FragDepth = ((gl_DepthRange.diff * ndcDepth) +
        gl_DepthRange.near + gl_DepthRange.far) / 2.0;

    // ambient
    float ambientStrength = 0.3;
    vec3 ambient = ambientStrength * lightColor;

    //vec3 lightPosView= vec3(view*vec4(lightPos,1.0));
    // diffuse
    vec3 norm = normalize(cameraNormal);
    vec3 lightDir = normalize(lightPos - cameraPos);
    float diff = max(dot(norm, lightDir), 0.0);
    vec3 diffuse = diff * lightColor;

    // specular
    float specularStrength = 1.0;
    vec3 viewDir = normalize(-cameraPos);
    vec3 reflectDir = reflect(-lightDir, norm);
    float spec = pow(max(dot(viewDir, reflectDir), 0.0), 32);
    vec3 specular = specularStrength * spec * lightColor;

    vec3 result = (ambient + diffuse + specular) * objectColor;
    FragColor = vec4(result, 1.0);
}
