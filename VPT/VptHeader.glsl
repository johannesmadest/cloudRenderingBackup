/**
 * MIT License
 *
 * Copyright (c) 2021-2022, Christoph Neuhauser, Timm Knörle, Ludwig Leonard
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

layout (local_size_x = LOCAL_SIZE_X, local_size_y = LOCAL_SIZE_Y, local_size_z = 1) in;

layout(push_constant) uniform Push
{
    ivec2 tileOffset;
    ivec2 tileSize;
} pc;

layout (binding = 0, rgba32f) uniform image2D resultImage;

#ifdef USE_NANOVDB
layout (binding = 1) readonly buffer NanoVdbBuffer {
    uint pnanovdb_buf_data[];
};
#ifdef USE_EMISSION
layout (binding = 2) readonly buffer EmissionNanoVdbBuffer {
    uint pnanovdb_emission_buf_data[];
};
#endif
#else // USE_NANOVDB
layout (binding = 1) uniform sampler3D gridImage;
#ifdef USE_EMISSION
layout (binding = 2) uniform sampler3D emissionImage;
#endif
#endif // USE_NANOVDB

layout (binding = 3) uniform Parameters {
    // Transform from normalized device coordinates to world space.
    mat4 inverseViewProjMatrix;
    mat4 previousViewProjMatrix;
    mat4 inverseTransposedViewMatrix;
    mat4 inverseViewMatrix;
    mat4 viewMatrix;

    // Cloud properties.
    vec3 boxMin; float voxelValueMin;
    vec3 boxMax; float voxelValueMax;
    vec3 gridMin; float minGradientVal;
    vec3 gridMax; float maxGradientVal;
    vec3 emissionBoxMin;
    vec3 emissionBoxMax;
    vec3 extinction; float tfScatteringAlbedoStrength;
    vec3 scatteringAlbedo;

    float phaseG;

    // Sky properties.
    vec3 sunDirection;
    vec3 sunIntensity;
    vec3 environmentMapIntensityFactor;
    mat3 envMapDirRot; //< Environment map sampling direction rotation matrix (mat3).
    mat3 invEnvMapDirRot; //< Inverse of matrix above.

    float emissionCap;
    float emissionStrength;
    int numFeatureMapSamplesPerFrame;

    // Whether to use linear RGB or sRGB.
    int useLinearRGB;

    // For residual ratio tracking and decomposition tracking.
    ivec3 superVoxelSize;
    ivec3 superVoxelGridSize;

    ivec3 gridResolution;
    uint isEnvMapBlack; //< Info for sampling environment map vs. single light sources.
    vec3 voxelTexelSize;
    float farDistance;

    // Isosurfaces.
    vec3 isoSurfaceColor;
    float isoValue;
    float isoStepWidth;
    float maxAoDist;
    int numAoSamples;

    // Clip plane
    int useClipPlane;
    vec3 clipPlaneNormal;
    float clipPlaneDistance;

    // Disney BRDF
    float subsurface;
    float metallic;
    float specular;
    float specularTint;
    float roughness;
    float anisotropic;
    float sheen;
    float sheenTint;
    vec3 camForward;
    float clearcoat;
    float clearcoatGloss;
} parameters;

layout (binding = 4) uniform FrameInfo {
    uint frameCount;
    // Either equivalent to frameNumber or a global frame ID not reset together with accumulation.
    uint globalFrameNumber;
    uvec2 other;
} frameInfo;

layout (binding = 5, rgba32f) uniform image2D accImage;

#ifdef WRITE_POSITION_MAP
layout (binding = 6, rgba32f) uniform image2D firstX;
#endif

#ifdef WRITE_FIRST_W_MAP
layout (binding = 7, rgba32f) uniform image2D firstW;
#endif

#if !defined(USE_NANOVDB) && (defined(USE_RESIDUAL_RATIO_TRACKING) || defined(USE_DECOMPOSITION_TRACKING))
layout (binding = 8) uniform sampler3D superVoxelGridImage;
layout (binding = 9) uniform usampler3D superVoxelGridOccupancyImage;
#endif

#ifdef USE_ENVIRONMENT_MAP_IMAGE
layout (binding = 10) uniform sampler2D environmentMapTexture;
layout (binding = 11) uniform sampler2D environmentMapOctahedralTexture;
#endif

#ifdef COMPUTE_PRIMARY_RAY_ABSORPTION_MOMENTS
layout (binding = 12, r32f) uniform image2DArray primaryRayAbsorptionMomentsImage;
#endif

#ifdef COMPUTE_SCATTER_RAY_ABSORPTION_MOMENTS
layout (binding = 13, r32f) uniform image2DArray scatterRayAbsorptionMomentsImage;
#endif

#ifdef WRITE_CLOUDONLY_MAP
layout (binding = 14, rgba32f) uniform image2D cloudOnlyImage;
#endif

#ifdef WRITE_DEPTH_MAP
layout (binding = 15, rg32f) uniform image2D depthImage;
#endif

#ifdef WRITE_DENSITY_MAP
layout (binding = 16, rg32f) uniform image2D densityImage;
#endif

#ifdef WRITE_BACKGROUND_MAP
layout (binding = 17, rg32f) uniform image2D backgroundImage;
#endif

#ifdef WRITE_REPROJ_UV_MAP
layout (binding = 18, rg32f) uniform image2D reprojUVImage;
#endif

#ifdef WRITE_NORMAL_MAP
layout (binding = 19, rgba32f) uniform image2D normalImage;
#endif

#ifdef WRITE_DEPTH_BLENDED_MAP
layout (binding = 20, rg32f) uniform image2D depthBlendedImage;
#endif

#ifdef WRITE_DEPTH_NEAREST_OPAQUE_MAP
layout (binding = 21, rg32f) uniform image2D depthNearestOpaqueImage;
#endif

#ifdef WRITE_FLOW_MAP
layout(binding = 22, rg32f) uniform image2D flowImage;
#endif

#ifdef WRITE_DEPTH_NABLA_MAP
layout (binding = 23, rg32f) uniform image2D depthNablaImage;
#endif

#ifdef WRITE_DEPTH_FWIDTH_MAP
layout (binding = 24, r32f) uniform image2D depthFwidthImage;
#endif

#ifdef WRITE_ALBEDO_MAP
layout (binding = 25, rgba32f) uniform image2D albedoImage;
#endif

#ifdef WRITE_TRANSMITTANCE_VOLUME
layout (binding = 26, r32ui) uniform uimage3D transmittanceVolumeImage;
#endif


/**
 * This code is part of an GLSL port of the HLSL code accompanying the paper "Moment-Based Order-Independent
 * Transparency" by Münstermann, Krumpen, Klein, and Peters (http://momentsingraphics.de/?page_id=210).
 * The original code was released in accordance to CC0 (https://creativecommons.org/publicdomain/zero/1.0/).
 *
 * This port is released under the terms of the MIT License.
 */
/*! This function implements complex multiplication.*/
layout(std140, binding = 27) uniform MomentUniformData {
    vec4 wrapping_zone_parameters;
    //float overestimation;
    //float moment_bias;
};
const float ABSORBANCE_MAX_VALUE = 10.0;

#ifdef USE_TRANSFER_FUNCTION
layout(binding = 28) uniform sampler1DArray transferFunctionTexture;
#endif

#if defined(ISOSURFACE_TYPE_GRADIENT) || (defined(ISOSURFACE_TYPE_DENSITY) && defined(ISOSURFACE_USE_TF) && defined(USE_TRANSFER_FUNCTION))
layout(binding = 29) uniform sampler3D gradientImage;
#endif

#ifdef USE_OCCUPANCY_GRID
// The occupancy grid can be used for empty space skipping with next event tracking transmittance rays.
layout(binding = 30, r8ui) uniform readonly uimage3D occupancyGridImage;
#endif

#if NUM_LIGHTS > 0
// Light::lightType
const uint LIGHT_TYPE_POINT = 0u;
const uint LIGHT_TYPE_SPOT = 1u;
const uint LIGHT_TYPE_DIRECTIONAL = 2u;
// Light::lightSpace
const uint LIGHT_SPACE_WORLD = 0u; // Positions & directions are in world space.
const uint LIGHT_SPACE_VIEW = 1u; // Positions & directions are in view space.
const uint LIGHT_SPACE_VIEW_ORIENTATION = 2u; // Positions are in world space, directions in view space.
struct Light {
    uint lightType;
    uint lightSpace; ///< All types; world space or view space position/direction?
    float spotTotalWidth; ///< SPOT
    float spotFalloffStart; ///< SPOT

    vec3 color; ///< All types
    float intensity; ///< All types

    // Point light & spotlight.
    vec3 position; ///< POINT & SPOT; distance for DIRECTIONAL.
    uint useDistance; ///< POINT & SPOT

    vec3 spotDirection;
    float padding;
};
layout (binding = 31) uniform LightsBuffer {
    Light lights[NUM_LIGHTS];
};
#endif

//buffers and structs for radiance caching
const uint MAX_LINEAR_ENTRIES = 512;
const uint MAX_RESULT_ENTRIES = 128;

const float SINGLE_MIN_RADIUS = 0.0001;
const float SINGLE_MAX_RADIUS = 0.01;

const float MULTI_MIN_RADIUS = 0.0001;
const float MULTI_MAX_RADIUS = 0.2;

const ivec2 debugPixel = ivec2(128, 128);



//layout for isotropic media
struct LinearNode
{
    vec4 location;
    vec4 radiance;
    vec4 gradientR;
    vec4 gradientG;
    vec4 gradientB;
    vec4 validRadius;
    uint hasValue;
    uint padding[3];
};

layout(std430, binding = 34) buffer LinearBuffer
{   
    LinearNode nodes[];
    
} linearBuffer;

layout(std430, binding = 35) buffer LinearMultiBuffer
{
    LinearNode nodes[];
} linearMultiBuffer;

layout(std430, binding = 36) buffer LinearDirectBuffer
{
    LinearNode nodes[];
} linearDirectBuffer;

//layout for anisotropic media

struct NodeSH
{
    vec4 locationAndRadius; //xyz is location, w is radius
    float coefficients[28]; // 9 coefficents per channel 27 + 1 padding
    vec4 gradients[27]; // 9 per channel
};

struct NodeSH2
{
    vec4 location;
    float coeffR[12];
    float coeffG[12]; //9 coefficients +3 padding for each channel
    float coeffB[12];
    vec4 gradient[3];
    vec4 validRadius;
};

layout(std430, binding = 37) buffer SingleSHBuffer
{   
    NodeSH nodes[];
} singleSHBuffer;

layout(std430, binding = 38) buffer MultiSHBuffer
{
    NodeSH nodes[];
} multiSHBuffer;



//counter buffers

layout(std430, binding = 39) buffer LinearCounter
{
uint nextFreeNode;
} linearCounter;




layout(std430, binding = 40) buffer LinearMultiCounter
{
uint nextFreeNode;
} linearMultiCounter;

layout(std430, binding = 41) buffer LinearDirectCounter
{
uint nextFreeNode;
} linearDirectCounter;


//debug buffers
struct Result{
    vec4 location;
    vec4 singleRadiance;
    vec4 multiRadiance;
    vec4 metaData; // x = step indicator 
};

layout(std430, binding = 42) buffer ResultBuffer
{
    Result results[];
} resultBuffer;

layout(std430, binding = 43) buffer ResultCounterBuffer
{
    uint nextFreeNode;
}resultCounterBuffer;



vec2 Multiply(vec2 LHS, vec2 RHS) {
    return vec2(LHS.x * RHS.x - LHS.y * RHS.y, LHS.x * RHS.y + LHS.y * RHS.x);
}
