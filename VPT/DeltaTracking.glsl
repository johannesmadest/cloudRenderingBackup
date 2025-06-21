/**
 * MIT License
 *
 * Copyright (c) 2021-2022, Christoph Neuhauser, Ludwig Leonard
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

// Pathtracing with Delta tracking and Spectral tracking.

#ifdef USE_SPECTRAL_DELTA_TRACKING
/**
 * For more details on spectral delta tracking, please refer to:
 * P. Kutz, R. Habel, Y. K. Li, and J. Novák. Spectral and decomposition tracking for rendering heterogeneous volumes.
 * ACM Trans. Graph., 36(4), Jul. 2017.
 */
vec3 deltaTrackingSpectral(vec3 x, vec3 w, inout ScatterEvent firstEvent) {
    float majorant = maxComponent(parameters.extinction);

    vec3 weights = vec3(1, 1, 1);
#ifdef USE_ISOSURFACES
    float lastScalarSign, currentScalarSign;
    bool isFirstPoint = true;
#endif
#ifdef CLOSE_ISOSURFACES
    bool isFirstPointFromOutside = true;
#endif

#if !defined(USE_TRANSFER_FUNCTION) || !defined(USE_TRANSFER_FUNCTION_SCATTERING_ALBEDO)
    vec3 absorptionAlbedo = vec3(1, 1, 1) - parameters.scatteringAlbedo;
    vec3 scatteringAlbedo = parameters.scatteringAlbedo;
    float PA = maxComponent(absorptionAlbedo * parameters.extinction);
    float PS = maxComponent(scatteringAlbedo * parameters.extinction);
#endif

    int i = 0;
    float tMin, tMax;
    if (rayBoxIntersect(parameters.boxMin, parameters.boxMax, x, w, tMin, tMax)) {
        x += w * tMin;
        float d = tMax - tMin;
        while (true) {
#ifdef USE_ISOSURFACES
            i++;
            if (i == 1000) {
                return vec3(0.0, 0.0, 0.0);
            }
#endif

            float t = -log(max(0.0000000001, 1 - random()))/majorant;

            if (t > d) {
                break;
            }

            vec3 xNew = x + w * t;

#ifdef USE_TRANSFER_FUNCTION
            vec4 densityEmission = sampleCloudDensityEmission(xNew);
            float density = densityEmission.a;
#ifdef USE_TRANSFER_FUNCTION_SCATTERING_ALBEDO
            vec3 scatteringAlbedo = mix(parameters.scatteringAlbedo, densityEmission.rgb, parameters.tfScatteringAlbedoStrength);
            vec3 absorptionAlbedo = vec3(1) - scatteringAlbedo;
            float PA = maxComponent(absorptionAlbedo * parameters.extinction);
            float PS = maxComponent(scatteringAlbedo * parameters.extinction);
#endif
#else
            float density = sampleCloud(xNew);
#endif

#include "CheckIsosurfaceHit.glsl"

            x = xNew;

            vec3 sigma_a = absorptionAlbedo * parameters.extinction * density;
            vec3 sigma_s = scatteringAlbedo * parameters.extinction * density;
            vec3 sigma_n = vec3(majorant) - parameters.extinction * density;

#if defined(MAX_BASED_PROBABILITY)
            float Pa = maxComponent(sigma_a);
            float Ps = maxComponent(sigma_s);
            float Pn = maxComponent(sigma_n);
#elif defined(AVG_BASED_PROBABILITY)
            float Pa = avgComponent(sigma_a);
            float Ps = avgComponent(sigma_s);
            float Pn = avgComponent(sigma_n);
#else // Path history average-based probability
            float Pa = avgComponent(sigma_a * weights);
            float Ps = avgComponent(sigma_s * weights);
            float Pn = avgComponent(sigma_n * weights);
#endif
            float C = Pa + Ps + Pn;
            Pa /= C;
            Ps /= C;
            Pn /= C;

            float xi = random();

            if (xi < Pa) {
                if (!firstEvent.hasValue) {
                    firstEvent.x = x;
                    firstEvent.pdf_x = 0; // TODO
                    firstEvent.w = vec3(0.);
                    firstEvent.pdf_w = 0;
                    firstEvent.hasValue = true;
                    firstEvent.density = density * maxComponent(parameters.extinction);
                    firstEvent.depth = tMax - d + t;
                }

#ifdef USE_TRANSFER_FUNCTION
                vec3 L_e = parameters.emissionStrength * densityEmission.rgb;
                return weights * sigma_a / (majorant * Pa) * L_e;
#else
                return vec3(0); // weights * sigma_a / (majorant * Pa) * L_e; // 0 - No emission
#endif
            }

            if (xi < Pa + Ps) { // scattering event
                float pdf_w;
                w = importanceSamplePhase(parameters.phaseG, w, pdf_w);

                if (!firstEvent.hasValue) {
                    firstEvent.x = x;
                    firstEvent.pdf_x = 0; // TODO
                    firstEvent.w = w;
                    firstEvent.pdf_w = pdf_w;
                    firstEvent.hasValue = true;
                    firstEvent.density = density * maxComponent(parameters.extinction);
                    firstEvent.depth = tMax - d + t;
                }

                if (rayBoxIntersect(parameters.boxMin, parameters.boxMax, x, w, tMin, tMax)) {
                    x += w*tMin;
                    d = tMax - tMin;
                }
                weights *= sigma_s / (majorant * Ps);
            } else {
                d -= t;
                weights *= sigma_n / (majorant * Pn);
            }
#if !defined(MAX_BASED_PROBABILITY) && !defined(AVG_BASED_PROBABILITY)
            weights = min(weights, vec3(100.0, 100.0, 100.0));
#endif
        }
    }

    return min(weights, vec3(100000, 100000, 100000)) * (sampleSkybox(w) + sampleLight(w));
}
#endif


#ifdef USE_DELTA_TRACKING

void insertResult(vec3 position, vec3 singleRadiance, vec3 multiRadiance, vec4 metaData)
    {
        
        if(resultCounterBuffer.nextFreeNode > 127) return;

        uint index = atomicAdd(resultCounterBuffer.nextFreeNode, 1);
        resultBuffer.results[index].location = vec4(position, 0);
        resultBuffer.results[index].singleRadiance = vec4(singleRadiance, 0);
        resultBuffer.results[index].multiRadiance = vec4(multiRadiance, 0);
        resultBuffer.results[index].metaData = metaData;
    }




vec3 deltaTracking(
        vec3 x, vec3 w, inout ScatterEvent firstEvent
#ifdef COMPUTE_SCATTER_RAY_ABSORPTION_MOMENTS
        , out float scatterRayAbsorptionMoments[NUM_SCATTER_RAY_ABSORPTION_MOMENTS + 1]
#endif
) {
#ifdef COMPUTE_SCATTER_RAY_ABSORPTION_MOMENTS
    for (int i = 0; i <= NUM_SCATTER_RAY_ABSORPTION_MOMENTS; i++) {
        scatterRayAbsorptionMoments[i] = 0.0;
    }
    float depth = 0.0;
#endif

//get variables from parameters
    float majorant = parameters.extinction.x;
    float absorptionAlbedo = 1.0 - parameters.scatteringAlbedo.x;
    float scatteringAlbedo = parameters.scatteringAlbedo.x;
    float PA = absorptionAlbedo * parameters.extinction.x;
    float PS = scatteringAlbedo * parameters.extinction.x;

#ifdef USE_ISOSURFACES
    vec3 weights = vec3(1, 1, 1);
    float lastScalarSign, octreeBuffer.nodes[nodeIndex]ScalarSign;
    bool isFirstPoint = true;
#endif
#ifdef CLOSE_ISOSURFACES
    bool isFirstPointFromOutside = true;
#endif

    int i = 0;
    float tMin, tMax;
    if (rayBoxIntersect(parameters.boxMin, parameters.boxMax, x, w, tMin, tMax)) {
        x += w * tMin;
#ifdef COMPUTE_SCATTER_RAY_ABSORPTION_MOMENTS
        //depth += tMin;
#endif
        float d = tMax - tMin;

        float pdf_x = 1;
        float transmittance = 1.0;

        while (true) {
#ifdef USE_ISOSURFACES
            i++;
            if (i == 1000) {
                return vec3(0.0, 0.0, 0.0);
            }
#endif
#ifdef COMPUTE_SCATTER_RAY_ABSORPTION_MOMENTS
            float absorbance = -log(transmittance);
            if (absorbance > ABSORBANCE_MAX_VALUE) {
                absorbance = ABSORBANCE_MAX_VALUE;
            }
#ifdef USE_POWER_MOMENTS_SCATTER_RAY
            for (int i = 0; i <= NUM_SCATTER_RAY_ABSORPTION_MOMENTS; i++) {
                scatterRayAbsorptionMoments[i] += absorbance * pow(depth, i);
            }
#else
            float phase = fma(depth, wrapping_zone_parameters.y, wrapping_zone_parameters.y);
            vec2 circlePoint = vec2(cos(phase), sin(phase));
            scatterRayAbsorptionMoments[0] = absorbance;
            scatterRayAbsorptionMoments[1] = absorbance * circlePoint.x;
            scatterRayAbsorptionMoments[2] = absorbance * circlePoint.y;
            vec2 circlePointNext = circlePoint;
            for (int i = 2; i <= NUM_SCATTER_RAY_ABSORPTION_MOMENTS / 2; i++) {
                circlePointNext = Multiply(circlePointNext, circlePoint);
                scatterRayAbsorptionMoments[i * 2] = absorbance * circlePointNext.x;
                scatterRayAbsorptionMoments[i * 2 + 1] = absorbance * circlePointNext.y;
            }
#endif
            transmittance = 1.0;
#endif
            float t = -log(max(0.0000000001, 1 - random())) / majorant;

            if (t > d) {
                break;
            }

            vec3 xNew = x + w * t;
#ifdef COMPUTE_SCATTER_RAY_ABSORPTION_MOMENTS
            depth += t;
#endif

#ifdef USE_TRANSFER_FUNCTION
            vec4 densityEmission = sampleCloudDensityEmission(xNew);
            float density = densityEmission.a;
#else
            float density = sampleCloud(xNew);
#endif

#include "CheckIsosurfaceHit.glsl"

            x = xNew;
            transmittance *= 1.0 - density;

            float sigma_a = PA * density;
            float sigma_s = PS * density;
            float sigma_n = majorant - parameters.extinction.x * density;

            float Pa = sigma_a / majorant;
            float Ps = sigma_s / majorant;
            float Pn = sigma_n / majorant;

            float xi = random();

            if (xi < Pa) { //absorption event -- ray ends and only emission is returned
                if (!firstEvent.hasValue) {
                    firstEvent.x = x;
                    firstEvent.pdf_x = 0;
                    firstEvent.w = vec3(0.);
                    firstEvent.pdf_w = 0;
                    firstEvent.hasValue = true;
                    firstEvent.density = density * parameters.extinction.x;
                    firstEvent.depth = tMax - d + t;
                }
#if defined(USE_EMISSION) || defined(USE_TRANSFER_FUNCTION)
                vec3 L_e = vec3(0.0);
#endif
#ifdef USE_EMISSION
                L_e += sampleEmission(x);
#ifdef USE_ISOSURFACES
                return weights * emission;
#else
                return emission;
#endif
#endif
#ifdef USE_TRANSFER_FUNCTION
                L_e += parameters.emissionStrength * densityEmission.rgb;
                //return weights * sigma_a / (majorant * Pa) * L_e;
#endif
#if defined(USE_EMISSION) || defined(USE_TRANSFER_FUNCTION)
#ifdef USE_ISOSURFACES
                return weights * L_e;
#else
                return L_e;
#endif
#else
                return vec3(0); // weights * sigma_a / (majorant * Pa) * L_e; // 0 - No emission
#endif
            }

            if (xi < 1 - Pn) // scattering event -- here i want to compute single scattering and multiple scattering radiances or use cache values to extrapolate
            {
                float pdf_w;
                w = importanceSamplePhase(parameters.phaseG, w, pdf_w);

                pdf_x *= exp(-majorant * t) * majorant * density;

                if (!firstEvent.hasValue) {
                    firstEvent.x = x;
                    firstEvent.pdf_x = sigma_s * pdf_x;
                    firstEvent.w = w;
                    firstEvent.pdf_w = pdf_w;
                    firstEvent.hasValue = true;
                    firstEvent.density = density * parameters.extinction.x;
                    firstEvent.depth = tMax - d + t;
                }

                if (rayBoxIntersect(parameters.boxMin, parameters.boxMax, x, w, tMin, tMax)) {
                    x += w*tMin;
#ifdef COMPUTE_SCATTER_RAY_ABSORPTION_MOMENTS
                    depth += tMin;
#endif
                    d = tMax - tMin;
                }
            } else { // null event 
                pdf_x *= exp(-majorant * t) * majorant * (1 - density);
                d -= t;
            }
        }
    }

#ifdef USE_ISOSURFACES
    return weights * (sampleSkybox(w) + sampleLight(w));
#else
    return sampleSkybox(w) + sampleLight(w);
#endif
}

/*
vec3 deltaTracking(
        vec3 x, vec3 w, inout ScatterEvent firstEvent
#ifdef COMPUTE_SCATTER_RAY_ABSORPTION_MOMENTS
        , out float scatterRayAbsorptionMoments[NUM_SCATTER_RAY_ABSORPTION_MOMENTS + 1]
#endif
) {
#ifdef COMPUTE_SCATTER_RAY_ABSORPTION_MOMENTS
    for (int i = 0; i <= NUM_SCATTER_RAY_ABSORPTION_MOMENTS; i++) {
        scatterRayAbsorptionMoments[i] = 0.0;
    }
    float depth = 0.0;
#endif

    vec3 radiance = vec3(0.0);
    float majorant = parameters.extinction.x;
    float absorptionAlbedo = 1.0 - parameters.scatteringAlbedo.x;
    float scatteringAlbedo = parameters.scatteringAlbedo.x;
    float PA = absorptionAlbedo * parameters.extinction.x;
    float PS = scatteringAlbedo * parameters.extinction.x;

#ifdef USE_ISOSURFACES
    vec3 weights = vec3(1, 1, 1);
    float lastScalarSign, currentScalarSign;
    bool isFirstPoint = true;
#endif
#ifdef CLOSE_ISOSURFACES
    bool isFirstPointFromOutside = true;
#endif

    int i = 0;
    float tMin, tMax;
    if (rayBoxIntersect(parameters.boxMin, parameters.boxMax, x, w, tMin, tMax)) {
        x += w * tMin;
#ifdef COMPUTE_SCATTER_RAY_ABSORPTION_MOMENTS
        //depth += tMin;
#endif
        float d = tMax - tMin;

        float pdf_x = 1;
        float transmittance = 1.0;

        while (true) {
#ifdef USE_ISOSURFACES
            i++;
            if (i == 1000) {
                return vec3(0.0, 0.0, 0.0);
            }
#endif
#ifdef COMPUTE_SCATTER_RAY_ABSORPTION_MOMENTS
            float absorbance = -log(transmittance);
            if (absorbance > ABSORBANCE_MAX_VALUE) {
                absorbance = ABSORBANCE_MAX_VALUE;
            }
#ifdef USE_POWER_MOMENTS_SCATTER_RAY
            for (int i = 0; i <= NUM_SCATTER_RAY_ABSORPTION_MOMENTS; i++) {
                scatterRayAbsorptionMoments[i] += absorbance * pow(depth, i);
            }
#else
            float phase = fma(depth, wrapping_zone_parameters.y, wrapping_zone_parameters.y);
            vec2 circlePoint = vec2(cos(phase), sin(phase));
            scatterRayAbsorptionMoments[0] = absorbance;
            scatterRayAbsorptionMoments[1] = absorbance * circlePoint.x;
            scatterRayAbsorptionMoments[2] = absorbance * circlePoint.y;
            vec2 circlePointNext = circlePoint;
            for (int i = 2; i <= NUM_SCATTER_RAY_ABSORPTION_MOMENTS / 2; i++) {
                circlePointNext = Multiply(circlePointNext, circlePoint);
                scatterRayAbsorptionMoments[i * 2] = absorbance * circlePointNext.x;
                scatterRayAbsorptionMoments[i * 2 + 1] = absorbance * circlePointNext.y;
            }
#endif
            transmittance = 1.0;
#endif
            float t = -log(max(0.0000000001, 1 - random())) / majorant;

            if (t > d) {
                break;
            }

            vec3 xNew = x + w * t;
#ifdef COMPUTE_SCATTER_RAY_ABSORPTION_MOMENTS
            depth += t;
#endif

#ifdef USE_TRANSFER_FUNCTION
            vec4 densityEmission = sampleCloudDensityEmission(xNew);
            float density = densityEmission.a;
#else
            float density = sampleCloud(xNew);
#endif

#include "CheckIsosurfaceHit.glsl"

            x = xNew;
            transmittance *= 1.0 - density;

            float sigma_a = PA * density;
            float sigma_s = PS * density;
            float sigma_n = majorant - parameters.extinction.x * density;

            float Pa = sigma_a / majorant;
            float Ps = sigma_s / majorant;
            float Pn = sigma_n / majorant;

            float xi = random();

            if (xi < Pa) {
                if (!firstEvent.hasValue) {
                    firstEvent.x = x;
                    firstEvent.pdf_x = 0;
                    firstEvent.w = vec3(0.);
                    firstEvent.pdf_w = 0;
                    firstEvent.hasValue = true;
                    firstEvent.density = density * parameters.extinction.x;
                    firstEvent.depth = tMax - d + t;
                }
#if defined(USE_EMISSION) || defined(USE_TRANSFER_FUNCTION)
                vec3 L_e = vec3(0.0);
#endif
#ifdef USE_EMISSION
                L_e += sampleEmission(x);
#ifdef USE_ISOSURFACES
                return weights * emission;
#else
                return emission;
#endif
#endif
#ifdef USE_TRANSFER_FUNCTION
                L_e += parameters.emissionStrength * densityEmission.rgb;
                //return weights * sigma_a / (majorant * Pa) * L_e;
#endif
#if defined(USE_EMISSION) || defined(USE_TRANSFER_FUNCTION)
#ifdef USE_ISOSURFACES
                return weights * L_e;
#else
                return L_e;
#endif
#else
                
               return vec3(0.0);// weights * sigma_a / (majorant * Pa) * L_e; // 0 - No emission
#endif
            }

            if (xi < 1 - Pn) // scattering event
            {
                float pdf_w;
                vec3 w_s = importanceSamplePhase(parameters.phaseG, w, pdf_w);

                pdf_x *= exp(-majorant * t) * majorant * density;

                if (!firstEvent.hasValue) {
                    firstEvent.x = x;
                    firstEvent.pdf_x = sigma_s * pdf_x;
                    firstEvent.w = w;
                    firstEvent.pdf_w = pdf_w;
                    firstEvent.hasValue = true;
                    firstEvent.density = density * parameters.extinction.x;
                    firstEvent.depth = tMax - d + t;
                }


                
                //new for reference
                vec3 x_s = x;
                float sigma_s = PS * density;
                vec3 L_env = sampleSkybox(w_s) + sampleLight(w_s);
                float neeP = evaluatePhase(parameters.phaseG, w, w_s);
                float Tr_nee = calculateTransmittanceTest(x_s, w_s);

                vec3 L1 = (sigma_s/majorant) * neeP * L_env * Tr_nee / pdf_w;

                vec3 lightDir = parameters.sunDirection;
                vec3 L_nee = sampleSkybox(lightDir) + sampleLight(lightDir);
                float neeP2 = evaluatePhase(parameters.phaseG, w, lightDir);
                float Tr_nee2 = calculateTransmittanceTest(xNew, lightDir);
                vec3 neeContribution = (sigma_s / majorant) * neeP2 * L_nee * Tr_nee2;
                return L1 + neeContribution;
                
                
                w = w_s;
                
                if (rayBoxIntersect(parameters.boxMin, parameters.boxMax, x, w, tMin, tMax)) {
                    x += w*tMin;
#ifdef COMPUTE_SCATTER_RAY_ABSORPTION_MOMENTS
                    depth += tMin;
#endif
                    d = tMax - tMin;
                }
            } else {
                pdf_x *= exp(-majorant * t) * majorant * (1 - density);
                d -= t;
            }
        }
    }

#ifdef USE_ISOSURFACES
    return weights * (sampleSkybox(w) + sampleLight(w));
#else
    //if(gl_GlobalInvocationID.xy == ivec2(127, 42)) insertResult(x, radiance, gl_GlobalInvocationID, vec4(0));
    return (sampleSkybox(w) + sampleLight(w));
#endif
}
*/
#endif






#ifdef USE_RAY_MARCHING

//small helper functions
float smoothCubic(float d)
{
    return (3*d*d) - (2*d*d*d);
}



//Functions to write to Buffers


void insertLinear(vec3 position, vec3 radiance, vec3 gradientR, vec3 gradientG, vec3 gradientB, float validRadius)
    {
        //read current state
        uint newIndex = atomicAdd(linearCounter.nextFreeNode, 0);
        
        if(newIndex >= MAX_LINEAR_ENTRIES) return;

        //only increase if allowed
        newIndex = atomicAdd(linearCounter.nextFreeNode, 1);
        

        linearBuffer.nodes[newIndex].location = vec4(position, 0);
        linearBuffer.nodes[newIndex].radiance = vec4(radiance, 0);
        linearBuffer.nodes[newIndex].gradientR = vec4(gradientR, 0);
        linearBuffer.nodes[newIndex].gradientG = vec4(gradientG, 0);
        linearBuffer.nodes[newIndex].gradientB = vec4(gradientB, 0);
        linearBuffer.nodes[newIndex].validRadius.x = validRadius;
        linearBuffer.nodes[newIndex].hasValue = 1;
    }

void insertLinearDirect(vec3 position, vec3 radiance, vec3 gradientR, vec3 gradientG, vec3 gradientB, float validRadius)
    {
        //read current state
        uint newIndex = atomicAdd(linearDirectCounter.nextFreeNode, 0);
        
        if(newIndex >= MAX_LINEAR_ENTRIES) return;

        //only increase if allowed
        newIndex = atomicAdd(linearDirectCounter.nextFreeNode, 1);
        

        linearDirectBuffer.nodes[newIndex].location = vec4(position, 0);
        linearDirectBuffer.nodes[newIndex].radiance = vec4(radiance, 0);
        linearDirectBuffer.nodes[newIndex].gradientR = vec4(gradientR, 0);
        linearDirectBuffer.nodes[newIndex].gradientG = vec4(gradientG, 0);
        linearDirectBuffer.nodes[newIndex].gradientB = vec4(gradientB, 0);
        linearDirectBuffer.nodes[newIndex].validRadius.x = validRadius;
        linearDirectBuffer.nodes[newIndex].hasValue = 1;
    }

void insertLinearMulti(vec3 position, vec3 radiance, vec3 gradientR, vec3 gradientG, vec3 gradientB, float validRadius)
    {
        if(linearMultiCounter.nextFreeNode >= MAX_LINEAR_ENTRIES) return;
        uint newIndex = atomicAdd(linearMultiCounter.nextFreeNode, 1);
        linearMultiBuffer.nodes[newIndex].location = vec4(position, 0);
        linearMultiBuffer.nodes[newIndex].radiance = vec4(radiance, 0);
        linearMultiBuffer.nodes[newIndex].gradientR = vec4(gradientR, 0);
        linearMultiBuffer.nodes[newIndex].gradientG = vec4(gradientG, 0);
        linearMultiBuffer.nodes[newIndex].gradientB = vec4(gradientB, 0);
        linearMultiBuffer.nodes[newIndex].validRadius.x = validRadius;
        linearMultiBuffer.nodes[newIndex].hasValue = 1;
    }


void insertSingleSH(vec3 position, float coefficients[27], vec3 gradients[27], float radius)
{   
    //read current state
    uint newIndex = atomicAdd(linearCounter.nextFreeNode, 0);
        
    if(newIndex >= MAX_LINEAR_ENTRIES) return;

    //only increase if allowed
    newIndex = atomicAdd(linearCounter.nextFreeNode, 1);
    singleSHBuffer.nodes[newIndex].locationAndRadius = vec4(position, radius);
    for(int i=0; i< 27; i++)
    {
        singleSHBuffer.nodes[newIndex].coefficients[i] = coefficients[i]; //r at 0 to 8
        singleSHBuffer.nodes[newIndex].gradients[i] = vec4(gradients[i], 0);
    }
   


    
}

void insertMultiSH(vec3 position, float coeffR[9], float coeffG[9], float coeffB[9], vec3 gradientR[9], vec3 gradientG[9], vec3 gradientB[9], float radius)
{   
    //read current state
    uint newIndex = atomicAdd(linearMultiCounter.nextFreeNode, 0);
        
    if(newIndex >= MAX_LINEAR_ENTRIES) return;

    //only increase if allowed
    newIndex = atomicAdd(linearMultiCounter.nextFreeNode, 1);
    multiSHBuffer.nodes[newIndex].locationAndRadius = vec4(position, radius);
    for(int i=0; i< 9; i++)
    {
        multiSHBuffer.nodes[newIndex].coefficients[i] = coeffR[i]; //r at 0 to 8
        multiSHBuffer.nodes[newIndex].coefficients[9 + i] = coeffG[i]; // g at 9 to 17
        multiSHBuffer.nodes[newIndex].coefficients[18 + i] = coeffB[i]; // b at 18 to 27

        multiSHBuffer.nodes[newIndex].gradients[i] = vec4(gradientR[i], 0);
        multiSHBuffer.nodes[newIndex].gradients[9 + i] = vec4(gradientG[i], 0);
        multiSHBuffer.nodes[newIndex].gradients[18 + i] = vec4(gradientB[i], 0);
    }

    
}



void insertResult(vec3 position, vec3 singleRadiance, vec3 multiRadiance, vec4 metaData)
    {
        
        if(resultCounterBuffer.nextFreeNode > 127) return;

        uint index = atomicAdd(resultCounterBuffer.nextFreeNode, 1);
        resultBuffer.results[index].location = vec4(position, 0);
        resultBuffer.results[index].singleRadiance = vec4(singleRadiance, 0);
        resultBuffer.results[index].multiRadiance = vec4(multiRadiance, 0);
        resultBuffer.results[index].metaData = metaData;
    }


//Query Functions

    int querySingle(vec3 position, out uint validEntries[MAX_LINEAR_ENTRIES], uint limit) {
        uint validEntryCounter=0;
        for(uint i=0; i< limit; i++)
        {
            
            if(distance(linearBuffer.nodes[i].location.xyz, position) < linearBuffer.nodes[i].validRadius.x)
            {
                validEntries[validEntryCounter] = i;
                validEntryCounter++;
            }
                      
        }
            

        return int(validEntryCounter);


    }

     int queryDirect(vec3 position, out uint validEntries[MAX_LINEAR_ENTRIES], uint limit) {
        uint validEntryCounter=0;
        for(uint i=0; i< limit; i++)
        {
            
            if(distance(linearDirectBuffer.nodes[i].location.xyz, position) < linearDirectBuffer.nodes[i].validRadius.x)
            {
                validEntries[validEntryCounter] = i;
                validEntryCounter++;
            }
                      
        }
            

        return int(validEntryCounter);


    }

    
    int querySingleSH(vec3 position, out uint validEntries[MAX_LINEAR_ENTRIES], uint limit) {
        uint validEntryCounter=0;
        for(uint i=0; i< limit; i++)
        {
            if(abs(parameters.phaseG) < 0.001)
            {
                if(distance(linearBuffer.nodes[i].location.xyz, position) < linearBuffer.nodes[i].validRadius.x)
                {
                    validEntries[validEntryCounter] = i;
                    validEntryCounter++;
                }
            }
            else
            {
                if(distance(singleSHBuffer.nodes[i].locationAndRadius.xyz, position) < singleSHBuffer.nodes[i].locationAndRadius.w)
                {
                    validEntries[validEntryCounter] = i;
                    validEntryCounter++;
                }
            }
            
        }
        //if(resultCounterBuffer.nextFreeNode < 128) insertResult(position, vec3(validEntryCounter), vec3(0), vec4(0));
            

        return int(validEntryCounter);


    }

     int queryMulti(vec3 position, out uint validEntries[MAX_LINEAR_ENTRIES], uint limit) {
        uint validEntryCounter=0;
        for(uint i=0; i< limit; i++)
        {
            
            if(distance(linearMultiBuffer.nodes[i].location.xyz, position) < linearMultiBuffer.nodes[i].validRadius.x)
            {
                validEntries[validEntryCounter] = i;
                validEntryCounter++;
            }                      
        }
        //if(resultCounterBuffer.nextFreeNode < 128) insertResult(position, vec3(validEntryCounter), vec3(0), vec4(0));
            

        return int(validEntryCounter);


    }

    int queryMultiSH(vec3 position, out uint validEntries[MAX_LINEAR_ENTRIES], uint limit) {
        uint validEntryCounter=0;
        for(uint i=0; i< limit; i++)
        {
            if(abs(parameters.phaseG) < 0.001)
            {
                if(distance(linearMultiBuffer.nodes[i].location.xyz, position) < linearMultiBuffer.nodes[i].validRadius.x)
                {
                    validEntries[validEntryCounter] = i;
                    validEntryCounter++;
                }
            }
            else
            {
                if(distance(multiSHBuffer.nodes[i].locationAndRadius.xyz, position) < multiSHBuffer.nodes[i].locationAndRadius.w)
                {
                    validEntries[validEntryCounter] = i;
                    validEntryCounter++;
                }
            }
            
        }
        
            

        return int(validEntryCounter);


    }



//spherical harmonics functions
const int M = 9;
const float inv4PI = 1.0 / (4.0 * PI);
const float g = 0.5;
const float A0 = (1.0) * inv4PI;
const float A1 = (3.0 * g) * inv4PI;
const float A2 = (5.0 * g * g) * inv4PI;



void evalRealSH(in vec3 w, out float shBasis[9])
{
    // normalization constants
    const float c0 = 0.28209479177387814;  //  sqrt 1/4pi
    const float c1 = 0.4886025119029199;   // sqrt 3/4pi
    const float c2 = 1.0925484305920792;   // sqrt 15/4pi
    const float c3 = 0.31539156525252005;  // sqrt 5/16pi
    const float c4 = 0.5462742152960396;   // sqrt 15/16pi

    float x = w.x;
    float y = w.y;
    float z = w.z;

    // l = 0, m = 0
    shBasis[0] =  c0;

    // l = 1
    shBasis[1] =  c1 * y;                      //  m = -1
    shBasis[2] =  c1 * z;                      // m = 0        
    shBasis[3] =  c1 * x;                      // m = +1 
    
    // l = 2
    shBasis[4] =  c2 * x * y;                  // m = -2
    shBasis[5] =  c2 * y * z;                  // m = -1
    shBasis[6] =  c3 * (3.0 * z * z - 1.0);     // m = 0
    shBasis[7] =  c2 * x * z;                  // m = +1
    shBasis[8] =  c4 * (x * x - y * y);        // m = +2
}



//extrapolate radiance function

vec3 extrapolateSingleRadiance(uint cacheIndices[MAX_LINEAR_ENTRIES], uint validEntries, vec3 position)
    {
        //vec3 numeratorSum = vec3(0);
        //vec3 denominatorSum = vec3(0);

        vec3 L = vec3(0);
        float weightSum = 0.0;
        float last_u;
        float lastWeight;

        float minLog = 1e+20;
        float maxLog = 1e-20;

        float gradientScalar = 0.3;
        for(int i=0; i< validEntries; i++)
        {
            
            //get entry
            LinearNode node = linearBuffer.nodes[cacheIndices[i]];
            
            vec3 d = position - node.location.xyz;
            float dist = length(d);

            float r = node.validRadius.x;
            //cubic kernel
            float u = clamp(1.0 - dist/r, 0.0, 1.0);
            last_u = u;
            float weight = 3.0 * u * u - 2.0 * u * u * u;  // 3u^2 - 2u^3
            

            //clamp stored radiance so log never blows up
            vec3 L_cached = max(node.radiance.xyz, vec3(1e-4));
            vec3 logL = log(L_cached);

            float gradR = dot(node.gradientR.xyz, d);
            float gradG = dot(node.gradientG.xyz, d);
            float gradB = dot(node.gradientB.xyz, d);


            vec3 Lc = vec3(logL.r + gradientScalar * gradR / L_cached.r, logL.g + gradientScalar * gradG / L_cached.g, logL.b + gradientScalar * gradB / L_cached.b);
            //vec3 Lc = vec3(logL.r, logL.g, logL.b);

            for(int c=0; c<3; c++)
            {
                minLog = min(minLog, Lc[c]);
                maxLog = max(maxLog, Lc[c]);
            }
            
                            

            L+= weight * Lc;
            weightSum += weight;
            lastWeight = weight;
            
            
        }

        vec3 avgLog = (weightSum > 0.0) ? (L / weightSum) : vec3(0.0);

        avgLog = clamp(avgLog, minLog, maxLog);
       
        vec3 result = exp(avgLog);
        
        return result;
        
        
    }

vec3 extrapolateSingleRadianceALT(uint cacheIndices[MAX_LINEAR_ENTRIES], uint validEntries, vec3 position)
    {
        //vec3 numeratorSum = vec3(0);
        //vec3 denominatorSum = vec3(0);

        vec3 L = vec3(0);
        float weightSum = 0.0;
       

        vec3 minLog = vec3(1e+20);
        vec3 maxLog = vec3(1e-20);

        float gradientScalar = 0.3;

        float w[MAX_LINEAR_ENTRIES];
        float W = 0.0;


        //first pass, store all raw wights and store in an array
        for(int i=0; i< validEntries; i++)
        {
            LinearNode node = linearBuffer.nodes[cacheIndices[i]];
            float dist = length(position - node.location.xyz);
            float u = clamp(1.0 - dist / node.validRadius.x, 0.0, 1.0);
            w[i] = 3.0 * u * u - 2.0 * u * u * u;
            W += w[i];
        }

        //if nothing to extrapolate, return
        if(W<= 0.0) return vec3(0);

        //second pass: build the normalized wighted log space sum

        for(int i=0; i< validEntries; i++)
        {
            
            //get entry
            LinearNode node = linearBuffer.nodes[cacheIndices[i]];
            //get weight
            float wi = w[i] / W;

            vec3 Lc = max(node.radiance.xyz, vec3(1e-4));
            vec3 logL = log(Lc);
            vec3 d = position - node.location.xyz;
            vec3 g = vec3( dot(node.gradientR.xyz, d),
                           dot(node.gradientG.xyz, d),
                           dot(node.gradientB.xyz, d)) 
                           * (gradientScalar / Lc);
            vec3 Lci = (logL + g);
            minLog = min(minLog, Lci);
            maxLog = max(maxLog, Lci);
            L += wi * Lci;
            
            
        }

        
        vec3 clampedL = clamp(L, minLog, maxLog);
        vec3 result = exp(clampedL);
        //insertResult(position, vec3(0), result, vec4(W, validEntries, 0, 0));
        return result;
        
        
    }


//with spherical harmonics THIS EXPECTS NORMALIZED DIRECTION VECTOR
void evalSH2(in vec3 dir, out float Y[9])
{
     //l = 0, m = 0
     Y[0] = 0.282095;           // 1/2 * sqrt(1/pi)

     //l = 1, m = -1, 0, 1
     Y[1] = 0.488603 * dir.y;  // sqrt(3/4π) * y
     Y[2] = 0.488603 * dir.z;  //    ”           ” * z
     Y[3] = 0.488603 * dir.x;  //    ”           ” * x

     //l = 2, m = -2, -1, 0, 1, 2
     Y[4] = 1.092548 * dir.x * dir.y;
     Y[5] = 1.092548 * dir.y * dir.z;
     Y[6] = 0.315392 * (3.0*dir.z*dir.z - 1.0);
     Y[7] = 1.092548 * dir.x * dir.z;
     Y[8] = 0.546274 * (dir.x*dir.x - dir.y*dir.y);
}
/*
vec3 extrapolateSingleSH(uint cacheIndices[MAX_LINEAR_ENTRIES], uint validEntries, vec3 position)
{
    
    vec3 L = vec3(0);
    float minLog = 1e+20;
    float maxLog = 1e-20;
    float R[9];
    float G[9];
    float B[9];


    for(int i=0; i< validEntries; i++)
    {
        LinearNode node = linearBuffer.nodes[cacheIndices[i]];
        float dist = length(position - node.location.xyz);
        float u = clamp(1.0 - dist / node.validRadius.x, 0.0, 1.0);
        w[i] = 3.0 * u * u - 2.0 * u * u * u;
        W += w[i];
    }

    for(int i =0; i < validEntries; i++)
    {
        NodeSH node = singleSHBuffer.nodes[cacheIndices[i]];
        vec3 d = position - node.location.xyz;
        float dist = length(d);

        float r = node.validRadius.x;
        //cubic kernel
        float wi = w[i] / W;

        //for first coefficient, use exponential
        float c0R = max(node.coeffR[0], 1e-4);
        vec3 g0R = node.gradient[0];
        float termR = log(c0R) * dot(g0R, d) / c0R; 
        R[0] = wi * termR;

        float c0G = max(node.coeffG[0], 1e-4);
        vec3 g0G = node.gradient[1];
        float termG = log(c0G) * dot(g0G, d) / c0G; 
        G[0] = wi * termG;

        float c0B = max(node.coeffB[0], 1e-4);
        vec3 g0B = node.gradient[2];
        float termB = log(c0B) * dot(g0B, d) / c0B; 
        B[0] = wi * termR;
        
    }
}
*/
vec3 extrapolateSingleSH_Version2(uint cacheIndices[MAX_LINEAR_ENTRIES], uint validEntries, vec3 position, vec3 outDir)
{
    
    vec3 L = vec3(0);
    vec3 minLog = vec3(1e+20);
    vec3 maxLog = vec3(1e-20);
    float weightSum = 0;
    float w[MAX_LINEAR_ENTRIES];
    float W = 0;

    float Y[9];
    evalSH2(outDir, Y);

   


    for(int i=0; i< validEntries; i++)
    {
        NodeSH node = singleSHBuffer.nodes[cacheIndices[i]];
        float dist = length(position - node.locationAndRadius.xyz);
        float u = clamp(1.0 - dist / node.locationAndRadius.w, 0.0, 1.0);
        w[i] = 3.0 * u * u - 2.0 * u * u * u;
        W += w[i];
    }

    for(int i=0; i< validEntries; i++)
    {
        NodeSH node = singleSHBuffer.nodes[cacheIndices[i]];
        vec3 d = position - node.locationAndRadius.xyz;
        vec3 radiance = vec3(0);
        vec3 gradientR = vec3(0);
        vec3 gradientG = vec3(0);
        vec3 gradientB = vec3(0);
        for(int k=0; k<9; k++)
        {
            radiance.r += node.coefficients[k] * Y[k];
            radiance.g += node.coefficients[9 + k] * Y[k];
            radiance.b += node.coefficients[18 + k] * Y[k];
            gradientR += node.gradients[k].xyz * Y[k];
            gradientG += node.gradients[9 + k].xyz * Y[k];
            gradientB += node.gradients[18 + k].xyz * Y[k];
        }
        radiance = max(radiance, vec3(1e-4));
         
        vec3 Lc = vec3(log(radiance.r) + dot(gradientR,d)/radiance.r, log(radiance.g) + dot(gradientG,d)/radiance.g, log(radiance.b) + dot(gradientB,d)/radiance.b);
        //vec3 Lc = log(radiance);
        L += (w[i]/W) * Lc;
        minLog = min(minLog, Lc);
        maxLog = max(maxLog, Lc);


    }
    vec3 clampedL = clamp(L, minLog, maxLog);
    vec3 result = exp(clampedL);
    //insertResult(position, vec3(0), result, vec4(W, validEntries, 0, 0));
    result = clamp(result, vec3(0.0), vec3(1.0));
    return result;

    
    
}

vec3 extrapolateDirectRadiance(uint cacheIndices[MAX_LINEAR_ENTRIES], uint validEntries, vec3 position)
    {
        //vec3 numeratorSum = vec3(0);
        //vec3 denominatorSum = vec3(0);

        vec3 L = vec3(0);
        float weightSum = 0.0;
        float last_u;
        float lastWeight;

        float minLog = 1e+20;
        float maxLog = 1e-20;

        float gradientScalar = 0.3;
        for(int i=0; i< validEntries; i++)
        {
            
            //get entry
            LinearNode node = linearDirectBuffer.nodes[cacheIndices[i]];
            
            vec3 d = position - node.location.xyz;
            float dist = length(d);

            float r = node.validRadius.x;
            //cubic kernel
            float u = clamp(1.0 - dist/r, 0.0, 1.0);
            last_u = u;
            float weight = 3.0 * u * u - 2.0 * u * u * u;  // 3u^2 - 2u^3
            //if(weight < 0.2) continue;
            //weight = pow(weight, 1.3);

            //clamp stored radiance so log never blows up
            vec3 L_cached = max(node.radiance.xyz, vec3(1e-4));
            vec3 logL = log(L_cached);

            float gradR = dot(node.gradientR.xyz, d);
            float gradG = dot(node.gradientG.xyz, d);
            float gradB = dot(node.gradientB.xyz, d);


            vec3 Lc = vec3(logL.r + gradientScalar * gradR / L_cached.r, logL.g + gradientScalar * gradG / L_cached.g, logL.b + gradientScalar * gradB / L_cached.b);
            //vec3 Lc = vec3(logL.r, logL.g, logL.b);

            for(int c=0; c<3; c++)
            {
                minLog = min(minLog, Lc[c]);
                maxLog = max(maxLog, Lc[c]);
            }
            
                            

            L+= weight * Lc;
            weightSum += weight;
            lastWeight = weight;
            
            
        }

        vec3 avgLog = (weightSum > 0.0) ? (L / weightSum) : vec3(0.0);

        avgLog = clamp(avgLog, minLog, maxLog);
       
        vec3 result = exp(avgLog);
       
        return result;
        
        
    }

vec3 extrapolateMultiRadiance(uint cacheIndices[MAX_LINEAR_ENTRIES], uint validEntries, vec3 position)
    {
        

         //vec3 numeratorSum = vec3(0);
        //vec3 denominatorSum = vec3(0);

        vec3 L = vec3(0);
        float weightSum = 0.0;
        float last_u;
        float lastWeight;

        float minLog = 1e+20;
        float maxLog = 1e-20;

        float gradientScalar = 0.3;
        for(int i=0; i< validEntries; i++)
        {
            
            //get entry
            LinearNode node = linearMultiBuffer.nodes[cacheIndices[i]];
            
            vec3 d = position - node.location.xyz;
            float dist = length(d);

            float r = node.validRadius.x;
            //cubic kernel
            float u = clamp(1.0 - dist/r, 0.0, 1.0);
            last_u = u;
            float weight = 3.0 * u * u - 2.0 * u * u * u;  // 3u^2 - 2u^3
            //if(weight < 0.2) continue;
            //weight = pow(weight, 1.3);

            //clamp stored radiance so log never blows up
            vec3 L_cached = max(node.radiance.xyz, vec3(1e-4));
            vec3 logL = log(L_cached);

            float gradR = dot(node.gradientR.xyz, d);
            float gradG = dot(node.gradientG.xyz, d);
            float gradB = dot(node.gradientB.xyz, d);


            vec3 Lc = vec3(logL.r + gradientScalar * gradR / L_cached.r, logL.g + gradientScalar * gradG / L_cached.g, logL.b + gradientScalar * gradB / L_cached.b);
            //vec3 Lc = vec3(logL.r, logL.g, logL.b);

            for(int c=0; c<3; c++)
            {
                minLog = min(minLog, Lc[c]);
                maxLog = max(maxLog, Lc[c]);
            }
            
                            

            L+= weight * Lc;
            weightSum += weight;
            lastWeight = weight;
            
            
        }

        vec3 avgLog = (weightSum > 0.0) ? (L / weightSum) : vec3(0.0);

        avgLog = clamp(avgLog, minLog, maxLog);
       
        vec3 result = exp(avgLog);
       
        return result;
        
    }

//cloud density gradient function
vec3 sampleCloudGradient(vec3 pos) {
    //domain size considerations
    vec3 domainMin = parameters.boxMin;
    vec3 domainMax = parameters.boxMax;
    float densityCenter = sampleCloud(pos);
    vec3 grad;
    vec3 voxelSize = (parameters.boxMax - parameters.boxMin) / parameters.gridResolution;

    // Choose a small epsilon value for finite differences
    float relEps = 0.005;
    float epsX = voxelSize.x;
    float epsY = voxelSize.y;
    float epsZ = voxelSize.z;


    // Compute finite differences with clamping to avoid sampling outside the domain
    vec3 posXPlus = clamp(pos + vec3(epsX, 0.0, 0.0), domainMin, domainMax);
    vec3 posXMinus = clamp(pos - vec3(epsX, 0.0, 0.0), domainMin, domainMax);
    vec3 posYPlus = clamp(pos + vec3(0.0, epsY, 0.0), domainMin, domainMax);
    vec3 posYMinus = clamp(pos - vec3(0.0, epsY, 0.0), domainMin, domainMax);
    vec3 posZPlus = clamp(pos + vec3(0.0, 0.0, epsZ), domainMin, domainMax);
    vec3 posZMinus = clamp(pos - vec3(0.0, 0.0, epsZ), domainMin, domainMax);
    
    //----------------- X AXIS ---------------------------
    float densityXPlus  = sampleCloud(posXPlus);
    float densityXMinus = sampleCloud(posXMinus);

    if(densityXPlus > 0.001 && densityXMinus > 0.001)
    {
        //both sides valid
        grad.x = (densityXPlus - densityXMinus) / (2.0 * epsX);
    }
    else if(densityXPlus > 0.001)
    {
        //+ side valid
        grad.x = (densityXPlus - densityCenter) / (epsX);
    }
    else if(densityXMinus > 0.001)
    {
        //- side valid
        grad.x = (densityCenter - densityXMinus) / (epsX);
    }
    else
    {
        //none valid
        grad.x = 0.0;
    }

    //----------------- Y AXIS ---------------------------
    float densityYPlus  = sampleCloud(posYPlus);
    float densityYMinus = sampleCloud(posYMinus);

    if(densityYPlus > 0.001 && densityYMinus > 0.001)
    {
        //both sides valid
        grad.y = (densityYPlus - densityYMinus) / (2.0 * epsY);
    }
    else if(densityYPlus > 0.001)
    {
        //+ side valid
        grad.y = (densityYPlus - densityCenter) / (epsY);
    }
    else if(densityYMinus > 0.001)
    {
        //- side valid
        grad.y = (densityCenter - densityYMinus) / (epsY);
    }
    else
    {
        //none valid
        grad.y = 0.0;
    }

    //----------------- Z AXIS ---------------------------
    float densityZPlus  = sampleCloud(posZPlus);
    float densityZMinus = sampleCloud(posZMinus);

     if(densityZPlus > 0.001 && densityZMinus > 0.001)
    {
        //both sides valid
        grad.z = (densityZPlus - densityZMinus) / (2.0 * epsZ);
    }
    else if(densityZPlus > 0.001)
    {
        //+ side valid
        grad.z = (densityZPlus - densityCenter) / (epsZ);
    }
    else if(densityZMinus > 0.001)
    {
        //- side valid
        grad.z = (densityCenter - densityZMinus) / (epsZ);
    }
    else
    {
        //none valid
        grad.z = 0.0;
    }
    
    // Compute the central difference for each axis
    
    
    
    
    
    return grad * voxelSize;
}

//Transmittance Calculation Functions

float opticalThicknessCalculation(vec3 x, vec3 x_i, out vec3 gradient, float offsetScalar)
{
    //grab parameters
    float majorant = parameters.extinction.x;
    float absorptionAlbedo = 1.0 - parameters.scatteringAlbedo.x;
    float scatteringAlbedo = parameters.scatteringAlbedo.x;
    float PA = absorptionAlbedo * parameters.extinction.x;
    float PS = scatteringAlbedo * parameters.extinction.x;
    float densityToExtinction = 3.0;
    float PA2 = absorptionAlbedo * densityToExtinction;
    float PS2 = scatteringAlbedo * densityToExtinction;



   //precompute
    vec3 dir = x_i - x;
    float L = length(dir);
    if(L < 1e-6)
    {
        gradient = vec3(0); return 0;
    }
    vec3 dir_n = dir / L;

    //marching setup
    float stepSize = 0.01;
    int N = max(1, int(L/stepSize));
    float invN = 1.0 / float(N);
    float dt = L / float(N);
    float offset = offsetScalar * dt;

    //integrate

    float tau = 0.0;
    gradient = vec3(0.0);

    for(int i=0; i<N; i++)
    {
        //get sampled position
        float tMid = offset + (float(i)+0.5) * dt;

        vec3 pos = x * dir_n * tMid;

        //get density at that positon
        float density = sampleCloud(pos);
        
        //get density gradient at that position using finite differences
        vec3 densityGradient = sampleCloudGradient(pos);
       

        //based on that density, calculate sigma_t (extinction coefficient)
        
        float sigma_t = (PA + PS) * density; // / majorant;

        //based on densityGradient, get delta_sigma_t
        vec3 delta_sigma_t = (PA + PS) * densityGradient;

      

        //accumulate
        tau += sigma_t * dt;
        gradient += delta_sigma_t * dt;
       
       


    }
    return tau;
}


float calculateTransmittance(vec3 x, vec3 w, out vec3 gradient) {
    float majorant = parameters.extinction.x;
    float absorptionAlbedo = 1.0 - parameters.scatteringAlbedo.x;
    float scatteringAlbedo = parameters.scatteringAlbedo.x;
    float PA = absorptionAlbedo * parameters.extinction.x;
    float PS = scatteringAlbedo * parameters.extinction.x;

    float transmittance = 1.0;
    float rr_factor = 1.0;
    int stepCounter = 0;

    gradient = vec3(0.0);
    ivec2 localID = ivec2(gl_GlobalInvocationID.xy);
    ivec2 imageCoord = pc.tileOffset + localID;

    float tMin, tMax;
    if (rayBoxIntersect(parameters.boxMin, parameters.boxMax, x, w, tMin, tMax)) {
        x += w * tMin;
        float d = tMax - tMin;

        float targetThroughput = 0.1;

        while (true) {
            stepCounter++;
            float t = -log(max(0.0000000001, 1 - random())) / majorant;

            if (t > d) {
                break;
            }


            vec3 xNew = x + w * t;


            float density = sampleCloud(xNew);
            vec3 densityGradient = sampleCloudGradient(xNew);
            if(density == 0.0) break; //assume to more cloud




            x = xNew;

            float sigma_a = PA * density;
            float sigma_s = PS * density;
            float sigma_n = majorant - parameters.extinction.x * density;

            float Pa = sigma_a / majorant;
            float Ps = sigma_s / majorant;
            float Pn = sigma_n / majorant;

            vec3 delta_sigma_t = (PA + PS) * densityGradient;

            if(density > 1e-4)
            {
                transmittance *= Pn;
            }
            
            gradient += delta_sigma_t * (tMax - tMin) / majorant;
            if(transmittance < 1e-5) break;
           
            

            d -= t;
        }
    }
    return transmittance;
}

float calculateTransmittanceDistance(vec3 x, vec3 w, float maxDist) {
    float majorant = parameters.extinction.x;
    float absorptionAlbedo = 1.0 - parameters.scatteringAlbedo.x;
    float scatteringAlbedo = parameters.scatteringAlbedo.x;
    float PA = absorptionAlbedo * parameters.extinction.x;
    float PS = scatteringAlbedo * parameters.extinction.x;

    float transmittance = 1.0;
    float rr_factor = 1.0;

    float tMin, tMax;
    float currDist = 0.0;
    if (rayBoxIntersect(parameters.boxMin, parameters.boxMax, x, w, tMin, tMax)) {
        x += w * tMin;
        currDist += tMin;
        float d = tMax - tMin;

        float targetThroughput = 0.1;

        while (true) {
            float t = -log(max(0.0000000001, 1 - random())) / majorant;
            currDist += t;

            if (t > d || currDist > maxDist) {
                break;
            }

           
            vec3 xNew = x + w * t;
            float density = sampleCloud(xNew);



            x = xNew;

            float sigma_a = PA * density;
            float sigma_s = PS * density;
            float sigma_n = majorant - parameters.extinction.x * density;

            float Pa = sigma_a / majorant;
            float Ps = sigma_s / majorant;
            float Pn = sigma_n / majorant;


            transmittance *= Pn;

            d -= t;
        }
    }
    return transmittance * rr_factor;
}


vec3 computeSingleScattering_Full(vec3 x, vec3 w, int N, out vec3 gradientR, out vec3 gradientG, out vec3 gradientB, out float radius)
    {
        //parameters
       
        float absorptionAlbedo = 1.0 - parameters.scatteringAlbedo.x;
        float scatteringAlbedo = parameters.scatteringAlbedo.x;
       
        float invN = 1.0/ float(N);

        //zero out accumulators for safety
        vec3 radiance = vec3(0.0);
        gradientR = vec3(0);
        gradientG = vec3(0);
        gradientB = vec3(0);
        radius = 0.0;

        
        float pdf_w = 0;

        float density = sampleCloud(x);
        float sigma_s = scatteringAlbedo * density * 1.0;
       

        vec3 radianceSums = vec3(0.0);
        vec3 gradientMagnitudes = vec3(0.0);
        float scalingFactor = 1;

        float offsetScalar = max(random(), 1e-8);

        float epsilon = 1.0;

        for(int i=0; i<N; i++)
        {
            //sample one incoming direction w' and its pdf
            vec3 wi = importanceSamplePhase(parameters.phaseG, w, pdf_w);
            
            //phase function value p (w' -> w)
            float p = evaluatePhase(parameters.phaseG, w, wi);


            //find how far to march to exit the medium
            float tMin, tMax;
            rayBoxIntersect(parameters.boxMin, parameters.boxMax, x, wi, tMin, tMax);

            //march to tMax

            //compute optical thickness
            vec3 delta_Tr = vec3(0.0);
            float Tr = calculateTransmittance(x, wi, delta_Tr);

            //sample light at the exit points
            vec3 L_env = sampleSkybox(wi) + sampleLight(wi);

            //accumulate weighted contribution
            //vec3 L = (Tr * L_env) / pdf_w;
            vec3 L = (p * Tr * L_env) / pdf_w;
            radiance += L;
            radianceSums += L;

            //delta p
            vec3 delta_p = evaluatePhaseGradient(parameters.phaseG, w, wi);

            //delta L
            
            //vec3 delta_Tr = -delta_tau * Tr;
            vec3 delta_L = L * delta_Tr;


            vec3 gradR = (p * L_env.r / pdf_w) * delta_Tr;
            vec3 gradG = (p * L_env.g / pdf_w) * delta_Tr;
            vec3 gradB = (p * L_env.b / pdf_w) * delta_Tr;

            gradientR += gradR;
            gradientG += gradG;
            gradientB += gradB;

            
            gradientMagnitudes.r += length(gradR);
            gradientMagnitudes.g += length(gradG);
            gradientMagnitudes.b += length(gradB);
           
            //insertResult(x, wi, vec3(p, pdf_w, Tr), vec4(calculatedGradient, 0));

                
               
                
               
        }

        //normalizing things
        radiance *= invN;
        gradientR *= invN;
        gradientG *= invN;
        gradientB *= invN;
        gradientMagnitudes *= invN;

        vec3 radiusVector = vec3(0.0);
        float radiusUnclamped = 0.0;
        if(gradientMagnitudes.x < 1e-6 || gradientMagnitudes.y < 1e-6 || gradientMagnitudes.z < 1e-6) 
        {
            radius= 0.0;
        }
        else
        {
             radiusVector = epsilon * scalingFactor * (radianceSums / gradientMagnitudes);
             radiusUnclamped = min(min(radiusVector.x, radiusVector.y), radiusVector.z);

             radius = clamp(radiusUnclamped, MIN_RADIUS, MAX_RADIUS);
        }
       

        //insertResult(radiance, gradientR, gradientMagnitudes, vec4(radiusUnclamped, 0, 0 , radius));
       
       
       return radiance;
    }


vec3 computeMultiScattering_Full(vec3 x, vec3 w, int N, out vec3 gradientR, out vec3 gradientG, out vec3 gradientB, out float radius)
{
     //parameters
        float majorant = parameters.extinction.x;
        float absorptionAlbedo = 1.0 - parameters.scatteringAlbedo.x;
        float scatteringAlbedo = parameters.scatteringAlbedo.x;
        float PA = absorptionAlbedo * majorant;
        float PS = scatteringAlbedo * majorant;
       
        //scattering coefficient
        float density = sampleCloud(x);
        vec3 densityGradient = sampleCloudGradient(x);
        float sigma_s = PS * density;
        float sigma_t = (PA+PS) * density;
        vec3 delta_sigma_t = (PA+PS) * densityGradient;
       
        float invN = 1.0/ float(N);

        //zero out accumulators for safety
        vec3 radiance = vec3(0.0);
        gradientR = vec3(0);
        gradientG = vec3(0);
        gradientB = vec3(0);
        radius = 0.0;

        vec3 radianceSums= vec3(0);
        vec3 gradientMagnitudes = vec3(0);
        float tMin;
        float tMax;
        int maxBounces = 20;

        float epsilon = 1.0;
        float scalingFactor = 1.0;


        for(int i=0; i < N; i++)
        {
            //first bounce: sample direction and distance
            float pdf_w0;
            vec3 w0 = importanceSamplePhase(parameters.phaseG, w, pdf_w0);

            if(!rayBoxIntersect(parameters.boxMin, parameters.boxMax, x, w0, tMin, tMax))
            {
                continue; //if no entry at all
            }

            float r0 = -log(max(0.000000001, 1.0 - random())) / sigma_t;
            float pdf_r0 = sigma_t * exp(-sigma_t * r0);

            //clamp to box exit
            float travel0 = min(r0, tMax);
            vec3 p0 = x + w0 * travel0;

            //initial throughput
            float offsetScalar = max(random(), 1e-6);
            vec3 delta_tau0;
            float tau0 = opticalThicknessCalculation(x, p0, delta_tau0, offsetScalar);
            float Tr0 = exp(-tau0);
            float p_val = evaluatePhase(parameters.phaseG, w, w0);
            float throughput = (sigma_s * p_val * Tr0) / (pdf_w0 * pdf_r0);
            throughput = clamp(throughput, 0.0, 1.0);
            //insertResult(p0, vec3(throughput, sigma_t, density), vec3(Tr0, r0, tMax), vec4(0.0));


            vec3 L_path = vec3(0.0);

            //if we already exceeded tmax, exit
            if(r0 >= tMax)
            {
                vec3 L_env = (sampleSkybox(w0) + sampleLight(w0));
                L_path = throughput * L_env;
                radiance += L_path;

                vec3 delta_Tr = -delta_tau0 * Tr0;
                vec3 delta_sigma_s = PS * densityGradient;
                vec3 delta_L_envR = delta_Tr * L_env.r;
                vec3 delta_L_envG = delta_Tr * L_env.g;
                vec3 delta_L_envB = delta_Tr * L_env.b;

                vec3 gradR = (p_val * (delta_Tr * sigma_s * L_env.r + Tr0 * delta_sigma_s * L_env.r + Tr0 * sigma_s * delta_L_envR)) / (pdf_w0 * pdf_r0); 
                vec3 gradG = (p_val * (delta_Tr * sigma_s * L_env.g + Tr0 * delta_sigma_s * L_env.g + Tr0 * sigma_s * delta_L_envG)) / (pdf_w0 * pdf_r0);
                vec3 gradB = (p_val * (delta_Tr * sigma_s * L_env.b + Tr0 * delta_sigma_s * L_env.b + Tr0 * sigma_s * delta_L_envB)) / (pdf_w0 * pdf_r0);
                gradientR += gradR;
                gradientG += gradG;
                gradientB += gradB;

                radianceSums += L_path;
                gradientMagnitudes.r += length(gradR);
                gradientMagnitudes.g += length(gradG);
                gradientMagnitudes.b += length(gradB);
                //insertResult(p0, vec3(throughput, sigma_t, density), radiance, vec4(0.0));
                continue;
            }

            //iterative random walk for subsequent bounces
            vec3 currPos = p0;
            vec3 currDir = w0;

            for(int b =0; b< maxBounces; b++)
            {
                //sample new direction
                float pdf_wi;
                vec3 wi = importanceSamplePhase(parameters.phaseG, currDir, pdf_wi);

                //sample distance within medium
                density = sampleCloud(currPos);
                densityGradient = sampleCloudGradient(currPos);

                //intersect box to get tMax
                if(!rayBoxIntersect(parameters.boxMin, parameters.boxMax, currPos, wi, tMin, tMax) || density < 0.01)
                {
                    //no further medium / this should never trigger
                    vec3 L_env = (sampleSkybox(wi) + sampleLight(wi));
                    L_path += throughput * L_env;

                    vec3 delta_Tr = -delta_tau0 * Tr0;
                    vec3 delta_sigma_s = PS * densityGradient;
                    vec3 delta_L_envR = delta_Tr * L_env.r;
                    vec3 delta_L_envG = delta_Tr * L_env.g;
                    vec3 delta_L_envB = delta_Tr * L_env.b;

                    vec3 gradR = (p_val * (delta_Tr * sigma_s * L_env.r + Tr0 * delta_sigma_s * L_env.r + Tr0 * sigma_s * delta_L_envR)) / (pdf_w0 * pdf_r0); 
                    vec3 gradG = (p_val * (delta_Tr * sigma_s * L_env.g + Tr0 * delta_sigma_s * L_env.g + Tr0 * sigma_s * delta_L_envG)) / (pdf_w0 * pdf_r0);
                    vec3 gradB = (p_val * (delta_Tr * sigma_s * L_env.b + Tr0 * delta_sigma_s * L_env.b + Tr0 * sigma_s * delta_L_envB)) / (pdf_w0 * pdf_r0);
                    gradientR += gradR;
                    gradientG += gradG;
                    gradientB += gradB;

                    radianceSums += L_path;
                    gradientMagnitudes.r += length(gradR);
                    gradientMagnitudes.g += length(gradG);
                    gradientMagnitudes.b += length(gradB);
                    //insertResult(currPos, vec3(throughput, sigma_t, density), radiance, vec4(1.0));
                    break;
                }

                
                float sigma_t_i = (PA + PS) * density;
                float sigma_s_i = PS * density;
                float r_i = -log(max(0.000000001, 1.0 - random())) / sigma_t_i;
                float pdf_ri = sigma_t_i * exp(-sigma_t_i * r_i);
                float travel_i = min(r_i, tMax);

                vec3 nextPos = currPos + wi * travel_i;

                //update throughput
                offsetScalar = max(random(), 1e-6);
                vec3 delta_tau_i;
                float tau_i = opticalThicknessCalculation(currPos, nextPos, delta_tau_i, offsetScalar);
                float Tr_i = exp(-tau_i);
                float p_i = evaluatePhase(parameters.phaseG, currDir, wi);
                float multiplicator = (sigma_s_i/sigma_t_i * p_i * Tr_i) / (pdf_wi * pdf_ri);
                //if(multiplicator > 1.0) insertResult(currPos, vec3(multiplicator, sigma_s_i, sigma_t_i), vec3(p_i, Tr_i, 0), vec4(pdf_wi, pdf_ri, 0, 0));
                multiplicator = clamp(multiplicator, 0.0, 1.0);

                throughput *= multiplicator;

                //if we escape medium, accumulate and break
                if(r_i >= tMax)
                {
                    vec3 L_env = (sampleSkybox(wi) + sampleLight(wi));
                    L_path += throughput * L_env;

                    vec3 delta_Tr = -delta_tau_i * Tr_i;
                    vec3 delta_sigma_s = PS * densityGradient;
                    vec3 delta_L_envR = delta_Tr * L_env.r;
                    vec3 delta_L_envG = delta_Tr * L_env.g;
                    vec3 delta_L_envB = delta_Tr * L_env.b;

                    vec3 gradR = (p_val * (delta_Tr * sigma_s_i * L_env.r + Tr_i * delta_sigma_s * L_env.r + Tr_i * sigma_s_i * delta_L_envR)) / (pdf_wi * pdf_ri); 
                    vec3 gradG = (p_val * (delta_Tr * sigma_s_i * L_env.g + Tr_i * delta_sigma_s * L_env.g + Tr_i * sigma_s_i * delta_L_envG)) / (pdf_wi * pdf_ri);
                    vec3 gradB = (p_val * (delta_Tr * sigma_s_i * L_env.b + Tr_i * delta_sigma_s * L_env.b + Tr_i * sigma_s_i * delta_L_envB)) / (pdf_wi * pdf_ri);
                    gradientR += gradR;
                    gradientG += gradG;
                    gradientB += gradB;

                    radianceSums += L_path;
                    gradientMagnitudes.r += length(gradR);
                    gradientMagnitudes.g += length(gradG);
                    gradientMagnitudes.b += length(gradB);


                    //insertResult(nextPos, vec3(throughput, sigma_t_i, density), radiance, vec4(3.0));

                    break;
                }

                //otherwise continue deeper
                currPos = nextPos;
                currDir = wi;

            }

            radiance += L_path;
        }


        //normalizing things
        radiance *= invN;
        gradientR *= invN;
        gradientG *= invN;
        gradientB *= invN;
        gradientMagnitudes *= invN;

        vec3 radiusVector = vec3(0.0);
        float radiusUnclamped = 0.0;
        if(gradientMagnitudes.x < 1e-6 || gradientMagnitudes.y < 1e-6 || gradientMagnitudes.z < 1e-6) 
        {
            radius= 0.0;
        }
        else
        {
             radiusVector = epsilon * scalingFactor * (radianceSums / gradientMagnitudes);
             radiusUnclamped = min(min(radiusVector.x, radiusVector.y), radiusVector.z);

             radius = clamp(radiusUnclamped, MIN_RADIUS, MAX_RADIUS);
        }
        return radiance;
}

vec3 deltaTrackingForMS(vec3 x, vec3 w, out vec3 w_exit) {

    vec3 radiance = vec3(0.0);
    float majorant = parameters.extinction.x;
    float absorptionAlbedo = 1.0 - parameters.scatteringAlbedo.x;
    float scatteringAlbedo = parameters.scatteringAlbedo.x;
    float PA = absorptionAlbedo * parameters.extinction.x;
    float PS = scatteringAlbedo * parameters.extinction.x;

    int i = 0;
    float tMin, tMax;
    float throughput = 1.0;
    if (rayBoxIntersect(parameters.boxMin, parameters.boxMax, x, w, tMin, tMax)) {
        x += w * tMin;

        float d = tMax - tMin;

        float pdf_x = 1;
        float transmittance = 1.0;

        while (true) {

            float t = -log(max(0.0000000001, 1 - random()))/majorant;

            if (t > d) {
                break;
            }

            if(throughput < 1e-4)
            {
                return vec3(0.0);
            }
            vec3 xNew = x + w * t;

            float density = sampleCloud(xNew);

            x = xNew;
            transmittance *= 1.0 - density;

            float sigma_a = PA * density;
            float sigma_s = PS * density;
            float sigma_n = majorant - parameters.extinction.x * density;

            float Pa = sigma_a / majorant;
            float Ps = sigma_s / majorant;
            float Pn = sigma_n / majorant;

            float xi = random();

            if (xi < Pa) {

                return vec3(0.0);// weights * sigma_a / (majorant * Pa) * L_e; // 0 - No emission
            }

            if (xi < 1 - Pn) // scattering event
            {
                float pdf_w;
                vec3 w_s = importanceSamplePhase(parameters.phaseG, w, pdf_w);

                float p = evaluatePhase(parameters.phaseG, w, w_s);
                w = w_s;
                w_exit = w;
                
                if (rayBoxIntersect(parameters.boxMin, parameters.boxMax, x, w, tMin, tMax)) {
                    x += w*tMin;
                    d = tMax - tMin;
                }
                throughput *= Ps * (p / pdf_w);
                if(throughput < 1e-4)
                {
                    return vec3(0.0);
                }
            } else {
                
                d -= t;               
            }
        }
    }


    //if(gl_GlobalInvocationID.xy == ivec2(127, 42)) insertResult(x, radiance, gl_GlobalInvocationID, vec4(0));
    return throughput * (sampleSkybox(w) + sampleLight(w));

}

vec3 MultiSimple(vec3 x, vec3 w, int N)
{
    float majorant = parameters.extinction.x;
    float absorptionAlbedo = 1.0 - parameters.scatteringAlbedo.x;
    float scatteringAlbedo = parameters.scatteringAlbedo.x;
    float PA = absorptionAlbedo * majorant;
    float PS = scatteringAlbedo * majorant;

    float tMin, tMax;
    
    float offsetScalar = max(random(), 1e-8);
    
    //init accumulators
    vec3 radiance = vec3(0.0);
    


    float invN = 1.0 / float(N);
    for(int i=0; i< N; i++)
    {
        //sample new direction

        float pdf_w;
        vec3 wi = importanceSamplePhase(parameters.phaseG, w, pdf_w);
        float p = evaluatePhase(parameters.phaseG, w, wi);

        //spherical harmonics projection
       
        vec3 w_exit = vec3(0.0);
        vec3 randomPathRadiance = deltaTrackingForMS(x, wi, w_exit);
        radiance += randomPathRadiance * (p / pdf_w);

    }

    
    radiance *= invN;
   

    return radiance;


}




void evalSHDebug(in vec3 dir, out float Y[9])
{
     //l = 0, m = 0
     Y[0] = 1.0;           // 1/2 * sqrt(1/pi)

     for(int i=1; i<9; i++)
     {
        Y[i] = 0.0;
     }
}


vec3 SingleSH(vec3 x, vec3 w, int N, out float coefficients[27], out vec3 gradients[27], out float radius)
{
    //g is 3 vec3s one for each color channel for the first coefficient
    //parameters 
    float majorant = parameters.extinction.x; // majorant is absorption + scattering + null coefficients
    float absorptionAlbedo = 1.0 - parameters.scatteringAlbedo.x; // absorption probability at scattering event
    float scatteringAlbedo = parameters.scatteringAlbedo.x; //scattering albedo [0,1] describes the probability of scattering (verus absorption) at a scattering evenet
    float PA = absorptionAlbedo * parameters.extinction.x;
    float PS = scatteringAlbedo * parameters.extinction.x;
    float density = sampleCloud(x);
    float sigma_s = PS * density;
    //float sigma_t = density * 3.0;
    //float sigma_s = scatteringAlbedo * sigma_t;
       
    float invN = 1.0 / float(N);
    float tMin, tMax;
    vec3 acc = vec3(0.0);

    vec3 radianceSums = vec3(0.0);
    vec3 gradientMagnitudes = vec3(0.0);

    float epsilon = 1;
    float scalingFactor = 1;

   


    float offsetScalar = max(random(), 1e-8);
    
    //init accumulators
    for(int k=0; k< 27; ++k)
    {
        coefficients[k] = 0.0;
        gradients[k] = vec3(0.0);
    }

    for(int i=0; i< N; ++i)
    {
        //sample new direction
        float pdf_w;
        vec3 wi = importanceSamplePhase(parameters.phaseG, w, pdf_w);
        float cosTH = dot(w, wi);

        float p = evaluatePhase(parameters.phaseG, w, wi);

        //compute tau and delta tau via ray marching
        float Tr = 0.0;
        vec3 delta_Tr = vec3(0.0);
        //get tMax where x + tMax * wi is the edge of the medium
        rayBoxIntersect(parameters.boxMin, parameters.boxMax, x, wi, tMin, tMax);

        Tr = calculateTransmittance(x, wi, delta_Tr);


        //sample the environment
        vec3 Le = sampleSkybox(wi) + sampleLight(wi);
        radianceSums += Le;
        acc += Tr * Le;

        //weight factors
        //sigma_s = min(sigma_s, sigma_s_max);
        float weightR = (p * Le.r * Tr) / pdf_w;
        float weightG = (p * Le.g * Tr) / pdf_w;
        float weightB = (p * Le.b * Tr) / pdf_w;

        vec3 delta_Le = delta_Tr * Le;
        vec3 delta_p = evaluatePhaseGradient(parameters.phaseG, w, wi);

        vec3 delta_weightR = (delta_p * Le.r + p * Le.r * Tr) / pdf_w;
        vec3 delta_weightG = (delta_p * Le.g + p * Le.g * Tr) / pdf_w; 
        vec3 delta_weightB = (delta_p * Le.b + p * Le.b * Tr) / pdf_w;
        
        gradientMagnitudes += length(delta_weightR) + length(delta_weightG) + length(delta_weightB);

        //project into sh
        float Y[9];
        vec3 dir = normalize(wi);
        evalSH2(dir, Y);
        


        for(int k=0; k < 9; ++k)
        {
            coefficients[k] += weightR * Y[k];
            coefficients[9 + k] += weightG * Y[k];
            coefficients[18+ k] += weightB * Y[k];

            gradients[k] += delta_weightR * Y[k];
            gradients[9 + k] += delta_weightG * Y[k];
            gradients[18 + k] += delta_weightB * Y[k];
        }
        
       

        //debug
        float ratio = p / pdf_w;
        
       
        
        //maxWeight = max(maxWeight, length(weight));
        //insertResult(x, vec3(p_debug, pdf_w, debugG), vec3(0.0), vec4(0));
        ivec2 localID = ivec2(gl_GlobalInvocationID.xy);
        
        ivec2 imageCoord = pc.tileOffset + localID;
        //if(imageCoord == ivec2(720, 365)) insertResult(x, vec3(p_debug), vec3(pdf_w), vec4(Tr, Le.r, coeffR[0], i));
        
    }

    //normalize
    for(int k=0; k < 27; ++k)
    {
        coefficients[k] *= invN;
        gradients[k] *= invN;
    }
   
    
    //evalute radiance
    vec3 radiance = vec3(0.0);
    float Y[9];
    vec3 outDir = normalize(w);
    evalSH2(outDir, Y);
    for(int k=0; k < 9; ++k)
    {
        radiance.r += coefficients[k] * Y[k];
        radiance.g += coefficients[9 + k] * Y[k];
        radiance.b += coefficients[18 + k] * Y[k];
    }

    vec3 radiusVector = vec3(0.0);
    float radiusUnclamped = 0.0;
    if(gradientMagnitudes.x < 1e-6 || gradientMagnitudes.y < 1e-6 || gradientMagnitudes.z < 1e-6) 
    {
        radius= 0.0;
    }
    else
    {
         radiusVector = epsilon * scalingFactor * (radianceSums / gradientMagnitudes);
         radiusUnclamped = min(min(radiusVector.x, radiusVector.y), radiusVector.z);
         radius = clamp(radiusUnclamped, MIN_RADIUS, MAX_RADIUS);
    }
    //float test = sigma_s_max / N * Y[0];
    //insertResult(x, radiance, acc, vec4(0));
    //if(radiance.r > 0.1 || radiance.g > 0.1 || radiance.b > 0.1) insertResult(vec3(coeffR[0], coeffG[0], coeffB[0]), g[0], g[1], vec4(g[2], radius));
    vec3 clampedRadiance = clamp(radiance, 0.0, 1.0);
    return clampedRadiance;

}

vec3 MultiSH(vec3 x, vec3 w, int N, out float coeffR[9], out float coeffG[9], out float coeffB[9], out vec3 g[9])
{
    //parameters
        float majorant = parameters.extinction.x;
        float absorptionAlbedo = 1.0 - parameters.scatteringAlbedo.x;
        float scatteringAlbedo = parameters.scatteringAlbedo.x;
        float PA = absorptionAlbedo * majorant;
        float PS = scatteringAlbedo * majorant;
       
        //scattering coefficient
        float density = sampleCloud(x);
        vec3 densityGradient = sampleCloudGradient(x);
        float sigma_s = PS * density;
        float sigma_t = (PA+PS) * density;
        vec3 delta_sigma_t = (PA+PS) * densityGradient;

       
        float invN = 1.0/ float(N);

        



        //replace with output variable
        float radius = 0.0;
        //init accumulators
        for(int k=0; k< 9; ++k)
        {
            coeffR[k] = 0.0;
            coeffG[k] = 0.0;
            coeffB[k] = 0.0;
            g[k] = vec3(0.0);
        }

        vec3 radianceSums= vec3(0);
        vec3 gradientMagnitudes = vec3(0);
        float tMin;
        float tMax;
        int maxBounces = 20;

        float epsilon = 1.0;
        float scalingFactor = 1.0;
        int n=0;

        for(int i=0; i < N; i++)
        {
            n++;
            //first bounce: sample direction and distance
            float pdf_w0;
            vec3 w0 = importanceSamplePhase(parameters.phaseG, w, pdf_w0);

            //get corresponding spherical harmonics part for first bounce direction
            float Y0[9];
            vec3 dir = normalize(w0);
            evalSH2(dir, Y0);

            if(!rayBoxIntersect(parameters.boxMin, parameters.boxMax, x, w0, tMin, tMax))
            {
                continue; //if no entry at all
            }

            float r0 = -log(max(0.0000000001, 1 - random())) / majorant;
            float pdf_r0 = majorant * exp(-majorant * r0);

            //clamp to box exit
            float travel0 = min(r0, tMax);
            vec3 p0 = x + w0 * travel0;

            //initial throughput
            float offsetScalar = max(random(), 1e-6);
            vec3 delta_tau0;
            //float tau0 = opticalThicknessCalculation(x, p0, delta_tau0, offsetScalar);
            //float Tr0 = exp(-tau0);
            float Tr0 = calculateTransmittanceDistance(x, w0, travel0);
            float p_val = evaluatePhase(parameters.phaseG, w, w0);
            float throughput = (sigma_s * p_val * Tr0) / (pdf_w0 * pdf_r0);
            //throughput = clamp(throughput, 0.0, 1.0);
            ivec2 localID = ivec2(gl_GlobalInvocationID.xy);
        
            ivec2 imageCoord = pc.tileOffset + localID;
            //if(imageCoord == ivec2(720, 365)) insertResult(p0, vec3(throughput, sigma_t, density), vec3(Tr0, r0, tMax), vec4(0, 0, 0, 0));
                


            vec3 L_path = vec3(0.0);

            //if we already exceeded tmax, exit
            if(r0 >= tMax)
            {
                //we have zero more bounces after the first -> skip


                /*
                vec3 L_env = (sampleSkybox(w0) + sampleLight(w0));
                L_path =throughput * L_env;
                //project L_path to spherical harmonics
                for(int k=0; k< 9; ++k)
                {
                    coeffR[k] += L_path.r * Y0[k];
                    coeffG[k] += L_path.g * Y0[k];
                    coeffB[k] += L_path.b * Y0[k];
                }

                vec3 delta_Tr = -delta_tau0 * Tr0;
                vec3 delta_sigma_s = PS * densityGradient;
                vec3 delta_L_envR = delta_Tr * L_env.r;
                vec3 delta_L_envG = delta_Tr * L_env.g;
                vec3 delta_L_envB = delta_Tr * L_env.b;

                vec3 gradR = (p_val * (delta_Tr * sigma_s * L_env.r + Tr0 * delta_sigma_s * L_env.r + Tr0 * sigma_s * delta_L_envR)) / (pdf_w0 * pdf_r0); 
                vec3 gradG = (p_val * (delta_Tr * sigma_s * L_env.g + Tr0 * delta_sigma_s * L_env.g + Tr0 * sigma_s * delta_L_envG)) / (pdf_w0 * pdf_r0);
                vec3 gradB = (p_val * (delta_Tr * sigma_s * L_env.b + Tr0 * delta_sigma_s * L_env.b + Tr0 * sigma_s * delta_L_envB)) / (pdf_w0 * pdf_r0);
                

                radianceSums += L_path;
                gradientMagnitudes.r += length(gradR);
                gradientMagnitudes.g += length(gradG);
                gradientMagnitudes.b += length(gradB);
                //insertResult(p0, vec3(throughput, sigma_t, density), radiance, vec4(0.0));
                
                */
                continue;
            }

            //iterative random walk for subsequent bounces
            vec3 currPos = p0;
            vec3 currDir = w0;

            for(int b =0; b< maxBounces; b++)
            {
                //sample new direction
                float pdf_wi;
                vec3 wi = importanceSamplePhase(parameters.phaseG, currDir, pdf_wi);

                //get corresponding spherical harmonics part for next bounce direction
                float Yi[9];
                dir = normalize(wi);
                evalSH2(dir, Yi);

                //sample distance within medium
                density = sampleCloud(currPos);
                densityGradient = sampleCloudGradient(currPos);
                if(density < 0.01) continue;

                //intersect box to get tMax
                if(!rayBoxIntersect(parameters.boxMin, parameters.boxMax, currPos, wi, tMin, tMax) || density < 0.01)
                {
                    
                    //no further medium / this should never trigger
                    vec3 L_env = (sampleSkybox(wi) + sampleLight(wi));
                    L_path = throughput * L_env;

                    for(int k=0; k< 9; ++k)
                    {
                        coeffR[k] += L_path.r * Yi[k];
                        coeffG[k] += L_path.g * Yi[k];
                        coeffB[k] += L_path.b * Yi[k];
                    }

                    vec3 delta_Tr = -delta_tau0 * Tr0;
                    vec3 delta_sigma_s = PS * densityGradient;
                    vec3 delta_L_envR = delta_Tr * L_env.r;
                    vec3 delta_L_envG = delta_Tr * L_env.g;
                    vec3 delta_L_envB = delta_Tr * L_env.b;

                    vec3 gradR = (p_val * (delta_Tr * sigma_s * L_env.r + Tr0 * delta_sigma_s * L_env.r + Tr0 * sigma_s * delta_L_envR)) / (pdf_w0 * pdf_r0); 
                    vec3 gradG = (p_val * (delta_Tr * sigma_s * L_env.g + Tr0 * delta_sigma_s * L_env.g + Tr0 * sigma_s * delta_L_envG)) / (pdf_w0 * pdf_r0);
                    vec3 gradB = (p_val * (delta_Tr * sigma_s * L_env.b + Tr0 * delta_sigma_s * L_env.b + Tr0 * sigma_s * delta_L_envB)) / (pdf_w0 * pdf_r0);
                    

                    radianceSums += L_path;
                    gradientMagnitudes.r += length(gradR);
                    gradientMagnitudes.g += length(gradG);
                    gradientMagnitudes.b += length(gradB);
                    //if(imageCoord == ivec2(720, 365)) insertResult(currPos, vec3(throughput, sigma_t, density), vec3(coeffR[0], coeffG[0], coeffB[0]), vec4(1.0));
                    break;
                    
                }

                
                float sigma_t_i = (PA + PS) * density;
                float sigma_s_i = PS * density;
                float r_i = -log(max(0.000000001, 1.0 - random())) / sigma_t_i;
                float pdf_ri = sigma_t_i * exp(-sigma_t_i * r_i);
                float travel_i = min(r_i, tMax);

               

                
                vec3 nextPos = currPos + wi * travel_i;
                //insertResult(currPos, nextPos, vec3(travel_i, tMin, tMax), vec4(travel_i, 0, 0, 0));

                //update throughput
                offsetScalar = max(random(), 1e-6);
                vec3 delta_tau_i;
                //float tau_i = opticalThicknessCalculation(currPos, nextPos, delta_tau_i, offsetScalar);
                //float Tr_i = exp(-tau_i);
                float Tr_i = calculateTransmittanceDistance(currPos, currDir, travel_i);
                float p_i = evaluatePhase(parameters.phaseG, currDir, wi);
                float multiplicator = (sigma_s_i * p_i * Tr_i) / (pdf_wi * pdf_ri);
                //insertResult(currPos, vec3(multiplicator, sigma_s_i, sigma_t_i), vec3(p_i, Tr_i, 0), vec4(pdf_wi, pdf_ri, 0, 0));
                //multiplicator = clamp(multiplicator, 0.0, 1.0);

                throughput *= multiplicator;

                //if we escape medium, accumulate and break
                if(r_i >= tMax)
                {
                    vec3 L_env = (sampleSkybox(wi) + sampleLight(wi));
                    L_path = throughput * L_env;

                    for(int k=0; k< 9; ++k)
                    {
                        coeffR[k] += L_path.r * Yi[k];
                        coeffG[k] += L_path.g * Yi[k];
                        coeffB[k] += L_path.b * Yi[k];
                    }

                    vec3 delta_Tr = -delta_tau_i * Tr_i;
                    vec3 delta_sigma_s = PS * densityGradient;
                    vec3 delta_L_envR = delta_Tr * L_env.r;
                    vec3 delta_L_envG = delta_Tr * L_env.g;
                    vec3 delta_L_envB = delta_Tr * L_env.b;

                    vec3 gradR = (p_val * (delta_Tr * sigma_s_i * L_env.r + Tr_i * delta_sigma_s * L_env.r + Tr_i * sigma_s_i * delta_L_envR)) / (pdf_wi * pdf_ri); 
                    vec3 gradG = (p_val * (delta_Tr * sigma_s_i * L_env.g + Tr_i * delta_sigma_s * L_env.g + Tr_i * sigma_s_i * delta_L_envG)) / (pdf_wi * pdf_ri);
                    vec3 gradB = (p_val * (delta_Tr * sigma_s_i * L_env.b + Tr_i * delta_sigma_s * L_env.b + Tr_i * sigma_s_i * delta_L_envB)) / (pdf_wi * pdf_ri);
                   

                    radianceSums += L_path;
                    gradientMagnitudes.r += length(gradR);
                    gradientMagnitudes.g += length(gradG);
                    gradientMagnitudes.b += length(gradB);

                    //if(imageCoord == ivec2(720, 365)) insertResult(nextPos, vec3(throughput, sigma_t_i, density), vec3(coeffR[0], coeffG[0], coeffB[0]), vec4(2.0));

                    break;
                }

                //otherwise continue deeper
                currPos = nextPos;
                currDir = wi;

            }

            //radiance += L_path;
        }
        


        //normalize
        for(int k=0; k < 9; ++k)
        {
            coeffR[k] *= invN;
            coeffG[k] *= invN;
            coeffB[k] *= invN;
            g[k] *= invN;
        }

        //reconstruct radiance
        vec3 radiance = vec3(0.0);
        float Y[9];
        vec3 outDir = normalize(w);
        evalSH2(outDir, Y);
        for(int k=0; k < 9; ++k)
        {
            radiance.r += coeffR[k] * Y[k];
            radiance.g += coeffG[k] * Y[k];
            radiance.b += coeffB[k] * Y[k];
        }
        vec3 clampedRadiance = clamp(radiance, 0.0, 1.0);
        


       

        vec3 radiusVector = vec3(0.0);
        float radiusUnclamped = 0.0;
        if(gradientMagnitudes.x < 1e-6 || gradientMagnitudes.y < 1e-6 || gradientMagnitudes.z < 1e-6) 
        {
            radius= 0.0;
        }
        else
        {
             radiusVector = epsilon * scalingFactor * (radianceSums / gradientMagnitudes);
             radiusUnclamped = min(min(radiusVector.x, radiusVector.y), radiusVector.z);

             radius = clamp(radiusUnclamped, MIN_RADIUS, MAX_RADIUS);
        }
        //insertResult(vec3(coeffR[0], coeffR[1], coeffR[2]), vec3(coeffG[0], coeffG[1], coeffG[2]), vec3(coeffB[0], coeffB[1], coeffB[2]), vec4(radiance, n));
        return clampedRadiance;
    
}




vec3 MultiSHSimple(vec3 x, vec3 w, int N, out float coeffR[9], out float coeffG[9], out float coeffB[9])
{
    float majorant = parameters.extinction.x;
    float absorptionAlbedo = 1.0 - parameters.scatteringAlbedo.x;
    float scatteringAlbedo = parameters.scatteringAlbedo.x;
    float PA = absorptionAlbedo * majorant;
    float PS = scatteringAlbedo * majorant;

    float tMin, tMax;
    
    float offsetScalar = max(random(), 1e-8);
    
    //init accumulators
    for(int k=0; k< 9; ++k)
    {
        coeffR[k] = 0.0;
        coeffG[k] = 0.0;
        coeffB[k] = 0.0;
    }


    float invN = 1.0 / float(N);
    for(int i=0; i< N; i++)
    {
        //sample new direction

        float pdf_w;
        vec3 wi = importanceSamplePhase(parameters.phaseG, w, pdf_w);

        //spherical harmonics projection
       
        vec3 w_exit = vec3(0.0);
        vec3 randomPathRadiance = deltaTrackingForMS(x, wi, w_exit);
        float Y[9];
        vec3 dir = normalize(w_exit);
        evalSH2(dir, Y);
        //float weight = p / pdf_w;
        for(int k=0; k< 9; ++k)
        {
            coeffR[k] += randomPathRadiance.r * Y[k];
            coeffG[k] += randomPathRadiance.g * Y[k];
            coeffB[k] += randomPathRadiance.b * Y[k];
        }

    }

     //normalize
    for(int k=0; k < 9; ++k)
    {
        coeffR[k] *= invN;
        coeffG[k] *= invN;
        coeffB[k] *= invN;
            
    }

    //spherical harmonics reprojection
    vec3 radiance = vec3(0);
    float Y[9];
    vec3 outDir = normalize(w);
    evalSH2(outDir, Y);
    for(int k=0; k < 9; ++k)
    {
        radiance.r += coeffR[k] * Y[k];
        radiance.g += coeffG[k] * Y[k];
        radiance.b += coeffB[k] * Y[k];
            
    }
    //insertResult(x, radiance, vec3(coeffR[0], coeffG[0], coeffB[0]), vec4(Y[0], Y[1], Y[3], 0));
    radiance = clamp(radiance, 0.0, 1.0);
    

    return radiance;

}


vec3 FirstEventCache(vec3 x, vec3 w, ScatterEvent firstEvent)
{
       
        float majorant = parameters.extinction.x;
        float absorptionAlbedo = 1.0 - parameters.scatteringAlbedo.x;
        float scatteringAlbedo = parameters.scatteringAlbedo.x;
        float PA = absorptionAlbedo * parameters.extinction.x;
        float PS = scatteringAlbedo * parameters.extinction.x;

        vec3 radiance = vec3(0.0);
        float throughput=1.0;

    


        int i = 0;
        float tMin, tMax;
        int contributions = 1;
        if (rayBoxIntersect(parameters.boxMin, parameters.boxMax, x, w, tMin, tMax)) {
            x += w * tMin;
   
            float d = tMax - tMin;

            float pdf_x = 1;
            float transmittance = 1.0;
            

            while (true) {
    
                //sample free flight distance
                float t = -log(max(0.0000000001, 1 - random())) / majorant;
                //if we exceed medium, exit the loop
                if (t > d) {
                    break;
                }
                //get new position
                vec3 xNew = x + w * t;
                
               
               
                x = xNew;

               

                float density = sampleCloud(xNew);

                float sigma_a = PA * density;
                float sigma_s = PS * density;
                float sigma_n = majorant - parameters.extinction.x * density;
                float sigma_t = sigma_a + sigma_s;

                float Pa = sigma_a / majorant;
                float Ps = sigma_s / majorant;
                float Pn = sigma_n / majorant;
                float Pt = sigma_t / majorant;

                float xi = random();
    

                if (xi < Pa) {
                   throughput = 0.0;
                    if (!firstEvent.hasValue)
                    {                   
                        firstEvent.hasValue = true;
                    }
                   break;
                }

                if (xi < 1 - Pn) // scattering event accumulate single and multi scattering, then continue along the ray
                {
                    if (!firstEvent.hasValue) 
                    {                   
                        firstEvent.hasValue = true;
                    }
                    contributions++;
                    vec3 gradient;
                    float radius;
                    vec3 singleScattering = vec3(0.0);
                    vec3 multiScattering = vec3(0.0);
                    vec3 neeContribution = vec3(0.0);
                    vec3 comparison = vec3(0.0);
                    vec3 singleResult=vec3(0.0);



                    #ifdef USE_CACHING
                        //for current position, sample both caches to find valid cache points

                        #ifdef USE_DIRECT_COMPONENT

                            uint snapShotDirect = atomicAdd(linearDirectCounter.nextFreeNode, 0);
                            uint directIndices[MAX_LINEAR_ENTRIES];
                            uint directCount = queryDirect(xNew, directIndices, snapShotDirect);
                            vec3 directGradientR = vec3(0.0);
                            vec3 directGradientG = vec3(0.0);
                            vec3 directGradientB = vec3(0.0);
                            float directRadius = 0.0;
                            float epsilon = 1.0;

                            if(directCount > 4)
                            {
                                 neeContribution = extrapolateDirectRadiance(directIndices, directCount, xNew);
                            }
                            else
                            {
                                float pdf_w;
                                vec3 lightDir = parameters.sunDirection;
                                vec3 L_nee = sampleSkybox(lightDir) + sampleLight(lightDir);
                                float neeP = evaluatePhase(parameters.phaseG, w, lightDir);
                                vec3 neeGradient;
                                float Tr_nee = calculateTransmittance(xNew, lightDir, neeGradient);
                                neeContribution = throughput * Ps * neeP * L_nee * Tr_nee;
                                //radiance += neeContribution;

                                //gradient. for isotropic, no delta p term
                                float wNee = throughput * Ps * neeP;

                                directGradientR = wNee * L_nee.r * neeGradient;
                                directGradientG = wNee * L_nee.g * neeGradient;
                                directGradientB = wNee * L_nee.b * neeGradient;
                                vec3 gradientMagnitudes = vec3(length(directGradientR), length(directGradientG), length(directGradientB));
                                vec3 radianceSums = neeContribution;
                                if(gradientMagnitudes.x < 1e-6 || gradientMagnitudes.y < 1e-6 || gradientMagnitudes.z < 1e-6) 
                                {
                                    directRadius= 0.0;
                                }
                                else
                                {
                                    vec3 radiusVector = epsilon * (radianceSums / gradientMagnitudes);
                                    float radiusUnclamped = min(min(radiusVector.x, radiusVector.y), radiusVector.z);

                                    directRadius = clamp(radiusUnclamped, MIN_RADIUS, MAX_RADIUS);
                                }

                                if(density > 1e-6 && random() < 0.1 && directRadius >= 0.0001) insertLinearDirect(xNew, neeContribution, directGradientR, directGradientG, directGradientB, directRadius);
                            }

                            

                   
                        #endif
                        

                        //single scattering component
                        #ifdef USE_SINGLESCATTERING
                            uint snapShotSingle = atomicAdd(linearCounter.nextFreeNode, 0);
                            uint singleIndices[MAX_LINEAR_ENTRIES];
                            uint singleCount = querySingleSH(xNew, singleIndices, snapShotSingle);
                            vec3 singleGradientR = vec3(0.0);
                            vec3 singleGradientG = vec3(0.0);
                            vec3 singleGradientB = vec3(0.0);
                            float singleRadius = 0.0;
                            

                            if(singleCount > 4)
                            {
                                //we have valid cache points  -> extrapolate

                                if(abs(parameters.phaseG) < 0.01)
                                {
                                    //isotropic case
                                    singleScattering = extrapolateSingleRadianceALT(singleIndices, singleCount, xNew);
                                }
                                else
                                {
                                    //anisotropic case
                                    singleScattering = extrapolateSingleSH_Version2(singleIndices, singleCount, xNew, -w);
                                }
                                //check this function again
                                
                                //comparison = computeSingleScattering_Full(xNew, w, 16, singleGradientR, singleGradientG, singleGradientB, singleRadius);
                                //insertResult(xNew, singleScattering, comparison, vec4(singleRadius, singleCount, 0, 0));
                            }
                            else
                            {
                                //we have no valid cache points -> compute as normal and save cache point
                                if(abs(parameters.phaseG) < 0.01)
                                {
                                    //isotropic case
                                    singleScattering = computeSingleScattering_Full(xNew, w, 16, singleGradientR, singleGradientG, singleGradientB, singleRadius);
                                     //only save cache point if we have reasonable density values
                                    if(density > 1e-6 && random() < 0.1 && singleRadius >= 0.0001)
                                    {
                                        insertLinear(xNew, singleScattering, singleGradientR, singleGradientG, singleGradientB, singleRadius);
                                    }
                                }
                                else
                                {
                                    //anisotropic case

                                    float coefficients[27];
                                    vec3 gradients[27];
                                    float radius;                                   
                                    singleScattering = SingleSH(xNew, w, 32, coefficients, gradients, radius);
                                    //insertResult(xNew, singleScattering, vec3(density, radius))
                                    if(density > 1e-6 && random() < 0.1 && radius >= 0.0001)
                                    {
                                        insertSingleSH(xNew, coefficients, gradients, radius);
                                    }
                                }

                                

                               
                                
                            }
                            
                        #endif

                        //multi scattering component
                        #ifdef USE_MULTISCATTERING
                            uint snapShotMulti = atomicAdd(linearMultiCounter.nextFreeNode, 0);
                            uint multiIndices[MAX_LINEAR_ENTRIES];
                            uint multiCount = queryMulti(x, multiIndices, snapShotMulti);
                            if(multiCount > 4)
                            {
                                //we have valid cache points  -> extrapolate
                                multiScattering = extrapolateMultiRadiance(multiIndices, multiCount, x);
                            }
                            else
                            {
                                //we have no valid cache points -> compute as normal and save cache point
                                vec3 multiGradientR;
                                vec3 multiGradientG;
                                vec3 multiGradientB;
                                float multiRadius;
                                multiScattering = computeMultiScattering_Full(xNew, w, 16, multiGradientR, multiGradientG, multiGradientB, multiRadius);
                                
                                //only save cache point if we have reasonable density values
                                if(density > 1e-6 && random() < 0.1 && multiRadius >= 0.0001)
                                {
                                    insertLinearMulti(x, multiScattering, multiGradientR, multiGradientG, multiGradientB, multiRadius);
                                }

                                
                            }
                           
                        #endif
                        


                    #else    //do not use caching       
                    
                        //nee contribution (optional)
                        #ifdef USE_DIRECT_COMPONENT
                            float pdf_w;
                            vec3 lightDir = parameters.sunDirection;
                            vec3 L_nee = sampleSkybox(lightDir) + sampleLight(lightDir);
                            float neeP = evaluatePhase(parameters.phaseG, w, lightDir);
                            vec3 neeGradient;
                            float Tr_nee = calculateTransmittance(xNew, lightDir, neeGradient);
                            neeContribution = throughput * Ps * neeP * L_nee * Tr_nee;
                            

                   
                        #endif

                        #ifdef USE_SINGLESCATTERING
                            
                            

                            if(abs(parameters.phaseG) < 0.001)
                            {
                                vec3 singleGradientR = vec3(0.0);
                                vec3 singleGradientG = vec3(0.0);
                                vec3 singleGradientB = vec3(0.0);
                                float singleRadius = 0.0;
                                singleScattering = computeSingleScattering_Full(xNew, w, 128, singleGradientR, singleGradientG, singleGradientB, singleRadius);
                            }
                            else
                            {
                                float coefficients[27];
                                vec3 gradients[27];
                                float radius;                                   
                                singleScattering = SingleSH(xNew, w, 32, coefficients, gradients, radius);
                                //insertResult(xNew, singleScattering, r[0], vec4(g[0], 0));
                            }
                            
                            //if(gl_GlobalInvocationID.xy == ivec2(127, 42)) insertResult(x, singleScattering, gl_GlobalInvocationID, vec4(1));
                            
                        #endif

                        #ifdef USE_MULTISCATTERING
                            if(abs(parameters.phaseG) < 0.001)
                            {
                                vec3 multiGradientR;
                                vec3 multiGradientG;
                                vec3 multiGradientB;
                                float multiRadius;
                                //multiScattering = computeMultiScattering_Full(xNew, w, 16, multiGradientR, multiGradientG, multiGradientB, multiRadius);
                                multiScattering = MultiSimple(xNew, w, 16);
                            }
                            else
                            {
                                float coeffR[9];
                                float coeffG[9];
                                float coeffB[9];
                                vec3 g[9];
                                multiScattering = MultiSHSimple(xNew, w, 32, coeffR, coeffG, coeffB);
                            }
                           
                        //multiScattering = multiSimple(xNew, w, 16);
                        #endif
                   
                  

                       
                   
                   //insertResult(xNew, singleScattering, multiScattering, vec4(neeContribution, 0));
                        


                   #endif // end of caching or no caching
                   vec3 Li = singleScattering + multiScattering + neeContribution;
                   return Li;

                   radiance += throughput * Ps * Li;
                   

                   /*
                   float pdf_w2;
                   vec3 wi = importanceSamplePhase(parameters.phaseG, w, pdf_w2);
                   float p = evaluatePhase(parameters.phaseG, w, wi);
                   //throughput *= (sigma_s / majorant) * (p / pdf_w2);
                   //bail if throughput too small
                   if(throughput < 1e-4)
                   {
                   return vec3(0.0);
                   }


                   w = wi; 

                   

                   if (rayBoxIntersect(parameters.boxMin, parameters.boxMax, x, w, tMin, tMax)) {
                       x += w*tMin;
                   
                       d = tMax - tMin;
                   }
                    
                    */
                   d -= t;
                }
                else //null event, continue along the path
                {
                    throughput *= Pn;
                    d -= t;
                }
                //throughput = clamp(throughput, 0.0, 1.0);
            }

            #ifdef EXIT_RADIANCE
                radiance += throughput * (sampleSkybox(w) + sampleLight(w));
            #endif

            #ifdef DIVIDE_BY_CONTRIBUTIONS
                return radiance * 1.0 / float(contributions);
            #else
                return radiance;
            #endif
        }
        else
        {
            radiance = (sampleSkybox(w) + sampleLight(w));
            return radiance;
        }

        
}


vec3 rayMarchDelta(vec3 x, vec3 w, ScatterEvent firstEvent)
{
       
        float majorant = parameters.extinction.x;
        float absorptionAlbedo = 1.0 - parameters.scatteringAlbedo.x;
        float scatteringAlbedo = parameters.scatteringAlbedo.x;
        float PA = absorptionAlbedo * parameters.extinction.x;
        float PS = scatteringAlbedo * parameters.extinction.x;

        vec3 radiance = vec3(0.0);
        float throughput=1.0;

    


        int i = 0;
        float tMin, tMax;
        int contributions = 1;
        if (rayBoxIntersect(parameters.boxMin, parameters.boxMax, x, w, tMin, tMax)) {
            x += w * tMin;
   
            float d = tMax - tMin;

            float pdf_x = 1;
            float transmittance = 1.0;
            

            while (true) {
    
                //sample free flight distance
                float t = -log(max(0.0000000001, 1 - random())) / majorant;
                //if we exceed medium, exit the loop
                if (t > d) {
                    break;
                }
                //get new position
                vec3 xNew = x + w * t;
                
               
               
                x = xNew;

               

                float density = sampleCloud(xNew);

                float sigma_a = PA * density;
                float sigma_s = PS * density;
                float sigma_n = majorant - parameters.extinction.x * density;
                float sigma_t = sigma_a + sigma_s;

                float Pa = sigma_a / majorant;
                float Ps = sigma_s / majorant;
                float Pn = sigma_n / majorant;
                float Pt = sigma_t / majorant;

                float xi = random();
    

                if (xi < Pa) {
                   throughput = 0.0;
                    if (!firstEvent.hasValue)
                    {                   
                        firstEvent.hasValue = true;
                    }
                   break;
                }

                if (xi < 1 - Pn) // scattering event accumulate single and multi scattering, then continue along the ray
                {
                    if (!firstEvent.hasValue) 
                    {                   
                        firstEvent.hasValue = true;
                    }
                    contributions++;
                    vec3 gradient;
                    float radius;
                    vec3 singleScattering = vec3(0.0);
                    vec3 multiScattering = vec3(0.0);
                    vec3 neeContribution = vec3(0.0);
                    vec3 comparison = vec3(0.0);
                    vec3 singleResult=vec3(0.0);



                    #ifdef USE_CACHING
                        //for current position, sample both caches to find valid cache points

                        #ifdef USE_DIRECT_COMPONENT

                            uint snapShotDirect = atomicAdd(linearDirectCounter.nextFreeNode, 0);
                            uint directIndices[MAX_LINEAR_ENTRIES];
                            uint directCount = queryDirect(xNew, directIndices, snapShotDirect);
                            vec3 directGradientR = vec3(0.0);
                            vec3 directGradientG = vec3(0.0);
                            vec3 directGradientB = vec3(0.0);
                            float directRadius = 0.0;
                            float epsilon = 1.0;

                            if(directCount > 4)
                            {
                                 neeContribution = extrapolateDirectRadiance(directIndices, directCount, xNew);
                            }
                            else
                            {
                                float pdf_w;
                                vec3 lightDir = parameters.sunDirection;
                                vec3 L_nee = sampleSkybox(lightDir) + sampleLight(lightDir);
                                float neeP = evaluatePhase(parameters.phaseG, w, lightDir);
                                vec3 neeGradient;
                                float Tr_nee = calculateTransmittance(xNew, lightDir, neeGradient);
                                neeContribution = throughput * Ps * neeP * L_nee * Tr_nee;
                                //radiance += neeContribution;

                                //gradient. for isotropic, no delta p term
                                float wNee = throughput * Ps * neeP;

                                directGradientR = wNee * L_nee.r * neeGradient;
                                directGradientG = wNee * L_nee.g * neeGradient;
                                directGradientB = wNee * L_nee.b * neeGradient;
                                vec3 gradientMagnitudes = vec3(length(directGradientR), length(directGradientG), length(directGradientB));
                                vec3 radianceSums = neeContribution;
                                if(gradientMagnitudes.x < 1e-6 || gradientMagnitudes.y < 1e-6 || gradientMagnitudes.z < 1e-6) 
                                {
                                    directRadius= 0.0;
                                }
                                else
                                {
                                    vec3 radiusVector = epsilon * (radianceSums / gradientMagnitudes);
                                    float radiusUnclamped = min(min(radiusVector.x, radiusVector.y), radiusVector.z);

                                    directRadius = clamp(radiusUnclamped, MIN_RADIUS, MAX_RADIUS);
                                }

                                if(density > 1e-6 && random() < 0.1 && directRadius >= 0.0001) insertLinearDirect(xNew, neeContribution, directGradientR, directGradientG, directGradientB, directRadius);
                            }

                            

                   
                        #endif
                        

                        //single scattering component
                        #ifdef USE_SINGLESCATTERING
                            uint snapShotSingle = atomicAdd(linearCounter.nextFreeNode, 0);
                            uint singleIndices[MAX_LINEAR_ENTRIES];
                            uint singleCount = querySingleSH(xNew, singleIndices, snapShotSingle);
                            vec3 singleGradientR = vec3(0.0);
                            vec3 singleGradientG = vec3(0.0);
                            vec3 singleGradientB = vec3(0.0);
                            float singleRadius = 0.0;
                            

                            if(singleCount > 4)
                            {
                                //we have valid cache points  -> extrapolate

                                if(abs(parameters.phaseG) < 0.01)
                                {
                                    //isotropic case
                                    singleScattering = extrapolateSingleRadianceALT(singleIndices, singleCount, xNew);
                                   
                                }
                                else
                                {
                                    singleScattering = extrapolateSingleSH_Version2(singleIndices, singleCount, xNew, -w);
                                    insertResult(xNew, -w, singleScattering, vec4(singleCount, 0, 0, 0));
                                    //anisotropic case
                                }
                                //check this function again
                                
                                //comparison = computeSingleScattering_Full(xNew, w, 16, singleGradientR, singleGradientG, singleGradientB, singleRadius);
                                //insertResult(xNew, singleScattering, comparison, vec4(singleRadius, singleCount, 0, 0));
                            }
                            else
                            {
                                //we have no valid cache points -> compute as normal and save cache point
                                if(abs(parameters.phaseG) < 0.01)
                                {
                                    //isotropic case
                                    singleScattering = computeSingleScattering_Full(xNew, w, 16, singleGradientR, singleGradientG, singleGradientB, singleRadius);
                                     //only save cache point if we have reasonable density values
                                    if(density > 1e-6 && random() < 0.1 && singleRadius >= 0.0001)
                                    {
                                        insertLinear(xNew, singleScattering, singleGradientR, singleGradientG, singleGradientB, singleRadius);
                                    }
                                }
                                else
                                {
                                    //anisotropic case
                                    float coefficients[27];
                                    vec3 gradients[27];
                                    float radius;                                   
                                    singleScattering = SingleSH(xNew, w, 32, coefficients, gradients, radius);
                                    if(density > 1e-6 && random() < 0.1 && radius >= 0.0001)
                                    {
                                        insertSingleSH(xNew, coefficients, gradients, radius);
                                    }
                                }

                                

                               
                                
                            }
                            
                        #endif

                        //multi scattering component
                        #ifdef USE_MULTISCATTERING
                            uint snapShotMulti = atomicAdd(linearMultiCounter.nextFreeNode, 0);
                            uint multiIndices[MAX_LINEAR_ENTRIES];
                            uint multiCount = queryMulti(x, multiIndices, snapShotMulti);
                            if(multiCount > 4)
                            {
                                //we have valid cache points  -> extrapolate
                                multiScattering = extrapolateMultiRadiance(multiIndices, multiCount, x);
                            }
                            else
                            {
                                //we have no valid cache points -> compute as normal and save cache point
                                vec3 multiGradientR;
                                vec3 multiGradientG;
                                vec3 multiGradientB;
                                float multiRadius;
                                multiScattering = computeMultiScattering_Full(xNew, w, 16, multiGradientR, multiGradientG, multiGradientB, multiRadius);
                                
                                //only save cache point if we have reasonable density values
                                if(density > 1e-6 && random() < 0.1 && multiRadius >= 0.0001)
                                {
                                    insertLinearMulti(x, multiScattering, multiGradientR, multiGradientG, multiGradientB, multiRadius);
                                }

                                
                            }
                           
                        #endif
                        


                    #else    //do not use caching       
                    
                        //nee contribution (optional)
                        #ifdef USE_DIRECT_COMPONENT
                            float pdf_w;
                            vec3 lightDir = parameters.sunDirection;
                            vec3 L_nee = sampleSkybox(lightDir) + sampleLight(lightDir);
                            float neeP = evaluatePhase(parameters.phaseG, w, lightDir);
                            vec3 neeGradient;
                            float Tr_nee = calculateTransmittance(xNew, lightDir, neeGradient);
                            neeContribution = throughput * Ps * neeP * L_nee * Tr_nee;
                            

                   
                        #endif

                        #ifdef USE_SINGLESCATTERING
                            
                            

                            if(abs(parameters.phaseG) < 0.001)
                            {
                                vec3 singleGradientR = vec3(0.0);
                                vec3 singleGradientG = vec3(0.0);
                                vec3 singleGradientB = vec3(0.0);
                                float singleRadius = 0.0;
                                singleScattering = computeSingleScattering_Full(xNew, w, 128, singleGradientR, singleGradientG, singleGradientB, singleRadius);
                            }
                            else
                            {
                                float coefficients[27];
                                vec3 gradients[27];
                                float radius;
                               
                                singleScattering = SingleSH(xNew, w, 32, coefficients, gradients, radius);
                                //insertResult(xNew, singleScattering, r[0], vec4(g[0], 0));
                            }
                            
                            //if(gl_GlobalInvocationID.xy == ivec2(127, 42)) insertResult(x, singleScattering, gl_GlobalInvocationID, vec4(1));
                            
                        #endif

                        #ifdef USE_MULTISCATTERING
                            if(abs(parameters.phaseG) < 0.001)
                            {
                                vec3 multiGradientR;
                                vec3 multiGradientG;
                                vec3 multiGradientB;
                                float multiRadius;
                                //multiScattering = computeMultiScattering_Full(xNew, w, 16, multiGradientR, multiGradientG, multiGradientB, multiRadius);
                                multiScattering = MultiSimple(xNew, w, 16);
                            }
                            else
                            {
                                float coeffR[9];
                                float coeffG[9];
                                float coeffB[9];
                                vec3 g[9];
                                multiScattering = MultiSHSimple(xNew, w, 32, coeffR, coeffG, coeffB);
                            }
                           
                        //multiScattering = multiSimple(xNew, w, 16);
                        #endif
                   
                  

                       
                   
                   //insertResult(xNew, singleScattering, multiScattering, vec4(neeContribution, 0));
                        


                   #endif // end of caching or no caching
                   vec3 Li = singleScattering + multiScattering + neeContribution;
                  
                   radiance += throughput * Ps * Li;
                   

                   
                   float pdf_w2;
                   vec3 wi = importanceSamplePhase(parameters.phaseG, w, pdf_w2);
                   float p = evaluatePhase(parameters.phaseG, w, wi);
                   throughput *= (sigma_s / majorant) * (p / pdf_w2);
                   //bail if throughput too small
                   if(throughput < 1e-4)
                   {
                        break;
                   }
                   /*
                   w = wi; 

                   

                   if (rayBoxIntersect(parameters.boxMin, parameters.boxMax, x, w, tMin, tMax)) {
                       x += w*tMin;
                   
                       d = tMax - tMin;
                   }
                   
                    */
                   d -= t;
                }
                else //null event, continue along the path
                {
                    throughput *= Pn;
                    d -= t;
                }
                //throughput = clamp(throughput, 0.0, 1.0);
            }

            #ifdef EXIT_RADIANCE
                radiance += throughput * (sampleSkybox(w) + sampleLight(w));
            #endif

            #ifdef DIVIDE_BY_CONTRIBUTIONS
                return radiance * 1.0 / float(contributions);
            #else
                return radiance;
            #endif
        }
        else
        {
            radiance = (sampleSkybox(w) + sampleLight(w));
            return radiance;
        }

        
}



vec3 probeDensity()
{
            vec3 dir = parameters.boxMax -parameters.boxMin;
            vec3 samplePos = parameters.boxMin + dir * 0.33;
            float density = sampleCloud(samplePos);
            vec3 densityGradient = sampleCloudGradient(samplePos);
            insertResult(samplePos, vec3(density), densityGradient, vec4(0));
   
    
    
    return vec3(0.0);
}
#endif
