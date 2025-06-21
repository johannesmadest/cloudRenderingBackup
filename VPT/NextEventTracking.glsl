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

/**
 * For more details on spectral delta tracking, please refer to:
 * P. Kutz, R. Habel, Y. K. Li, and J. Novák. Spectral and decomposition tracking for rendering heterogeneous volumes.
 * ACM Trans. Graph., 36(4), Jul. 2017.
 */
#ifdef USE_NEXT_EVENT_TRACKING_SPECTRAL
vec3 nextEventTrackingSpectral(vec3 x, vec3 w, inout ScatterEvent firstEvent, bool onlyFirstEvent) {
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

    vec3 color = vec3(0);
    float bw_phase = 1.;

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
            vec4 densityEmission = sampleCloudDensityEmission(x);
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
                //return color;
#ifdef USE_EMISSION
                vec3 emission = sampleEmission(x);
#ifdef USE_ISOSURFACES
                return color + weights * emission;
#else
                return color + emission;
#endif
#else
#ifdef USE_TRANSFER_FUNCTION
#ifdef USE_ISOSURFACES
                return weights * densityEmission.rgb * parameters.emissionStrength;
#else
                return densityEmission.rgb * parameters.emissionStrength;
#endif
#endif
                return color; // weights * sigma_a / (majorant * Pa) * L_e; // 0 - No emission
#endif
            }

            if (xi < Pa + Ps) { // scattering event
                float pdf_w;
                vec3 next_w = importanceSamplePhase(parameters.phaseG, w, pdf_w);

                if (!firstEvent.hasValue) {
                    firstEvent.x = x;
                    firstEvent.pdf_x = 0; // TODO
                    firstEvent.w = next_w;
                    firstEvent.pdf_w = pdf_w;
                    firstEvent.hasValue = true;
                    firstEvent.density = density * maxComponent(parameters.extinction);
                    firstEvent.depth = tMax - d + t;
                }
                if (onlyFirstEvent){
                    return vec3(0);
                }

                float pdfLightNee; // only used for skybox.
                vec3 dirLightNee;
#if defined(USE_HEADLIGHT) || NUM_LIGHTS > 0
                // We are sampling the environment map or headlight with 50/50 chance.
                bool isSamplingHeadlight = (parameters.isEnvMapBlack != 0u) ? true : (random() > 0.5);
                float lightProbabilityFactor = parameters.isEnvMapBlack != 0u ? 1.0 : 2.0;
                float lightDistance = 0.0;
                uint lightIdx = 0;
                if (isSamplingHeadlight) {
                    dirLightNee = getHeadlightDirection(x, lightIdx, lightProbabilityFactor, lightDistance);
                } else {
#endif
                    dirLightNee = importanceSampleSkybox(pdfLightNee);
#if defined(USE_HEADLIGHT) || NUM_LIGHTS > 0
                }
#endif

                float pdf_nee_phase = evaluatePhase(parameters.phaseG, w, dirLightNee);
                w = next_w;

                weights *= sigma_s / (majorant * Ps);

#if defined(USE_HEADLIGHT) || NUM_LIGHTS > 0
                vec3 commonFactor = (lightProbabilityFactor * pdf_nee_phase) * min(weights, vec3(100000, 100000, 100000));
                if (isSamplingHeadlight) {
                    color +=
                            commonFactor * calculateTransmittanceDistance(x, dirLightNee, lightDistance)
                            * sampleHeadlight(x, lightIdx);
                } else {
                    color +=
                            commonFactor / pdfLightNee * calculateTransmittance(x, dirLightNee)
                            * (sampleSkybox(dirLightNee) + sampleLight(dirLightNee));
                }
#else
                // Normal NEE.
                color +=
                        (pdf_nee_phase / pdfLightNee * calculateTransmittance(x, dirLightNee))
                        * min(weights, vec3(100000, 100000, 100000)) * (sampleSkybox(dirLightNee) + sampleLight(dirLightNee));
#endif

                if (rayBoxIntersect(parameters.boxMin, parameters.boxMax, x, w, tMin, tMax)) {
                    x += w*tMin;
                    d = tMax - tMin;
                }
            } else {
                d -= t;
                weights *= sigma_n / (majorant * Pn);
            }
#if !defined(MAX_BASED_PROBABILITY) && !defined(AVG_BASED_PROBABILITY)
            weights = min(weights, vec3(100.0, 100.0, 100.0));
#endif
        }
    }

    if (!firstEvent.hasValue){
        color += bw_phase * min(weights, vec3(100000, 100000, 100000)) * (sampleSkybox(w) + sampleLight(w));
    }
    return color;
}
#endif

#ifdef USE_NEXT_EVENT_TRACKING
vec3 nextEventTracking(vec3 x, vec3 w, inout ScatterEvent firstEvent, bool onlyFirstEvent) {
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

    vec3 color = vec3(0.);

    int i = 0;
    float tMin, tMax;
    if (rayBoxIntersect(parameters.boxMin, parameters.boxMax, x, w, tMin, tMax)) {
        x += w * tMin;
        float d = tMax - tMin;
        float pdf_x = 1;

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
#else
            float density = sampleCloud(xNew);
#endif

#include "CheckIsosurfaceHit.glsl"

            x = xNew;

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
                    firstEvent.pdf_x = sigma_s * pdf_x;
                    firstEvent.w = w;
                    firstEvent.pdf_w = 0;
                    firstEvent.hasValue = true;
                    firstEvent.density = density * maxComponent(parameters.extinction);
                    firstEvent.depth = tMax - d + t;
                }
#ifdef USE_EMISSION
                vec3 emission = sampleEmission(x);
#ifdef USE_ISOSURFACES
                return color + weights * emission;
#else
                return color + emission;
#endif
#else
#ifdef USE_TRANSFER_FUNCTION
#ifdef USE_ISOSURFACES
                return color + weights * parameters.emissionStrength * densityEmission.rgb;
#else
                return color + parameters.emissionStrength * densityEmission.rgb;
#endif
#endif
                return color; // weights * sigma_a / (majorant * Pa) * L_e; // 0 - No emission
#endif
            }

            if (xi < 1 - Pn)// scattering event
            {
                float pdf_w;
                vec3 next_w = importanceSamplePhase(parameters.phaseG, w, pdf_w);

                if (!firstEvent.hasValue) {
                    firstEvent.x = x;
                    firstEvent.pdf_x = sigma_s * pdf_x;
                    firstEvent.w = next_w;
                    firstEvent.pdf_w = pdf_w;
                    firstEvent.hasValue = true;
                    firstEvent.density = density * maxComponent(parameters.extinction);
                    firstEvent.depth = tMax - d + t;
                }
                if (onlyFirstEvent){
                    return vec3(0);
                }

                float pdfLightNee; // only used for skybox.
                vec3 dirLightNee;
#if defined(USE_HEADLIGHT) || NUM_LIGHTS > 0
                // We are sampling the environment map or headlight with 50/50 chance.
                bool isSamplingHeadlight = (parameters.isEnvMapBlack != 0u) ? true : (random() > 0.5);
                float lightProbabilityFactor = parameters.isEnvMapBlack != 0u ? 1.0 : 2.0;
                float lightDistance = 0.0;
                uint lightIdx = 0;
                if (isSamplingHeadlight) {
                    dirLightNee = getHeadlightDirection(x, lightIdx, lightProbabilityFactor, lightDistance);
                } else {
#endif
                    dirLightNee = importanceSampleSkybox(pdfLightNee);
#if defined(USE_HEADLIGHT) || NUM_LIGHTS > 0
                }
#endif

                float pdf_nee_phase = evaluatePhase(parameters.phaseG, w, dirLightNee);
                w = next_w;

#if defined(USE_HEADLIGHT) || NUM_LIGHTS > 0
                float commonFactor = lightProbabilityFactor * pdf_nee_phase;
                vec3 colorNew;
                if (isSamplingHeadlight) {
                    colorNew =
                            (commonFactor * calculateTransmittanceDistance(x, dirLightNee, lightDistance))
                            * sampleHeadlight(x, lightIdx);
                } else {
                    colorNew =
                            (commonFactor / pdfLightNee * calculateTransmittance(x, dirLightNee))
                            * (sampleSkybox(dirLightNee) + sampleLight(dirLightNee));
                }
#else
                // Normal NEE.
                vec3 colorNew =
                        (pdf_nee_phase / pdfLightNee * calculateTransmittance(x, dirLightNee))
                        * (sampleSkybox(dirLightNee) + sampleLight(dirLightNee));
#endif

#ifdef USE_ISOSURFACES
                colorNew *= weights;
#endif
                color += colorNew;
                pdf_x *= exp(-majorant * t) * majorant * density;

                if (rayBoxIntersect(parameters.boxMin, parameters.boxMax, x, w, tMin, tMax)) {
                    x += w*tMin;
                    d = tMax - tMin;
                }
            } else {
                pdf_x *= exp(-majorant * t) * majorant * (1 - density);
                d -= t;
            }
        }
    }

    if (!firstEvent.hasValue){
        color += sampleSkybox(w) + sampleLight(w);
    }
    return color;
}
#endif



#ifdef CUSTOM
//spherical harmonics helper function
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




//------------------------------extrapolation functions--------------------------------------------------------------
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

vec3 extrapolateMultiRadianceALT(uint cacheIndices[MAX_LINEAR_ENTRIES], uint validEntries, vec3 position)
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
            LinearNode node = linearMultiBuffer.nodes[cacheIndices[i]];
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
            LinearNode node = linearMultiBuffer.nodes[cacheIndices[i]];
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


//-------------------------------query functions------------------------------------------------------------------------
int querySingleSH(vec3 position, out uint validEntries[MAX_LINEAR_ENTRIES], uint limit)
{
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

int queryMultiSH(vec3 position, out uint validEntries[MAX_LINEAR_ENTRIES], uint limit)
{
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
        //if(resultCounterBuffer.nextFreeNode < 128) insertResult(position, vec3(validEntryCounter), vec3(0), vec4(0));
            

        return int(validEntryCounter);


}




//----------------------buffer insertion functions---------------------------------------------------------------------

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

void insertResult(vec3 position, vec3 singleRadiance, vec3 multiRadiance, vec4 metaData)
    {
        
        if(resultCounterBuffer.nextFreeNode > 127) return;

        uint index = atomicAdd(resultCounterBuffer.nextFreeNode, 1);
        resultBuffer.results[index].location = vec4(position, 0);
        resultBuffer.results[index].singleRadiance = vec4(singleRadiance, 0);
        resultBuffer.results[index].multiRadiance = vec4(multiRadiance, 0);
        resultBuffer.results[index].metaData = metaData;
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

//transmittance calculation functions
float calculateTransmittance(vec3 x, vec3 w, out vec3 gradient) {
    float majorant = parameters.extinction.x;
    float absorptionAlbedo = 1.0 - parameters.scatteringAlbedo.x;
    float scatteringAlbedo = parameters.scatteringAlbedo.x;
    float PA = absorptionAlbedo * parameters.extinction.x;
    float PS = scatteringAlbedo * parameters.extinction.x;

    float transmittance = 1.0;
    float rr_factor = 1.0;
    int stepCounter = 0;

    vec3 gradTau = vec3(0.0);
    

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

            //insertResult(x, densityGradient, vec3(density), vec4(0));
            if(density <= 0.0) break; //assume no more cloud




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
            
            gradTau += delta_sigma_t / majorant;
            if(transmittance < 1e-5) break;
           
            

            d -= t;
        }
    }
    gradient = -transmittance * gradTau;
    //insertResult(x, w, gradient, vec4(transmittance, 0, 0, 0));
    return transmittance;
}


//------------------SINGLE SCATTERING FUNCTIONS-----------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------
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

             radius = clamp(radiusUnclamped, SINGLE_MIN_RADIUS, SINGLE_MAX_RADIUS);
        }
       

        //insertResult(gradientR, gradientG, gradientB, vec4(radiusUnclamped, 0, 0 , radius));
       
       
       return radiance;
    }


vec3 singleScattering(vec3 x, vec3 w, out IsoGradients gradients)
{
    vec3 Ls = vec3(0.0);
    vec3 gradR = vec3(0.0);
    vec3 gradG = vec3(0.0);
    vec3 gradB = vec3(0.0);

    gradients.R = vec3(0.0);
    gradients.G = vec3(0.0);
    gradients.B = vec3(0.0);

    //sun contribution
    vec3 lightDir = parameters.sunDirection;
    float phase_sun = evaluatePhase(parameters.phaseG, w, lightDir);
    vec3 delta_Tr_sun;
    float Tr_sun = calculateTransmittance(x, lightDir, delta_Tr_sun);
    vec3 L_sun = sampleLight(lightDir);
    Ls += phase_sun * Tr_sun * L_sun;

    //gradient
    //vec3 delta_phase_sun = evaluatePhaseGradient(parameters.phaseG, w, lightDir);
    vec3 delta_L_sun = L_sun * delta_Tr_sun;

    gradR += phase_sun * delta_L_sun.r;
    gradG += phase_sun * delta_L_sun.g;
    gradB += phase_sun * delta_L_sun.b;

    //skybox contribution
    float pdf_w;
    vec3 skyDir = importanceSamplePhase(parameters.phaseG, w, pdf_w);
    float phase_sky = evaluatePhase(parameters.phaseG, w, skyDir);
    vec3 delta_Tr_sky;
    float Tr_sky = calculateTransmittance(x, skyDir, delta_Tr_sky);

    vec3 L_sky = sampleSkybox(skyDir);
    Ls += (phase_sky * Tr_sky * L_sky) / pdf_w;

    //gradient
    vec3 delta_L_sky = L_sky * delta_Tr_sky;

    gradR += phase_sky * delta_L_sky.r / pdf_w;
    gradG += phase_sky * delta_L_sky.g / pdf_w;
    gradB += phase_sky * delta_L_sky.b / pdf_w;

    gradients.R = gradR;
    gradients.G = gradG;
    gradients.B = gradB;

    //insertResult(x, delta_Tr_sun, delta_Tr_sky, vec4(Tr_sun, Tr_sky, 1, 1));


    return Ls;
}

vec3 cacheSingleScattering(vec3 x, vec3 w, int N)
{
    float invN = 1.0 / float(N);
    IsoGradients gradients;
    gradients.R = vec3(0.0);
    gradients.G = vec3(0.0);
    gradients.B = vec3(0.0);
    vec3 radiance = vec3(0.0);
    vec3 radianceSums = vec3(0.0);
    vec3 gradientMagnitudes = vec3(0.0);
    float epsilon = 1;
    ivec2 localID = ivec2(gl_GlobalInvocationID.xy);
    ivec2 imageCoord = pc.tileOffset + localID;

    vec3 voxelSize = (parameters.boxMax - parameters.boxMin) / parameters.gridResolution;

    for(int i = 0; i< N; i++)
    {
        IsoGradients resultGradients;
        resultGradients.R = vec3(0.0);
        resultGradients.G = vec3(0.0);
        resultGradients.B = vec3(0.0);
        vec3 result = singleScattering(x, w, resultGradients);
        radiance += result;
        radianceSums += result;

        gradients.R += resultGradients.R;
        gradients.G += resultGradients.G;
        gradients.B += resultGradients.B;

        gradientMagnitudes.r += length(resultGradients.R);
        gradientMagnitudes.g += length(resultGradients.G);
        gradientMagnitudes.b += length(resultGradients.B);
                
    }
    //normalizing things
    radiance *= invN;
    

    gradients.R *= invN;
    gradients.G *= invN;
    gradients.B *= invN;

    //radius calculation
    vec3 radiusVector = vec3(0.0);
    float radiusUnclamped = 0.0;
    float radius = 0.0;
    if(gradientMagnitudes.x < 1e-6 || gradientMagnitudes.y < 1e-6 || gradientMagnitudes.z < 1e-6) 
    {
        radius= 0.0;
    }
    else
    {
         radiusVector = epsilon * (radianceSums / gradientMagnitudes) * voxelSize;
         radiusUnclamped = min(min(radiusVector.x, radiusVector.y), radiusVector.z);
         radius = clamp(radiusUnclamped, SINGLE_MIN_RADIUS, SINGLE_MAX_RADIUS);
    }
    float highest = 0.0;
    for(int i=0; i<3; i++)
    {
        highest = max(abs(gradients.R[i]), highest);
        highest = max(abs(gradients.G[i]), highest);
        highest = max(abs(gradients.B[i]), highest);
    }
    insertResult(gradients.R, gradients.G, gradients.B, vec4(gradientMagnitudes, highest));
    //insert data into cache point
    if(radius > 0.0001)
    {
        insertLinear(x, radiance, gradients.R, gradients.G, gradients.B, radius);
    }

    return radiance;
    

}

vec3 singleScatteringFromCache(vec3 x, vec3 w)
{
    // 0. Initialization
    vec3 Ls = vec3(0.0);

    // 1. get snapshot of how many cache points are there
    uint snapShotSingle = atomicAdd(linearCounter.nextFreeNode, 0);
                            
    // 2. query for available cache points using radius
    uint singleIndices[MAX_LINEAR_ENTRIES];
    uint singleCount = querySingleSH(x, singleIndices, snapShotSingle);

    // 3a. If no cache point is found, create using cacheSingleScattering and return its result
    if(singleCount < 2)
    {
        Ls = cacheSingleScattering(x, w, 16);
    }
    // 3b.If enough cache points are found, extrapolate and return the result
    else
    {
        Ls = extrapolateSingleRadianceALT(singleIndices, singleCount, x);
    }
    //insertResult(x, w, Ls, vec4(singleCount, 0, 0, 0));
    return Ls;
}

vec3 singleScatteringSH(vec3 x, vec3 w, out Anisotropic data)
{
    vec3 Ls = vec3(0.0);

    float weightR = 0;
    float weightG = 0;
    float weightB = 0;

    vec3 gradR = vec3(0.0);
    vec3 gradG = vec3(0.0);
    vec3 gradB = vec3(0.0);

    for(int i=0; i< 9; i++)
    {
        data.R[i] = 0;
        data.G[i] = 0;
        data.B[i] = 0;
        data.gradientR[i] = vec3(0);
        data.gradientG[i] = vec3(0);
        data.gradientB[i] = vec3(0);
    }

    //project into sh
    float Y[9];
    vec3 dir = normalize(wi);
    evalSH2(dir, Y);
    

    //sun contribution
    vec3 lightDir = parameters.sunDirection;
    float phase_sun = evaluatePhase(parameters.phaseG, w, lightDir);
    vec3 delta_Tr_sun;
    float Tr_sun = calculateTransmittance(x, lightDir, delta_Tr_sun);
    vec3 L_sun = sampleLight(lightDir);
    Ls += phase_sun * Tr_sun * L_sun;

    //gradient
    //vec3 delta_phase_sun = evaluatePhaseGradient(parameters.phaseG, w, lightDir);
    vec3 delta_L_sun = L_sun * delta_Tr_sun;

    gradR += phase_sun * delta_L_sun.r;
    gradG += phase_sun * delta_L_sun.g;
    gradB += phase_sun * delta_L_sun.b;

    //skybox contribution
    float pdf_w;
    vec3 skyDir = importanceSamplePhase(parameters.phaseG, w, pdf_w);
    float phase_sky = evaluatePhase(parameters.phaseG, w, skyDir);
    vec3 delta_Tr_sky;
    float Tr_sky = calculateTransmittance(x, skyDir, delta_Tr_sky);

    vec3 L_sky = sampleSkybox(skyDir);
    Ls += (phase_sky * Tr_sky * L_sky) / pdf_w;

    //gradient
    vec3 delta_L_sky = L_sky * delta_Tr_sky;

    gradR += phase_sky * delta_L_sky.r / pdf_w;
    gradG += phase_sky * delta_L_sky.g / pdf_w;
    gradB += phase_sky * delta_L_sky.b / pdf_w;

    gradients.R = gradR;
    gradients.G = gradG;
    gradients.B = gradB;

    //insertResult(x, delta_Tr_sun, delta_Tr_sky, vec4(Tr_sun, Tr_sky, 1, 1));


    return Ls;
}

//------------------ MULTI SCATTERING FUNCTIONS ----------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------


vec3 multiScattering(vec3 x, vec3 w, out IsoGradients gradients)
{
    //get variables from parameters
    float majorant = parameters.extinction.x;
    float absorptionAlbedo = 1.0 - parameters.scatteringAlbedo.x;
    float scatteringAlbedo = parameters.scatteringAlbedo.x;
    float PA = absorptionAlbedo * parameters.extinction.x;
    float PS = scatteringAlbedo * parameters.extinction.x;


    int i = 0;
    float tMin, tMax;
    vec3 throughput = vec3(1.0);
    float transmittance = 1.0;

    vec3 radiance = vec3(0.0);

    gradients.R = vec3(0.0);
    gradients.G = vec3(0.0);
    gradients.B = vec3(0.0);

    vec3 gradTau = vec3(0);


    if (rayBoxIntersect(parameters.boxMin, parameters.boxMax, x, w, tMin, tMax)) {
        x += w * tMin;

        float d = tMax - tMin;

        float pdf_x = 1;
       

        while (true) {

            float t = -log(max(0.0000000001, 1 - random())) / majorant;
            

            if (t > d) {
                radiance += throughput * (sampleSkybox(w) + sampleLight(w));
                break;
            }

            vec3 xNew = x + w * t;

            float density = sampleCloud(xNew);

            vec3 densityGradient = sampleCloudGradient(xNew);



            x = xNew;

            float sigma_a = PA * density;
            float sigma_s = PS * density;
            float sigma_t = sigma_a + sigma_s;
            float sigma_n = majorant - parameters.extinction.x * density;

            vec3 delta_sigma_s = PS * densityGradient;
            vec3 delta_sigma_t = (PA + PS) * densityGradient;

            gradTau += delta_sigma_t / majorant;

            if(sigma_t <= 0.0)
            {
                radiance += throughput * (sampleSkybox(w) + sampleLight(w));
                break;
            }

            float Pa = sigma_a / majorant;
            float Ps = sigma_s / majorant;
            float Pn = sigma_n / majorant;

            float xi = random();

            if (xi < Pa) { //absorption event -- ray ends and only emission is returned
                //return vec3(0); // weights * sigma_a / (majorant * Pa) * L_e; // 0 - No emission
                break;
            }

            if (xi < 1 - Pn) // scattering event -- here i want to compute single scattering and multiple scattering radiances or use cache values to extrapolate
            {
                i++;
                throughput *= Ps;
                transmittance *= Ps;

                float pdf_w;
                vec3 wi = importanceSamplePhase(parameters.phaseG, w, pdf_w);
                float p = evaluatePhase(parameters.phaseG, w, wi);
                
                throughput *= (p / pdf_w);
                IsoGradients gradLi;
                vec3 Li = singleScattering(x, w, gradLi);
                radiance += throughput * Li;


                //gradient terms
                vec3 delta_Tr = -transmittance * gradTau;
                
                //R channel
                gradients.R += p * delta_Tr * sigma_s * Li.r
                            +  p * transmittance * delta_sigma_s * Li.r
                            +  p * transmittance * gradLi.R;

                gradients.G += p * delta_Tr * sigma_s * Li.g
                            +  p * transmittance * delta_sigma_s * Li.g
                            +  p * transmittance * gradLi.G;

                gradients.B += p * delta_Tr * sigma_s * Li.b
                            +  p * transmittance * delta_sigma_s * Li.b
                            +  p * transmittance * gradLi.B;

                //insertResult(delta_Tr, delta_sigma_s, Li, vec4(p, transmittance, sigma_s, 0));


                w = wi;

                if (rayBoxIntersect(parameters.boxMin, parameters.boxMax, x, w, tMin, tMax)) {
                    x += w*tMin;

                    d = tMax - tMin;
                }
                else
                {
                    radiance += throughput * (sampleSkybox(w) + sampleLight(w));
                    break;
                }
          
                
            } else { // null event

                d -= t;
            }

        }
    }
    //insertResult(radiance, gradients.R, gradients.G, vec4(gradients.B, 0));
    return radiance;

}

vec3 cacheMultiScattering(vec3 x, vec3 w, int N)
{
    float invN = 1.0 / float(N);
    IsoGradients gradients;
    gradients.R = vec3(0.0);
    gradients.G = vec3(0.0);
    gradients.B = vec3(0.0);
    vec3 radiance = vec3(0.0);
    vec3 radianceSums = vec3(0.0);
    vec3 gradientMagnitudes = vec3(0.0);
    float epsilon = 1;
    ivec2 localID = ivec2(gl_GlobalInvocationID.xy);
    ivec2 imageCoord = pc.tileOffset + localID;

    vec3 voxelSize = (parameters.boxMax - parameters.boxMin) / parameters.gridResolution;

    for(int i = 0; i< N; i++)
    {
        IsoGradients resultGradients;
        resultGradients.R = vec3(0.0);
        resultGradients.G = vec3(0.0);
        resultGradients.B = vec3(0.0);
        vec3 result = multiScattering(x, w, resultGradients);
        radiance += result;
        radianceSums += result;

        gradients.R += resultGradients.R;
        gradients.G += resultGradients.G;
        gradients.B += resultGradients.B;

        gradientMagnitudes.r += length(resultGradients.R);
        gradientMagnitudes.g += length(resultGradients.G);
        gradientMagnitudes.b += length(resultGradients.B);
                
    }
    //normalizing things
    radiance *= invN;

    gradients.R *= invN;
    gradients.G *= invN;
    gradients.B *= invN;

    //radius calculation
    vec3 radiusVector = vec3(0.0);
    float radiusUnclamped = 0.0;
    float radius = 0.0;
    if(gradientMagnitudes.x < 1e-6 || gradientMagnitudes.y < 1e-6 || gradientMagnitudes.z < 1e-6) 
    {
        radius= 0.0;
    }
    else
    {
         radiusVector = epsilon * (radianceSums / gradientMagnitudes) * voxelSize;
         radiusUnclamped = min(min(radiusVector.x, radiusVector.y), radiusVector.z);
         radius = clamp(radiusUnclamped, MULTI_MIN_RADIUS, MULTI_MAX_RADIUS);
    }
    float highest = 0.0;
    for(int i=0; i<3; i++)
    {
        highest = max(abs(gradients.R[i]), highest);
        highest = max(abs(gradients.G[i]), highest);
        highest = max(abs(gradients.B[i]), highest);
    }
    //insertResult(gradients.R, gradients.G, gradients.B, vec4(radius, highest, 1, 1));
    //insert data into cache point
    if(radius > 0.0001)
    {
        insertLinearMulti(x, radiance, gradients.R, gradients.G, gradients.B, radius);
    }

    return radiance;
}

vec3 multiScatteringFromCache(vec3 x, vec3 w)
{
    // 0. Initialization
    vec3 Lm = vec3(0.0);

    // 1. get snapshot of how many cache points are there
    uint snapShotMulti = atomicAdd(linearMultiCounter.nextFreeNode, 0);
                            
    // 2. query for available cache points using radius
    uint multiIndices[MAX_LINEAR_ENTRIES];
    uint multiCount = queryMultiSH(x, multiIndices, snapShotMulti);

    // 3a. If no cache point is found, create using cacheMultiScattering and return its result
    if(multiCount < 2)
    {
        Lm = cacheMultiScattering(x, w, 8);
    }
    // 3b.If enough cache points are found, extrapolate and return the result
    else
    {
        Lm = extrapolateMultiRadianceALT(multiIndices, multiCount, x);
    }
    //insertResult(x, w, Ls, vec4(multiCount, 0, 0, 0));
    return Lm;
}




//main function
vec3 deltaTracer(vec3 x, vec3 w, ScatterEvent event)
{
   

    //get variables from parameters
    float majorant = parameters.extinction.x;
    float absorptionAlbedo = 1.0 - parameters.scatteringAlbedo.x;
    float scatteringAlbedo = parameters.scatteringAlbedo.x;
    float PA = absorptionAlbedo * parameters.extinction.x;
    float PS = scatteringAlbedo * parameters.extinction.x;

    ivec2 localID = ivec2(gl_GlobalInvocationID.xy);
    ivec2 imageCoord = pc.tileOffset + localID;


    int i = 0;
    float tMin, tMax;
    float transmittance = 1.0;
    vec3 radiance = vec3(0.0);
    bool earlyBreak = false;
    if (rayBoxIntersect(parameters.boxMin, parameters.boxMax, x, w, tMin, tMax)) {
        x += w * tMin;

        float d = tMax - tMin;

        float pdf_x = 1;
       

        while (true) {

            float t = -log(max(0.0000000001, 1 - random())) / majorant;
            

            if (t > d) {
                break;
            }

            vec3 xNew = x + w * t;

            float density = sampleCloud(xNew);

            x = xNew;

            float sigma_a = PA * density;
            float sigma_s = PS * density;
            float sigma_t = sigma_a + sigma_s;
            float sigma_n = majorant - parameters.extinction.x * density;

            float Pa = sigma_a / majorant;
            float Ps = sigma_s / majorant;
            float Pn = sigma_n / majorant;

            float xi = random();

            if (xi < Pa) { //absorption event -- ray ends and only emission is returned
                return vec3(0); // weights * sigma_a / (majorant * Pa) * L_e; // 0 - No emission
            }

            if (xi < 1 - Pn) // scattering event -- here i want to compute single scattering and multiple scattering radiances or use cache values to extrapolate
            {
                vec3 Ls = vec3(0.0);
                vec3 Lm = vec3(0.0);

                //single scattering only at first bounce
                #ifdef USE_SINGLESCATTERING
                if(event.bounces >= 0)
                {
                    #ifdef USE_CACHING
                        Ls = singleScatteringFromCache(x, w);
                    #else
                        IsoGradients g;
                        Ls = singleScattering(x, w, g);
                    #endif
                }
                #endif

                #ifdef USE_MULTISCATTERING
                    #ifdef USE_CACHING
                        Lm = multiScatteringFromCache(x, w);
                    #else
                        IsoGradients g_multi;
                        Lm = multiScattering(x, w, g_multi);
                    #endif
                    
                    
                    
                #endif
                radiance += transmittance * (Ls + Lm);
                //if((Ls != vec3(0.0) && Lm != vec3(0.0)))insertResult(x, Ls, Lm, vec4(event.bounces, transmittance, imageCoord));
                #if defined(USE_SINGLESCATTERING) && defined(USE_MULTISCATTERING)
                //if(imageCoord == ivec2(660, 340)) insertResult(x, Ls, Lm, vec4(event.bounces, transmittance, imageCoord));
                #endif

                event.bounces++;
                
                float pdf_w;
                w = importanceSamplePhase(parameters.phaseG, w, pdf_w);              

                if (rayBoxIntersect(parameters.boxMin, parameters.boxMax, x, w, tMin, tMax)) {
                    x += w*tMin;

                    d = tMax - tMin;
                }
                
            } else { // null event

                d -= t;
            }
            transmittance *= exp(-sigma_t * t);
            if(transmittance< 1e-6)
            {
                earlyBreak = true;
                break;
            }

        }
    }


    if(!earlyBreak) radiance += transmittance * (sampleSkybox(w) + sampleLight(w));

    return radiance;

}

#endif
