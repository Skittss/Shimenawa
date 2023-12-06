// Simple bloom filter. 
// Do this step in HDR space to avoid blooming non-bright objects. Combine 4 mipmap levels for blur, and apply.
// There is some artifacting due to the gaussian size.

#define THRESH vec3(EXPOSURE * BLOOM_THRESHOLD)

vec3 makeBloom(float lod, vec2 offset, vec2 uv){
    
    vec2 pixelSize = 1.0 / vec2(iResolution.x, iResolution.y);

    offset += pixelSize;

    float lodFactor = exp2(lod);

    vec3 bloom = vec3(0.0);
    vec2 scale = lodFactor * pixelSize;

    vec2 coord = (uv.xy  -offset) * lodFactor;
    float totalWeight = 0.0;

    if (any(greaterThanEqual(abs(coord - 0.5), scale + 0.5)))
        return vec3(0.0);

    for (int i = -5; i < 5; i++) {
        for (int j = -5; j < 5; j++) {

            float wg = pow(1.0-length(vec2(i,j)) * 0.125,6.0); // Apply pseudo-gaussian weights

            bloom = max(vec3(0.0), (HDR_MAX_COL * texture(iChannel0, vec2(i,j) * scale + lodFactor * pixelSize + coord, lod).rgb - THRESH)) * wg + bloom;
            totalWeight += wg;

        }
    }

    bloom /= totalWeight;

    return bloom;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 uv = fragCoord / iResolution.xy;
    
	vec3 blur = makeBloom(2.0, vec2(0.0,0.0), uv);
		blur += makeBloom(3.0, vec2(0.3,0.0), uv);
		blur += makeBloom(4.0, vec2(0.0,0.3), uv);
		blur += makeBloom(5.0, vec2(0.1,0.3), uv);
		blur += makeBloom(6.0, vec2(0.2,0.3), uv);

    fragColor = vec4(blur, 1.0);
    //fragColor = vec4(pow(blur, vec3(1.0 / 2.2)),1.0);
}