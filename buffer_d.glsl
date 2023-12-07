// Simple bloom filter. 
// Do this step in HDR space to avoid blooming non-bright objects. 
//  Combine 4 mipmap levels for progressive gaussian blur, and apply.
//  There is some artifacting due to the gaussian size.

#define THRESH vec3(EXPOSURE * BLOOM_THRESHOLD)

vec3 mipmap(float mipmap_exp, vec2 offset, vec2 uv)
{
    vec2 pixel_size = 1.0 / vec2(iResolution.x, iResolution.y);
    float ds_factor = exp2(mipmap_exp);

    offset.x += 1.0 - 2.0 / ds_factor; // Simplification of 1 - 1/(2^(n-1))
    vec2 scale = ds_factor * pixel_size;
    vec2 coord = (uv - offset) * ds_factor;

    // Don't render any more pixels outside the mipmap rect (otherwise boundary stretches).
    if (any(greaterThanEqual(abs(coord - 0.5), scale + 0.5)))
        return vec3(0.0);
    
    float tot_weight = 0.0;
    vec3 bloom = vec3(0.0);

    // TODO: There is some colour bleeding at mipmap edges because of the tight packing of the mipmap levels
    //         This is circumvented by diagonally offsetting mipmap levels, but needs fixing if more than one mipmap is used
    //         in this pass.
    for (int i = -5; i < 5; i++) 
    for (int j = -5; j < 5; j++) 
    {
        float wg = pow(1.0 - length(vec2(i,j)) * 0.125, 6.0); // Apply pseudo-gaussian weights, TODO: this is inaccurate and kind of slow.

        vec2 blur_offset = vec2(i,j) * scale + ds_factor * pixel_size;
        bloom = max(vec3(0.0), (HDR_MAX_COL * texture(iChannel0, blur_offset + coord, mipmap_exp + 1.0).rgb - THRESH)) * wg + bloom;
        tot_weight += wg;
    }

    bloom /= tot_weight;

    return bloom;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 uv = fragCoord / iResolution.xy;
    
    vec3 blur = mipmap(1.0, vec2(0.0, 0.0   ), uv);
        blur += mipmap(2.0, vec2(0.0, 0.5   ), uv);
        blur += mipmap(3.0, vec2(0.0, 0.75  ), uv);
        blur += mipmap(4.0, vec2(0.0, 0.875 ), uv);
        blur += mipmap(5.0, vec2(0.0, 0.9375), uv);
    
    fragColor = vec4(blur, 1.0);
}