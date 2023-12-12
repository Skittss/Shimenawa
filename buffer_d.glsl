// Simple buffer pass bloom. 
// Do this step in HDR space to avoid blooming non-bright objects. 
//  Combine 4 mipmap levels for progressive gaussian blur, and apply.
//  There is some artifacting due to the gaussian size.

// I honestly think this is a pretty bad way of doing this on shadertoy in hindsight
//  Would benefit greatly from individual buffers for each level + billinear sampling ...Too bad!

#define THRESH vec3(EXPOSURE * BLOOM_THRESHOLD)

vec3 mipmap(float mipmap_exp, vec2 offset, vec2 uv)
{
    vec2 pixel_size = 1.0 / vec2(iResolution.x, iResolution.y);
    float ds_factor = exp2(mipmap_exp);

    offset.x += mipmap_exp * BLOOM_MIPMAP_NUDGE + 1.0 - 2.0 / ds_factor; // Simplification of 1 - 1/(2^(n-1))
    vec2 scale = ds_factor * pixel_size;
    vec2 coord = (uv - offset) * ds_factor;
        
    // Toggle this line if you use another mipmap on top
    #if 0
    if (coord.x > 1.0 || coord.y > 1.0 || coord.x < 0.0 || coord.y < 0.0)
    #else
    if (coord.x > 1.0 || coord.x < 0.0)
    #endif
        return vec3(0.0);
    
    float tot_weight = 0.0;
    vec3 bloom = vec3(0.0);

    for (int i = -5; i < 5; i++) 
    for (int j = -5; j < 5; j++) 
    {
        float wg = pow(1.0 - length(vec2(i,j)) * 0.125, 6.0); // Apply pseudo-gaussian weights, TODO: this is inaccurate and kind of slow.

        vec2 blur_coord = vec2(i,j) * scale + ds_factor * pixel_size + coord;
                
        bloom = max(vec3(0.0), (HDR_MAX_COL * texture(iChannel0, blur_coord, mipmap_exp + 1.0).rgb - THRESH)) * wg + bloom;
        tot_weight += wg;
    }

    bloom /= tot_weight;

    return bloom;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 uv = fragCoord / iResolution.xy;
    
    vec3 blur = mipmap(1.0, vec2(0.0), uv);
        blur += mipmap(2.0, vec2(0.0), uv);
        blur += mipmap(3.0, vec2(0.0), uv);
        blur += mipmap(4.0, vec2(0.0), uv);
        blur += mipmap(5.0, vec2(0.0), uv);
    
    fragColor = vec4(blur, 1.0);
}