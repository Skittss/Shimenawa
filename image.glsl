// Post processing

#define colorRange 1.0

// This curve is kinda trash icl
vec3 aces(vec3 x) {
  const float a = 2.51;
  const float b = 0.03;
  const float c = 2.43;
  const float d = 0.59;
  const float e = 0.14;
  return clamp((x * (a * x + b)) / (x * (c * x + d) + e), 0.0, 1.0);
}

// Really useful visualisation of the following:
// https://www.desmos.com/calculator/gslcdxvipg
vec3 gt_tonemap(vec3 x) {
    const vec3 P = vec3(GT_MAX_BRIGHTNESS           );
	const vec3 a = vec3(GT_CONTRAST                 );
	const vec3 m = vec3(GT_LINEAR_OFFSET            );
	const vec3 l = vec3(GT_LINEAR_LENGTH            );
	const vec3 c = vec3(GT_BLACK_TIGHTNESS_CURVATURE);
	const vec3 b = vec3(GT_BLACK_TIGHTNESS_OFFSET   );
    
	vec3 l0 = ((P - m) * l) / a;
    vec3 S0 = m + l0;
    vec3 S1 = m + a * l0;
    vec3 C2 = (a * P) / (P - S1);
    vec3 CP = -C2 / P;

    vec3 w0 = 1.0f - smoothstep(x, vec3(0.0f), m);
    vec3 w2 = smoothstep(x, m + l0, m + l0);
    vec3 w1 = 1.0f - w0 - w2;

    vec3 T = m * pow(x / m, c) + b;
    vec3 S = P - (P - S1) * exp(CP * (x - S0));
    vec3 L = m + a * (x - m);

    // Should probably clamp this
    return T * w0 + L * w1 + S * w2;
}

vec3 bloom_mipmap(float mipmap_exp, vec2 offset, vec2 uv) 
{
    offset.x += 1.0 - 2.0 / exp2(mipmap_exp); // Simplification of 1 - 1/(2^(n-1))
    
    return texture(iChannel1, uv * exp2(-mipmap_exp) + offset).rgb;
}

vec3 get_bloom(vec2 uv)
{
    vec3 blur = bloom_mipmap(1.0, vec2(0.0, 0.0   ), uv);
        blur += bloom_mipmap(2.0, vec2(0.0, 0.5   ), uv);
        blur += bloom_mipmap(3.0, vec2(0.0, 0.75  ), uv);
        blur += bloom_mipmap(4.0, vec2(0.0, 0.875 ), uv);
        blur += bloom_mipmap(5.0, vec2(0.0, 0.9375), uv);
    
    return blur * colorRange;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 uv = fragCoord.xy / iResolution.xy;
    
    vec3 col = HDR_MAX_COL * texture(iChannel0, uv).rgb;
    vec3 bloom_map = HDR_MAX_COL * texture(iChannel1, uv).rgb;
    vec3 bloom = get_bloom(uv);
    
    #if POSTPROCESS
    {
        // HDR Post-processing
        col = col + BLOOM_INTENSITY * bloom; // Bloom
        col = gt_tonemap(col);               // Tonemap
        col = pow(col, vec3(GAMMA));         // Gamma

        // LDR Post-processing
        col = CONTRAST * (col - vec3(0.5)) + vec3(0.5);  // Contrast
        col = BRIGHTNESS + col;                          // Brightness
    }
    #else
    {
        col = pow(col, vec3(GAMMA));
    }
    #endif

    fragColor = vec4(col, 1.0);
}