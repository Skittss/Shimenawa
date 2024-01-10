/*                                  

                                                 ████                                               
                   ██                            █████                              █████████       
       █████      ██████             ███         ████████              ████   ██  █████   █████     
      ███████     ███████           ██████  ██████████                 █████  ██  ████████████      
        █████    ███████             █████    █ ██████████            ██████████  ███████████       
                ██  █████            █████   ██ █████  █████         █████    ███ ██████████        
              ███   ████████               █  ██ ███████████        ████████  ███████████           
     █████  ████████████████        ██████ ██  █████████████       ██████████    ██████ █████       
   █████████████████████          ████████ ████ ██████████            ██████    ███████ ██████      
    ████████   ██ ███████          █████     ███████  ████           ███████ █████████ ██████       
      █  ██     ██████████            ████    █  ████████████         ██████████ ██████████         
     █  ██   ███████████               ██████████████████████          ███ ███████████   ██         
     █ ███  ███   ████             ███████████   ███  █               █ ███ █   ██████     ██       
    █ ███       ███████ █████      ████████      ███  █              █  ███      █████      ███     
    ████   ████████████████████      ██ ████████ ██        ██         ████      ████████████████    
    ███  ████████          ████               █████████████████         ████           █████████    
                                                                                                    

    注連縄 (Shimenawa) by Henry A. 
    
    --LAYOUT------------------------------------------------------------------------------------
    
      COMMON: Colour schemes, Rendering settings, SDF primitives, Intersectors, and util funcs.
       
    BUFFER B: Perlin-Worley texture atlas generation for clouds. 
              (Source: https://www.shadertoy.com/view/3sffzj)
              
    BUFFER C: Scene rendering
   
    BUFFER D: HDR Bloom pass
    
       IMAGE: LDR Post processing (Tonemap, Gamma correction, etc.)
       
    --------------------------------------------------------------------------------------------
    
*/

// Post processing

#define colorRange 1.0

// ACES for that sweet film look ;)
vec3 aces(vec3 x) {
  const float a = 2.51;
  const float b = 0.03;
  const float c = 2.43;
  const float d = 0.59;
  const float e = 0.14;
  return clamp((x * (a * x + b)) / (x * (c * x + d) + e), 0.0, 1.0);
}

// GT Tonemapper from "HDR theory and practice", Hajime Uchimura:
//  Curve Visualisation: https://www.desmos.com/calculator/gslcdxvipg
//  GT tonemap is more linear in LDR colour space so can be easier to work with.
vec3 gt_tonemap(vec3 x) {
    const vec3 P = vec3(GT_MAX_BRIGHTNESS           );
	const vec3 a = vec3(GT_CONTRAST                 );
	const vec3 m = vec3(GT_LINEAR_OFFSET            );
	const vec3 l = vec3(GT_LINEAR_LENGTH            );
	const vec3 c = vec3(GT_BLACK_TIGHTNESS_CURVATURE);
	const vec3 b = vec3(GT_BLACK_TIGHTNESS_OFFSET   );
    
    vec3 l0 = ((P - m) * l) / a;
    vec3 L0 = m - m / a;
    vec3 L1 = m + (1.0 - m) / a;
    vec3 S0 = m + l0;
    vec3 S1 = m + a * l0;
    vec3 C2 = (a * P) / (P - S1);
    vec3 CP = -C2 / P;

    vec3 w0 = 1.0 - smoothstep(vec3(0.0), m, x);
    vec3 w2 = step(m + l0, x);
    vec3 w1 = 1.0 - w0 - w2;

    vec3 T = m * pow(x / m, c) + b;
    vec3 S = P - (P - S1) * exp(CP * (x - S0));
    vec3 L = m + a * (x - m);

    return T * w0 + L * w1 + S * w2;
}

vec3 bloom_mipmap(float mipmap_exp, vec2 offset, vec2 uv) 
{
    // This reverse mapping is inaccurate, and causes slight edge bleeding T^T
    float ds_factor = exp2(mipmap_exp);
    offset.x += mipmap_exp * BLOOM_MIPMAP_NUDGE + 1.0 - 2.0 / ds_factor; // Simplification of 1 - 1/(2^(n-1))
    
    return texture(iChannel1, uv / ds_factor + offset).rgb;
}

vec3 get_bloom(vec2 uv)
{ 
    vec3 blur  = bloom_mipmap(1.0, vec2(0.0), uv);
         blur += bloom_mipmap(2.0, vec2(0.0), uv);
         blur += bloom_mipmap(3.0, vec2(0.0), uv);
         blur += bloom_mipmap(4.0, vec2(0.0), uv);
         blur += bloom_mipmap(5.0, vec2(0.0), uv);
    
    return blur * colorRange;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 uv = fragCoord.xy / iResolution.xy;
    
    vec3 col = HDR_MAX_COL * texture(iChannel0, uv).rgb;
    vec3 bloom_map = HDR_MAX_COL * texture(iChannel1, uv).rgb;
    vec3 bloom = get_bloom(uv);
    
    // HDR Post-processing
    #if POSTPROCESS
    col = col + BLOOM_INTENSITY * bloom; // Bloom
    #endif
    
    // Tonemap
    #ifdef USE_ACES_TONEMAP
    col = aces(col);
    #elifdef USE_GT_TONEMAP
    col = gt_tonemap(col);
    #endif
    
    col = pow(col, vec3(GAMMA)); // Gamma

    // LDR Post-processing
    #if POSTPROCESS
    col = CONTRAST * (col - vec3(0.5)) + vec3(0.5);  // Contrast
    col = BRIGHTNESS + col;                          // Brightness
    
    #ifdef VIGNETTE
    col *= 0.8 + 0.3*sqrt( 16.0*uv.x*uv.y*(1.0-uv.x)*(1.0-uv.y) ); // Vignette
    #endif
    
    #endif
    
    fragColor = vec4(col, 1.0);
}