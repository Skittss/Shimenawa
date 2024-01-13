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
                                                                                                    

    注連縄 (Shimenawa) by Henry A.                 (Self: https://www.shadertoy.com/view/clVyzW)
    
    --LAYOUT------------------------------------------------------------------------------------
    
      COMMON: Colour schemes, Rendering settings, SDF primitives, Intersectors, and util funcs.
       
    BUFFER B: Perlin-Worley texture atlas generation for clouds. 
              (Source: https://www.shadertoy.com/view/3sffzj)
              
    BUFFER C: Scene rendering
   
    BUFFER D: HDR Bloom pass
    
       IMAGE: LDR Post processing (Tonemap, Gamma correction, etc.)
       
    --------------------------------------------------------------------------------------------
    
*/

//========================================================
// COLOUR SCHEMES :)
//   0 = day, 1 = sunset, 2 = night
#define COLOUR_SCHEME 2
//========================================================
// RENDERING PARAMS

// Use stars on certain colour schemes
#define STARS

// Use orbit controls (debug cam - click and drag.)
//#define USE_ORBIT_CAMERA

// Increasing max steps can help render objects at oblique angles
#define RAYMARCH_MAX_STEPS 128
#define BRIDGE_LOD_DIST 800.0;

// Width of SSAA square (e.g. AA = 2 corresponds to 4x SSAA)
#define AA 1
#define AO_SAMPLES 64.0

// Cloud settings
//   Fast Raymarch
#define CAMERA_RAY_STEPS 32
#define  LIGHT_RAY_STEPS 8
//   Slow Raymarch - Less edge noise.
//#define CAMERA_RAY_STEPS 64
//#define  LIGHT_RAY_STEPS 12

// Dithering (removes artifacts from discrete ray sampling as a result of large-ish step size)
#define CLOUD_BLUE_NOISE

// Cast shadows from the background -> the clouds. Very pretty!!! but extremely slow. 
//#define CLOUD_SHADOW_CAST

//========================================================
// POST-PROCESSING PARAMS

// Toggle - (Gamma and tonemapper always applied)
#define POSTPROCESS 1

// HDR
#define EXPOSURE 1.0
#define BLOOM_INTENSITY 1.1
#define BLOOM_THRESHOLD 1.3

// Tonemappers
//#define USE_GT_TONEMAPPER
#define USE_ACES_TONEMAPPER

// GT Tonemap Params
#define GT_MAX_BRIGHTNESS 1.00
#define GT_CONTRAST 1.00
#define GT_LINEAR_OFFSET 0.22
#define GT_LINEAR_LENGTH 0.40
#define GT_BLACK_TIGHTNESS_CURVATURE 1.33
#define GT_BLACK_TIGHTNESS_OFFSET 0.00

// LDR (post-HDR effects)
#define GAMMA 0.4545
#define BRIGHTNESS 0.0
#define CONTRAST 1.15
#define TINT vec3(1.0, 1.0, 1.0)
//#define VIGNETTE

//========================================================
// MESS WITH THESE AT YOUR OWN PERIL

// Shadertoy does not have arbitrary float buffers, so we have to store HDR values in the range [0, 1].
//  Hacky solution, but we can divide through by a max HDR value before to compensate with HDR calculations.
#define HDR_MAX_COL 25.0

#define BLOOM_MIPMAP_NUDGE 0.005

//========================================================
// Util

// https://iquilezles.org/articles/smin
float smin(float a, float b, float k)
{
    float h = max(k-abs(a-b),0.0);
    return min(a, b) - h*h*0.25/k;
}

mat2 rot(in float a) { float c = cos(a); float s = sin(a); return mat2(c, s, -s, c); }

// https://www.shadertoy.com/view/3sffzj
// For Perlin-Worley noise
float saturate(float x)
{
	return clamp(x, 0.0, 1.0);
}

float remap(float x, float low1, float high1, float low2, float high2)
{
	return low2 + (x - low1) * (high2 - low2) / (high1 - low1);
}

//========================================================
// SDFs

float sdSphere( vec3 p, float s )
{
    return length(p)-s;
}

// https://iquilezles.org/articles/distfunctions
// a and b are start points, with radii provided
float sdCone( vec3 p, vec3 a, vec3 b, float ra, float rb )
{
    float rba  = rb-ra;
    float baba = dot(b-a,b-a);
    float papa = dot(p-a,p-a);
    float paba = dot(p-a,b-a)/baba;

    float x = sqrt( papa - paba*paba*baba );

    float cax = max(0.0,x-((paba<0.5)?ra:rb));
    float cay = abs(paba-0.5)-0.5;

    float k = rba*rba + baba;
    float f = clamp( (rba*(x-ra)+paba*baba)/k, 0.0, 1.0 );

    float cbx = x-ra - f*rba;
    float cby = paba - f;
    
    float s = (cbx < 0.0 && cay < 0.0) ? -1.0 : 1.0;
    
    return s*sqrt( min(cax*cax + cay*cay*baba,
                       cbx*cbx + cby*cby*baba) );
}

// https://iquilezles.org/articles/distfunctions
float sdCappedCylinder(vec3 p, float h, float r)
{
    vec2 d = abs(vec2(length(p.xz),p.y)) - vec2(r,h);
    return min(max(d.x,d.y),0.0) + length(max(d,0.0));
}

float sdInfVerticalCylinder( vec3 p, vec3 c )
{
  return length(p.xz - c.xy) - c.z;
}

// https://iquilezles.org/articles/distfunctions
float sdBox(vec3 p, vec3 b)
{
  vec3 q = abs(p) - b;
  return length(max(q,0.0)) + min(max(q.x,max(q.y,q.z)),0.0);
}

// https://iquilezles.org/articles/distfunctions
float sdBoxFrame(vec3 p, vec3 b, float e)
{
    p = abs(p  )-b;
    vec3 q = abs(p+e)-e;
    return min(min(
    length(max(vec3(p.x,q.y,q.z),0.0))+min(max(p.x,max(q.y,q.z)),0.0),
    length(max(vec3(q.x,p.y,q.z),0.0))+min(max(q.x,max(p.y,q.z)),0.0)),
    length(max(vec3(q.x,q.y,p.z),0.0))+min(max(q.x,max(q.y,p.z)),0.0));
}

// https://iquilezles.org/articles/distfunctions
float sdVerticalCapsule(vec3 p, float r, float h)
{
  p.x -= clamp( p.x, 0.0, h );
  return length( p ) - r;
}

float sdCircleXZ(in vec3 p, in float r, out float t) 
{
    float closest_p_theta = atan(p.z, p.x) + 3.14159265;
    t = closest_p_theta / (2.0 * 3.14159265);
    
    return length(vec2(p.x, p.z))-r;
}

// https://www.shadertoy.com/view/msVBzy
float sdPrism(vec3 position, float halfWidth1, float halfWidth2, float halfHeight, float halfDepth) {
    position.x = abs(position.x);
    position.x -= 0.5 * (halfWidth2 + halfWidth1);
    vec2 e = vec2(0.5 * (halfWidth2 - halfWidth1), halfHeight);
    vec2 q = position.xy - e * clamp(dot(position.xy, e) / dot(e, e), -1.0, 1.0);
    float d1 = length(q);
    if (q.x < 0.0) {
        d1 = max(-d1, abs(position.y) - halfHeight);
    }
    float d2 = abs(position.z) - halfDepth;
    return length(max(vec2(d1, d2), 0.0)) + min(max(d1, d2), 0.0);
}

// https://www.shadertoy.com/view/wsSGDG
float sdOctahedron(vec3 p, float s)
{
    p = abs(p);
    float m = p.x + p.y + p.z - s;
    vec3 r = 3.0*p - m;
    vec3 o = min(r, 0.0);
    o = max(r*2.0 - o*3.0 + (o.x+o.y+o.z), 0.0);
    return length(p - s*o/(o.x+o.y+o.z));
}

// https://www.shadertoy.com/view/Ntd3DX
// TODO: this is probably too expensive for its use case
float sdPyramid(vec3 position, float halfWidth, float halfDepth, float halfHeight) {
    position.y += halfHeight;
    position.xz = abs(position.xz);
    vec3 d1 = vec3(max(position.x - halfWidth, 0.0), position.y, max(position.z - halfDepth, 0.0));
    vec3 n1 = vec3(0.0, halfDepth, 2.0 * halfHeight);
    float k1 = dot(n1, n1);
    float h1 = dot(position - vec3(halfWidth, 0.0, halfDepth), n1) / k1;
    vec3 n2 = vec3(k1, 2.0 * halfHeight * halfWidth, -halfDepth * halfWidth);
    float m1 = dot(position - vec3(halfWidth, 0.0, halfDepth), n2) / dot(n2, n2);
    vec3 d2 = position - clamp(position - n1 * h1 - n2 * max(m1, 0.0), vec3(0.0), vec3(halfWidth, 2.0 * halfHeight, halfDepth));
    vec3 n3 = vec3(2.0 * halfHeight, halfWidth, 0.0);
    float k2 = dot(n3, n3);
    float h2 = dot(position - vec3(halfWidth, 0.0, halfDepth), n3) / k2;
    vec3 n4 = vec3(-halfWidth * halfDepth, 2.0 * halfHeight * halfDepth, k2);
    float m2 = dot(position - vec3(halfWidth, 0.0, halfDepth), n4) / dot(n4, n4);    
    vec3 d3 = position - clamp(position - n3 * h2 - n4 * max(m2, 0.0), vec3(0.0), vec3(halfWidth, 2.0 * halfHeight, halfDepth));
    float d = sqrt(min(min(dot(d1, d1), dot(d2, d2)), dot(d3, d3)));
    return max(max(h1, h2), -position.y) < 0.0 ? -d : d;
}

// https://iquilezles.org/articles/distfunctions
float sdLink(in vec3 p, in float le, in float r1, in float r2)
{
    vec3 q = vec3( p.x, max(abs(p.y)-le,0.0), p.z );
    return length(vec2(length(q.xy)-r1,q.z)) - r2;
}

// https://iquilezles.org/articles/distfunctions
float sdTorus(in vec3 p, in vec2 t)
{
    vec2 q = vec2(length(p.xz)-t.x,p.y);
    return length(q)-t.y;
}

//========================================================
// Intersectors

// https://iquilezles.org/articles/intersectors
vec2 iSphere(in vec3 ro, in vec3 rd, in float rad)
{
	float b = dot( ro, rd );
	float c = dot( ro, ro ) - rad*rad;
	float h = b*b - c;
	if( h<0.0 ) return vec2(-1.0); // no intersection
    h = sqrt(h);
	return vec2( -b-h, -b+h );
}

// https://www.shadertoy.com/view/3sffzj
vec2 intersectAABB(vec3 rayOrigin, vec3 rayDir, vec3 boxMin, vec3 boxMax) 
{
    vec3 tMin = (boxMin - rayOrigin) / rayDir;
    vec3 tMax = (boxMax - rayOrigin) / rayDir;
    vec3 t1 = min(tMin, tMax);
    vec3 t2 = max(tMin, tMax);
    float tNear = max(max(t1.x, t1.y), t1.z);
    float tFar = min(min(t2.x, t2.y), t2.z);
    return vec2(tNear, tFar);
}

bool insideAABB(vec3 p, vec3 min_corner, vec3 max_corner)
{
    float eps = 1e-4;
	return (p.x > min_corner.x-eps) && (p.y > min_corner.y-eps) && (p.z > min_corner.z-eps) && 
		   (p.x < max_corner.x+eps) && (p.y < max_corner.y+eps) && (p.z < max_corner.z+eps);
}
