// Shadertoy does not support the use of proper float buffers, so we have to store HDR values in the range [0, 1].
//  Hacky solution, but we can divide through by a max HDR value before to compensate with HDR calculations.
#define HDR_MAX_COL 25.0

//========================================================
// POST-PROCESSING PARAMS

// Toggle - Gamma applied either way
#define POSTPROCESS 1

// HDR
#define EXPOSURE 1.0
#define BLOOM_INTENSITY 2.0
#define BLOOM_THRESHOLD 1.3

// LDR
#define GAMMA 2.2
#define BRIGHTNESS 0.0
#define CONTRAST 1.15

// GT TONEMAP
#define GT_MAX_BRIGHTNESS 1.00
#define GT_CONTRAST 1.00
#define GT_LINEAR_OFFSET 0.22
#define GT_LINEAR_LENGTH 0.40
#define GT_BLACK_TIGHTNESS_CURVATURE 1.33
#define GT_BLACK_TIGHTNESS_OFFSET 0.00

// RENDERING PARAMS
#define AA 1
#define AO_SAMPLES 64.0

// MESS WITH THESE AT YOUR OWN PERIL
#define BLOOM_MIPMAP_NUDGE 0.005

//========================================================
// Util

// https://iquilezles.org/articles/smin
float smin( float a, float b, float k )
{
    float h = max(k-abs(a-b),0.0);
    return min(a, b) - h*h*0.25/k;
}

//========================================================
// SDFs

// https://iquilezles.org/articles/distfunctions
// a and b are start points, with radii provided
float sdCone(vec3 p, vec3 a, vec3 b, float ra, float rb)
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
float sdCappedCylinder( vec3 p, float h, float r )
{
    vec2 d = abs(vec2(length(p.xz),p.y)) - vec2(r,h);
    return min(max(d.x,d.y),0.0) + length(max(d,0.0));
}

// https://iquilezles.org/articles/distfunctions
float sdBox( vec3 p, vec3 b )
{
  vec3 q = abs(p) - b;
  return length(max(q,0.0)) + min(max(q.x,max(q.y,q.z)),0.0);
}

// https://iquilezles.org/articles/distfunctions
float sdVerticalCapsule( vec3 p, float r, float h )
{
  p.x -= clamp( p.x, 0.0, h );
  return length( p ) - r;
}

float sdVertCylinderInf( in vec3 p, in float r )
{
    return length(p.xz)-r;
}

float sdCircleXZ( in vec3 p, in float r, out float t) 
{
    float closest_p_theta = atan(p.z, p.x) + 3.14159265;
    t = closest_p_theta / (2.0 * 3.14159265);
    
    return length(vec2(p.x, p.z))-r;
}

// t-valued ring defined in xz so that periodic functions can be applied around the tree
vec2 sdTorusRing( in vec3 p, in float r, in float c)
{
    vec3 ring_origin = p - vec3(r, 0.0, 0.0); // use point in x axis as the origin
    vec2 p_xz = normalize(p).xz;
    vec3 closest_p = r * vec3(p_xz.x, 0.0, p_xz.y);
    float t = dot(normalize(ring_origin.xz), p_xz) / 6.283185;
    
    vec2 q = vec2(length(p.xz)-r,p.y);
    
    return vec2(t, length(q)-c);
}

// https://iquilezles.org/articles/distfunctions
float sdTorus( in vec3 p, in vec2 t )
{
    vec2 q = vec2(length(p.xz)-t.x,p.y);
    return length(q)-t.y;
}

mat2 rot(in float a) { float c = cos(a); float s = sin(a); return mat2(c, s, -s, c); }

//========================================================
// Intersectors

// https://iquilezles.org/articles/intersectors
vec2 iSphere( in vec3 ro, in vec3 rd, in float rad )
{
	float b = dot( ro, rd );
	float c = dot( ro, ro ) - rad*rad;
	float h = b*b - c;
	if( h<0.0 ) return vec2(-1.0);
    h = sqrt(h);
	return vec2(-b-h, -b+h );
}
