#define PHI 1.61803399
#define PI  3.14159265
#define TAU 6.28318533

#define USE_DEBUG_CAMERA 1

#define ZERO (min(iFrame,0))

// Sky

// Rope Params
const vec3 _ShideWindParams = vec3(0.10, 4.0, 1.75); // Wave amplitude, distance modifier (stiffness), anim speed
const vec3 _ShideWindParams_s = vec3(0.013, 72.0, 3.5);
const vec2 _KiraretanawaWindParamsYZ = vec2(PI / 20.0, 2.3); // Max rot, anim speed
const vec2 _KiraretanawaWindParamsYX = vec2(PI / 30.0, 2.0);

// Materials
#define MAT_ROPE 1.0
#define MAT_SHIDE 2.0
#define MAT_SHIDE_SECONDARY 3.0
const vec3 _MatRope = vec3(0.95, 0.89, 0.74);
//const vec3 _RopeTerminatorLineCol = 0.8*vec3(0.49, 0.329, 1.0);
const vec3 _RopeTerminatorLineCol = 0.8 * vec3(1.0, 0.329, 0.518);
const vec3 _MatShide = vec3(1.0, 1.0, 1.0);
const vec3 _MatShideSecondary = 0.8*vec3(1.0, 0.412, 0.412);

// Illumination
const vec3  _SunPos  = vec3(30.0, 20.0, 30.0);
//const vec3 _SunPos = vec3(0.2, 56, -40.1);

const float _SunSize = 3.5;
const float _SunBrightness = 1.4;
const vec3  _SunCol  = vec3(0.98, 0.72, 0.31);
//const vec3  _SunCol  = vec3(0.51471, 0.79919, 1.0);

const vec3 _ZenithCol = 0.5 * vec3(0.37, 0.14, 0.43);
const vec3 _HorizonCol = 1.5 * vec3(0.36, 0.17, 0.66);
//const vec3 _HorizonCol = vec3(0.0);
const vec3 _NadirCol = vec3(0.35, 0.32, 0.8);

const float _ZenithAttenuation = 1.8;
const float _NadirAttenuation  = 1.2;
const float _HorizonOffset = 0.0;

const vec3 _LightDir = normalize(_SunPos);
//const vec3 _LightDir = normalize(vec3(2.0, 1.0, 2.0));

const float _RopeExtraShadowBrightness = 0.025; // Set to zero for flat shaded look

// SSS Sun Outline
const vec3  _OutlineCol = _SunCol;
const float _OutlineDotThreshold = 0.95;
const float _OutlineRadialAttenuation = 2.0;
const float _OutlineMaxDist = 0.05;
const float _OutlineAttenuation = 16.0;
const float _OutlineExtraThickness = 0.0035;


// Util
#define CMP_MAT_LT(a, b) a < (b + 0.5)
#define CMP_MAT_GT(a, b) a > (b - 0.5)

//==HELPER FUNC===================================================================================================================================

vec2 MIN_MAT(vec2 a, vec2 b)
{
    return (a.x < b.x) ? a : b;
}

//util functions for AO
// https://www.shadertoy.com/view/ld3Gz2
float hash1(float n)
{
    return fract(sin(n)*43758.5453123);
}

vec3 forwardSF(float i, float n) 
{
    float phi = 2.0*PI*fract(i/PHI);
    float zi = 1.0 - (2.0*i+1.0)/n;
    float sinTheta = sqrt( 1.0 - zi*zi);
    return vec3( cos(phi)*sinTheta, sin(phi)*sinTheta, zi);
}

//==SDF===========================================================================================================================================

float sdShide( in vec3 p, in int s_n, in float sec_id, out float seg_id )
{
    // 紙垂 multiple zig-zag boxes
    // 紙で作られているから、ちょっとSSS必要あり。
    
    
    // 全物体を４５C回る
    p.yz = (p.yz + vec2(p.z, -p.y))*sqrt(0.5); // Shortcut for 45-degrees rotation

    // 大きさ
    vec3 dim = 0.02 * vec3(0.03, 1.0, 2.0); // 2:1:0.1 ratio
    vec3 smallDim = dim;
    smallDim.y /= 2.0;
    smallDim.z *= 0.9;
    
    //vec3 shidePos = p - vec3(0.0, 1.5 * dim.y, 0.0);
        
    //shidePos += 0.1* sin(shidePos * 3.0 + iTime) * distToConnector;
       
    // 風から動きの準備
    // 垂直動きが欲しいので、回りの後でする。
    // 接続点から距離を使って動きの強さが決まる
    vec3 connectorPos = vec3(0., 2.0 * dim.y, -1.5 * dim.z);
    float distToConnector = length(connectorPos - p);

    // Thinner box, mount to rope    
    
    vec3 thinShidePos = p + vec3(.0, -1.5 * dim.y, 0.9 * dim.z);
    // 動かす - Sinusoidal motion
    thinShidePos.x += _ShideWindParams.x * sin(distToConnector * _ShideWindParams.y + _ShideWindParams.z * iTime + 53.0 * sec_id) * distToConnector;
    float d = sdBox(thinShidePos, smallDim);
    
    // Zig-zag - domain repetition is a viable option here but I'm not sure how to
    //   ensure the sdf remains correct when applying a stepwise shift in the domain
    
    seg_id = 0.0;
    for (int i=0; i<s_n; i++)
    {
        vec3 shidePos = p + vec3(-5.0 * dim.x * mod(float(i+1), 2.0), 2.0 * dim.y * float(i), -dim.z * float(i));
        shidePos.x += _ShideWindParams.x * sin(distToConnector * _ShideWindParams.y + _ShideWindParams.z * iTime + 53.0 * sec_id) * distToConnector;
        shidePos.x += _ShideWindParams_s.x * sin(distToConnector * _ShideWindParams_s.y + _ShideWindParams_s.z * iTime + 13.0 * sec_id) * distToConnector;
        //d = min(d, sdBox(shidePos, dim));
        float d2 = sdBox(shidePos, dim);
        if (d2 < d) {d = d2; seg_id = float(i+1);}
    }
    
    return d;
}

// p is position of top, s is scale
float sdKiraretanawa( in vec3 p, in float s, in float connectorOffset, in float id )
{
    // Offset rotation to make object sway
    p.yz *= rot(_KiraretanawaWindParamsYZ.x * sin(_KiraretanawaWindParamsYZ.y * 0.89 * iTime + 71.0 * id));
    p.yx *= rot(_KiraretanawaWindParamsYX.x * cos(_KiraretanawaWindParamsYX.y * 0.57 * iTime + 31.0 * id));
    p.y += connectorOffset;

    // TODO: is diving multiplying d by s at the end equivalent to this scaling??
        
    float d = sdCone(p + s * vec3(0., 0.6, 0.), vec3(0.), s * vec3(0.,0.5,0.), s * 0.32, s * 0.1) - s * 0.05;
    d = min(d, sdCappedCylinder(p, s * 0.06, s * 0.12) - s * 0.015);
    
    {
    // 縛り縄 - Binding rope, 3 for cost of 1
    vec3 q = p;
    q.y = abs(p.y + s * 0.02) - s * 0.0165; // XZ plane symmetry for duplication, and translate downwards
    q.y = abs(q.y) - s * 0.0165;
    d = min(d, sdTorus(q, vec2(s * 0.15, s * 0.013)));
    }
    
    // 接続縄 - Connecting rope
    d = min(d, sdCappedCylinder(p - vec3(0.0, connectorOffset/2.0, 0.0), connectorOffset, s * 0.013));
    
    return d;
}


float sdSwirl(vec3 p, in float r, in float c, in float f, out float id)
{
    // Rotate and twist a capsule to make the rope along the x axis
    float l = TAU * r; // circumference rope length
    p.yz*=rot(p.x*PI*f); // twisting; 45-deg per x
    //p.yz = abs(p.yz)-0.25; // 4 for the price of one
    //p.yz = (p.yz + vec2(p.z, -p.y))*sqrt(0.5); // Shortcut for 45-degrees rotation
    vec2 sector = step(0.0,p.yz);
    id = sector.x + 2.0 * sector.y;
    p.y = abs(p.y) - 0.02;
    //p.yz = abs(p.yz)-0.05;
    //p.yz*=rot(p.x*PI*4.0);
    //p.yz = abs(p.yz)-0.02;
    float d = sdVerticalCapsule(p,c,l); //Potentially better solution exists with capped torus
    return d;
}

vec2 sdShimenawa(in vec3 p, in float r, in float c, in float f, out float id)
{
    vec2 res = vec2(1e10);

    float ring_t = 0.0;
    float dRing = sdCircleXZ(p, r, ring_t);
    ring_t = 1.0 - ring_t; // reflect rope direction
    
    //p.y = p.y + 0.06*sin(2.0*2.0*PI*ring_t+1.5*time);
    
    #if 1
    if (dRing < 0.8) // approx bounding ring
    {    
        // I'm fairly confident this local-UV transformation majorly messes up the SDF for isolines > 0.0 but... oh well.
        //vec3 ring_uv = vec3(ring_t * TAU * r, p.y + 0.06*sin(2.0*2.0*PI*ring_t+1.5*time), dRing);
        vec3 ring_uv = vec3(ring_t * TAU * r, p.y, dRing);
        res = vec2(sdSwirl(ring_uv, r, c, f, id), MAT_ROPE);
    }
    #endif
    
    #if 1
    {
    vec3 q = p;
    const float an = TAU/7.0;
    float sector = round(atan(p.z,p.x)/an);
    float angrot = sector*an;
    q.xz *= rot(angrot);
    float seg_id;
    float d = sdShide(q - vec3(r + 0.4*c, -2.3*c, 0.0), 4, sector+1.0, seg_id);
    float shide_mat = (mod(seg_id, 2.0) == 0.0) ? MAT_SHIDE : MAT_SHIDE + 1.0;
    res = MIN_MAT(res, vec2(d, shide_mat));
    }
    #endif
    
    #if 1
    {
    // I think this rotation could be done in one go, but my brain is a mess thinking about the domain repetition here -.-
    vec3 q_s = p;
    
    // Prerotate offset
    float an = TAU/14.0;
    float co = cos(an), si = sin(an);
    q_s.xz *= rot(an);
    
    an = TAU/7.0;
    float sector = round(atan(q_s.z,q_s.x)/an);
    float angrot = sector*an;
    q_s.xz *= rot(angrot);
    float d = sdKiraretanawa(q_s - vec3(r, -1.4*c, 0.0), 0.14, 0.05, sector+1.0);
    res = MIN_MAT(res, vec2(d, MAT_ROPE));
    }
    #endif
    
    return res;
}

float sdBark( in vec3 p, in float h, in float r, in float d, in float w, in float n, in float phase)
{
    float angle = TAU / n;
    float sector = round(atan(p.z, p.x)/angle);
    vec3 q = p;
    float an = sector * angle + phase;
    q.xz = mat2(cos(an), -sin(an),
                sin(an), cos(an)) * q.xz;
                
    return sdBox(q-vec3(r, 0.0, 0.01*cos(50.0*q.y)), vec3(w, h, d)) - 0.001;
}

vec2 sdTree( in vec3 p, in float h, in float r )
{
    float m = -1.0;

    float db = sdBark(p, h, r, 0.0025, 0.005, 24.0, 0.0);
    db = min(db, sdBark(p, h, r, 0.001, 0.005, 12.0, 0.17));
    float dTree = sdCappedCylinder(p, h, r) - 0.001;
    
    return vec2(min(dTree, db), m);
}

//==ILLUMINATION================================================================================================================================

vec3 sky( in vec3 ro, in vec3 rd ) 
{
    float sunDist = length(_SunPos);
    float dist = length(_SunPos - (ro + rd * sunDist)) - _SunSize;
    dist = dist < 0.0 ? 0.0 : dist;
    
    float ry = _HorizonOffset + rd.y;
    
    //float haloDist = clamp((dist - 4.0) * 0.015, 0.0, 1.0);
    //float halo = pow(haloDist, 0.5);
    //vec3 sun = halo * _SunCol;
    
    // Atmosphere
    float zenith  = 1.0 - pow(min(1.0, 1.0 - ry), _ZenithAttenuation);
    float nadir   = 1.0 - pow(min(1.0, 1.0 + ry), _NadirAttenuation);
    float horizon = 1.0 - zenith - nadir;
    
    vec3 skycol = zenith * _ZenithCol + nadir * _NadirCol + horizon * _HorizonCol;
    
    // Clouds
    //  Generate 3 cloud levels - Close, mid, far, which change the shape of the clouds
    
    //  Also add a few wispy streaks for extra detail.
    
    
    float halo = pow((3.0/dist), 0.9);
    vec3 sun = halo * _SunCol;
    skycol = 1.0 - exp(-(skycol + sun));
    skycol *= _SunBrightness;
    
    return skycol;
}

//==RENDERING===================================================================================================================================

vec2 map( in vec3 p )
{    
    vec2 res = vec2(1e10); // (Distance, Material)
    
    #if 0
    res = sdTree(p, 0.40, 0.40);
    #endif
    
    #if 1
    float r_id;
    res = MIN_MAT(res, sdShimenawa(p, 0.4615, 0.04, 10.0, r_id));
    #endif
           
    return res; // returns (distance, material) pair.
}

// https://iquilezles.org/articles/normalsSDF
vec3 calcNormal( in vec3 pos )
{
#if 0
    vec2 e = vec2(1.0,-1.0)*0.5773;
    const float eps = 0.00025;
    return normalize( e.xyy*map( pos + e.xyy*eps ).x + 
					  e.yyx*map( pos + e.yyx*eps ).x + 
					  e.yxy*map( pos + e.yxy*eps ).x + 
					  e.xxx*map( pos + e.xxx*eps ).x );
#else
    // klems's trick to prevent the compiler from inlining map() 4 times
    vec3 n = vec3(0.0);
    for( int i=ZERO; i<4; i++ )
    {
        vec3 e = 0.5773*(2.0*vec3((((i+3)>>1)&1),((i>>1)&1),(i&1))-1.0);
        n += e*map(pos+0.0005*e).x;
    }
    return normalize(n);
#endif    
}

float calcSSS( in vec3 pos, in vec3 nor )
{
	float occ = 0.0;
    for( int i=ZERO; i<8; i++ )
    {
        float h = 0.002 + 0.11*float(i)/7.0;
        vec3 dir = normalize( sin( float(i)*13.0 + vec3(0.0,2.1,4.2) ) );
        dir *= sign(dot(dir,nor));
        occ += (h-map(pos-h*dir).x);
    }
    occ = clamp( 1.0 - 11.0*occ/8.0, 0.0, 1.0 );    
    return occ*occ;
}

// https://iquilezles.org/articles/rmshadows
// Soft Shadows with no backtracking
// Generally slower and less accurate than with backtracking, but with more predictable behavior for distorted interior SDFs
float softShadow( in vec3 ro, in vec3 rd, float mint, float maxt, float w )
{
	float res = 1.0;
    float t = mint;
    float ph = 1e10; // big, such that y = 0 on the first iteration
    
    for( int i=0; i<32; i++ )
    {
        // TODO: Mapping of big rope is problematic, separate out conservative step.
		float h = map( ro + rd*t ).x;
        float y = h*h/(2.0*ph);
        float d = sqrt(h*h-y*y);
        res = min( res, d/(w*max(0.0,t-y)) );
        ph = h;
        
        // TODO: Apply this conservative step ONLY to the rope
        t += 0.75*h;
        
        if( res < 0.0001 || t > maxt ) break;
        
    }
    res = clamp( res, 0.0, 1.0 );
    return res*res*(3.0-2.0*res);
}

// https://iquilezles.org/articles/rmshadows
// Soft Shadows with backtracking.
float softShadowBacktrack( in vec3 ro, in vec3 rd, float k )
{
    float res = 1.0;
    float t = 0.01;
    for( int i=ZERO; i<32; i++ )
    {
        float h = map(ro + rd*t).x;
        res = min( res, smoothstep(0.0,1.0,k*h/t) );
        //t += clamp( h, 0.0, 0.1 );
        t += h;
		//if( res<0.01 ) break;
		if( abs(res)<0.01 ) break;
    }
    return clamp(res,0.0,1.0);
}

// https://www.shadertoy.com/view/ld3Gz2
//  AO by (pseudo-random) sampling a number of distances to surfaces in a hemisphere about the normal
float calcAO( in vec3 pos, in vec3 nor )
{
	float ao = 0.0;
    for( int i=ZERO; i<int(AO_SAMPLES); i++ )
    {
        vec3 ap = forwardSF( float(i), AO_SAMPLES );
        float h = hash1(float(i));
		ap *= sign( dot(ap,nor) ) * h*0.1;
        ao += clamp( map(pos + nor*0.01 + ap).x*3.0, 0.0, 1.0 );
    }
	ao /= AO_SAMPLES;
	
    return clamp( ao*6.0, 0.0, 1.0 );
}

vec3 sunSSSOutline( in vec3 ro, in vec3 rd, in float d )
{
    // Create an outline around objects within a certain radius around the sun.
    //    The view direction dot product can be used to determine the radius.
    //    Raymarching near missesare used to generate the outline. 
    //   Then, a power curve is used to decay the outline based on distance from the sun.
    float rd_dot_sun = dot(rd, normalize(_SunPos - ro));
    
    float angle_factor = 1.0 / (1.0 - _OutlineDotThreshold) * (max(0.0, rd_dot_sun - _OutlineDotThreshold));
    angle_factor = pow(angle_factor, _OutlineRadialAttenuation);
    
    // It might be faster to check the dot threshold here when computing the outline if the attenuation power curve is expensive.
    float a = step(0.001, d); // Stencil for non-edges
    float b = clamp((d - _OutlineExtraThickness) / _OutlineMaxDist, 0.0, 1.0);
    
    return a * angle_factor * _OutlineCol * pow(1.0 - b, _OutlineAttenuation);
}

vec3 intersect( in vec3 ro, in vec3 rd )
{
    vec3 res = vec3(-1.0, 1e10, -1.0);
    
    // bounding sphere
    vec2 tminmax = iSphere( ro, rd, 1.0 );
	if( tminmax.y>0.0 )
    {
        // raymarch
        float t = max(tminmax.x,0.001);
        for( int i=0; i<128 && t<tminmax.y; i++ )
        {
            vec2 h = map(ro+t*rd);
            res.y = max(min(res.y, h.x), 0.0); // Track near misses for strong light outine
            if( abs(h.x)<0.001 ) { res=vec3(t, res.y, h.y); break; }
            t += 0.8 * h.x; // Coeff here is to try and avoid overshooting at the cost of performance
                            //  TODO: There should be a much smarter solution than this. The problem only arises with large domain distortions.
        }
    }
    
    return res; // t, nearest, mat_id
}

vec3 shade( in vec3 ro, in vec3 rd, in float t, in float m ) 
{
    vec3 pos = ro + t*rd;
    vec3 nor = calcNormal(pos);
    float shadow = softShadowBacktrack(pos - 0.01*rd, _LightDir, 2.0);
    float occ = calcAO(pos, nor);
    shadow = pow(occ, 2.0) * (shadow + _RopeExtraShadowBrightness) / (1.0 + _RopeExtraShadowBrightness);
    //shadow = (shadow + occ) / 2.0;
    //float shadow = softShadow(pos - 0.01*rd, _LightDir, 0.002, 1.0, 0.4);

    if (CMP_MAT_LT(m, MAT_ROPE)) 
    {
        vec3 base_shadow = mix(0.6*_MatRope, _HorizonCol, 0.2);
        // TODO: I think this multiplier of the shadow coeff changes the base shadow colour too.
        vec3 sss_style_mix = mix(base_shadow, _RopeTerminatorLineCol, min(1.0, 4.0 * shadow));
        //vec3 sss_style_mix = mix(base_shadow, _RopeTerminatorLineCol, shadow);
        
        vec3 albedo = mix(sss_style_mix, _MatRope, min(1.0, 2.0 * shadow));
        //vec3 albedo = mix(sss_style_mix, _MatRope, shadow);
               
        return albedo;
    }
    else if (CMP_MAT_LT(m, MAT_SHIDE_SECONDARY)) // handle both shide mat variatons here
    {
        vec3 mat = (m == MAT_SHIDE) ? _MatShide : _MatShideSecondary;
        // This simple SSS approximation is good enough for sun -> paper.
        float tr_range = t / 5.0;
        float view_bias = abs(dot(normalize(ro - pos), normalize(pos - _SunPos)));
        float sun_transmission = map(pos + _LightDir * tr_range).x / tr_range;
        vec3 sss = 0.3*_SunCol * smoothstep(0.0, 1.0, sun_transmission);
        
        vec3 base_shadow = mix(0.6*mat, _HorizonCol, 0.2);

        //sss = sss + fre + (0.5+0.5*fre)*pow(abs(t-0.2),1.0);
        
        //return vec3(occ);
        return sss + mix(base_shadow, mat, shadow);
    }
    
    vec3 col = vec3(0.0);
    col = 0.5 + 0.5*nor;
    col = mix(vec3(0.0), col, shadow);

    return col;
}

vec3 render ( in vec3 ro, in vec3 rd ) 
{
    // background
    //vec3 col = vec3(1.0+rd.y)*0.03;
    vec3 col = sky(ro, rd);

    // raymarch geometry
    vec3 tm = intersect( ro, rd );
    if( tm.x>0.0 )
    {
        col = shade(ro, rd, tm.x, tm.z);
    }
    
    // TODO: This should only be applied to certain materials that SSS (and to varying extents?)
    col += sunSSSOutline(ro, rd, tm.y);
    
    return col;
}

mat3 setCamera( in vec3 ro, in vec3 ta, float cr )
{
	vec3 cw = normalize(ta-ro);
	vec3 cp = vec3(sin(cr), cos(cr),0.0);
	vec3 cu = normalize( cross(cw,cp) );
	vec3 cv =          ( cross(cu,cw) );
    return mat3( cu, cv, cw );
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec3 tot = vec3(0.0);

    #if AA>1
    for( int m=ZERO; m<AA; m++ )
    for( int n=ZERO; n<AA; n++ )
    {
        // pixel coordinates
        vec2 o = vec2(float(m),float(n)) / float(AA) - 0.5;
        vec2 p = (2.0*(fragCoord+o)-iResolution.xy)/iResolution.y;
        float d = 0.5*sin(fragCoord.x*147.0)*sin(fragCoord.y*131.0);
        #else    
        vec2 p = (2.0*fragCoord-iResolution.xy)/iResolution.y;
        #endif

	    // camera & movement
        float an = TAU*1.4*iTime/40.0;
        vec3 ta = vec3( 0.0, 0.0, 0.0 );
        
        vec2 m = iMouse.xy / iResolution.xy-.5;
        vec3 ro;
        #if USE_DEBUG_CAMERA
        if ( iMouse.x >= 0.0 ) 
        {
            vec3 target = vec3(0.0);
            ro = vec3(0.8, 0.3, 0.8);
            ro.yz *= rot(m.y*PI);
            ro.xz *= rot(-m.x*TAU*2.0);
        }
        #else
        //ro = ta + vec3( 1.5*cos(an), 0.3, 1.5*sin(an) );
        ro = ta + vec3( cos(an), -0.3, sin(an) );
        //ro = ta + vec3( 0.001, 2.0, 0.0 );
        #endif
        
        // camera-to-world transformation
        mat3 ca = setCamera( ro, ta, 0.0 );
        
        // ray direction
        float fl = 2.0;
        vec3 rd = ca * normalize( vec3(p,fl) );
        
        vec3 col = render( ro, rd );
        
        // gamma        
	    tot += col;
    #if AA>1
    }
    tot /= float(AA*AA);
    #endif
    
    tot *= EXPOSURE;

    // cheap dithering
    tot += sin(fragCoord.x*114.0)*sin(fragCoord.y*211.1)/512.0;
    tot /= HDR_MAX_COL;
    
    fragColor = vec4( tot, 1.0 );
}