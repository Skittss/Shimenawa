#define AA 2
#define PI 3.14159265
#define TAU 6.2831853

#define USE_DEBUG_CAMERA 1

const vec3 _ShideWindParams = vec3(0.10, 4.0, 1.75);
const vec3 _ShideWindParams_s = vec3(0.013, 72.0, 3.5);

float layeredPerlin1D()
{
    return 1.0;
}

float sdShide( in vec3 p, in int s_n, in float id )
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
    thinShidePos.x += _ShideWindParams.x * sin(distToConnector * _ShideWindParams.y + _ShideWindParams.z * iTime + 53.0 * id) * distToConnector;
    float d = sdBox(thinShidePos, smallDim);
    
    // Zig-zag - domain repetition is a viable option here but I'm not sure how to
    //   ensure the sdf remains correct when applying a stepwise shift in the domain
    
    for (int i=0; i<s_n; i++)
    {
        vec3 shidePos = p + vec3(-2.0 * dim.x * mod(float(i+1), 2.0), 2.0 * dim.y * float(i), -dim.z * float(i));
        shidePos.x += _ShideWindParams.x * sin(distToConnector * _ShideWindParams.y + _ShideWindParams.z * iTime + 53.0 * id) * distToConnector;
        shidePos.x += _ShideWindParams_s.x * sin(distToConnector * _ShideWindParams_s.y + _ShideWindParams_s.z * iTime + 13.0 * id) * distToConnector;
        d = min(d, sdBox(shidePos, dim));
    }
    
    return d;
}

// p is position of top, s is scale
float sdKiraretanawa( in vec3 p, in float s )
{
    // TODO: is diving multiplying d by s at the end equivalent of this scaling??
    
    float d = sdCone(p + s * vec3(0., 0.6, 0.), vec3(0.), s * vec3(0.,0.5,0.), s * 0.32, s * 0.1) - s * 0.05;
    d = min(d, sdCappedCylinder(p, s * 0.06, s * 0.12) - s * 0.015);
    
    {
    // 縛り縄 - Binding rope, 3 for cost of 1
    vec3 q = p;
    q.y = abs(p.y + s * 0.02) - s * 0.0165; // XZ plane symmetry for duplication, and translate downwards
    q.y = abs(q.y) - s * 0.0165;
    d = min(d, sdTorus(q, vec2(s * 0.15, s * 0.013)));
    }
    
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

float sdShimenawa(in vec3 p, in float r, in float c, in float f, in float time, out float id)
{
    float ring_t = 0.0;
    float dRing = sdCircleXZ(p, r, ring_t);
    
    //p.y = p.y + 0.06*sin(2.0*2.0*PI*ring_t+1.5*time);
    
    float d = dRing;
    // I'm fairly confident this local-UV transformation majorly messes up the SDF for isolines > 0.0 but... oh well.
    if (dRing < 0.8) // approx bounding ring
    {    
        //vec3 ring_uv = vec3(ring_t * TAU * r, p.y + 0.06*sin(2.0*2.0*PI*ring_t+1.5*time), dRing);
        vec3 ring_uv = vec3(ring_t * TAU * r, p.y, dRing);
        d = sdSwirl(ring_uv, r, c, f, id);
    }
    
    {
    vec3 q = p;
    const float an = TAU/7.0;
    float sector = round(atan(p.z,p.x)/an);
    float angrot = sector*an;
    q.xz *= rot(angrot);
    d = min(d, sdShide(q - vec3(r + c, -2.*c, 0.0), 4, sector+1.0));
    }
    
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
    d = min(d, sdKiraretanawa(q_s - vec3(r, -1.8*c, 0.0), 0.14));
    }
    
    return d;
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

float sdTree( in vec3 p, in float h, in float r)
{
    float db = sdBark(p, h, r, 0.0025, 0.005, 24.0, 0.0);
    db = min(db, sdBark(p, h, r, 0.001, 0.005, 12.0, 0.17));
    float dTree = sdCappedCylinder(p, h, r) - 0.001;
    
    return min(dTree, db);
}

vec4 map( in vec3 p, float time )
{
    float t = sdTree(p, 0.40, 0.40);
    float r_id;
    //float r = sdShimenawa(p + vec3(.0, .1, .0), 0.4615, 0.04, 10.0, time, r_id);
    float r = sdShimenawa(p, 0.4615, 0.04, 10.0, time, r_id);
        
    float d = min(t, r);
    
    return vec4(r, p);
}

#define ZERO min(iFrame,0)

// https://iquilezles.org/articles/normalsSDF
vec3 calcNormal( in vec3 pos, in float time )
{
#if 0
    vec2 e = vec2(1.0,-1.0)*0.5773;
    const float eps = 0.00025;
    return normalize( e.xyy*map( pos + e.xyy*eps, time ).x + 
					  e.yyx*map( pos + e.yyx*eps, time ).x + 
					  e.yxy*map( pos + e.yxy*eps, time ).x + 
					  e.xxx*map( pos + e.xxx*eps, time ).x );
#else
    // klems's trick to prevent the compiler from inlining map() 4 times
    vec3 n = vec3(0.0);
    for( int i=ZERO; i<4; i++ )
    {
        vec3 e = 0.5773*(2.0*vec3((((i+3)>>1)&1),((i>>1)&1),(i&1))-1.0);
        n += e*map(pos+0.0005*e,time).x;
    }
    return normalize(n);
#endif    
}

vec4 intersect( in vec3 ro, in vec3 rd, in float time )
{
    vec4 res = vec4(-1.0);
    
    // bounding sphere
    vec2 tminmax = iSphere( ro, rd, 1000.0 );
	if( tminmax.y>0.0 )
    {
        // raymarch
        float t = max(tminmax.x,0.001);
        for( int i=0; i<128 && t<tminmax.y; i++ )
        {
            vec4 h = map(ro+t*rd,time);
            if( h.x<0.001 ) { res=vec4(t,h.yzw); break; }
            t += 0.9 * h.x; // Coeff here is to try and avoid overshooting at the cost of performance
        }
    }
    
    return res;
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
        float time = iTime;
        #else    
        vec2 p = (2.0*fragCoord-iResolution.xy)/iResolution.y;
        float time = 1.0*iTime;
        #endif

	    // camera & movement
        float an = TAU*time/40.0;
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
        ro = ta + vec3( 0.56*cos(an), -0.3, 0.56*sin(an) );
        //ro = ta + vec3( 0.001, 2.0, 0.0 );
        #endif
        
        // camera-to-world transformation
        mat3 ca = setCamera( ro, ta, 0.0 );
        
        // ray direction
        float fl = 2.0;
        vec3 rd = ca * normalize( vec3(p,fl) );

        // background
        vec3 col = vec3(1.0+rd.y)*0.03;
        
        // raymarch geometry
        vec4 tuvw = intersect( ro, rd, time );
        if( tuvw.x>0.0 )
        {
            // shading/lighting	
            vec3 pos = ro + tuvw.x*rd;
            vec3 nor = calcNormal(pos, time);
                        
            col = 0.5 + 0.5*nor;
        }
        
        
        // gamma        
	    tot += pow(col,vec3(0.45) );
    #if AA>1
    }
    tot /= float(AA*AA);
    #endif

    // cheap dithering
    tot += sin(fragCoord.x*114.0)*sin(fragCoord.y*211.1)/512.0;

    fragColor = vec4( tot, 1.0 );
}