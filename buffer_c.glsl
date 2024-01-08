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
    
    This shader has served as my learning playground for raymarching.
      It is by no means perfect, but I hope you enjoy it as much as I enjoyed making it.
    このシェーダーを作ることでレイマーチングの基本を学びました。
      決して完璧じゃありませんが、見て楽しんでください。
    
    Also, apologies for the long compile time! I tried to make it as quick as possible.
    長いコンパイル時間にすみません。なるべく短くしようとしました。
    
*/
//================================================================================================================================================
//--Consts and Debug---------------------------------------------------------------------------------------------------------------
#define ZERO (min(iFrame,0))

// Mathematical Consts
#define PHI 1.61803399
#define PI  3.14159265
#define TAU 6.28318533

// Debug (Orbital Mouse Controls) camera settings
#define USE_DEBUG_CAMERA
#define DEBUG_CAMERA_DIST 1.0
#define DEBUG_CAMERA_TARGET vec3(0.0, 0.0, 0.0)

// Cinematic camera settings
#define CAMERA_DIST 1.6
#define CAMERA_TARGET vec3(0.0, -.1, 0.0)

// Debug Rendering settings
//#define RENDER_ROPE
//#define RENDER_PILLARS
//#define RENDER_BRIDGES

// These Slow flags are almost solely responsible for the long compile time. They are important for ensuring domain repetition SDF
//   correctness but I surmise they cause the compiler to unwind a bunch of complex SDF code, causing the slowdown.
//   TODO: Look into ways to fix (or alternative approaches).

// Fast animated bridges (& pillars) do not correct discontinuous domain rep SDF and will cause artifacts, but are considerably faster.
//  TODO: I'm fairly confident this can be worked around by just doing 3 SDF evaluations. I end up doing 3 anyway to ensure correctness???
#define SLOWER_BRIDGES 0
// Fast pillars are not reccommended (very broken).
#define SLOWER_PILLARS 1

//--Bridge Construction Params-----------------------------------------------------------------------------------------------------
const float _BelowCloudBottom = 5.0;
const float _CloudBottomExtra = 150.0; // For blending into the clouds.

const float _BridgeStrutInterval = 3.50;
const float _BridgeNRep = 3.0;
const float _BridgeStrutOffset = 0.8;
const float _BridgeStrutThickness = 0.2;
const float _BridgeStrutRoundness = 0.01;

const float _BridgeTopWidth = 0.35;
const float _BridgeTopThickness = 0.25;
const float _BridgeTopBevelExtrusion = 0.03;
const float _BridgeTopBevelThickness = 0.05;
const float _BridgeTopRoundness = 0.005;

const float _BridgeWedgeHeight = 0.3;
const float _BridgeWedgeTopWidth = 2.2 * _BridgeStrutThickness;
const float _BridgeWedgeBevelHeight = 0.05;
const float _BridgeWedgeBevelExtrusion = 0.05;
const float _BridgeWedgeBevel2Ratio = 0.618; // height of second bevel along wedge

const float _BridgeStrutBoxFrameExtrusion = 0.04;
const float _BridgeStrutBoxFrameThickness = 0.03;

const float _BridgeBoxFrameBevelExtrusion = 0.03;
const float _BridgeBoxFrameBevelHeight = 0.04;
const float _BridgeBoxFrameBevelSep = 0.4;
const float _BridgeBoxFrameBevelBottomOffset = 1.618;

const float _BridgeLinkTopSpacing = 0.4;
const float _BridgeLinkThickness = 0.025;

const float _BridgeMiniStrutInterval = _BridgeStrutInterval / 4.0;
const float _BridgeMiniStrutNRep = 13.0;
const float _BridgeMiniStrutOffset = _BridgeStrutOffset - 0.5 * _BridgeMiniStrutInterval;
const float _BridgeMiniStrutHeight = 0.35;
const float _BridgeMiniStrutThickness = 0.075;
const float _BridgeMiniStrut_Y_Extrusion = 0.05;
const float _BridgeMiniStrut_Z_Extrusion = 0.07;

//--Infinite Bridge Construction Params--------------------------------------------------------------------------------------------
const float _BelowCloudBottomPillar = 10.0;

const float _InfBridgeAnimSpeed = 0.3;
const float _InfBridgeLowOffset = 8.0;
const float _InfBridgeAnimSegLen = 6.0;

//--Pillar Construction Params-----------------------------------------------------------------------------------------------------
const float _PillarRoundness = 0.001;
const float _PillarBevelRoundness = 0.01;
const float _PillarBevelExtrusion = 0.1;
const float _PillarBevelHeight = 0.075;
const float _PillarBevelNRep = 4.0;

const float _PillarVBevelExtrusion = 0.05;
const float _PillarVBevelThickness = 0.05;
const float _PillarVBevelRoundness = 0.03;
const float _PillarVBevel_N = 16.0; // TODO: High n looks bad at far distances due to aliasing, so LOD this.

const float _PillarCapHeight  = 0.47;
const float _PillarCapExtrusion = 0.45;

const float _PillarSegHeight = 4.3;

const float _LittlePillar_H = 1.0;
const float _LittlePillar_R = 0.125;
const float _LittlePillar_N = 3.0; // smaller numbers (generally odd) give better silhouettes across multiple viewing angles

//--Rope Params--------------------------------------------------------------------------------------------------------------------
const vec3 _ShideWindParams = vec3(0.10, 4.0, 1.75); // Wave amplitude, distance modifier (stiffness), anim speed
const vec3 _ShideWindParams_s = vec3(0.013, 72.0, 3.5); // Smaller sub-waves
const vec2 _KiraretanawaWindParamsYZ = vec2(PI / 20.0, 2.3); // Max rot, anim speed
const vec2 _KiraretanawaWindParamsYX = vec2(PI / 30.0, 2.0);

//--Cloud Params-------------------------------------------------------------------------------------------------------------------
const vec3 _CloudSigmaS = vec3(1.0); // Scattering coeff
//const vec3 _CloudSigmaS = vec3(1.0, 0.7, 0.7); // Scattering coeff
const vec3 _CloudSigmaA = vec3(0.0); // Absorption coeff
const vec3 _CloudSigmaE = max(_CloudSigmaS + _CloudSigmaA, vec3(1e-6)); // Extinction coeff

const float _CloudWidth = 1500.0;
const float _CloudBot = -0.5*_CloudWidth;
const float _CloudTop = 10.0;
const vec3  _CloudMinCorner = vec3(-_CloudWidth, _CloudBot, -_CloudWidth);
const vec3  _CloudMaxCorner = vec3(_CloudWidth, _CloudTop, _CloudWidth);

const float _CloudDensityMultiplier = 0.075;

const float _CloudShapeSpeed  = -5.0;
const float _CloudShapeSize  = 0.05;
const float _CloudShapeStrength  = 0.7;

const float _CloudDetailSpeed = -10.0;
const float _CloudDetailSize = 0.3;
const float _CloudDetailStrength = 0.2;

//--Materials----------------------------------------------------------------------------------------------------------------------
#define MAT_ROPE 1.0
#define MAT_SHIDE 2.0
#define MAT_SHIDE_SECONDARY 3.0
#define MAT_BRIDGE_STONE 4.0
#define MAT_BRIDGE_BRASS 5.0
#define MAT_PILLAR_STONE 6.0
#define MAT_PILLAR_GOLD  7.0
#define MAT_DEBUG 10000.0

//const vec3 _MatRope = vec3(0.95, 0.89, 0.74);

//const vec3 _RopeTerminatorLineCol = 0.8*vec3(0.49, 0.329, 1.0);

//const vec3 _MatBridgeStone = 2.0*vec3(0.361, 0.329, 0.370);
const vec3 _MatBridgeStone = 0.5*vec3(0.361, 0.329, 0.370);

//const vec3 _MatBridgeBrass = vec3(0.940, 0.841, 0.517);

//const vec3 _MatBridgeBrass = vec3(0.840, 0.730, 0.370);
const vec3 _MatBridgeBrass = 0.5*vec3(0.840, 0.730, 0.370);

const vec3 _MatBridgeBrassSpe = 0.5*vec3(0.370, 0.840, 0.832);
//const vec3 _MatBridgeBrassSpe = vec3(0.370, 0.840, 0.832);


//--Illumination-------------------------------------------------------------------------------------------------------------------
const vec3  _SunPos  = vec3(30.0, 20.0, 30.0);
//const vec3 _SunPos = vec3(0.2, 56, -40.1);

const vec3  _LightDir = normalize(_SunPos);

//const vec3  _SunCol  = vec3(0.51471, 0.79919, 1.0);

const float _ZenithAttenuation = 1.8;
const float _NadirAttenuation  = 1.2;
const float _HorizonOffset = 0.0;

const float _FogFrontRadius = 100.0;
const float _RopeExtraShadowBrightness = 0.025; // Set to zero for flat shaded look

//--Colour Schemes-----------------------------------------------------------------------------------------------------------------
//----0: Sunset----------------------------------------------------------------------------------------------------
#if COLOUR_SCHEME == 0
// Sun
const float _SunSize = 3.5;
const vec3  _SunCol = vec3(0.98, 0.72, 0.31);
const float _SunBrightness = 1.4;
const float _SunHaloAttenuation = 1.4;
const float _SunHaloRadius = 3.0;

// Atmosphere
const vec3 _ZenithCol = 0.5 * vec3(0.37, 0.14, 0.43);
const vec3 _HorizonCol = 1.5 * vec3(0.36, 0.17, 0.66);
const vec3 _NadirCol = vec3(0.35, 0.32, 0.8);
const vec3 _FogColour = _HorizonCol;

// Materials
const vec3 _MatRope = 0.60*vec3(0.95, 0.89, 0.74);
const vec3 _RopeTerminatorLineCol = 0.8 * vec3(1.0, 0.329, 0.518);

const vec3 _MatShide = vec3(1.0, 1.0, 1.0);
const vec3 _MatShideSecondary = 0.8*vec3(1.0, 0.412, 0.412);

const vec3 _MatPillarStone = 0.8*vec3(0.969, 0.961, 0.918);
const vec3 _MatPillarStoneFre = vec3(0.471, 0.737, 0.941);
const vec3 _MatPillarStoneHardlight = vec3(1.0, 0.957, 0.882);
#endif

//----1: Night-----------------------------------------------------------------------------------------------------
#if COLOUR_SCHEME == 1
// Sun
const float _SunSize = 3.0;
const vec3  _SunCol = vec3(0.65, 0.8, 1.0);
const float _SunBrightness = 1.0;
const float _SunHaloAttenuation = 2.0;
const float _SunHaloRadius = 2.0;

// Atmosphere
const vec3 _ZenithCol = 0.02 * vec3(0.32, 0.65, 1.0);
const vec3 _HorizonCol = 0.04 * vec3(0.32, 0.65, 1.0);
const vec3 _NadirCol = 0.02 * vec3(0.32, 0.65, 1.0);
const vec3 _FogColour = _HorizonCol;

// Materials
const vec3 _MatRope = 0.80*vec3(0.89, 0.631, 0.341);
const vec3 _RopeTerminatorLineCol = 0.5 * vec3(0.859, 0.761, 0.596);

const vec3 _MatShide = vec3(1.0, 1.0, 1.0);
const vec3 _MatShideSecondary = vec3(1.0, 1.0, 1.0);

const vec3 _MatPillarStone = 0.15*vec3(0.969, 0.961, 0.918);
const vec3 _MatPillarStoneFre = vec3(0.471, 0.737, 0.941);
const vec3 _MatPillarStoneHardlight = vec3(1.0, 0.957, 0.882);
#endif
//-----------------------------------------------------------------------------------------------------------------

// SSS Sun Outline
const vec3  _OutlineCol = _SunCol;
const float _OutlineDotThreshold = 0.95;
const float _OutlineRadialAttenuation = 2.0;
const float _OutlineMaxDist = 0.05;
const float _OutlineAttenuation = 16.0;
const float _OutlineExtraThickness = 0.0035;

//==HELPER FUNCS==================================================================================================================================

// Material Selection Helper Funcs
#define CMP_MAT_LT(a, b) a < (b + 0.5)
#define CMP_MAT_GT(a, b) a > (b - 0.5)

// material comparisons
vec2 MIN_MAT(vec2 a, vec2 b) { return (a.x < b.x) ? a : b; }

// util functions for AO
// https://www.shadertoy.com/view/ld3Gz2
float hash1(float n) { return fract(sin(n)*43758.5453123); }
float hash1( vec2 n) { return fract(43758.5453123*sin(dot(n,vec2(1.0,113.0)))); }
float hash2( vec2 n) { return fract(74193.9957812*sin(dot(n,vec2(99.0,33.0)))); }

vec3 forwardSF(float i, float n) 
{
    float phi = 2.0*PI*fract(i/PHI);
    float zi = 1.0 - (2.0*i+1.0)/n;
    float sinTheta = sqrt( 1.0 - zi*zi);
    return vec3( cos(phi)*sinTheta, sin(phi)*sinTheta, zi);
}

// Phase function for cloud isotropic scattering
//  See PBR: From Theory to Implementation, Section 11.3
float henyeyGreenstein(float g, float costh) 
{ 
    return (1.0 / (4.0 * PI)) * ((1.0 - g * g) / pow(1.0 + g*g - 2.0*g*costh, 1.5)); 
}

//==SDF===========================================================================================================================================

//--Shimenawa----------------------------------------------------------------------------------------------------------------------
float sdShide( in vec3 p, in int s_n, in float sec_id, out float seg_id )
{
    // 紙垂 (Shide) Multiple zig-zag boxes    
    
    // Rotate the whole thing 45c
    p.yz = (p.yz + vec2(p.z, -p.y))*sqrt(0.5); // Shortcut for 45-degrees rotation

    // Dimensions
    vec3 dim = 0.02 * vec3(0.03, 1.0, 2.0); // 2:1:0.1 ratio
    vec3 smallDim = dim;
    smallDim.y /= 2.0;
    smallDim.z *= 0.9;
           
    // 風から動きの準備
    // 垂直動きが欲しいので、回りの後でする。
    // 接続点から距離を使って動きの強さが決まる
    vec3 connectorPos = vec3(0., 2.0 * dim.y, -1.5 * dim.z);
    float distToConnector = length(connectorPos - p);

    // Thinner box, mount to rope    
    vec3 thinShidePos = p + vec3(.0, -1.5 * dim.y, 0.9 * dim.z);
    // Sinusoidal motion for wind
    thinShidePos.x += _ShideWindParams.x * sin(distToConnector * _ShideWindParams.y + _ShideWindParams.z * iTime + 53.0 * sec_id) * distToConnector;
    float d = sdBox(thinShidePos, smallDim);
    
    // Rest of the zig-zag
    seg_id = 0.0;
    for (int i=0; i<s_n; i++)
    {
        vec3 shidePos = p + vec3(-5.0 * dim.x * mod(float(i+1), 2.0), 2.0 * dim.y * float(i), -dim.z * float(i));
        shidePos.x += _ShideWindParams.x * sin(distToConnector * _ShideWindParams.y + _ShideWindParams.z * iTime + 53.0 * sec_id) * distToConnector;
        shidePos.x += _ShideWindParams_s.x * sin(distToConnector * _ShideWindParams_s.y + _ShideWindParams_s.z * iTime + 13.0 * sec_id) * distToConnector;
        float d2 = sdBox(shidePos, dim);
        if (d2 < d) {d = d2; seg_id = float(i+1);}
    }
    
    return d;
}

// p is position of top, s is scale
float sdKiraretanawa( in vec3 p, in float s, in float connectorOffset, in float id )
{
    // TODO: Scaling is done in the SDF calls (i.e. by scaling dimensions by s), 
    //   but could be done by a simple pre-scaling of the domain (i.e. some q = p / s).
    //   As the same scale is applied to all in the domain repetition, this prescaling wouldn't affect the SDF continuity.
    
    // Offset rotation to make object sway
    p.yz *= rot(_KiraretanawaWindParamsYZ.x * sin(_KiraretanawaWindParamsYZ.y * 0.89 * iTime + 71.0 * id));
    p.yx *= rot(_KiraretanawaWindParamsYX.x * cos(_KiraretanawaWindParamsYX.y * 0.57 * iTime + 31.0 * id));
    p.y += connectorOffset;
        
    float d = sdCone(p + s * vec3(0., 0.6, 0.), vec3(0.), s * vec3(0.,0.5,0.), s * 0.32, s * 0.1) - s * 0.05;
    d = min(d, sdCappedCylinder(p, s * 0.06, s * 0.12) - s * 0.015);
    
    {
    // Binding rope, 3 for cost of 1
    vec3 q = p;
    q.y = abs(p.y + s * 0.02) - s * 0.0165; // XZ plane symmetry for duplication, and translate downwards
    q.y = abs(q.y) - s * 0.0165;
    d = min(d, sdTorus(q, vec2(s * 0.15, s * 0.013)));
    }
    
    // Connecting rope
    d = min(d, sdCappedCylinder(p - vec3(0.0, connectorOffset/2.0, 0.0), connectorOffset, s * 0.013));
    
    return d;
}


float sdSwirl(vec3 p, in float r, in float c, in float f, out float id)
{
    // Rotate and twist a capsule to make the rope along the x axis
    float l = TAU * r; // circumference rope length
    
    p.yz*=rot(p.x*PI*f); // twisting; 45-deg per x
    vec2 sector = step(0.0,p.yz);
    id = sector.x + 2.0 * sector.y;
    p.y = abs(p.y) - 0.02;
    
    float d = sdVerticalCapsule(p,c,l); //Potentially better solution exists with capped torus
    return d;
}

vec2 sdShimenawa(in vec3 p, in float r, in float c, in float f, out float id)
{
    vec2 res = vec2(1e10);
    
    // Bounding sphere
    float d = sdSphere(p, r + 0.2);
    if (d > 0.5) return vec2(d, 1e10);

    float ring_t = 0.0;
    float dRing = sdCircleXZ(p, r, ring_t);
    ring_t = 1.0 - ring_t; // reflect rope direction
    
    //p.y = p.y + 0.06*sin(2.0*2.0*PI*ring_t+1.5*time);
    
    // Swirly Rope
    if (dRing < 0.8) // approx bounding ring
    {    
        // I'm fairly confident this local-UV transformation majorly messes up the SDF for isolines > 0.0 but... oh well.
        //vec3 ring_uv = vec3(ring_t * TAU * r, p.y + 0.06*sin(2.0*2.0*PI*ring_t+1.5*time), dRing);
        vec3 ring_uv = vec3(ring_t * TAU * r, p.y, dRing);
        res = vec2(sdSwirl(ring_uv, r, c, f, id), MAT_ROPE);
    }
    
    // Shide
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
    
    // Cut Ropes
    {
    // TODO: This rotation is horrendously jank - I *should* fix it.
    vec3 q_s = p;
    
    // Prerotate offset
    float an = TAU/14.0;
    q_s.xz *= rot(an);
    
    an = TAU/7.0;
    float sector = round(atan(q_s.z,q_s.x)/an);
    float angrot = sector*an;
    q_s.xz *= rot(angrot);
    float d = sdKiraretanawa(q_s - vec3(r, -1.4*c, 0.0), 0.14, 0.05, sector+1.0);
    res = MIN_MAT(res, vec2(d, MAT_ROPE));
    }
    
    res.x *= 0.8; // Deal with messed up SDF T^T
    return res;
}

//--Pillars------------------------------------------------------------------------------------------------------------------------
vec2 sdPillarSeg( 
    in vec3 p, in float r, in float seg_h,
    in float seg_n_sep, in float seg_n_bevels, in float n_small_pillars )
{
    // Radial scale - 1.0 is base scale, can be changed.
    //   Want to keep nice ratio between main pillar radius and the size of objects which comprise the gap,
    //   so scale these objects by a factor in comparison to a base radius which has a aesthetically pleasing ratio.
    float rs = r / 1.0; 
    
    // I am so sorry if you happen to unfortunately read this
    float rs_PillarCapHeight = rs*_PillarCapHeight;
    float rs_PillarCapExtrusion = rs*_PillarCapExtrusion;
    float rs_LittlePillar_H = rs*_LittlePillar_H;
    float rs_LittlePillar_R = rs*_LittlePillar_R;
    float rs_PillarBevelExtrusion = rs*_PillarBevelExtrusion;
    float rs_PillarVBevelExtrusion = rs*_PillarVBevelExtrusion;
    float rs_PillarVBevelThickness = rs*_PillarVBevelThickness;
    
    // TODO: Could LOD this pillar segment.
    vec2 res = vec2(1e10);
    
    float pillar_height = seg_h - 0.5*rs_PillarCapHeight;
    
    // Bounding Volume (Exact Cylinder)
    float d = sdCappedCylinder(p, 2.0*(seg_h + rs_LittlePillar_H) + rs_PillarCapHeight + _PillarBevelHeight, r + rs_PillarCapExtrusion);
    if (d > 1.0) return vec2(d, 1e10);
    
    // big cylinder
    vec3 cyl_base = p + vec3(0.0, 0.5*rs_PillarCapHeight, 0.0);
    res = MIN_MAT(res, vec2(
        sdCappedCylinder(cyl_base, pillar_height, r) - _PillarRoundness,
        MAT_PILLAR_STONE
    ));
    
    // horizontal bevels
    {
    vec3 q = cyl_base + vec3(0.0, pillar_height, 0.0);
    float interval = 2.0*pillar_height / seg_n_sep;
    float rep_id = clamp(round((q.y / interval)), 0.0, seg_n_sep);
    q.y -= rep_id * interval;
    res = MIN_MAT(res, vec2(
        sdCappedCylinder(q, _PillarBevelHeight, r + rs_PillarBevelExtrusion) - _PillarBevelRoundness,
        MAT_PILLAR_GOLD
    ));
    }
    
    // vertical bevels (to look like pillar indents)
    {
    vec3 q = cyl_base;
    
    //   rotatational domain repetition
    float an = TAU/seg_n_bevels;
    float sector = round(atan(q.z,q.x)/an);
    float angrot = an*(sector);
    q.xz *= rot(angrot);
    
    res = MIN_MAT(res, vec2(
        sdBox(q - vec3(r, 0.0, 0.0), vec3(rs_PillarVBevelExtrusion, pillar_height, rs_PillarVBevelThickness)) - _PillarVBevelRoundness,
        MAT_PILLAR_STONE
    ));
    }
    
    // small radial sub-pillars
    {
    vec3 q = p;
    q.y -= seg_h + rs_LittlePillar_H;
    
    //   rotatational domain repetition
    float an = TAU/n_small_pillars;
    float sector = round(atan(q.z,q.x)/an);
    float angrot = an*(sector);
    q.xz *= rot(angrot);

    res = MIN_MAT(res, vec2(
        sdCappedCylinder(q - vec3(0.7*r, 0.0, 0.0), rs_LittlePillar_H, rs_LittlePillar_R),
        MAT_PILLAR_STONE
    ));
    }
    
    // Cone caps
    {
    vec3 q = p;
    q.y -= seg_h - rs_PillarCapHeight;
    float interval = 2.0*rs_LittlePillar_H + rs_PillarCapHeight;
    float rep_id = clamp(round((q.y / interval)), 0.0, 1.0);
    q.y -= rep_id * interval;
    res = MIN_MAT(res, vec2(
        sdCone(q, vec3(0.0), vec3(0.0, rs_PillarCapHeight, 0.0), r, r + rs_PillarCapExtrusion) - _PillarRoundness,
        MAT_PILLAR_STONE
    ));
    }
        
    return res;
}

vec2 sdPillar( 
    in vec3 p, in float n_rep, in float r, in float seg_h, 
    in float seg_n_sep, in float seg_n_bevels, in float n_small_pillars )
{    
    float rs = r / 1.0; // Factor in radial scale of pillar ( see comment in sdPillarSeg() )
    float full_pillar_height = 2.0*(seg_h + rs*_LittlePillar_H) + rs*_PillarCapHeight;
    p.y += _BelowCloudBottomPillar - 0.5 * full_pillar_height;
    
    // check i = n'th repetition for current rep's SDF
    float rep_id = clamp(round((p.y / full_pillar_height)), 0.0, n_rep);
    p.y -= rep_id * full_pillar_height;

    vec2 res = sdPillarSeg(p, r, seg_h, seg_n_sep, seg_n_bevels, n_small_pillars);
    
    if (rep_id < 1.0) return res;
    
    // check i = n-1'th repetition to ensure SDF correctness.
    p.y += full_pillar_height;
    res = MIN_MAT(res, sdPillarSeg(p, r, seg_h, seg_n_sep, seg_n_bevels, n_small_pillars));

    return res;
}

vec2 sdPillars( in vec3 p )
{
    // TODO: Make pillar field interesting by modifying certain parameters 
    //    - Bottom offset between 0 -> 1 * h (so pillars have 'phase')
    //    - N small pillars, probably 2 -> 6
    //    - Big pillar radius
    //    - N repetitions (further pillars can be taller)
    //    - N H and V bevels (further pillars should have fewer)
    //    
    
    // TODO: Pillar variations - would like a whiter, but thinner variation through materials.
        
    // TODO: BB on wedge between cloud bottom and max pillar height.
    
    const vec2 spacing = vec2(145.0);
    const float len_spacing = length(spacing);

    // Maximum and minimum radial extents of the domain repetition defined by 2 params (r, t)
    //  TODO: I would like this to be independent of spacing, but it messes with the SDF so i cba
    const float ring_rad = 2.5;
    const float ring_thickness = 6.0;
    const float ext_min = ring_rad - ring_thickness;
    const float ext_max = ring_rad + ring_thickness;
    
    vec2 res = vec2(1e10);
    
    //return sdPillar(p, 0.0, 3.0, 4.3, 4.0, 16.0, 3.0 );
        
    vec2 base_id = round(p.xz / spacing);
    float len_base_id = length(base_id);
    
    // Early exit on maximum extent for opt (no more pillars beyond this point so raymarch accuracy not affected)
    if (len_base_id > ext_max) return vec2(1e10);
    
    #if SLOWER_PILLARS
    for (int i=-1; i<2; i++)
    for (int j=-1; j<2; j++)
    {
    #else
    {
    int i = 0; int j = 0;
    #endif 
    
    vec2 id = base_id + vec2(i, j);
    float len_id = length(id);
    float id_hash1 = hash1(id);
    float id_hash2 = hash2(id);
    
    // Near and far heights - for a dynamic scene, make further away pillars *considerably* larger.
    float far = step(ring_rad, len_id);
    
    float height_dev = (far == 1.0) ? 5.0 : 2.0;
    float height_avg = (far == 1.0) ? 8.0: 5.0;
    float height_bias = 2.0 * (smoothstep(ext_min, ext_max, len_id) * id_hash1) - 1.0; // generate taller pillars when further away
    float height = height_avg + height_dev * height_bias;
    height *= 1.2;

    float nseg = 2.0 + floor(2.0*id_hash2);
    
    float quick_hash3 = fract(id_hash1 * id_hash2);
    float n_pillar = 2.0 + 2.0*floor(2.0 * quick_hash3) + 1.0; // 3 or 5 pillars for better silhouette (even is bad)
    float h_offset = len_id * 2.0 * quick_hash3;
    
    vec3 q = p;
    q.xz = p.xz - spacing * id;
    q.xz *= rot(id_hash1 * 2.0 * PI);
    q.xz += 0.53 * spacing * (sin(131.5 * id.x + 13.8 * id.y + vec2(1.5, 0.0)));
    q.y += h_offset;

    //res = MIN_MAT( res, sdPillar(q, 2.0, height, 4.3, 4.0, 16.0, 3.0 ) );
    res = MIN_MAT( res, sdPillar(q, nseg, 0.8*height, 3.0*height, 4.0, 16.0, n_pillar ) );
    //res = MIN_MAT(res, vec2(sdCappedCylinder(q+vec3(0.0, _BelowCloudBottomPillar - height, 0.0), height, height / 4.3), MAT_PILLAR_STONE));
    
    }
    
    res.x = max(res.x, abs(len_base_id - ring_rad) - ring_thickness);
    
    return res;

    // Smaller pillars closer to origin, some really large ones in the distance.
    //return sdPillar(p + vec3(0.0, 0.0, 0.0), 3.0, 1.0, 1.0);
}

//--Bridges------------------------------------------------------------------------------------------------------------------------
vec2 sdBridgeStrut( in vec3 strut_base, in float h) 
{
    vec2 res = vec2(1e10);
    
    // strut beam
    vec3 q = strut_base;
    res = MIN_MAT(res, vec2(
        sdBox(q, vec3(_BridgeStrutThickness, h, _BridgeStrutThickness)) - _BridgeStrutRoundness, 
        MAT_BRIDGE_STONE
    ));

    // strut-bridge connector wedge
    q = strut_base - vec3(0.0, h - _BridgeWedgeHeight, 0.0);
    res = MIN_MAT(res, vec2(
        sdPrism(q, _BridgeStrutThickness, _BridgeWedgeTopWidth, _BridgeWedgeHeight, _BridgeStrutThickness) - _BridgeStrutRoundness, 
        MAT_BRIDGE_STONE
    ));

    // wedge bevels
    //   first wedge
    q = strut_base - vec3(0.0, h, 0.0);
    res = MIN_MAT(res, vec2(
        sdBox(q, vec3(
            _BridgeWedgeTopWidth + _BridgeWedgeBevelExtrusion, 
            _BridgeWedgeBevelHeight, 
            _BridgeStrutThickness + _BridgeWedgeBevelExtrusion)
        ) - _BridgeStrutRoundness,
        MAT_BRIDGE_BRASS
    ));
    //   second wedge
    q.y += 2.0 * _BridgeWedgeHeight * (1.0 - _BridgeWedgeBevel2Ratio);
    res = MIN_MAT(res, vec2( 
        sdBox(q, vec3(
            mix(_BridgeStrutThickness, _BridgeWedgeTopWidth, _BridgeWedgeBevel2Ratio) + _BridgeWedgeBevelExtrusion, 
            _BridgeWedgeBevelHeight / 2.0, 
            _BridgeStrutThickness + _BridgeWedgeBevelExtrusion)
        ) - _BridgeStrutRoundness,
        MAT_BRIDGE_BRASS
    ));

    // box frame
    q = strut_base; q.y += _BridgeWedgeHeight;
    float box_frame_height = h - _BridgeWedgeHeight;
    res = MIN_MAT(res, vec2(
        sdBoxFrame(q, vec3(
            _BridgeStrutThickness + _BridgeStrutBoxFrameExtrusion, 
            box_frame_height, 
            _BridgeStrutThickness + _BridgeStrutBoxFrameExtrusion
        ), _BridgeStrutBoxFrameThickness) - _BridgeStrutRoundness,
        MAT_BRIDGE_STONE
    ));

    // box frame bevels
    float box_frame_top_offset = _BridgeWedgeHeight - box_frame_height + _BridgeBoxFrameBevelHeight;
    // top
    q = strut_base;
    q.y += box_frame_top_offset; // translate up along strut
    q.y = abs(q.y + 0.5 * _BridgeBoxFrameBevelSep) - 0.5 * _BridgeBoxFrameBevelSep; // domain repetition (reflection in xz plane), and center about top bevel.
    float bevel_width = _BridgeStrutThickness + _BridgeStrutBoxFrameExtrusion + _BridgeBoxFrameBevelExtrusion;
    res = MIN_MAT(res, vec2(
        sdBox(q, vec3(
            bevel_width,
            _BridgeBoxFrameBevelHeight,
            bevel_width)
        ) - _BridgeStrutRoundness,
        MAT_BRIDGE_BRASS
    ));

    // bottom
    q = strut_base;
    q.y += box_frame_top_offset + _BridgeBoxFrameBevelSep + _BridgeBoxFrameBevelBottomOffset; // translate up along strut
    res = MIN_MAT(res, vec2(
        sdBox(q, vec3(
            bevel_width,
            _BridgeBoxFrameBevelHeight,
            bevel_width)
        ) - _BridgeStrutRoundness,
        MAT_BRIDGE_BRASS
    ));


    // link
    q.y -= 0.5 * _BridgeBoxFrameBevelBottomOffset;
    q.z -= _BridgeStrutThickness + _BridgeLinkThickness;
    res = MIN_MAT(res, vec2(
        sdLink(q, 0.65*0.5*_BridgeBoxFrameBevelBottomOffset, 0.10, _BridgeLinkThickness),
        MAT_BRIDGE_BRASS
    ));
    
    return res;
}

vec2 sdBridgeTop( in vec3 p, in float h, in float l ) 
{
    vec2 res = vec2(1e10);
    
    vec3 q = p - vec3(l, _BridgeTopThickness, 0.0);
    res = MIN_MAT(res, vec2(
        sdBox(q, vec3(l, _BridgeTopThickness, _BridgeTopWidth)),
        MAT_BRIDGE_STONE
    ));
    // bridge top bevels
    q.y = abs(q.y) - _BridgeTopThickness + _BridgeTopBevelThickness; // domain repetition
    res = MIN_MAT(res, vec2(
        sdBox(q, vec3(l, _BridgeTopBevelThickness, _BridgeTopWidth + _BridgeTopBevelExtrusion)),
        MAT_BRIDGE_BRASS
    ));
    
    return res;
}

vec2 sdBridgeTopStruts( in vec3 p, in float h, in float l )
{
    vec2 res = vec2(1e10);
    
    // mini struts        
    //   block
    vec3 q = p;
    //      domain repetition (reflection in xz plane), and center about top bevel.
    q.z = abs(q.z) - _BridgeTopWidth + _BridgeMiniStrutThickness - _BridgeMiniStrut_Z_Extrusion; 
    res = MIN_MAT(res, vec2(
        sdBox(q, vec3(_BridgeMiniStrutThickness, _BridgeMiniStrutHeight, _BridgeMiniStrutThickness)),
        MAT_BRIDGE_STONE
    ));

    //   spike
    q.y += _BridgeMiniStrutHeight + _BridgeMiniStrutThickness;
    q.y = -q.y;
    res = MIN_MAT(res, vec2(
        sdPyramid(q, _BridgeMiniStrutThickness, _BridgeMiniStrutThickness, _BridgeMiniStrutThickness),
        MAT_BRIDGE_BRASS
    ));
    
    return res;
}

vec2 sdBridgeSegmentLOD( in vec3 p, in float h, in float l )
{
    // TODO: If we LOD the bridges, we should be able to a few more.
    return vec2(1e10);
}

vec2 sdBridgeSegment( in vec3 p, in float h, in float l, in float s )
{
    // TODO: These arches most likely should only cast shadows on eachother, not the scene foreground.
    // TODO: I made this quite detailed... might have to cut back if optimisation is not enough with BB, etc.    
    vec2 res = vec2(1e10);
            
    // Transform to correct origin plane in y
    //p.y += _BelowCloudBottom - h;
    
    // Bounding box
    // TODO: If there were a way to move this to before the raymarch, to avoid doing any raymarching at all,
    //         this would be a lot faster.
    // TODO: Could create a more exact bounding box here if needed.
    float d = sdBox(p, vec3(2.0*l, h + _BridgeTopThickness + 0.3, _BridgeTopWidth + 0.05));
    if (d > 1.0) return vec2(d, 1e10);
    
    // Domain repetition for struts - render as many struts as we like for the price of one.
    vec3 strut_base = p - vec3(_BridgeStrutOffset, 0.0, 0.0);
    float rep_id = clamp(round((strut_base.x / _BridgeStrutInterval)), 0.0, _BridgeNRep);
    strut_base.x -= rep_id * _BridgeStrutInterval;
    
    res = MIN_MAT(res, sdBridgeStrut(strut_base, h));
    
    vec3 top_base = p - vec3(0.0, h, 0.0);
    res = MIN_MAT(res, sdBridgeTop(top_base, h, l));

    // Domain repetition for struts - render as many struts as we like for the price of one.
    vec3 mini_strut_base = p - vec3(_BridgeMiniStrutOffset, h + _BridgeMiniStrutHeight - _BridgeMiniStrut_Y_Extrusion, 0.0);
    rep_id = clamp(round((mini_strut_base.x / _BridgeMiniStrutInterval)), 0.0, _BridgeMiniStrutNRep);
    mini_strut_base.x -= rep_id * _BridgeMiniStrutInterval;

    res = MIN_MAT(res, sdBridgeTopStruts(mini_strut_base, h, l));
    
    return res;
}

vec2 sdInfiniteBridge( in vec3 p, in vec3 o, in float an, in float h, in float scale )
{
    // TODO: BB on wedge between cloud bottom and bridge height (top).

    // Map bridge to line in XZ plane with line origin (o) and direction in XZ defined by angle from x basis vector (an).
    // Translate to point on bridge line
    p -= o;
    
    // Rotate domain to fit bridge to desired line
    p.xz *= rot(an);

    float seg_l = 12.5; 
    p.x = p.x - 2.0*seg_l * round((p.x - seg_l) / (2.0*seg_l));;
    p.y += _BelowCloudBottom - h;

    p /= scale;
    h /= scale;
    vec2 res = sdBridgeSegment(p, h, seg_l, scale);
    res.x *= scale;
    
    return res;
}

vec2 sdInfiniteBridgeAnimated( 
    in vec3 p, in vec3 o, in float an, in float h,
    const float sep, const float phase, const float n_seg, const float min_bound, const float max_bound
) {
    // TODO: I'm 99% sure it will be faster (for compilation, mostly) to just do 3 SDF calls: domain rep segment, falling, and rising sections.
    vec2 res = vec2(1e10);
    
    // Potentially 3x the work.. but the SDF is correct.
    //  The inclusion of this for loop also makes compilation very slow. Not exactly sure why but I don't think its due to 
    //     Suggestions for fixes are welcome :) (if even possible - it is a lot of work being done here.)
    #if SLOWER_BRIDGES
    for (int i=-1; i<2; i++)
    {
    #else
    int i = 0;
    #endif
    
    vec3 q = p;
    q -= o;
    q.xz *= rot(an);
    
    float seg_l = 6.0;
    float rep_id = float(i) + (round((q.x - seg_l) / (2.0*seg_l)));

    q.x = q.x - 2.0*seg_l * rep_id; // domain repetition

    // which repetition id should be at its animation zenith
    float high_center_id = mod(_InfBridgeAnimSpeed*iTime, sep) + phase; 
    // make repetitions periodic
    float mod_rep_id = mod(rep_id, sep);
    
    // Get circular difference between two periodic ids (i.e. 19 - 0  = 1 (mod 20), not 19).
    //  This determines how far we are from the bridge arc zenith.
    //  as it turns out, this is the same as the distance between two elements of a ring 0 -> n-1: min(|i - j|, n - |i - j|).
    float diff = min(abs(mod_rep_id - high_center_id), sep - abs(mod_rep_id - high_center_id));
    
    // This optimisation is terrible and bugged as hell, but somehow it kind of works.
    if (diff > n_seg + 6.0) return vec2(1e10);
    
    // Clamp signal to: linear 0->1 if on edge, 1 if low, 0 if high
    float edge_diff = clamp(diff, n_seg, n_seg + 1.0) - n_seg;
    float diff_plus_one = clamp(diff, n_seg, n_seg + 2.0) - (n_seg + 1.0);
    edge_diff = pow(edge_diff, 3.5); // ease in edges with power curve
    edge_diff = (rep_id < min_bound || rep_id > max_bound) ? 1.0 : edge_diff; // Limit periodic function to certain range.
    
    q.y += _InfBridgeLowOffset * edge_diff;
    
    #if SLOWER_BRIDGES
    res = MIN_MAT(res, sdBridgeSegment(q, h, seg_l));
    }
    #else
    res = sdBridgeSegment(q, h, seg_l, 1.0);
    #endif
    
    // Clip low bridges
    float clip_dist = p.y + _BelowCloudBottom;
    res.x = max(-clip_dist, res.x); // Boolean subtraction

    return res;
}

vec2 sdCurvedBridge( in vec3 p, in float h, in float l, in float r )
{
    // Create bridge on z axis line, then wrap it around an arc with circle uvs
    float ring_t = 0.0;
    float dRing = sdCircleXZ(p, r, ring_t);
        
    // I'm fairly confident this local-UV transformation majorly messes up the SDF for isolines > 0.0 but... oh well.
    vec3 ring_uv = vec3(ring_t * TAU * r, p.y, dRing);
    vec2 res = sdBridgeSegment(ring_uv / 4.0, h, l, 1.0);
    res.x *= 0.8;
    
    return res;
}

//==ILLUMINATION================================================================================================================================
vec3 sky( in vec3 ro, in vec3 rd ) 
{
    // Make sun always appear as if viewed from a certain point to deal with the fact it
    //  is infact very close geometrically.
    vec3 sun_pov_origin = vec3(0.0);
    float sunDist = length(_SunPos);
    float dist = length(_SunPos - (sun_pov_origin + rd * sunDist)) - _SunSize;
    dist = dist < 0.0 ? 0.0 : dist;
    
    float ry = _HorizonOffset + rd.y;
        
    // Atmosphere
    float zenith  = 1.0 - pow(min(1.0, 1.0 - ry), _ZenithAttenuation);
    float nadir   = 1.0 - pow(min(1.0, 1.0 + ry), _NadirAttenuation);
    float horizon = 1.0 - zenith - nadir;
    
    vec3 skycol = zenith * _ZenithCol + nadir * _NadirCol + horizon * _HorizonCol;
    
    // Skybox Clouds?
    //  Also add a few wispy streaks for extra detail?
    
    float halo = pow((_SunHaloRadius/dist), _SunHaloAttenuation);
    vec3 sun = halo * _SunCol;
    skycol = 1.0 - exp(-(skycol + sun));
    skycol *= _SunBrightness;
    
    return skycol;
}

vec3 fog( in vec3 col, in float t, in vec3 rd )
{
    float amt = 1.0 - exp(-t/2500.0);
    float sun_factor = max(dot(rd, _LightDir), 0.0);
    
    //return col;
    return mix(col, _FogColour, amt);
}

//==RENDERING===================================================================================================================================

//--GEOMETRY-----------------------------------------------------------------------------------------------------------------------
vec2 mapBg( in vec3 p )
{
    // TODO: Arches are unlikely to interact with the rest of the scene significantly, so give them their own map for optimisation
    return vec2(-1.0);
}

vec2 map( in vec3 p )
{    
    vec2 res = vec2(1e10); // (Distance, Material)

    // Rope
    #ifdef RENDER_ROPE
    float r_id;
    res = MIN_MAT(res, sdShimenawa(p, 0.4615, 0.04, 10.0, r_id));
    #endif
    
    // Pillars
    #ifdef RENDER_PILLARS
    res = MIN_MAT(res, sdPillars(p + vec3(0.0, _CloudBottomExtra, 0.0)));
    #endif
    
    // Bridge Segment
    #if 0
    res = MIN_MAT(res, sdBridgeSegment(p, 2.5, 6.0));
    #endif

    // Static Inf Bridges
    #ifdef RENDER_BRIDGES
    {
    vec3 q = p + vec3(0.0, _CloudBottomExtra, 0.0);
    res = MIN_MAT(res, sdInfiniteBridge(q, 
        4.0*vec3(-25.0, 0.0, -25.0), 
        -PI / 4.0, 
        _CloudBottomExtra/2.0 + 2.5, 
        7.5));
        
    res = MIN_MAT(res, sdInfiniteBridge(q,
        4.0*vec3(-55.0, 0.0, -50.0),
        -7.0*PI / 12.0,
        _CloudBottomExtra/2.0 + 18.5,
        7.5));
    }
    #endif
    
    // Animated Inf Bridges
    #if 0
    {
    res = MIN_MAT(res, sdInfiniteBridgeAnimated(
        p, vec3(60.0, 0.0, 60.0), PI / 2.0, 2.5,
        40.0, 0.0, 1.0, -20.0, 20.0
    ));
    
    // TODO: There is a bug where the high bridges do not update correctly with the following.
    /*
    vec3 q = p;
    q.xz *= rot(PI);
    res = MIN_MAT(res, sdInfiniteBridgeAnimated(
        q, vec3(-35.0, 0.0, -35.0), PI / 2.0, 1.5,
        40.0, 2.0*TAU, 1.0, -20.0, 20.0, 1.0
    ));
    */
    }
    #endif
    
    // Curved bridge segment
    #if 0
    res = MIN_MAT(res, sdCurvedBridge(p, 2.5, 6.0, 15.0));
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

//--LIGHTING-----------------------------------------------------------------------------------------------------------------------
float calcSSS( in vec3 pos, in vec3 nor )
{
    const int N_SAMPLES = 9;
	float occ = 0.0;
    for( int i=ZERO; i<N_SAMPLES; i++ )
    {
        float h = 0.002 + 0.11*float(i)/7.0;
        vec3 dir = normalize( sin( float(i)*13.0 + vec3(0.0,2.1,4.2) ) );
        dir *= sign(dot(dir,nor));
        occ += (h-map(pos-h*dir).x);
    }
    occ = clamp( 1.0 - 11.0*occ/float(N_SAMPLES), 0.0, 1.0 );    
    return occ*occ;
}

// https://iquilezles.org/articles/rmshadows
// Soft Shadows with backtracking.
float softShadow( in vec3 ro, in vec3 rd, float k )
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

// TODO: I don't think this object scales well with distance / object size.
vec3 sunSSSOutline( in vec3 ro, in vec3 rd, in float d )
{
    // Create an outline around objects within a certain radius around the sun.
    //    The view direction dot product can be used to determine the radius.
    //    Raymarching near misses are used to generate the outline. 
    //   Then, a power curve is used to decay the outline based on distance from the edge of the sun.
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
    
    // TODO: One bounding volume might not be enough for this whole scene.
    //        Ideally could find a way to stop raymarching for specific objects without having to pass ro, rd into map()...
    vec2 tminmax = iSphere( ro, rd, 1000000000.0 );
	if( tminmax.y>0.0 )
    {
        // raymarch
        float t = max(tminmax.x,0.001);
        for( int i=0; i<RAYMARCH_MAX_STEPS && t<tminmax.y; i++ )
        {
            vec2 h = map(ro+t*rd);
            res.y = max(min(res.y, h.x), 0.0); // Track near misses for strong light outine
            if( abs(h.x)<0.001 ) { res=vec3(t, res.y, h.y); break; }
            t += h.x; // Coeff here is to try and avoid overshooting at the cost of performance
                            //  TODO: There should be a much smarter solution than this. The problem only arises with large domain distortions.
        }
    }
    
    return res; // t, nearest, mat_id
}

//--MATERIALS----------------------------------------------------------------------------------------------------------------------
vec3 shade( in vec3 ro, in vec3 rd, in float t, in float m ) 
{
    // TODO: Bridges and pillars should cast (relatively) sharp shadows on one another, so they need a separate map func.
    vec3 pos = ro + t*rd;
    vec3 nor = calcNormal(pos);
    //float shadow = softShadow(pos - 0.01*rd, _LightDir, 2.0); // Soft
    float shadow = softShadow(pos - 0.01*rd, _LightDir, 10.0); // Sharp
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
        //vec3 albedo = mix(base_shadow, _MatRope, shadow);
               
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
    else if (CMP_MAT_LT(m, MAT_BRIDGE_STONE))
    {
        float fre = clamp(1.0 + dot(nor, rd), 0.0, 1.0 );
        vec3 base_shadow = mix(0.6*_MatBridgeStone, _HorizonCol, 0.15);
        vec3 albedo = _MatBridgeStone + fre *_SunCol * _SunBrightness;
        return mix(base_shadow, albedo, shadow);
    }
    else if (CMP_MAT_LT(m, MAT_BRIDGE_BRASS))
    {
        float fre = clamp(1.0 + dot(nor, rd), 0.0, 1.0 );
        float ref = dot(reflect(rd, nor), normalize(_SunPos - pos));
        ref = smoothstep(0.7, 0.8, ref);
        vec3 base_shadow = mix(0.6*_MatBridgeBrass, _HorizonCol, 0.15);
        vec3 albedo = _MatBridgeBrass + (fre + ref) * (_SunCol * _SunBrightness);
        //vec3 albedo = _MatBridgeBrass + fre * (2.0 * _SunCol * _SunBrightness) + ref * _MatBridgeBrassSpe;
        
        return mix(base_shadow, albedo, shadow);
    }
    else if (CMP_MAT_LT(m, MAT_PILLAR_STONE))
    {
        float fre = clamp(1.0 + dot(nor, rd), 0.0, 1.0 );
        vec3 base_shadow = mix(0.6*_MatPillarStone, _HorizonCol, 0.15);
        vec3 albedo = _MatPillarStone + fre *_SunCol * _SunBrightness;
        return mix(base_shadow, albedo, shadow);
    }
    else if (CMP_MAT_LT(m, MAT_PILLAR_GOLD))
    {
        float fre = clamp(1.0 + dot(nor, rd), 0.0, 1.0 );
        float ref = dot(reflect(rd, nor), normalize(_SunPos - pos));
        ref = smoothstep(0.7, 0.8, ref);
        vec3 base_shadow = mix(0.6*_MatBridgeBrass, _HorizonCol, 0.15);
        vec3 albedo = _MatBridgeBrass + (fre + ref) * (_SunCol * _SunBrightness);
        //vec3 albedo = _MatBridgeBrass + fre * (2.0 * _SunCol * _SunBrightness) + ref * _MatBridgeBrassSpe;
        
        return mix(base_shadow, albedo, shadow);
    }


    
    vec3 col = vec3(0.0);
    col = 0.5 + 0.5*nor;
    //col = mix(vec3(0.0), col, shadow);

    return 1.2*col;
}

//--CLOUDS-------------------------------------------------------------------------------------------------------------------------
vec3 multipleOctaveScattering(in float extinction, in float mu, in float step_size)
{
    vec3 li = vec3(0.0);
    const float OCTAVES = 4.0;
    
    float a = 1.0;
    float b = 1.0;
    float c = 1.0;
    
    float phase;
    
    for (float i = 0.0; i<OCTAVES; i++)
    {
        phase = mix(henyeyGreenstein(-0.1 * c, mu), henyeyGreenstein(0.3 * c, mu), 0.7);
        li += b * phase * exp(-step_size * extinction * a);
        a*=0.25;
        b*=0.5;
        c*=0.5;
    }
    
    return li;
}

bool intersectClouds(in vec3 ro, in vec3 rd, out float t_near, out float isect_dist)
{
    // This is the bounding volume for the cloud.
    vec2 isect = intersectAABB(ro, rd, _CloudMinCorner, _CloudMaxCorner);
    
    if (insideAABB(ro, _CloudMinCorner, _CloudMaxCorner)) isect.x = 0.0001;
        
    t_near = isect.x;
    isect_dist = isect.y - isect.x;
    
    return isect.x > 0.0 && isect.x < isect.y;
}

// https://www.shadertoy.com/view/3sffzj
float getPerlinWorleyNoise(vec3 pos)
{
    // The cloud shape texture is an atlas of 6*6 tiles (36). 
    // Each tile is 32*32 with a 1 pixel wide boundary.
    // Per tile:		32 + 2 = 34.
    // Atlas width:	6 * 34 = 204.
    // The rest of the texture is black.
    // The 3D texture the atlas represents has dimensions 32 * 32 * 36.
    // The green channel is the data of the red channel shifted by one tile.
    // (tex.g is the data one level above tex.r). 
    // To get the necessary data only requires a single texture fetch.
    const float dataWidth = 204.0;
    const float tileRows = 6.0;
    const vec3 atlasDimensions = vec3(32.0, 32.0, 36.0);

    // Change from Y being height to Z being height.
    vec3 p = pos.xzy;

    // Pixel coordinates of point in the 3D data.
    vec3 coord = vec3(mod(p, atlasDimensions));
    float f = fract(coord.z);  
    float level = floor(coord.z);
    float tileY = floor(level/tileRows); 
    float tileX = level - tileY * tileRows;

    // The data coordinates are offset by the x and y tile, the two boundary cells 
    // between each tile pair and the initial boundary cell on the first row/column.
    vec2 offset = atlasDimensions.x * vec2(tileX, tileY) + 2.0 * vec2(tileX, tileY) + 1.0;
    vec2 pixel = coord.xy + offset;
    vec2 data = texture(iChannel0, mod(pixel, dataWidth)/iChannelResolution[0].xy).xy;
    return mix(data.x, data.y, f);
}

float mapCloud(vec3 p) { return 1.0; }

// Modified from https://www.shadertoy.com/view/3sffzj
float mapCloudDensity(in vec3 p, out float cloudHeight)
{  
    if(!insideAABB(p, _CloudMinCorner, _CloudMaxCorner)) return 0.0;
    
    cloudHeight = saturate((p.y - _CloudBot)/(_CloudTop-_CloudBot));
    float cloud = mapCloud(p);

    //If there are no clouds, exit early.
    if(cloud <= 0.0) return 0.0;

    //Sample texture which determines how high clouds reach.
    float height = cloud;
    
    //Round the top of the cloud. From "Real-time rendering of volumetric clouds". 
    cloud *= saturate(remap(cloudHeight, 0.8*height, height, 1.0, 0.0));
    //cloud *= saturate(remap(cloudHeight, 0.8*height, height, 1.0, 0.0));
    //cloud *= saturate(remap(cloudHeight, 0.0, 0.25 * (1.0-cloud), 0.0, 1.0))
    //      * saturate(remap(cloudHeight, 0.75 * height, height, 1.0, 0.0));


    //Animate main shape.
    p += vec3(_CloudShapeSpeed * iTime);
    
    //Get main shape noise, invert and scale it.
    float shape = 1.0-getPerlinWorleyNoise(_CloudShapeSize * p);
    shape *= _CloudShapeStrength;

    //Carve away density from cloud based on noise.
    cloud = saturate(remap(cloud, shape, 1.0, 0.0, 1.0));

    //Early exit from empty space
    if(cloud <= 0.0) return 0.0;
    
    //Animate details.
    p += vec3(_CloudDetailSpeed * iTime, 0.0, 0.5 * _CloudDetailSpeed * iTime);
    
    float detail = getPerlinWorleyNoise(_CloudDetailSize * p);
	detail *= _CloudDetailStrength;
    
	//Carve away detail based on the noise
	cloud = saturate(remap(cloud, detail, 1.0, 0.0, 1.0));
    return _CloudDensityMultiplier * cloud;
}

vec3 getCloudLi( in vec3 ro, in vec3 p, in float mu, in vec3 wi ) 
{
    // Sample density of volume between sample point p and the light source (wi = incident light dir).
    
    float isect_length = _CloudWidth*1.5;
    float tnear = 0.0;
    
    intersectClouds(p, wi, tnear, isect_length);

    // Similar approach here as with the raymarch into the cloud, we do a fixed number of steps out of the cloud.
    float ray_density = 0.0;
    
    float step_size = isect_length/float(LIGHT_RAY_STEPS);
    
    float h_placeholder;
    
    for (int i = 0; i<LIGHT_RAY_STEPS; i++)
    {
        vec3 q = p + float(i) * wi * step_size;
        ray_density += mix(1.0, 0.75, mu) * mapCloudDensity(q, h_placeholder);
    }
    
    // The light that reaches the point is then affected by Beer-Powder curve.
    vec3 beers_law = multipleOctaveScattering(ray_density, mu, step_size);
    vec3 powder = beers_law * 2.0 * (1.0 - (exp(-step_size * ray_density * 2.0 * _CloudSigmaE)));
    
    return mix(powder, beers_law, 0.5 + 0.5 * mu);
}

vec3 renderClouds( in vec3 ro, in vec3 rd, in float ray_offset, out vec3 ray_transmittance, out float start_t )
{
    /*
    REFERENCE MATERIALS:
    
    (VOLUMETRIC CLOUD) https://www.shadertoy.com/view/3sffzj
    
    */
    
    // We do a raycast through the cloud volume, at each point, figuring out the attenuation of light.
    
    ray_transmittance = vec3(1.0);
    vec3 col = vec3(0.0);
    
    float mu = dot(rd, _LightDir);
    // Forward and backwards scattering
    float phase_function = mix(henyeyGreenstein(-0.3, mu), henyeyGreenstein(0.3, mu), 0.7);
    
    // As we have an SDF for the cloud, we can start the raymarch at the cloud surface,
    //   instead of doing some sort of variable step size approach.
    float t = 0.0;
    float isect_length = 0.0; // length of cloud-ray intersection.
    bool isect = intersectClouds(ro, rd, t, isect_length);
    
    if (!isect) return col;
    
    // TODO: Early exit if no intersection.
    
    float step_size = isect_length / float(CAMERA_RAY_STEPS);
    t += step_size * ray_offset;
    start_t = -1.0;

    //const vec3 light_col = 100.0 * mix(vec3(0.65, 0.8, 1.0), _SunCol, 0.5);
    const vec3 light_col = vec3(0.65, 0.8, 1.0) * 100.0;
    
    // raymarch
    for(int i = 0; i < CAMERA_RAY_STEPS; i++)
    {    
        vec3 p = ro + t*rd;
        
        float height;
        float density = mapCloudDensity(p, height);
        vec3 sigma_s = density*_CloudSigmaS;
        vec3 sigma_e = density*_CloudSigmaE;
        
        if (density > 0.0)
        {
            if (start_t < 0.0) { start_t = t; }
            // Amount of light that reaches the sample point (Li = incident radiance, from rendering eq.)
            
            vec3 ambient = mix(4.0, 8.0, height)*_HorizonCol;
            //vec3 ambient = vec3(1.0);
            
            // Shadow casting onto the clouds is possible, but extremely expensive.
            //   This is even when approximating, as an exact shadow cast should compute the shadow on each light ray sample.
            #ifdef CLOUD_SHADOW_CAST
            float shadow = softShadow(p, _LightDir, 10.0);
            #else
            float shadow = 1.0;
            #endif

            vec3 li = ambient + 
                      shadow * light_col * phase_function * getCloudLi(ro, p, mu, _LightDir);
                        
            li *= sigma_s; // Which is proportional to cloud density at the point
            
            // Beer-Lambert
            vec3 transmittance = exp(-sigma_e * step_size);
            
            // integrate incident radiance for reflected radiance of the ray (pixel colour).
            col += ray_transmittance * (li - li * transmittance) / sigma_e;
            // update transmittance for the entire raycast
            ray_transmittance *= transmittance;
            
            // When the transmittance is almost zero, there is no more cloud occluding the sample pt, so early exit the raymarch.
            if (length(ray_transmittance) <= 0.001) { ray_transmittance = vec3(0.0); break; }
        }
        
        t += step_size;
    }
        
    return col;
}

//--CAMERA & RAYMARCHING-----------------------------------------------------------------------------------------------------------
vec3 render( in vec3 ro, in vec3 rd, in vec2 fragCoord ) 
{
    // --SKY----------------------------------------------------------------------
    vec3 col = sky(ro, rd);

    // --GEOMETRY-----------------------------------------------------------------
    vec3 tm = intersect( ro, rd );
    if( tm.x>0.0 )
    {
        col = shade(ro, rd, tm.x, tm.z);
    }
    col += sunSSSOutline(ro, rd, tm.y);
    
    // --CLOUDS-------------------------------------------------------------------
    vec3 cloud_transmittance = vec3(1.0);    
    float ray_offset = 0.0;
    
    #ifdef CLOUD_BLUE_NOISE
    // Wait for texture to load
    if(iChannelResolution[1].xy == vec2(1024.0))
    {
        float blue_noise = texture(iChannel1, fragCoord / 1024.0).r;
        ray_offset = fract(blue_noise + float(iFrame%32) * PHI);
        
    }
    #endif
    
    float cloud_t;
    //float min_t = tm.x; // track smallest t for fog
    vec3 cloud_col = 0.5 * renderClouds( ro, rd, ray_offset, cloud_transmittance, cloud_t ); 
    
    if ( tm.x < 0.0 || cloud_t < tm.x ) {
        //min_t = cloud_t;
        col = cloud_col + col * cloud_transmittance;
        //col = mix(cloud_col, col, cloud_transmittance);
    }
    
    //col = fog(col, min_t, rd);
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

    // SSAA
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
        #ifdef USE_DEBUG_CAMERA
        vec3 ta = DEBUG_CAMERA_TARGET;
        #else
        vec3 ta = CAMERA_TARGET;
        #endif
        
        vec2 m = iMouse.xy / iResolution.xy-.5;
        vec3 ro; // ray origin
        #ifdef USE_DEBUG_CAMERA
        if ( iMouse.x >= 0.0 ) 
        {
            ro = DEBUG_CAMERA_DIST*vec3(0.8, 0.3, 0.8);
            ro.yz *= rot(m.y*PI);
            ro.xz *= rot(-m.x*TAU*2.0);
            ro += ta;
        }
        #else
        //ro = ta + vec3( cos(an), -0.3, sin(an) );
        ro = ta + CAMERA_DIST*vec3( cos(an), 0.3*sin(2.0*an), sin(an) );
        #endif
        
        // camera-to-world transformation
        mat3 ca = setCamera( ro, ta, 0.0 );
        
        // ray direction
        float fl = 2.0;
        vec3 rd = ca * normalize( vec3(p,fl) );
        
        vec3 col = render( ro, rd, fragCoord );
        
	    tot += col;
    #if AA>1
    }
    tot /= float(AA*AA);
    #endif
    
    tot *= EXPOSURE;

    // fast dithering to remove banding artifacts, particularly from the sun.
    tot += sin(fragCoord.x*114.0)*sin(fragCoord.y*211.1)/512.0;
    
    // Divide through by HDR max so HDR values can be stored in 0->1 float buffer
    //   There is probably a loss of numerical precision, but oh well.
    tot /= HDR_MAX_COL;
        
    fragColor = vec4( tot, 1.0 );
}