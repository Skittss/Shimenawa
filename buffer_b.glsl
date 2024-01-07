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
    
      COMMON: Colour schemes, Rendering settings, SDF primitives, intersectors, and util funcs.
       
    BUFFER B: Perlin-Worley texture atlas generation for clouds. 
              (Source: https://www.shadertoy.com/view/3sffzj)
              
    BUFFER C: Scene rendering
   
    BUFFER D: HDR Bloom pass
    
       IMAGE: LDR Post processing (Tonemap, Gamma correction, etc.)
       
    --------------------------------------------------------------------------------------------
    
*/

//==============================================================================================
//  SOURCE: https://www.shadertoy.com/view/3sffzj)
//  Thank you alro for the noise implementation for the clouds here!
//==============================================================================================

//Create a Perlin-Worley texture atlas for cloud shape carving.
//Runs only once in the first frame.
//Based on https://github.com/sebh/TileableVolumeNoise/blob/master/main.cpp

//The atlas is a 6*6 grid of 32*32 tiles with a single layer of halo cells around each tile. 

//TODO: Assumes a size of at least 204 * 204. Make it work with any reasonable resolution.

#define PERLIN_WORLEY 0
#define WORLEY 1

#define SIZE 8.0
#define NUM_CELLS 2.0

vec3 modulo(vec3 m, float n) { return mod(mod(m, n) + n, n); }

// 5th order polynomial interpolation
vec3 fade(vec3 t) { return (t * t * t) * (t * (t * 6.0 - 15.0) + 10.0); }

// https://www.shadertoy.com/view/4djSRW
vec3 hash(vec3 p3)
{
    p3 = modulo(p3, SIZE);
    p3 = fract(p3 * vec3(0.1031, 0.1030, 0.0973));
    p3 += dot(p3, p3.yxz + 33.33);
    return 2.0 * fract((p3.xxy + p3.yxx) * p3.zyx) - 1.0;
}

float gradientNoise(vec3 p)
{
    vec3 i = floor(p);
    vec3 f = fract(p);
	
	vec3 u = fade(f);
    
    /*
    * For 1D, the gradient of slope g at vertex u has the form h(x) = g * (x - u), where u 
    * is an integer and g is in [-1, 1]. This is the equation for a line with slope g which 
    * intersects the x-axis at u.
    * For N dimensional noise, use dot product instead of multiplication, and do 
    * component-wise interpolation (for 3D, trilinear)
    */
    return mix( mix( mix( dot( hash(i + vec3(0.0,0.0,0.0)), f - vec3(0.0,0.0,0.0)), 
              dot( hash(i + vec3(1.0,0.0,0.0)), f - vec3(1.0,0.0,0.0)), u.x),
         mix( dot( hash(i + vec3(0.0,1.0,0.0)), f - vec3(0.0,1.0,0.0)), 
              dot( hash(i + vec3(1.0,1.0,0.0)), f - vec3(1.0,1.0,0.0)), u.x), u.y),
    mix( mix( dot( hash(i + vec3(0.0,0.0,1.0)), f - vec3(0.0,0.0,1.0)), 
              dot( hash(i + vec3(1.0,0.0,1.0)), f - vec3(1.0,0.0,1.0)), u.x),
         mix( dot( hash(i + vec3(0.0,1.0,1.0)), f - vec3(0.0,1.0,1.0)), 
              dot( hash(i + vec3(1.0,1.0,1.0)), f - vec3(1.0,1.0,1.0)), u.x), u.y), u.z );
}


float getPerlinNoise(vec3 pos, float frequency)
{
	//Compute the sum for each octave.
	float sum = 0.0;
	float weightSum = 0.0;
	float weight = 1.0;

	for(int oct = 0; oct < 3; oct++){

        vec3 p = pos * frequency;
        float val = 0.5 + 0.5 * gradientNoise(p);
        sum += val * weight;
        weightSum += weight;

        weight *= 0.5;
        frequency *= 2.0;
	}

	return saturate(sum / weightSum);
}

float worley(vec3 pos, float numCells)
{
	vec3 p = pos * numCells;
	float d = 1.0e10;
	for (int x = -1; x <= 1; x++){
		for (int y = -1; y <= 1; y++){
			for (int z = -1; z <= 1; z++){
                vec3 tp = floor(p) + vec3(x, y, z);
                tp = p - tp - (0.5 + 0.5 * hash(mod(tp, numCells)));
                d = min(d, dot(tp, tp));
            }
        }
    }
	return 1.0 - saturate(d);
}

//Return the 3D coordinate corresponding to the 2D atlas uv coordinate.
vec3 get3Dfrom2D(vec2 uv, float tileRows)
{
    vec2 tile = floor(uv);
    float z = floor(tileRows * tile.y + tile.x);
    return vec3(fract(uv), z);
}

float getTextureForPoint(vec3 p, int type)
{
	float res;
    if(type == PERLIN_WORLEY){
        
        //Perlin-Worley.
        float perlinNoise = getPerlinNoise(p, SIZE);
        res = perlinNoise;

        //Special weights from example code.
        float worley0 = worley(p, NUM_CELLS * 2.0);
        float worley1 = worley(p, NUM_CELLS * 8.0);
        float worley2 = worley(p, NUM_CELLS * 14.0);

        float worleyFBM = worley0 * 0.625 + worley1 * 0.25 + worley2 * 0.125;
        res = remap(perlinNoise, 0.0, 1.0, worleyFBM, 1.0);
        
	}else{

        //Worley
        float worley0 = worley(p, NUM_CELLS);
        float worley1 = worley(p, NUM_CELLS * 2.0);
        float worley2 = worley(p, NUM_CELLS * 4.0);
        float worley3 = worley(p, NUM_CELLS * 8.0);

        float FBM0 = worley0 * 0.625 + worley1 * 0.25 + worley2 * 0.125;
		float FBM1 = worley1 * 0.625 + worley2 * 0.25 + worley3 * 0.125;
		float FBM2 = worley2 * 0.75 + worley3 * 0.25;

        res = FBM0 * 0.625 + FBM1 * 0.25 + FBM2 * 0.125;
	}
    
	return res;
}

void mainImage(out vec4 fragColor, in vec2 fragCoord)
{
    if(iFrame < 1 || length(texelFetch(iChannel0, ivec2(0), 0).rgba) == 0.0){
        vec3 col = vec3(0);
        //32 with 1 pixel on either side.
        float tileSize = 34.0;
        float padWidth = 1.0;
        float coreSize = tileSize - 2.0 * padWidth;
        float tileRows = 6.0;
        float tileCount = tileRows * tileRows;
        vec2 tile = floor((fragCoord.xy - 0.5) / tileSize);

        bool padCell = false;
        if(mod(fragCoord.x, tileSize) == 0.5 || mod(fragCoord.x, tileSize) == tileSize - 0.5){
            padCell = true;
        }
        if(mod(fragCoord.y, tileSize) == 0.5 || mod(fragCoord.y, tileSize) == tileSize - 0.5){
            padCell = true;
        }

        bool startPadX = false;
        bool endPadX = false;
        bool startPadY = false;
        bool endPadY = false;

        if(fragCoord.x == tile.x * tileSize + 0.5){
            startPadX = true;
        }
        if(fragCoord.y == tile.y * tileSize + 0.5){
            startPadY = true;
        }
        if(fragCoord.x == (tile.x + 1.0) * tileSize - 0.5){
            endPadX = true;
        }
        if(fragCoord.y == (tile.y + 1.0) * tileSize - 0.5){
            endPadY = true;
        }

        vec2 padding = vec2(2.0 * padWidth) * tile;
        vec2 pixel;
        vec2 uv;
        
        if(!padCell){
            pixel = fragCoord.xy - padWidth - padding;
            uv = vec2(pixel.xy/coreSize);
        }else{
            pixel = fragCoord.xy - padWidth - padding;
            if(startPadX){
                pixel.x += coreSize;	
            }
            if(startPadY){
                pixel.y += coreSize;	
            }
            if(endPadX){
                pixel.x -= coreSize;	
            }
            if(endPadY){
                pixel.y -= coreSize;	
            }
            uv = vec2(pixel.xy/coreSize);
        }
        
        vec3 p_ = get3Dfrom2D(uv, tileRows);
        vec3 p = p_;
        p.z /= (tileRows*tileRows);

        // Get Perlin-Worley noise for level l
        float worleyPerlinNoise = getTextureForPoint(p, PERLIN_WORLEY);

        // Get Worley noise for level l
        float worleyNoise = getTextureForPoint(p, WORLEY);
        col.r = saturate(remap(worleyPerlinNoise, worleyNoise, 1.0, 0.0, 1.0));

        p_ = mod(p_ + 1.0, tileRows * tileRows);
        p = p_;
        p.z /= (tileRows*tileRows);

        // Get Perlin-Worley noise for level l+1
        worleyPerlinNoise = getTextureForPoint(p, PERLIN_WORLEY);

        // Get Worley noise for level l+1
        worleyNoise = getTextureForPoint(p, WORLEY);
        col.g = saturate(remap(worleyPerlinNoise, worleyNoise, 1.0, 0.0, 1.0));

        // Unused cells
        if(gl_FragCoord.x > tileRows * tileSize || gl_FragCoord.y > tileRows * tileSize){
            col = vec3(0);
        }
        
    	fragColor = vec4(col,1.0);
        
    } else {
        
    	fragColor = texelFetch(iChannel0, ivec2(fragCoord - 0.5), 0).rgba;
        
    }

}