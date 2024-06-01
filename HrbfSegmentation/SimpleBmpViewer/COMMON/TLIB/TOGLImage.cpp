#include "StdAfx.h"
#include "TOGLImage.h"


#include <atlimage.h>

#define MAX_VOLUME_SIZE 512 * 512 * 256

inline static double getLuminance( GLubyte* color )
{
	return ( 0.298912 * ( color[0] / 255.0 ) + 
		     0.586611 * ( color[1] / 255.0 ) + 
			 0.114478 * ( color[2] / 255.0 ) );
}

static void t_unbindTexture( const TOGL* ogl, unsigned int &texName)
{
	if( ogl != 0 && !ogl->IsDrawing() ) ogl->MakeOpenGLCurrent();
	if( texName != -1 && glIsTexture( texName ) ) glDeleteTextures(1,&texName ) ; 
	texName = -1;
	if( ogl != 0 && !ogl->IsDrawing() ) wglMakeCurrent( NULL, NULL);
}




//////////////////////////////////////////////////////////////////////////////
//volume volume volume////////////////////////////////////////////////////////
TOGL3DImage4::TOGL3DImage4(void)
{
	m_height = m_width = m_depth = 0;
	m_RGBA        =  0;
	m_textureName = -1;
	m_maxGradientValue = 0;
	m_DoInterpolation = true;
}

TOGL3DImage4::~TOGL3DImage4(void){ clear(0); }

void TOGL3DImage4::clear(const TOGL* ogl)
{
	t_unbindTexture( ogl, m_textureName);
	if( m_RGBA != 0) delete[] m_RGBA;

	m_height = m_width = m_depth = 0;
	m_RGBA        =  0;
	m_textureName = -1;
	m_maxGradientValue = 0;
}

void TOGL3DImage4::unbind(TOGL* ogl)
{
	t_unbindTexture( ogl, m_textureName);
}

void TOGL3DImage4::bind(TOGL* ogl)
{
	if( m_RGBA == 0 ) return;

	if( ogl != 0 && !ogl->IsDrawing() ) ogl->MakeOpenGLCurrent();
	if( m_textureName == -1 || !glIsTexture( m_textureName  ) )
	{
		//textureをoglに送る//
		glPixelStorei( GL_UNPACK_ALIGNMENT,4);
		glGenTextures( 1,&m_textureName ) ;

	
		glBindTexture( GL_TEXTURE_3D, m_textureName ) ;
		glTexImage3D( GL_TEXTURE_3D, 0 ,GL_RGBA, m_width, m_height, m_depth, 0, 
		                                GL_RGBA, GL_UNSIGNED_BYTE, m_RGBA);	
		glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE  );
		glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE  );
		glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE  );
		glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
		glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
		glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_BORDER);
		if( m_DoInterpolation ) {
			glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR) ;
			glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR) ;
		}else{
			glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
			glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		}


		//glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
	} else {	
		glBindTexture  (GL_TEXTURE_3D, m_textureName) ;//bindの後にsetする!!	
		if( m_DoInterpolation ) {
			glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR) ;
			glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR) ;
		}else{
			glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
			glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		}
	}
	if( ogl != 0 && !ogl->IsDrawing() ) wglMakeCurrent( NULL, NULL);
}


static inline int log_2( int n)
{
	if( n == 0 ) return -1;

	int i = 0;
	while( true )
	{
		if( n == 1 ) return i;
		n /= 2;
		++i;
	}
}



void TOGL3DImage4::allocateImage( unsigned int width, unsigned int  height, unsigned int depth, TOGL* ogl )
{
	clear(ogl);
	m_width  = width ;
	m_height = height;
	m_depth  = depth ;

	fprintf( stderr, "system allocates color volume image size is %d %d %d\n", width, height, depth);

	m_RGBA = new GLubyte[ width * height * depth * 4 ];
	memset( m_RGBA, 0 , sizeof( GLubyte ) * width * height * depth * 4 );
}

void TOGL3DImage4::allocateImage( const TOGL3DImage4 &src, TOGL* ogl )
{
	clear(ogl);
	m_width  = src.m_width ;
	m_height = src.m_height;
	m_depth  = src.m_depth ;
	fprintf( stderr, "system allocates color volume image size is %d %d %d\n", m_width, m_height, m_depth);
	m_RGBA = new GLubyte[ m_width * m_height * m_depth * 4 ];
	memcpy( m_RGBA, src.m_RGBA, sizeof(GLubyte) * m_width * m_height * m_depth * 4);
}

void TOGL3DImage4::allocateSphereVolume( int sizeR, TOGL *ogl )
{
	allocateImage( sizeR, sizeR, sizeR, ogl);
	TOGL3DImage1 tmp;
	tmp.allocateSphereVolume(sizeR, ogl);
	
	for( int i=0; i<sizeR*sizeR*sizeR; ++i)
	{
		m_RGBA[4*i+0] = m_RGBA[4*i+1] = 
		m_RGBA[4*i+2] = m_RGBA[4*i+3] = tmp.m_img[i];
	}
}


/* 元画像をサイズが同じになるようにresamplingする*/
void TOGL3DImage4::setVolumeColor_pack( const Image3D &img, CHANNEL_ID id)
{
	if( id == CHANNEL_ALL ){
		for( unsigned int z = 0; z < m_depth ; ++z)
		for( unsigned int y = 0; y < m_height; ++y)
		for( unsigned int x = 0; x < m_width ; ++x)
		{
			int    idx = 4 * ( x + y * m_width + z * m_width * m_height);
			double imgX = x / (double) m_width  * (double) img.sx; 
			double imgY = y / (double) m_height * (double) img.sy; 
			double imgZ = z / (double) m_depth  * (double) img.sz; 

			double c = img.img[ (int)imgZ ][ (int)imgY ][ (int)imgX ];
			if( c < 0 ) c = 0  ;
			if( c >255) c = 255;

			m_RGBA[idx + 0] = m_RGBA[idx + 1] = m_RGBA[idx + 2] = m_RGBA[idx + 3] = (GLubyte) c;
		}
	}
	else 
	{
		int channelI = ( id == CHANNEL_RED   ) ? 0 :
					   ( id == CHANNEL_GREEN ) ? 1 :
					   ( id == CHANNEL_BLUE  ) ? 2 : 3;
		for( unsigned int z = 0; z < m_depth ; ++z)
		for( unsigned int y = 0; y < m_height; ++y)
		for( unsigned int x = 0; x < m_width ; ++x)
		{
			int    idx = 4 * ( x + y * m_width + z * m_width * m_height);
			double imgX = x / (double) m_width  * (double) img.sx; 
			double imgY = y / (double) m_height * (double) img.sy; 
			double imgZ = z / (double) m_depth  * (double) img.sz; 
			
			double c = img.img[ (int)imgZ ][ (int)imgY ][ (int)imgX ];
			if( c < 0 ) c = 0  ;
			if( c >255) c = 255;
			m_RGBA[ idx + channelI ] = (GLubyte) c;
		}
	}
}


/* 元画像をサイズが同じになるようにresamplingする*/
void TOGL3DImage4::setVolumeColor_pack( const vector<TOGL2DImage*> &images, CHANNEL_ID id)
{
	if( images.size() ==  0){
		fprintf( stderr, "imput images contain error\n"); return;
	}
	
	int channelI =	( id == CHANNEL_RED   ) ? 0 :
					( id == CHANNEL_GREEN ) ? 1 :
					( id == CHANNEL_BLUE  ) ? 2 : 3;

	for( unsigned int z = 0; z < m_depth; ++z)
	{
		int    imgI0, imgI1;

		double t;
		double zVal = (double)z / (double)(m_depth-1);

		if     ( zVal  == 0  ){ imgI0 = imgI1 = 0                   ; t = 0;}
		else if( zVal  >= 1.0){ imgI0 = imgI1 = (int)images.size()-1; t = 0;}
		else{
			imgI0 = (int) ( zVal * (images.size()-1) );
			imgI1 = imgI0+1;
			t =  zVal * (images.size()-1)  -  imgI0;
		}

		const TOGL2DImage &img0 = *(images[ imgI0 ]);
		const TOGL2DImage &img1 = *(images[ imgI1 ]);

		//image0, image1, volumeすべてのwidth heightが異なることも想定
		double img0Xrate = (double) img0.m_width  / (double) m_width;
		double img0Yrate = (double) img0.m_height / (double) m_height;
		double img1Xrate = (double) img1.m_width  / (double) m_width;
		double img1Yrate = (double) img1.m_height / (double) m_height;

		for(unsigned int y = 0; y < m_height; ++y)
		for(unsigned int x = 0; x < m_width ; ++x)
		{
			int idx = 4 * (x + y * m_width + z * m_width * m_height);

			int idx2D0 = 4 * ( ((int)(x * img0Xrate)) + 
				               ((int)(y * img0Yrate)) * img0.m_width);
			int idx2D1 = 4 * ( ((int)(x * img1Xrate)) + 
				               ((int)(y * img1Yrate)) * img1.m_width);

			if( id == CHANNEL_ALL ) for( int i=0;i<4;++i)
				m_RGBA[idx +  i      ] = (GLubyte)( (1.0-t) * img0.m_RGBA[idx2D0 +i ] + 
												         t  * img1.m_RGBA[idx2D1 +i ] );
			else
				m_RGBA[idx + channelI] = (GLubyte)( (1.0-t) * img0.m_RGBA[idx2D0 +channelI ] + 
								                         t  * img1.m_RGBA[idx2D1 +channelI ] );
		}
	}
}

	
void TOGL3DImage4::setVolumeColor( const Image3D &img, CHANNEL_ID id)
{

	int channelI = ( id == CHANNEL_RED   ) ? 0 :
				   ( id == CHANNEL_GREEN ) ? 1 :
				   ( id == CHANNEL_BLUE  ) ? 2 : 3;
	//まずsetしたいchannelに0を入れる
	if( id == CHANNEL_ALL ) memset( m_RGBA, 0, sizeof(GLubyte) * 4 * m_depth * m_height * m_width );
	else
		for(unsigned int i= 0 ; i < m_width * m_height * m_depth * 4; i += 4) m_RGBA[ i + channelI ] = 0;


	if( id == CHANNEL_ALL ){
		for( unsigned int z = 0; z < m_depth && z < (unsigned int)img.sz; ++z) 
		for( unsigned int y = 0; y < m_height&& y < (unsigned int)img.sy; ++y) 
		for( unsigned int x = 0; x < m_width && x < (unsigned int)img.sx; ++x) 
		{
			double c = img.img[z][y][x]; 
			if( c < 0 ) c = 0; 
			if( c >255) c = 255;

			int    idx = 4 * ( x + y * m_width + z * m_width * m_height);
			m_RGBA[idx + 0] = m_RGBA[idx + 1] =  m_RGBA[idx + 2] = m_RGBA[idx + 3] =  (GLubyte) c;
		}
	}else {
		for( unsigned int z = 0; z < m_depth && z < (unsigned int)img.sz; ++z) 
		for( unsigned int y = 0; y < m_height&& y < (unsigned int)img.sy; ++y) 
		for( unsigned int x = 0; x < m_width && x < (unsigned int)img.sx; ++x) 
		{
			double c = img.img[z][y][x];   if( c < 0 ) c = 0;  if( c >255) c = 255;

			int    idx = 4 * ( x + y * m_width + z * m_width * m_height);
			m_RGBA[idx + channelI] =(GLubyte) c;
		}
	}
}




void TOGL3DImage4::setVolumeColor( const vector<TOGL2DImage*> &images, CHANNEL_ID id)
{
	if( images.size() ==  0){
		fprintf( stderr, "imput images contain error\n"); return;
	}
	
	int channelI =	( id == CHANNEL_RED   ) ? 0 :
					( id == CHANNEL_GREEN ) ? 1 :
					( id == CHANNEL_BLUE  ) ? 2 : 3;
	//まずsetしたいchannelに0を入れる
	if( id == CHANNEL_ALL ) memset( m_RGBA, 0, sizeof(GLubyte) * 4 * m_depth * m_height * m_width );
	else
		for(unsigned int i= 0 ; i < m_width * m_height * m_depth * 4; i += 4) m_RGBA[ i + channelI ] = 0;


	for( unsigned int z = 0; z < m_depth && z < (unsigned int) images.size(); ++z)
	{
		const TOGL2DImage &img = *images[ z ];
		
		for(unsigned int y = 0; y < m_height && y < img.m_height; ++y)
		for(unsigned int x = 0; x < m_width  && x < img.m_width ; ++x)
		{
			int idx   = 4 * (x + y * m_width + z * m_width * m_height);
			int idx2D = 4 * (y * img.m_width + x);

			if( id == CHANNEL_ALL ) for( int i=0;i<4;++i) m_RGBA[idx +  i      ] = img.m_RGBA[ idx2D + i        ]; 
			else                                          m_RGBA[idx + channelI] = img.m_RGBA[ idx2D + channelI ]; 
		}
	}
}




void TOGL3DImage4::convertToImage3D( Image3D& trgtImg, const TVector3 &cubeSize, CHANNEL_ID id ) const
{
	if( id == CHANNEL_ALL ) fprintf( stderr, "system cant convert all channgel of TOGL3DImage4 to simple image 3D!\n");

	//必要ならimgを解放
	if( trgtImg.img != 0 )
	{
		for(int i = 0; i < trgtImg.sz ; ++i)
		{
			for(int j =0 ; j < trgtImg.sy; ++j) delete[] trgtImg.img[i][j];
			delete[] trgtImg.img[i];
		}
    delete[] trgtImg.img;

	}

	// size
	trgtImg.sx = m_width ;
	trgtImg.sy = m_height;
	trgtImg.sz = m_depth ;

	trgtImg.px = cubeSize.data[0] / (double) m_width ;
	trgtImg.py = cubeSize.data[1] / (double) m_height;
	trgtImg.pz = cubeSize.data[2] / (double) m_depth ;

	trgtImg.maxgray   = 0;
	trgtImg.graylevel = 0;
	
	//img
    trgtImg.img = new double**[ trgtImg.sz ];
   
	for(int i = 0; i < trgtImg.sz; ++i)
	{
      trgtImg.img[i] = new double* [trgtImg.sy];
	  for(int j = 0; j < trgtImg.sy; ++j) trgtImg.img[i][j] = new double[trgtImg.sx]; 
	}

	int offset = ( id == CHANNEL_RED   ) ? 0 : 
				 ( id == CHANNEL_GREEN ) ? 1 : 
				 ( id == CHANNEL_BLUE  ) ? 2 : 3;
	for( int z = 0; z < (int) m_depth ; ++z )
	for( int y = 0; y < (int) m_height; ++y )
	for( int x = 0; x < (int) m_width ; ++x )
	{
		int idx4 = getVoxelIdx4( x, y, z);
		trgtImg.img[z][y][x] = (double) m_RGBA[idx4 + offset] / 255.0;
	}
}



void TOGL3DImage4::setIntensityToAlphaChannel()
{
	for( unsigned int z = 0; z < m_depth ; ++z)
	for( unsigned int y = 0; y < m_height; ++y)
	for( unsigned int x = 0; x < m_width ; ++x)
	{		
		int idx = 4 * (x + y * m_width + z * m_width * m_height);
		m_RGBA[idx+3] = max3( m_RGBA[idx+0], m_RGBA[idx+1], m_RGBA[idx+2]); //intensity (明暗度 HSVのV)
	}
}




/*
以下のような値が入る
a = length
2 * ( (r, g, b) - (0.5, 0.5, 0.5) )= (x, y, z)
void TOGL3DImage4::allocateGradientImage( const TOGL3DImage4 &trgtVol, TOGL* ogl)
{	
	allocateImage( trgtVol.m_width, trgtVol.m_height, trgtVol.m_depth, ogl);

	double rt_3 = sqrt( 3.0 );

	TVector3 grad;
	for(unsigned int z = 0; z < m_depth ; ++z )
	for(unsigned int y = 0; y < m_height; ++y )
	for(unsigned int x = 0; x < m_width ; ++x )
	{
		int pivIdx = 4 * ( x + y * m_width + z * m_width * m_height);

		if( x == 0          || y == 0           || z == 0          || 
			x == m_width -1 || y == m_height- 1 || z == m_depth - 1  )	
		{
			m_RGBA[pivIdx  ] = (GLubyte)127;
			m_RGBA[pivIdx+1] = (GLubyte)127;
			m_RGBA[pivIdx+2] = (GLubyte)127;
			m_RGBA[pivIdx+3] = (GLubyte)127;
			continue;
		}

		double lumi      = getLuminance( &trgtVol.m_RGBA[pivIdx                         ] );

		double lumi_preX = getLuminance( &trgtVol.m_RGBA[pivIdx - 4                     ] );
		double lumi_nexX = getLuminance( &trgtVol.m_RGBA[pivIdx + 4                     ] );
		
		double lumi_preY = getLuminance( &trgtVol.m_RGBA[pivIdx - 4 * m_width           ] );
		double lumi_nexY = getLuminance( &trgtVol.m_RGBA[pivIdx + 4 * m_width           ] );
		
		double lumi_preZ = getLuminance( &trgtVol.m_RGBA[pivIdx - 4 * m_width * m_height] );
		double lumi_nexZ = getLuminance( &trgtVol.m_RGBA[pivIdx + 4 * m_width * m_height] );

		grad.data[0] = lumi_nexX - lumi_preX;
		grad.data[1] = lumi_nexY - lumi_preY;
		grad.data[2] = lumi_nexZ - lumi_preZ;
		
		double l = grad.Length();
		if( l >  0.01 )
		{		
			m_RGBA[pivIdx + 0] = (GLubyte) (( 0.5 * grad.data[0] + 0.5 ) * 255.0 );
			m_RGBA[pivIdx + 1] = (GLubyte) (( 0.5 * grad.data[1] + 0.5 ) * 255.0 );
			m_RGBA[pivIdx + 2] = (GLubyte) (( 0.5 * grad.data[2] + 0.5 ) * 255.0 );
			m_RGBA[pivIdx + 3] = (GLubyte) ( min( 1, grad.Length() / rt_3 ) * 255.0 );
		}
		else
		{
			m_RGBA[pivIdx + 0] = 127;
			m_RGBA[pivIdx + 1] = 127;
			m_RGBA[pivIdx + 2] = 127;
			m_RGBA[pivIdx + 3] = 0;
		}
	}
}
*/

void TOGL3DImage4::allocateGradientImage( const TOGL3DImage1 &trgtVol, double gradMagCoef,  TOGL* ogl)
{	
	const int W = trgtVol.m_width ;
	const int H = trgtVol.m_height;
	const int D = trgtVol.m_depth ;
	const int WH = W*H;
	allocateImage( W, H, D, ogl);

	const double inv_255 = 1.0/255.0   ;

	TVector3 grad;
	int idx = 0;
	for(int z = 0; z < D; ++z       )
	for(int y = 0; y < H; ++y       )
	for(int x = 0; x < W; ++x, ++idx)
	{
		int idx4 = 4*idx;
		m_RGBA[idx4  ] = m_RGBA[idx4+1] = m_RGBA[idx4+2] = (GLubyte)127; m_RGBA[idx4+3] =0;
		if( x == 0   || y == 0    || z == 0   || x == W-1 || y == H- 1 || z == D- 1  ) continue;

		grad.data[0] = trgtVol.m_img[idx + 1  ]-  trgtVol.m_img[idx - 1 ]; 
		grad.data[1] = trgtVol.m_img[idx + W  ]-  trgtVol.m_img[idx - W ];
		grad.data[2] = trgtVol.m_img[idx + WH ]-  trgtVol.m_img[idx - WH];
		grad *= inv_255;
		
		double l = grad.Length();
		if( l >  0.01 )
		{		
			m_RGBA[idx4    ] = (GLubyte) (127.5 * grad.data[0] + 127.5);//(GLubyte) (( 0.5 * grad.data[0] + 0.5 ) * 255.0 );
			m_RGBA[idx4 + 1] = (GLubyte) (127.5 * grad.data[1] + 127.5);//(GLubyte) (( 0.5 * grad.data[1] + 0.5 ) * 255.0 );
			m_RGBA[idx4 + 2] = (GLubyte) (127.5 * grad.data[2] + 127.5);//(GLubyte) (( 0.5 * grad.data[2] + 0.5 ) * 255.0 );
			m_RGBA[idx4 + 3] = (GLubyte) ( min( 1, gradMagCoef * l) * 255.0 );
		}
	}
}







void TOGL3DImage4::simpleSmoothing(int n)
{
	clock_t t = clock();
	for( int i = 0; i < n; ++i) simpleSmoothing();
	fprintf( stderr, "smoothing %d   %f\n", n,(double)(clock()-t) / CLOCKS_PER_SEC);
}

void TOGL3DImage4::simpleSmoothing()
{
	
	GLubyte *tmpRGBA = new GLubyte[ m_width * m_height * m_depth * 4];

	int yOffset = 4 * m_width;
	int zOffset = 4 * m_width * m_height;

	int idx = 0;
	for(unsigned int z = 0; z < m_depth ; ++z )
	for(unsigned int y = 0; y < m_height; ++y )
	for(unsigned int x = 0; x < m_width ; ++x, idx += 4 )
	{
		if( x == 0          || y == 0            || z == 0          || 
			x == m_width -1 || y == m_height - 1 || z == m_depth -1  )	
		{
			memcpy( &tmpRGBA[ idx ], &m_RGBA[ idx ], sizeof( GLubyte ) * 4 );
		}
		else
			smoothing( &tmpRGBA[ idx],   &m_RGBA[idx		  ],
										 &m_RGBA[idx + 4	  ], &m_RGBA[idx - 4      ],
										 &m_RGBA[idx + yOffset], &m_RGBA[idx - yOffset],
										 &m_RGBA[idx + zOffset], &m_RGBA[idx - zOffset]);
	}

	swap( tmpRGBA, m_RGBA );
	delete[] tmpRGBA;
}


void TOGL3DImage4::halfSize( TOGL *ogl )
{
	unsigned int newWidth  = m_width  / 2;
	unsigned int newHeight = m_height / 2;
	unsigned int newDepth  = m_depth  / 2;
	GLubyte *newRGBA = new GLubyte[ newWidth * newHeight * newDepth * 4];

	for(unsigned int z = 0; z < newDepth  ; ++z )
	for(unsigned int y = 0; y < newHeight ; ++y )
	for(unsigned int x = 0; x < newWidth  ; ++x )
	{
		int newIdx = 4 * (     x +      y * newWidth  +   z * newHeight * newWidth);
		int oldIdx = 4 * (   2*x +    2*y * m_width   + 2*z * m_width * m_height);
		for( int i = 0; i < 4; ++i)	
			newRGBA[newIdx+i] = m_RGBA[oldIdx + i];
	}

	allocateImage( newWidth, newHeight, newDepth, ogl);
	memcpy( m_RGBA, newRGBA, sizeof( GLubyte ) * m_width * m_height * m_depth * 4);
	delete[] newRGBA;
}



void TOGL3DImage4::resize( int newWidth, int newHeight, int newDepth, TOGL *ogl)
{
	if( !isAllocated() ) return;
	GLubyte *newRGBA = new GLubyte[ newWidth * newHeight * newDepth * 4];

	for(int z = 0; z < newDepth  ; ++z )
	for(int y = 0; y < newHeight ; ++y )
	for(int x = 0; x < newWidth  ; ++x )
	{	
		int oldX = (int)( (double) x / (double) newWidth  * (double)m_width  );
		int oldY = (int)( (double) y / (double) newHeight * (double)m_height );
		int oldZ = (int)( (double) z / (double) newDepth  * (double)m_depth  );

		int newIdx = 4 * (    x +     y * newWidth  +   z  * newHeight * newWidth);
		int oldIdx = 4 * ( oldX +  oldY * m_width   + oldZ * m_height  * m_width );
	
		for( int i = 0; i < 4; ++i)	newRGBA[newIdx+i] = m_RGBA[oldIdx + i];
	}

	allocateImage( newWidth, newHeight, newDepth, ogl);
	memcpy( m_RGBA, newRGBA, sizeof( GLubyte ) * m_width * m_height * m_depth * 4);
	delete[] newRGBA;
}




void TOGL3DImage4::clipImage( int minX, int minY, int minZ, int maxX, int maxY, int maxZ, TOGL *ogl )
{
	if( !isAllocated() ) return;
	
	int newWidth = maxX - minX + 1;
	int newHeight= maxY - minY + 1;
	int newDepth = maxZ - minZ + 1;
	fprintf( stderr, " clip image new size is = %d %d %d\n", newWidth, newHeight, newDepth );
	if(     minX <      0       ||    minY   <      0       || minZ<0             || 
		newWidth > (int)m_width || newHeight > (int)m_height|| newDepth > (int)m_depth) return;


	GLubyte *newRGBA = new GLubyte[ newWidth * newHeight * newDepth * 4];

	for(int z = minZ; z <= maxZ  ; ++z )
	for(int y = minY; y <= maxY  ; ++y )
	for(int x = minX; x <= maxX  ; ++x )
	{	
		int newX = x - minX;
		int newY = y - minY;
		int newZ = z - minZ;
		int oldIdx = 4 * (   x +     y * m_width  +    z  * m_width * m_height  );
		int newIdx = 4 * (newX +  newY * newWidth + newZ  * newWidth * newHeight);
	
		for( int i = 0; i < 4; ++i)	newRGBA[newIdx+i] = m_RGBA[oldIdx + i];
	}

	allocateImage( newWidth, newHeight, newDepth, ogl);
	memcpy( m_RGBA, newRGBA, sizeof( GLubyte ) * m_width * m_height * m_depth * 4);
	delete[] newRGBA;
}


void TOGL3DImage4::normalizeGradientImage()
{
	//obtain max Gradient Magnitude
	TVector3 grad;
	m_maxGradientValue = 0;

	const int vidSize4= 4*(m_depth * m_height * m_width);
	for(int i = 0; i < vidSize4; i += 4 )
	{
		grad.data[0] = ( m_RGBA[i+ 0] / 255.0  - 0.5) * 2;
		grad.data[1] = ( m_RGBA[i+ 1] / 255.0  - 0.5) * 2;
		grad.data[2] = ( m_RGBA[i+ 2] / 255.0  - 0.5) * 2;
		m_maxGradientValue = max( m_maxGradientValue , grad.Length() );
	}

	for(int i = 0; i < vidSize4; i += 4 )
	{
		grad.data[0] = ( m_RGBA[i+ 0] / 255.0  - 0.5) * 2;
		grad.data[1] = ( m_RGBA[i+ 1] / 255.0  - 0.5) * 2;
		grad.data[2] = ( m_RGBA[i+ 2] / 255.0  - 0.5) * 2;
		double l = grad.Length();
		if( l  >  0.01)
		{	
			grad.Normalize_Self();
			m_RGBA[i + 0] = (GLubyte) min(255.0, max(0.0, (( 0.5 * grad.data[0] + 0.5 ) * 255.0 ) ));
			m_RGBA[i + 1] = (GLubyte) min(255.0, max(0.0, (( 0.5 * grad.data[1] + 0.5 ) * 255.0 ) ));
			m_RGBA[i + 2] = (GLubyte) min(255.0, max(0.0, (( 0.5 * grad.data[2] + 0.5 ) * 255.0 ) ));
			m_RGBA[i + 3] = (GLubyte) ( l / m_maxGradientValue * 255.0) ;
		}
		else
		{
			m_RGBA[i + 0] = 127;
			m_RGBA[i + 1] = 127;
			m_RGBA[i + 2] = 127;
			m_RGBA[i + 3] = 0;
		}
	}
}



void TOGL3DImage4::gainContrast( GLubyte bottom, GLubyte top)
{
	if( bottom >= top ) return;

	for(unsigned int z = 0; z < m_depth ; ++z )
	for(unsigned int y = 0; y < m_height; ++y )
	for(unsigned int x = 0; x < m_width ; ++x )
	{
		int idx = 4 * ( x + y * m_width + z * m_width * m_height);

		for( int rgb = 0; rgb < 4; ++rgb)
		{
			if     ( m_RGBA[ idx + rgb ] < bottom ) m_RGBA[idx + rgb] = 0  ;
			else if( m_RGBA[ idx + rgb ] > top    ) m_RGBA[idx + rgb] = 255;
			else
			{
				double val = (255.0 /( top - bottom )) * ( m_RGBA[idx + rgb] - bottom);
				if( val > 255 ) val = 255; 
				m_RGBA[ idx + rgb ] = (GLubyte) val; 
			}
		}

	}
}





/////////////////////////////////////////////////////////////////////////////////////
//2D image //////////////////////////////////////////////////////////////////////////
TOGL2DImage::~TOGL2DImage(){
	clear(0);
}

TOGL2DImage::TOGL2DImage ()
{
	m_RGBA = 0;
	m_width = m_height = 0;
	m_textureName = -1;
	m_DoInterpolation = true;

}

void TOGL2DImage::clear(TOGL* ogl)
{
	t_unbindTexture( ogl, m_textureName );

	if( m_RGBA != 0 ) delete[] m_RGBA ; 
	m_RGBA = 0;	
	m_width = m_height = 0;
}


void TOGL2DImage::unbind(TOGL* ogl){
	t_unbindTexture( ogl, m_textureName );
}


void TOGL2DImage::bind(TOGL* ogl)
{
	if( ogl != 0 && !ogl->IsDrawing() ) ogl->MakeOpenGLCurrent();
	if( m_textureName == -1 || !glIsTexture( m_textureName  ) )
	{
		//fprintf( stderr, "ogl 2D image bind  %d --> ", m_textureName );
		//textureをoglに送る//
		glPixelStorei( GL_UNPACK_ALIGNMENT,4);
		glGenTextures( 1,&m_textureName ) ;

	
		glBindTexture( GL_TEXTURE_2D, m_textureName ) ;
		glTexImage2D( GL_TEXTURE_2D, 0 ,GL_RGBA, m_width, m_height, 0, GL_RGBA, GL_UNSIGNED_BYTE, m_RGBA);		

		//fprintf( stderr, " %d \n", m_textureName );


		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		if( m_DoInterpolation ) {
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR) ;
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR) ;
		}else{
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		}
	} else {
		

		glBindTexture  (GL_TEXTURE_2D, m_textureName) ;
		if( m_DoInterpolation ) {
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR) ;
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR) ;
		}else{
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		}

	}
	if( ogl != 0 && !ogl->IsDrawing() ) wglMakeCurrent( NULL, NULL);
}



void TOGL2DImage::allocateImage( unsigned int width, unsigned int height, TOGL* ogl)
{
	clear(ogl);
	m_width  = width ;
	m_height = height;
	m_RGBA = new GLubyte[ width * height * 4 ];
	memset( m_RGBA, 0 , sizeof( GLubyte ) * width * height * 4 );
}


void TOGL2DImage::allocateImage( const TOGL2DImage &src , TOGL* ogl)
{
	clear(ogl);
	m_width  = src.m_width ;
	m_height = src.m_height;
	m_RGBA = new GLubyte[ m_width * m_height * 4 ];
	memcpy( m_RGBA, src.m_RGBA, sizeof( GLubyte ) * m_width * m_height * 4 );
}


/*-------------------------------------------------------------------------------
CImageを利用した画像読み込み
BMP、GIF、JPEG、PNG、および TIFF に対応している(らしい)

const char* fname :file name
bool &inverted    :負なら原点が画像の左下隅 (画素自体は左上から並ぶので逆に向ける必要がある)
TOGL* ogl         :既にビデオメモリに乗っている場合は、unbindする必要がある。その対象のoglを指定する
--------------------------------------------------------------------------------*/
bool TOGL2DImage::allocateFromFile( const char *fname, bool &inverted, TOGL* ogl)
{
	CImage pic;
	if( !SUCCEEDED( pic.Load( _T(fname) ) ) ) return false;
	
	//メージのピッチ(行方向のバイト数)を返します。
	//戻り値が負の値の場合、ビットマップは左下隅を起点とする逆方向 (下から上) の DIB です。
	//戻り値が正の値の場合、ビットマップは左上隅を起点とする順方向 (上から下の向き) の DIB です
	int pitch   = pic.GetPitch();
	if( pitch < 0) { inverted = true ; pitch *= -1; }
	else             inverted = false;
	int W       = pic.GetWidth ();
	int H       = pic.GetHeight();
	int byteNum = pitch / W;

	allocateImage( W, H, ogl);

	if( byteNum == 1 )
	{
		for( int y=0, i=0; y<H; ++y     )
		for( int x=0     ; x<W; ++x, ++i)
		{
			//RGB(24bit)に対するアドレスをbyteにキャストしてるので, R G B順の一番下位のBを指すものになる
			byte *c = (byte*)pic.GetPixelAddress( x, y );
			m_RGBA[ i*4 + 3] = 128;
			m_RGBA[ i*4 + 2] = *c ; 
			m_RGBA[ i*4 + 1] = *c ; 
			m_RGBA[ i*4 + 0] = *c ; 
		}
	}
	else if( byteNum == 3 ) //24bit color 
	{
		for( int y=0, i=0; y<H; ++y     )
		for( int x=0     ; x<W; ++x, ++i)
		{
			//RGB(24bit)に対するアドレスをbyteにキャストしてるので, R G B順の一番下位のBを指すものになる
			byte *c = (byte*)pic.GetPixelAddress( x, y );
			m_RGBA[ i*4 + 3] = 128;
			m_RGBA[ i*4 + 2] = *c ; ++c;
			m_RGBA[ i*4 + 1] = *c ; ++c;
			m_RGBA[ i*4 + 0] = *c ; 
		}
	}
	else //その他(ピクセルごとにビデオメモリをポーリングする遅い実装)
	{
		for( int y=0, i=0; y<H; ++y     )
		for( int x=0     ; x<W; ++x, ++i)
		{
			COLORREF c = pic.GetPixel( x, y );
			m_RGBA[ i*4 + 3] = 128;
			m_RGBA[ i*4 + 0] = GetRValue(c);
			m_RGBA[ i*4 + 1] = GetGValue(c);
			m_RGBA[ i*4 + 2] = GetBValue(c);
		}
	}	
	return true;
}

bool TOGL2DImage::saveAsFile( const char *fname, int flg_BmpJpgPngTiff)
{
	CImage pic;
	pic.Create(m_width, m_height, 24, 0 );
	for( unsigned int y = 0, i=0; y < m_height; ++y     )
	for( unsigned int x = 0     ; x < m_width ; ++x, ++i)
	{
		byte *c = (byte*)pic.GetPixelAddress( x, y );
		*c = m_RGBA[ i*4 + 2]; ++c;
		*c = m_RGBA[ i*4 + 1]; ++c;
		*c = m_RGBA[ i*4 + 0]; 
	}
	if(      flg_BmpJpgPngTiff == 0 )  pic.Save( _T( fname ), Gdiplus::ImageFormatBMP  );
	else if( flg_BmpJpgPngTiff == 1 )  pic.Save( _T( fname ), Gdiplus::ImageFormatJPEG );
	else if( flg_BmpJpgPngTiff == 2 )  pic.Save( _T( fname ), Gdiplus::ImageFormatPNG  );
	else if( flg_BmpJpgPngTiff == 3 )  pic.Save( _T( fname ), Gdiplus::ImageFormatTIFF );	
	return true;
}






void TOGL2DImage::setAllAlphaCh( byte alpha)
{
	if( !isAllocated() ) return;
	for( int i = 0; i < (int)m_height; ++i )
	for( int j = 0; j < (int)m_width ; ++j ) m_RGBA[ 4*(i * m_width + j) + 3 ] = alpha;
}

void TOGL2DImage::gaussianFilter33()
{
	GLubyte *tmp = new GLubyte[ 4*m_width * m_height];
	static double f[3][3]= { {1,2,1}, {2,4,2}, {1,2,1} };
	for( unsigned int y=0,i=0; y<m_height; ++y     )
	for( unsigned int x=0    ; x<m_width ; ++x, ++i)
	{
		double r = 0, g = 0, b = 0, a = 0, s=0;
		for( int yy=-1; yy<2; ++yy) if( 0 <= y+yy && y+yy < m_height )
		for( int xx=-1; xx<2; ++xx) if( 0 <= x+xx && x+xx < m_width  )
		{
			int idx  = i + xx + yy*m_width;
			int idx4 = 4*idx;
			s += f[yy+1][xx+1];
			r += f[yy+1][xx+1] * m_RGBA[ idx4+0];
			g += f[yy+1][xx+1] * m_RGBA[ idx4+1];
			b += f[yy+1][xx+1] * m_RGBA[ idx4+2];
			a += f[yy+1][xx+1] * m_RGBA[ idx4+3];
		}
		tmp[ 4*i +0 ] = (GLubyte) (r/s);
		tmp[ 4*i +1 ] = (GLubyte) (g/s);
		tmp[ 4*i +2 ] = (GLubyte) (b/s);
		tmp[ 4*i +3 ] = (GLubyte) (a/s);
	}

	memcpy( m_RGBA, tmp, sizeof( GLubyte) * m_width * m_height *4 );
	delete[] tmp;

}


/////////////////////////////////////////////////////////////////////////////////////
//1D image 1D image 1D image////////////////////////////////////////////////////////

TOGL1DImage::~TOGL1DImage(){
	clear(0);
}

TOGL1DImage::TOGL1DImage ()
{
	m_width       =  0;
	m_RGBA        =  0;
	m_textureName = -1;
	m_DoInterpolation = true;
}


void TOGL1DImage::clear(TOGL *ogl)
{
	t_unbindTexture( ogl, m_textureName );

	if( m_RGBA != 0 ) delete[] m_RGBA ; 
	m_RGBA  = 0;
	m_width = 0;
}


void TOGL1DImage::unbind(TOGL *ogl){
	t_unbindTexture( ogl, m_textureName );
}


void TOGL1DImage::bind(TOGL *ogl)
{
	if( m_RGBA == 0 ) return;

	if( ogl != 0 && !ogl->IsDrawing() ) ogl->MakeOpenGLCurrent();
	if( m_textureName == -1 || !glIsTexture( m_textureName  ) )
	{
		fprintf( stderr, "Here TOGL Image 1 genTex\n");
		//textureをoglに送る//
		glPixelStorei( GL_UNPACK_ALIGNMENT,4);
		glGenTextures( 1, &m_textureName ) ;

		glBindTexture( GL_TEXTURE_1D, m_textureName ) ;
		glTexImage1D( GL_TEXTURE_1D, 0 ,GL_RGBA, m_width, 0, GL_RGBA, GL_UNSIGNED_BYTE, m_RGBA);		

		glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		if( m_DoInterpolation ) {
			glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR) ;
			glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR) ;
		}else{
			glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
			glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		}
		//glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
	} else {

		glBindTexture(GL_TEXTURE_1D, m_textureName) ;//より後のものが影響する
		
		if( m_DoInterpolation )
		{
			glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR) ;
			glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR) ;
		}else{
			glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
			glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		}
	}
	if( ogl != 0 && !ogl->IsDrawing() ) wglMakeCurrent( NULL, NULL);
}




void TOGL1DImage::allocateImage( unsigned int width, TOGL* ogl)
{
	clear( ogl );
	m_width  = width ;
	m_RGBA = new GLubyte[ width * 4 ];
	memset( m_RGBA, 0 , sizeof( GLubyte ) * width * 4 );
	m_DoInterpolation = true;
}


void TOGL1DImage::allocateHeuColorBar( TOGL* ogl )
{
	clear(ogl);
	m_width = 256;
	m_RGBA  = new GLubyte[ m_width * 4 ];
	memset( m_RGBA, 0 , sizeof( GLubyte ) * m_width * 4 );

	int texSize   = 256;
	int texSize_2 = texSize / 2;
	for (int i = 0; i < texSize; ++i) 
	{
		int r = 512 - 1024 * i / texSize;
		int g = 512 - 512 * abs(i - texSize_2) / texSize_2;
		int b = -512 + 1024 * i / texSize;
		GLubyte r_ = (GLubyte)max(min(r, 255), 0);
		GLubyte g_ = (GLubyte)max(min(g, 255), 0);
		GLubyte b_ = (GLubyte)max(min(b, 255), 0);
		m_RGBA[4 * i    ] = r_;
		m_RGBA[4 * i + 1] = g_;
		m_RGBA[4 * i + 2] = b_;
		m_RGBA[4 * i + 3] = 255;
	}
}



///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
//Signel channel image//////////////////////////////////////////////////////////////////////////////////////
TOGL3DImage1::~TOGL3DImage1(void){
	clear(0);
}
TOGL3DImage1::TOGL3DImage1(void)
{
	m_height = m_width = m_depth = 0;
	m_img    =  0;
	m_textureName = -1;
	m_DoInterpolation = true;
}



void TOGL3DImage1::clear(const TOGL *ogl)
{
	t_unbindTexture( ogl, m_textureName );
	if( m_img != 0 ) delete[] m_img; 
	m_height = m_width = m_depth = 0;
	m_img    =  0;
	m_textureName = -1;	
}

void TOGL3DImage1::unbind(const  TOGL *ogl ){
	t_unbindTexture( ogl, m_textureName );
}

void TOGL3DImage1::bind(const TOGL *ogl )
{
	if( m_img == 0 ) return;

	if( ogl != 0 && !ogl->IsDrawing() ) ogl->MakeOpenGLCurrent();
	if( m_textureName == -1 || !glIsTexture( m_textureName  ) )
	{
		//fprintf( stderr, "bindTexture %d(これが-1でないときはおかしい) %d\n ", m_textureName, glIsTexture( m_textureName  ) );

		//textureをoglに送る//
		glPixelStorei( GL_UNPACK_ALIGNMENT , 1);
		glGenTextures( 1,&m_textureName ) ;

		glBindTexture( GL_TEXTURE_3D, m_textureName ) ;
		glTexImage3D( GL_TEXTURE_3D, 0 ,GL_LUMINANCE8, m_width, m_height, m_depth, 0, 
		                                GL_LUMINANCE, GL_UNSIGNED_BYTE, m_img );
		//glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		//glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		//glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
		glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
		glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_BORDER);
				
		if( m_DoInterpolation ){
			glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR) ;
			glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR) ;
		}else{
			glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_NEAREST) ;
			glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_NEAREST) ;
		}
		//glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
	} else {
		glBindTexture  (GL_TEXTURE_3D, m_textureName) ;
		if( m_DoInterpolation ){
			glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR) ;
			glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR) ;
		}else{
			glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_NEAREST) ;
			glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_NEAREST) ;
		}
	}
	if( ogl != 0 && !ogl->IsDrawing() ) wglMakeCurrent( NULL, NULL);
}




//////image allocation////////////////////////////////////////////////////////////////////////////
void TOGL3DImage1::allocateImage( unsigned int width, unsigned int  height, unsigned int depth, const TOGL *ogl )
{
	clear(ogl);
	if( width == 0 || height == 0 || depth == 0) return;
	m_width  = width ;
	m_height = height;
	m_depth  = depth ;
	fprintf( stderr, "system allocates intensity volume image size is %d %d %d\n", width, height, depth);

	m_img = new GLubyte[ width * height * depth ];
	memset( m_img, 0 , sizeof( GLubyte ) * width * height * depth );
}



void TOGL3DImage1::allocateImage( unsigned int width, unsigned int  height, unsigned int depth ,
	                           const float* img, const double minV, const double maxV, const TOGL* ogl)
{
	clear(ogl);
	if( width == 0 || height == 0 || depth == 0) return;
	m_width  = width ;
	m_height = height;
	m_depth  = depth ;
	fprintf( stderr, "system allocates intensity volume image size is %d %d %d\n", width, height, depth);

	m_img = new GLubyte[ width * height * depth ];

	const int WH = width * height;

#pragma omp parallel for
	for( int z = 0; z < (int) depth; ++z){
		for( int y = 0; y < (int) height; ++y)
		for( int x = 0; x < (int) width; ++x)
		{
			int idx = x + y * width + z *WH;
			double d = max(0, min( 1.0,  (img[idx]-minV) / (maxV-minV) ) );
			m_img[ idx ] = (byte)( 255.0 * d );
		}
	}
}



void TOGL3DImage1::allocateImage( const TOGL3DImage1 &srcImg, const TOGL *ogl  ){
	if( !srcImg.isAllocated() ) return;
	allocateImage( srcImg.m_width, srcImg.m_height, srcImg.m_depth, ogl);
	copy( srcImg ); 
}

void TOGL3DImage1::allocateImageFromAlphaChannel(const TOGL3DImage4 &vol, const TOGL *ogl){
	allocateImage( vol.m_width, vol.m_height, vol.m_depth, ogl);
	for( int i=0; i<(int) (vol.m_width * vol.m_height * vol.m_depth) ; ++i) m_img[i] = vol.m_RGBA[i*4+3];
}


void TOGL3DImage1::allocateSphereVolume( int sizeR, const TOGL *ogl  )
{
	clear(ogl);
	allocateImage( sizeR, sizeR, sizeR, ogl);
	double cubeR   = 10;
	double sphereR =  3;

	TVector3 pos, rad;
	TVector3 center( cubeR/2.0 , cubeR/2.0, cubeR/2.0);
	for( int z = 0; z < sizeR; ++z)
	for( int y = 0; y < sizeR; ++y)
	for( int x = 0; x < sizeR; ++x)
	{
		pos.Set( x / (double) sizeR * cubeR, 
			     y / (double) sizeR * cubeR, 
				 z / (double) sizeR * cubeR);
		pos -= center;

		// val =  - ( x-r ) ^ 2 + 255
		double r = pos.Length();
		double val =  - 10 * ( r - sphereR ) * ( r - sphereR ) * 
			                 ( r - sphereR ) * ( r - sphereR ) * 
			                 ( r - sphereR ) * ( r - sphereR ) + 255;
		
		if( val < 0) val = 0;
		m_img[ x + y*sizeR + z*sizeR*sizeR ] = (GLubyte) val; 
	}
}

/* もとの実装(2012/4/7更新)
static double sblCoef = 1/18.0;
static double sblX[3][3][3]={ {{-1,  0, +1},
	                           {-3,  0, +3},
							   {-1,  0, +1}},
						      {{-3,  0,  3},
							   {-6,  0,  6},
							   {-3,  0,  3}},
							  {{-1,  0,  1},
							   {-3,  0,  3},
							   {-1,  0,  1}} };

static double sblY[3][3][3]={ {{-1, -3, -1},
	                           { 0,  0,  0},
							   { 1,  3,  1}},
						      {{-3, -6, -3},
							   { 0,  0,  0},
							   { 3,  6,  3}},
							  {{-1, -3, -1},
							   { 0,  0,  0},
							   { 1,  3,  1}} };

static double sblZ[3][3][3]={ {{-1, -3, -1},
	                           {-1, -6, -1},
							   {-1, -3, -1}},
						      {{ 0,  0,  0},
							   { 0,  0,  0},
							   { 0,  0,  0}},
							  {{ 1,  3,  1},
							   { 1,  6,  1},
							   { 1,  3,  1}} };
*/
static double sblCoef      = 1/16.0/2.0;
static double sblCoef_edge = 1/12.0/2.0;
static double sblX[3][3][3]={ {{-1,  0, +1},
	                           {-2,  0, +2},
							   {-1,  0, +1}},   {{-2,  0,  2},
												 {-4,  0,  4},
												 {-2,  0,  2}},  {{-1,  0,  1},
																  {-2,  0,  2},
																  {-1,  0,  1}} };
static double sblY[3][3][3]={ {{-1, -2, -1},
	                           { 0,  0,  0},
							   { 1,  2,  1}},
						      {{-2, -4, -2},
							   { 0,  0,  0},
							   { 2,  4,  2}},
							  {{-1, -2, -1},
							   { 0,  0,  0},
							   { 1,  2,  1}} };

static double sblZ[3][3][3]={ {{-1, -2, -1},
	                           {-2, -4, -2},
							   {-1, -2, -1}},
						      {{ 0,  0,  0},
							   { 0,  0,  0},
							   { 0,  0,  0}},
							  {{ 1,  2,  1},
							   { 2,  4,  2},
							   { 1,  2,  1}} };

void TOGL3DImage1::allocateGradMagImage(const TOGL3DImage1 &vol, const TOGL* ogl, double scale )
{
	if( !vol.isAllocated() ) return;
	const int W  = vol.m_width ;
	const int H  = vol.m_height;
	const int D  = vol.m_depth ;
	const int WH = W*H;
	this->allocateImage( W, H, D, ogl);
	fprintf( stderr, "gradient magnitude !!!! %d %d %d\n", m_width,m_height,m_depth );
	
	//id[z][y][x]
	const int id[3][3][3]={ {{-1-W-WH, -W-WH, 1-W-WH},
						     {-1  -WH,   -WH, 1  -WH},
						     {-1+W-WH, +W-WH, 1+W-WH}},
						    {{-1-W   , -W   , 1-W   },
						     {-1     ,  0   , 1     },
						     {-1+W   , +W   , 1+W   }},
						    {{-1-W+WH, -W+WH, 1-W+WH},
						     {-1  +WH,   +WH, 1  +WH},
						     {-1+W+WH, +W+WH, 1+W+WH}}};
#pragma omp parallel for
	for(int z = 0; z < D; ++z  ){
		for(int y = 0; y < H; ++y  )
		for(int x = 0; x < W; ++x  )
		{
			int idx = x + y*W + z*WH;
			double gx=0, gy=0, gz=0;
			for( int zz=-1;zz<2;++zz) if( 0<=z+zz && z+zz<D)
			for( int yy=-1;yy<2;++yy) if( 0<=y+yy && y+yy<H)
			for( int xx=-1;xx<2;++xx) if( 0<=x+xx && x+xx<W)
			{
				double c = vol[ idx +id[zz+1][yy+1][xx+1] ];
				gx += sblX[zz+1][yy+1][xx+1] * c;
				gy += sblY[zz+1][yy+1][xx+1] * c;
				gz += sblZ[zz+1][yy+1][xx+1] * c;
			}

			//端の場合はその方向の微分は計算しない
			double coef = sblCoef * scale;
			if( x==0 || x==W-1) { gx = 0; coef = sblCoef_edge;}
			if( y==0 || y==H-1) { gy = 0; coef = sblCoef_edge;}
			if( z==0 || z==D-1) { gz = 0; coef = sblCoef_edge;}
			gx *= coef;
			gy *= coef;
			gz *= coef;
			m_img[ idx ] = (GLubyte) ( min( 255, sqrt( gx*gx + gy*gy + gz*gz) ) );
		}
	}
}
void TOGL3DImage1::allocateGradMagImage( const TOGL3DImage1 &r, 
									     const TOGL3DImage1 &g, const TOGL* ogl, double scale )
{
	const int W = r.m_width ;
	const int H = r.m_height;
	const int D = r.m_depth ;
	const int WH = W*H;
	this->allocateImage( W, H, D, ogl);
	
	//id[z][y][x]
	int id[3][3][3]={ {{-1-W-WH, -W-WH, 1-W-WH},
	                   {-1  -WH,   -WH, 1  -WH},
	                   {-1+W-WH, +W-WH, 1+W-WH}},
					  {{-1-W   , -W   , 1-W   },
	                   {-1     ,  0   , 1     },
	                   {-1+W   , +W   , 1+W   }},
					  {{-1-W+WH, -W+WH, 1-W+WH},
	                   {-1  +WH,   +WH, 1  +WH},
	                   {-1+W+WH, +W+WH, 1+W+WH}}};

	TVector3 grad;
	for(int z = 0; z < D; ++z  )
	for(int y = 0; y < H; ++y  )
	for(int x = 0; x < W; ++x  )
	{
		int idx = x + y*W + z*WH;

		grad.Set(0,0,0);
		for( int zz=-1;zz<2;++zz) if( 0<=z+zz && z+zz<D)
		for( int yy=-1;yy<2;++yy) if( 0<=y+yy && y+yy<H)
		for( int xx=-1;xx<2;++xx) if( 0<=x+xx && x+xx<W)
		{
			double c = r[ idx +id[zz+1][yy+1][xx+1] ] + g[ idx +id[zz+1][yy+1][xx+1] ];
			grad.data[0] +=  sblX[zz+1][yy+1][xx+1] * c;
			grad.data[1] +=  sblY[zz+1][yy+1][xx+1] * c;
			grad.data[2] +=  sblZ[zz+1][yy+1][xx+1] * c;
		}

		//橋の場合はその方向の微分は計算しない
		double coef = sblCoef* scale;
		if( x==0 || x==W-1) { grad.data[0] = 0; coef = sblCoef_edge;}
		if( y==0 || y==H-1) { grad.data[1] = 0; coef = sblCoef_edge;}
		if( z==0 || z==D-1) { grad.data[2] = 0; coef = sblCoef_edge;}
		grad *= coef/2.0;
		m_img[ idx ] = (GLubyte) ( min( 255, grad.Length())  );
	}
}

void TOGL3DImage1::allocateGradMagImage( const TOGL3DImage1 &r, 
									     const TOGL3DImage1 &g, 
									     const TOGL3DImage1 &b, const TOGL* ogl, double scale )
{
	const int W = r.m_width ;
	const int H = r.m_height;
	const int D = r.m_depth ;
	const int WH = W*H;
	this->allocateImage( W, H, D, ogl);
	
	//id[z][y][x]
	int id[3][3][3]={ {{-1-W-WH, -W-WH, 1-W-WH},
	                   {-1  -WH,   -WH, 1  -WH},
	                   {-1+W-WH, +W-WH, 1+W-WH}},
					  {{-1-W   , -W   , 1-W   },
	                   {-1     ,  0   , 1     },
	                   {-1+W   , +W   , 1+W   }},
					  {{-1-W+WH, -W+WH, 1-W+WH},
	                   {-1  +WH,   +WH, 1  +WH},
	                   {-1+W+WH, +W+WH, 1+W+WH}}};

	TVector3 grad;
	for(int z = 0; z < D; ++z  )
	for(int y = 0; y < H; ++y  )
	for(int x = 0; x < W; ++x  )
	{
		int idx = x + y*W + z*WH;

		grad.Set(0,0,0);
		for( int zz=-1;zz<2;++zz) if( 0<=z+zz && z+zz<D)
		for( int yy=-1;yy<2;++yy) if( 0<=y+yy && y+yy<H)
		for( int xx=-1;xx<2;++xx) if( 0<=x+xx && x+xx<W)
		{
			double c = r[ idx +id[zz+1][yy+1][xx+1] ] + g[ idx +id[zz+1][yy+1][xx+1] ] + b[ idx +id[zz+1][yy+1][xx+1] ];
			grad.data[0] += sblX[zz+1][yy+1][xx+1] * c;
			grad.data[1] += sblY[zz+1][yy+1][xx+1] * c;
			grad.data[2] += sblZ[zz+1][yy+1][xx+1] * c;
		}

		//端の場合はその方向の微分は計算しない
		double coef = sblCoef * scale;
		if( x==0 || x==W-1) { grad.data[0] = 0; coef = sblCoef_edge;}
		if( y==0 || y==H-1) { grad.data[1] = 0; coef = sblCoef_edge;}
		if( z==0 || z==D-1) { grad.data[2] = 0; coef = sblCoef_edge;}
		grad *= coef/3.0;
		m_img[ idx ] = (GLubyte) min( 255, grad.Length());
	}
}

void TOGL3DImage1::calcGrad33( const int x, const int y, const int z, TVector3 &grad ) const
{
	const int W = m_width ;
	const int H = m_height;
	const int D = m_depth ;
	const int WH = W*H;
	grad.Set(0,0,0);
	int idx = x + y*W + z*WH;

	for( int zz=-1;zz<2;++zz) if( 0<=z+zz && z+zz<D)
	for( int yy=-1;yy<2;++yy) if( 0<=y+yy && y+yy<H)
	for( int xx=-1;xx<2;++xx) if( 0<=x+xx && x+xx<W)
	{
		double c = m_img[ idx + xx + yy*W + zz*WH];
		grad.data[0] += sblX[zz+1][yy+1][xx+1] * c;
		grad.data[1] += sblY[zz+1][yy+1][xx+1] * c;
		grad.data[2] += sblZ[zz+1][yy+1][xx+1] * c;
	}
	double coef = sblCoef;
	if( x==0 || x==W-1) { grad.data[0] = 0; coef = sblCoef_edge;}
	if( y==0 || y==H-1) { grad.data[1] = 0; coef = sblCoef_edge;}
	if( z==0 || z==D-1) { grad.data[2] = 0; coef = sblCoef_edge;}
	grad *= coef;

}

//x方向 sbl5[y][z][x]
//y方向 sbl5[z][x][y]
//z方向 sbl5[x][y][z]
static double sbl5[5][5][5] = 
{ { { -1, -2, -4, -2, -1}, 
	{ -2, -6, -8, -6, -2}, 
	{ -4, -8,-16, -8, -4} , 
	{ -2, -6, -8, -6, -2} , 
	{ -1, -2, -4, -2, -1} },

  { { -2, -4, -8, -4, -2}, 
	{ -4,-12,-16,-12, -4}, 
	{ -8,-16,-32,-16, -8} , 
	{ -4,-12,-16,-12, -4} , 
	{ -2, -4, -8, -4, -2} },
								 
	{ { 0,0,0,0,0}, 
      { 0,0,0,0,0}, 
      { 0,0,0,0,0}, 
      { 0,0,0,0,0}, 
      { 0,0,0,0,0} }, 
 {  {  2,  4,  8,  4,  2} ,
    {  4, 12, 16, 12,  4} , 
	{  8, 16, 32, 16,  8} , 
	{  4, 12, 16, 12,  4} , 
	{  2,  4,  8,  4,  2} },

  { {  1,  2,  4,  2,  1},        
	{  2,  6,  8,  6,  2}, 
	{  4,  8, 16,  8,  4} , 
	{  2,  6,  8,  6,  2} , 
    {  1,  2,  4,  2,  1} } };
void TOGL3DImage1::calcGrad55( const int x, const int y, const int z, TVector3 &grad ) const
{
	int W = m_width, H = m_height, D = m_depth, WH = W*H;
	int idx = x + y*m_width + z * WH;

	grad.Set( 0,0,0 );
	for( int zz=-2;zz<3;++zz) if( 0<=z+zz && z+zz<D)
	for( int yy=-2;yy<3;++yy) if( 0<=y+yy && y+yy<H)
	for( int xx=-2;xx<3;++xx) if( 0<=x+xx && x+xx<W)
	{
		double c = m_img[ idx + xx + yy*W + zz*WH];
		grad[0] += sbl5[xx+2][zz+2][yy+2] * c;
		grad[1] += sbl5[yy+2][xx+2][zz+2] * c;
		grad[2] += sbl5[zz+2][yy+2][xx+2] * c;
	}
	grad /= (324 * 2.0);
	if( x<2 || W-2 < x ) {grad[0] = 0; fprintf( stderr, "a");}
	if( y<2 || H-2 < y ) {grad[1] = 0; fprintf( stderr, "b");}
	if( z<2 || D-2 < z ) {grad[2] = 0; fprintf( stderr, "c");}
}


static double gSig = 1.0 / 2;
static double gauss33[3][3][3] = 
{ { {exp(-3/gSig), exp(-2/gSig), exp(-3/gSig)}, 
    {exp(-2/gSig), exp(-1/gSig), exp(-2/gSig)}, 
    {exp(-3/gSig), exp(-2/gSig), exp(-3/gSig)} },

  { {exp(-2/gSig), exp(-1/gSig), exp(-2/gSig)}, 
    {exp(-1/gSig), exp(-0/gSig), exp(-1/gSig)}, 
    {exp(-2/gSig), exp(-1/gSig), exp(-2/gSig)} },

  { {exp(-3/gSig), exp(-2/gSig), exp(-3/gSig)}, 
    {exp(-2/gSig), exp(-1/gSig), exp(-2/gSig)}, 
    {exp(-3/gSig), exp(-2/gSig), exp(-3/gSig)} },
  };

static double gauss55[5][5][5] = 
{ {{exp(-12/gSig), exp(-9/gSig), exp(-8/gSig), exp(-9/gSig),exp(-12/gSig)}, 
    {exp(-9/gSig), exp(-6/gSig), exp(-5/gSig), exp(-6/gSig), exp(-9/gSig)}, 
    {exp(-8/gSig), exp(-5/gSig), exp(-4/gSig), exp(-5/gSig), exp(-8/gSig)}, 
    {exp(-9/gSig), exp(-6/gSig), exp(-5/gSig), exp(-6/gSig), exp(-9/gSig)}, 
   {exp(-12/gSig), exp(-9/gSig), exp(-8/gSig), exp(-9/gSig),exp(-12/gSig)} },
	
  {	{exp(-9/gSig), exp(-6/gSig), exp(-5/gSig), exp(-6/gSig), exp(-9/gSig)}, 
    {exp(-6/gSig), exp(-3/gSig), exp(-2/gSig), exp(-3/gSig), exp(-6/gSig)}, 
    {exp(-5/gSig), exp(-2/gSig), exp(-1/gSig), exp(-2/gSig), exp(-5/gSig)}, 
    {exp(-6/gSig), exp(-3/gSig), exp(-2/gSig), exp(-3/gSig), exp(-6/gSig)}, 
    {exp(-9/gSig), exp(-6/gSig), exp(-5/gSig), exp(-6/gSig), exp(-9/gSig)} },

  {	{exp(-8/gSig), exp(-5/gSig), exp(-4/gSig), exp(-5/gSig), exp(-8/gSig)}, 
    {exp(-5/gSig), exp(-2/gSig), exp(-1/gSig), exp(-2/gSig), exp(-5/gSig)}, 
    {exp(-4/gSig), exp(-1/gSig), exp(-0/gSig), exp(-1/gSig), exp(-4/gSig)}, 
    {exp(-5/gSig), exp(-2/gSig), exp(-1/gSig), exp(-2/gSig), exp(-5/gSig)}, 
    {exp(-8/gSig), exp(-5/gSig), exp(-4/gSig), exp(-5/gSig), exp(-8/gSig)} },

  {	{exp(-9/gSig), exp(-6/gSig), exp(-5/gSig), exp(-6/gSig), exp(-9/gSig)}, 
    {exp(-6/gSig), exp(-3/gSig), exp(-2/gSig), exp(-3/gSig), exp(-6/gSig)}, 
    {exp(-5/gSig), exp(-2/gSig), exp(-1/gSig), exp(-2/gSig), exp(-5/gSig)}, 
    {exp(-6/gSig), exp(-3/gSig), exp(-2/gSig), exp(-3/gSig), exp(-6/gSig)}, 
    {exp(-9/gSig), exp(-6/gSig), exp(-5/gSig), exp(-6/gSig), exp(-9/gSig)} },
 
  {{exp(-12/gSig), exp(-9/gSig), exp(-8/gSig), exp(-9/gSig),exp(-12/gSig)}, 
    {exp(-9/gSig), exp(-6/gSig), exp(-5/gSig), exp(-6/gSig), exp(-9/gSig)}, 
    {exp(-8/gSig), exp(-5/gSig), exp(-4/gSig), exp(-5/gSig), exp(-8/gSig)}, 
    {exp(-9/gSig), exp(-6/gSig), exp(-5/gSig), exp(-6/gSig), exp(-9/gSig)}, 
   {exp(-12/gSig), exp(-9/gSig), exp(-8/gSig), exp(-9/gSig),exp(-12/gSig)} } };


void TOGL3DImage1::gaussianSmoothing33(int n)
{
	GLubyte *tmpImg = new GLubyte[ m_width * m_height * m_depth ];
	const int W = m_width, H = m_height, D = m_depth, WH = W*H;

	for( int kkk=0; kkk<n; ++kkk)
	{
#pragma omp parallel for
		for( int z = 0; z < D ; ++z )
		for( int y = 0; y < H ; ++y )
		for( int x = 0; x < W ; ++x )
		{
			int idx = x + y*W + z * WH;
			double coef = 0, v = 0;
			for( int zz=-1;zz<2;++zz) if( 0<=z+zz && z+zz<D)
			for( int yy=-1;yy<2;++yy) if( 0<=y+yy && y+yy<H)
			for( int xx=-1;xx<2;++xx) if( 0<=x+xx && x+xx<W)
			{
				coef += gauss33[zz+1][yy+1][xx+1];
				v    += gauss33[zz+1][yy+1][xx+1] * m_img[ idx + xx + yy*W + zz*WH];
			}
			tmpImg[idx] = (byte)(v / coef);
		}
		memcpy( m_img, tmpImg, sizeof( byte ) * W*H*D);
	}
	delete[] tmpImg;
}

void TOGL3DImage1::gaussianSmoothing55(int n)
{
	GLubyte *tmpImg = new GLubyte[ m_width * m_height * m_depth ];
	const int W = m_width, H = m_height, D = m_depth, WH = W*H;

	for( int kkk=0; kkk<n; ++kkk)
	{
#pragma omp parallel for
		for( int z = 0; z < D ; ++z )
		{
			for( int y = 0; y < H ; ++y )
			for( int x = 0; x < W ; ++x )
			{
				int idx = x + y*W + z * WH;
				double coef = 0, v = 0;
				for( int zz=-2;zz<3;++zz) if( 0<=z+zz && z+zz<D)
				for( int yy=-2;yy<3;++yy) if( 0<=y+yy && y+yy<H)
				for( int xx=-2;xx<3;++xx) if( 0<=x+xx && x+xx<W)
				{
					coef += gauss55[zz+2][yy+2][xx+2];
					v    += gauss55[zz+2][yy+2][xx+2] * m_img[ idx + xx + yy*W + zz*WH];
				}
				tmpImg[idx] = (byte)(v / coef);
			}
		}
		memcpy( m_img, tmpImg, sizeof( byte ) * W*H*D);
		fprintf( stderr, "smoothing done %d / %d\n", kkk, n);
	}
	delete[] tmpImg;
}


void TOGL3DImage1::gaussianSmoothing55(int n, int sX,int sY,int sZ, int eX,int eY,int eZ)
{
	GLubyte *tmpImg = new GLubyte[ m_width * m_height * m_depth ];
	const int W = m_width, H = m_height, D = m_depth, WH = W*H;

	for( int kkk=0; kkk<n; ++kkk)
	{
		memcpy( tmpImg, m_img, sizeof( byte ) * W*H*D);
#pragma omp parallel for
		for( int z = sZ; z < eZ ; ++z )
		{
			for( int y = sY; y < eY ; ++y )
			for( int x = sX; x < eX ; ++x )
			{
				int idx = x + y*W + z * WH;
				double coef = 0, v = 0;
				for( int zz=-2;zz<3;++zz) if( 0<=z+zz && z+zz<D)
				for( int yy=-2;yy<3;++yy) if( 0<=y+yy && y+yy<H)
				for( int xx=-2;xx<3;++xx) if( 0<=x+xx && x+xx<W)
				{
					coef += gauss55[zz+2][yy+2][xx+2];
					v    += gauss55[zz+2][yy+2][xx+2] * m_img[ idx + xx + yy*W + zz*WH];
				}
				tmpImg[idx] = (byte)(v / coef);
			}
		}
		memcpy( m_img, tmpImg, sizeof( byte ) * W*H*D);
		fprintf( stderr, "smoothing done %d / %d\n", kkk, n);
	}
	delete[] tmpImg;
}




void TOGL3DImage1::simpleSmoothing(int n)
{
	clock_t t = clock();
	for( int i = 0; i < n; ++i) simpleSmoothing();
	fprintf( stderr, "smoothing %d   %f\n", n,(double)(clock()-t) / CLOCKS_PER_SEC);
}

void TOGL3DImage1::simpleSmoothing()
{
	GLubyte *tmpImg = new GLubyte[ m_width * m_height * m_depth ];

	int yOffset = m_width;
	int zOffset = m_width * m_height;

	int idx = 0;
	for(unsigned int z = 0; z < m_depth ; ++z )
	for(unsigned int y = 0; y < m_height; ++y )
	for(unsigned int x = 0; x < m_width ; ++x, ++idx )
	{
		if( x == 0          || y == 0            || z == 0          || 
			x == m_width -1 || y == m_height - 1 || z == m_depth -1  )	
		{
			tmpImg[ idx ] = m_img[ idx ];
		}
		else
		{
			double d = 0.1 * (4.0 * m_img[idx] + m_img[idx + 1	    ] +  m_img[idx - 1      ]+
											     m_img[idx + yOffset] +  m_img[idx - yOffset]+
											     m_img[idx + zOffset] +  m_img[idx - zOffset] );				
			tmpImg[ idx ]= ( d > 255.0 )? 255 : (GLubyte) d;
		}
	}

	swap( tmpImg, m_img);
	delete[] tmpImg;
}






void TOGL3DImage1::gainContrast( GLubyte bottom, GLubyte top)
{
	if( bottom >= top ) return;

	for( int idx = 0; idx<getVolIdxSize();++idx)
	{
		if     ( m_img[ idx ] < bottom ) m_img[ idx ] = 0  ;
		else if( m_img[ idx ] > top    ) m_img[ idx ] = 255;
		else{
			double val = (255.0 /( top - bottom )) * ( m_img[ idx ] - bottom);
			if( val > 255 ) val = 255; 
			m_img[ idx ] = (GLubyte) val; 
		}
	}
}

void TOGL3DImage1::resize(int newWidth, int newHeight, int newDepth, TOGL *ogl)
{
	if( !isAllocated() ) return;
	GLubyte *newImg = new GLubyte[ newWidth * newHeight * newDepth ];

	for(int z = 0; z < newDepth  ; ++z )
	for(int y = 0; y < newHeight ; ++y )
	for(int x = 0; x < newWidth  ; ++x )
	{	
		int oldX = (int)( (double) x / (double) newWidth  * (double)m_width  );
		int oldY = (int)( (double) y / (double) newHeight * (double)m_height );
		int oldZ = (int)( (double) z / (double) newDepth  * (double)m_depth  );

		int newIdx =    x +    y * newWidth  +   z  * newHeight * newWidth ;
		int oldIdx = oldX + oldY * m_width   + oldZ * m_height  * m_width  ;
		newImg[ newIdx ] = m_img[ oldIdx ];
	}

	allocateImage( newWidth, newHeight, newDepth, ogl );
	copy( newImg );
	delete[] newImg;
}

void TOGL3DImage1::setImage(const Image3D &img)
{
	setImgZero();
	for( unsigned int z = 0; z < m_depth  && z < (unsigned int)img.sz; ++z) 
	for( unsigned int y = 0; y < m_height && y < (unsigned int)img.sy; ++y) 
	for( unsigned int x = 0; x < m_width  && x < (unsigned int)img.sx; ++x) 
	{
		double c = img.img[z][y][x]; 
		if( c < 0 ) c = 0  ; 
		if( c >255) c = 255;
		m_img[ ( x + y * m_width + z * m_width * m_height) ] = (byte)c;
	}
}



void TOGL3DImage1::setImage_fitSize(const TOGL3DImage1 &img)
{
	setImgZero();
	double xRate = (double)(img.m_width -1)/ (double) (m_width -1);
	double yRate = (double)(img.m_height-1)/ (double) (m_height-1);
	double zRate = (double)(img.m_depth -1)/ (double) (m_depth -1);
	
	int idx = 0;
	for( unsigned int z = 0; z < m_depth  ; ++z) 
	for( unsigned int y = 0; y < m_height ; ++y) 
	for( unsigned int x = 0; x < m_width  ; ++x, ++idx) 
		m_img[idx] = img.m_img[ img.getVoxelIdx( (int)( x * xRate), 
			                                     (int)( y * yRate), 
												 (int)( z * zRate)) ];
}


void TOGL3DImage1::setImageIntens(const Image3D &Rimg, const Image3D &Gimg, const Image3D &Bimg)
{
	setImgZero();
	for( unsigned int z = 0; z < m_depth  && z < (unsigned int)Rimg.sz; ++z) 
	for( unsigned int y = 0; y < m_height && y < (unsigned int)Rimg.sy; ++y) 
	for( unsigned int x = 0; x < m_width  && x < (unsigned int)Rimg.sx; ++x) 
	{
		double c = (Rimg.img[z][y][x] + Gimg.img[z][y][x] + Bimg.img[z][y][x]) / 3.0; 
		if( c < 0 ) c = 0  ; 
		if( c >255) c = 255;
		m_img[ ( x + y * m_width + z * m_width * m_height) ] = (byte)c;
	}
}

//static byte calcLuminance(byte R, byte G, byte B)
//{
//	return (byte)( 0.3 * R + 0.59 * G + 0.11 * B);
//}
void TOGL3DImage1::setImageIntens(const vector<TOGL2DImage> &images)
{
	if( images.size() ==  0){ fprintf( stderr, "imput images contain error\n"); return; }

	setImgZero();
	for( unsigned int z = 0; z < m_depth && z < (unsigned int) images.size(); ++z)
	{
		const TOGL2DImage &img = images[ z ];
		
		for(unsigned int y = 0; y < m_height && y < img.m_height; ++y)
		for(unsigned int x = 0; x < m_width  && x < img.m_width ; ++x)
		{
			int idx   =   (x + y * m_width + z * m_width * m_height);
			int idx2D = 4*(x + y * img.m_width                     );
			//m_img[ idx ] = max( max( img.m_RGBA[idx2D+0], img.m_RGBA[idx2D+1]), img.m_RGBA[idx2D+2]);
			//m_img[ idx ] = calcLuminance( img.m_RGBA[idx2D+0], img.m_RGBA[idx2D+1], img.m_RGBA[idx2D+2]);
			m_img[ idx ] = (byte)( (img.m_RGBA[idx2D+0] + img.m_RGBA[idx2D+1] + img.m_RGBA[idx2D+2]) / 3.0);			
		}
	}
}


void TOGL3DImage1::setImage(const vector<TOGL2DImage > &images, int channel/*0:r 1:g 2:b*/)
{
	if( images.size() ==  0){ fprintf( stderr, "imput images contain error\n"); return; }

	setImgZero();
	for( unsigned int z = 0; z < m_depth && z < (unsigned int) images.size(); ++z)
	{
		const TOGL2DImage &img = images[ z ];
		for(unsigned int y = 0; y < m_height && y < img.m_height; ++y)
		for(unsigned int x = 0; x < m_width  && x < img.m_width ; ++x)
		{
			int idx   =   (x + y * m_width + z * m_width * m_height);
			int idx2D = 4*(x + y * img.m_width                     );
			m_img[ idx ] = img.m_RGBA[idx2D+channel];
		}
	}
}

//trgtIという値を持つ連結画素のラべリング
//初期値1で最大ラベルIDを返す
//int *lavel既にallocated
int TOGL3DImage1::labeling26ConnectedRegion(byte trgtI, int *label) const 
{
	struct voxelInfo{
		int x,y,z,idx;
		voxelInfo(){}
		voxelInfo(int _x,int _y,int _z,int _idx){ x=_x; y=_y; z=_z; idx=_idx;}
	};

	const int W = m_width ;
	const int H = m_height;
	const int D = m_depth , WH= W*H, WHD = W*H*D;
	memset( label, 0, sizeof( int ) * WHD );

	int searchI = 0;
	int labelId = 0;
	while( true )
	{
		//pivot検索
		int pivIdx = -1;
		for( ; searchI <WHD; ++searchI) if( label[searchI] == 0 && m_img[searchI] == trgtI ){
			pivIdx = searchI;
			break;
		}
		if( pivIdx == -1 ) break;

		//region growing 26近傍利用 fronteerを利用する
		labelId++;
		fprintf( stderr, "start region growing!! %d\n", labelId);
		list<voxelInfo> frontier; 
		int _x,_y,_z; getXYZIdx( pivIdx, _x,_y,_z );

		frontier.push_back( voxelInfo( _x,_y,_z, pivIdx ) );
		while( !frontier.empty() )
		{
			voxelInfo v = frontier.front(); frontier.pop_front();
			for( int z = v.z-1; z < v.z+2; ++z) if(0 <= z && z < D) 
			for( int y = v.y-1; y < v.y+2; ++y) if(0 <= y && y < H)
			for( int x = v.x-1; x < v.x+2; ++x) if(0 <= x && x < W)
			{
				int idx = z * WH + y * W + x;
				if( label[ idx ] != 0 || m_img[idx] != trgtI) continue;
				label[ idx ] = labelId;
				frontier.push_back( voxelInfo( x,y,z, idx) );
			}
		}
	}
	return labelId;
}


double TOGL3DImage1::triLinearSampling( const TVector3 &cSize, const TVector3 &pos ) const
{
	const int W = m_width ;
	const int H = m_height;
	const int D = m_depth , WH= W*H, WHD = W*H*D;
	double px = cSize[0] / W;
	double py = cSize[1] / H;
	double pz = cSize[2] / D;
	//voxelの中心に値が埋め込まれているとする
	int xi = (int)( (pos[0] - 0.5* px) / px); double tx = (pos[0] - (xi +0.5)*px) / px;
	int yi = (int)( (pos[1] - 0.5* py) / py); double ty = (pos[1] - (yi +0.5)*py) / py;
	int zi = (int)( (pos[2] - 0.5* pz) / pz); double tz = (pos[2] - (zi +0.5)*pz) / pz;
	if( xi < 0 || W-1 < xi) return 0;
	if( yi < 0 || H-1 < yi) return 0;
	if( zi < 0 || D-1 < zi) return 0;

	const int idx = xi + yi * W + zi * WH;
	if( xi== 0 || W-1== xi) return m_img[idx];
	if( yi== 0 || H-1== yi) return m_img[idx];;
	if( zi== 0 || D-1== zi) return m_img[idx];;

	double v[8] = { m_img[idx   ], m_img[idx+1   ], m_img[idx+1+W   ], m_img[idx+W   ], 
		            m_img[idx+WH], m_img[idx+1+WH], m_img[idx+1+W+WH], m_img[idx+W+WH] };

	double p0132 = (1-tx)*(1-ty) * v[0] +  tx *(1-ty) * v[1] +   
		           (1-tx)*  ty   * v[3] +  tx *  ty   * v[2]; 
	double p4576 = (1-tx)*(1-ty) * v[4] +  tx *(1-ty) * v[5] +   
		           (1-tx)*  ty   * v[7] +  tx *  ty   * v[6]; 
	return (1-tz) * p0132 + tz * p4576;
}



///////////////////////////////////////////////////////////////////
//calc Histram for intensity, gardient magnitude and multiDimension//
void TI_calcLogScaleHistgramesOfIntensity( const TOGL3DImage4 &volume, 
										   const TOGL3DImage4 &gradient, 
										   vector<       double > &volHist ,
										   vector<       double > &gradHist,
										   vector<vector<double>> &MultHist )
{
	const int binSize = 256;
	volHist .clear(); volHist .resize( binSize, 0 );
	gradHist.clear(); gradHist.resize( binSize, 0 );
	MultHist.clear(); MultHist.resize( binSize, vector<double>( binSize, 0) );

	//gradientのサイズがvolumeとは違う可能性を考慮
	double gXrate = (gradient.m_width - 1) /( volume.m_width - 1), 
		   gYrate = (gradient.m_height- 1) / (volume.m_height- 1), 
		   gZrate = (gradient.m_depth - 1) / (volume.m_depth - 1);

	for(unsigned int z = 0; z < volume.m_depth ; ++z )
	for(unsigned int y = 0; y < volume.m_height; ++y )
	for(unsigned int x = 0; x < volume.m_width ; ++x )
	{
		int gx = (int)( x * gXrate), 
			gy = (int)( y * gYrate), 
			gz = (int)( z * gZrate);
		int Gidx = 4 * ( gx + gy * gradient.m_width  +  gz * gradient.m_width * gradient.m_height) ;
		int Vidx = 4 * ( x  +  y * volume  .m_width  +   z * volume  .m_width  *  volume.m_height) ;
		
		//double l = calcLuminance( &vol.m_RGBA[idx] );//[0, 1]
		double l = max3( volume.m_RGBA[Vidx+0], volume.m_RGBA[Vidx+1], volume.m_RGBA[Vidx+2] ) / 255.0; //[0, 1]
		double g = gradient.m_RGBA[Gidx + 3] / 255.0;     //[0, 1]

		volHist [ (int)( l * ( binSize - 1.0)) ] += 1.0; //histgram of intensity
		gradHist[ (int)( g * ( binSize - 1.0)) ] += 1.0; //histram  of gradient
		MultHist[ (int)( l * ( binSize - 1.0)) ]
		        [ (int)( g * ( binSize - 1.0)) ] += 1.0;
	}

	//log and normalize
	for( int i = 0; i < binSize; ++i){
		volHist [i] = log( volHist [i] );
		gradHist[i] = log( gradHist[i] );
		for( int j = 0; j < binSize; ++j ) MultHist[i][j] = log( MultHist[i][j]);
	}
	
	double maxForV = 0, maxForG = 0, maxForMult = 0;
	for( int i = 1; i < binSize-1; ++i)
	{
		maxForV = max( maxForV, volHist [i] );
		maxForG = max( maxForG, gradHist[i] );
		for( int j = 1; j < binSize-1; ++j ) maxForMult = max( maxForMult, MultHist[i][j] );
	}
	for( int i = 0; i < binSize; ++i)
	{
		volHist [i] /= maxForV;
		gradHist[i] /= maxForG;
		for( int j = 0; j < binSize; ++j ) MultHist[i][j] /= maxForMult;
	}
}



void TI_calcLogScaleHistgramesOfRGBA(const TOGL3DImage4 &volume, vector<double > &hist_R, vector<double > &hist_G, 
																vector<double > &hist_B, vector<double > &hist_A)
{
	const int binSize = 256;
	hist_R.clear();  hist_R.resize( binSize, 0);
	hist_G.clear();  hist_G.resize( binSize, 0);
	hist_B.clear();  hist_B.resize( binSize, 0);
	hist_A.clear();  hist_A.resize( binSize, 0);

	for(unsigned int z = 0; z < volume.m_depth ; ++z )
	for(unsigned int y = 0; y < volume.m_height; ++y )
	for(unsigned int x = 0; x < volume.m_width ; ++x )
	{
		int Vidx = 4 * ( x  +  y * volume.m_width  +   z * volume.m_width*volume.m_height) ;
		
		//double l = calcLuminance( &vol.m_RGBA[idx] );//[0, 1]
		double r = volume.m_RGBA[Vidx + 0] / 255.0;
		double g = volume.m_RGBA[Vidx + 1] / 255.0;
		double b = volume.m_RGBA[Vidx + 2] / 255.0;
		double a = volume.m_RGBA[Vidx + 3] / 255.0;
		hist_R[ (int)( r * ( binSize - 1.0)) ] += 1.0;
		hist_G[ (int)( g * ( binSize - 1.0)) ] += 1.0;
		hist_B[ (int)( b * ( binSize - 1.0)) ] += 1.0;
		hist_A[ (int)( a * ( binSize - 1.0)) ] += 1.0;
	}

	//log and normalize//
	for( int i = 0; i < binSize; ++i)
	{
		hist_R[i] = (hist_R[i]==0) ? 0 : log( hist_R[i] );
		hist_G[i] = (hist_G[i]==0) ? 0 : log( hist_G[i] );
		hist_B[i] = (hist_B[i]==0) ? 0 : log( hist_B[i] );
		hist_A[i] = (hist_A[i]==0) ? 0 : log( hist_A[i] );
	}
	
	double maxR = 0, maxG = 0, maxB = 0, maxA = 0;
	for( int i = 1; i < binSize-1; ++i){
		maxR = max( maxR, hist_R[i] );
		maxG = max( maxG, hist_G[i] );
		maxB = max( maxB, hist_B[i] );
		maxA = max( maxA, hist_A[i] );
	}
	for( int i = 0; i < binSize; ++i)
	{
		hist_R[i]  /= maxR;
		hist_G[i]  /= maxG;
		hist_B[i]  /= maxB;
		hist_A[i]  /= maxA;
	}
}




/*
inline HBITMAP LoadBitmapFile(CString szFileName)
{
	HANDLE	hFile;
	HDC	hDC;
	HBITMAP	hBmp;
	DWORD	dwFileSize;
	DWORD	dwHeaderSize;
	DWORD	dwScanDataSize;
	DWORD	dw;
	BITMAPFILEHEADER	bmpFileHeader;
	BYTE	*pHeaderBuffer;
	BYTE	*pScanDataBuffer;
	BITMAPINFOHEADER	*pBmpInfoHdr;

	hFile = CreateFile( szFileName, GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);

	dwFileSize = GetFileSize(hFile, NULL);

	ReadFile(hFile, &bmpFileHeader, sizeof(BITMAPFILEHEADER), &dw, NULL);
	dwHeaderSize = bmpFileHeader.bfOffBits - sizeof(BITMAPFILEHEADER);
	dwScanDataSize = dwFileSize - bmpFileHeader.bfOffBits;

	pHeaderBuffer = new BYTE[dwHeaderSize];
	pScanDataBuffer = new BYTE[dwScanDataSize];

	ReadFile(hFile, pHeaderBuffer, dwHeaderSize, &dw, NULL);
	ReadFile(hFile, pScanDataBuffer, dwScanDataSize, &dw, NULL);

	CloseHandle(hFile);

	hDC = GetDC(NULL);
	pBmpInfoHdr = (BITMAPINFOHEADER*)pHeaderBuffer;

	hBmp = CreateDIBitmap(hDC, pBmpInfoHdr, CBM_INIT, pScanDataBuffer, (BITMAPINFO*)pBmpInfoHdr, DIB_RGB_COLORS);

	ReleaseDC(NULL, hDC);
	delete[] pScanDataBuffer;
	delete[] pHeaderBuffer;
	return hBmp;
}

bool TOGL2DImage::loadFromBMP( const char* filename, TOGL* ogl)
{
	FILE* fp = fopen(filename,"r") ;
	if( !fp ){ fprintf( stderr, "cant open file\n" ); return false; }
	fclose(fp) ;

	fprintf( stderr, "a");
	HBITMAP  hBmp = LoadBitmapFile( (CString) filename );
	fprintf( stderr, "b");

	BITMAP    bmp; 
	GetObject( hBmp, sizeof( BITMAP ), &bmp);
	CBitmap *cbmp = CBitmap::FromHandle( hBmp );

	if( bmp.bmBitsPixel != 32 ){ return false;}
	allocateImage( (unsigned int) bmp.bmWidth, (unsigned int) bmp.bmHeight, ogl);

	BYTE *pBuffer     = new BYTE[ bmp.bmWidthBytes * bmp.bmHeight + 1];
	cbmp->GetBitmapBits( bmp.bmWidthBytes * bmp.bmHeight , pBuffer);
		
	for(int j = 0; j < (int) m_height; ++j)
	for(int i = 0; i < (int) m_width * 4; i += 4)
	{
		int idx = j * m_width * 4 + i; 
		m_RGBA[ idx + 0 ] = pBuffer[bmp.bmWidthBytes * j + i + 2];//r
		m_RGBA[ idx + 1 ] = pBuffer[bmp.bmWidthBytes * j + i + 1];//g
		m_RGBA[ idx + 2 ] = pBuffer[bmp.bmWidthBytes * j + i + 0];//b
		m_RGBA[ idx + 3 ] = 0;
	}
	DeleteObject( hBmp );
	delete[] pBuffer;
	return true;
}
*/

/*
static bool t_loadImg1(const char *fname, int &W, int &H, byte* rgba)
{
	CImage pic;
	if( !SUCCEEDED( pic.Load( _T(fname) ) ) ) return false;//bmp jpeg png tiff画像の読み込み(bmp jpegのみ動作確認済み)
	
	W    = pic.GetWidth ();
	H    = pic.GetHeight();
	rgba = new byte[W*H*4];

	for( int y=0, i=0; y<H; ++y     )
	for( int x=0     ; x<W; ++x, ++i)
	{
		COLORREF c = pic.GetPixel( x, y );//画素値get (非常に遅い)
		rgba [ i*4 + 0] = GetRValue(c); //画素値コピー
		rgba [ i*4 + 1] = GetGValue(c); //画素値コピー
		rgba [ i*4 + 2] = GetBValue(c); //画素値コピー
		rgba [ i*4 + 3] = 128;
	}
	return true;
}

static bool t_loadImg2(const char *fname, int &W, int &H, byte* rgba)
{
	CImage pic;
	if( !SUCCEEDED( pic.Load( _T(fname) ) ) ) return false;//bmp jpeg png tiff画像の読み込み(bmp jpegのみ動作確認済み)
	
	//メージのピッチ(行方向のバイト数)を返します。
	//戻り値が負の値の場合、ビットマップは左下隅を起点とする逆方向 (下から上) の DIB です。
	//戻り値が正の値の場合、ビットマップは左上隅を起点とする順方向 (上から下の向き) の DIB です
	int pitch   = pic.GetPitch();
	if( pitch < 0) pitch *= -1;

	W    = pic.GetWidth ();
	H    = pic.GetHeight();
	rgba = new byte[W*H*4];

	int byteNum = pitch / W; //これが1pixelあたりのバイト数


	if( byteNum == 3 ) //24bit color 
	{
		for( int y=0, i=0; y<H; ++y     )
		for( int x=0     ; x<W; ++x, ++i)
		{
			//RGB(24bit)に対するアドレスをbyteにキャストしてるので, R G B順の一番下位のBを指すものになる
			byte *c = (byte*)pic.GetPixelAddress( x, y );
			rgba[ i*4 + 3] = 128;
			rgba[ i*4 + 2] = *c ; ++c;
			rgba[ i*4 + 1] = *c ; ++c;
			rgba[ i*4 + 0] = *c ; 
		}
	}
	else //その他
	{
		for( int y=0, i=0; y<H; ++y     )
		for( int x=0     ; x<W; ++x, ++i)
		{
			COLORREF c = pic.GetPixel( x, y );
			rgba[ i*4 + 3] = 128;
			rgba[ i*4 + 0] = GetRValue(c);
			rgba[ i*4 + 1] = GetGValue(c);
			rgba[ i*4 + 2] = GetBValue(c);
		}
	}	
	return true;
}
*/

/*
BOOL SaveBitmapFile2( const HBITMAP hBmp, LPCTSTR p_pchFileName )
{ 
 //   HBITMAP hBmp = (HBITMAP)GetCurrentObject( p_hDC, OBJ_BITMAP );
	HDC p_hDC = GetDC( NULL );
    BITMAPINFO  stBmpInfo; 
    stBmpInfo.bmiHeader.biSize = sizeof( stBmpInfo.bmiHeader );
    stBmpInfo.bmiHeader.biBitCount = 0;
    GetDIBits( p_hDC, hBmp, 0, 0, NULL, &stBmpInfo, DIB_RGB_COLORS );
 
    ULONG iBmpInfoSize; 
    switch( stBmpInfo.bmiHeader.biBitCount )
    { 
        case 24:
            iBmpInfoSize = sizeof(BITMAPINFOHEADER); 
            break;
        case 16: 
        case 32:
            iBmpInfoSize = sizeof(BITMAPINFOHEADER)+sizeof(DWORD)*3; 
            break;
        default:
            iBmpInfoSize = sizeof(BITMAPINFOHEADER)
                    + sizeof(RGBQUAD)
                    * ( 1 << stBmpInfo.bmiHeader.biBitCount ); 
            break;
    } 
 
    PBITMAPINFO pstBmpInfo;
    if( iBmpInfoSize != sizeof(BITMAPINFOHEADER) ) 
    { 
        pstBmpInfo = (PBITMAPINFO)GlobalAlloc
                    ( GMEM_FIXED | GMEM_ZEROINIT, iBmpInfoSize );
        PBYTE pbtBmpInfoDest
            = (PBYTE)pstBmpInfo; 
        PBYTE pbtBmpInfoSrc
            = (PBYTE)&stBmpInfo; 
        ULONG iSizeTmp
            = sizeof( BITMAPINFOHEADER );
 
        while( iSizeTmp-- ) 
        { 
            *( ( pbtBmpInfoDest )++ ) = *( ( pbtBmpInfoSrc )++ ); 
        } 
    } 

    HANDLE hFile
        = CreateFile
                ( p_pchFileName, GENERIC_WRITE, 0
                , NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_ARCHIVE 
                , NULL );
 
    BITMAPFILEHEADER stBmpFileHder; 
    stBmpFileHder.bfType = 0x4D42; // 'BM' 
    stBmpFileHder.bfSize
        = sizeof(BITMAPFILEHEADER)
        + sizeof(BITMAPINFOHEADER)
        + iBmpInfoSize
        + pstBmpInfo->bmiHeader.biSizeImage; 
    stBmpFileHder.bfReserved1 = 0;
    stBmpFileHder.bfReserved2 = 0; 
    stBmpFileHder.bfOffBits
        = sizeof(BITMAPFILEHEADER) + iBmpInfoSize; 
 
    DWORD dRet;
    WriteFile
        ( hFile, (LPCVOID)&stBmpFileHder
        , sizeof(BITMAPFILEHEADER), &dRet, NULL );
 
    PBYTE pBits
        = (PBYTE)GlobalAlloc
                ( GMEM_FIXED | GMEM_ZEROINIT
                , stBmpInfo.bmiHeader.biSizeImage );

    HBITMAP hBmpOld; 
    HBITMAP hTmpBmp
        = CreateCompatibleBitmap
            ( p_hDC
            , pstBmpInfo->bmiHeader.biWidth
            , pstBmpInfo->bmiHeader.biHeight );
    hBmpOld = (HBITMAP)SelectObject( p_hDC, hTmpBmp ); 
    GetDIBits
        ( p_hDC, hBmp, 0, pstBmpInfo->bmiHeader.biHeight
        , (LPSTR)pBits, pstBmpInfo, DIB_RGB_COLORS );

    WriteFile
        ( hFile, (LPCVOID)pstBmpInfo
        , iBmpInfoSize, &dRet, NULL );
 
    WriteFile( hFile, (LPCVOID)pBits
        , pstBmpInfo->bmiHeader.biSizeImage
        , &dRet, NULL );
 
    SelectObject( p_hDC, hBmpOld ); 
    DeleteObject( hTmpBmp ); 
    CloseHandle( hFile ); 
    GlobalFree( pstBmpInfo ); 
    GlobalFree( pBits ); 
    return TRUE;
} 
 



inline
bool SaveBitmapFile(const HBITMAP hBmp, const char *szFileName, const int iColors) 
{
	HANDLE           hFile;  
	HDC              hDC;  
	DWORD            dw; 
	DWORD            dwHeaderSize; 
	DWORD            dwScanDataSize; 
	BITMAPFILEHEADER bmpFileHeader;
	BITMAPINFOHEADER *pBmpInfoHdr;  
	BYTE			 *pHeaderBuffer; 
	BYTE			 *pScanDataBuffer; 
	BITMAP			 bmp;
	GetObject(hBmp, sizeof(BITMAP), &bmp);

	dwHeaderSize  = sizeof(BITMAPINFOHEADER) + sizeof(RGBQUAD) * iColors;
	pHeaderBuffer = new BYTE[dwHeaderSize]; 
	memset(pHeaderBuffer, 0, dwHeaderSize);

	pBmpInfoHdr = (BITMAPINFOHEADER*)pHeaderBuffer;  
	switch (iColors) 
	{
	case 2:
		pBmpInfoHdr->biBitCount = 1;
		break;
	case 16:
		pBmpInfoHdr->biBitCount = 4;
		break;
	case 256:
		pBmpInfoHdr->biBitCount = 8;
		break;
	case 0:
		pBmpInfoHdr->biBitCount = 24;
		break;
	}

	pBmpInfoHdr->biSize   = sizeof(BITMAPINFOHEADER); 
	pBmpInfoHdr->biWidth  = bmp.bmWidth;  
	pBmpInfoHdr->biHeight = bmp.bmHeight;
	pBmpInfoHdr->biPlanes = 1;

	pBmpInfoHdr->biCompression   = BI_RGB; 
	pBmpInfoHdr->biSizeImage     = 0;
	pBmpInfoHdr->biXPelsPerMeter = 0;
	pBmpInfoHdr->biYPelsPerMeter = 0; 
	pBmpInfoHdr->biClrUsed       = 0;
	pBmpInfoHdr->biClrImportant  = 0;

	hDC = GetDC(NULL);
	GetDIBits(hDC, hBmp, 0, bmp.bmHeight, NULL, (LPBITMAPINFO)pBmpInfoHdr, DIB_RGB_COLORS); 
	dwScanDataSize  = pBmpInfoHdr->biSizeImage; 
	pScanDataBuffer = new BYTE[dwScanDataSize];  
	GetDIBits(hDC, hBmp, 0, bmp.bmHeight, pScanDataBuffer, (LPBITMAPINFO)pBmpInfoHdr, DIB_RGB_COLORS); 
	ReleaseDC(NULL, hDC);

	bmpFileHeader.bfType      = 0x4d42;  //"BM" 
	bmpFileHeader.bfReserved1 = 0;
	bmpFileHeader.bfReserved2 = 0;

	bmpFileHeader.bfOffBits = sizeof(BITMAPFILEHEADER) + dwHeaderSize; 
	bmpFileHeader.bfSize = dwScanDataSize + bmpFileHeader.bfOffBits;

	hFile = CreateFile(szFileName, GENERIC_WRITE, FILE_SHARE_WRITE, NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
	WriteFile(hFile, &bmpFileHeader, sizeof(BITMAPFILEHEADER), &dw, NULL);  WriteFile(hFile, pHeaderBuffer, dwHeaderSize, &dw, NULL);  WriteFile(hFile, pScanDataBuffer, dwScanDataSize, &dw, NULL);  CloseHandle(hFile);

	delete[] pScanDataBuffer;
	delete[] pHeaderBuffer;

	return true;
}


void TOGL2DImage::saveAsBMP( const char* fName )
{
	CDC *pDC = CDC::FromHandle( ::GetDC( NULL ) ) ; 
	//bmp生成
	CBitmap bitmap;
	bitmap.CreateCompatibleBitmap( pDC , m_width, m_height );

	//bitmapに色を書き込む
	CDC bitmapDC;
	bitmapDC.CreateCompatibleDC( pDC );
	bitmapDC.SelectObject( bitmap );

	for(int j = 0; j < (int)m_height; ++j)
	for(int i = 0; i < (int)m_width ; ++i)
	{
		int idx = 4*(i + j*m_width);
		bitmapDC.SetPixel( i, j, RGB( m_RGBA[idx], m_RGBA[idx+1], m_RGBA[idx+2] ) );
	}
	SaveBitmapFile2((HBITMAP) bitmap, fName); 
	//SaveBitmapFile((HBITMAP) bitmap, fName, 24); 
}
*/



