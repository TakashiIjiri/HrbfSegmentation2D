#pragma once

#include "TOGL.h"
#include <vector>
#include "../VIO/SimpleImage.h"

#define INIT_TWODIM_CONST_R 0.1

enum CHANNEL_ID
{
	CHANNEL_ALL,
	CHANNEL_RED,
	CHANNEL_GREEN,
	CHANNEL_BLUE,
	CHANNEL_ALPHA
};


struct TVoxelInfo
{
	int x,y,z,idx;
	void Set(   int _x, int _y, int _z, int _idx){ x = _x   ;  y = _y   ;  z = _z   ;  idx =    _idx;}
	TVoxelInfo( int _x, int _y, int _z, int _idx){ x = _x   ;  y = _y   ;  z = _z   ;  idx =    _idx;}
	TVoxelInfo(           const TVoxelInfo &src){  x = src.x;  y = src.y;  z = src.z;  idx = src.idx;}
	TVoxelInfo& operator=(const TVoxelInfo &src){  x = src.x;  y = src.y;  z = src.z;  idx = src.idx; return *this;}
};

struct TPixelInfo
{
	int x,y,idx;
	void Set(   int _x, int _y, int _idx){ x = _x   ;  y = _y   ;  idx =    _idx;}
	TPixelInfo( int _x, int _y, int _idx){ x = _x   ;  y = _y   ;  idx =    _idx;}
	TPixelInfo(           const TVoxelInfo &src){  x = src.x;  y = src.y; idx = src.idx;}
	TPixelInfo& operator=(const TVoxelInfo &src){  x = src.x;  y = src.y; idx = src.idx; return *this;}
};


class TOGL1DImage
{
public:
	GLubyte* m_RGBA;
	unsigned int m_textureName;
	unsigned int m_width;

	bool m_DoInterpolation;
	~TOGL1DImage();
	TOGL1DImage ();

	TOGL1DImage(const TOGL1DImage & src)
	{
		m_RGBA          = 0;
		m_width			= src.m_width ;
		m_textureName	= src.m_textureName;
		m_DoInterpolation = src.m_DoInterpolation;
		if( src.m_RGBA != 0){
			m_RGBA = new GLubyte[ m_width  * 4 ] ;
			memcpy(m_RGBA, src.m_RGBA, sizeof( GLubyte) * m_width * 4) ;
		}
	}
	//////////////////////////////////
	//make openGl currentÇÇ∑ÇÈ ogl=0ÇÃèÍçáwglMakeCurrentÇµÇ»Ç¢
	void clear ( TOGL* ogl ); 
	void bind  ( TOGL* ogl );
	void unbind( TOGL* ogl );
	void allocateImage( unsigned int width, TOGL* ogl);//oglÇÕclearÇ…óòópÇ≥ÇÍÇÈ
	void allocateHeuColorBar(               TOGL* ogl);
};


class TOGL2DImage
{
public:
	GLubyte* m_RGBA;
	unsigned int m_textureName;
	unsigned int m_width, m_height;

	bool m_DoInterpolation;

	~TOGL2DImage();
	TOGL2DImage ();

	//ÉRÉsÅ[ÉRÉìÉXÉgÉâÉNÉ^ : Ç±ÇÍÇ™Ç»Ç¢Ç∆ vector< TOGL2DImage> Ç™Ç§Ç‹Ç≠ìÆÇ©Ç»Ç¢
	//( allocateÇ™ë´ÇËÇ»Ç≠Ç»Ç¡ÇΩÇÁé©ìÆÇ≈ÉRÉsÅ[Ç∑ÇÈÇ©ÇÁ)
	TOGL2DImage(const TOGL2DImage & src)
	{
		m_RGBA          = 0;
		m_width			= src.m_width ;
		m_height		= src.m_height;
		m_textureName	= src.m_textureName;
		m_DoInterpolation = src.m_DoInterpolation;
		
		if( src.m_RGBA != 0)
		{
			m_RGBA = new GLubyte[ m_width * m_height * 4 ] ;
			memcpy(m_RGBA, src.m_RGBA, sizeof( GLubyte) * m_width * m_height * 4) ;
		}
	}

	//////////////////////////////////
	//make openGl currentÇÇ∑ÇÈ ogl=0ÇÃèÍçáwglMakeCurrentÇµÇ»Ç¢
	void clear (TOGL* ogl); 
	void bind  (TOGL* ogl);
	void unbind(TOGL* ogl);

	//////////////////////////////////
	void allocateImage( unsigned int width, unsigned int height, TOGL* ogl);
	void allocateImage( const TOGL2DImage &src                 , TOGL* ogl);
	bool allocateFromFile( const char *fname, bool &inverted, TOGL* ogl);
	bool saveAsFile      ( const char *fname, int flg_BmpJpgPngTiff = 0 );

	void gaussianFilter33();

	inline void flipImageInY()
	{
		byte *tmp = new byte[ m_width * m_height * 4 ];
		for( unsigned int y = 0; y < m_height; ++y)
		for( unsigned int x = 0; x < m_width ; ++x)
		{
			int iold = 4*(       y        * m_width + x) ;
			int inew = 4*( (m_height-1-y) * m_width + x) ;
			memcpy( &tmp[ inew ], &m_RGBA[ iold ], sizeof( byte ) * 4 ); 
		}
		memcpy( m_RGBA, tmp, sizeof( byte ) * m_width* m_height * 4 );
		delete[] tmp;
	}

	inline bool isAllocated(){ return m_RGBA != 0; }

	void setAllAlphaCh( byte alpha);

	inline void getGrad3( int x, int y, TVector2 &gra)const{
		int idx = 4*(x + y*m_width);

		double pix[3][3] = {{0,0,0}, {0,0,0}, {0,0,0}};
		for( int yy = 0; yy<3; ++yy) if( 1<= y + yy && y+yy <= (int)m_height )
		for( int xx = 0; xx<3; ++xx) if( 1<= x + xx && x+xx <= (int)m_width  ){
			int i = idx + 4*(xx-1 + (yy-1) * m_width);
			pix[yy][xx]	= (m_RGBA[ i] + m_RGBA[ i + 1] + m_RGBA[ i + 2] ) * 0.33333;
		}

		gra[0] = gra[1] = 0;
		gra[0] +=        pix[0][2] - pix[0][0];   gra[1] +=        pix[2][0] - pix[0][0];
		gra[0] += 2.0 * (pix[1][2] - pix[1][0]);  gra[1] += 2.0 * (pix[2][1] - pix[0][1]);
		gra[0] +=        pix[2][2] - pix[2][0];   gra[1] +=        pix[2][2] - pix[0][2];
		gra *= 1/8.0;
	}
	inline void getGrad5( int x, int y, TVector2 &gra){
		int idx = 4*(x + y*m_width);

		double pix[5][5] = { {0,0,0,0,0}, {0,0,0,0,0}, {0,0,0,0,0} , {0,0,0,0,0} , {0,0,0,0,0} };
		for( int yy = 0; yy<5; ++yy) if( 2<= y + yy && y+yy < (int)m_height+2 )
		for( int xx = 0; xx<5; ++xx) if( 2<= x + xx && x+xx < (int)m_width +2 ){
			int i = idx + 4*(xx-2 + (yy-2) * m_width);
			pix[yy][xx]	= (m_RGBA[ i] + m_RGBA[ i + 1] + m_RGBA[ i + 2] ) * 0.33333;
		}

		gra[0] = gra[1] = 0;
		gra[0] += 1.0 * (pix[0][4] + pix[4][4] - pix[0][0]-pix[4][0]);
		gra[0] += 4.0 * (pix[1][4] + pix[3][4] - pix[1][0]-pix[3][0]);
		gra[0] += 6.0 * (pix[2][4]             - pix[2][0]          );
		gra[0] += 2.0 * (pix[0][3] + pix[4][3] - pix[0][1]-pix[4][1]);
		gra[0] += 8.0 * (pix[1][3] + pix[3][3] - pix[1][1]-pix[3][1]);
		gra[0] +=12.0 * (pix[2][3]             - pix[2][1]          );

		gra[1] += 1.0 * (pix[4][0] + pix[4][4] - pix[0][0]-pix[0][4]);
		gra[1] += 4.0 * (pix[4][1] + pix[4][3] - pix[0][1]-pix[0][3]);
		gra[1] += 6.0 * (pix[4][2]             - pix[0][2]          );
		gra[1] += 2.0 * (pix[3][0] + pix[3][4] - pix[1][0]-pix[1][4]);
		gra[1] += 8.0 * (pix[3][1] + pix[3][3] - pix[1][1]-pix[1][3]);
		gra[1] +=12.0 * (pix[3][2]             - pix[1][2]          );

		gra *= 1/96.0;//Ç±Ç±Ç≈48*2=96Ç≈äÑÇ¡ÇƒÇ¢ÇÈÇÃÇΩÇﬂÅAå„Ç≈/2Ç∑ÇÈïKóvÇ»Ç¢
	}

	inline int getPixelIdx ( int x, int y) const{ return  x + y * m_width;}
	inline int getPixelIdx4( int x, int y) const{ return 4 * getPixelIdx( x,y);}

	inline bool getIdxByPos( double xPos, double yPos, double xSize, double ySize, int &xIdx, int &yIdx, int &idx )const
	{
		if( xPos < 0 || xSize < xPos ) return false;
		if( yPos < 0 || ySize < yPos ) return false;

		xIdx = (int) (xPos / xSize * (m_width -1) );
		yIdx = (int) (yPos / ySize * (m_height-1) );
		idx = getPixelIdx( xIdx,  yIdx);
		return true;
	}
	
	inline void setPix( int idx     , byte* rgba){ memcpy( &m_RGBA[idx], rgba, sizeof( byte ) * 4); }
	inline void setPix( int x, int y, byte* rgba){ int idx = 4*(x + y*m_width);
		                                           memcpy( &m_RGBA[idx], rgba, sizeof( byte ) * 4 );}

	inline void setPix( int x, int y, byte r, byte g, byte b, byte a){
		int idx = 4*(x + y*m_width);
		m_RGBA[idx+0] = r;
		m_RGBA[idx+1] = g;
		m_RGBA[idx+2] = b;
		m_RGBA[idx+3] = a;
	}
	inline void setPix( int id4, byte r, byte g, byte b, byte a){
		m_RGBA[id4+0] = r;
		m_RGBA[id4+1] = g;
		m_RGBA[id4+2] = b;
		m_RGBA[id4+3] = a;
	}

	inline double getIntens( int id4 )const{ return (m_RGBA[ id4+0] + m_RGBA[ id4+1] + m_RGBA[ id4+2] ) * 0.333333333;}

};




class TOGL3DImage1;

class TOGL3DImage4
{
public:
	GLubyte*     m_RGBA       ;
	unsigned int m_textureName;
	unsigned int m_width, m_height, m_depth;// x, y, z

	double       m_maxGradientValue; //only for gradient volume image
	bool         m_DoInterpolation ; //bindëOÇ…ÉZÉbÉgÇ∑ÇÈ

public:
	~TOGL3DImage4();
	TOGL3DImage4 ();
	TOGL3DImage4(const TOGL3DImage4 & src)
	{
		m_RGBA             = 0;
		m_width			   = src.m_width ;
		m_height		   = src.m_height;
		m_depth			   = src.m_depth;
		m_textureName	   = src.m_textureName;
		m_maxGradientValue = src.m_maxGradientValue;
		m_DoInterpolation  = src.m_DoInterpolation;
		
		if( src.m_RGBA != 0){
			m_RGBA = new GLubyte[ m_width * m_height * m_depth * 4 ] ;
			memcpy(m_RGBA, src.m_RGBA, sizeof( GLubyte) * m_width * m_height * m_depth * 4) ;
		}
	}

	//////////////////////////////////
	void clear (const TOGL* ogl); 
	void bind  (TOGL* ogl);		
	void unbind(TOGL* ogl);

	void allocateImage        ( unsigned int width, unsigned int  height, unsigned int depth, TOGL* ogl);
	void allocateImage        ( const TOGL3DImage4 &src                                     , TOGL* ogl);
	void allocateGradientImage( const TOGL3DImage1 &trgtVol, double gradMagCoef, TOGL* ogl);

	void allocateSphereVolume( int sizeR , TOGL *ogl );

	inline bool isAllocated(){ return m_RGBA != 0; }

	void setVolumeColor(      const Image3D              &img   , CHANNEL_ID id);
	void setVolumeColor(      const vector<TOGL2DImage*> &images, CHANNEL_ID id);
	void setVolumeColor_pack( const Image3D              &img   , CHANNEL_ID id);
	void setVolumeColor_pack( const vector<TOGL2DImage*> &images, CHANNEL_ID id);

	void convertToImage3D( Image3D& trgtImg, const TVector3 &cubeSize, CHANNEL_ID id) const;

	void simpleSmoothing(int n);
	void simpleSmoothing();
	
	void halfSize(                                                              TOGL *ogl );
	void resize( int  newWidth, int newHeight, int newDepth                   , TOGL *ogl );
	void clipImage( int minX, int minY, int minZ, int maxX, int maxY, int maxZ, TOGL *ogl );

	void normalizeGradientImage();
	void gainContrast( GLubyte bottom, GLubyte top);
	void setIntensityToAlphaChannel();


	inline int getVoxelIdx( int x, int y, int z ) const{
		return  x + y * m_width + z * m_width * m_height;
	}

	inline int getVoxelIdx4( int x, int y, int z ) const{
		return 4 * getVoxelIdx( x,y,z);
	}

	inline bool getIdxByPos( const TVector3 &p, const TVector3 &volSize, int &idx ) const
	{
		if( p.data[0] < 0 || volSize.data[0] < p.data[0] ) return false;
		if( p.data[1] < 0 || volSize.data[1] < p.data[1] ) return false;
		if( p.data[2] < 0 || volSize.data[2] < p.data[2] ) return false;

		unsigned int xIdx = (int) (p.data[0] / volSize.data[0] * (m_width -1) );
		unsigned int yIdx = (int) (p.data[1] / volSize.data[1] * (m_height-1) );
		unsigned int zIdx = (int) (p.data[2] / volSize.data[2] * (m_depth -1) );

		idx = getVoxelIdx( xIdx,  yIdx,  zIdx);
		return true;
	}

	inline bool getIdxByPos( const TVector3 &p, const TVector3 &volSize, int &xIdx, int &yIdx, int &zIdx, int &idx )const
	{
		if( p.data[0] < 0 || volSize.data[0] < p.data[0] ) return false;
		if( p.data[1] < 0 || volSize.data[1] < p.data[1] ) return false;
		if( p.data[2] < 0 || volSize.data[2] < p.data[2] ) return false;

		xIdx = (int) (p.data[0] / volSize.data[0] * (m_width -1) );
		yIdx = (int) (p.data[1] / volSize.data[1] * (m_height-1) );
		zIdx = (int) (p.data[2] / volSize.data[2] * (m_depth -1) );
		idx = getVoxelIdx( xIdx,  yIdx,  zIdx);
		return true;
	}
};




class TOGL3DImage1
{
public:
	GLubyte*     m_img        ;
	unsigned int m_textureName;
	unsigned int m_width, m_height, m_depth;// x, y, z
	bool         m_DoInterpolation;

	~TOGL3DImage1();
	TOGL3DImage1();
	TOGL3DImage1(const TOGL3DImage1 &src){
		m_img = 0;
		m_width           = src.m_width ;
		m_height          = src.m_height;
		m_depth           = src.m_depth ;
		m_textureName     = src.m_textureName;
		m_DoInterpolation = src.m_DoInterpolation;
		if( src.m_img != 0){
			m_img = new GLubyte[ m_width * m_height * m_depth ];
			memcpy( m_img, src.m_img , sizeof( GLubyte ) * m_width * m_height * m_depth );
		} 
	}

	void clear (const TOGL* ogl); 
	void bind  (const TOGL* ogl);
	void unbind(const TOGL* ogl);

	//allocation////////////////////////////////////////////////////////////////////////////////////////
	void allocateImage( unsigned int width, unsigned int  height, unsigned int depth , const TOGL* ogl);
	void allocateImage( unsigned int width, unsigned int  height, unsigned int depth , const float* img, const double minV, const double maxV, const TOGL* ogl);
	void allocateImage                (const TOGL3DImage1 &srcImg                    , const TOGL* ogl);
	void allocateImageFromAlphaChannel(const TOGL3DImage4 &vol                       , const TOGL* ogl);
	void allocateSphereVolume         ( int sizeR                                    , const TOGL* ogl);

	void allocateGradMagImage( const TOGL3DImage1 &r,                                               const TOGL* ogl, double scale );
	void allocateGradMagImage( const TOGL3DImage1 &r, const TOGL3DImage1 &g,                        const TOGL* ogl, double scale );
	void allocateGradMagImage( const TOGL3DImage1 &r, const TOGL3DImage1 &g, const TOGL3DImage1 &b, const TOGL* ogl, double scale );
	void calcGrad33          ( const int x, const int y, const int z, TVector3 &grad ) const;
	void calcGrad55          ( const int x, const int y, const int z, TVector3 &grad ) const;

	//standard manipulation//
	void simpleSmoothing(int n);
	void simpleSmoothing();
	void gaussianSmoothing33(int n=1);
	void gaussianSmoothing55(int n=1);
	void gaussianSmoothing55(int n, int sX,int sY,int sZ, int eX,int eY,int eZ);


	void gainContrast( GLubyte bottom, GLubyte top);
	void resize(int newWidth, int newHeight, int newDepth, TOGL *ogl);

	void setImage        (const Image3D      &img);
	void setImage_fitSize(const TOGL3DImage1 &img);
	void setImage        (const vector<TOGL2DImage > &images, int channel/*0:r 1:g 2:b*/);
	void setImageIntens  (const Image3D &Rimg, const Image3D &Gimg, const Image3D &Bimg);
	void setImageIntens  (const vector<TOGL2DImage > &images);

	inline int getVolIdxSize()const{return (int)(m_width*m_height*m_depth);}
	inline void copy(const byte         *img){ if( img != 0 ) memcpy(m_img, img, sizeof(GLbyte)*getVolIdxSize());}
	inline void copy(const TOGL3DImage1 &img){ if( img.isAllocated() ) copy(img.m_img );};
	inline bool isAllocated()const{ return m_img != 0; }
	inline void setImgZero() { memset( m_img, 0, sizeof(GLubyte)*getVolIdxSize());}

	int labeling26ConnectedRegion(byte trgtI, int *label/*should be already allocated*/) const ;

	inline int getVoxelIdx( int x, int y, int z ) const{
		return  x + y * m_width + z * m_width * m_height;
	}
	inline void getXYZIdx( const int idx, int &x, int &y, int &z ) const{
		z = idx / (m_width*m_height);
		y = ( idx - z * m_width*m_height ) / m_width;
		x = idx - z * m_width*m_height - y * m_width;
	}

	inline bool getIdxByPos( const TVector3 &p, const TVector3 &volSize, int &xIdx, int &yIdx, int &zIdx, int &idx )const
	{
		if( p.data[0] < 0 || volSize.data[0] < p.data[0] ) return false;
		if( p.data[1] < 0 || volSize.data[1] < p.data[1] ) return false;
		if( p.data[2] < 0 || volSize.data[2] < p.data[2] ) return false;
		double pichX = volSize.data[0] / m_width ;
		double pichY = volSize.data[1] / m_height;
		double pichZ = volSize.data[2] / m_depth ;
		xIdx = (int)( (p.data[0] ) / pichX ); if(xIdx<0)xIdx=0; if(xIdx>(int)m_width -1) xIdx=m_width -1;
		yIdx = (int)( (p.data[1] ) / pichY ); if(yIdx<0)yIdx=0; if(yIdx>(int)m_height-1) yIdx=m_height-1;
		zIdx = (int)( (p.data[2] ) / pichZ ); if(zIdx<0)zIdx=0; if(zIdx>(int)m_depth -1) zIdx=m_depth -1;
		idx = getVoxelIdx( xIdx,  yIdx,  zIdx);
		return true;
	}
	inline bool getIdxByPos( const TVector3 &p, const TVector3 &volSize, int &idx ) const{
		int x,y,z;
		return getIdxByPos(p, volSize, x,y,z,idx);
	}

	inline byte getMaxVal(){
		byte maxV=0;
		for( int i=0, size=getVolIdxSize(); i<size; ++i) maxV = max( maxV, m_img[i]);
		return maxV;
	}

	inline GLubyte& operator[](int id)       { return m_img[id];}
	inline GLubyte  operator[](int id) const { return m_img[id];}
	
	inline void flipImageInY()
	{
		if( !isAllocated() ) return;
		const int WH = m_width*m_height, WHD = m_width * m_height * m_depth;
		byte *tmp = new byte[ m_width * m_height * m_depth ];
		for( unsigned int z = 0; z < m_depth ; ++z)
		for( unsigned int y = 0; y < m_height; ++y)
		for( unsigned int x = 0; x < m_width ; ++x)
			tmp[ z*WH + (m_height-1-y) * m_width + x  ] = m_img[ z * WH + y* m_width + x  ];

		memcpy( m_img, tmp, sizeof( byte ) * WHD );
		delete[] tmp;
	}
	double triLinearSampling( const TVector3 &cSize, const TVector3 &pos ) const;
};


void TI_calcLogScaleHistgramesOfIntensity( const TOGL3DImage4 &volume, 
										   const TOGL3DImage4 &gradient, 
										   vector<       double > &volHist ,
										   vector<       double > &gradHist,
										   vector<vector<double>> &MultHist ); 

void TI_calcLogScaleHistgramesOfRGBA(const TOGL3DImage4 &gradient, 
									   vector<double > &hist_R,
									   vector<double > &hist_G,
									   vector<double > &hsit_B,
									   vector<double > &hist_a);

inline void smoothing( GLubyte *trgt, GLubyte *src0, 
								      GLubyte *src1, GLubyte *src2, GLubyte *src3,
								      GLubyte *src4, GLubyte *src5, GLubyte *src6)
{
	//x
	double d = 0.1 * (4 * src0[0] + src1[0] + src2[0] + src3[0] + src4[0] + src5[0] + src6[0]); trgt[0]  = ( d > 255.0 )? 255 : (GLubyte) d;
	d        = 0.1 * (4 * src0[1] + src1[1] + src2[1] + src3[1] + src4[1] + src5[1] + src6[1]); trgt[1]  = ( d > 255.0 )? 255 : (GLubyte) d;
	d        = 0.1 * (4 * src0[2] + src1[2] + src2[2] + src3[2] + src4[2] + src5[2] + src6[2]); trgt[2]  = ( d > 255.0 )? 255 : (GLubyte) d;
}


/*///////////////////////////////////////////////////////////////////////////////////////////
imgInOutÇ…imageÇÃëOåiîwåièÓïÒÇ™äiî[Ç≥ÇÍÇƒÇ¢ÇÈ 0 - out   255 - in
Ç±ÇÍÇErodeé˚èkÇ≥ÇπÇÈèàóùÇçsÇ§
W,HÇÕÇªÇÍÇºÇÍwidth height
/////////////////////////////////////////////////////////////////////////////////////////////*/
inline void t_imgMorpho_Erode( byte *imgInOut, const int W, const int H)
{
	const int WH = W*H;

#pragma omp parallel for 
	for( int y = 0; y < H; ++y){
		for( int x = 0; x < W; ++x) 
		{
			int idx = x + y * W ;
			if( imgInOut[idx] != 255 ) continue;

			//é©ï™Ç™ã´äEÇ≈Ç†ÇÍÇŒerode flag (2ÇÇΩÇƒÇƒÇ®Ç≠)
			if( x == 0 || y == 0 || x == W-1 || y == H-1 ) imgInOut[idx] = 2;
			else if( !imgInOut[ idx - 1] || !imgInOut[ idx - W] ||
				     !imgInOut[ idx + 1] || !imgInOut[ idx + W] ) imgInOut[idx] = 2;
		}
	}
#pragma omp parallel for 
	for( int i = 0; i < WH; ++i) if( imgInOut[i] == 2) imgInOut[i] = 0;
}

inline void t_imgMorpho_Dilate(byte *imgInOut, const int W, const int H)
{
	const int WH = W*H;

#pragma omp parallel for 
	for( int y = 0; y < H; ++y){
		for( int x = 0; x < W; ++x) 
		{
			int idx = x + y * W;
			if( imgInOut[idx] != 0) continue;

			//é©ï™Ç™ã´äEÇ≈Ç†ÇÍÇŒerode flag (2ÇÇΩÇƒÇƒÇ®Ç≠)
			if(  ( 0 <= x-1 && imgInOut[ idx - 1]==255)  ||  ( 0 <= y-1 && imgInOut[ idx - W]==255)  ||
				 ( W >  x+1 && imgInOut[ idx + 1]==255)  ||  ( H >  y+1 && imgInOut[ idx + W]==255)  ) imgInOut[idx] = 2;
		}
	}
#pragma omp parallel for 
	for( int i = 0; i < WH; ++i) if( imgInOut[i] == 2) imgInOut[i] = 255;
}




/*///////////////////////////////////////////////////////////////////////////////////////////
voLInOutÇ…volumeÇÃëOåiîwåièÓïÒÇ™äiî[Ç≥ÇÍÇƒÇ¢ÇÈ 0 - out   255 - in
Ç±ÇÍÇErodeé˚èkÇ≥ÇπÇÈèàóùÇçsÇ§
W,H,DÇÕÇªÇÍÇºÇÍwidth height depth
/////////////////////////////////////////////////////////////////////////////////////////////*/
inline void t_volMorpho_Erode(byte *volInOut, const int W, const int H, const int D, const int minIdx[3], const int maxIdx[3])
{
	const int WH = W*H, WHD = WH*D;

#pragma omp parallel for 
	for( int z = minIdx[2]; z <  maxIdx[2]; ++z)
	{
		for( int y = minIdx[1]; y <  maxIdx[1]; ++y)
		for( int x = minIdx[0]; x <  maxIdx[0]; ++x) 
		{
			int idx = x + y * W + z*WH;
			if( volInOut[idx] != 255 ) continue;

			//é©ï™Ç™ã´äEÇ≈Ç†ÇÍÇŒerode flag (2ÇÇΩÇƒÇƒÇ®Ç≠)
			if( x == 0 || y == 0 || z == 0 || x == W-1 || y == H-1 || z == D-1 ) volInOut[idx] = 2;
			else if( !volInOut[ idx - 1] || !volInOut[ idx - W] || !volInOut[ idx - WH] ||
				     !volInOut[ idx + 1] || !volInOut[ idx + W] || !volInOut[ idx + WH] ) volInOut[idx] = 2;
		}
	}
#pragma omp parallel for 
	for( int i = 0; i < WHD; ++i) if( volInOut[i] == 2) volInOut[i] = 0;

}

inline void t_volMorpho_Erode(byte *voLInOut, const int W, const int H, const int D)
{
	int minIdx[3] = {0,0,0};
	int maxIdx[3] = {W,H,D};
	t_volMorpho_Erode(voLInOut, W, H, D, minIdx, maxIdx);
}


inline void t_volMorpho_Dilate(byte *voLInOut, const int W, const int H, const int D, const int minIdx[3], const int maxIdx[3])
{
	const int WH = W*H, WHD = WH*D;

#pragma omp parallel for 
	for( int z = minIdx[2]; z <  maxIdx[2]; ++z){
	for( int y = minIdx[1]; y <  maxIdx[1]; ++y)
	for( int x = minIdx[0]; x <  maxIdx[0]; ++x) 
	{
		int idx = x + y * W + z*WH;
		if( voLInOut[idx] != 0) continue;

		//é©ï™Ç™ã´äEÇ≈Ç†ÇÍÇŒerode flag (2ÇÇΩÇƒÇƒÇ®Ç≠)
		if(  ( 0 <= x-1 && voLInOut[ idx - 1]==255)  ||  ( 0 <= y-1 && voLInOut[ idx - W]==255)  ||  ( 0 <= z-1 && voLInOut[ idx - WH]==255) ||
			 ( W >  x+1 && voLInOut[ idx + 1]==255)  ||  ( H >  y+1 && voLInOut[ idx + W]==255)  ||  ( D >  z+1 && voLInOut[ idx + WH]==255) ) voLInOut[idx] = 2;
	}}
#pragma omp parallel for 
	for( int i = 0; i < WHD; ++i) if( voLInOut[i] == 2) voLInOut[i] = 255;
}

inline void t_volMorpho_Dilate(byte *voLInOut, const int W, const int H, const int D)
{
	int minIdx[3] = {0,0,0};
	int maxIdx[3] = {W,H,D};
	t_volMorpho_Dilate(voLInOut, W, H, D, minIdx, maxIdx);
}


