#ifndef __ILMATH_H_INCLUDED__
#define __ILMATH_H_INCLUDED__

/*tmath.h//////////////////////////////////////////////////////////////////////
written by takashi ijiri @ Riken
2011/10/7
defined class
TVector3  TMatrix16  TMatrix9  TMatrix4
/////////////////////////////////////////////////////////////////////////////*/


#pragma warning (disable:4786)
#pragma warning (disable:4996)

#include <memory.h>
#include <cfloat>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <vector>
#include <list>
using namespace std ;

//#define for if(0);else for

#ifndef VSN
#define VSN 0.00000000001
#endif

#ifndef byte
typedef unsigned char byte;
#endif

#ifndef MIN
#define MIN(a,b)	((a)<(b)?(a):(b))
#define MAX(a,b)	((a)>(b)?(a):(b))
#endif //MIN

#ifndef MIN3
#define MIN3(  a,b,c)	((a)<(b)?((a)<(c)?(a):(c)):((b)<(c)?(b):(c)))
#define MAX3(  a,b,c)	((a)>(b)?((a)>(c)?(a):(c)):((b)>(c)?(b):(c)))
#define MIN3ID(a,b,c)	((a)<(b)?((a)<(c)?(0):(2)):((b)<(c)?(1):(2)))
#define MAX3ID(a,b,c)	((a)>(b)?((a)>(c)?(0):(2)):((b)>(c)?(1):(2)))
template <class T> T   max3  (T a,T b,T c){return MAX3(  a,b,c);}
template <class T> int max3id(T a,T b,T c){return MAX3ID(a,b,c);}
template <class T> T   min3  (T a,T b,T c){return MIN3(  a,b,c);}
template <class T> int min3id(T a,T b,T c){return MIN3ID(a,b,c);}
#endif

#ifndef MIN4
#define MIN4(a,b,c,d) MIN(MIN(a,b),MIN(c,d))
#endif
#ifndef MAX4
#define MAX4(a,b,c,d) MAX(MAX(a,b),MAX(c,d))
#endif

#ifndef CLIPUCHAR
#define CLIPUCHAR(a)	((a)<0?0:((a)>255?255:(unsigned char)(a)))
inline unsigned char ClipUchar(double       a){ return CLIPUCHAR(a); }
inline unsigned char ClipUchar(float        a){ return CLIPUCHAR(a); }
inline unsigned char ClipUchar(int          a){ return CLIPUCHAR(a); }
inline unsigned char ClipUchar(unsigned int a){ return CLIPUCHAR(a); }
#endif



#ifndef M_PI
#define M_PI 3.14159265358979
#endif


#ifndef _WIN32
class CRect {
public:
	int top,bottom,right,left;
};
class CPoint {
public:
	CPoint(int cx=0,int cy=0):x(cx),y(cy){}
	int x,y;
};
#endif //_WIN32


template <class T>  bool SquareMat_CalcInvMatrix( T* src,T* targ,int dim ){
	if( dim <= 0 )	return false ;
	if( dim == 1 ){	if( *src==0 ) return false ; *targ = 1/(*src) ; return true ; }
	if( dim == 2 ){
		T det = src[0]*src[3] - src[1]*src[2] ;
		if( det==0 )	return false ;
		T det_inv = 1/det ;
		targ[0] = src[3]*det_inv ;		targ[1] = -src[1]*det_inv ;
		targ[2] = -src[2]*det_inv ;		targ[3] = src[0]*det_inv ;
		return true ;
	}

	T* adjugate = new T[ dim*dim ] ;	// Adjugate matrix = 余因子行列
	int dim_s = dim-1 ;
	T* submat = new T[ dim_s*dim_s ] ;
	T det = 0 ;
	int i,j,ii,jj ;
	{	// Calc j == 0 ( includes det calculation )
		for( i=0;i<dim;i++ ){
			for( ii=0;ii<dim_s;ii++ ){
				if( ii < i )	memcpy( &submat[ii*dim_s],&src[ii*dim+1],sizeof(T)*dim_s ) ;
				else			memcpy( &submat[ii*dim_s],&src[(ii+1)*dim+1],sizeof(T)*dim_s ) ;
			}
			// cofactor = 余因子
			T cofactor = pow(-1.0,i) * SquareMat_CalcDeterminant( submat,dim_s ) ;
			adjugate[i] = cofactor ;
			det += src[i*dim] * cofactor ;
		}
	}

	if( det == 0 ){
		delete[] adjugate ;
		delete[] submat ;
		return false ;
	}

	for( i=0;i<dim;i++ )for( j=1;j<dim;j++ ){
		for( ii=0;ii<dim_s;ii++ )for( jj=0;jj<dim_s;jj++ ){
			if( ii < i && jj < j)	submat[ii*dim_s + jj] = src[ii*dim + jj] ;
			else if( ii >= i && jj < j)	submat[ii*dim_s + jj] = src[(ii+1)*dim + jj] ;
			else if( ii < i && jj >= j)	submat[ii*dim_s + jj] = src[ii*dim + (jj+1)] ;
			else /*if( ii >= i && jj >= j)*/ submat[ii*dim_s + jj] = src[(ii+1)*dim + (jj+1)] ;
		}
		adjugate[ j*dim + i ] = pow(-1.0,i+j) * SquareMat_CalcDeterminant( submat,dim_s ) ;
	}

	T det_inv = 1/det ;
	for( i=0;i<dim;i++ )for( j=0;j<dim;j++ ){
		targ[i*dim+j] = det_inv*adjugate[i*dim+j] ;
	}

	delete[] adjugate ;
	delete[] submat ;
	return true ;

}
inline static double SquareMat_CalcDeterminant( double* src,int dim ){
	if( dim <= 0 )	return 0 ;
	if( dim == 1 )	return *src ;
	if( dim == 2 )	return src[0]*src[3] - src[1]*src[2] ;
	if( dim == 3 )	return src[0]*src[4]*src[8]+src[3]*src[7]*src[2]+src[1]*src[5]*src[6]
		- src[2]*src[4]*src[6] - src[1]*src[3]*src[8] - src[0]*src[5]*src[7] ;

	int dim_s = dim-1 ;
	double* tmp = new double[ dim_s*dim_s ] ;
	double result = 0 ;
	int i,ii ;
	for( i=0;i<dim;i++ ){
		for( ii=0;ii<dim_s;ii++ ){
			if( ii < i )	memcpy( &tmp[ii*dim_s],&src[ii*dim+1],sizeof(double)*dim_s ) ;
			else			memcpy( &tmp[ii*dim_s],&src[(ii+1)*dim+1],sizeof(double)*dim_s ) ;
		}
		result += src[i*dim] * pow(-1.0,i) * SquareMat_CalcDeterminant( tmp,dim_s ) ;
	}
	delete[] tmp ;

	return result ;
}

inline static float SquareMat_CalcDeterminant( float* src,int dim ){
	if( dim <= 0 )	return 0 ;
	if( dim == 1 )	return *src ;
	if( dim == 2 )	return src[0]*src[3] - src[1]*src[2] ;
	if( dim == 3 )	return src[0]*src[4]*src[8]+src[3]*src[7]*src[2]+src[1]*src[5]*src[6]
		- src[2]*src[4]*src[6] - src[1]*src[3]*src[8] - src[0]*src[5]*src[7] ;

	int dim_s = dim-1 ;
	float* tmp = new float[ dim_s*dim_s ] ;
	float result = 0 ;
	int i,ii ;
	for( i=0;i<dim;i++ ){
		for( ii=0;ii<dim_s;ii++ ){
			if( ii < i )	memcpy( &tmp[ii*dim_s],&src[ii*dim+1],sizeof(float)*dim_s ) ;
			else			memcpy( &tmp[ii*dim_s],&src[(ii+1)*dim+1],sizeof(float)*dim_s ) ;
		}
		result += (float) (src[i*dim] * pow(-1.0,i) * SquareMat_CalcDeterminant( tmp,dim_s )) ;
	}
	delete[] tmp ;

	return result ;
}













class TMatrix4 ;
class TMatrix9 ;
class TMatrix16;

class TVector3 {
public:
	double data[3];
	//Constructors////////////////////////////////////////////////////////////////
	TVector3()                             { data[0]=data[1]=data[2]=0; }
	TVector3(double f1,double f2,double f3){ data[0]=f1;data[1]=f2;data[2]=f3; }
	TVector3(double f[3]){ if( f ){ memcpy(data,f,sizeof(double)*3); } }	

	//TMatrix16& operator=(const TMatrix16 src){memcpy(data,src.data,sizeof(double)*16);return*this;}

	inline void Set(double f1=0.0,double f2=0.0,double f3=0.0){ data[0]=f1;data[1]=f2;data[2]=f3; }
	inline void Set(const double* f){ memcpy(data,f,sizeof(double)*3); }
	inline void Set(const TVector3 &v){ memcpy(data, v.data,sizeof(double)*3);}

	inline double x(){return data[0];}
	inline double y(){return data[1];}
	inline double z(){return data[2];}
	inline   double& operator[](int id)       { return data[id]; }
	inline   double  operator[](int id) const { return data[id]; }
	operator double*(){ return data; }

	operator void  *(){ return this; }
	void Trace() const { fprintf(stderr,"%f %f %f\n",data[0],data[1],data[2]); }
	

	inline double getVectorAngle2D() const{
		double l = Length2D();
		return atan2( data[1] / l, data[0] / l); 
	}

	inline bool Normalize_Self2D(){
		double l = Length2D(); if(l == 0.0) return false;
		data[0] /= l;
		data[1] /= l;
		return true;
	}
	inline bool Normalize_Self(){
		double l = Length(); if(l==0.0f) return false;
		data[0] /= l;
		data[1] /= l;
		data[2] /= l;
		return true;
	}
	inline TVector3 Normalize(){
		double l = Length(); 
		if(l==0.0f) return TVector3(1,0,0);
		else return        TVector3( data[0] / l, data[1] / l, data[2] /l);
	}

	inline double Length_Square()   const { return data[0]*data[0] + data[1]*data[1] + data[2]*data[2]; }
	inline double Length2D_Square() const { return data[0]*data[0] + data[1]*data[1]                  ; }	
	inline double Length()          const { return sqrt( Length_Square  () ); }
	inline double Length2D()        const { return sqrt( Length2D_Square() ); }

	//////operators + & - //////////////////////////////////////////////////////
	inline TVector3& operator=(const TVector3& vec1) {
		memcpy( (void*)data,vec1.data,sizeof(double)*3 ) ;
		return *this ;
	}
	inline TVector3 operator+(const TVector3& vec1) const {
		return TVector3( data[0] + vec1.data[0],
						 data[1] + vec1.data[1],
						 data[2] + vec1.data[2]);
	}
	inline TVector3& operator+=(const TVector3& vec1){
		data[0] += vec1.data[0];
		data[1] += vec1.data[1];
		data[2] += vec1.data[2];
		return *this ;
	}
	inline TVector3 operator-(const TVector3& vec1) const {
		return TVector3( data[0]-vec1.data[0],
						 data[1]-vec1.data[1],
						 data[2]-vec1.data[2]);
	}
	inline TVector3& operator-=(const TVector3& vec1){
		data[0] -= vec1.data[0];
		data[1] -= vec1.data[1];
		data[2] -= vec1.data[2];
		return *this ;
	}
	inline TVector3 operator-() const { return (-1)**this; }
	//multiply//////////////////////////////////////////////////////////////////////////////////////////////
	inline        TVector3  operator* (double f) const{ return TVector3( f*data[0], f*data[1],f*data[2]); }
	inline friend TVector3  operator* (double f, const TVector3& vec1){ return vec1*f; }
	inline        TVector3& operator*=(double f){ data[0] *= f;
												  data[1] *= f;
												  data[2] *= f;
												  return *this ;}
	//divide/////////////////////////////////////////////////////////////////////////////////////////////////
	inline        TVector3  operator/(double f                      ) const{ return TVector3( data[0]/f,data[1]/f,data[2]/f); }
	inline friend TVector3  operator/(double f,const TVector3& vec1){ return vec1 / f; }
	inline        TVector3& operator/=(double f){ data[0]/=f;
												  data[1]/=f;
												  data[2]/=f;
												  return *this; }
	//dot product///////////////////////////////////////////////////////////////////////////////////////////
	inline double operator*(const TVector3& v) const { return data[0]*v.data[0] + data[1]*v.data[1] + data[2]*v.data[2];}
	//cross product/////////////////////////////////////////////////////////////////////////////////////////
	inline TVector3 operator^(const TVector3& v) const {
		return TVector3( data[1]*v.data[2] - data[2]*v.data[1],
						 data[2]*v.data[0] - data[0]*v.data[2],
						 data[0]*v.data[1] - data[1]*v.data[0]);
	}
	inline TVector3& operator^=(const TVector3& v){
		double res[3];
		res[0] = data[1]*v.data[2] - data[2]*v.data[1];
		res[1] = data[2]*v.data[0] - data[0]*v.data[2];
		res[2] = data[0]*v.data[1] - data[1]*v.data[0];
		data[0]=res[0]; data[1]=res[1]; data[2]=res[2];
		return *this ;
	}

	////fast setters////////////////////////////////////////////////////////////////////
	inline void SetSubtract(const TVector3& v1, const TVector3& v2){    data[0] = v1.data[0] - v2.data[0];
																		data[1] = v1.data[1] - v2.data[1];
																		data[2] = v1.data[2] - v2.data[2]; }

	inline void AddSubtract(const TVector3& v1, const TVector3& v2){    data[0] += v1.data[0] - v2.data[0];
																		data[1] += v1.data[1] - v2.data[1];
																		data[2] += v1.data[2] - v2.data[2]; }

	inline void SetAddition(const TVector3& v1, const TVector3& v2){    data[0]  = v1.data[0] + v2.data[0];
																		data[1]  = v1.data[1] + v2.data[1];
																		data[2]  = v1.data[2] + v2.data[2];}

	inline void Set_CoefMultVec(double coef, const TVector3& v){		data[0] = coef * v.data[0];
																		data[1] = coef * v.data[1];
																		data[2] = coef * v.data[2]; }

	inline void Add_CoefMultVec(double coef, const TVector3& v){        data[0] += coef * v.data[0];
																		data[1] += coef * v.data[1];
																		data[2] += coef * v.data[2]; }

	inline void Set_V1_Add_CoefMultV2(const TVector3 &v1, double coef, const TVector3& v2){
		data[0] = v1.data[0] + coef * v2.data[0];
		data[1] = v1.data[1] + coef * v2.data[1];
		data[2] = v1.data[2] + coef * v2.data[2];
	}

	inline void SetCrossProd(const TVector3& v1, const TVector3& v2){
		data[0] = v1.data[1] * v2.data[2] - v1.data[2] * v2.data[1];
		data[1] = v1.data[2] * v2.data[0] - v1.data[0] * v2.data[2];
		data[2] = v1.data[0] * v2.data[1] - v1.data[1] * v2.data[0];
	}
	inline void AddCrossProd(const TVector3& v1, const TVector3& v2){
		data[0] += v1.data[1] * v2.data[2] - v1.data[2] * v2.data[1];
		data[1] += v1.data[2] * v2.data[0] - v1.data[0] * v2.data[2];
		data[2] += v1.data[0] * v2.data[1] - v1.data[1] * v2.data[0];
	}
	inline void AddCrossProdWithCoef(double c, const TVector3& v1, const TVector3& v2){
		data[0] += c* (v1.data[1] * v2.data[2] - v1.data[2] * v2.data[1]);
		data[1] += c* (v1.data[2] * v2.data[0] - v1.data[0] * v2.data[2]);
		data[2] += c* (v1.data[0] * v2.data[1] - v1.data[1] * v2.data[0]);
	}
	
	inline void SetGravityCenter(const TVector3 &v0, const  TVector3 &v1, const TVector3 &v2){
		data[0] = ( v0.data[0] + v1.data[0] + v2.data[0]) * 0.33333333333;
		data[1] = ( v0.data[1] + v1.data[1] + v2.data[1]) * 0.33333333333;
		data[2] = ( v0.data[2] + v1.data[2] + v2.data[2]) * 0.33333333333;
	}
	
	inline void SetGravityCenter(const TVector3 &v0, const  TVector3 &v1, const TVector3 &v2, const TVector3 &v3){
		data[0] = ( v0.data[0] + v1.data[0] + v2.data[0] + v3.data[0]) * 0.25;
		data[1] = ( v0.data[1] + v1.data[1] + v2.data[1] + v3.data[1]) * 0.25;
		data[2] = ( v0.data[2] + v1.data[2] + v2.data[2] + v3.data[2]) * 0.25;
	}

	inline void SetAdditionWithCoef(double c1, const TVector3 &v1, double c2, const TVector3 &v2){
		data[0] = c1 * v1.data[0] + c2 * v2.data[0];
		data[1] = c1 * v1.data[1] + c2 * v2.data[1];
		data[2] = c1 * v1.data[2] + c2 * v2.data[2];
	}
	inline void SetAdditionWithCoef(double c1, const TVector3 &v1, 
		                            double c2, const TVector3 &v2,
									double c3, const TVector3 &v3){
		data[0] = c1 * v1.data[0] + c2 * v2.data[0] + c3 * v3.data[0];
		data[1] = c1 * v1.data[1] + c2 * v2.data[1] + c3 * v3.data[1];
		data[2] = c1 * v1.data[2] + c2 * v2.data[2] + c3 * v3.data[2];
	}
	inline void SetAdditionWithCoef(double c1, const TVector3 &v1, 
		                            double c2, const TVector3 &v2,
									double c3, const TVector3 &v3,
									double c4, const TVector3 &v4){
		data[0] = c1 * v1.data[0] + c2 * v2.data[0] + c3 * v3.data[0] + c4 * v4.data[0];
		data[1] = c1 * v1.data[1] + c2 * v2.data[1] + c3 * v3.data[1] + c4 * v4.data[1];
		data[2] = c1 * v1.data[2] + c2 * v2.data[2] + c3 * v3.data[2] + c4 * v4.data[2];
	}
	inline void Set_V1subtV2_crsprd_V3subtV4( const TVector3 &v1, const TVector3 &v2, 
											  const TVector3 &v3, const TVector3 &v4)
	{
		data[0] = (v1.data[1] - v2.data[1]) * (v3.data[2] - v4.data[2])   -   (v1.data[2] - v2.data[2]) * (v3.data[1] - v4.data[1]);
		data[1] = (v1.data[2] - v2.data[2]) * (v3.data[0] - v4.data[0])   -   (v1.data[0] - v2.data[0]) * (v3.data[2] - v4.data[2]);
		data[2] = (v1.data[0] - v2.data[0]) * (v3.data[1] - v4.data[1])   -   (v1.data[1] - v2.data[1]) * (v3.data[0] - v4.data[0]);	
	}
	inline void t_swap( TVector3 &trgt){ std::swap( data[0], trgt.data[0]);
										 std::swap( data[1], trgt.data[1]);
										 std::swap( data[2], trgt.data[2]); }
	///////////FILE MANIPULATION/////////////////////////////////////////////////////////////////////////////////
	inline void writeToFile ( FILE *fp )const { fprintf( fp, "%f %f %f\n" ,  data[0],  data[1],  data[2]); }
	inline void readFromFile( FILE *fp )      { fscanf(  fp, "%lf%lf%lf\n", &data[0], &data[1], &data[2]); }

	inline void writeToFile_normalize( FILE *fp, double sizeX, double sizeY, double sizeZ ){
		fprintf( fp, "%f %f %f\n", data[0]/sizeX, data[1]/sizeY, data[2]/sizeZ); 
	}
	inline void readFromFile_stretch( FILE *fp, double sizeX, double sizeY, double sizeZ  ){
		fscanf( fp, "%lf%lf%lf\n", &data[0], &data[1], &data[2]); 
		data[0] *= sizeX;
		data[1] *= sizeY;
		data[2] *= sizeZ;
	}
	inline void writeToFile_b( FILE *fp)const{
		fwrite( &data[0], sizeof(double), 1, fp );
		fwrite( &data[1], sizeof(double), 1, fp );
		fwrite( &data[2], sizeof(double), 1, fp );
	}
	inline void readFromFile_b( FILE *fp){ fread( data, sizeof(double), 3, fp ); }

	inline void writeToFile_b_normalize( FILE *fp, double sizeX, double sizeY, double sizeZ )const{
		double d;
		d = data[0] / sizeX; fwrite( &d, sizeof(double), 1, fp );
		d = data[1] / sizeY; fwrite( &d, sizeof(double), 1, fp );
		d = data[2] / sizeZ; fwrite( &d, sizeof(double), 1, fp );
	}
	inline void readFromFile_b_stretch( FILE *fp, double sizeX, double sizeY, double sizeZ  ){
		fread( data, sizeof(double), 3, fp );
		data[0] *= sizeX;
		data[1] *= sizeY;
		data[2] *= sizeZ;
	}
};



class TVector2 
{
public:
	double data[2];
	//Constructors////////////////////////////////////////////////////////////////
	TVector2(double f1=0,double f2=0){ Set(f1,f2); }
	TVector2(double f[2]            ){ Set(f)    ; }	
	TVector2(const TVector2 &src    ){ Set(src)  ; }	

	inline void Set(double f1=0.0,double f2=0.0){ data[0]=f1;data[1]=f2; }
	inline void Set(const double* f            ){ if( f ) memcpy(data, f     ,sizeof(double)*2); }
	inline void Set(const TVector2 &v          ){         memcpy(data, v.data,sizeof(double)*2); }

	inline   double& operator[](int id)       { return data[id]; }
	inline   double  operator[](int id) const { return data[id]; }
	inline operator double*(){ return data; }
	inline operator void  *(){ return this; }
	inline TVector2& operator=(const TVector2& src) { Set( src ); return *this; }

	inline double Length_sq()const { return data[0]*data[0] + data[1]*data[1]; }
	inline double Length   ()const { return sqrt( Length_sq() );               }

	inline bool Normalize_Self(){
		double l = Length(); if(l==0.0f) return false;
		data[0] /= l;
		data[1] /= l;
		return true;
	}

	inline void Trace() const { fprintf(stderr,"%f %f\n",data[0],data[1]); }
	inline        TVector2 operator* (double f) const{ return TVector2( f*data[0], f*data[1] ); }
	inline friend TVector2 operator* (double f, const TVector2& vec1){ return vec1*f; }
	inline        TVector2 operator+ (const TVector2& vec1) const {
		return TVector2( data[0] + vec1.data[0], data[1] + vec1.data[1]);
	}
	inline        TVector2 operator- (const TVector2& vec1) const {
		return TVector2( data[0] - vec1.data[0], data[1] - vec1.data[1]);
	}
	inline        TVector2& operator*=(double f){ data[0] *= f;
												  data[1] *= f;
												  return *this ;}
	inline double operator*(const TVector2& v) const { return data[0]*v.data[0] + data[1]*v.data[1];}
	inline TVector2& operator+=(const TVector2& vec1){
		data[0] += vec1.data[0];
		data[1] += vec1.data[1];
		return *this ;
	}


};



/*////////////////////////////////////////////////////////////////////////////////
data[16] spec (Same as OpenGL)

matrix =	data[0]  data[4]  data[8]  data[12]
			data[1]  data[5]  data[9]  data[13]
			data[2]  data[6]  data[10] data[14]
			data[3]  data[7]  data[11] data[15]
////////////////////////////////////////////////////////////////////////////////*/

class TMatrix16 {
public:
	double data[16];

	TMatrix16(){SetIdentity();}
	TMatrix16(double* cmat){Set(cmat);}
	
	//Setters////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	inline void Set(const double* cmat){ memcpy(data,cmat,sizeof(double)*16);}
	inline void Set(const float*  cmat){ for( int i=0;i<16;++i ) data[i]=cmat[i];}
	
	inline void SetZeroMatrix(){memset(data,0,sizeof(double)*16);}
	inline void SetIdentity  (){SetZeroMatrix(); data[0] = data[5] = data[10] = data[15] = 1.0;}
	inline void SetAsTranslateMatrix(const TVector3& trans){ SetIdentity(); memcpy(&data[12],trans.data,sizeof(double)*3); }
	inline void SetAsZoomMatrix(double xzoom,double yzoom,double zzoom){ SetZeroMatrix(); data[0]=xzoom; data[5]=yzoom; data[10]=zzoom; data[15]=1;}
	
	//Operator "=", "*","+" //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	inline TMatrix16& operator=(const TMatrix16 src){memcpy(data, src.data, sizeof(double)*16 ); return *this;}
	inline TMatrix16& operator=(const double*   src){memcpy(data, src     , sizeof(double)*16 ); return *this;}

	inline void operator*=(const TMatrix16& mat){
		double tmp[16]; memset(tmp,0,sizeof(double)*16) ;
		for(int i=0;i<4;i++)for(int j=0;j<4;j++)for(int k=0;k<4;k++)  tmp[j*4+i] += data[k*4+i] * mat.data[j*4+k];
		memcpy( data,tmp,sizeof(double)*16 ) ;
	}
	inline void operator*=(const double mat[16]){
		double tmp[16]; memset(tmp,0,sizeof(double)*16) ;
		for(int i=0;i<4;i++)for(int j=0;j<4;j++)for(int k=0;k<4;k++)  tmp[j*4+i] += data[k*4+i] * mat[j*4+k];
		memcpy( data,tmp,sizeof(double)*16 ) ;
	}
	inline TMatrix16 operator*(const TMatrix16& mat){
		TMatrix16 tmp ; memset(tmp.data,0,sizeof(double)*16) ;
		for(int i=0;i<4;i++)for(int j=0;j<4;j++)for(int k=0;k<4;k++) tmp[i*4+j]+=(*this)[k*4+j]*mat.data[i*4+k];
		return tmp ;
	}

	inline TVector3 operator*(const TVector3& v) const {
		return TVector3( data[0]*v.data[0] + data[4]*v.data[1] + data[ 8]*v.data[2] + data[12],
						 data[1]*v.data[0] + data[5]*v.data[1] + data[ 9]*v.data[2] + data[13],
						 data[2]*v.data[0] + data[6]*v.data[1] + data[10]*v.data[2] + data[14] ) ;
	}

	inline void operator+=(TMatrix16& mat){ for(int i=0;i<16;i++) data[i] += mat.data[i]; }
	inline TMatrix16 operator+(const TMatrix16& mat){
		TMatrix16 tmp ;
		for(int i=0;i<16;i++) tmp[i]= data[i] + mat.data[i];
		return tmp ;
	}

	//Operator others//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	inline double& operator[](int i){ return data[i]; }
	inline operator double*(){        return data; }

	TMatrix16 GetTransposeMatrix(){
		static const int map_mat[16]={
		0,4,8,12 ,1,5,9,13 ,2,6,10,14 ,3,7,11,15 } ;
		TMatrix16 ret ; for( int i=0;i<16;i++ )ret[i]=data[map_mat[i]]; return ret ; 
	}

	//inversions////////////////////////////////////////////////////////////////////////////////////////////////////
	TMatrix16 GetInvMatrix(){
		if( fabs(SquareMat_CalcDeterminant(data,4))<VSN )	return false ;
		TMatrix16 ret ;
		SquareMat_CalcInvMatrix(GetTransposeMatrix().data,ret.data,4) ;
		return ret.GetTransposeMatrix() ;
	}

	bool GetInvMatrix_Self(){
		if( fabs(SquareMat_CalcDeterminant(data,4))<VSN )	return false ;
		TMatrix16 tmp ;
		SquareMat_CalcInvMatrix(GetTransposeMatrix().data,tmp.data,4) ;
		memcpy( data,tmp.GetTransposeMatrix().data,sizeof(double)*16 ) ;
		return true ;
	}
	//rotations////////////////////////////////////////////////////////////////////////////////////////////////////
	void SetAsRotateX(double angle){
		SetIdentity();
		double c=(double)cos((double)angle);
		double s=(double)sin((double)angle);
		data[5]=c; data[9]=-s;
		data[6]=s; data[10]=c;
	}
	void SetAsRotateY(double angle){
		SetIdentity();
		double c=(double)cos((double)angle);
		double s=(double)sin((double)angle);
		data[0]= c;  data[8]=s;	
		data[2]=-s; data[10]=c;
	}
	void SetAsRotateZ(double angle){
		SetIdentity();
		double c=(double)cos((double)angle);
		double s=(double)sin((double)angle);
		data[0]=c;  data[4]=-s;	
		data[1]=s;	data[5]= c;
	}
	bool RotateAlongArbitraryAxis(const TVector3& axis, double cosangle,double sinangle)
	{
		if(cosangle==1.0f && sinangle==0.0f){ SetIdentity(); return true ;}
		TVector3 ax=axis;
		if(!ax.Normalize_Self()            ){ SetIdentity(); return false;}

		double u,v,w;
		u=ax.data[0];
		v=ax.data[1];
		w=ax.data[2];

		SetIdentity();
		data[0]=u*u+(1-u*u)*cosangle;
		data[4]=u*v*(1-cosangle)-w*sinangle;
		data[8]=u*w*(1-cosangle)+v*sinangle;
		data[1]=u*v*(1-cosangle)+w*sinangle;
		data[5]=v*v+(1-v*v)*cosangle;
		data[9]=v*w*(1-cosangle)-u*sinangle;
		data[2]=u*w*(1-cosangle)-v*sinangle;
		data[6]=v*w*(1-cosangle)+u*sinangle;
		data[10]=w*w+(1-w*w)*cosangle;
		return true;
	}
	bool RotateAlongArbitraryAxis(const TVector3& axis, double angle)                         { return RotateAlongArbitraryAxis(axis       ,cos(angle), sin(angle)); }
	bool RotateAlongArbitraryAxis(double x,double y,double z, double cosangle,double sinangle){ return RotateAlongArbitraryAxis(TVector3(x,y,z), cosangle,sinangle); }
	bool RotateAlongArbitraryAxis(double x,double y,double z,double angle){                     return RotateAlongArbitraryAxis(TVector3(x,y,z), angle);             }
	void RotateAlongArbitraryAxis_WithShift(TVector3& axis_ori, TVector3& axis_pos, double angle){
		TMatrix16 mat1,rot,mat2 ;
		mat1[12] = -axis_pos[0] ;   mat2[12] =  axis_pos[0] ;
		mat1[13] = -axis_pos[1] ;   mat2[13] =  axis_pos[1] ;
		mat1[14] = -axis_pos[2] ;   mat2[14] =  axis_pos[2] ;
		rot.RotateAlongArbitraryAxis(axis_ori, angle) ;
		*this = mat2 * rot * mat1 ;
	}
	bool RotateFromVectorToVector(TVector3 src,TVector3 targ){
		src .Normalize_Self() ; 
		targ.Normalize_Self() ;
		TVector3 cp = src^targ ;
		return RotateAlongArbitraryAxis( cp, src*targ, cp.Length()) ;
	}

	//Get rotation axis and rotation angle from a Rotation matrix
	//This function works well only when the target matrix is rotation matrix (A == -At)
	void GetRotCenterAxisAndAngle(TVector3& cent,double& angle){
		double Q[4] ;	// quaternion
		{
			TMatrix16 tp = GetTransposeMatrix() ;
			double T = tp[0] + tp[5] + tp[10] + 1 ;	// Trace
		
			if( T > 0 ){
				double S = 0.5 / sqrt(T) ;
				Q[0] = ( tp[9] - tp[6] ) * S ;
				Q[1] = ( tp[2] - tp[8] ) * S ;
				Q[2] = ( tp[4] - tp[1] ) * S ;
				Q[3] = 0.25 / S ;
			} else {
				switch(MAX3ID(tp[0],tp[5],tp[10])){
				case 0:
					{
						double S = sqrt( 1.0 + tp[0] - tp[5] - tp[10] ) * 2;
						Q[0] = 0.5 / S;
						Q[1] = (tp[1] + tp[4] ) / S;
						Q[2] = (tp[2] + tp[8] ) / S;
						Q[3] = (tp[6] + tp[9] ) / S;
					} break ;
				case 1:
					{
						double S  = sqrt( 1.0 + tp[5] - tp[0] - tp[10] ) * 2;
						Q[0] = (tp[1] + tp[4] ) / S;
						Q[1] = 0.5 / S;
						Q[2] = (tp[6] + tp[9] ) / S;
						Q[3] = (tp[2] + tp[8] ) / S;
					} break ;
				case 2:
					{
						double S  = sqrt( 1.0 + tp[10] - tp[0] - tp[5] ) * 2;
						Q[0] = (tp[2] + tp[8] ) / S;
						Q[1] = (tp[6] + tp[9] ) / S;
						Q[2] = 0.5 / S;
						Q[3] = (tp[1] + tp[4] ) / S;
					} break ;
				}
			}
		}
		angle = acos(Q[3])*2 ;
		cent.Set(Q[0],Q[1],Q[2]) ;
		cent.Normalize_Self() ;
	}

	void ComposeFromRotZoomTranslate_byNaiveMethod(TVector3& rotAxis,double rotAngle,double zoomFactor,TVector3& translate){
		TMatrix16 transMat , zoomMat , rotateMat ;
		transMat.SetAsTranslateMatrix(translate) ;
		zoomMat.SetAsZoomMatrix(zoomFactor,zoomFactor,zoomFactor) ;
		rotateMat.RotateAlongArbitraryAxis(rotAxis,rotAngle) ;
		*this = transMat * zoomMat * rotateMat ;
	}
	
	inline void t_GetPosition( TVector3 &pos ) const{
		pos.data[0] = data[12];
		pos.data[1] = data[13];
		pos.data[2] = data[14];
	}

	//Mat * TVector3( 0, l, 0)
	inline void t_GetPositionOnY( double l, TVector3 &pos ) const{
		pos.data[0] = l * data[4] + data[12];
		pos.data[1] = l * data[5] + data[13];
		pos.data[2] = l * data[6] + data[14];
	};
	//M*TVector3(0,1,0) 
	inline void t_GetYDirection( TVector3 &yDir ) const{
		yDir.data[0] = data[4];	
		yDir.data[1] = data[5];	
		yDir.data[2] = data[6];	
	}
	//return M * ILVector(0,0,1)
	inline void t_GetZDirection( TVector3 &zDir ) const {
		zDir.data[0] = data[8 ];
		zDir.data[1] = data[9 ];
		zDir.data[2] = data[10];
	}
	inline void SetTransTerm( const TVector3 &v){ data[12] = v.data[0];  data[13] = v.data[1]; data[14] = v.data[2];}
	inline void SetTranslateValueZero(){ data[12] = data[13] = data[14] = 0; }

	void Trace() const {
		for( int i=0;i<4;i++ ){
			for( int j=0;j<4;j++ ) fprintf(stderr, "%f ",data[j*4+i]) ;
			fprintf(stderr,"\n") ;
		}
	}
};

/*/////////////////////////////////////////////////////////////////////////
data[9] spec (Same as OpenGL)
matrix =data[0]  data[3]  data[6]
		data[1]  data[4]  data[7]
		data[2]  data[5]  data[8] 
/////////////////////////////////////////////////////////////////////////*/


class TMatrix9 {
public:
	double data[9];
	//Constructers/////////////////////////////////////////////////
	TMatrix9(){SetIdentity();}
	TMatrix9(double* cmat){Set(cmat);}
	TMatrix9(double m11, double m12, double m13,
			 double m21, double m22, double m23,
			 double m31, double m32, double m33){
		Set( m11, m12, m13, 
			 m21, m22, m23, 
			 m31, m32, m33);
	}

	//Setters//////////////////////////////////////////////////////
	inline void SetZero()    { memset(data,0,sizeof(double)*9);}
	inline void SetIdentity(){ memset(data,0,sizeof(double)*9); data[0]=1.0f; data[4]=1.0f; data[8]=1.0f; }

	inline void Set(double* cmat){memcpy(data,cmat,sizeof(double)*9);}
	inline void Set(float*  cmat){for( int i=0;i<9;++i ) data[i]=cmat[i];}

	inline void Set(const double &m11, const double &m12, const double &m13,
			        const double &m21, const double &m22, const double &m23,
			        const double &m31, const double &m32, const double &m33)
	{
		data[0] = m11;
		data[1] = m21;
		data[2] = m31;

		data[3] = m12;
		data[4] = m22;
		data[5] = m32;
		
		data[6] = m13;
		data[7] = m23;
		data[8] = m33;
	}

	////Transpose/////////////////////////////////////////////////////
	inline TMatrix9 getTransposed(){
		return TMatrix9(data[0], data[1], data[2],
						data[3], data[4], data[5],
						data[6], data[7], data[8]);
	}
	inline void setTransposed(const TMatrix9& M){
		data[0] = M.data[0]; data[3] = M.data[1], data[6] = M.data[2],
		data[1] = M.data[3]; data[4] = M.data[4], data[7] = M.data[5],
		data[2] = M.data[6]; data[5] = M.data[7], data[8] = M.data[8];
	}
	void transpose_Self(){
		std::swap(data[1], data[3]);
		std::swap(data[2], data[6]);
		std::swap(data[5], data[7]);
	}


	void multScholar(double a){ for(int i = 0; i < 9; ++i) data[i] *= a; }
	

	inline       double& operator[](int i)       { return data[i]; }
	inline const double  operator[](int i) const { return data[i]; }
	operator double*(){ return data; }

	void operator*=(const TMatrix9& mat){
		double tmp[9];
		tmp[0] = data[0] * mat.data[0] + data[3] * mat.data[1] + data[6] * mat.data[2];
		tmp[1] = data[1] * mat.data[0] + data[4] * mat.data[1] + data[7] * mat.data[2];
		tmp[2] = data[2] * mat.data[0] + data[5] * mat.data[1] + data[8] * mat.data[2];

		tmp[3] = data[0] * mat.data[3] + data[3] * mat.data[4] + data[6] * mat.data[5];
		tmp[4] = data[1] * mat.data[3] + data[4] * mat.data[4] + data[7] * mat.data[5];
		tmp[5] = data[2] * mat.data[3] + data[5] * mat.data[4] + data[8] * mat.data[5];
		
		tmp[6] = data[0] * mat.data[6] + data[3] * mat.data[7] + data[6] * mat.data[8];
		tmp[7] = data[1] * mat.data[6] + data[4] * mat.data[7] + data[7] * mat.data[8];
		tmp[8] = data[2] * mat.data[6] + data[5] * mat.data[7] + data[8] * mat.data[8];
		memcpy( data,tmp,sizeof(double)*9 ) ;
	}

	TMatrix9 operator*(const TMatrix9& mat){
		double tmp[9] ;
		tmp[0] = data[0] * mat.data[0] + data[3] * mat.data[1] + data[6] * mat.data[2];
		tmp[1] = data[1] * mat.data[0] + data[4] * mat.data[1] + data[7] * mat.data[2];
		tmp[2] = data[2] * mat.data[0] + data[5] * mat.data[1] + data[8] * mat.data[2];

		tmp[3] = data[0] * mat.data[3] + data[3] * mat.data[4] + data[6] * mat.data[5];
		tmp[4] = data[1] * mat.data[3] + data[4] * mat.data[4] + data[7] * mat.data[5];
		tmp[5] = data[2] * mat.data[3] + data[5] * mat.data[4] + data[8] * mat.data[5];
		
		tmp[6] = data[0] * mat.data[6] + data[3] * mat.data[7] + data[6] * mat.data[8];
		tmp[7] = data[1] * mat.data[6] + data[4] * mat.data[7] + data[7] * mat.data[8];
		tmp[8] = data[2] * mat.data[6] + data[5] * mat.data[7] + data[8] * mat.data[8];
		return TMatrix9(tmp) ;
	}

	TVector3 operator*(const TVector3& v)const{
		return TVector3(data[0] * v.data[0] + data[3] * v.data[1] + data[6] * v.data[2],
						data[1] * v.data[0] + data[4] * v.data[1] + data[7] * v.data[2],
						data[2] * v.data[0] + data[5] * v.data[1] + data[8] * v.data[2]) ;
	}
	TMatrix9 operator+(const TMatrix9& mat){
		return TMatrix9( data[0] + mat.data[0], data[3] + mat.data[3], data[6] + mat.data[6],
						 data[1] + mat.data[1], data[4] + mat.data[4], data[7] + mat.data[7],
						 data[2] + mat.data[2], data[5] + mat.data[5], data[8] + mat.data[8]);
	}
	void operator+=(TMatrix9& mat){
		data[0]+=mat.data[0]; data[3]+=mat.data[3]; data[6]+=mat.data[6];
		data[1]+=mat.data[1]; data[4]+=mat.data[4]; data[7]+=mat.data[7];
		data[2]+=mat.data[2]; data[5]+=mat.data[5]; data[8]+=mat.data[8];
	}

	inline friend TMatrix9 operator*(double f,const TMatrix9& mat1){
		TMatrix9 tmp;
		for(int i = 0; i < 9;++i) tmp[i] = mat1.data[i] * f;
		return tmp;
	}
	
	TMatrix9& operator=(const TMatrix9& src){memcpy(data,src.data,sizeof(double)*9);return*this;}
	TMatrix9& operator=(const double* src  ){memcpy(data,src     ,sizeof(double)*9);return*this;}

	void Trace(){
		for( int i=0;i<3;++i ){
			for( int j=0;j<3;++j ) fprintf(stderr,"%f ",data[j*3+i]) ;	
			fprintf(stderr, "\n") ;
		}
	}

	// M = l * l.transpose
	inline void setMultVectors(const TVector3& l, const TVector3& r){
		data[0] = l.data[0] * r.data[0];   data[3] = l.data[0] * r.data[1];   data[6] = l.data[0] * r.data[2];
		data[1] = l.data[1] * r.data[0];   data[4] = l.data[1] * r.data[1];   data[7] = l.data[1] * r.data[2];
		data[2] = l.data[2] * r.data[0];   data[5] = l.data[2] * r.data[1];   data[8] = l.data[2] * r.data[2];
	}
	// M += l * r.transpose
	inline void addMultVectors(const TVector3& l, const TVector3& r){
		data[0] += l.data[0] * r.data[0];   data[3] += l.data[0] * r.data[1];   data[6] += l.data[0] * r.data[2];
		data[1] += l.data[1] * r.data[0];   data[4] += l.data[1] * r.data[1];   data[7] += l.data[1] * r.data[2];
		data[2] += l.data[2] * r.data[0];   data[5] += l.data[2] * r.data[1];   data[8] += l.data[2] * r.data[2];
	}
	// M += l * r.transpose
	inline void subtMultVectors(const TVector3& l, const TVector3& r){
		data[0] -= l.data[0] * r.data[0];   data[3] -= l.data[0] * r.data[1];   data[6] -= l.data[0] * r.data[2];
		data[1] -= l.data[1] * r.data[0];   data[4] -= l.data[1] * r.data[1];   data[7] -= l.data[1] * r.data[2];
		data[2] -= l.data[2] * r.data[0];   data[5] -= l.data[2] * r.data[1];   data[8] -= l.data[2] * r.data[2];
	}

	// M = coef * (l * r.transpose)
	inline void setMultVectors_withCoef(double c, const TVector3& l, const TVector3& r){
		data[0] = c * l.data[0] * r.data[0];   data[3] = c * l.data[0] * r.data[1];   data[6] = c * l.data[0] * r.data[2];
		data[1] = c * l.data[1] * r.data[0];   data[4] = c * l.data[1] * r.data[1];   data[7] = c * l.data[1] * r.data[2];
		data[2] = c * l.data[2] * r.data[0];   data[5] = c * l.data[2] * r.data[1];   data[8] = c * l.data[2] * r.data[2];
	}
	// M += coef * (l * r.transpose)
	inline void addMultVectors_withCoef(double c, const TVector3& l, const TVector3& r){
		data[0] += c * l.data[0] * r.data[0];   data[3] += c * l.data[0] * r.data[1];   data[6] += c * l.data[0] * r.data[2];
		data[1] += c * l.data[1] * r.data[0];   data[4] += c * l.data[1] * r.data[1];   data[7] += c * l.data[1] * r.data[2];
		data[2] += c * l.data[2] * r.data[0];   data[5] += c * l.data[2] * r.data[1];   data[8] += c * l.data[2] * r.data[2];
	}

	// M += d * M 
	inline void setMatMultScholar(double d, const TMatrix9& M){
		data[0] = d * M.data[0];  data[3] = d * M.data[3];  data[6] = d * M.data[6];
		data[1] = d * M.data[1];  data[4] = d * M.data[4];  data[7] = d * M.data[7];
		data[2] = d * M.data[2];  data[5] = d * M.data[5];  data[8] = d * M.data[8];
	}
	inline void addMatMultScholar(double d, const TMatrix9& M){
		data[0] += d * M.data[0];  data[3] += d * M.data[3];  data[6] += d * M.data[6];
		data[1] += d * M.data[1];  data[4] += d * M.data[4];  data[7] += d * M.data[7];
		data[2] += d * M.data[2];  data[5] += d * M.data[5];  data[8] += d * M.data[8];
	}
	//M = M1 * M2
	inline void setMultM1M2(const TMatrix9& M1, const TMatrix9& M2)
	{
		data[0] = M1.data[0] * M2.data[0]  +  M1.data[3] * M2.data[1]  +  M1.data[6] * M2.data[2];
		data[1] = M1.data[1] * M2.data[0]  +  M1.data[4] * M2.data[1]  +  M1.data[7] * M2.data[2];
		data[2] = M1.data[2] * M2.data[0]  +  M1.data[5] * M2.data[1]  +  M1.data[8] * M2.data[2];

		data[3] = M1.data[0] * M2.data[3]  +  M1.data[3] * M2.data[4]  +  M1.data[6] * M2.data[5];
		data[4] = M1.data[1] * M2.data[3]  +  M1.data[4] * M2.data[4]  +  M1.data[7] * M2.data[5];
		data[5] = M1.data[2] * M2.data[3]  +  M1.data[5] * M2.data[4]  +  M1.data[8] * M2.data[5];
		
		data[6] = M1.data[0] * M2.data[6]  +  M1.data[3] * M2.data[7]  +  M1.data[6] * M2.data[8];
		data[7] = M1.data[1] * M2.data[6]  +  M1.data[4] * M2.data[7]  +  M1.data[7] * M2.data[8];
		data[8] = M1.data[2] * M2.data[6]  +  M1.data[5] * M2.data[7]  +  M1.data[8] * M2.data[8];
	}
	//M = M1 * M2^t
	inline void setMultM1M2t(const TMatrix9& M1, const TMatrix9& M2t)
	{
		data[0] = M1.data[0] * M2t.data[0]  +  M1.data[3] * M2t.data[3]  +  M1.data[6] * M2t.data[6];
		data[1] = M1.data[1] * M2t.data[0]  +  M1.data[4] * M2t.data[3]  +  M1.data[7] * M2t.data[6];
		data[2] = M1.data[2] * M2t.data[0]  +  M1.data[5] * M2t.data[3]  +  M1.data[8] * M2t.data[6];

		data[3] = M1.data[0] * M2t.data[1]  +  M1.data[3] * M2t.data[4]  +  M1.data[6] * M2t.data[7];
		data[4] = M1.data[1] * M2t.data[1]  +  M1.data[4] * M2t.data[4]  +  M1.data[7] * M2t.data[7];
		data[5] = M1.data[2] * M2t.data[1]  +  M1.data[5] * M2t.data[4]  +  M1.data[8] * M2t.data[7];
		
		data[6] = M1.data[0] * M2t.data[2]  +  M1.data[3] * M2t.data[5]  +  M1.data[6] * M2t.data[8];
		data[7] = M1.data[1] * M2t.data[2]  +  M1.data[4] * M2t.data[5]  +  M1.data[7] * M2t.data[8];
		data[8] = M1.data[2] * M2t.data[2]  +  M1.data[5] * M2t.data[5]  +  M1.data[8] * M2t.data[8];
	}
	
	//M = M1_transpose * M2
	inline void setMultM1tM2(const TMatrix9& M1t, const TMatrix9& M2)
	{
		data[0] = M1t.data[0] * M2.data[0]  +  M1t.data[1] * M2.data[1]  +  M1t.data[2] * M2.data[2];
		data[1] = M1t.data[3] * M2.data[0]  +  M1t.data[4] * M2.data[1]  +  M1t.data[5] * M2.data[2];
		data[2] = M1t.data[6] * M2.data[0]  +  M1t.data[7] * M2.data[1]  +  M1t.data[8] * M2.data[2];

		data[3] = M1t.data[0] * M2.data[3]  +  M1t.data[1] * M2.data[4]  +  M1t.data[2] * M2.data[5];
		data[4] = M1t.data[3] * M2.data[3]  +  M1t.data[4] * M2.data[4]  +  M1t.data[5] * M2.data[5];
		data[5] = M1t.data[6] * M2.data[3]  +  M1t.data[7] * M2.data[4]  +  M1t.data[8] * M2.data[5];
		
		data[6] = M1t.data[0] * M2.data[6]  +  M1t.data[1] * M2.data[7]  +  M1t.data[2] * M2.data[8];
		data[7] = M1t.data[3] * M2.data[6]  +  M1t.data[4] * M2.data[7]  +  M1t.data[5] * M2.data[8];
		data[8] = M1t.data[6] * M2.data[6]  +  M1t.data[7] * M2.data[7]  +  M1t.data[8] * M2.data[8];
	}
	//M = M1_transpose * M1	
	inline void setMultM1tM1(const TMatrix9& M1)
	{
		data[0] = M1.data[0] * M1.data[0]  +  M1.data[1] * M1.data[1]  +  M1.data[2] * M1.data[2];
		data[1] = M1.data[3] * M1.data[0]  +  M1.data[4] * M1.data[1]  +  M1.data[5] * M1.data[2];
		data[2] = M1.data[6] * M1.data[0]  +  M1.data[7] * M1.data[1]  +  M1.data[8] * M1.data[2];
		data[3] = data[1];
		data[4] = M1.data[3] * M1.data[3]  +  M1.data[4] * M1.data[4]  +  M1.data[5] * M1.data[5];
		data[5] = M1.data[6] * M1.data[3]  +  M1.data[7] * M1.data[4]  +  M1.data[8] * M1.data[5];
		data[6] = data[2];
		data[7] = data[5];
		data[8] = M1.data[6] * M1.data[6]  +  M1.data[7] * M1.data[7]  +  M1.data[8] * M1.data[8];
	}
	//M = M1 * M2 * M3_transpose
	inline void setMultM1_M2_M3t(const TMatrix9& L, const TMatrix9& C, const TMatrix9& R)
	{
		double T[9];
		T[0] = C.data[0]*R.data[0] + C.data[3]*R.data[3] + C.data[6]*R.data[6];
		T[1] = C.data[1]*R.data[0] + C.data[4]*R.data[3] + C.data[7]*R.data[6];
		T[2] = C.data[2]*R.data[0] + C.data[5]*R.data[3] + C.data[8]*R.data[6];
		T[3] = C.data[0]*R.data[1] + C.data[3]*R.data[4] + C.data[6]*R.data[7];
		T[4] = C.data[1]*R.data[1] + C.data[4]*R.data[4] + C.data[7]*R.data[7];
		T[5] = C.data[2]*R.data[1] + C.data[5]*R.data[4] + C.data[8]*R.data[7];
		T[6] = C.data[0]*R.data[2] + C.data[3]*R.data[5] + C.data[6]*R.data[8];
		T[7] = C.data[1]*R.data[2] + C.data[4]*R.data[5] + C.data[7]*R.data[8];
		T[8] = C.data[2]*R.data[2] + C.data[5]*R.data[5] + C.data[8]*R.data[8];

		data[0] = L.data[0]*T[0] + L.data[3]*T[1] + L.data[6]*T[2];
		data[1] = L.data[1]*T[0] + L.data[4]*T[1] + L.data[7]*T[2];
		data[2] = L.data[2]*T[0] + L.data[5]*T[1] + L.data[8]*T[2];

		data[3] = L.data[0]*T[3] + L.data[3]*T[4] + L.data[6]*T[5];
		data[4] = L.data[1]*T[3] + L.data[4]*T[4] + L.data[7]*T[5];
		data[5] = L.data[2]*T[3] + L.data[5]*T[4] + L.data[8]*T[5];

		data[6] = L.data[0]*T[6] + L.data[3]*T[7] + L.data[6]*T[8];
		data[7] = L.data[1]*T[6] + L.data[4]*T[7] + L.data[7]*T[8];
		data[8] = L.data[2]*T[6] + L.data[5]*T[7] + L.data[8]*T[8];
	}

	//M = V * diag(d0,d1,d2) * V_transpose
	inline void setMultV_diag_Vt(const TMatrix9& V, const double &d0, const double &d1, const double &d2)
	{
		double T[9];
		T[0] = d0 * V.data[0];
		T[1] = d1 * V.data[3];
		T[2] = d2 * V.data[6];
		T[3] = d0 * V.data[1];
		T[4] = d1 * V.data[4];
		T[5] = d2 * V.data[7];
		T[6] = d0 * V.data[2];
		T[7] = d1 * V.data[5];
		T[8] = d2 * V.data[8];

		data[0] = V.data[0]*T[0] + V.data[3]*T[1] + V.data[6]*T[2];
		data[1] = V.data[1]*T[0] + V.data[4]*T[1] + V.data[7]*T[2];
		data[2] = V.data[2]*T[0] + V.data[5]*T[1] + V.data[8]*T[2];

		data[3] = data[1];
		data[4] = V.data[1]*T[3] + V.data[4]*T[4] + V.data[7]*T[5];
		data[5] = V.data[2]*T[3] + V.data[5]*T[4] + V.data[8]*T[5];

		data[6] = data[2];
		data[7] = data[5];
		data[8] = V.data[2]*T[6] + V.data[5]*T[7] + V.data[8]*T[8];
	}

	// S(v1) * S(v2).getTransposeを一度に計算
	inline void setS1multS2transpose(TVector3& v1, TVector3& v2)
	{
		data[0] =   v1.data[2] * v2.data[2]  +  v1.data[1] * v2.data[1];
		data[1] = - v1.data[0] * v2.data[1];
		data[2] = - v1.data[0] * v2.data[2];

		data[3] = - v1.data[1] * v2.data[0];
		data[4] =   v1.data[2] * v2.data[2]  +  v1.data[0] * v2.data[0];
		data[5] = - v1.data[1] * v2.data[2];
		
		data[6] = - v1.data[2] * v2.data[0];
		data[7] = - v1.data[2] * v2.data[1];
		data[8] =   v1.data[1] * v2.data[1]  +  v1.data[0] * v2.data[0];
	}

	inline void addS1multS2transpose(TVector3& v1, TVector3& v2)
	{
		data[0] +=   v1.data[2] * v2.data[2]  +  v1.data[1] * v2.data[1];
		data[1] -=   v1.data[0] * v2.data[1];
		data[2] -=   v1.data[0] * v2.data[2];

		data[3] -=  v1.data[1] * v2.data[0];
		data[4] +=  v1.data[2] * v2.data[2]  +  v1.data[0] * v2.data[0];
		data[5] -=  v1.data[1] * v2.data[2];
		
		data[6] -=  v1.data[2] * v2.data[0];
		data[7] -=  v1.data[2] * v2.data[1];
		data[8] +=  v1.data[1] * v2.data[1]  +  v1.data[0] * v2.data[0];
	}

	inline void addMultM1M2(TMatrix9& M1, TMatrix9& M2){
		data[0] += M1.data[0] * M2.data[0]  +  M1.data[3] * M2.data[1]  +  M1.data[6] * M2.data[2];
		data[1] += M1.data[1] * M2.data[0]  +  M1.data[4] * M2.data[1]  +  M1.data[7] * M2.data[2];
		data[2] += M1.data[2] * M2.data[0]  +  M1.data[5] * M2.data[1]  +  M1.data[8] * M2.data[2];

		data[3] += M1.data[0] * M2.data[3]  +  M1.data[3] * M2.data[4]  +  M1.data[6] * M2.data[5];
		data[4] += M1.data[1] * M2.data[3]  +  M1.data[4] * M2.data[4]  +  M1.data[7] * M2.data[5];
		data[5] += M1.data[2] * M2.data[3]  +  M1.data[5] * M2.data[4]  +  M1.data[8] * M2.data[5];
		
		data[6] += M1.data[0] * M2.data[6]  +  M1.data[3] * M2.data[7]  +  M1.data[6] * M2.data[8];
		data[7] += M1.data[1] * M2.data[6]  +  M1.data[4] * M2.data[7]  +  M1.data[7] * M2.data[8];
		data[8] += M1.data[2] * M2.data[6]  +  M1.data[5] * M2.data[7]  +  M1.data[8] * M2.data[8];
	}


	inline double calcDet9() const {
		return   + data[0] * data[4] * data[8] 
			     + data[1] * data[5] * data[6] 
				 + data[2] * data[3] * data[7] 
				  
				 - data[0] * data[5] * data[7] 
				 - data[2] * data[4] * data[6] 
				 - data[1] * data[3] * data[8];
	}

	inline static double calcDet9(const TMatrix9 &Mat){
		return   + Mat.data[0] * Mat.data[4] * Mat.data[8] 
			     + Mat.data[1] * Mat.data[5] * Mat.data[6] 
				 + Mat.data[2] * Mat.data[3] * Mat.data[7] 
				  
				 - Mat.data[0] * Mat.data[5] * Mat.data[7] 
				 - Mat.data[2] * Mat.data[4] * Mat.data[6] 
				 - Mat.data[1] * Mat.data[3] * Mat.data[8];
	}

	inline bool getInvertSelf(){
		double det = calcDet9();
		if(abs(det) < 0.0000001) return false;
		double detInv = 1.0 / det;
		double tmp[9];
		tmp[0] = detInv * ( data[4] * data[8] - data[7] * data[5]);
		tmp[1] = detInv * ( data[7] * data[2] - data[1] * data[8]);
		tmp[2] = detInv * ( data[1] * data[5] - data[4] * data[2]);

		tmp[3] = detInv * ( data[6] * data[5] - data[3] * data[8]);
		tmp[4] = detInv * ( data[0] * data[8] - data[6] * data[2]);
		tmp[5] = detInv * ( data[3] * data[2] - data[0] * data[5]);

		tmp[6] = detInv * ( data[3] * data[7] - data[6] * data[4]);
		tmp[7] = detInv * ( data[6] * data[1] - data[0] * data[7]);
		tmp[8] = detInv * ( data[0] * data[4] - data[3] * data[1]);
		memcpy(data,tmp,sizeof(double)*9);
		return true;	
	}

	//自分の逆行列をtrgtにセット
	inline bool getInvert(TMatrix9 &trgt) const
	{
		double det = calcDet9();

		if(abs(det) < 0.0000001){
			trgt.SetIdentity();
			return false;
		}
		double detInv = 1.0 / det;
		
		trgt.data[0] = detInv * ( data[4] * data[8] - data[7] * data[5]);
		trgt.data[1] = detInv * ( data[7] * data[2] - data[1] * data[8]);
		trgt.data[2] = detInv * ( data[1] * data[5] - data[4] * data[2]);

		trgt.data[3] = detInv * ( data[6] * data[5] - data[3] * data[8]);
		trgt.data[4] = detInv * ( data[0] * data[8] - data[6] * data[2]);
		trgt.data[5] = detInv * ( data[3] * data[2] - data[0] * data[5]);

		trgt.data[6] = detInv * ( data[3] * data[7] - data[6] * data[4]);
		trgt.data[7] = detInv * ( data[6] * data[1] - data[0] * data[7]);
		trgt.data[8] = detInv * ( data[0] * data[4] - data[3] * data[1]);
		return true;
	}

	
	bool RotateAlongArbitraryAxis(double axisX, double axisY, double axisZ, double cosangle,double sinangle){
		double L = sqrt( axisX*axisX + axisY*axisY + axisZ*axisZ );
		if(cosangle==1.0 && sinangle==0.0 ){ SetIdentity(); return true ;}
		if( L == 0.0                      ){ SetIdentity(); return false;}
		axisX/= L;
		axisY/= L;
		axisZ/= L;
		double u = axisX,v=axisY,w=axisZ;
		data[0]=u*u+(1-u*u)*cosangle;        data[3]=u*v*(1-cosangle)-w*sinangle; data[6]=u*w*(1-cosangle)+v*sinangle;
		data[1]=u*v*(1-cosangle)+w*sinangle; data[4]=v*v+(1-v*v)*cosangle;        data[7]=v*w*(1-cosangle)-u*sinangle;
		data[2]=u*w*(1-cosangle)-v*sinangle; data[5]=v*w*(1-cosangle)+u*sinangle; data[8]=w*w+(1-w*w)*cosangle;
		return true;
	}
	bool RotateAlongArbitraryAxis(double axisX, double axisY, double axisZ, double angle){ return RotateAlongArbitraryAxis( axisX  , axisY  , axisZ  , cos(angle), sin(angle) );}
	bool RotateAlongArbitraryAxis(const TVector3 &axis                    , double angle){ return RotateAlongArbitraryAxis( axis[0], axis[1], axis[2], cos(angle), sin(angle) );}
	void RotateAlongX(double angle){
		double c=cos(angle);
		double s=sin(angle);
		data[0]= 1; data[3]= 0; data[6]=0;
		data[1]= 0; data[4]= c; data[7]=-s;
		data[2]= 0;	data[5]= s; data[8]=c;
	}
	void RotateAlongY(double angle){
		double c=cos(angle);
		double s=sin(angle);
		data[0]= c; data[3]= 0; data[6]=s;
		data[1]= 0; data[4]= 1; data[7]=0;
		data[2]=-s;	data[5]= 0; data[8]=c;
	}
	void RotateAlongZ(double angle){
		double c=cos(angle);
		double s=sin(angle);
		data[0]= c; data[3]=-s; data[6]=0;
		data[1]= s; data[4]= c; data[7]=0;
		data[2]= 0;	data[5]= 0; data[8]=1;
	}
	void writeToFile( FILE* fp){
		fprintf( fp, "%f %f %f\n", data[0], data[3], data[6]);
		fprintf( fp, "%f %f %f\n", data[1], data[4], data[7]);
		fprintf( fp, "%f %f %f\n", data[2], data[5], data[8]);
	}
};	

/*////////////////////////////////////////////////////////////////////////////////////////////////////////
data[4] spec (Same as OpenGL)
matrix =
	data[0]  data[2]
	data[1]  data[3]
////////////////////////////////////////////////////////////////////////////////////////////////////////*/
class TMatrix4 {
public:
	double data[4];
	//Constucters/////////////////////////
	TMatrix4(){SetIdentity();}
	TMatrix4(double* cmat){Set(cmat);}
	TMatrix4(double m11, double m12, 
			 double m21, double m22){ Set( m11, m12, m21, m22 );}
	inline void Set(const double &m11, const double &m12, 
			        const double &m21, const double &m22){ data[0] = m11; data[2] = m12;
		                                                   data[1] = m21; data[3] = m22;}
	//Setters/////////////////////////////
	inline void SetZero()    { memset(data, 0  , sizeof(double)*4);                            }
	inline void SetIdentity(){ memset(data, 0  , sizeof(double)*4); data[0]=1.0f; data[3]=1.0f;}
	void Set( double* cmat)  { memcpy(data,cmat, sizeof(double)*4);}
	void Set( float*  cmat) {for( int i=0;i<4;++i ) data[i]=cmat[i];}

	inline void transpose_Self(){ std::swap(data[1], data[2]); }
	void multScholar(double a){ for(int i = 0; i < 4; ++i) data[i] *= a; }
	
	inline double&      operator[](int i){        return data[i]; }
	inline const double operator[](int i) const { return data[i]; }
	operator double*(){ return data; }
	
	inline TVector3 operator*(const TVector3& v){
		return TVector3( data[0] * v.data[0] + data[2] * v.data[1],
			             data[1] * v.data[0] + data[3] * v.data[1], 0);
	}
	inline TMatrix4 operator+(const TMatrix4& mat){
		return TMatrix4( data[0] + mat.data[0], data[2] + mat.data[2],
						 data[1] + mat.data[1], data[3] + mat.data[3] );
	}
	inline TMatrix4 operator-(const TMatrix4& mat){
		return TMatrix4( data[0] - mat.data[0],   data[2] - mat.data[2],
						 data[1] - mat.data[1],   data[3] - mat.data[3] );
	}
	inline TMatrix4 operator*(const TMatrix4& M){
		return TMatrix4( data[0]*M[0] + data[2]*M[1],   data[0]*M[2] + data[2]*M[3],
			             data[1]*M[0] + data[3]*M[1],   data[1]*M[2] + data[3]*M[3] );
	}
	inline TMatrix4 operator-() const { return TMatrix4(-data[0], -data[2], 
		                                                -data[1], -data[3]);}
	TMatrix4& operator=(const TMatrix4& src){ memcpy(data,src.data, sizeof(double)*4); return*this;}
	TMatrix4& operator=(const double*   src){ memcpy(data,src     , sizeof(double)*4); return*this;}

	inline void operator*=(const TMatrix4& mat)
	{
		double tmp[4]; memset(tmp, 0 , sizeof(double)*4);
		tmp[0] = data[0]*mat.data[0] + data[2]*mat.data[1]; tmp[2] = data[0]*mat.data[2] + data[2]*mat.data[3];
		tmp[1] = data[1]*mat.data[0] + data[3]*mat.data[1]; tmp[3] = data[1]*mat.data[2] + data[3]*mat.data[3];
		memcpy( data,tmp,sizeof(double)*4 ) ;
	}
	void operator+=(TMatrix4& mat){
		data[0] += mat.data[0];
		data[1] += mat.data[1];
		data[2] += mat.data[2];
		data[3] += mat.data[3];
	}
	inline void operator/=(double a){
		data[0] /= a;
		data[1] /= a;
		data[2] /= a;
		data[3] /= a;
	}

	//functions//////////////////////////////////////////////////////////
	inline void setSubtM1M2( const TMatrix4 &L, const TMatrix4 &R){
		data[0] = L[0]-R[0];
		data[1] = L[1]-R[1];
		data[2] = L[2]-R[2];
		data[3] = L[3]-R[3];
	}
	inline void setMultM1M2( const TMatrix4 &L, const TMatrix4 &R){
		data[0] = L[0]*R[0] + L[2]*R[1]; data[2] = L[0]*R[2] + L[2]*R[3];
		data[1] = L[1]*R[0] + L[3]*R[1]; data[3] = L[1]*R[2] + L[3]*R[3];	
	}
	inline void setMultM1M2t( const TMatrix4 &L, const TMatrix4 &R){
		data[0] = L[0]*R[0] + L[2]*R[2];  data[2] = L[0]*R[1] + L[2]*R[3];
		data[1] = L[1]*R[0] + L[3]*R[2];  data[3] = L[1]*R[1] + L[3]*R[3];
	}
	inline void setMultM1tM2( const TMatrix4 &L, const TMatrix4 &R){
		data[0] = L[0]*R[0] + L[1]*R[1];  data[2] = L[0]*R[2] + L[1]*R[3];
		data[1] = L[2]*R[0] + L[3]*R[1];  data[3] = L[2]*R[2] + L[3]*R[3];	
	}
	inline void setMultM1tM2M3( const TMatrix4 &M1, const TMatrix4 &M2, const TMatrix4 &M3){
		data[0] = M1.data[0] * (M2.data[0]*M3.data[0] + M2.data[2]*M3.data[1]) + M1.data[1] * (M2.data[1]*M3.data[0] + M2.data[3]*M3.data[1]);  
		data[2] = M1.data[0] * (M2.data[0]*M3.data[2] + M2.data[2]*M3.data[3]) + M1.data[1] * (M2.data[1]*M3.data[2] + M2.data[3]*M3.data[3]); 
		data[1] = M1.data[2] * (M2.data[0]*M3.data[0] + M2.data[2]*M3.data[1]) + M1.data[3] * (M2.data[1]*M3.data[0] + M2.data[3]*M3.data[1]);  
		data[3] = M1.data[2] * (M2.data[0]*M3.data[2] + M2.data[2]*M3.data[3]) + M1.data[3] * (M2.data[1]*M3.data[2] + M2.data[3]*M3.data[3]); 
	}
	inline void setMultM1M2M3t( const TMatrix4 &M1, const TMatrix4 &M2, const TMatrix4 &M3){
		data[0] = (M1.data[0]*M2.data[0] + M1.data[2]*M2.data[1]) * M1.data[0]  +  (M1.data[0]*M2.data[2] + M1.data[2]*M2.data[3]) * M1.data[2];  
		data[2] = (M1.data[0]*M2.data[0] + M1.data[2]*M2.data[1]) * M1.data[1]  +  (M1.data[0]*M2.data[2] + M1.data[2]*M2.data[3]) * M1.data[3];
		data[1] = (M1.data[1]*M2.data[0] + M1.data[3]*M2.data[1]) * M1.data[0]  +  (M1.data[1]*M2.data[2] + M1.data[3]*M2.data[3]) * M1.data[2];
		data[3] = (M1.data[1]*M2.data[0] + M1.data[3]*M2.data[1]) * M1.data[1]  +  (M1.data[1]*M2.data[2] + M1.data[3]*M2.data[3]) * M1.data[3];	
	}
	inline void setV1multV2t( const TVector3 &v1, const TVector3 &v2){
		data[0] = v1.data[0] * v2.data[0];
		data[1] = v1.data[1] * v2.data[0];
		data[2] = v1.data[0] * v2.data[1];
		data[3] = v1.data[1] * v2.data[1];
	}
	inline void addV1multV2t( const TVector3 &v1, const TVector3 &v2){
		data[0] += v1.data[0] * v2.data[0];
		data[1] += v1.data[1] * v2.data[0];
		data[2] += v1.data[0] * v2.data[1];
		data[3] += v1.data[1] * v2.data[1];
	}
	inline void addM1tmultM2( const TMatrix4 &L, const TMatrix4 &R){
		data[0] += L[0]*R[0] + L[1]*R[1];  data[2] += L[0]*R[2] + L[1]*R[3];
		data[1] += L[2]*R[0] + L[3]*R[1];  data[3] += L[2]*R[2] + L[3]*R[3];
	}
	inline bool getInvertSelf(){
		double tmp[4]; memset(tmp, 0 , sizeof(double)*4);
		double det = data[0] * data[3] - data[2] * data[1];
		if( fabs(det) == 0 ) return false;
		det = 1.0 / det;
		tmp[0] = det * data[3];   tmp[2] =-det * data[2];
		tmp[1] =-det * data[1];   tmp[3] = det * data[0];
		memcpy( data,tmp,sizeof(double)*4 ) ;

		return true;
	}

	inline void SetAsRotation(double theta){
		data[0] =  data[3] = cos(theta);
		data[1] = sin(theta);
		data[2] = -data[1];
	}
	inline void Trace(){ fprintf( stderr, "%f %f\n", data[0], data[2]);
		                 fprintf( stderr, "%f %f\n", data[1], data[3]);}
};

////////////////////////////////////////////////////////////////////////////////////////////
//Mathmatical functions (independent to ILVector/Matrix)////////////////////////////////////////////
inline double t_pow(double d,int n){
	if(n < 0){ d = 1/d; n = -n; }
	double result = 1.0;
	for(int i = 0;i < n;++i) result *= d;
	return result;
}


//	  | a b | |s|    w1
//    | c d | |t|  = w2
inline bool t_solve2by2LinearEquation(const double a,  const double b, 
									  const double c,  const double d,  const double w1, const double w2, 
									                                          double &s,       double &t)
{
	double det = (a*d - b*c); 
	if(det == 0) return false;
	det = 1.0 / det;
	s = (  d*w1 - b*w2 ) * det;
	t = ( -c*w1 + a*w2 ) * det;
	return true;
}
inline bool t_isPowerOfTwo(unsigned int val)
{
	unsigned int k = 1;
	for( int i = 0; i < 1024; ++i)
	{
		if( val == k ) return true;
		if( val <  k ) return false;
		k *= 2;
	}
	return false;
}
inline int t_minPowerOfTwo( int n ){
	int val = 1;
	for( int i = 1; i < 63; ++i){
		if( n <= val) return val;
		val *= 2;
	}
	return val;
}
inline double t_cubert( const double &a)
{
	if( a >= 0) return  pow( a, 1.0/3.0);
	else        return -pow(-a, 1.0/3.0);
}

inline void t_solveSquareEq( const double &a, const double &b, const double &c, 
	                                  int &numOfRealSolusion, 
								   double &s1, double &s2)
{
	double D = b*b - 4*a*c;
	if     ( D       < 0            ) numOfRealSolusion = 0;
	else if( fabs(D) <= 0.0000000001){
		numOfRealSolusion = 1;
		s1 = -b / (2*a);
	}else{
		numOfRealSolusion = 2;
		s1 = (-b + sqrt(D) ) / (2*a); //fprintf( stderr, "aaa1   %f\n", a*s1*s1 + b*s1 + c);
		s2 = (-b - sqrt(D) ) / (2*a); //fprintf( stderr, "aaa2   %f\n", a*s2*s2 + b*s2 + c);
	}
}

//solve x^3 + A2 x^2 + A1 x + A0 = 0 
inline void t_solveCubeEq(const double &A2, const double &A1, const double &A0, 
	                      int &NumOfRealSolusion, 
						  double &s1, double &s2, double &s3)
{
	const static double pi_2_3 = M_PI * 2 / 3;
	double p = A1 - A2*A2/3;
	double q = A0 - A1*A2/3 + 2 * A2*A2*A2 / 27;

	double D = q*q/4 + p*p*p / 27;

	if(      D < 0) 
	{
		NumOfRealSolusion = 3;//実数解3個

		//x =  (cos( 60k)+i sin( 60k)) ( a + b i )^(1/3) + (cos(-60k)+i sin(-60k)) ( a - b i )^(1/3); k=1,2,3; 
		double a = -q/2    ; //cos
		double b = sqrt(-D); //sin

		double aabb = sqrt( a*a+b*b );
		double t    = atan2( b/aabb, a/aabb);//cos(t) = a/(aabb)   sin(t) = b/aabb     stan2(sinT, cosT)
		
		double a1 = t_cubert( aabb );
		s1 = a1 * 2 * cos( t/3.0          ) - A2/3.0; //fprintf( stderr, "33res1 %f\n", s1*s1*s1 + A2*s1*s1 + A1*s1 + A0);
		s2 = a1 * 2 * cos( t/3.0 + pi_2_3 ) - A2/3.0; //fprintf( stderr, "33res2 %f\n", s2*s2*s2 + A2*s2*s2 + A1*s2 + A0);
		s3 = a1 * 2 * cos( t/3.0 - pi_2_3 ) - A2/3.0; //fprintf( stderr, "33res3 %f\n", s3*s3*s3 + A2*s3*s3 + A1*s3 + A0);
	}
	else if( fabs(D) <= 0.00000000001)//重解
	{
		if( q = 0 ) 
		{
			NumOfRealSolusion = 1;//3重解
			s1 = -A2 / 3.0;                           //fprintf( stderr, "1res1 %f\n", s1*s1*s1 + A2*s1*s1 + A1*s1 + A0);
		}
		else
		{
			NumOfRealSolusion = 2;//重解
			double a1 = t_cubert( -q/2 );
			s1 = a1 * 2                 - A2/3.0;     //fprintf( stderr, "22res1 %f\n", s1*s1*s1 + A2*s1*s1 + A1*s1 + A0);
			s2 = a1 * 2 * cos( pi_2_3 ) - A2/3.0;     //fprintf( stderr, "22res2 %f\n", s2*s2*s2 + A2*s2*s2 + A1*s2 + A0);
		}
	}
	else if( D > 0)//実解一個
	{
		NumOfRealSolusion = 1;
		s1 = t_cubert( -q/2 + sqrt(D) ) + 
			 t_cubert( -q/2 - sqrt(D) ) - A2/3.0;     //fprintf( stderr, "11res1 %f\n", s1*s1*s1 + A2*s1*s1 + A1*s1 + A0);
	}
}

//平均は1, 分散widthのガウス分布の乱数ジェネレータ
inline double t_getRandomGaussian(double width = 1)
{
	//Box-Muller transformによりガウス分布に
	double a = ( (double) rand() / (double) RAND_MAX);
	double b = ( (double) rand() / (double) RAND_MAX);
	a = sqrt( -2 * log(a)) * sin(M_PI * 2 * b);//ガウス分布平均0 分散1
	a = a * width + 1.0;
	if( a < 0  ) return 0;
	if( a > 2.0) return 2.0;
	return a;//ガウス分布 平均1.0　分散width
}
inline void t_solveLinearEquationByCG(const int N, double** A, const double* b, double* result)
{
	size_t double_N = sizeof( double ) * N;

	double *cg_r = new double[ N ];
	double *cg_d = new double[ N ];
	double *cg_q = new double[ N ];

	memset( result, 0, double_N );
	int iteration = 0;

	//r = b - Ax (  A is symmetric!!  )
	memcpy( cg_r, b, double_N );
	for( int i = 0; i < N; ++i) 
	for( int j = 0; j < N; ++j) cg_r[ i ] -= A[i][j] * result[ j ];
		
	//d = r
	memcpy( cg_d, cg_r, double_N );

	//deltaNew = r_t * r
	double deltaNew = 0;
	for( int i = 0; i < N; ++i) deltaNew += cg_r[i] * cg_r[i];

	while( iteration < N && deltaNew > 0.000000001)
	{
		//q = Ad
		memset( cg_q, 0, double_N );
		for( int i = 0; i < N; ++i) 
		for( int j = 0; j < N; ++j) cg_q[i] += A[i][j] * cg_d[j];

		//alpha = deltaNew / (d_t * q )
		double alpha = 0;
		for( int i = 0; i < N ; ++i ) alpha += cg_d[i] * cg_q[i];
		alpha = deltaNew / alpha;

		// x = x + alpha * d
		for(int i = 0; i < N; ++i) result[i] += alpha * cg_d[i];

		if( iteration % 30 == 0 ){//r = b - Ax
			memcpy( cg_r, b, double_N );
			for( int i = 0; i < N; ++i) 
			for( int j = 0; j < N; ++j) cg_r[i] -= A[i][j] * result[j];

		}else{//r = r - alpha * q
			for( int i = 0; i < N; ++i) cg_r[i] -= alpha * cg_q[i];		
		}

		//deltaOld = deltaNew
		double deltaOld = deltaNew;

		//deltaNew = r*r
		deltaNew = 0;
		for( int i = 0; i < N; ++i) deltaNew += cg_r[i] * cg_r[i];

		//beta = deltaNew / deltaOld
		double beta = deltaNew / deltaOld;
	
		//d = r + beta + d
		for( int i = 0; i < N; ++i) cg_d[i] = cg_r[i] + beta * cg_d[i];
		++iteration;
	}

	delete[] cg_r;
	delete[] cg_d;
	delete[] cg_q;
}


inline bool t_solveLinearEquationByCG(const int N, 
	                                  const double* Ax,
							          const int*    Ap, 
									  const int*    Ai, 
									  const double* b, 
											      double* result, 
												  double  threshold)
{
	size_t double_N = sizeof( double ) * N;	
	double *m_cg_r = new double[ N ];
	double *m_cg_d = new double[ N ];
	double *m_cg_q = new double[ N ];
	const int    *m_Ap_c = Ap;
	const int    *m_Ai_c = Ai;
	const double *m_Ax_c = Ax;

	int iteration = 0;

	//r = b - Ax (  A is symmetric!!  )
	memcpy( m_cg_r, b, double_N );
	memset( result, 0, double_N );

	for( int i = 0        ; i < N          ; ++i)
	for( int j = m_Ap_c[i]; j < m_Ap_c[i+1]; ++j)
			m_cg_r[ i ] -= m_Ax_c[ j ] * result[ m_Ai_c[j] ];
		
	//d = r
	memcpy( m_cg_d, m_cg_r, double_N );

	//deltaNew = r_t * r
	double deltaNew = 0;
	for( int i = 0; i < N; ++i) deltaNew += m_cg_r[i] * m_cg_r[i];


	//deltaZero = deltaNew
	double deltaZero = deltaNew;

	while( ( iteration < N && deltaNew > threshold) )
	{
		//q = Ad
		memset( m_cg_q, 0, double_N );
		for( int i = 0        ; i < N          ; ++i)
		for( int j = m_Ap_c[i]; j < m_Ap_c[i+1]; ++j)
			m_cg_q[ i ] += m_Ax_c[ j ] * m_cg_d[ m_Ai_c[j] ];

	
		//alpha = deltaNew / (d_t * q )
		double alpha = 0;
		for( int i = 0; i < N ; ++i ) alpha += m_cg_d[i] * m_cg_q[i];
		alpha = deltaNew / alpha;


		// x = x + alpha * d
		for(int i = 0; i < N; ++i) result[i] += alpha * m_cg_d[i];

		if( iteration % 20 == 0 ){
			//r = b - Ax
			memcpy( m_cg_r, b, double_N );
			for( int i = 0        ; i < N          ; ++i)
			for( int j = m_Ap_c[i]; j < m_Ap_c[i+1]; ++j)
					m_cg_r[ i ] -= m_Ax_c[ j ] * result[ m_Ai_c[j] ];

		}else{	
			//r = r - alpha * q
			for( int i = 0; i < N; ++i) m_cg_r[i] -= alpha * m_cg_q[i];		
		}

		//deltaOld = deltaNew
		double deltaOld = deltaNew;

		//deltaNew = r*r
		deltaNew = 0;
		for( int i = 0; i < N; ++i) deltaNew += m_cg_r[i] * m_cg_r[i];

		//beta = deltaNew / deltaOld
		double beta = deltaNew / deltaOld;
	
		//d = r + beta + d
		for( int i = 0; i < N; ++i) m_cg_d[i] = m_cg_r[i] + beta * m_cg_d[i];
		
		++iteration;
	}

	delete[] m_cg_r ;
	delete[] m_cg_d ;
	delete[] m_cg_q ;
	return true;
}
template <class T> void t_copy(vector<T>& src,list<T>& targ){
	targ.clear() ;
	for( vector<T>::iterator vit = src.begin() ; vit != src.end() ; vit++ )
		targ.push_back(*vit) ;
}
template <class T> void t_copy(list<T>& src,vector<T>& targ){
	targ.clear() ;
	targ.resize(src.size()) ;
	int idx=0 ;
	for( list<T>::iterator lit = src.begin() ; lit != src.end() ; lit++,idx++ )
		targ[idx] = *lit ;
}
inline void t_crop( double &t, const double &minVal, const double &maxVal){
	if( t < minVal ) t = minVal;
	if( t > maxVal ) t = maxVal;
}
inline void t_cropI( int &t, const int &minVal, const  int &maxVal)
{
	if( t < minVal ) t = minVal;
	if( t > maxVal ) t = maxVal;
}
//double *M, double *Uは TMatrix9と同じ並び(縦方向の行列)
inline void t_Jacobi_usePreU(double *M9, double *U9, double threshold =0.0000001, int maxIteration = 100)
{
	double d[9];
	//M = Utmp * M * U
	//d = M*U
	d[0] = M9[0] * U9[0]  +  M9[3] * U9[1]  +  M9[6] * U9[2];
	d[1] = M9[1] * U9[0]  +  M9[4] * U9[1]  +  M9[7] * U9[2];
	d[2] = M9[2] * U9[0]  +  M9[5] * U9[1]  +  M9[8] * U9[2];
	d[3] = M9[0] * U9[3]  +  M9[3] * U9[4]  +  M9[6] * U9[5];
	d[4] = M9[1] * U9[3]  +  M9[4] * U9[4]  +  M9[7] * U9[5];
	d[5] = M9[2] * U9[3]  +  M9[5] * U9[4]  +  M9[8] * U9[5];
	d[6] = M9[0] * U9[6]  +  M9[3] * U9[7]  +  M9[6] * U9[8];
	d[7] = M9[1] * U9[6]  +  M9[4] * U9[7]  +  M9[7] * U9[8];
	d[8] = M9[2] * U9[6]  +  M9[5] * U9[7]  +  M9[8] * U9[8];

	M9[0] = U9[0] * d[0]  +  U9[1] * d[1]  +  U9[2] * d[2];
	M9[1] = U9[3] * d[0]  +  U9[4] * d[1]  +  U9[5] * d[2];
	M9[2] = U9[6] * d[0]  +  U9[7] * d[1]  +  U9[8] * d[2];

	M9[3] = U9[0] * d[3]  +  U9[1] * d[4]  +  U9[2] * d[5];
	M9[4] = U9[3] * d[3]  +  U9[4] * d[4]  +  U9[5] * d[5];
	M9[5] = U9[6] * d[3]  +  U9[7] * d[4]  +  U9[8] * d[5];

	M9[6] = U9[0] * d[6]  +  U9[1] * d[7]  +  U9[2] * d[8];
	M9[7] = U9[3] * d[6]  +  U9[4] * d[7]  +  U9[5] * d[8];
	M9[8] = U9[6] * d[6]  +  U9[7] * d[7]  +  U9[8] * d[8];


	for(int count = 0; count < maxIteration; ++count)
	{
		//非対角成分をまわる
		if     ( fabs( M9[3]) >= fabs( M9[6]) && fabs(M9[3]) >= fabs( M9[7]) ) 
		{ 
			if(fabs( M9[3]) <= threshold) break;
			double theta =   (fabs( M9[ 0 ] -   M9[ 4 ]) >0.000001) ? 0.5  * atan(- 2.0 * M9[ 3 ] / ( M9[ 0 ] - M9[ 4 ]) ) : 0.25 * M_PI ;

			double u00 = cos(theta ), u01 = sin(theta );
			double u10 = -u01;
			double u11 = u00;

			//U *= tmpU
			d[0] = U9[0] * u00 + U9[3] * u10;    d[3] = U9[0] * u01 + U9[3] * u11;    d[6] = U9[6] ;
			d[1] = U9[1] * u00 + U9[4] * u10;    d[4] = U9[1] * u01 + U9[4] * u11;    d[7] = U9[7] ;
			d[2] = U9[2] * u00 + U9[5] * u10;    d[5] = U9[2] * u01 + U9[5] * u11;    d[8] = U9[8] ;
			memcpy( U9, d, sizeof(double)*9 ) ;

			//M = Utmp.getTransposed() * M * Utmp
			//d = M * Utmp
			d[0] = M9[0] * u00 + M9[3] * u10;   d[3] = M9[0] * u01 + M9[3] * u11;    d[6] = M9[6] ;
			d[1] = M9[1] * u00 + M9[4] * u10;   d[4] = M9[1] * u01 + M9[4] * u11;    d[7] = M9[7] ;
			d[2] = M9[2] * u00 + M9[5] * u10;   d[5] = M9[2] * u01 + M9[5] * u11;    d[8] = M9[8] ;

			//M = Utmp.trans * d
			M9[0] = u00 * d[0]  +  u10 * d[1]  ;    M9[3] = u00 * d[3]  +  u10 * d[4];   M9[6] = u00 * d[6]  +  u10 * d[7];
			M9[1] = u01 * d[0]  +  u11 * d[1]  ;    M9[4] = u01 * d[3]  +  u11 * d[4];   M9[7] = u01 * d[6]  +  u11 * d[7];
			M9[2] = d[2];                           M9[5] = d[5];                        M9[8] = d[8];
		}
		else if( fabs( M9[6]) >= fabs( M9[7]) )                                        
		{ 
			if(fabs( M9[6]) <= threshold) break; 
			double theta =   (fabs( M9[ 0 ] - M9[ 8 ]) >0.000001) ?  0.5  * atan(- 2.0 * M9[ 6 ] / ( M9[ 0 ] - M9[ 8 ]) ) : 0.25 * M_PI ;
			double u00 = cos( theta ), u02 = sin( theta );
			double u20 = -u02;
			double u22 =  u00;
			//U *= Utmp;
			d[0] = U9[0] * u00 + U9[6] * u20;    d[3] = U9[3] ;    d[6] = U9[0] * u02 + U9[6] * u22;
			d[1] = U9[1] * u00 + U9[7] * u20;    d[4] = U9[4] ;    d[7] = U9[1] * u02 + U9[7] * u22;
			d[2] = U9[2] * u00 + U9[8] * u20;    d[5] = U9[5] ;    d[8] = U9[2] * u02 + U9[8] * u22;
			memcpy( U9, d, sizeof(double)*9 ) ;
			
			//M = Utmp.getTransposed() * M * Utmp;
			//d = M * Utmp
			d[0] = M9[0] * u00 + M9[6] * u20;    d[3] = M9[3];     d[6] = M9[0] * u02 + M9[6] * u22;
			d[1] = M9[1] * u00 + M9[7] * u20;    d[4] = M9[4];     d[7] = M9[1] * u02 + M9[7] * u22;
			d[2] = M9[2] * u00 + M9[8] * u20;    d[5] = M9[5];     d[8] = M9[2] * u02 + M9[8] * u22;

			//M = Utmp.trans * d
			M9[0] = u00 * d[0]  +  u20 * d[2];     M9[3] = u00 * d[3]  +  u20 * d[5];    M9[6] = u00 * d[6]  +  u20 * d[8];
			M9[1] =       d[1];                    M9[4] =       d[4];                   M9[7] =       d[7]; 
 			M9[2] = u02 * d[0]  +  u22 * d[2];     M9[5] = u02 * d[3]  +  u22 * d[5];    M9[8] = u02 * d[6]  +  u22 * d[8];
		}
		else
		{ 
			if(fabs( M9[7]) <= threshold) break; 
			double theta =   (fabs( M9[ 4 ] - M9[ 8 ]) >0.000001) ? 0.5  * atan(- 2.0 * M9[ 7 ] / ( M9[ 4 ] - M9[ 8 ]) ) : 0.25 * M_PI ;

			double u11 = cos( theta ), u12 = sin( theta );
			double u22 = u11;
			double u21 = -u12;

			//U *= Utmp;
			d[0] = U9[0] ;    d[3] = U9[3] * u11 + U9[6] * u21;    d[6] = U9[3] * u12 + U9[6] * u22;
			d[1] = U9[1] ;    d[4] = U9[4] * u11 + U9[7] * u21;    d[7] = U9[4] * u12 + U9[7] * u22;
			d[2] = U9[2] ;    d[5] = U9[5] * u11 + U9[8] * u21;    d[8] = U9[5] * u12 + U9[8] * u22;
			memcpy( U9, d,sizeof(double)*9 ) ;

			//d = M*Utmp
			d[0] = M9[0] ;    d[3] = M9[3] * u11 + M9[6] * u21;    d[6] = M9[3] * u12 + M9[6] * u22;
			d[1] = M9[1] ;    d[4] = M9[4] * u11 + M9[7] * u21;    d[7] = M9[4] * u12 + M9[7] * u22;
			d[2] = M9[2] ;    d[5] = M9[5] * u11 + M9[8] * u21;    d[8] = M9[5] * u12 + M9[8] * u22;

			//M = Utmp.trans * d
			M9[0] = d[0];                         M9[3] = d[3];                         M9[6] = d[6];
			M9[1] = u11 * d[1]  +  u21 * d[2];    M9[4] = u11 * d[4]  +  u21 * d[5];    M9[7] = u11 * d[7]  +  u21 * d[8];
			M9[2] = u12 * d[1]  +  u22 * d[2];    M9[5] = u12 * d[4]  +  u22 * d[5];    M9[8] = u12 * d[7]  +  u22 * d[8];			
		}
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//Vector matrix operations///////////////////////////////////////////////////////////////////////////
//trgt = M * trgt
inline void t_MatMultVec(TVector3 &trgt, const TMatrix16 &M){
	double x = trgt.data[0] * M.data[0]   + trgt.data[1] * M.data[4]   + trgt.data[2] * M.data[ 8]   + M.data[12]; 
	double y = trgt.data[0] * M.data[1]   + trgt.data[1] * M.data[5]   + trgt.data[2] * M.data[ 9]   + M.data[13]; 
	double z = trgt.data[0] * M.data[2]   + trgt.data[1] * M.data[6]   + trgt.data[2] * M.data[10]   + M.data[14]; 
	trgt.Set(x,y,z);
}
inline void t_MatMultVec(TVector3 &trgt, const TMatrix9 &M){
	double x = trgt.data[0] * M.data[0]   + trgt.data[1] * M.data[3]   + trgt.data[2] * M.data[ 6];
	double y = trgt.data[0] * M.data[1]   + trgt.data[1] * M.data[4]   + trgt.data[2] * M.data[ 7];
	double z = trgt.data[0] * M.data[2]   + trgt.data[1] * M.data[5]   + trgt.data[2] * M.data[ 8];
	trgt.Set(x,y,z);
}
inline void t_MatMultVecOnlyRotation(TVector3 &trgt, const TMatrix16 &M){
	double x = trgt.data[0] * M.data[0]   + trgt.data[1] * M.data[4]   + trgt.data[2] * M.data[ 8]; 
	double y = trgt.data[0] * M.data[1]   + trgt.data[1] * M.data[5]   + trgt.data[2] * M.data[ 9]; 
	double z = trgt.data[0] * M.data[2]   + trgt.data[1] * M.data[6]   + trgt.data[2] * M.data[10]; 
	trgt.Set(x,y,z);
}
//trgt = M * V1
inline void t_MatMultVec(const TMatrix16 &M, const TVector3 &V1, TVector3 &trgt){
	trgt.data[0] = V1.data[0] * M.data[0]   + V1.data[1] * M.data[4]   + V1.data[2] * M.data[8]  + M.data[12]; 
	trgt.data[1] = V1.data[0] * M.data[1]   + V1.data[1] * M.data[5]   + V1.data[2] * M.data[9]  + M.data[13]; 
	trgt.data[2] = V1.data[0] * M.data[2]   + V1.data[1] * M.data[6]   + V1.data[2] * M.data[10] + M.data[14]; 
}
inline void t_MatMultVec(const TMatrix9  &M, const TVector3 &V1, TVector3 &trgt){
	trgt.data[0] = V1.data[0] * M.data[0]   + V1.data[1] * M.data[3]   + V1.data[2] * M.data[6]; 
	trgt.data[1] = V1.data[0] * M.data[1]   + V1.data[1] * M.data[4]   + V1.data[2] * M.data[7]; 
	trgt.data[2] = V1.data[0] * M.data[2]   + V1.data[1] * M.data[5]   + V1.data[2] * M.data[8]; 
}
//trgt = V1 - M*V2
inline void t_MatV1minusMatMultV2(const TVector3 &V1, const TMatrix9 &M, const TVector3 &V2, TVector3 &trgt)
{
	trgt.data[0] = V1.data[0] - V2.data[0] * M.data[0] - V2.data[1] * M.data[3] - V2.data[2] * M.data[6]; 
	trgt.data[1] = V1.data[1] - V2.data[0] * M.data[1] - V2.data[1] * M.data[4] - V2.data[2] * M.data[7]; 
	trgt.data[2] = V1.data[2] - V2.data[0] * M.data[2] - V2.data[1] * M.data[5] - V2.data[2] * M.data[8]; 
}
//trgt += M*V1 + V2
inline void t_addMatV1plusV2(const TMatrix9 &M, const TVector3 &V1, const TVector3 &V2, TVector3 &trgt)
{
	trgt.data[0] += V1.data[0] * M.data[0]   + V1.data[1] * M.data[3]   + V1.data[2] * M.data[6] + V2.data[0]; 
	trgt.data[1] += V1.data[0] * M.data[1]   + V1.data[1] * M.data[4]   + V1.data[2] * M.data[7] + V2.data[1]; 
	trgt.data[2] += V1.data[0] * M.data[2]   + V1.data[1] * M.data[5]   + V1.data[2] * M.data[8] + V2.data[2];
}
//trgt = M*V1 + V2
inline void t_setMatV1plusV2(const TMatrix9 &M, const TVector3 &V1, const TVector3 &V2, TVector3 &trgt)
{
	trgt.data[0] = V1.data[0] * M.data[0]   + V1.data[1] * M.data[3]   + V1.data[2] * M.data[6] + V2.data[0]; 
	trgt.data[1] = V1.data[0] * M.data[1]   + V1.data[1] * M.data[4]   + V1.data[2] * M.data[7] + V2.data[1]; 
	trgt.data[2] = V1.data[0] * M.data[2]   + V1.data[1] * M.data[5]   + V1.data[2] * M.data[8] + V2.data[2];
}
inline void t_MatMultVec(TVector3 &tgt, const TMatrix4 &m, const TVector3 &v){
	tgt.data[0] = m.data[0] * v.data[0] + m.data[2] * v.data[1];
	tgt.data[1] = m.data[1] * v.data[0] + m.data[3] * v.data[1];
}
inline void t_MatMultVec2D( const TMatrix16& M, TVector3 &v)
{
	double  x = M.data[0] * v.data[0] + M.data[4] * v.data[1] + M.data[12];
	double  y = M.data[1] * v.data[0] + M.data[5] * v.data[1] + M.data[13];
	v.Set(x,y,0);
}
inline void t_MatMultVec2D(TVector3 &trgt, const TMatrix16& M, const TVector3 &v)
{
	trgt.data[0] = M.data[0] * v.data[0] + M.data[4] * v.data[1] + M.data[12];
	trgt.data[1] = M.data[1] * v.data[0] + M.data[5] * v.data[1] + M.data[13];
}
//trgt =  (v1-v2)^(v3-v4)
inline void t_V1subtV2_cros_V3subtV4(
	const TVector3 &v1, const TVector3 &v2, 
	const TVector3 &v3, const TVector3 &v4,
	      TVector3 &trgt)
{
	trgt.data[0] = (v1.data[1] - v2.data[1]) * (v3.data[2] - v4.data[2])   -   (v1.data[2] - v2.data[2]) * (v3.data[1] - v4.data[1]);
	trgt.data[1] = (v1.data[2] - v2.data[2]) * (v3.data[0] - v4.data[0])   -   (v1.data[0] - v2.data[0]) * (v3.data[2] - v4.data[2]);
	trgt.data[2] = (v1.data[0] - v2.data[0]) * (v3.data[1] - v4.data[1])   -   (v1.data[1] - v2.data[1]) * (v3.data[0] - v4.data[0]);	
}

//trgt =  (v1-v2)^(v3-v4)
inline void t_add_V1subtV2_cros_V3subtV4( const TVector3 &v1, const TVector3 &v2, 
										  const TVector3 &v3, const TVector3 &v4,
												TVector3 &trgt)
{
	trgt.data[0] += (v1.data[1] - v2.data[1]) * (v3.data[2] - v4.data[2])   -   (v1.data[2] - v2.data[2]) * (v3.data[1] - v4.data[1]);
	trgt.data[1] += (v1.data[2] - v2.data[2]) * (v3.data[0] - v4.data[0])   -   (v1.data[0] - v2.data[0]) * (v3.data[2] - v4.data[2]);
	trgt.data[2] += (v1.data[0] - v2.data[0]) * (v3.data[1] - v4.data[1])   -   (v1.data[1] - v2.data[1]) * (v3.data[0] - v4.data[0]);	
}
// ((v1-v2)^(v3-v4))  * V5
inline double t_V1subtV2_cros_V3subtV4_multV5(const TVector3 &v1, const TVector3 &v2,
											  const TVector3 &v3, const TVector3 &v4, const TVector3 &v5)
{
	return  ( (v1.data[1]-v2.data[1]) * (v3.data[2]-v4.data[2])  -  (v1.data[2]-v2.data[2]) * (v3.data[1]-v4.data[1]) ) * v5.data[0]+
			( (v1.data[2]-v2.data[2]) * (v3.data[0]-v4.data[0])  -  (v1.data[0]-v2.data[0]) * (v3.data[2]-v4.data[2]) ) * v5.data[1]+
			( (v1.data[0]-v2.data[0]) * (v3.data[1]-v4.data[1])  -  (v1.data[1]-v2.data[1]) * (v3.data[0]-v4.data[0]) ) * v5.data[2];	
}
//(V1^V2) * V3
inline double t_V1crosV2_multV3(const TVector3& V1, const TVector3& V2, const TVector3& V3)
{
	return	( V1.data[1] * V2.data[2] - V1.data[2] * V2.data[1] ) * V3.data[0]
		   +( V1.data[2] * V2.data[0] - V1.data[0] * V2.data[2] ) * V3.data[1]
		   +( V1.data[0] * V2.data[1] - V1.data[1] * V2.data[0] ) * V3.data[2];
}
//(V1-V2) * V3
inline double t_V1subtV2_multV3(const TVector3& V1, const TVector3& V2, const TVector3& V3){
	return	(V1.data[0] - V2.data[0]) * V3.data[0] + 
			(V1.data[1] - V2.data[1]) * V3.data[1] + 
			(V1.data[2] - V2.data[2]) * V3.data[2] ;
}
//trgt = M1 * M2
inline void t_MatMultMat(const TMatrix16 &M1, const TMatrix16 &M2, TMatrix16 &trgt )
{
	trgt.data[0]  = M1.data[0]*M2.data[0]  +  M1.data[4]*M2.data[1 ] + M1.data[8 ]*M2.data[2 ] + M1.data[12]*M2.data[3];
	trgt.data[1]  = M1.data[1]*M2.data[0]  +  M1.data[5]*M2.data[1 ] + M1.data[9 ]*M2.data[2 ] + M1.data[13]*M2.data[3]; 
	trgt.data[2]  = M1.data[2]*M2.data[0]  +  M1.data[6]*M2.data[1 ] + M1.data[10]*M2.data[2 ] + M1.data[14]*M2.data[3]; 
	trgt.data[3]  = M1.data[3]*M2.data[0]  +  M1.data[7]*M2.data[1 ] + M1.data[11]*M2.data[2 ] + M1.data[15]*M2.data[3]; 

	trgt.data[4]  = M1.data[0]*M2.data[4]  +  M1.data[4]*M2.data[5 ] + M1.data[8 ]*M2.data[6 ] + M1.data[12]*M2.data[7];
	trgt.data[5]  = M1.data[1]*M2.data[4]  +  M1.data[5]*M2.data[5 ] + M1.data[9 ]*M2.data[6 ] + M1.data[13]*M2.data[7]; 
	trgt.data[6]  = M1.data[2]*M2.data[4]  +  M1.data[6]*M2.data[5 ] + M1.data[10]*M2.data[6 ] + M1.data[14]*M2.data[7]; 
	trgt.data[7]  = M1.data[3]*M2.data[4]  +  M1.data[7]*M2.data[5 ] + M1.data[11]*M2.data[6 ] + M1.data[15]*M2.data[7]; 

	trgt.data[8]  = M1.data[0]*M2.data[8]  +  M1.data[4]*M2.data[9 ] + M1.data[8 ]*M2.data[10] + M1.data[12]*M2.data[11];
	trgt.data[9]  = M1.data[1]*M2.data[8]  +  M1.data[5]*M2.data[9 ] + M1.data[9 ]*M2.data[10] + M1.data[13]*M2.data[11]; 
	trgt.data[10] = M1.data[2]*M2.data[8]  +  M1.data[6]*M2.data[9 ] + M1.data[10]*M2.data[10] + M1.data[14]*M2.data[11]; 
	trgt.data[11] = M1.data[3]*M2.data[8]  +  M1.data[7]*M2.data[9 ] + M1.data[11]*M2.data[10] + M1.data[15]*M2.data[11]; 

	trgt.data[12] = M1.data[0]*M2.data[12] +  M1.data[4]*M2.data[13] + M1.data[8 ]*M2.data[14] + M1.data[12]*M2.data[15];
	trgt.data[13] = M1.data[1]*M2.data[12] +  M1.data[5]*M2.data[13] + M1.data[9 ]*M2.data[14] + M1.data[13]*M2.data[15]; 
	trgt.data[14] = M1.data[2]*M2.data[12] +  M1.data[6]*M2.data[13] + M1.data[10]*M2.data[14] + M1.data[14]*M2.data[15]; 
	trgt.data[15] = M1.data[3]*M2.data[12] +  M1.data[7]*M2.data[13] + M1.data[11]*M2.data[14] + M1.data[15]*M2.data[15]; 
}
//Vector matrix operations///////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
//Geomtric operations////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
//Distances/////////////////////////////////////////////////////////////////////////
inline double t_distance_sq(const TVector3 &v1, const TVector3 &v2){
	return  (v1.data[0] - v2.data[0]) * (v1.data[0] - v2.data[0]) + 
		    (v1.data[1] - v2.data[1]) * (v1.data[1] - v2.data[1]) + 
		    (v1.data[2] - v2.data[2]) * (v1.data[2] - v2.data[2]) ;		         
}
inline double t_distance2D_sq(const TVector3 &x1, const TVector3 &x2){
	return (x1.data[0] - x2.data[0]) * (x1.data[0] - x2.data[0]) + 
		   (x1.data[1] - x2.data[1]) * (x1.data[1] - x2.data[1]);
}
inline double t_distance   (const TVector3 &v1, const TVector3 &v2){ return sqrt( t_distance_sq  (v1, v2) ); }
inline double t_distance2D( const TVector3& x1, const TVector3 &x2){ return sqrt( t_distance2D_sq(x1, x2) );}

inline double t_distance2D_sq( double x1, double y1, double x2, double y2 ){ 
	return (x1-x2) * (x1-x2) + (y1-y2) * (y1-y2);
}
inline double t_distance2D( double x1, double y1, double x2, double y2 ){
	return sqrt( t_distance2D_sq( x1, y1, x2, y2));
}
//rot*v1 and  v2//
inline double t_distance2D_sq(const TMatrix4 &rot, const TVector3 &v1, const TVector3 &v2)
{
	return ( v2.data[0]  -  (rot.data[0] * v1.data[0] + rot.data[2] * v1.data[1]) ) * 
		   ( v2.data[0]  -  (rot.data[0] * v1.data[0] + rot.data[2] * v1.data[1]) ) 
		   + 
		   ( v2.data[1]  -  (rot.data[1] * v1.data[0] + rot.data[3] * v1.data[1]) ) * 
		   ( v2.data[1]  -  (rot.data[1] * v1.data[0] + rot.data[3] * v1.data[1]) ); 
}

inline double t_distance_sq(const TVector2 &v1, const TVector2 &v2){
	return  (v1.data[0] - v2.data[0]) * (v1.data[0] - v2.data[0]) + 
		    (v1.data[1] - v2.data[1]) * (v1.data[1] - v2.data[1]);
}
inline double t_distance(const TVector2 &v1, const TVector2 &v2){
	return  sqrt( t_distance_sq( v1,v2) );
}
/*---------------------------------------------------------------------------------
H  = lineP0 + t * lineDir
PH = lineP0 + t * lineDir - P
PH * lineDir = 0

lineDir * (lineP0 + t * lineDir - P) = 0
t = (P-lineP0) * lineDir / (lineDir*lineDir)

---------------------------------------------------------------------------------*/
inline double t_distPointToLine_sq( const TVector3 &p, const TVector3 &linePos, const TVector3 &lineDir, double &t)
{
	t =  (p.data[0] - linePos.data[0]) * lineDir.data[0] + 
		 (p.data[1] - linePos.data[1]) * lineDir.data[1] +
		 (p.data[2] - linePos.data[2]) * lineDir.data[2];
	t = t / lineDir.Length_Square();

	double x = linePos.data[0]  + t * lineDir.data[0] - p.data[0];
	double y = linePos.data[1]  + t * lineDir.data[1] - p.data[1];
	double z = linePos.data[2]  + t * lineDir.data[2] - p.data[2];
	return x * x + y*y + z * z;
}
inline double t_distPointToLine_sq( const TVector3 &p, const TVector3 &linePos, const TVector3 &lineDir){
	double t = 0; 
	return t_distPointToLine_sq( p, linePos, lineDir, t) ;
}
inline double t_distPointToLine( const TVector3 &p, const TVector3 &linePos, const TVector3 &lineDir, double &t){
	return sqrt( t_distPointToLine_sq( p, linePos, lineDir, t) );
}
inline double t_distPointToLine( const TVector3 &p, const TVector3 &linePos, const TVector3 &lineDir){
	double t = 0;
	return sqrt( t_distPointToLine_sq( p, linePos, lineDir, t) );
}
inline double t_distPointToLineSegment_sq( const TVector3 &p , 
								           const TVector3 &lineP0, 
								           const TVector3 &lineP1)
{
	double t =  (p.data[0] - lineP0.data[0]) * (lineP1.data[0] - lineP0.data[0]) + 
		        (p.data[1] - lineP0.data[1]) * (lineP1.data[1] - lineP0.data[1]) + 
				(p.data[2] - lineP0.data[2]) * (lineP1.data[2] - lineP0.data[2]);
	t /= t_distance_sq( lineP1, lineP0 );

	if( t < 0 ) return t_distance_sq( p, lineP0 );
	if( t > 1 ) return t_distance_sq( p, lineP1 );

	double x = lineP0.data[0]  + t * (lineP1.data[0]-lineP0.data[0]) - p.data[0];
	double y = lineP0.data[1]  + t * (lineP1.data[1]-lineP0.data[1]) - p.data[1];
	double z = lineP0.data[2]  + t * (lineP1.data[2]-lineP0.data[2]) - p.data[2];
	return x*x + y*y + z*z;
}
inline double t_distPointToLineSegment( const TVector3 &p, const TVector3 &lineP0, const TVector3 &lineP1){
	return sqrt( t_distPointToLineSegment_sq( p, lineP0, lineP1 ) );
}
inline double t_distLineToLineSegment_sq( const TVector3 &lineP0, const TVector3 &lineP1, 
								          const TVector3 &segmP0, const TVector3 &segmP1,
										        bool     &isOnSegm /*最近傍線の垂線の足がline segm上に乗っていたらtrue*/)
{
	//2直線の最短距離を結ぶ垂線の足をそれぞれ 
	//x0 = lineP0     + t0 * (lineP1    - lineP0    )
	//x1 = lineSegmP0 + t1 * (lineSegmP1- lineSegmP0)としてt0,t1を計算する
	isOnSegm = false;
	TVector3 d0, d1;
	d0.SetSubtract( lineP1, lineP0 );
	d1.SetSubtract( segmP1, segmP0 );

	double d0d0 = d0*d0; if( d0d0 <0.00000001 ) return 0;
	double d1d1 = d1*d1; if( d1d1 <0.00000001 ) return 0;
	double d0d1 = d0*d1;
	double b0   = d0 * (segmP0-lineP0);
	double b1   = d1 * (segmP0-lineP0);
	double t0,t1;
	t_solve2by2LinearEquation( d0d0, -d0d1, 
		                       d0d1, -d1d1, b0, b1, t0,t1);

	if( 0<=t1 && t1 <= 1) isOnSegm = true;

	return t_distance_sq( lineP0 + t0 * (lineP1-lineP0), 
		                  segmP0 + t1 * (segmP1-segmP0) );
}


// ph = p0 + t * d - p; ph * d = 0//
inline double t_distPointToLineSegment_2D(  const TVector3 &p , 
											const TVector3 &lineP0, 
											const TVector3 &lineP1)
{
	double t =  (p.data[0] - lineP0.data[0]) * (lineP1.data[0] - lineP0.data[0]) + 
		        (p.data[1] - lineP0.data[1]) * (lineP1.data[1] - lineP0.data[1]);
	t /= t_distance2D_sq( lineP1, lineP0 );

	if( t < 0 ) return t_distance2D( p, lineP0 );
	if( t > 1 ) return t_distance2D( p, lineP1 );

	double x = lineP0.data[0]  + t * (lineP1.data[0]-lineP0.data[0]) - p.data[0];
	double y = lineP0.data[1]  + t * (lineP1.data[1]-lineP0.data[1]) - p.data[1];

	return sqrt( x*x + y*y);
}
//////////////////////////////////////////////////////////////////////////////////////////////
//Window check////////////////////////////////////////////////////////////////////////////////
inline bool t_isInWindow2D( const CPoint &p, const CPoint &minXY, const CPoint &maxXY)
{	
	if( minXY.x <= p.x && p.x <= maxXY.x &&
		minXY.y <= p.y && p.y <= maxXY.y ) return true;
	return false;
}
inline bool t_isInWindow2D( const CPoint &p, const CPoint &minXY, int width, int height)
{	
	if( minXY.x <= p.x && p.x <= minXY.x + width &&
		minXY.y <= p.y && p.y <= minXY.y +height) return true;
	return false;
}
inline bool t_isInWindow2D( const TVector3 &piv, double windowR, const TVector3 &p)
{	
	if( p.data[0] < piv.data[0] - windowR) return false;
	if( p.data[1] < piv.data[1] - windowR) return false;
	if( p.data[0] > piv.data[0] + windowR) return false;
	if( p.data[1] > piv.data[1] + windowR) return false;
	return true;
}
inline bool t_isInWindow2D( const TVector3 &p0, const TVector3 &p1, const TVector3 &p, double offset = 0)
{
	double xmin = MIN( p0.data[0], p1.data[0] );
	double ymin = MIN( p0.data[1], p1.data[1] );
	double xmax = MAX( p0.data[0], p1.data[0] );
	double ymax = MAX( p0.data[1], p1.data[1] );
	if( xmin - offset <= p.data[0] && p.data[0] <= xmax + offset &&
		ymin - offset <= p.data[1] && p.data[1] <= ymax + offset ) return true;
	return false;
}
inline bool t_isInWindow3D( const TVector3 &p0, const TVector3 &p1, const TVector3 &p, double offset = 0)
{
	double xmin = MIN( p0.data[0], p1.data[0] ), xmax = MAX( p0.data[0], p1.data[0] );
	double ymin = MIN( p0.data[1], p1.data[1] ), ymax = MAX( p0.data[1], p1.data[1] );
	double zmin = MIN( p0.data[2], p1.data[2] ), zmax = MAX( p0.data[2], p1.data[2] );
	if( xmin - offset <= p.data[0] && p.data[0] <= xmax + offset &&
		ymin - offset <= p.data[1] && p.data[1] <= ymax + offset && 
		zmin - offset <= p.data[2] && p.data[2] <= zmax + offset  ) return true;
	return false;
}
inline bool t_isInWindow3D_MinMax( const TVector3 &minP, const TVector3 &maxP, const TVector3 &p, double offset = 0)
{
	if( minP.data[0] - offset <= p.data[0] && p.data[0] <= maxP.data[0] + offset &&
		minP.data[1] - offset <= p.data[1] && p.data[1] <= maxP.data[1] + offset && 
		minP.data[2] - offset <= p.data[2] && p.data[2] <= maxP.data[2] + offset  ) return true;
	return false;
}
//////////////////////////////////////////////////////////////////////////////////////
//calcAngle and rotation matrix///////////////////////////////////////////////////////
inline double t_cosine(const TVector3 &v0, const TVector3 &v1)
{
	double l = v0.Length() * v1.Length(); 
	if(l < 0.000001) return 0;
	else             return (v0 * v1) / l;
}
//cot = 1 / tan
inline double t_cotangent(const TVector3 &v0, const TVector3 &v1)
{
	double cosine = t_cosine(v0,v1);
	if(cosine <= 0.0001) return 0;      //cot90 == 0
	if(cosine > 0.99999) return 1000.0; //cot90 >> 1000.0

	//cot = cos / sin = cos / sqrt( 1 - cos*cos)   (sineは確実に正 theta < 180)
	double sine =  sqrt(1.0 - cosine * cosine);
	if( sine <= 0.00001) return 1000.0;
	else				 return cosine / sine;
}
/*-------------------------------------------------
//return angle between v1 and v2.
//axis is the axis of rotation from v1 to v2.
//	axisを軸に、v1からv2までの角度を求める
//	軸方向に右ねじなら正が帰ってくる
-------------------------------------------------*/
inline double t_getAngle_normalized_FixedAxis(const TVector3 &v1_normarized, 
											  const TVector3 &v2_normarized, const TVector3 &axis)
{
	double cosT = v1_normarized * v2_normarized;
	if(cosT >  1.0) cosT =  1.0;
	if(cosT < -1.0) cosT = -1.0;
	double angle = acos(cosT);
	if( (v1_normarized ^ v2_normarized) * axis >= 0) return  angle;
	else											 return -angle;
}

inline double t_getAngle_FixedAxis(const TVector3 &v1, const TVector3 &v2, const TVector3 &axis)
{
	double l = v1.Length() * v2.Length();
	if( l == 0 ) return 0;
	double cosT = ( v1 * v2 ) / l;
	if(cosT >  1.0) cosT =  1.0;
	if(cosT < -1.0) cosT = -1.0;
	double angle = acos(cosT); // -1< cosT < 1 に対し0 <= acos( ) <= pi
	if( (v1 ^ v2) * axis >= 0) return  angle;
	else			           return -angle;
}

//angle between v1 and v2
inline double t_getAngle2D(const TVector3 &v1, const TVector3 &v2 )
{
	double l = v1.Length2D() * v2.Length2D();
	if( l == 0 ) return 0;
	double cosT = ( v1 * v2 ) / l;
	if(cosT >  1.0) cosT =  1.0;
	if(cosT < -1.0) cosT = -1.0;
	double angle = acos(cosT); // -1< cosT < 1 に対し0 <= acos( ) <= pi
	if( v1.data[0] * v2.data[1] - v1.data[1] * v2.data[0] >= 0) return  angle;
	else					                                    return -angle;

}
//ベクトルv1をベクトルv2に沿わせるような観点ベクトル・回転軸を求める
inline void t_getRotMat_V1ToV2_normalized( const TVector3 &v1, const TVector3 &v2, TMatrix16 &rotationMat)
{
	TVector3 rotAxis;

	rotAxis.SetCrossProd(v1, v2);
	double l = rotAxis.Length();
	if( l < 0.000001 ) { rotationMat.SetIdentity(); return; }
	
	rotAxis/= l ;
	double cosT = v1 * v2; 
	if(cosT > 1.0) cosT = 1.0; if(cosT < -1.0) cosT = -1.0;
	rotationMat.RotateAlongArbitraryAxis(rotAxis,acos(cosT));
}
inline void t_getRotMat_V1ToV2    ( TVector3 &v1, TVector3 &v2, TMatrix16 &rotationMat)
{
	v1.Normalize_Self(); 
	v2.Normalize_Self();
	t_getRotMat_V1ToV2_normalized( v1, v2, rotationMat);
}
inline bool t_getRotMat_YToV1( TVector3 &v1, TVector3 &axis, TMatrix16 &rotationMat)
{
	if( !v1.Normalize_Self() ) { rotationMat.SetIdentity(); return false;}

	if( v1.data[2] * v1.data[2] + v1.data[0] * v1.data[0] < 0.001 ) return false;
	axis.Set( v1.data[2], 0, - v1.data[0]);
	rotationMat.RotateAlongArbitraryAxis( axis, acos(v1.data[1]) ); //double angle = acos ( (0,1,0) * v1 )
	return true;
}
inline void t_getRotMat_ZToV1( TVector3 &v1, TVector3 &axis, TMatrix16 &rotationMat)
{
	if( !v1.Normalize_Self() ) { rotationMat.SetIdentity(); return; }
	axis.Set(  -v1.data[1], v1.data[0], 0);
	rotationMat.RotateAlongArbitraryAxis( axis, acos(v1.data[2]) );//double angle = acos ( (0,0,1) * v1 )
}

inline void t_getRotMat_ZToV1_normalized ( const TVector3 &v1_normalized, TMatrix16 &rotationMat)
{
	TVector3 axis; 
	axis.Set(  -v1_normalized.data[1], v1_normalized.data[0], 0);
	rotationMat.RotateAlongArbitraryAxis( axis, acos(v1_normalized.data[2]) ); //double angle = acos ( (0,0,1) * v1 )
}

//local X Y (Z) directionからなる局所座標系に変換する回転行列を生成
inline bool t_getRotMat_LocalXYdir( TVector3 &localX, TVector3 &localY, TMatrix9  &rotationMat)
{
	if( !localX.Normalize_Self() ) return false;
	localY.Add_CoefMultVec( -(localX*localY), localX );
	if( !localY.Normalize_Self() ) return false;

	//localZ = localX ^ localY;
	double localZx = localX.data[1] * localY.data[2] - localX.data[2] * localY.data[1];
	double localZy = localX.data[2] * localY.data[0] - localX.data[0] * localY.data[2];
	double localZz = localX.data[0] * localY.data[1] - localX.data[1] * localY.data[0];

	rotationMat.Set( localX.data[0], localY.data[0], localZx, 
		             localX.data[1], localY.data[1], localZy, 
					 localX.data[2], localY.data[2], localZz);
	return true;
}
////////////////////////////////////////////////////////////////////////////////////////
//BARYCENTRIC COORDINATE////////////////////////////////////////////////////////////////
/* ---------------------------------------------------------------------
面に平行なdirectioにのみ利用可能
 * returns a new Vector3d(alpha, beta, gamma) which satisfies:
 * 
 * T = alpha * p0 + beta * p1 + gamma * p2
 * alpha + beta + gamma = 0
	At	p0.x, p0.y, p0.z, 1,
		p1.x, p1.y, p1.z, 1,
		p2.x, p2.y, p2.z, 1

	A 	p0.x, p1.x, p2.x,
		p0.y, p1.y, p2.y,
		p0.z, p1.z, p2.z,
		1,    1,    1
Pos = a1 * ( p1 - p0) + a2 * (p2 - p0)とおいた

Pos = a0 * p0 + a1 * p1 + a2 * p2
 a0 + a1 + a2 = 0となる
 ---------------------------------------------------------------------*/
inline bool t_calcBaryCoordOfDirectionOnPoly(const TVector3 &T, const TVector3 &p0, const TVector3 &p1, const TVector3 &p2, TVector3 &trgt)
{
	TVector3 d1, d2;
	TMatrix4  AtA4;
	d1.SetSubtract( p1,p0 );
	d2.SetSubtract( p2,p0 );

	AtA4.Set( d1 * d1, d1 * d2, 
		      d2 * d1, d2 * d2 );
	if(!AtA4.getInvertSelf() ){ fprintf( stderr, "error!!!!\n" ); trgt.Set( 0,0,0 ); return false; }

	double b1 = d1 * T,
		   b2 = d2 * T;
	trgt.data[1] = AtA4.data[0] * b1 + AtA4.data[2] * b2; 
	trgt.data[2] = AtA4.data[1] * b1 + AtA4.data[3] * b2;
	trgt.data[0] = - trgt.data[1] - trgt.data[2]; 
	return true;
}
/*---------------------------------------------------------------------
三角形と同一平面に乗る点のみに利用可能 (基本的にはa1 a2しか動かしていない)
 * pos = alpha * p0 + beta * p1 + gamma * p2
 * alpha + beta + gamma = 1

 Pos = p0 + a1 * (p1 -p0) + a2 * (p2 - p0);とおくと
 Pos = a0 * p0 + a1 * p1 + a2 * p2
 a0 + a1 + a2 = 1となる

 0 <= a0, a1, a2 <= 1　ならばPosは平面内にある 
--------------------------------------------------------------------- */
inline bool t_calcBaryCoordOfPositionOnPoly(const TVector3 &pos, const TVector3 &p0, const TVector3 &p1, const TVector3 &p2, TVector3 &trgt)
{
	TMatrix9  AtA;
	TVector3 Atb;
	AtA.Set( p0 * p0 + 1,  p0 * p1 + 1,  p0 * p2 + 1,
		     p1 * p0 + 1,  p1 * p1 + 1,  p1 * p2 + 1,
		     p2 * p0 + 1,  p2 * p1 + 1,  p2 * p2 + 1);
	Atb.Set( p0 * pos + 1, p1 * pos + 1, p2 * pos + 1 );

	if( !AtA.getInvertSelf() ){trgt.Set(0,0,0); fprintf( stderr, "error1111222\n" ); return false;}
	
	t_MatMultVec( AtA, Atb, trgt );
	return true;
}
// pos = alpha * p0 + beta * p1 + gamma * p2
inline void t_calcBaryCoord(const TVector3 &pos, const TVector3 &p0, const TVector3 &p1, const TVector3 &p2, TVector3 &trgt)
{
	TMatrix9  A;
	TVector3 b;
	A.Set(   p0.data[0],  p1.data[0],  p2.data[0],
		     p0.data[1],  p1.data[1],  p2.data[1],
		     p0.data[2],  p1.data[2],  p2.data[2]);
	b.Set( pos);

	if( !A.getInvertSelf() ){trgt.Set(0,0,0); fprintf( stderr, "error1111\n" ); return;}
	t_MatMultVec( A, b, trgt );
}
/* tetra頂点の重心座標系を計算
    pos = x0  +  a1 * v1  +  a2 * v2  +  a3 * v3;

	pos = a0 x0 + a1 x1 + a2 x2 + a3 x3
	      a0 + a1 + a2 + a3 = 1

    a1
    a2  = A(-1) ( pos - x0)
    a3
	を計算してから a0 = 1- a1 - a2 - a3

	0 <= a0,a1,a2,a3 <= 1ならばこの点はtetraの内部に含まれる
*/
inline bool t_calcBaryCoordTetra( const TVector3 &pos, 
								  const TVector3 &x0, const TVector3 &x1, 
								  const TVector3 &x2, const TVector3 &x3, double trgt[4])
{
	TMatrix9 M( x1.data[0] - x0.data[0], x2.data[0] - x0.data[0], x3.data[0] - x0.data[0], 
				 x1.data[1] - x0.data[1], x2.data[1] - x0.data[1], x3.data[1] - x0.data[1], 
				 x1.data[2] - x0.data[2], x2.data[2] - x0.data[2], x3.data[2] - x0.data[2]);

	if( !M.getInvertSelf() ) return false;
	trgt[1] = M.data[0] * (pos.data[0]-x0.data[0])  +  M.data[3] * (pos.data[1]-x0.data[1])  +  M.data[6] * (pos.data[2]-x0.data[2]);
	trgt[2] = M.data[1] * (pos.data[0]-x0.data[0])  +  M.data[4] * (pos.data[1]-x0.data[1])  +  M.data[7] * (pos.data[2]-x0.data[2]);
	trgt[3] = M.data[2] * (pos.data[0]-x0.data[0])  +  M.data[5] * (pos.data[1]-x0.data[1])  +  M.data[8] * (pos.data[2]-x0.data[2]);
	trgt[0] = 1 - (trgt[1] + trgt[2] + trgt[3]);
	return true;
}
/*  dir  = a1 * v1  +  a2 * v2  +  a3 * v3;とおく
	pos = a0 x0 + a1 x1 + a2 x2 + a3 x3
	  0 = a0 + a1 + a2 + a3 
	tetra頂点の平行移動に依存しない方向vectorに利用する*/
inline bool t_calcBaryCoordOfDirectionOnTetra( const TVector3 &pos, 
											    const TVector3 &x0, const TVector3 &x1, 
											    const TVector3 &x2, const TVector3 &x3, double trgt[4])
{
	TMatrix9 M( x1.data[0] - x0.data[0], x2.data[0] - x0.data[0], x3.data[0] - x0.data[0], 
				 x1.data[1] - x0.data[1], x2.data[1] - x0.data[1], x3.data[1] - x0.data[1], 
				 x1.data[2] - x0.data[2], x2.data[2] - x0.data[2], x3.data[2] - x0.data[2]);

	if( !M.getInvertSelf() ) return false;
	trgt[1] = M.data[0] * pos.data[0] +  M.data[3] * pos.data[1] +  M.data[6] * pos.data[2] ;
	trgt[2] = M.data[1] * pos.data[0] +  M.data[4] * pos.data[1] +  M.data[7] * pos.data[2] ;
	trgt[3] = M.data[2] * pos.data[0] +  M.data[5] * pos.data[1] +  M.data[8] * pos.data[2] ;
	trgt[0] = -(trgt[1] + trgt[2] + trgt[3]);
	return true;
}
//Barycentric coordinate////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////





inline void t_internalDivision( double val0, const TVector3 &p0, 
								double val1, const TVector3 &p1, 
								double val,        TVector3 &pos ){
	double u = (val - val0) / (val1 - val0);
	pos.SetAdditionWithCoef( 1-u, p0, u, p1 );
}
inline double t_calcTriangleArea( const TVector3 &v0, const TVector3 &v1)
{
	//s = 0.5 * sqrt(v0*v0*v1*v1 - v0v1*v0v1)
	double v0v1 = v0*v1;
	return 0.5 * sqrt( v0.Length_Square() * v1.Length_Square() - v0v1*v0v1  );
}
inline void t_getRandDirByGaussianRandom( TVector3 &trgt )
{
	trgt.data[0] = t_getRandomGaussian( 1 )-1; 
	trgt.data[1] = t_getRandomGaussian( 1 )-1; 
	trgt.data[2] = t_getRandomGaussian( 1 )-1; 
	trgt.Normalize_Self();
}


///////////////////////////////////////////////////////////////////////////////////////
////////////////SHAPE MATCHING/////////////////////////////////////////////////////////
/*
Mをjacobi法により対角化する
対角が後の行列がMに入り,Uには対角化するための行列が入る
!!Mは対称である必要がある!!
U-1 * M * U                  = diag( r0, r1, r2 )
U   * diag(r0, r1, r2) * U-1 = A
http://case.f7.ems.okayama-u.ac.jp/statedu/eigen/node5.html参照
*/
inline void t_Jacobi(TMatrix9 &M, TMatrix9& U, double threshold=0.00001 , int maxIteration = 100)
{
	U.SetIdentity();
	int dim = 3;

	TMatrix9 Utmp,tmp;
	for(int count = 0; count < maxIteration; ++count)
	{
		int pivI, pivJ;
		double maxVal = -1;
		//非対角成分をまわる
		if     ( fabs( M.data[3]) >= fabs( M.data[6]) && fabs(M.data[3]) >= fabs( M.data[7]) ) { pivI = 0; pivJ = 1; maxVal = fabs( M.data[3]); }
		else if( fabs( M.data[6]) >= fabs( M.data[7]) )                                        { pivI = 0; pivJ = 2; maxVal = fabs( M.data[6]); }
		else                                                                                   { pivI = 1; pivJ = 2; maxVal = fabs( M.data[7]); }

		if(maxVal <= threshold) break;

		double theta = M_PI / 4.0;
		double aii   = M.data[ pivI + 3*pivI ];
		double aij   = M.data[ pivI + 3*pivJ ];
		double ajj   = M.data[ pivJ + 3*pivJ ];

		if( fabs(aii - ajj) >0.0001)
			theta = 0.5 * atan(- 2.0 * aij / ( aii - ajj) );
		
		Utmp.SetIdentity();
		Utmp.data[pivI + 3*pivI ] =  cos(theta);
		Utmp.data[pivJ + 3*pivJ ] =  cos(theta);
		Utmp.data[pivI + 3*pivJ ] =  sin(theta);
		Utmp.data[pivJ + 3*pivI ] = -sin(theta);

		U *= Utmp;

		//M = Utmp.getTransposed() * M * Utmp;
		tmp.setMultM1M2(M,Utmp);
		Utmp.transpose_Self();
		M.setMultM1M2(Utmp,tmp);
	}
}


/*
通常は 対象行列 Aをもらって 
U-1 * A * U = diag( r0, r1, r2 )
というUとdiagを返す
*/
#if 0
void TI_Jacobi_usePreU(TMatrix9 &M, TMatrix9& U, double threshold , int maxIteration)
{
	TMatrix9 Utmp,tmp;
	
	Utmp.setTransposed( U );
	//M = Utmp * M * U
	tmp.setMultM1M2( M, U);
	M.setMultM1M2(Utmp,tmp);

	for(int count = 0; count < maxIteration; ++count)
	{

		int pivI, pivJ;
		double maxVal = -1;
		//非対角成分をまわる
		if     ( fabs( M.data[3]) >= fabs( M.data[6]) && fabs(M.data[3]) >= fabs( M.data[7]) ) { pivI = 0; pivJ = 1; maxVal = fabs( M.data[3]); }
		else if( fabs( M.data[6]) >= fabs( M.data[7]) )                                        { pivI = 0; pivJ = 2; maxVal = fabs( M.data[6]); }
		else                                                                                   { pivI = 1; pivJ = 2; maxVal = fabs( M.data[7]); }
		if(maxVal <= threshold) break;

		double theta = (fabs( M.data[ pivI + 3*pivI ] - M.data[ pivJ + 3*pivJ ]) >0.0001) ? 0.5 * atan(- 2.0 * M.data[ pivI + 3*pivJ ] / ( M.data[ pivI + 3*pivI ] - M.data[ pivJ + 3*pivJ ]) ) :
			                                                                                0.25 * M_PI ;
		Utmp.SetIdentity();
		Utmp.data[pivI + 3*pivI ] =   Utmp.data[pivJ + 3*pivJ ] = cos(theta);
		Utmp.data[pivI + 3*pivJ ] =   sin(theta);
		Utmp.data[pivJ + 3*pivI ] = - Utmp.data[pivI + 3*pivJ ];

		U *= Utmp;
		//M = Utmp.getTransposed() * M * Utmp;重い・・・
		tmp.setMultM1M2(M,Utmp);
		Utmp.transpose_Self();
		M.setMultM1M2(Utmp,tmp);
	}
}
#else
inline void t_Jacobi_usePreU(TMatrix9 &M, TMatrix9& U, double threshold =0.0000001, int maxIteration = 100){
	t_Jacobi_usePreU(M.data, U.data, threshold, maxIteration );
}
#endif

inline void t_Jacobi(TMatrix4 &M, TMatrix4& U, double threshold = 0.000001, int maxIteration = 100)
{
	U.SetIdentity();

	TMatrix4 Utmp,tmp;
	for(int count = 0; count < maxIteration; ++count){
		if(fabs(M.data[2]) <= threshold) break;

		double theta = M_PI / 4.0;
		double aii   = M.data[ 0 ];
		double aij   = M.data[ 2 ];
		double ajj   = M.data[ 3 ];

		if( fabs(aii - ajj) >0.0001) theta = 0.5 * atan(- 2.0 * aij / ( aii - ajj) );
		
		Utmp.data[ 0 ] =  cos(theta);
		Utmp.data[ 3 ] =  cos(theta);
		Utmp.data[ 2 ] =  sin(theta);
		Utmp.data[ 1 ] = -sin(theta);

		U *= Utmp;

		//M = Utmp_t * M * Utmp
		tmp.setMultM1M2(M,Utmp);
		Utmp.transpose_Self();
		M.setMultM1M2(Utmp,tmp);
	}
}

inline void t_Jacobi1(TMatrix4 &M, TMatrix4& U)
{
	double a = M.data[0], b = M.data[2];
	double c = M.data[1], d = M.data[3];

	if( fabs(b) < 0.000001 && fabs(c) < 0.0000001){ 
		U.SetIdentity(); return; 
	}

	double b2_4ac = (a+d)*(a+d) - 4*(a*d-b*c);
	if( b2_4ac < 0) {M.SetIdentity();return;}
	b2_4ac = sqrt( b2_4ac );
	
	double lambda1 = 0.5 * ( (a+d) + b2_4ac);
	double lambda2 = 0.5 * ( (a+d) - b2_4ac);

	//if( fabs( lambda1 - lambda2) < 0.0000001){ M.dataU.SetIdentity(); return; }
	
	U.data[1] = U.data[3] = 1;
	U.data[0] = ( fabs( a-lambda1 ) > fabs( c ) ) ? ( - b / (a-lambda1) ) : ((lambda1-d) / c);
	U.data[2] = ( fabs( a-lambda2 ) > fabs( c ) ) ? ( - b / (a-lambda2) ) : ((lambda2-d) / c);

	double l1 = sqrt( U.data[0]*U.data[0] + U.data[1]*U.data[1] );
	double l2 = sqrt( U.data[2]*U.data[2] + U.data[3]*U.data[3] );
	U.data[0] /= l1;
	U.data[1] /= l1;
	U.data[2] /= l2;
	U.data[3] /= l2;

	M.data[0] = lambda1;
	M.data[1] = 0;
	M.data[2] = 0;
	M.data[3] = lambda2;
}

inline void t_extractRotationTermFromMatrix(const TMatrix9 &M, TMatrix9 &R)
{
	TMatrix9 MtM, U, MtMrt;
	MtM.setMultM1tM2( M, M );

	t_Jacobi( MtM, U, 0.000001, 8);
	
    MtM.data[0] = sqrt( abs( MtM.data[0]) );//普通はabsいらない(0付近になるときに誤差の関係に負になることがある)
    MtM.data[4] = sqrt( abs( MtM.data[4]) );
    MtM.data[8] = sqrt( abs( MtM.data[8]) );
	MtM.data[1] = MtM.data[2] = MtM.data[3] =
	MtM.data[5] = MtM.data[6] = MtM.data[7] = 0;

	MtMrt = U * MtM * U.getTransposed();
	
	MtMrt.getInvertSelf();
	R.setMultM1M2( M, MtMrt );
}

inline void t_extractRotationAndScalingTermFromMatrix(const TMatrix9 &M, TMatrix9 &R)
{
	TMatrix9 MtM, U, MtMrt;
	MtM.setMultM1tM2( M, M );

	t_Jacobi( MtM, U, 0.000001, 8);

	double scale = MtM.data[0] * MtM.data[0] + MtM.data[1] * MtM.data[1] + MtM.data[2] * MtM.data[2];
	scale = sqrt( scale );

    MtM.data[0] = sqrt( abs( MtM.data[0]) );//普通はabsいらない(0付近になるときに誤差の関係に負になることがある)
    MtM.data[4] = sqrt( abs( MtM.data[4]) );
    MtM.data[8] = sqrt( abs( MtM.data[8]) );
	MtM.data[1] = MtM.data[2] = MtM.data[3] =
	MtM.data[5] = MtM.data[6] = MtM.data[7] = 0;

	MtMrt = U * MtM * U.getTransposed();
	
	MtMrt.getInvertSelf();
	R.setMultM1M2( M, MtMrt );
	for( int i = 0; i < 9 ; ++i) R.data[i] *= scale;
}

inline void t_PCA_3D(const vector<TVector3> &points, TVector3 &v1, TVector3 &v2, TVector3 &v3)
{
	TVector3 gCenter;
	for(int i=0; i<(int)points.size(); ++i) gCenter += points[i];
	gCenter /= (int) points.size();

	TVector3 x;
	TMatrix9 covM, U; // Σ x * x_t
	for(int i=0; i<(int)points.size(); ++i) {
		x.SetSubtract( points[i], gCenter );
		covM.addMultVectors( x, x );
	}

	t_Jacobi( covM, U );
	if( covM.data[0] >= covM.data[4] && covM.data[4] >= covM.data[8]      ){//123
		v1.Set(U.data[0], U.data[1], U.data[2]);//1
		v2.Set(U.data[3], U.data[4], U.data[5]);//2
		v3.Set(U.data[6], U.data[7], U.data[8]);//3
	}else if( covM.data[0] >= covM.data[8] && covM.data[8] >= covM.data[4]){//132
		v1.Set(U.data[0], U.data[1], U.data[2]);//1
		v2.Set(U.data[6], U.data[7], U.data[8]);//3
		v3.Set(U.data[3], U.data[4], U.data[5]);//2
	}else if( covM.data[4] >= covM.data[0] && covM.data[0] >= covM.data[8]){//213
		v1.Set(U.data[3], U.data[4], U.data[5]);//2
		v2.Set(U.data[0], U.data[1], U.data[2]);//1
		v3.Set(U.data[6], U.data[7], U.data[8]);//3
	}else if( covM.data[4] >= covM.data[8] && covM.data[8] >= covM.data[0]){//231
		v1.Set(U.data[3], U.data[4], U.data[5]);//2
		v2.Set(U.data[6], U.data[7], U.data[8]);//3
		v3.Set(U.data[0], U.data[1], U.data[2]);//1
	}else if( covM.data[8] >= covM.data[0] && covM.data[0] >= covM.data[4]){//312
		v1.Set(U.data[6], U.data[7], U.data[8]);//3
		v2.Set(U.data[0], U.data[1], U.data[2]);//1
		v3.Set(U.data[3], U.data[4], U.data[5]);//2
	}else/* if( covM.data[8] >= covM.data[4] && covM.data[4] >= covM.data[0])*/{//321
		v1.Set(U.data[6], U.data[7], U.data[8]);//3
		v2.Set(U.data[3], U.data[4], U.data[5]);//2
		v3.Set(U.data[0], U.data[1], U.data[2]);//1
	}
}

//2辺の垂直二等分線の交点が中心　外接円
inline void t_getCircumCircle2D( const TVector3 &x0, const TVector3 &x1, const TVector3 &x2, TVector3 &center, double &r)
{
	double a  = x0.data[0] - x1.data[0];
	double b  = x0.data[1] - x1.data[1];
	double c  = x1.data[0] - x2.data[0];
	double d  = x1.data[1] - x2.data[1];

	double x1Len = x1.Length_Square();
	double b0 = 0.5 * ( x0.Length_Square() - x1Len              );
	double b1 = 0.5 * ( x1Len              - x2.Length_Square() );

	double x,y;
	t_solve2by2LinearEquation( a,b,c,d, b0, b1, x,y);

	center.data[0] = x;
	center.data[1] = y;
	center.data[2] = 0;
	
	r = t_distance( x0, center );
}
/////////////////////////////////////////////////////////////////////////////////////////////
//intersections//////////////////////////////////////////////////////////////////////////////
inline bool t_intersect2D_Lines_preCheck( const TVector3 &x0, const TVector3 &x1,
							              const TVector3 &p0, const TVector3 &p1 )
{
	//x 
	if( x0.data[0] < x1.data[0] ){
		if( p0.data[0] <= x0.data[0] && p1.data[0] <= x0.data[0] ) return false;
		if( p0.data[0] >= x1.data[0] && p1.data[0] >= x1.data[0] ) return false;
	}else{
		if( p0.data[0] <= x1.data[0] && p1.data[0] <= x1.data[0] ) return false;
		if( p0.data[0] >= x0.data[0] && p1.data[0] >= x0.data[0] ) return false;
	}
	//y
	if( x0.data[1] < x1.data[1] ){
		if( p0.data[1] <= x0.data[1] && p1.data[1] <= x0.data[1] ) return false;
		if( p0.data[1] >= x1.data[1] && p1.data[1] >= x1.data[1] ) return false;
	}else{
		if( p0.data[1] <= x1.data[1] && p1.data[1] <= x1.data[1] ) return false;
		if( p0.data[1] >= x0.data[1] && p1.data[1] >= x0.data[1] ) return false;
	}
	return true;
}

inline bool t_intersect2D_Lines( const TVector3 &x0, const TVector3 &d0, const TVector3 &x1, const TVector3 &d1,TVector3 &result )
{
	double t0, t1;
	double a = d0.data[0] , b = -d1.data[0];
	double c = d0.data[1] , d = -d1.data[1];
	double u = x1.data[0] - x0.data[0];
	double v = x1.data[1] - x0.data[1];

	if( !t_solve2by2LinearEquation(a,b,c,d, u,v, t0, t1) ) return false;
	
	result.data[0] = x0.data[0] + t0 * d0.data[0];
	result.data[1] = x0.data[1] + t0 * d0.data[1];
	result.data[2] = 0; 
	return true;
}
//method return true if there may be a intersection
inline bool t_intersect2D_LineSegments( const TVector3 &x0, const TVector3 &d0, 
										const TVector3 &x1, const TVector3 &d1,
											  TVector3 &result )
{
	double t0, t1;
	double a = d0.data[0] , b = -d1.data[0];
	double c = d0.data[1] , d = -d1.data[1];
	double u = x1.data[0] - x0.data[0];
	double v = x1.data[1] - x0.data[1];

	if( !t_solve2by2LinearEquation(a,b,c,d, u,v, t0, t1) ) return false;
	
	if( 0 <= t0 && t0 <= 1 && 
		0 <= t1 && t1 <= 1 )
	{
		result.data[0] = x0.data[0] + t0 * d0.data[0];
		result.data[1] = x0.data[1] + t0 * d0.data[1];
		result.data[2] = 0; 

		fprintf( stderr, "t=%f %f\n", t0, t1);
		return true;
	}
	return false;
}

//strokeのbounding boxとv0-v1の交点の有無をチェックした方が速い(todo)
inline bool t_intersect2D_strokeToLineSegment( const vector<TVector3> &str, const TVector3 &v0, const TVector3 &v1, TVector3 &result)
{
	TVector3 dir = v1-v0, dx;
	for( int i=1,s=(int)str.size(); i<s; ++i)
	{
		const TVector3 &x  = str[i-1];
		dx.SetSubtract( str[ i ], str[i-1]);
		if( t_intersect2D_LineSegments( x,dx, v0,dir, result) ) return true;
	}
	return false;
}


//sphere center: s0, sphere radius: sR  Ray: p0 + t*norm_d(normalized), pos0,pos1 are the intersection points
inline bool t_intersect_RayAndSphere( const TVector3 &s0, double sR, 
	                                  const TVector3 &rayP, const TVector3 &rayD_norm,
							         TVector3 &pos0, TVector3 &pos1)
{
	TVector3 P; 

	P.SetSubtract( rayP, s0 );
	double a = rayD_norm * rayD_norm;
	double b = 2 * ( rayD_norm * P);
	double c = P * P - sR * sR;

	//a x*x + bx + c = 0を得
	double D = b * b - 4 * a * c;
	if( D < 0 ) return false;
	
	double t0 = (-b + sqrt( D )) / (2 * a);
	double t1 = (-b - sqrt( D )) / (2 * a);
	
	pos0.Set_V1_Add_CoefMultV2( rayP, t0, rayD_norm );
	pos1.Set_V1_Add_CoefMultV2( rayP, t1, rayD_norm );
	return true;
}
inline void t_intersect_RayAndTriangle( const TVector3 &lineP, const TVector3 &lineD,
							            const TVector3 &x0   , const TVector3 &x1   , const TVector3 &x2, TVector3 &trgt )
{
	TMatrix9 M;
	TVector3 b, tsu;
	M.Set(  x1.data[0] - x0.data[0],   x2.data[0] - x0.data[0],   - lineD.data[0], 
			x1.data[1] - x0.data[1],   x2.data[1] - x0.data[1],   - lineD.data[1], 
			x1.data[2] - x0.data[2],   x2.data[2] - x0.data[2],   - lineD.data[2]);
	if( !M.getInvertSelf() ){fprintf( stderr, "error12332\n" ); trgt.Set(0.333,0.333,0.333); return;}
	b.SetSubtract( lineP, x0 );
	t_MatMultVec( M, b, tsu );
	trgt.Set_V1_Add_CoefMultV2( lineP, tsu.data[2], lineD);
}
//1線分 p1(x1,y1) p2(x2,y2)にポイントが乗っていたらtrue
inline bool t_isOnLine2D(const TVector3 &p1, const TVector3 &p2, const TVector3 &point)
{
	double x1 = p1   .data[0], y1 = p1   .data[1];
	double x2 = p2   .data[0], y2 = p2   .data[1];
	double x  = point.data[0], y  = point.data[1];
	if(x1 == x2 && y1 == y2){ 
		if(x == x1 && y == y1) return true ;		
		else                   return false;
	}
	if( ((x - x1) * (x - x2) <= 0) && ((y - y1) * (y - y2) <= 0)){
		double a = y2-y1;
		double b = -(x2-x1);
		double c = y1 * x2 - x1 * y2;
		double p = abs(a * x + b * y + c);
		p = p / sqrt(a * a + b * b);
		
		if( p <= 0.001) return true ;//本来は== 0とするべきだが。		
		else            return false;
	}else
		return false;
}

inline bool t_intersect_LineSegmentAndTriangle( const TVector3 &x0, const TVector3 &x1, const TVector3 &x2,
										        const TVector3 &p0, const TVector3 &p1,
										       TVector3 &h)
										
{
	TVector3 v1, v2, d, P;
	TMatrix9 M;

	v1.SetSubtract(x1,x0);
	v2.SetSubtract(x2,x0);
	d .SetSubtract(p1,p0);
	P .SetSubtract(p0,x0);
	
	M.Set(v1.data[0], v2.data[0], -d.data[0],
		  v1.data[1], v2.data[1], -d.data[1],
          v1.data[2], v2.data[2], -d.data[2]);
	
	if(!M.getInvertSelf()) return false;

	double a = M.data[0] * P.data[0] + M.data[3] * P.data[1] + M.data[6] * P.data[2];//abt = Minv * p0
	double b = M.data[1] * P.data[0] + M.data[4] * P.data[1] + M.data[7] * P.data[2];
	double t = M.data[2] * P.data[0] + M.data[5] * P.data[1] + M.data[8] * P.data[2];

	if( a >=   0 && b >=   0 && t >= 0   &&
		a <= 1.0 && b <= 1.0 && t <= 1.0 &&
		a + b <= 1.0)
	{
		h.data[0] = p0.data[0] + t * d.data[0];
		h.data[1] = p0.data[1] + t * d.data[1];
		h.data[2] = p0.data[2] + t * d.data[2];
		return true;
	}
	return false;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//vertices & strokes (verticesは連続性の必要ないものに対して/strokeは連続性のひつようなもの)/////////////////////////////////////////////////////////////////////////////////////////////////
inline void t_verts_getNearestPoint(const vector< TVector3 > &verts, const TVector3 &pos, int &idx, double &dist)
{
	dist = DBL_MAX;
	idx  = -1;
	for( int i = 0; i < (int) verts.size(); ++i)
	{
		double d = t_distance_sq( pos, verts[i] );
		if( d < dist ) { dist = d; idx  = i;}
	}
}
inline void t_verts_getNearestPoint(const int vSize, const TVector3 *verts, const TVector3 &pos, int &idx, double &dist)
{
	dist = DBL_MAX;
	idx  = -1;
	for( int i = 0; i < vSize; ++i){
		double d = t_distance_sq( pos, verts[i] );
		if( d < dist ) { dist = d; idx  = i;}
	}
}
inline void t_verts_calcBoundary(const vector< TVector3 > &verts,  TVector3 &minima, TVector3 &maxima)
{
	minima.Set( DBL_MAX, DBL_MAX, DBL_MAX );
	maxima.Set(-DBL_MAX,-DBL_MAX,-DBL_MAX );
	for( vector< TVector3 >::const_iterator it = verts.begin(); it != verts.end(); ++it)
	{
		for( int i = 0; i < 3; ++i){
			minima.data[i] = MIN( minima.data[i], it->data[i] );
			maxima.data[i] = MAX( maxima.data[i], it->data[i] );
		}
	}
}
inline void t_verts_calcBoundary( int vSize, const TVector3* verts, TVector3 &minB, TVector3 &maxB)
{
	minB.Set( DBL_MAX, DBL_MAX, DBL_MAX );
	maxB.Set(-DBL_MAX,-DBL_MAX,-DBL_MAX );
	for( int i=0; i<vSize; ++i ){
		minB.data[0] = min( minB.data[0], verts[i].data[0] );
		minB.data[1] = min( minB.data[1], verts[i].data[1] );
		minB.data[2] = min( minB.data[2], verts[i].data[2] );
		maxB.data[0] = max( maxB.data[0], verts[i].data[0] );
		maxB.data[1] = max( maxB.data[1], verts[i].data[1] );
		maxB.data[2] = max( maxB.data[2], verts[i].data[2] );
	}
}

inline void t_scaleVertices( vector< TVector3 > &vertices, double scale          ){
	for( vector< TVector3 >::iterator it = vertices.begin(); it != vertices.end(); ++it)
		(*it) *= scale;
}
inline void t_translateVertices( vector< TVector3 > &vertices, const TVector3 &translate ){
	for( vector< TVector3 >::iterator it = vertices.begin(); it != vertices.end(); ++it)
		(*it) += translate;
}

inline void t_fitVertsIntoBox( vector< TVector3 > &verts, const TVector3 &center, double edgeLength)
{
	if( verts.size() <= 1 ) return;

	TVector3 minPos, maxPos, gCenter; 
	t_verts_calcBoundary( verts, minPos, maxPos );

	for( int i=0;i<(int) verts.size(); ++i) gCenter += verts[i]; 
	gCenter /= (double) verts.size();
	
	double scale = 0.5 * edgeLength / MAX3( maxPos.data[0] - minPos.data[0],
									        maxPos.data[1] - minPos.data[1], 
									        maxPos.data[2] - minPos.data[2] );
	t_translateVertices( verts, gCenter );
	t_scaleVertices    ( verts, scale   );
	t_translateVertices( verts, center  );
}
////////////////////////////////////////////////////////////////////////////////////////////////////
////STROKE//////////////////////////////////////////////////////////////////////////////////////////
inline double t_stroke_Length(const vector<TVector3> &str, bool bClosed = false)
{
	double d = 0;
	if( bClosed && str.size() >= 2) d += t_distance( str.back(), str.front() );
	for( int i = 1; i < (int) str.size(); ++i ) d += t_distance( str[i], str[i-1] );
	return d;
}
inline bool t_stroke2D_closedStrokeHasSelfCollision( const vector<TVector3> &stroke2D )
{
	TVector3 dv,dk,res;
	for( int i=0  ; i<(int)stroke2D.size()-1;++i)
	for( int j=i+2; j<(int)stroke2D.size()  ;++j)if( !(i==0 && j== stroke2D.size()-1) )
	{
		const TVector3 &v0 = stroke2D[ i ]; 
		const TVector3 &v1 = stroke2D[i+1];
		const TVector3 &k0 = stroke2D[ j ]; 
		const TVector3 &k1 = stroke2D[(j==stroke2D.size()-1) ? 0:j+1]; 
		dv.SetSubtract(v1,v0);
		dk.SetSubtract(k1,k0);
		if( t_intersect2D_LineSegments( v0,dv, k0,dk, res) ) return true;
	}
	return false;
}
inline bool t_stroke2D_pointInsideStroke( const TVector3 &pos, const vector<TVector3> &stroke, double &sumOfAngles)
{
	TVector3 v0, v1;
	sumOfAngles = 0;
	for( int i=0; i<(int) stroke.size() -1; ++i)
	{
		v0.SetSubtract( stroke[ i ], pos );
		v1.SetSubtract( stroke[i+1], pos );

		sumOfAngles += t_getAngle2D( v0, v1 );
	}

	v0.SetSubtract( stroke.back() , pos );
	v1.SetSubtract( stroke.front(), pos );
	sumOfAngles += t_getAngle2D( v0, v1 );

	if     ( sumOfAngles > M_PI * 1.5 ) return true ;
	else if( sumOfAngles <-M_PI * 1.5 ) return true ;
	else                                return false;
}


/* pos    : 2D position 
   stroke : 2D stroke
   this function return true when pos is inside of the closed poly line (stroke )and 
                       false when pos is outside of the stroke */
inline bool t_stroke2D_pointInsideStroke( const TVector3 &pos, const vector<TVector3> &stroke){
	double sum = 0;
	return t_stroke2D_pointInsideStroke( pos, stroke, sum);
}
/*--------------------------------------------------------------------------------
//devide the stroke into "n" section so that they are equaly spaced along itself.
//This method return n+1 points;
//Strokeをｎ等分する。帰ってくる点列はn+1。
--------------------------------------------------------------------------------*/
inline void t_stroke_devideEquals(int n, const vector<TVector3> &stroke ,vector<TVector3> &result)
{
	if( stroke.size() < 2 ) {result.resize( n + 1); return;}

	double strokeLength = t_stroke_Length( stroke  ); 
	double stepD        = strokeLength / ((double) n);
	
	if(stepD == 0){	 result.clear(); result.resize(n+1); return; }

	result.clear();

	double distance  = 0;
	result.push_back( stroke[0]);
	TVector3 pivot = stroke[0];

	TVector3 vec, point;
	for( int index = 1 ; index < (int)stroke.size();)
	{
		const TVector3 &newP = stroke[index];
		distance += t_distance( newP, pivot );

		if( distance >= stepD )//踏み越え
		{
			vec.SetSubtract( pivot, newP ); 
			vec.Normalize_Self();
			vec *= (distance - stepD);//踏み越え分.;

			point.SetAddition(newP, vec);
			result.push_back(point);
			//さらに次に向けて更新.indexはインクリメントしない
			distance = 0;
			pivot.Set(point);
		}
		else                   //踏み越えていない
		{
			pivot.Set( newP );
			++index;
		}
	}
	//最後の要素が誤差の関係で入れられてないので、入れる
	if(result.size() != n + 1) result.push_back( stroke.back() );
}
/*--------------------------------------------------------------------------------
//devide the stroke into "n" section so that they are equaly spaced along itself.
//This method return n+1 points;
//closed Strokeをｎ等分する。帰ってくる点列はn
--------------------------------------------------------------------------------*/
inline void t_stroke_devideEquals_closed(int n, const vector<TVector3> &stroke ,vector<TVector3> &result)
{
	if( stroke.size() < 2 ) {result.resize( n ); return;}

	double strokeLength = t_stroke_Length( stroke  ) + t_distance(stroke.front(), stroke.back()); 
	double stepD        = strokeLength / (double) n;
	result.clear();
	
	if(stepD == 0){	result.resize(n); return; }


	double distance  = 0;
	result.push_back( stroke[0]);
	TVector3 pivot = stroke[0];

	TVector3 vec, point;
	for( int index = 1, sSize = (int)stroke.size(); index <= sSize;)
	{
		const TVector3 &newP = (index == sSize) ? stroke.front() : stroke[index];
		distance += t_distance( newP, pivot );

		if( distance >= stepD )//踏み越え
		{
			vec.SetSubtract( pivot, newP ); 
			vec.Normalize_Self();
			vec *= (distance - stepD);//踏み越え分.;

			point.SetAddition(newP, vec);
			result.push_back(point);
			//さらに次に向けて更新.indexはインクリメントしない
			distance = 0;
			pivot.Set(point);
			if( result.size() == n ) return;
		}
		else                   //踏み越えていない
		{
			pivot.Set( newP );
			++index;
		}
	}
	//最後の要素が誤差の関係で入れられてないので、入れる
	if(result.size() != n ) { fprintf( stderr, "something wrong error 14231\n");}
}

inline void t_stroke_calcTangent( const vector<TVector3> &str, const int &idx, const bool &bClosed, TVector3 &dir)
{
	dir.Set( 0,0,0 ); 
	const int s = (int)str.size();
	if( bClosed ){ if( idx == 0   ) dir.AddSubtract( str[ 0   ], str[s-1  ] );
				   else             dir.AddSubtract( str[ idx ], str[idx-1] );
				   if( idx == s-1 ) dir.AddSubtract( str[ 0   ], str[idx  ] );
				   else             dir.AddSubtract( str[idx+1], str[idx  ] );}
	else{          if( idx != 0   ) dir.AddSubtract( str[ idx ], str[idx-1] );
			       if( idx != s-1 ) dir.AddSubtract( str[idx+1], str[ idx ] );}
	dir.Normalize_Self();
}

//Smoothing stroke using simple algorithm
inline void t_stroke_Smoothing( vector< TVector3 > &stroke )
{
	vector<TVector3> result( stroke.size() );
	result[0] = stroke.front();
	for (int i = 1; i < (int)stroke.size() - 1; i++)  result[i].SetAdditionWithCoef( 0.5 , stroke[ i ], 
																				     0.25, stroke[i-1], 
																				     0.25, stroke[i+1] );
	result.back() = stroke.back();
	stroke = result;
}
//Smoothing stroke using simple algorithm
inline void t_Smoothing( vector<double> &verts)
{
	vector<double> result( verts.size() );
	result[0] = verts.front();
	for (int i = 1; i < (int)verts.size() - 1; i++)  verts[i] =  0.5 * verts[ i ] + 0.25* verts[i-1] + 
																				    0.25* verts[i+1] ;
	result.back() = verts.back();
	verts = result;
}
inline void t_stroke_Smoothing( vector<CPoint> &points){
	vector< CPoint > result( points.size() );
	result.front() = points.front();
	result.back()  = points.back ();
	for(int i = 1; i < (int)points.size() - 1; ++i)
	{
		result[i].x = (int)( ( points[ i ].x * 2.0 + points[i-1].x + points[i+1].x) / 4.0);
		result[i].y = (int)( ( points[ i ].y * 2.0 + points[i-1].y + points[i+1].y) / 4.0);
	}
	points.swap( result);
}




/*-----------------------------------------------------------
                   p
                   |
                   |
        *----------*---------------*
     p[idx1]                    p[idx2]      
	p = p[idx1] + t * d
	d = p[idx2] - p[idx1]
	
TVector3 &pos : found point 
-----------------------------------------------------------*/
inline bool t_stroke_getNearestPoint( const TVector3 &p, const vector< TVector3 > &str, 
	                                   TVector3 &pos, double &mindist_sq, int &pivIdx1, int &pivIdx2, double &tVal/*[0,1]*/, double &lenAlongStr)
{
	bool found = false;
	TVector3 d, tmpPos;
	mindist_sq = DBL_MAX;
	int tmpID1, tmpID2;
	double totalLen = 0;
	for( int i = 1; i < (int) str.size(); ++i)
	{
		const TVector3 &x0 = str[ i-1 ];
		const TVector3 &x1 = str[  i  ];
		tmpID1 = i-1;
		tmpID2 =  i ;		

		d.SetSubtract( x1,x0);
		double l       = d.Length_Square();
		double length  = sqrt( l );
		
		if( l == 0) continue;
		double t = ( 1 / l ) * (  ( p.data[0] - x0.data[0] ) * d.data[0] + 
			                      ( p.data[1] - x0.data[1] ) * d.data[1] + 
								  ( p.data[2] - x0.data[2] ) * d.data[2] ) ;
		double distSq = 100000;
		if(      t < 0 ){ distSq = t_distance_sq( x0, p ); tmpPos.Set( x0 ); t = 0;}
		else if( t > 1 ){ distSq = t_distance_sq( x1, p ); tmpPos.Set( x1 ); t = 1;}
		else{
			tmpPos.data[0] = x0.data[0] + t * d.data[0];
			tmpPos.data[1] = x0.data[1] + t * d.data[1];
			tmpPos.data[2] = x0.data[2] + t * d.data[2];
			distSq = t_distance_sq(tmpPos , p);
		}

		if( distSq < mindist_sq )
		{
			pos.Set( tmpPos );
			mindist_sq  = distSq ;
			found       = true;
			pivIdx1     = tmpID1;
			pivIdx2     = tmpID2;
			tVal        = t     ;
			lenAlongStr = totalLen + t * length;
		}
		totalLen += length;
	}
	return found;
}


inline void t_points_getNearestPoint( const TVector3 &p, const vector< TVector3 > &points, int &idx, double &distSq)
{
	idx    = -1;
	distSq = DBL_MAX;
	for( int i=0, s = (int) points.size(); i < s; ++i)
	{
		double d = t_distance( p, points[i]);
		if( d < distSq )
		{
			distSq = d;
			idx    = i;
		}
	}
}


inline bool t_stroke_getNearestPoint( const TVector3 &p, const vector< TVector3 > &str , TVector3 &pos, double &mindist_sq)
{
	int idxtmp1, idxtmp2;
	double t, lenAlongStr;
	return t_stroke_getNearestPoint( p, str , pos, mindist_sq, idxtmp1, idxtmp2, t, lenAlongStr);
}

inline void t_stroke_getNearestPoint( const TVector3 &p, const vector< TVector3 > &str , int &idx, double &distance)
{
	idx = -1;

	distance = 1000000000;
	for( int i=0; i<(int) str.size(); ++i)
	{
		double d = t_distance( p, str[i]);
		if( d < distance )
		{
			distance = d;
			idx = i;
		}
	}
}

inline void t_stroke2D_BoundBox( const vector<TVector3> &str, double b[4])
{
	b[0] =  100000; b[1] =  100000;
	b[2] = -100000; b[3] = -100000;
	for( vector< TVector3 >::const_iterator it = str.begin(); it != str.end(); ++it){
		b[0] = MIN( b[0], it->data[0]);  b[1] = MIN( b[1], it->data[1]);
		b[2] = MAX( b[2], it->data[0]);  b[3] = MAX( b[3], it->data[1]);
	}
}
inline void t_stroke_BoundBox( const vector<TVector3> &str, TVector3 &minPos, TVector3 &maxPos)
{
	minPos.data[0] =  100000;  minPos.data[1] =  100000;  minPos.data[2] =  100000; 
	maxPos.data[0] = -100000;  maxPos.data[1] = -100000;  maxPos.data[2] = -100000;

	for( vector< TVector3 >::const_iterator it = str.begin(); it != str.end(); ++it){
		minPos.data[0] = min( minPos.data[0], it->data[0]);  
		minPos.data[1] = min( minPos.data[1], it->data[1]);  
		minPos.data[2] = min( minPos.data[2], it->data[2]);  
		
		maxPos.data[0] = max( maxPos.data[0], it->data[0]);  
		maxPos.data[1] = max( maxPos.data[1], it->data[1]);  
		maxPos.data[2] = max( maxPos.data[2], it->data[2]); 
	}
}

inline void t_stroke_multMat(const TMatrix16 &M, vector<TVector3> &stroke){	
	for( vector<TVector3>::iterator it = stroke.begin(); it != stroke.end(); ++it) t_MatMultVec( *it, M );
}


//the stroke should be closed form( the function automatically connect str.back() to str.front()
inline bool t_stroke_isClosedStrokeUnClockwise2D(const vector<TVector3> &str)
{
	double yval = -1000000;
	double xval = -1000000;
	int    piv  = -1;
	//get top point
	for( int i = 0; i < (int) str.size(); ++i) 
		if(      yval < str[i].data[1] )
		{
			xval = str[i].data[0];
			yval = str[i].data[1];
			piv  = i;
		}
		else if( yval == str[i].data[1] && xval < str[i].data[0] )
		{
			xval = str[i].data[0];
			yval = str[i].data[1];
			piv  = i;
		}
	TVector3 d1, d2;
	d1.SetSubtract( str[ (piv ==  0                 ) ? (int)str.size() -1 : piv-1 ] ,str[piv]);
	d2.SetSubtract( str[ (piv == (int)str.size() -1 ) ?        0           : piv+1 ] ,str[piv]);
	double angle = t_getAngle2D( d1,d2 );
	if( angle >= 0) return false;
	else            return true ;
}

inline void t_stroke_invert( vector<TVector3> &v)
{
	vector< TVector3 > tmpV;
	for( int i = (int) v.size() -1; i >= 0; --i)
		tmpV.push_back( v[i] );
	v = tmpV;
}

inline void t_stroke_invert(  vector<int> &v)
{
	vector< int > tmpV;
	for( int i = (int) v.size() -1; i >= 0; --i) tmpV.push_back( v[i] );
	v = tmpV;
}
inline double t_stroke_calcSumStrokeAngleAroundPoint( const TVector3 &point, const vector< TVector3 > &str )
{
	static TVector3 dir1, dir2;
	
	double sumOfAngle = 0;

	for( int i = 0; i < (int) str.size(); ++i)
	{
		dir1.SetSubtract( str[ i                                          ], point );
		dir2.SetSubtract( str[ ( i == 0 ) ? (int) str.size() - 1 : i - 1  ], point );
		sumOfAngle += t_getAngle2D( dir1, dir2 );
	}
	return sumOfAngle;
}

inline bool t_stroke_isPointInStroke2D( const TVector3 &point, const vector< TVector3 > &str )
{
	double sumOfAngle = t_stroke_calcSumStrokeAngleAroundPoint( point, str );
	if( ( 2 * M_PI - 0.5 < sumOfAngle && sumOfAngle < 2 * M_PI + 0.5 ) || 
		(-2 * M_PI - 0.5 < sumOfAngle && sumOfAngle <-2 * M_PI + 0.5 ) ) return true;
	return                                                                     false;
}

inline void t_stroke_getGravityCenter( const vector<TVector3> &str, TVector3 &gCenter)
{
	gCenter.Set(0,0,0);
	for( int i=0, s = (int)str.size();i<s;++i) gCenter += str[i];
	gCenter /= (double) str.size();
}

//2D strokeとscreenに射影されたfiberの交差する点を検索
inline void t_stroke2D_crossingPt_toFiber( const vector< TVector3 > &stroke2D, 
										   const vector< TVector3 > &fiber2D ,
										   bool                        fib_closed ,
										  vector<pair<int,double>> &crossPointOnFiber)
{
	crossPointOnFiber.clear();

	TVector3 xDir, pDir;
	for( int i1 = 0; i1< (int)stroke2D.size()-1                           ; ++i1 )
	for( int i2 = 0; i2< (int)(fib_closed?fiber2D.size():fiber2D.size()-1); ++i2 )
	{
		const TVector3 &x0 = stroke2D[i1  ],  &p0 = fiber2D[                                 i2  ];
		const TVector3 &x1 = stroke2D[i1+1],  &p1 = fiber2D[ (i2==(int)fiber2D.size()-1) ? 0:i2+1];

		if( !t_intersect2D_Lines_preCheck(x0, x1, p0, p1) ) continue;
		xDir.SetSubtract( x1, x0);
		pDir.SetSubtract( p1, p0);
		
		double s, t;
		//(xDir -pDir)(t,s) = -x0 + p0    交点h = x0 + t*xDir = p0 + sP0
		if( t_solve2by2LinearEquation(  xDir.data[0], -pDir.data[0], xDir.data[1], -pDir.data[1],
									   -x0.data[0] + p0.data[0], -x0.data[1] + p0.data[1], t,s) && 
									    0<=s && s<1 && 0<=t && t<1)
		{
			crossPointOnFiber.push_back( pair<int,double>( i2,s ) );
		}
	}

	//全く同じsetが発見されていたら削除
	while( true ) 
	{
		bool found = false;
		for( int piv1=   0  ; piv1<(int) crossPointOnFiber.size(); ++piv1)
		for( int piv2=piv1+1; piv2<(int) crossPointOnFiber.size(); ++piv2)
		if( crossPointOnFiber[piv1].first  == crossPointOnFiber[piv2].first )
		{
			fprintf( stderr, "here is a strange crossPoint" );
			crossPointOnFiber.erase( crossPointOnFiber.begin() + piv2 );
			found = true;break;
		}
		if( !found ) break;
	}
	
}
//Stroke Manipulation//////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//save and load/////////////////////////////////////////////////////////////////
inline bool t_SaveArrayIntToBin(const char*fname, vector<int> &intValues)
{
	FILE *fp = fopen( fname, "wb");

	unsigned int sizeOfArray = (unsigned int)intValues.size();
	fwrite( &sizeOfArray, sizeof(unsigned int),1,fp );
	int *p = new int[sizeOfArray];
	for( unsigned int i=0;i<sizeOfArray;++i) p[i] = intValues[i];	
	fwrite( p, sizeof( int )*sizeOfArray, 1, fp ) ;

	delete[] p;
	fclose( fp );
	return true;
}
inline bool t_LoadArrayIntToBin(const char*fname, vector<int> &intValues)
{
	FILE *fp = fopen( fname, "rb");
	if( fp == 0 ){ return false;}

	intValues.clear();
	unsigned int sizeOfArray = 0;
	fread( &sizeOfArray, sizeof(unsigned int),1,fp ) ;

	if( sizeOfArray <= 0 ) { fclose(fp); return false; }
	int *p = new int[ sizeOfArray ];
	fread( p, sizeof(int)*sizeOfArray,1 , fp ) ;

	intValues.resize(sizeOfArray);
	for(unsigned int i=0; i<sizeOfArray; ++i) intValues[i] = p[i];	

	delete[] p;
	fclose(fp); 
	return true;
}

inline void t_save3DVectors        ( FILE *fp, const vector<TVector3> &vec)
{
	fprintf( fp, "3DVectors %d\n", (int) vec.size() );
	for( int i=0; i<(int)vec.size(); ++i) fprintf( fp, "%f %f %f\n", vec[i].data[0], vec[i].data[1], vec[i].data[2]);
}
inline void t_load3DVectors        ( FILE *fp,       vector<TVector3> &vec)
{
	char buf[256];
	int  vSize;
	fscanf( fp, "%s%d", buf, &vSize);   fprintf( stderr, "%s %d\n", buf, vSize);
	vec.resize( vSize );
	for( int i=0; i<vSize; ++i) fscanf( fp, "%lf%lf%lf", &vec[i].data[0], &vec[i].data[1], &vec[i].data[2]);
}


// ( A * qi - pi )^2を最小化するようなaffine Aを求める　
//qiが動かしたいvertex(cloth simulationの結果)
//piがfitしたいvertex (Initial shape)
inline void t_MesshlessDefomationsRotationMatrix( const vector<TVector3> &vertices1, const vector<TVector3> &vertices2, 
	                                             TMatrix16 &R, TVector3 &X0cm, TVector3 &Xcm)
{
	//fprintf(stderr, "start calculation for meshless deformation's rotation!!\n"	);
	X0cm.Set(0,0,0);
	Xcm. Set(0,0,0);
	for(int i = 0; i < (int) vertices1.size();++i){ X0cm  += vertices1[i];
		                                            Xcm   += vertices2[i]; }
	double  numInv = 1.0 /(double) vertices1.size();
	X0cm *= numInv;
	Xcm  *= numInv;	
	//fprintf(stderr, "xcm  ") ;Xcm.Trace();
	//fprintf(stderr, "x0cm  ");X0cm.Trace();
	TMatrix9 Apq(0,0,0,0,0,0,0,0,0);
	TVector3 p, q;
	for(int i = 0; i < (int) vertices1.size();++i){
		q.SetSubtract( vertices1[i], X0cm);
		p.SetSubtract( vertices2[i], Xcm);
		//fprintf(stderr, "p  ");p.Trace();
		//fprintf(stderr, "q  ");q.Trace();
		Apq.addMultVectors(p , q);
	}
	Apq = numInv * Apq;
	//fprintf(stderr, "-------------------------\nApq\n");
	TMatrix9 Apq_T(Apq.data), ApqT_Apq;
	Apq_T.transpose_Self();
	//Apq.Trace();
	//Apq_T.Trace();
	ApqT_Apq.setMultM1M2(Apq_T, Apq);
	//fprintf(stderr, "-------------------------\nAt * A\n");
	//ApqT_Apq.Trace();
	TMatrix9 U, tmp;
	t_Jacobi(ApqT_Apq,U, 0.000001, 5);
	//fprintf(stderr, "対角化後At * A\n");
	//ApqT_Apq.Trace();
    ApqT_Apq.data[0] = sqrt( fabs(ApqT_Apq.data[0]) );//普通はabsいらない(0付近になるときに誤差の関係に負になることがある)
    ApqT_Apq.data[4] = sqrt( fabs(ApqT_Apq.data[4]) );
    ApqT_Apq.data[8] = sqrt( fabs(ApqT_Apq.data[8]) );
	ApqT_Apq.data[1] = 0;  ApqT_Apq.data[2] = 0;  ApqT_Apq.data[3] = 0;
	ApqT_Apq.data[5] = 0;  ApqT_Apq.data[6] = 0;  ApqT_Apq.data[7] = 0;
	//これが √At * A
	TMatrix9 rootAt_A = U * ApqT_Apq * U.getTransposed();	
	//fprintf(stderr, "rtAt * rtAt \n");
	//(rootAt_A * rootAt_A).Trace();
	if(!rootAt_A.getInvertSelf()){fprintf(stderr, "system couldn't invert matrix\n");}
	TMatrix9 Rtmp = Apq * rootAt_A;
	//fprintf(stderr, "R  det %f\n", Rtmp.calcDet9());
	//Rtmp.Trace();
	R.data[0] = Rtmp.data[0];  R.data[4] = Rtmp.data[3];	R.data[ 8] = Rtmp.data[6];
	R.data[1] = Rtmp.data[1];  R.data[5] = Rtmp.data[4];	R.data[ 9] = Rtmp.data[7];
	R.data[2] = Rtmp.data[2];  R.data[6] = Rtmp.data[5];	R.data[10] = Rtmp.data[8];
}


inline float t_sqrt( const float& x )
{
    float xHalf = 0.5f * x;
    int tmp = 0x5F3759DF - ( *(int*)&x >> 1 );//initial guess
    float xRes = *(float*)&tmp;

 	xRes *= ( 1.5f - ( xHalf * xRes * xRes ) );
	//xRes *= ( 1.5f - ( xHalf * xRes * xRes ) );//コメントアウトを外すと精度が上がる
    return xRes * x;
}

inline double t_sqrt( const double &x) 
{
	double xHalf = 0.5 * x;
	long long int tmp = 0x5FE6EB50C7B537AAl - ( *(long long int*)&x >> 1);//initial guess
	
	double xRes = * (double*)&tmp;
	xRes *= ( 1.5 - ( xHalf * xRes * xRes ) );
	xRes *= ( 1.5 - ( xHalf * xRes * xRes ) );//コメントアウトを外すと精度が上がる
	return xRes * x;
}


/*
以下利用するなら書き直しが必要だが参考にはなりそうなmethods

//uvに対してGrid状に並んだ, spline surface頂点を生成する
//定義域は CP0 - CPn-1で t [1, n-2]
//両短点に2点同じＣＰを加えたことにすると端点をなるべく通るものが出来る
inline void t_createBSplineSurfaceGrid(int cpVert, int cpHori, const vector< vector< TVector3 > > &CPs, 
							    int   vert, int   hori,       vector< TVector3 > &trgtVertex)
{
	double stepU = (cpVert-1) / ( vert -1.0);//本来は u[-1, n+1]とする
	double stepV = (cpHori-1) / ( hori -1.0);//これは u[ 1, n  ]の実装

	for( int vertI = 0; vertI < vert; ++vertI )
	for( int horiI = 0; horiI < hori; ++horiI )
	{
		double u = stepU * vertI;
		double v = stepV * horiI;

		int trgtIdx = hori * vertI + horiI;
		trgtVertex[ trgtIdx ].Set(0,0,0);

		for(int i = -2 ; i <= cpVert -1 + 2 ; ++i) // -2 to n+2
		for(int j = -2 ; j <= cpHori -1 + 2 ; ++j)
		{
			int indexI = (i <= 0)? 0 : (i >= cpVert-1)? cpVert - 1: i;
			int indexJ = (j <= 0)? 0 : (j >= cpHori-1)? cpHori - 1: j;

			const TVector3 &nextP = CPs[indexI][indexJ];
			double g = bSpline(u - i), 
				   l = bSpline(v - j);
			trgtVertex[ trgtIdx ].data[0] += g * l * (nextP.data[0]);
			trgtVertex[ trgtIdx ].data[1] += g * l * (nextP.data[1]);
			trgtVertex[ trgtIdx ].data[2] += g * l * (nextP.data[2]);
		}	
	}
}
*/

//todo 1 順番整理
//todo 3 TMathEx.hとかにしてもいいかも

#endif	// __ILMATH_H_INCLUDED__