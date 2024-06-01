#pragma once


#include "math.h"
#include <vector>
#include <map>
#include <set>
#include <algorithm>

using namespace std;

/*
a graph for graph cut algorithm
by takashi ijiri at 17/6/2009
*/

class colorRGB
{
public:
	byte m_rgb[3];//[0,1] * [0,1] * [0,1];

	~colorRGB(){}
	colorRGB( byte r = 0, byte g = 0, byte b = 0){ m_rgb[0] = r;  m_rgb[1] = g;  m_rgb[2] = b; }
	colorRGB( const colorRGB &c ){
		memcpy( m_rgb, c.m_rgb, sizeof( byte ) * 3 );
	}

	inline double distance2( const colorRGB &c ) const{
		return (double)( (m_rgb[0] - c.m_rgb[0]) * (m_rgb[0] - c.m_rgb[0]) + 
					     (m_rgb[1] - c.m_rgb[1]) * (m_rgb[1] - c.m_rgb[1]) + 
					     (m_rgb[2] - c.m_rgb[2]) * (m_rgb[2] - c.m_rgb[2]) );
	}

	inline void Set( const colorRGB &c ){
		memcpy( m_rgb, c.m_rgb, sizeof( byte ) * 3 );
	}

	inline void Set( const byte* rgba ){
		m_rgb[0] = rgba[0]; 
		m_rgb[1] = rgba[1]; 
		m_rgb[2] = rgba[2]; 
	}
};




class TGCnode
{
	double m_totalR;
	double m_totalG;
	double m_totalB;
public:
	colorRGB    m_color   ; //平均pixel color[0, 1] 
	vector<int> m_pixelIDs;
	bool        m_isFore  ;
	bool        m_isEnable;//partialのみで利用

public:
	~TGCnode(){}
	TGCnode (){
		m_totalR = m_totalG = m_totalB = 0;
		m_isFore   = false;
		m_isEnable = true;
		m_pixelIDs.clear();
	}

	TGCnode( const TGCnode &node)
	{
		m_color .Set( node.m_color  );
		m_totalR   = node.m_totalR;
		m_totalG   = node.m_totalG;
		m_totalB   = node.m_totalB;
		m_pixelIDs = node.m_pixelIDs;
		m_isFore   = node.m_isFore;
		m_isEnable = node.m_isEnable;//partialのみで利用
	}

	inline void addPixel(int pixelIdx, byte R, byte G, byte B)
	{
		m_pixelIDs.push_back( pixelIdx );
		m_totalR += R; 
		m_totalG += G; 
		m_totalB += B;
		m_color.m_rgb[0] = (byte) (m_totalR / (double) m_pixelIDs.size());
		m_color.m_rgb[1] = (byte) (m_totalG / (double) m_pixelIDs.size());
		m_color.m_rgb[2] = (byte) (m_totalB / (double) m_pixelIDs.size());
	}
};


class TGCpixelNode
{
public:
	int      m_pixelID;
	colorRGB m_color;
	inline TGCpixelNode(int pixelIdx, const byte* rgb) : m_color(rgb[0], rgb[1], rgb[2])
	{
		m_pixelID = pixelIdx;
	}
};

class TGCpixelNode_onlyId
{
public:
	int      m_pixelID;
	inline TGCpixelNode_onlyId(int pixelIdx)
	{
		m_pixelID = pixelIdx;
	}
};





/*----------------------------------------------------------------------------------------------------------------
1, watershed による over segmentationからgraph cutを計算できる
2, over segmentation graph cutのseam部分のみをもう一度graph cutかけることができる

以下呼ぶ手順
 createByWsdImage2D()               over segmentationから荒いグラフ構造を構築 (!準備としてWsd後一度だけ!)

 calcMinCut()              荒いグラフ上でgraph cutを計算
 calcMinCut_PixelLebel2D() 荒いグラフで計算後、その境界から幅widthの帯状の部分において pixel levelでgraph cut
----------------------------------------------------------------------------------------------------------------*/
class TGCgraph
{
public:
	static double m_lambda;
	vector<int             >  m_WsdLabels;//pixel         idx --> node (region) idx 

protected:
	vector<map<int, double>>  m_WsdEdges ;//edge i jがある時には, m_edges[i][a].first = j, m_edges[i][?].second = error = E2 
	vector<TGCnode         >  m_WsdNodes ;//node (region) idx --> pixel         idx  ( m_pixelIDs )

public:
	~TGCgraph(void){}
	TGCgraph (void){}
	void testMinCut() const;
	
	//label値は、正整数, 最小=1, 連続 
	void createByWsdImage2D(int width, int height,            const vector<int> &wsd_labels, const byte* rgba );
	void createByWsdImage3D(int width, int height, int depth, const vector<int> &wsd_labels, const byte* rgba );

	void calcMinCut(const set<int> &forePixelIdxs, const set<int> &backPixelIdxs, vector<int> &isForeBack);

	void calcMinCut_PixelLebel2D(int width, int height,            byte* rgba,const set<int> &forePixelIdxs, const set<int> &backPixelIdxs, vector<int> &isForeBack, int boundWidth);
	void calcMinCut_PixelLebel3D(int width, int height, int depth, byte* rgba,const set<int> &foreVoxIdxs  , const set<int> &backVoxIdxs  , vector<int> &isForeBack, int boundWidth);

	void calcMinCur_PixelLevel3D_onlyBoundary( int width, int height, int depth, const byte *rgba, short* volumeInOut, int boundWidth);
	void clear(){m_WsdEdges.clear(); m_WsdNodes.clear(); m_WsdLabels.clear();}

protected:
	inline void addEdge( int nodeIdx1, int nodeIdx2)
	{
		if( nodeIdx1 == nodeIdx2 ) return;
		if( nodeIdx1 >  nodeIdx2 ) std::swap( nodeIdx1, nodeIdx2 );
		m_WsdEdges[nodeIdx1].insert( pair<int,double>( nodeIdx2, 0 ) );
	}
	
	void preCheckAndDivideInconsistentNodes(const set<int> &forePixelIdxs, const set<int> &backPixelIdxs, int sizeW, int sizeH, int sizeD, const byte *rgba);
	void divideInconsistentNodes( int foreVoxIdx, int backVoxIdx, const int sizeW, const int sizeH, const int sizeD, const byte *rgba);
};



/*----------------------------------------------------------------------------------------------------------------

img.RGBAのalphaが0の部分を無視してgraph cut
volumeの一部分だけのgraph cutが可能

 createByWsdImage3D()           (親class)over segmentationから荒いグラフ構造を構築 (!準備としてWsd後一度だけ!)

 updateEnabeWsdNode ()          計算すべきWsdNodeとそうでないWsdNodeをセット
 calcMinCut_part()              荒いグラフ上の一部分のみにおいてgraph cutを計算
 calcMinCut_PixelLebel3D_part() 荒いグラフで計算後、その境界から幅widthの帯状の部分においてpixel levelでgraph cut
----------------------------------------------------------------------------------------------------------------*/

class TGCgraph_partial : public TGCgraph
{
public:
	~TGCgraph_partial(void){}
	TGCgraph_partial (void){}

	void updateEnableWsdNode(byte* rgba);
	void updateEnableWsdNode(byte* rgba, int sizeW, int sizeH, int sizeD, int minX, int minY, int minZ, int maxX, int maxY, int maxZ);
	void calcMinCut_part(const set<int> &forePixelIdxs, const set<int> &backPixelIdxs, vector<int> &isForeBack ,
		                 int width, int height, int depth, const byte *rgba);
	void calcMinCut_PixelLebel3D_part(int width, int height, int depth, byte* rgba,const set<int> &foreVoxIdxs  , const set<int> &backVoxIdxs  , vector<int> &isForeBack, int boundWidth);
};










//error E1_0 E1_1, E2　The segmentation results are very sensitive to coeficients.
inline double calcE2( const TGCnode &p0, const TGCnode &p1 ){
	return  TGCgraph::m_lambda / ( 1.0 +  p0.m_color.distance2( p1.m_color ) );
}

inline double calcE2( const colorRGB &c1, const colorRGB &c2){
	return  TGCgraph::m_lambda / ( 1.0 +  c1.distance2( c2 ) );
}

inline void calcE1(  const vector<colorRGB> &f_Cols,
					 const vector<colorRGB> &b_Cols,
					 const colorRGB &c,
					 double &ef,//fore
					 double &eb) //back
{
	double df = 1000000000;
	double db = 1000000000;
	for(vector<colorRGB>::const_iterator it = f_Cols.begin(); it != f_Cols.end(); ++it)  df = min(df, it->distance2( c ) );
	for(vector<colorRGB>::const_iterator it = b_Cols.begin(); it != b_Cols.end(); ++it)  db = min(db, it->distance2( c ) );
	df = sqrt( df );
	db = sqrt( db );

	if( df == 0 && db == 0 ) 
	{
		ef = eb = 0.5;
	}else{
		ef =  db / ( df + db );
		eb =  df / ( df + db );
	}
}



