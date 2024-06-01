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
	double m_rgb[3];//[0,1] * [0,1] * [0,1];

	~colorRGB(){}
	colorRGB( double r = 0, double g = 0, double b = 0)
	{
		m_rgb[0] = r; m_rgb[1] = g; m_rgb[2] = b;
	}
	colorRGB( const colorRGB &c ){
		memcpy( m_rgb, c.m_rgb, sizeof( double ) * 3 );
	}

	inline double distance2( const colorRGB &c ) const{
		return (m_rgb[0] - c.m_rgb[0]) * (m_rgb[0] - c.m_rgb[0]) + 
               (m_rgb[1] - c.m_rgb[1]) * (m_rgb[1] - c.m_rgb[1]) + 
			   (m_rgb[2] - c.m_rgb[2]) * (m_rgb[2] - c.m_rgb[2]) ;
	}

	inline void Set( const colorRGB &c ){
		memcpy( m_rgb, c.m_rgb, sizeof( double ) * 3 );
	}
	inline void SetByByte( const byte* rgba ){
		m_rgb[0] = rgba[0]/255.0; 
		m_rgb[1] = rgba[1]/255.0; 
		m_rgb[2] = rgba[2]/255.0; 
	}
};


/*

class colorRGB
{
public:
	double m_rgb[3];//[0,1] * [0,1] * [0,1];

	~colorRGB(){}
	colorRGB( double r = 0, double g = 0, double b = 0)
	{
		m_rgb[0] = r; m_rgb[1] = g; m_rgb[2] = b;
	}
	colorRGB( const colorRGB &c ){
		memcpy( m_rgb, c.m_rgb, sizeof( double ) * 3 );
	}

	inline double distance2( const colorRGB &c ) const{
		return (m_rgb[0] - c.m_rgb[0]) * (m_rgb[0] - c.m_rgb[0]) + 
               (m_rgb[1] - c.m_rgb[1]) * (m_rgb[1] - c.m_rgb[1]) + 
			   (m_rgb[2] - c.m_rgb[2]) * (m_rgb[2] - c.m_rgb[2]) ;
	}

	inline void Set( const colorRGB &c ){
		memcpy( m_rgb, c.m_rgb, sizeof( double ) * 3 );
	}
	inline void SetByByte( const byte* rgba ){
		m_rgb[0] = rgba[0]/255.0; 
		m_rgb[1] = rgba[1]/255.0; 
		m_rgb[2] = rgba[2]/255.0; 
	}
};



*/




class TGCnode
{
	colorRGB m_colorT;
public:
	colorRGB    m_color   ; //平均pixel color[0, 1] 
	double      m_E1_f, m_E1_b;
	vector<int> m_pixelIDs;
	bool        m_isFore  ;

public:
	~TGCnode(){}
	TGCnode (){
		m_E1_f = m_E1_b = 0;
		m_isFore = false;
		m_pixelIDs.clear();
	}

	TGCnode( const TGCnode &node)
	{
		m_E1_f = node.m_E1_f; 
		m_E1_b = node.m_E1_b;
		m_color .Set( node.m_color  );
		m_colorT.Set( node.m_colorT );
		for(int i=0; i<(int)m_pixelIDs.size(); ++i) m_pixelIDs.push_back( m_pixelIDs[i] );
	}

	inline void addPixel(int pixelIdx, double R, double G, double B)
	{
		m_pixelIDs.push_back( pixelIdx );
		m_colorT.m_rgb[0] += R; 
		m_colorT.m_rgb[1] += G; 
		m_colorT.m_rgb[2] += B;
		m_color.m_rgb[0] = m_colorT.m_rgb[0] / (double) m_pixelIDs.size();
		m_color.m_rgb[1] = m_colorT.m_rgb[1] / (double) m_pixelIDs.size();
		m_color.m_rgb[2] = m_colorT.m_rgb[2] / (double) m_pixelIDs.size();
	}
};


class TGCpixelNode
{
public:
	int      m_pixelID;
	colorRGB m_color;
	inline TGCpixelNode(int pixelIdx, byte* rgb)
		:m_color(rgb[0] / 255.0, rgb[1] / 255.0, rgb[2] / 255.0 )
	{
		m_pixelID = pixelIdx;
	}
};





/*
1, watershed による over segmentationからgraph cutを計算できる


2, over segmentation graph cutのseam部分のみをもう一度graph cutかけることができる


*/
class TGCgraph
{
private:
	static double m_lambda;
	vector<map<int, double>>  m_WsdEdges ;//edge i jがある時には, m_edges[i][a].first = j, m_edges[i][?].second = error = E2 
	vector<TGCnode         >  m_WsdNodes ;//node (region) idx --> pixel         idx  ( m_pixelIDs )
	vector<int             >  m_WsdLabels;//pixel         idx --> node (region) idx 

private:
	inline void addPixelToNode( int nodeIdx, int pixelIdx, double R, double G, double B){
		m_WsdNodes[nodeIdx].addPixel( pixelIdx, R, G, B);
	}

	inline void addEdge( int nodeIdx1, int nodeIdx2)
	{
		if( nodeIdx1 == nodeIdx2 ) return;
		if( nodeIdx1 >  nodeIdx2 ) std::swap( nodeIdx1, nodeIdx2 );
		m_WsdEdges[nodeIdx1].insert( pair<int,double>( nodeIdx2, 0 ) );
	}

public:
	~TGCgraph(void){}
	TGCgraph (void){}
	void testMinCut() const;
	void clear(){m_WsdEdges.clear(); m_WsdNodes.clear(); m_WsdLabels.clear();}

	//label値は、正整数, 最小=1, 連続 
	void createByWsdImage2D(int width, int height,            const vector<int> &wsd_labels, const byte* rgba );
	void createByWsdImage3D(int width, int height, int depth, const vector<int> &wsd_labels, const byte* rgba );

	//graphcut segmentation in wsd node level
	void calcMinimumGraphCut(const set<int> &forePixelIdxs, 
		                     const set<int> &backPixelIdxs, vector<int> &isForeBack);

	//graphcut segmentation in pixel level
	void calcMinimumGraphCut_PixelLebel2D(int width, int height, byte* rgba,
										  const set<int> &forePixelIdxs, 
		                                  const set<int> &backPixelIdxs, vector<int> &isForeBack);
	void calcMinimumGraphCut_PixelLebel3D(int width, int height, int depth, byte* rgba,
										  const set<int> &foreVoxIdxs, 
		                                  const set<int> &backVoxIdxs, vector<int> &isForeBack);


	//error E1_0 E1_1, E2　The segmentation results are very sensitive to coeficients.
	inline double calcE2( const TGCnode &p0, const TGCnode &p1 ){
		static double b = 255 * 255;
		return  m_lambda / ( 1.0 +  b * p0.m_color.distance2( p1.m_color ) );
	}

	inline double calcE2( const colorRGB &c1, const colorRGB &c2){
		static double b = 255 * 255;
		return  m_lambda / ( 1.0 +  b * c1.distance2( c2 ) );
	}


	inline void calcE1(  const vector<colorRGB> &f_Cols,
						 const vector<colorRGB> &b_Cols,
						 const colorRGB &c,
						 double &ef,//fore
						 double &eb) //back
	{
		static double a = 255 * 255;

		double df = 1000000000;
		double db = 1000000000;
		for(vector<colorRGB>::const_iterator it = f_Cols.begin(); it != f_Cols.end(); ++it)  df = min(df, it->distance2( c ) );
		for(vector<colorRGB>::const_iterator it = b_Cols.begin(); it != b_Cols.end(); ++it)  db = min(db, it->distance2( c ) );
		df = sqrt( a * df );
		db = sqrt( a * db );

		if( df == 0 && db == 0 ) 
		{
			ef = eb = 0.5;
		}else{
			ef =  db / ( df + db );
			eb =  df / ( df + db );
		}
	}
};
