#pragma once
#include "math.h"
#include <vector>
#include <map>
#include <set>
#include <algorithm>

using namespace std;
#ifndef byte
typedef unsigned char byte;
#endif

struct TGCColor
{
	double r, g, b;
	TGCColor(){r=g=b=0;}
	void Set( byte _r, byte _g, byte _b)
	{
		r = _r;
		g = _g;
		b = _b;
	}
};



/*------------------------------------------------------------
2010/10/23 coded by takashi ijiri (RIKEN)
2010/9/25  revised by Takashi Ijiri @ RIKEN
Graph Cut Image segmentation Algorithm

this class supports 
 - graph cut segmentation 
 - for volumetric 
 - single/double/triple channel image
 - 
------------------------------------------------------------*/

//Graph Cut node   contains multiple pixels  &  have average color
class TGCnode
{
public:
	vector<int> m_pixelIDs ;
	TGCColor    m_total    ;
	TGCColor    m_c        ;
	bool        m_isFore   ;
	bool        m_isEnable ;//segmentationÇÃëŒè€Ç©Ç«Ç§Ç©
	
	~TGCnode(){}
	TGCnode (){
		m_isFore   = false;
		m_isEnable = true ;
		m_pixelIDs.clear();
	}
	TGCnode(const TGCnode &n){
		m_total    = n.m_total;
		m_c        = n.m_c;
		m_isFore   = n.m_isFore  ;
		m_isEnable = n.m_isEnable;
		m_pixelIDs = n.m_pixelIDs;
	}
	inline void addPixel(int pixelIdx, byte r, byte g, byte b){
		m_pixelIDs.push_back( pixelIdx );
		m_total.r += r;   m_c.r = m_total.r / (double) m_pixelIDs.size();
		m_total.g += g;   m_c.g = m_total.g / (double) m_pixelIDs.size();
		m_total.b += b;   m_c.b = m_total.b / (double) m_pixelIDs.size();
	}
};


/*----------------------------------------------------------------------------------------------------------------
How to use...
 1. createByWsdImage2D()     : initialize graph shape from watershad segmentation results (call only onece)
 2. calcMinCut()             : compute graph cut segmentation in the watershad level
 3. calcMinCut_PixelLebel2D(): compute graph cut segmentation in the pixel level. This computes graph cut only on thin belt shaped boundary region of the results of 2
----------------------------------------------------------------------------------------------------------------*/
class TGCgraph3D
{
public:
	static double    m_lambda   ;
	vector<int     > m_pixNodeID;//pixel idx --> node (region) idx
	vector<set<int>> m_n_edges  ;//node  idx --> edge starting from the node
	vector<TGCnode > m_n_nodes  ;//node  idx --> node info (idxs of pixels / colors);

private:

public:
	~TGCgraph3D(void){}
	TGCgraph3D (void){}
	inline void clear(){m_pixNodeID.clear(); m_n_edges.clear(); m_n_nodes.clear();}

	//////Initialization and preparation/////////////////////////////////////////////////////////////////////////////////////////////
	//label ID is  sequential integer start from "1". i.e. 1,2,3,4... imgG, imgB can be "0" if the image has smaller channels////////
	void createByWsdImage   (int W, int H, int D, const byte* imgR, const byte* imgG, const byte* imgB, const vector<int> &wsd_labels ); 
	void updateEnableWsdNode(byte* volInOut);

	//imgG, imgB can be "0" if the image has smaller channels!!!
	void calcMinCut(              int W, int H, int D, const byte *imgR, const byte *imgG, const byte *imgB, const set<int> &forePixIds,
															                                                 const set<int> &backPixIds, 
																												vector<int> &isForeBack);
	void calcMinCut_voxelLv_bound(int W, int H, int D, const byte *imgR, const byte *imgG, const byte *imgB, byte *bVoxInOut, int boundWidth);

private:
	inline void addEdge( int nodeIdx1, int nodeIdx2)
	{
		if( nodeIdx1 == nodeIdx2 ) return;
		if( nodeIdx1 >  nodeIdx2 ) std::swap( nodeIdx1, nodeIdx2 );
		m_n_edges[nodeIdx1].insert( nodeIdx2 );
	}

	//imgR, imgG, imgB can be "0" if the image has smaller channels
	void preCheck_ConflictWsdNodes(const int W, const int H, const int D, const byte *imgR, const byte *imgG, const byte *imgB, const set<int> &forePixelIdxs,  const set<int> &backPixelIdxs);
	void divideConflictWsdNodes   (const int W, const int H, const int D, const byte *imgR, const byte *imgG, const byte *imgB, int foreVoxIdx, int backVoxIdx);
};



class TGCgraph2D
{
public:
	static double m_lambda;
	vector<int     > m_pixNodeID;//pixel idx --> node (region) idx
	vector<set<int>> m_n_edges  ;//node  idx --> edge starting from the node
	vector<TGCnode > m_n_nodes  ;//node  idx --> node info (idxs of pixels / colors);

public:

	~TGCgraph2D(void){}
	TGCgraph2D (void){}
	inline void clear(){m_pixNodeID.clear(); m_n_edges.clear(); m_n_nodes.clear();}

	//////Initialization and preparation/////////////////////////////////////////////////////////////////////////////////////////////
	//label ID is  sequential integer start from "1". i.e. 1,2,3,4... imgG, imgB can be "0" if the image has smaller channels////////
	void createByWsdImage      (int W, int H,        const byte* rgba,  const vector<int> &wsd_labels ); 
	void calcMinCut            (const set<int> &forePixIds, const set<int> &backPixIds, vector<int> &isForeBack);
	void calcMinCut_pixLv_bound(int W, int H, byte* rgba,
										const set<int> &forePixIdxs, 
										const set<int> &backPixIdxs, vector<int> &isForeBack, int boundWidth);

private:
	inline void addEdge( int nodeIdx1, int nodeIdx2){
		if( nodeIdx1 == nodeIdx2 ) return;
		if( nodeIdx1 >  nodeIdx2 ) std::swap( nodeIdx1, nodeIdx2 );
		m_n_edges[nodeIdx1].insert( nodeIdx2 );
	}
};












