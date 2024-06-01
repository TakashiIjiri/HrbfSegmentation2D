#include "StdAfx.h"
#include "TGCgraph.h"
#include "graphcut/adjacency_list/graph.h"
#include <list>
double TGCgraph::m_lambda = 1000;
#define FOREBACK_MAX 1000000000.0

//wsd_label値は、正整数、連続, 最小 = 1
void TGCgraph::createByWsdImage2D( int width, int height, const vector<int> &wsd_labels, const byte* rgba )
{
	clear();
	m_WsdLabels.resize( wsd_labels.size() );
	int maxLabel = 0;
	for( int i= 0; i < (int) m_WsdLabels.size(); ++i )
	{
		m_WsdLabels[i] = wsd_labels[i] - 1;
		maxLabel = max( maxLabel, wsd_labels[i] );
	}
	
	fprintf( stderr, "On creating TGraphCut Graph node num = %d\n", maxLabel );

	m_WsdNodes.clear();  m_WsdNodes.resize( maxLabel );
	m_WsdEdges.clear();  m_WsdEdges.resize( maxLabel );
	
	for( int y = 0; y < height; ++y)
	for( int x = 0; x < width ; ++x)
	{
		int idx = y * width + x;

		//pixelを挿入
		m_WsdNodes[ m_WsdLabels[idx] ].addPixel( idx, rgba[ idx*4 + 0 ] , rgba[ idx*4 + 1 ] , rgba[ idx*4 + 2 ] );
		//edgeを挿入. 右と下のpixel
		if( x < width  - 1 ) addEdge( m_WsdLabels[idx], m_WsdLabels[idx + 1    ] );
		if( y < height - 1 ) addEdge( m_WsdLabels[idx], m_WsdLabels[idx + width] );
	}

	
	//各edgeのenergy E2を計算
	for( int idx0 = 0; idx0 < (int) m_WsdEdges.size(); ++idx0 )
	for( map<int,double>::iterator itr = m_WsdEdges[idx0].begin(); itr != m_WsdEdges[idx0].end(); ++itr)
	{
		itr->second = calcE2( m_WsdNodes[ idx0 ], m_WsdNodes[ itr->first ]) ;
	}
}




void TGCgraph::createByWsdImage3D(int width, int height, int depth, const vector<int> &wsd_labels, const byte* rgba )
{
	clear();

	m_WsdLabels.resize( wsd_labels.size() );
	
	int maxLabel = 0;
	for( int i= 0; i < (int) m_WsdLabels.size(); ++i ){
		m_WsdLabels[i] = wsd_labels[i] - 1;
		maxLabel       = max( maxLabel, wsd_labels[i] );
	}
	
	fprintf( stderr, "On creating TGraphCut Graph node num = %d\n", maxLabel );

	m_WsdNodes.clear();  m_WsdNodes.resize( maxLabel );
	m_WsdEdges.clear();  m_WsdEdges.resize( maxLabel );
	
	int idx =0;
	for( int z = 0; z < depth ; ++z)
	for( int y = 0; y < height; ++y)
	for( int x = 0; x < width ; ++x, ++idx)
	{
		//pixelを挿入
		m_WsdNodes[ m_WsdLabels[idx] ].addPixel( idx, rgba[ idx*4 + 0 ], rgba[ idx*4 + 1 ], rgba[ idx*4 + 2 ] );

		//edgeを挿入. 右と下のpixel
		if( x < width  - 1 ) addEdge( m_WsdLabels[idx], m_WsdLabels[idx + 1           ] );
		if( y < height - 1 ) addEdge( m_WsdLabels[idx], m_WsdLabels[idx + width       ] );
		if( z < depth  - 1 ) addEdge( m_WsdLabels[idx], m_WsdLabels[idx + width*height] );
	}
	
	//各edgeのenergy E2を計算
	for( int idx0 = 0; idx0 < (int) m_WsdEdges.size(); ++idx0 )
	for( map<int,double>::iterator itr = m_WsdEdges[idx0].begin(); itr != m_WsdEdges[idx0].end(); ++itr)
	{
		itr->second = calcE2( m_WsdNodes[ idx0 ], m_WsdNodes[ itr->first ]) ;
	}

	fprintf( stderr, "On creating TGraphCut Graph node num = %d\n", maxLabel );
}






void TGCgraph::calcMinCut(const set<int> &forePixIds, 
								   const set<int> &backPixIds, vector<int> &isForeBack)
{

	//fore/back pixels--> fore/back nodes
	set< int > foreLs, backLs;
	set<int>::const_iterator it;
	for( it = forePixIds.begin(); it != forePixIds.end(); ++it ) foreLs.insert( m_WsdLabels[ *it ] );//y * width + xが入っている
	for( it = backPixIds.begin(); it != backPixIds.end(); ++it ) backLs.insert( m_WsdLabels[ *it ] );//y * width + xが入っている

	//node --> color 抽出
	vector< colorRGB > foreColors, backColors;
	foreColors.reserve(foreLs.size());
	backColors.reserve(backLs.size());
	for( it = foreLs.begin(); it != foreLs.end(); ++it ) foreColors.push_back( m_WsdNodes[ *it].m_color );
	for( it = backLs.begin(); it != backLs.end(); ++it ) backColors.push_back( m_WsdNodes[ *it].m_color );


	Graph::node_id *nodes = new Graph::node_id[ m_WsdNodes.size() ];
	Graph *g = new Graph();

	//nodeを挿入
	for(int i = 0; i < (int)m_WsdNodes.size(); ++i )
	{
		nodes[i] = g->add_node();
			
		double e1_f, e1_b;
		calcE1( foreColors, backColors, m_WsdNodes[i].m_color, e1_f, e1_b);
		g->set_tweights( nodes[i], e1_f, e1_b);
	}

	//前景背景にinfを入れる
	for( it = foreLs.begin(); it != foreLs.end(); ++it ){ g->set_tweights( nodes[*it], FOREBACK_MAX, 0); }
	for( it = backLs.begin(); it != backLs.end(); ++it ){ g->set_tweights( nodes[*it], 0, FOREBACK_MAX); }


	//edgeを挿入
	for(int i = 0; i < (int)m_WsdNodes.size(); ++i)
	{
		map<int, double>::const_iterator edgeIt;
		for( edgeIt = m_WsdEdges[i].begin(); edgeIt != m_WsdEdges[i].end(); ++edgeIt )if( i < edgeIt->first )
		{
			//e2計算しなおす
			double e2 = calcE2( m_WsdNodes[ i ], m_WsdNodes[ edgeIt->first ]);
			g->add_edge( nodes[      i       ], 
				         nodes[ edgeIt->first], e2,e2); //edgeIt->second, edgeIt->second );
		}
	}


	//graph cut//
	Graph::flowtype flow = g->maxflow();  fprintf( stderr, "flow = %f\n", flow );
	
	for(int i = 0; i < (int)m_WsdNodes.size(); ++i) 
		m_WsdNodes[i].m_isFore = ( g->what_segment( nodes[i] ) == Graph::SOURCE);

	delete g;
	delete [] nodes;

	//各pixelがforeかbackかを入れる
	isForeBack.resize( m_WsdLabels.size() );
	for( int i = 0; i < (int)m_WsdLabels.size(); ++i)
		isForeBack[ i ] = m_WsdNodes[ m_WsdLabels[i] ].m_isFore ? 1 : 0;
}
















/*
境界からbound width分の帯状の部分についてのみ計算する
*/
void TGCgraph::calcMinCut_PixelLebel2D(int width, int height, byte* rgba,
												  const set<int> &forePixIdxs, 
												  const set<int> &backPixIdxs, vector<int> &isForeBack, int boundWidth)
{
	const int imageSize = width * height;
	
	/////////////計算すべきpixelを抽出//////////////////////////////////////
	int   *pixId2nodeId = new int  [ imageSize ]; //node以外はfore=-1 back=-2, nodeにはnodeIdxが入る
	short *pixelStates  = new short[ imageSize ]; //0: is not in the bound,  1: is in the bound;

	for( int i = 0; i < (int)m_WsdLabels.size(); ++i)
		pixId2nodeId[ i ] = m_WsdNodes[ m_WsdLabels[i] ].m_isFore ? -1 : -2;

	memset( pixelStates, 0, sizeof( short )* imageSize );

	for( int y = 1; y < height; ++y)
	for( int x = 1; x < width ; ++x) //手前との差をみるから1からスタート
	{
		int idx = x + y * width;

		if( pixId2nodeId[ idx ] != pixId2nodeId[ idx-1     ] || 
			pixId2nodeId[ idx ] != pixId2nodeId[ idx-width ])
		{
			//flagを立てる
			for( int y_off = -boundWidth; y_off <= boundWidth; ++y_off)
			for( int x_off = -boundWidth; x_off <= boundWidth; ++x_off)
			if( 0 <= x+x_off && x+x_off < width && 0 <= y+y_off && y+y_off < height)
			{
				pixelStates[idx + x_off + y_off * width] = 1;
			}
		}
	}
	
	//isForeBack.resize( m_WsdLabels.size() );
	//for( int i = 0; i < imageSize; ++i) isForeBack[i] = ( pixelStates[i]==1) ? 1:0; 
	//return;

	////////////計算すべきpixelのindexを計算//////////////////////////////////////
	vector< TGCpixelNode* > pixNodes;  
	pixNodes.reserve( width * height / 4 );

	for( int i = 0; i < imageSize; ++i) if( pixelStates[i] == 1)
	{
		pixId2nodeId[ i ] = (int) pixNodes.size();
		pixNodes.push_back( new TGCpixelNode( i, &rgba[4*i])  );
	}


	////////////////////////graph構築//////////////////////////////////////
	Graph::node_id *nodes = new Graph::node_id[ pixNodes.size() ];
	Graph *g = new Graph();
	for(int i = 0; i < (int)pixNodes.size(); ++i ) nodes[i] = g->add_node();	

	//pixel --> color 抽出
	vector< colorRGB > foreColors(forePixIdxs.size()),
		               backColors(backPixIdxs.size());
	{
		int i; set<int>::const_iterator it;
		for(i=0, it = forePixIdxs.begin(); it != forePixIdxs.end(); ++it,++i ) foreColors[i].Set( &rgba[4*(*it)] );
		for(i=0, it = backPixIdxs.begin(); it != backPixIdxs.end(); ++it,++i ) backColors[i].Set( &rgba[4*(*it)] );
	}


	for( int nodeI = 0; nodeI < (int) pixNodes.size(); ++nodeI)
	{
		int pixIdx = pixNodes[nodeI]->m_pixelID;
		int pixX   = pixIdx  %  width;
		int pixY   = pixIdx  /  width;
		bool edgeToRight = false;
		bool edgeToBelow = false;
		
		int pixCondition = 0;// -1 backpixelが隣接  1 forepixelが隣接
		//上下左右をみてfore/backに隣接していたら、E1f/E1b = 無限
		if( pixX != 0       ) {
			if     ( pixId2nodeId[ pixIdx-1 ] == -1 ) pixCondition =  1;
			else if( pixId2nodeId[ pixIdx-1 ] == -2 ) pixCondition = -1;
		}
		if( pixX != width-1 ){
			if     ( pixId2nodeId[ pixIdx+1 ] == -1 ) pixCondition =  1;
			else if( pixId2nodeId[ pixIdx+1 ] == -2 ) pixCondition = -1;
			else edgeToRight = true;
		}
		if( pixY != 0       ){
			if     ( pixId2nodeId[ pixIdx-width ] == -1 ) pixCondition =  1;
			else if( pixId2nodeId[ pixIdx-width ] == -2 ) pixCondition = -1;
		}
		if( pixY != height-1 ){
			if     ( pixId2nodeId[ pixIdx+width ] == -1 ) pixCondition =  1;
			else if( pixId2nodeId[ pixIdx+width ] == -2 ) pixCondition = -1;
			else edgeToBelow = true;

		}

		double E1_f,E1_b;
		if(     pixCondition == 1 ){ E1_f = FOREBACK_MAX; E1_b = 0; }
		else if(pixCondition ==-1 ){ E1_b = FOREBACK_MAX; E1_f = 0; }
		else {
			calcE1( foreColors, backColors, pixNodes[nodeI]->m_color, E1_f, E1_b);
		}
		g->set_tweights( nodes[nodeI], E1_f, E1_b);
		
		//edgeを追加
		if( edgeToRight ){
			int nodeJ = pixId2nodeId[ pixIdx+1 ];
			double E2 = calcE2( pixNodes[nodeI]->m_color, pixNodes[nodeJ]->m_color ); 
			g->add_edge( nodes[ nodeI ], nodes[ nodeJ ], E2, E2);
		}
		if( edgeToBelow ){
			int nodeJ = pixId2nodeId[ pixIdx+width ];
			double E2 = calcE2( pixNodes[nodeI]->m_color, pixNodes[nodeJ]->m_color ); 
			g->add_edge( nodes[ nodeI ], nodes[ nodeJ ], E2, E2);
		}	
	}
	

	//graph cut//
	fprintf( stderr, "node size == %d  ", (int)pixNodes.size() );
	Graph::flowtype flow = g->maxflow();  fprintf( stderr, "flow = %f\n", flow );

	isForeBack.resize( m_WsdLabels.size() );
	for( int i = 0; i < imageSize; ++i)
	{
		isForeBack[i] = ( pixId2nodeId[i]==-1) ?  1 :
			            ( pixId2nodeId[i]==-2) ?  0 : 
						( g->what_segment( nodes[ pixId2nodeId[i] ] ) == Graph::SOURCE) ? 1 : 0;
	}

	delete g;
	delete[] nodes;
	delete[] pixelStates;
	delete[] pixId2nodeId;
	for( int i= 0; i < (int) pixNodes.size(); ++i) delete pixNodes[i]; pixNodes.clear();
}





/*
境界からbound width分の帯状の部分についてのみ計算する
*/
void TGCgraph::calcMinCut_PixelLebel3D(int width, int height, int depth, byte* rgba,
												  const set<int> &forePixIdxs, 
												  const set<int> &backPixIdxs, vector<int> &isForeBack, int boundWidth)
{
	const int volumeSize = width * height * depth;
	
	/////////////計算すべきpixelを抽出//////////////////////////////////////
	int   *voxId2nodeId = new int  [ volumeSize ]; //node以外はfore=-1 back=-2, nodeにはnodeIdxが入る
	short *voxelStates  = new short[ volumeSize ]; //0: 境界帯外,  1: 境界帯内;

	for( int i = 0; i < volumeSize; ++i) voxId2nodeId[i] =  m_WsdNodes[m_WsdLabels[i]].m_isFore ? -1 : -2;
	memset( voxelStates, 0, sizeof( short )* volumeSize );


	for( int z = 1; z < depth ; ++z)
	for( int y = 1; y < height; ++y)
	for( int x = 1; x < width ; ++x) //手前との差をみるから1からスタート
	{
		int idx = x + y * width + z * width*height;

		if( voxId2nodeId[ idx ] != voxId2nodeId[ idx-1            ] || 
			voxId2nodeId[ idx ] != voxId2nodeId[ idx-width        ] ||
			voxId2nodeId[ idx ] != voxId2nodeId[ idx-width*height ] )
		{
			//flagを立てる
			for( int z_off = -boundWidth; z_off <= boundWidth; ++z_off) if( 0 <= z+z_off && z+z_off < depth  )
			for( int y_off = -boundWidth; y_off <= boundWidth; ++y_off) if( 0 <= y+y_off && y+y_off < height )
			for( int x_off = -boundWidth; x_off <= boundWidth; ++x_off) if( 0 <= x+x_off && x+x_off < width  )
			{
				voxelStates[idx + x_off 
					            + y_off * width 
				                + z_off * width*height] = 1;
			}
		}
	}
	
	//isForeBack.resize( m_WsdLabels.size() );
	//for( int i = 0; i < volumeSize; ++i) isForeBack[i] = ( voxelStates[i]==1) ? 1:0; 
	//return;

	////////////計算すべきpixelのindexを計算//////////////////////////////////////
	vector< TGCpixelNode* > voxNodes;  
	voxNodes.reserve( volumeSize / 4 );

	for( int i = 0; i < volumeSize; ++i) if( voxelStates[i] == 1)
	{
		voxId2nodeId[ i ] = (int) voxNodes.size();
		voxNodes.push_back( new TGCpixelNode( i, &rgba[4*i])  );
	}


	////////////////////////graph構築//////////////////////////////////////
	Graph::node_id *nodes = new Graph::node_id[ voxNodes.size() ];
	Graph *g = new Graph();
	for(int i = 0; i < (int)voxNodes.size(); ++i ) nodes[i] = g->add_node();	

	//pixel --> color 抽出
	vector< colorRGB > foreColors(forePixIdxs.size()),
		               backColors(backPixIdxs.size());
	{
		int i; set<int>::const_iterator it;
		for(i=0, it = forePixIdxs.begin(); it != forePixIdxs.end(); ++it,++i ) foreColors[i].Set( &rgba[4*(*it)] );
		for(i=0, it = backPixIdxs.begin(); it != backPixIdxs.end(); ++it,++i ) backColors[i].Set( &rgba[4*(*it)] );
	}


	for( int nodeI = 0; nodeI < (int) voxNodes.size(); ++nodeI)
	{
		int voxIdx = voxNodes[nodeI]->m_pixelID;
		int voxZ   =  voxIdx         /  (width * height);
		int voxY   = (voxIdx - voxZ) /  width;
		int voxX   =  voxIdx - voxZ - voxY;
		bool edgeToNexX = false;
		bool edgeToNexY = false;
		bool edgeToNexZ = false;
		
		int voxCondition = 0;// -1 backpixelが隣接  1 forepixelが隣接
		//左右前後上下をみてfore/backに隣接していたら、E1f/E1b = 無限
		if( voxX != 0       ) {
			if     ( voxId2nodeId[ voxIdx-1		] == -1 ) voxCondition =  1;
			else if( voxId2nodeId[ voxIdx-1		] == -2 ) voxCondition = -1;
		}
		if( voxX != width-1 ){
			if     ( voxId2nodeId[ voxIdx+1		] == -1 ) voxCondition =  1;
			else if( voxId2nodeId[ voxIdx+1		] == -2 ) voxCondition = -1;
			else edgeToNexX = true;
		}
		if( voxY != 0       ){
			if     ( voxId2nodeId[ voxIdx-width ] == -1 ) voxCondition =  1;
			else if( voxId2nodeId[ voxIdx-width ] == -2 ) voxCondition = -1;
		}
		if( voxY != height-1 ){
			if     ( voxId2nodeId[ voxIdx+width ] == -1 ) voxCondition =  1;
			else if( voxId2nodeId[ voxIdx+width ] == -2 ) voxCondition = -1;
			else edgeToNexY = true;
		}
		if( voxZ != 0       ){
			if     ( voxId2nodeId[ voxIdx-width*height ] == -1 ) voxCondition =  1;
			else if( voxId2nodeId[ voxIdx-width*height ] == -2 ) voxCondition = -1;
		}
		if( voxZ != depth-1 ){
			if     ( voxId2nodeId[ voxIdx+width*height ] == -1 ) voxCondition =  1;
			else if( voxId2nodeId[ voxIdx+width*height ] == -2 ) voxCondition = -1;
			else edgeToNexZ = true;
		}

		double E1_f,E1_b;
		if(     voxCondition == 1 ){ E1_f = FOREBACK_MAX; E1_b = 0; }
		else if(voxCondition ==-1 ){ E1_b = FOREBACK_MAX; E1_f = 0; }
		else {
			calcE1( foreColors, backColors, voxNodes[nodeI]->m_color, E1_f, E1_b);
		}
		g->set_tweights( nodes[nodeI], E1_f, E1_b);
		
		//edgeを追加
		if( edgeToNexX){
			int nodeJ = voxId2nodeId[ voxIdx+1 ];
			double E2 = calcE2( voxNodes[nodeI]->m_color, voxNodes[nodeJ]->m_color ); 
			g->add_edge( nodes[ nodeI ], nodes[ nodeJ ], E2, E2);
		}
		if( edgeToNexY ){
			int nodeJ = voxId2nodeId[ voxIdx+width ];
			double E2 = calcE2( voxNodes[nodeI]->m_color, voxNodes[nodeJ]->m_color ); 
			g->add_edge( nodes[ nodeI ], nodes[ nodeJ ], E2, E2);
		}
		if( edgeToNexZ ){
			int nodeJ = voxId2nodeId[ voxIdx+width*height ];
			double E2 = calcE2( voxNodes[nodeI]->m_color, voxNodes[nodeJ]->m_color ); 
			g->add_edge( nodes[ nodeI ], nodes[ nodeJ ], E2, E2);
		}	
	}
	
	
	//graph cut//
	fprintf( stderr, "node size == %d  ", (int)voxNodes.size() );
	Graph::flowtype flow = g->maxflow();  fprintf( stderr, "flow = %f\n", flow );

	isForeBack.resize( m_WsdLabels.size() );
	for( int i = 0; i < volumeSize; ++i)
	{
		isForeBack[i] = ( voxId2nodeId[i]==-1) ?  1 :
			            ( voxId2nodeId[i]==-2) ?  0 : 
						( g->what_segment( nodes[ voxId2nodeId[i] ] ) == Graph::SOURCE) ? 1 : 0;
	}


	delete g;
	delete[] nodes;
	delete[] voxelStates;
	delete[] voxId2nodeId;
	for( int i= 0; i < (int) voxNodes.size(); ++i) delete voxNodes[i]; voxNodes.clear();
}





//制約点による、色の類似度を考慮しない
void TGCgraph::calcMinCur_PixelLevel3D_onlyBoundary( int width, int height, int depth, const byte *rgba, short* volumeInOut, int boundWidth)
{
	const int volumeSize = width * height * depth;
	const int widhei     = width * height;
	/////////////計算すべきpixelを抽出////////////////////////////////////////////////////////////////////
	int   *voxId2nodeId = new int  [ volumeSize ]; //node以外はfore=-1 back=-2, nodeにはnodeIdxが入る
	for( int i = 0; i < volumeSize; ++i) voxId2nodeId[i] =  volumeInOut[i] ? -1 : -2;//fore : -1  back : -2

	short *voxelStates  = new short[ volumeSize ]; //0: 境界帯外,  1: 境界帯内;
	memset( voxelStates, 0, sizeof( short )* volumeSize );


	for( int z = 1; z < depth ; ++z)
	for( int y = 1; y < height; ++y)
	for( int x = 1; x < width ; ++x) //手前との差をみるから1からスタート
	{
		int idx = x + y * width + z * width*height;

		if( voxId2nodeId[ idx ] != voxId2nodeId[ idx-1      ] || 
			voxId2nodeId[ idx ] != voxId2nodeId[ idx-width  ] ||
			voxId2nodeId[ idx ] != voxId2nodeId[ idx-widhei ] )
		{
			//flagを立てる
			for( int z_off = -boundWidth; z_off <= boundWidth; ++z_off) if( 0 <= z+z_off && z+z_off < depth  )
			for( int y_off = -boundWidth; y_off <= boundWidth; ++y_off) if( 0 <= y+y_off && y+y_off < height )
			for( int x_off = -boundWidth; x_off <= boundWidth; ++x_off) if( 0 <= x+x_off && x+x_off < width  )
			{
				voxelStates[idx + x_off 
					            + y_off * width 
				                + z_off * widhei ] = 1;
			}
		}
	}
	
	////////////計算すべきpixelのindexを計算//////////////////////////////////////
	vector< TGCpixelNode* > voxNodes;  
	voxNodes.reserve( volumeSize / 4 );

	for( int i = 0; i < volumeSize; ++i) if( voxelStates[i] == 1){
		voxId2nodeId[ i ] = (int) voxNodes.size();
		voxNodes.push_back( new TGCpixelNode( i, &rgba[4*i])  );
	}

	////////////////////////graph構築//////////////////////////////////////
	Graph::node_id *nodes = new Graph::node_id[ voxNodes.size() ];
	Graph *g = new Graph();
	for(int i = 0; i < (int)voxNodes.size(); ++i ) nodes[i] = g->add_node();	

	for( int nodeI = 0; nodeI < (int) voxNodes.size(); ++nodeI)
	{
		int voxIdx = voxNodes[nodeI]->m_pixelID;
		int voxZ   =  voxIdx         /  (width * height);
		int voxY   = (voxIdx - voxZ) /  width;
		int voxX   =  voxIdx - voxZ - voxY;
		bool edgeToNexX = false;
		bool edgeToNexY = false;
		bool edgeToNexZ = false;
		
		int voxCondition = 0;// -1 backpixelが隣接  1 forepixelが隣接
		//左右前後上下をみてfore/backに隣接していたら、E1f/E1b = 無限
		if( voxX != 0       ) {
			if     ( voxId2nodeId[ voxIdx-1		] == -1 ) voxCondition =  1;
			else if( voxId2nodeId[ voxIdx-1		] == -2 ) voxCondition = -1;
		}
		if( voxX != width-1 ){
			if     ( voxId2nodeId[ voxIdx+1		] == -1 ) voxCondition =  1;
			else if( voxId2nodeId[ voxIdx+1		] == -2 ) voxCondition = -1;
			else edgeToNexX = true;
		}
		if( voxY != 0       ){
			if     ( voxId2nodeId[ voxIdx-width ] == -1 ) voxCondition =  1;
			else if( voxId2nodeId[ voxIdx-width ] == -2 ) voxCondition = -1;
		}
		if( voxY != height-1 ){
			if     ( voxId2nodeId[ voxIdx+width ] == -1 ) voxCondition =  1;
			else if( voxId2nodeId[ voxIdx+width ] == -2 ) voxCondition = -1;
			else edgeToNexY = true;
		}
		if( voxZ != 0       ){
			if     ( voxId2nodeId[ voxIdx-width*height ] == -1 ) voxCondition =  1;
			else if( voxId2nodeId[ voxIdx-width*height ] == -2 ) voxCondition = -1;
		}
		if( voxZ != depth-1 ){
			if     ( voxId2nodeId[ voxIdx+width*height ] == -1 ) voxCondition =  1;
			else if( voxId2nodeId[ voxIdx+width*height ] == -2 ) voxCondition = -1;
			else edgeToNexZ = true;
		}

		double E1_f,E1_b;
		if(     voxCondition == 1 ){ E1_f = FOREBACK_MAX; E1_b = 0; }
		else if(voxCondition ==-1 ){ E1_b = FOREBACK_MAX; E1_f = 0; }
		else {                       E1_f = E1_b = 1;               }

		g->set_tweights( nodes[nodeI], E1_f, E1_b);
		
		//edgeを追加
		if( edgeToNexX){
			int nodeJ = voxId2nodeId[ voxIdx+1 ];
			double E2 = calcE2( voxNodes[nodeI]->m_color, voxNodes[nodeJ]->m_color ); 
			g->add_edge( nodes[ nodeI ], nodes[ nodeJ ], E2, E2);
		}
		if( edgeToNexY ){
			int nodeJ = voxId2nodeId[ voxIdx+width ];
			double E2 = calcE2( voxNodes[nodeI]->m_color, voxNodes[nodeJ]->m_color ); 
			g->add_edge( nodes[ nodeI ], nodes[ nodeJ ], E2, E2);
		}
		if( edgeToNexZ ){
			int nodeJ = voxId2nodeId[ voxIdx+width*height ];
			double E2 = calcE2( voxNodes[nodeI]->m_color, voxNodes[nodeJ]->m_color ); 
			g->add_edge( nodes[ nodeI ], nodes[ nodeJ ], E2, E2);
		}	
	}
	
	
	//graph cut//
	fprintf( stderr, "node size == %d  ", (int)voxNodes.size() );
	Graph::flowtype flow = g->maxflow();  fprintf( stderr, "flow = %f\n", flow );

	for( int i = 0; i < volumeSize; ++i)
	{
		volumeInOut[i] = ( voxId2nodeId[i]==-1) ?  1 :
			             ( voxId2nodeId[i]==-2) ?  0 : 
						  ( g->what_segment( nodes[ voxId2nodeId[i] ] ) == Graph::SOURCE) ? 1 : 0;
	}


	delete g;
	delete[] nodes;
	delete[] voxelStates;
	delete[] voxId2nodeId;
	for( int i= 0; i < (int) voxNodes.size(); ++i) delete voxNodes[i]; voxNodes.clear();
}








void TGCgraph::testMinCut() const
{
	Graph::node_id nodes[2];
	Graph *g = new Graph();

	nodes[0] = g -> add_node();
	nodes[1] = g -> add_node();
	g -> set_tweights(nodes[0], 5, 2);
	g -> set_tweights(nodes[1], 2, 4.2);
	g -> add_edge(nodes[0], nodes[1], 0.5, 0.5);

	Graph::flowtype flow = g -> maxflow();

	fprintf(stderr, "Flow = %f\n", flow);
	fprintf(stderr, "Minimum cut:\n");
	if (g->what_segment(nodes[0]) == Graph::SOURCE)
		fprintf(stderr, "node0 is in the SOURCE set\n");
	else
		fprintf(stderr, "node0 is in the SINK set\n");


	if (g->what_segment(nodes[1]) == Graph::SOURCE)
		fprintf(stderr, "node1 is in the SOURCE set\n");
	else
		fprintf(stderr, "node1 is in the SINK set\n");


	delete g;
}



//////////////////////////////////////////////////////////////////////////////////////////
/////GRAPH CUT PARTIAL///////////////////////////////////////////////////////////////////

/*
含んでいるすべてのpixelのalpha値が0 のwsd nodeをdisableする
*/
void TGCgraph_partial::updateEnableWsdNode( byte* rgba)
{
	if( m_WsdNodes.size() == 0 ){fprintf( stderr ,"wsd node should be generated before calling this !\n"); return;}

	for( int i=0; i<(int) m_WsdNodes.size(); ++i)
	{
		TGCnode &n = m_WsdNodes[i];
		n.m_isEnable = false;
		for( int k = 0; k < (int) n.m_pixelIDs.size(); ++k) 
			if( rgba[ 4 * n.m_pixelIDs[k] + 3 ] != 0 )
			{
				n.m_isEnable = true;
				break;
			}
	}
}






void TGCgraph::preCheckAndDivideInconsistentNodes(const set<int> &forePixelIdxs, const set<int> &backPixelIdxs, int sizeW, int sizeH, int sizeD, const byte *rgba )
{
	for( set<int>::const_iterator fit = forePixelIdxs.begin(); fit != forePixelIdxs.end(); ++fit) 
	for( set<int>::const_iterator bit = backPixelIdxs.begin(); bit != backPixelIdxs.end(); ++bit)
	{
		if( m_WsdLabels[*fit] == m_WsdLabels[*bit] && *fit != *bit )
		{
			divideInconsistentNodes( *fit, *bit, sizeW, sizeH, sizeD, rgba ); 
		}
	}
}

void calcVoxelXYZfromIdx( const int idx, 
								 const int W, 
								 const int H, 
								 const int D, 
								 int &x, int &y, int &z)
{
	z =     idx       / (W * H);
	y = (idx - z*W*H) / W;
	x =  idx - z*W*H - y*W;
}
/*
変更する物
vector<map<int, double>>  m_WsdEdges ;//edge i jがある時には, m_edges[i][a].first = j, m_edges[i][?].second = error = E2 
vector<TGCnode         >  m_WsdNodes ;//node (region) idx --> pixel         idx  ( m_pixelIDs )
vector<int             >  m_WsdLabels;//pixel         idx --> node (region) idx 
*/
void TGCgraph::divideInconsistentNodes(int foreVoxIdx, int backVoxIdx, const int W, const int H, const int D, const byte *rgba)
{

	const int nId = m_WsdLabels[ foreVoxIdx ];
	const int WH = W*H;

	list<int> fBoundVoxIds; set<int> fVox;  
	list<int> bBoundVoxIds; set<int> bVox;  
	fBoundVoxIds.push_back( foreVoxIdx ); fVox.insert(foreVoxIdx );
	bBoundVoxIds.push_back( backVoxIdx ); bVox.insert(backVoxIdx );

	int minX= 100000, minY= 100000, minZ= 100000;
	int maxX=-100000, maxY=-100000, maxZ=-100000;
	//growth
	while( !bBoundVoxIds.empty() || !fBoundVoxIds.empty())
	{
		//fore成長
		int x,y,z;

		if( !bBoundVoxIds.empty() )
		{
			int idx = bBoundVoxIds.front(); bBoundVoxIds.pop_front( );
			calcVoxelXYZfromIdx( idx, W,H,D, x, y, z);

			if( x>0  && m_WsdLabels[idx-1 ]==nId && fVox.find(idx-1 )==fVox.end()&&bVox.find(idx-1 )==bVox.end()){ bVox.insert(idx-1 ); bBoundVoxIds.push_back(idx-1 ); }
			if( x<W-1&& m_WsdLabels[idx+1 ]==nId && fVox.find(idx+1 )==fVox.end()&&bVox.find(idx+1 )==bVox.end()){ bVox.insert(idx+1 ); bBoundVoxIds.push_back(idx+1 ); }
			if( y>0  && m_WsdLabels[idx-W ]==nId && fVox.find(idx-W )==fVox.end()&&bVox.find(idx-W )==bVox.end()){ bVox.insert(idx-W ); bBoundVoxIds.push_back(idx-W ); }
			if( y<H-1&& m_WsdLabels[idx+W ]==nId && fVox.find(idx+W )==fVox.end()&&bVox.find(idx+W )==bVox.end()){ bVox.insert(idx+W ); bBoundVoxIds.push_back(idx+W ); }
			if( z>0  && m_WsdLabels[idx-WH]==nId && fVox.find(idx-WH)==fVox.end()&&bVox.find(idx-WH)==bVox.end()){ bVox.insert(idx-WH); bBoundVoxIds.push_back(idx-WH); }
			if( z<D-1&& m_WsdLabels[idx+WH]==nId && fVox.find(idx+WH)==fVox.end()&&bVox.find(idx+WH)==bVox.end()){ bVox.insert(idx+WH); bBoundVoxIds.push_back(idx+WH); }
			minX=min(minX,x); minY=min(minY,y); minZ=min(minZ,z);
			maxX=max(maxX,x); maxY=max(maxY,y); maxZ=max(maxZ,z);
		}

		if( !fBoundVoxIds.empty() )
		{
			int idx = fBoundVoxIds.front(); fBoundVoxIds.pop_front();
			calcVoxelXYZfromIdx( idx, W,H,D, x, y, z);

			if( x>0  && m_WsdLabels[idx-1 ]==nId && fVox.find(idx-1 )==fVox.end()&&bVox.find(idx-1 )==bVox.end()){ fVox.insert(idx-1 ); fBoundVoxIds.push_back(idx-1 ); }
			if( x<W-1&& m_WsdLabels[idx+1 ]==nId && fVox.find(idx+1 )==fVox.end()&&bVox.find(idx+1 )==bVox.end()){ fVox.insert(idx+1 ); fBoundVoxIds.push_back(idx+1 ); }
			if( y>0  && m_WsdLabels[idx-W ]==nId && fVox.find(idx-W )==fVox.end()&&bVox.find(idx-W )==bVox.end()){ fVox.insert(idx-W ); fBoundVoxIds.push_back(idx-W ); }
			if( y<H-1&& m_WsdLabels[idx+W ]==nId && fVox.find(idx+W )==fVox.end()&&bVox.find(idx+W )==bVox.end()){ fVox.insert(idx+W ); fBoundVoxIds.push_back(idx+W ); }
			if( z>0  && m_WsdLabels[idx-WH]==nId && fVox.find(idx-WH)==fVox.end()&&bVox.find(idx-WH)==bVox.end()){ fVox.insert(idx-WH); fBoundVoxIds.push_back(idx-WH); }
			if( z<D-1&& m_WsdLabels[idx+WH]==nId && fVox.find(idx+WH)==fVox.end()&&bVox.find(idx+WH)==bVox.end()){ fVox.insert(idx+WH); fBoundVoxIds.push_back(idx+WH); }
			minX=min(minX,x); minY=min(minY,y); minZ=min(minZ,z);
			maxX=max(maxX,x); maxY=max(maxY,y); maxZ=max(maxZ,z);
		}
	}

	//back//
	m_WsdNodes.push_back( TGCnode() );
	m_WsdEdges.push_back( map<int,double>() );

	TGCnode &backN = m_WsdNodes.back();
	m_WsdNodes[nId]= TGCnode();
	TGCnode &foreN = m_WsdNodes[nId] ;

	const int fNodeId = nId;
	const int bNodeId = (int)m_WsdNodes.size()-1;


	for( set<int>::iterator it = fVox.begin(); it != fVox.end(); ++it){
		m_WsdLabels[ *it ] = fNodeId;
		m_WsdNodes[ fNodeId ].addPixel( (*it),	rgba[ (*it)*4 + 0 ], rgba[ (*it)*4 + 1 ], rgba[ (*it)*4 + 2 ] );
	}
	for( set<int>::iterator it = bVox.begin(); it != bVox.end(); ++it){
		m_WsdLabels[ *it ] = bNodeId;
		m_WsdNodes[ bNodeId ].addPixel( (*it),	rgba[ (*it)*4 + 0 ],  rgba[ (*it)*4 + 1 ], rgba[ (*it)*4 + 2 ] );
	}

	//edgeの付け替えを行う//
	
	//i)nId関係のedgeを削除
	m_WsdEdges[nId].clear();
	for( int i=0; i< (int)m_WsdNodes.size(); ++i) if( i != nId )
	{
		map<int,double>::iterator it = m_WsdEdges[i].find( nId );
		if( it != m_WsdEdges[i].end()) m_WsdEdges[i].erase( it );
	}

	//ii)新たにedgeを張る
	int idx =0;
	for( int z = minZ; z < maxZ; ++z)
	for( int y = minY; y < maxY; ++y)
	for( int x = minX; x < maxX; ++x)
	{
		int idx = x + y*W + z*WH;
		//edgeを挿入. 右と下のpixel
		if( x < W - 1 ) addEdge( m_WsdLabels[idx], m_WsdLabels[idx + 1  ] );
		if( y < H - 1 ) addEdge( m_WsdLabels[idx], m_WsdLabels[idx + W  ] );
		if( z < D - 1 ) addEdge( m_WsdLabels[idx], m_WsdLabels[idx + WH ] );
	}
}






void TGCgraph_partial::updateEnableWsdNode(byte* rgba, int sizeW, int sizeH, int sizeD, 
										   int minX, int minY, int minZ, int maxX, int maxY, int maxZ)

{
	if( m_WsdNodes.size() == 0 ){fprintf( stderr ,"wsd node should be generated before calling this !\n"); return;}


	int counter = 0;
	for( int i=0; i<(int) m_WsdNodes.size(); ++i) m_WsdNodes[i].m_isEnable = false;

	for( int z = 0; z<sizeD; ++z) if( minZ <= z && z <= maxZ)
	for( int y = 0; y<sizeH; ++y) if( minY <= y && y <= maxY)
	for( int x = 0; x<sizeW; ++x) if( minX <= x && x <= maxX)
	{
		int idx = x + y * sizeW + z*sizeW*sizeH;
		if( rgba[ 4 * idx+ 3 ] != 0 ) {m_WsdNodes[ m_WsdLabels[idx] ].m_isEnable = true; }
	}

	//fprintf( stderr, "%d %d %d %d %d %d\n", minX, minY, minZ, maxX, maxY, maxZ);
}




void TGCgraph_partial::calcMinCut_part(const set<int> &forePixIds, 
									   const set<int> &backPixIds, 
									   vector   <int> &isForeBack, int width, int height, int depth, const byte *rgba)
{
	preCheckAndDivideInconsistentNodes( forePixIds, backPixIds, width, height, depth, rgba);

	//fore/back pixels--> fore/back nodes
	set< int > foreLs, backLs;
	set<int>::const_iterator it;
	for( it = forePixIds.begin(); it != forePixIds.end(); ++it ) foreLs.insert( m_WsdLabels[ *it ] );
	for( it = backPixIds.begin(); it != backPixIds.end(); ++it ) backLs.insert( m_WsdLabels[ *it ] );

/*
	fprintf( stderr, "\n\n  %d \n", m_WsdNodes.size() );
	for( int i=0; i<m_WsdNodes.size(); ++i) fprintf( stderr, "%d %d \n", i, m_WsdNodes[i].m_pixelIDs.size() );
	for( it = foreLs.begin(); it != foreLs.end(); ++it ) fprintf( stderr, "fore %d\n", *it);
	for( it = backLs.begin(); it != backLs.end(); ++it ) fprintf( stderr, "back %d\n", *it);
*/




	//node --> color 抽出
	vector< colorRGB > foreColors, backColors; 
	foreColors.reserve(foreLs.size());
	backColors.reserve(backLs.size());
	for( it = foreLs.begin(); it != foreLs.end(); ++it ) foreColors.push_back( m_WsdNodes[ *it].m_color );
	for( it = backLs.begin(); it != backLs.end(); ++it ) backColors.push_back( m_WsdNodes[ *it].m_color );


	Graph *g = new Graph();
	Graph::node_id *nodes = new Graph::node_id[ m_WsdNodes.size() ];

	//以下は node.m_isEnable = falseのnodeも挿入はするが, edgeはつながない実装

	//nodeを挿入
	for(int i = 0; i < (int)m_WsdNodes.size(); ++i )
	{
		nodes[i] = g->add_node();
		
		if( m_WsdNodes[i].m_isEnable ){
			double e1_f, e1_b;
			calcE1( foreColors, backColors, m_WsdNodes[i].m_color, e1_f, e1_b);
			g->set_tweights( nodes[i], e1_f, e1_b);
		}
		else
			g->set_tweights( nodes[i], 0, FOREBACK_MAX);
	}

	//前景背景にinfを入れる (EnableのNodeのみ)
	for( it=foreLs.begin(); it != foreLs.end(); ++it) if(m_WsdNodes[*it].m_isEnable){ g->set_tweights( nodes[*it], FOREBACK_MAX, 0); }
	for( it=backLs.begin(); it != backLs.end(); ++it) if(m_WsdNodes[*it].m_isEnable){ g->set_tweights( nodes[*it], 0, FOREBACK_MAX); }

	//edgeを挿入(EnableのNodeどうしのみ)
	for(int i = 0; i < (int)m_WsdNodes.size(); ++i) if(m_WsdNodes[i].m_isEnable)
	{
		map<int, double>::const_iterator edgeIt;
		for( edgeIt = m_WsdEdges[i].begin(); edgeIt != m_WsdEdges[i].end(); ++edgeIt )
		if( i < edgeIt->first && m_WsdNodes[ edgeIt->first ].m_isEnable )
		{
			//e2計算しなおす
			double e2 = calcE2( m_WsdNodes[ i ], m_WsdNodes[ edgeIt->first ]);
			g->add_edge( nodes[      i       ], 
				         nodes[ edgeIt->first], e2,e2); //edgeIt->second, edgeIt->second );
		}
	}

	//graph cut//
	Graph::flowtype flow = g->maxflow();  fprintf( stderr, "flow = %f\n", flow );
	
	for(int i = 0; i < (int)m_WsdNodes.size(); ++i) 
		m_WsdNodes[i].m_isFore = ( g->what_segment( nodes[i] ) == Graph::SOURCE);


	//各pixelがforeかbackかを入れる
	isForeBack.resize( m_WsdLabels.size() );
	for( int i = 0; i < (int)m_WsdLabels.size(); ++i) 
		isForeBack[ i ] = m_WsdNodes[ m_WsdLabels[i] ].m_isFore ? 1 : 0;

	delete g;
	delete [] nodes;
}


void TGCgraph_partial::calcMinCut_PixelLebel3D_part(int width, 
													int height, 
													int depth, 
													byte* rgba,
													const set<int> &foreVoxIdxs  , 
													const set<int> &backVoxIdxs  , 
													vector   <int> &isForeBack   ,
													int boundWidth)
{





}











