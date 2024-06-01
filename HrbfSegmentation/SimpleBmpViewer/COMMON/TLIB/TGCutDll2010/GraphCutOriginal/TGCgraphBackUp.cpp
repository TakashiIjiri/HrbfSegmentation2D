#include "StdAfx.h"
#include "TGCgraph.h"
#include "graphcut/adjacency_list/graph.h"

double TGCgraph::m_lambda = 1000;


//wsd_label値は、正整数、連続, 最小 = 1
void TGCgraph::createByWsdImage2D( int width, int height, const vector<int> &wsd_labels, const byte* rgba )
{
	int maxLabel = 0;
	m_WsdLabels.resize( wsd_labels.size() );
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
		m_WsdNodes[ m_WsdLabels[idx] ].addPixel( idx, rgba[ idx*4 + 0 ] / 255.0, 
													  rgba[ idx*4 + 1 ] / 255.0, 
								   				 	  rgba[ idx*4 + 2 ] / 255.0);
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
	int maxLabel = 0;
	m_WsdLabels.resize( wsd_labels.size() );
	
	for( int i= 0; i < (int) m_WsdLabels.size(); ++i ){
		m_WsdLabels[i] = wsd_labels[i] - 1;
		maxLabel    = max( maxLabel, wsd_labels[i] );
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
		m_WsdNodes[ m_WsdLabels[idx] ].addPixel( idx, rgba[ idx*4 + 0 ] / 255.0, 
													  rgba[ idx*4 + 1 ] / 255.0, 
									   				  rgba[ idx*4 + 2 ] / 255.0);

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






void TGCgraph::calcMinimumGraphCut(const set<int> &forePixIds, 
								   const set<int> &backPixIds, vector<int> &isForeBack)
{

	//fore/back pixels--> fore/back nodes
	set< int > foreLs, backLs;
	set<int>::const_iterator it;
	for( it = forePixIds.begin(); it != forePixIds.end(); ++it ) foreLs.insert( m_WsdLabels[ *it ] );//y * width + xが入っている
	for( it = backPixIds.begin(); it != backPixIds.end(); ++it ) backLs.insert( m_WsdLabels[ *it ] );//y * width + xが入っている

	//node --> color 抽出
	vector< colorRGB > foreColors, backColors;
	for( it = foreLs.begin(); it != foreLs.end(); ++it ) foreColors.push_back( m_WsdNodes[ *it].m_color );
	for( it = backLs.begin(); it != backLs.end(); ++it ) backColors.push_back( m_WsdNodes[ *it].m_color );


	//////////////////////////////////////////////////////
	//各nodeにerror E1_fore E1_backを挿入/////////////////
	for( int i= 0; i < (int) m_WsdNodes.size(); ++i){
		TGCnode &p = m_WsdNodes[i];
		calcE1( foreColors, backColors, p.m_color, p.m_E1_f, p.m_E1_b );
	}

	//前景背景にinfを入れる
	for( it = foreLs.begin(); it != foreLs.end(); ++it ){ m_WsdNodes[ *it ].m_E1_f = 10000000000000000.0;
														  m_WsdNodes[ *it ].m_E1_b = 0;                   }
	for( it = backLs.begin(); it != backLs.end(); ++it ){ m_WsdNodes[ *it ].m_E1_f = 0;
														  m_WsdNodes[ *it ].m_E1_b = 10000000000000000.0; }

	Graph::node_id *nodes = new Graph::node_id[ m_WsdNodes.size() ];
	Graph *g = new Graph();

	//nodeを挿入
	for(int i = 0; i < (int)m_WsdNodes.size(); ++i ) nodes[i] = g->add_node();
	for(int i = 0; i < (int)m_WsdNodes.size(); ++i ) g->set_tweights( nodes[i], m_WsdNodes[i].m_E1_f, m_WsdNodes[i].m_E1_b);
	
	//edgeを挿入
	for(int i = 0; i < (int)m_WsdNodes.size(); ++i)
	{
		map<int, double>::const_iterator edgeIt;
		for( edgeIt = m_WsdEdges[i].begin(); edgeIt != m_WsdEdges[i].end(); ++edgeIt )if( i < edgeIt->first )
			g->add_edge( nodes[      i       ], 
				         nodes[ edgeIt->first], edgeIt->second, edgeIt->second );
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












void TGCgraph::calcMinimumGraphCut_PixelLebel2D(int width, int height, byte* rgba,
												  const set<int> &forePixIdxs, 
												  const set<int> &backPixIdxs, vector<int> &isForeBack)
{

	//pixel --> color 抽出
	vector< colorRGB > foreColors(forePixIdxs.size()),
		               backColors(backPixIdxs.size());
	{
		int i; set<int>::const_iterator it;
		for(i=0, it = forePixIdxs.begin(); it != forePixIdxs.end(); ++it,++i ) foreColors[i].SetByByte( &rgba[4*(*it)] );
		for(i=0, it = backPixIdxs.begin(); it != backPixIdxs.end(); ++it,++i ) backColors[i].SetByByte( &rgba[4*(*it)] );
	}

	vector<int> pixId2nodeId( width*height, -1 );//node以外はfore=-1 back=-1, nodeにはnodeIdxが入る
	for( int i = 0; i < (int)m_WsdLabels.size(); ++i)
		pixId2nodeId[ i ] = m_WsdNodes[ m_WsdLabels[i] ].m_isFore ? -1 : -2;

	
	vector< TGCpixelNode* > pixNodes; 
	pixNodes.reserve( width * height / 4 );

	//計算すべきpixelを抽出
	for( int i=0; i<(int) m_WsdEdges.size(); ++i)
	for( map<int,double>::const_iterator it = m_WsdEdges[i].begin(); it != m_WsdEdges[i].end(); ++it)
	{
		TGCnode &ni = m_WsdNodes[ i        ];
		TGCnode &nj = m_WsdNodes[ it->first];
		if( ni.m_isFore != nj.m_isFore )
		{
			//add pixel Nodes 
			for( int kk = 0; kk < (int) ni.m_pixelIDs.size(); ++kk)
			{
				int pixIdx = ni.m_pixelIDs[kk];
				pixId2nodeId[pixIdx] = (int) pixNodes.size();
				pixNodes.push_back( new TGCpixelNode(pixIdx, &rgba[4*pixIdx])  );
			}
			for( int kk = 0; kk < (int)nj.m_pixelIDs.size(); ++kk)
			{
				int pixIdx = nj.m_pixelIDs[kk];
				pixId2nodeId[pixIdx] = (int) pixNodes.size();
				pixNodes.push_back( new TGCpixelNode(pixIdx, &rgba[4*pixIdx])  );
			}
		}
	}


	Graph::node_id *nodes = new Graph::node_id[ pixNodes.size() ];
	Graph *g = new Graph();
	for(int i = 0; i < (int)pixNodes.size(); ++i ) nodes[i] = g->add_node();

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
		if(     pixCondition == 1 ){ E1_f = 10000000000000000.0; E1_b = 0; }
		else if(pixCondition ==-1 ){ E1_b = 10000000000000000.0; E1_f = 0; }
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
	for( int i = 0; i < (int)pixId2nodeId.size(); ++i)
	{
		isForeBack[i] = ( pixId2nodeId[i]==-1) ?  1 :
			            ( pixId2nodeId[i]==-2) ?  0 : 
						( g->what_segment( nodes[ pixId2nodeId[i] ] ) == Graph::SOURCE) ? 1 : 0;
	}

	delete g;
	delete [] nodes;
	for( int i= 0; i < (int) pixNodes.size(); ++i) delete pixNodes[i]; pixNodes.clear();
}



//graphcut segmentation in pixel level
void TGCgraph::calcMinimumGraphCut_PixelLebel3D(int width, 
												int height, 
												int depth, byte* rgba,
												const set<int> &foreVoxIdxs, 
												const set<int> &backVoxIdxs, vector<int> &isForeBack)
{
	//pixel --> color 抽出
	vector< colorRGB > foreColors(foreVoxIdxs.size()),
		               backColors(backVoxIdxs.size());
	{
		int i; set<int>::const_iterator it;
		for(i=0, it = foreVoxIdxs.begin(); it != foreVoxIdxs.end(); ++it,++i ) foreColors[i].SetByByte( &rgba[4*(*it)] );
		for(i=0, it = backVoxIdxs.begin(); it != backVoxIdxs.end(); ++it,++i ) backColors[i].SetByByte( &rgba[4*(*it)] );
	}

	vector<int> voxId2nodeId( width*height*depth, -1 );//node以外はfore=-1 back=-1, nodeにはnodeIdxが入る
	for( int i = 0; i < (int)m_WsdLabels.size(); ++i)
		voxId2nodeId[ i ] = m_WsdNodes[ m_WsdLabels[i] ].m_isFore ? -1 : -2;

	
	vector< TGCpixelNode* > voxNodes; 
	voxNodes.reserve( width * height * depth / 4 );

	//計算すべきpixelを抽出
	for( int i=0; i<(int) m_WsdEdges.size(); ++i)
	for( map<int,double>::const_iterator it = m_WsdEdges[i].begin(); it != m_WsdEdges[i].end(); ++it)
	{
		TGCnode &ni = m_WsdNodes[ i        ];
		TGCnode &nj = m_WsdNodes[ it->first];
		if( ni.m_pixelIDs.size() > 1000 || nj.m_pixelIDs.size() > 1000 ) continue; //あまりに多いノードは確実に背景の可能性あり
		if( ni.m_isFore != nj.m_isFore )
		{
			//add pixel Nodes 
			for( int kk = 0; kk < (int) ni.m_pixelIDs.size(); ++kk)
			{
				int voxIdx = ni.m_pixelIDs[kk];
				voxId2nodeId[voxIdx] = (int) voxNodes.size();
				voxNodes.push_back( new TGCpixelNode(voxIdx, &rgba[4*voxIdx])  );
			}
			for( int kk = 0; kk < (int)nj.m_pixelIDs.size(); ++kk)
			{
				int voxIdx = nj.m_pixelIDs[kk];
				voxId2nodeId[voxIdx] = (int) voxNodes.size();
				voxNodes.push_back( new TGCpixelNode(voxIdx, &rgba[4*voxIdx])  );
			}
		}
	}


	Graph::node_id *nodes = new Graph::node_id[ voxNodes.size() ];
	Graph *g = new Graph();
	for(int i = 0; i < (int)voxNodes.size(); ++i ) nodes[i] = g->add_node();

	fprintf( stderr, "node size == %d  ", (int)voxNodes.size() );


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
		if(     voxCondition == 1 ){ E1_f = 10000000000000000.0; E1_b = 0; }
		else if(voxCondition ==-1 ){ E1_b = 10000000000000000.0; E1_f = 0; }
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
	Graph::flowtype flow = g->maxflow();  fprintf( stderr, "flow = %f\n", flow );

	isForeBack.resize( m_WsdLabels.size() );
	for( int i = 0; i < (int)voxId2nodeId.size(); ++i)
	{
		isForeBack[i] = ( voxId2nodeId[i]==-1) ?  1 :
			            ( voxId2nodeId[i]==-2) ?  0 : 
						( g->what_segment( nodes[ voxId2nodeId[i] ] ) == Graph::SOURCE) ? 1 : 0;
	}

	delete g;
	delete [] nodes;
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
