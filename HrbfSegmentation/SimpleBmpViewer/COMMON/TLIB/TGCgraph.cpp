#include "StdAfx.h"
#include "TGCGraph.h"

#ifndef TGCUT_STATIC_LINK
#define TGCUT_STATIC_LINK
#endif 

#ifndef FOREBACK_MAX 
#define FOREBACK_MAX 1000000000.0
#endif 

#include "./TGCutDll2010/TGcutStdHead.h" //include after defining the "TGCUT_STATIC_LINK"必要あり
#include <list>


double TGCgraph3D::m_lambda = 1000;
double TGCgraph2D::m_lambda = 1000;



static inline void calcVoxelXYZfromIdx( const int &idx, const int &W, const int &H, const int &D, int &x, int &y, int &z)
{
	z =     idx       / (W * H);
	y = (idx - z*W*H) / W;
	x =  idx - z*W*H - y*W;
}


static inline double distSq( const TGCColor &c1, const TGCColor &c2){
	return (c1.r-c2.r) * (c1.r-c2.r) + 
		   (c1.g-c2.g) * (c1.g-c2.g) + 
		   (c1.b-c2.b) * (c1.b-c2.b);
}

static inline double calcE2( const TGCColor &c0, const TGCColor &c1, double lambda )
{
	return  lambda / ( 1.0 +  distSq( c0, c1) );
}

static inline void calcE1(  const vector<TGCColor> &f_Cs,
							const vector<TGCColor> &b_Cs,
							const TGCColor         &c,
							double &ef,//fore energy
							double &eb)//back energy
{
	double df = 1000000000;
	double db = 1000000000;
	for( int i=0, s=(int)f_Cs.size(); i<s; ++i)  df = min(df, distSq( f_Cs[i], c ) );
	for( int i=0, s=(int)b_Cs.size(); i<s; ++i)  db = min(db, distSq( b_Cs[i], c ) );
	df = sqrt( df );
	db = sqrt( db );

	if( df == 0 && db == 0 ) ef = eb = 0.5;
	else{
		ef =  db / ( df + db );
		eb =  df / ( df + db );
	}
}




/*-----------------------------------------------------------------------------
//water shad segmentationからTCGgraphを作成する
//int W, int H, int D           --> voxel width, height, depth
//const vector<int> &wsd_labels --> 各voxelのwatershad segmentationのregion id
//const byte* volImg            --> volume image (intensity image)
// input label値は、正整数, 最小=1, 連続 
-----------------------------------------------------------------------------*/
void TGCgraph3D::createByWsdImage(int W, int H, int D, const byte* imgR, const byte* imgG, const byte* imgB, const vector<int> &wsd_labels )
{
	clear();
	//create nodes//////////////////////////////////////////////////
	m_pixNodeID.resize( wsd_labels.size() );
	
	int maxLabel = 0;
	for( int i= 0, s = (int) m_pixNodeID.size(); i < s; ++i ){
		m_pixNodeID[i] = wsd_labels[i] - 1;
		maxLabel       = max( maxLabel, wsd_labels[i] );
	}

	m_n_nodes.resize( maxLabel );
	m_n_edges.resize( maxLabel );
	
	for( int z = 0, idx=0; z < D ; ++z)
	for( int y = 0       ; y < H ; ++y)
	for( int x = 0       ; x < W ; ++x, ++idx)
	{
		m_n_nodes[ m_pixNodeID[idx] ].addPixel( idx, imgR==0 ? 0 : imgR[ idx ],    
			                                         imgG==0 ? 0 : imgG[ idx ],
													 imgB==0 ? 0 : imgB[ idx ]);
		//edgeを挿入. 右と下のpixel   E2はまだ計算しない　     この辺を確認中 
		if( x < W - 1 ) addEdge( m_pixNodeID[idx], m_pixNodeID[idx + 1  ] );
		if( y < H - 1 ) addEdge( m_pixNodeID[idx], m_pixNodeID[idx + W  ] );
		if( z < D - 1 ) addEdge( m_pixNodeID[idx], m_pixNodeID[idx + W*H] );
	}
	fprintf( stderr, "On creating TGraphCut Graph node num = %d \n", maxLabel );
	fprintf( stderr, "GCnode memory size  node = %d  %d  edge %d %d\n", sizeof(TGCnode         ), sizeof(TGCnode         ) * maxLabel, 
		                                                                sizeof(pair<int,double>), sizeof(pair<int,double>) * m_n_edges.size()  );
}



/*-------------------------------------------------------
volInOut: voxelをsegmentationの対象にするかどうか
 0      : out of the segmentaton target
 1-255  : segmentation target
-------------------------------------------------------*/
void TGCgraph3D::updateEnableWsdNode( byte* volInOut )
{
	if( m_n_nodes.size() == 0 ){fprintf( stderr ,"wsd node should be generated before calling this !\n"); return;}

	for( int i=0, s = (int) m_n_nodes.size(); i < s; ++i)
	{
		TGCnode &n   = m_n_nodes[i];
		n.m_isEnable = false;
		for( int k = 0; k < (int) n.m_pixelIDs.size(); ++k) if( volInOut[ n.m_pixelIDs[k]] != 0 )
		{
			n.m_isEnable = true;
			break;
		}
	}
}



/*------------------------------------------------------------
check the conflict in the foreground constraints and background constraints.
If one(same) node has fore\back consts, we divide the node 
-------------------------------------------------------------*/
void TGCgraph3D::preCheck_ConflictWsdNodes( const int W, const int H, const int D, 
											const byte *imgR,
											const byte *imgG,
											const byte *imgB,
											const set<int> &forePixIds,
											const set<int> &backPixIds)
{
	bool modified = false;
	for( set<int>::const_iterator fit = forePixIds.begin(); fit != forePixIds.end(); ++fit) 
	for( set<int>::const_iterator bit = backPixIds.begin(); bit != backPixIds.end(); ++bit)
	{
		if( m_pixNodeID[*fit] == m_pixNodeID[*bit] && *fit != *bit )
		{
			divideConflictWsdNodes( W, H, D, imgR, imgG, imgB, *fit, *bit);
			modified = true;
		}
	}

	if( !modified ) return;
	//reconstruct edges
	m_n_edges.clear();
	m_n_edges.resize( m_n_nodes.size() );
	for( int z = 0, idx=0; z < D ; ++z)
	for( int y = 0       ; y < H ; ++y)
	for( int x = 0       ; x < W ; ++x, ++idx)
	{
		if( x < W - 1 ) addEdge( m_pixNodeID[idx], m_pixNodeID[idx + 1  ] );
		if( y < H - 1 ) addEdge( m_pixNodeID[idx], m_pixNodeID[idx + W  ] );
		if( z < D - 1 ) addEdge( m_pixNodeID[idx], m_pixNodeID[idx + W*H] );
	}
}





/*---------------------------------------------------------------
一つのnodeにfore制約(foreVoxIdx)とback制約(backVoxIdx)が存在するとき
このWsdNodeを分割する

変更する物
vector<map<int, double>>  m_WsdEdges ;//edge i jがある時には, m_edges[i][a].first = j, m_edges[i][?].second = error = E2 
vector<TGCnode         >  m_WsdNodes ;//node (region) idx --> pixel         idx  ( m_pixelIDs )
vector<int             >  m_WsdLabels;//pixel         idx --> node (region) idx 

Node num will be incremented.
This method do not modify edges. 
----------------------------------------------------------------*/
void TGCgraph3D::divideConflictWsdNodes(  const int W, const int H, const int D, 
										const byte *imgR, 
										const byte *imgG, 
										const byte *imgB,
										int foreVoxIdx, 
										int backVoxIdx )
{
	if( m_pixNodeID[ foreVoxIdx ] != m_pixNodeID[ backVoxIdx ] ) return;

	const int nId = m_pixNodeID[ foreVoxIdx ];
	const int WH = W*H;

	list<int> fBoundVoxIds; fBoundVoxIds.push_back( foreVoxIdx ); 
	list<int> bBoundVoxIds; bBoundVoxIds.push_back( backVoxIdx ); 
	set<int> fVox; fVox.insert(foreVoxIdx );
	set<int> bVox; bVox.insert(backVoxIdx );

	int minX= 100000, minY= 100000, minZ= 100000;
	int maxX=-100000, maxY=-100000, maxZ=-100000;
	//growth
	while( !bBoundVoxIds.empty() || !fBoundVoxIds.empty())
	{
		int x,y,z;

		if( !bBoundVoxIds.empty() )//fore
		{
			int i = bBoundVoxIds.front(); bBoundVoxIds.pop_front( );
			calcVoxelXYZfromIdx( i, W,H,D, x, y, z);

			if( x>0  && m_pixNodeID[i-1 ]==nId && fVox.find(i-1 )==fVox.end()&&bVox.find(i-1 )==bVox.end()){ bVox.insert(i-1 ); bBoundVoxIds.push_back(i-1 ); }
			if( x<W-1&& m_pixNodeID[i+1 ]==nId && fVox.find(i+1 )==fVox.end()&&bVox.find(i+1 )==bVox.end()){ bVox.insert(i+1 ); bBoundVoxIds.push_back(i+1 ); }
			if( y>0  && m_pixNodeID[i-W ]==nId && fVox.find(i-W )==fVox.end()&&bVox.find(i-W )==bVox.end()){ bVox.insert(i-W ); bBoundVoxIds.push_back(i-W ); }
			if( y<H-1&& m_pixNodeID[i+W ]==nId && fVox.find(i+W )==fVox.end()&&bVox.find(i+W )==bVox.end()){ bVox.insert(i+W ); bBoundVoxIds.push_back(i+W ); }
			if( z>0  && m_pixNodeID[i-WH]==nId && fVox.find(i-WH)==fVox.end()&&bVox.find(i-WH)==bVox.end()){ bVox.insert(i-WH); bBoundVoxIds.push_back(i-WH); }
			if( z<D-1&& m_pixNodeID[i+WH]==nId && fVox.find(i+WH)==fVox.end()&&bVox.find(i+WH)==bVox.end()){ bVox.insert(i+WH); bBoundVoxIds.push_back(i+WH); }
			minX=min(minX,x); minY=min(minY,y); minZ=min(minZ,z);
			maxX=max(maxX,x); maxY=max(maxY,y); maxZ=max(maxZ,z);
		}

		if( !fBoundVoxIds.empty() )
		{
			int i = fBoundVoxIds.front(); fBoundVoxIds.pop_front();
			calcVoxelXYZfromIdx( i, W,H,D, x, y, z);

			if( x>0  && m_pixNodeID[i-1 ]==nId && fVox.find(i-1 )==fVox.end()&&bVox.find(i-1 )==bVox.end()){ fVox.insert(i-1 ); fBoundVoxIds.push_back(i-1 ); }
			if( x<W-1&& m_pixNodeID[i+1 ]==nId && fVox.find(i+1 )==fVox.end()&&bVox.find(i+1 )==bVox.end()){ fVox.insert(i+1 ); fBoundVoxIds.push_back(i+1 ); }
			if( y>0  && m_pixNodeID[i-W ]==nId && fVox.find(i-W )==fVox.end()&&bVox.find(i-W )==bVox.end()){ fVox.insert(i-W ); fBoundVoxIds.push_back(i-W ); }
			if( y<H-1&& m_pixNodeID[i+W ]==nId && fVox.find(i+W )==fVox.end()&&bVox.find(i+W )==bVox.end()){ fVox.insert(i+W ); fBoundVoxIds.push_back(i+W ); }
			if( z>0  && m_pixNodeID[i-WH]==nId && fVox.find(i-WH)==fVox.end()&&bVox.find(i-WH)==bVox.end()){ fVox.insert(i-WH); fBoundVoxIds.push_back(i-WH); }
			if( z<D-1&& m_pixNodeID[i+WH]==nId && fVox.find(i+WH)==fVox.end()&&bVox.find(i+WH)==bVox.end()){ fVox.insert(i+WH); fBoundVoxIds.push_back(i+WH); }
			minX=min(minX,x); minY=min(minY,y); minZ=min(minZ,z);
			maxX=max(maxX,x); maxY=max(maxY,y); maxZ=max(maxZ,z);
		}
	}

	//update m_pixNodeID/m_n_edges/m_n_nodes
	m_n_nodes.push_back( TGCnode()  ); //pushback single node 
	m_n_edges.push_back( set<int>() );

	const int fNodeId = nId;
	const int bNodeId = (int)m_n_nodes.size()-1;
	TGCnode &foreN = m_n_nodes[fNodeId];
	TGCnode &backN = m_n_nodes[bNodeId];

	foreN = TGCnode();
	for( set<int>::iterator it = fVox.begin(); it != fVox.end(); ++it){
		m_pixNodeID[ *it ] = fNodeId;
		m_n_nodes[ fNodeId ].addPixel( (*it), imgR==0 ? 0 : imgR[ *it ],    
			                                  imgG==0 ? 0 : imgG[ *it ],
											  imgB==0 ? 0 : imgB[ *it ]);
	}
	for( set<int>::iterator it = bVox.begin(); it != bVox.end(); ++it){
		m_pixNodeID[ *it ] = bNodeId;
		m_n_nodes[ bNodeId ].addPixel( (*it), imgR==0 ? 0 : imgR[ *it ],    
			                                  imgG==0 ? 0 : imgG[ *it ],
											  imgB==0 ? 0 : imgB[ *it ]);
	}
}




//imgG, imgB can be "0" if the image has smaller channels
void TGCgraph3D::calcMinCut(int W, int H, int D, const byte *imgR, 
											   const byte *imgG, 
											   const byte *imgB, 
											   const set<int> &forePixIds,
											   const set<int> &backPixIds, 
											   vector<int> &isForeBack)
{
	preCheck_ConflictWsdNodes( W, H, D, imgR, imgG, imgB, forePixIds, backPixIds );

	//fore/back pixels--> fore/back nodes
	set< int > foreNs, backNs;
	set<int>::const_iterator it;
	for( it = forePixIds.begin(); it != forePixIds.end(); ++it ) foreNs.insert( m_pixNodeID[ *it ] );
	for( it = backPixIds.begin(); it != backPixIds.end(); ++it ) backNs.insert( m_pixNodeID[ *it ] );

	//fore/back node --> core/back color
	vector<TGCColor> foreCs, backCs; 
	foreCs.reserve(foreNs.size());
	backCs.reserve(backNs.size());
	for( it = foreNs.begin(); it != foreNs.end(); ++it ) foreCs.push_back( m_n_nodes[ *it].m_c);
	for( it = backNs.begin(); it != backNs.end(); ++it ) backCs.push_back( m_n_nodes[ *it].m_c);

	//construct graph
	const int nSize = (int)m_n_nodes.size();
	Graph *g = new Graph();
	Graph::node_id *nodes = new Graph::node_id[ nSize ];

	// insert node and compute E1/////////////////////////////////////////////////////////////////////////
	for(int i = 0; i < nSize; ++i )
	{
		nodes[i] = g->add_node();
		
		if( m_n_nodes[i].m_isEnable ){
			double e1_f, e1_b;
			calcE1( foreCs, backCs, m_n_nodes[i].m_c, e1_f, e1_b);
			g->set_tweights( nodes[i], e1_f, e1_b);
		}
		else
			g->set_tweights( nodes[i], 0, FOREBACK_MAX);
	}

	// insert inf into fore/back constrained nodes
	for( it=foreNs.begin(); it != foreNs.end(); ++it) if(m_n_nodes[*it].m_isEnable){ g->set_tweights( nodes[*it], FOREBACK_MAX, 0); }
	for( it=backNs.begin(); it != backNs.end(); ++it) if(m_n_nodes[*it].m_isEnable){ g->set_tweights( nodes[*it], 0, FOREBACK_MAX); }

	// insert edge and compute E2/////////////////////////////////////////////////////////////////////////
	for(int i = 0; i < nSize; ++i) if(m_n_nodes[i].m_isEnable)
	{
		for( set<int>::const_iterator eIt = m_n_edges[i].begin(); eIt != m_n_edges[i].end(); ++eIt)
		if( i < *eIt && m_n_nodes[ *eIt ].m_isEnable )
		{
			double e2 = calcE2( m_n_nodes[ i ].m_c, m_n_nodes[ *eIt ].m_c, m_lambda);
			g->add_edge( nodes[ i ], nodes[ *eIt ], e2,e2); 
		}
	}

	//graph cut//
	Graph::flowtype flow = g->maxflow();  fprintf( stderr, "flow = %f\n", flow );
	
	//update fore/back nodes --> fore/back pixels
	for(int i = 0; i < nSize; ++i) m_n_nodes[i].m_isFore = ( g->what_segment( nodes[i] ) == Graph::SOURCE);

	isForeBack.resize( m_pixNodeID.size() );
	for( int i = 0, s = (int)m_pixNodeID.size(); i < s ; ++i) 
		isForeBack[ i ] = m_n_nodes[ m_pixNodeID[i] ].m_isFore ? 1 : 0;

	delete g;
	delete [] nodes;
}



struct TGCpixelNode
{
	int      m_pixID;
	TGCColor m_c;
	inline TGCpixelNode(int pixelIdx, const byte r, const byte g, const byte b)
	{
		m_pixID = pixelIdx;
		m_c.r=r;
		m_c.g=g;
		m_c.b=b;
	}
};

/*-------------------------------------------------------------------------------
arguments
 bVoxInOut 0 : out of segmentation target
           1 : segmentation target & initial back ground
       other : segmentation target & initial FORE ground

result
 bVoxInOut 0 : out of segmentation target
           1 : back ground
       other : FORE ground

境界からbound width分の帯状の部分についてのみ計算する
色の情報を考慮せず、E2のみを利用
--------------------------------------------------------------------------------*/
void TGCgraph3D::calcMinCut_voxelLv_bound(int W, int H, int D, 
										const byte *imgR, 
										const byte *imgG, 
										const byte *imgB, byte *bVoxInOut, int boundWidth)
{
	const int vSize = W * H * D;
	const int WH    = W * H;

	/////////////extract target pixel////////////////////////////////////////////////////////////////////
	int   *voxId2nodeId = new int[ vSize ]; // target node: pixNodeIdx  /  fore=-1  /   back=-2
	for( int i = 0; i < vSize; ++i) voxId2nodeId[i] =  (bVoxInOut[i] ==0 || bVoxInOut[i] == 1) ? -2 : -1;

	short *voxelStates  = new short[ vSize ]        ; //0: out of Bound,  1: in bound belt;
	memset( voxelStates, 0, sizeof( short )* vSize );


	for( int z = 1; z < D ; ++z)
	for( int y = 1; y < H ; ++y)
	for( int x = 1; x < W ; ++x)
	{
		int idx = x + y * W + z * WH;

		if( voxId2nodeId[ idx ] != voxId2nodeId[ idx-1  ] || 
			voxId2nodeId[ idx ] != voxId2nodeId[ idx-W  ] ||
			voxId2nodeId[ idx ] != voxId2nodeId[ idx-WH ] )
		{
			//flagを立てる
			for( int z_off = -boundWidth; z_off <= boundWidth; ++z_off) if( 0 <= z+z_off && z+z_off < D )
			for( int y_off = -boundWidth; y_off <= boundWidth; ++y_off) if( 0 <= y+y_off && y+y_off < H )
			for( int x_off = -boundWidth; x_off <= boundWidth; ++x_off) if( 0 <= x+x_off && x+x_off < W )
				voxelStates[idx + x_off 
					            + y_off * W 
				                + z_off * WH ] = 1;
		}
	}
	
	////////////create target pixel nodes
	vector< TGCpixelNode* > voxNodes;  
	voxNodes.reserve( vSize / 4 );

	for( int i = 0; i < vSize; ++i) if( voxelStates[i] == 1){
		voxId2nodeId[ i ] = (int) voxNodes.size();
		voxNodes.push_back( new TGCpixelNode( i, imgR==0?0:imgR[i], 
			                                     imgG==0?0:imgG[i],
												 imgB==0?0:imgB[i] )  );
	}

	//construct graph
	Graph::node_id *nodes = new Graph::node_id[ voxNodes.size() ];
	Graph *g = new Graph();
	for(int i = 0, s = (int)voxNodes.size(); i < s; ++i ) nodes[i] = g->add_node();	

	for( int nodeI = 0, s = (int)voxNodes.size(); nodeI < s; ++nodeI)
	{
		int voxIdx = voxNodes[nodeI]->m_pixID;
		int voxX, voxY, voxZ;
		calcVoxelXYZfromIdx( voxIdx, W, H, D, voxX, voxY, voxZ);

		bool edgeToNexX = false;
		bool edgeToNexY = false;
		bool edgeToNexZ = false;

		int voxCondition = 0;// -1 : adjacent to back pix   / 1:adjucent to forepixel
		if( voxX != 0   ){ if     ( voxId2nodeId[ voxIdx-1 ] == -1 ) voxCondition =  1;
			               else if( voxId2nodeId[ voxIdx-1 ] == -2 ) voxCondition = -1;}
		if( voxY != 0   ){ if     ( voxId2nodeId[ voxIdx-W ] == -1 ) voxCondition =  1;
			               else if( voxId2nodeId[ voxIdx-W ] == -2 ) voxCondition = -1;}
		if( voxZ != 0   ){ if     ( voxId2nodeId[ voxIdx-WH] == -1 ) voxCondition =  1;
			               else if( voxId2nodeId[ voxIdx-WH] == -2 ) voxCondition = -1;}

		if( voxX != W-1 ){ if     ( voxId2nodeId[ voxIdx+1 ] == -1 ) voxCondition =  1;
						   else if( voxId2nodeId[ voxIdx+1 ] == -2 ) voxCondition = -1;
						   else edgeToNexX = true;                                     }
		if( voxZ != D-1 ){ if     ( voxId2nodeId[ voxIdx+WH] == -1 ) voxCondition =  1;
						   else if( voxId2nodeId[ voxIdx+WH] == -2 ) voxCondition = -1;
						   else edgeToNexZ = true;                                     }
		if( voxY != H-1 ){ if     ( voxId2nodeId[ voxIdx+W ] == -1 ) voxCondition =  1;
						   else if( voxId2nodeId[ voxIdx+W ] == -2 ) voxCondition = -1;
						   else edgeToNexY = true;                                     }

		double E1_f,E1_b;
		if(     voxCondition == 1 ){ E1_f = FOREBACK_MAX; E1_b = 0; }
		else if(voxCondition ==-1 ){ E1_b = FOREBACK_MAX; E1_f = 0; }
		else {                       E1_f = E1_b = 1;               }
		g->set_tweights( nodes[nodeI], E1_f, E1_b);
		
		//add edges
		if( edgeToNexX){
			int nodeJ = voxId2nodeId[ voxIdx+1 ];
			double E2 = calcE2( voxNodes[nodeI]->m_c, voxNodes[nodeJ]->m_c, m_lambda); 
			g->add_edge( nodes[ nodeI ], nodes[ nodeJ ], E2, E2);
		}
		if( edgeToNexY ){
			int nodeJ = voxId2nodeId[ voxIdx+W ];
			double E2 = calcE2( voxNodes[nodeI]->m_c, voxNodes[nodeJ]->m_c, m_lambda); 
			g->add_edge( nodes[ nodeI ], nodes[ nodeJ ], E2, E2);
		}
		if( edgeToNexZ ){
			int nodeJ = voxId2nodeId[ voxIdx+WH ];
			double E2 = calcE2( voxNodes[nodeI]->m_c, voxNodes[nodeJ]->m_c, m_lambda); 
			g->add_edge( nodes[ nodeI ], nodes[ nodeJ ], E2, E2);
		}	
	}
	
	
	//graph cut//
	Graph::flowtype flow = g->maxflow();  fprintf( stderr, "flow = %f\n", flow );

	for( int i = 0; i < vSize; ++i) if( bVoxInOut[i] != 0 )
	{
		bVoxInOut[i] = ( voxId2nodeId[i]==-1) ?  255 :
			           ( voxId2nodeId[i]==-2) ?  1 : 
					   ( g->what_segment( nodes[ voxId2nodeId[i] ] ) == Graph::SOURCE) ? 255 : 1;
	}

	delete g;
	delete[] nodes;
	delete[] voxelStates;
	delete[] voxId2nodeId;
	for( int i= 0; i < (int) voxNodes.size(); ++i) delete voxNodes[i]; voxNodes.clear();
}





//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//2D//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void TGCgraph2D::createByWsdImage (int W, int H, const byte* rgba, const vector<int> &wsd_labels )
{
	clear();
	//create nodes//////////////////////////////////////////////////
	m_pixNodeID.resize( wsd_labels.size() );
	
	int maxLabel = 0;
	for( int i= 0, s = (int) m_pixNodeID.size(); i < s; ++i ){
		m_pixNodeID[i] = wsd_labels[i] - 1;
		maxLabel       = max( maxLabel, wsd_labels[i] );
	}

	m_n_nodes.resize( maxLabel );
	m_n_edges.resize( maxLabel );
	
	for( int y = 0, idx=0; y < H ; ++y)
	for( int x = 0       ; x < W ; ++x, ++idx)
	{
		m_n_nodes[ m_pixNodeID[idx] ].addPixel( idx, rgba[idx*4+0],
			                                         rgba[idx*4+1],
													 rgba[idx*4+2]); 
		//edgeを挿入. 右と下のpixel   E2はまだ計算しない　     この辺を確認中 
		if( x < W - 1 ) addEdge( m_pixNodeID[idx], m_pixNodeID[idx + 1  ] );
		if( y < H - 1 ) addEdge( m_pixNodeID[idx], m_pixNodeID[idx + W  ] );
	}
	fprintf( stderr, "On creating TGraphCut Graph node num = %d \n", maxLabel );
	fprintf( stderr, "GCnode memory size  node = %d  %d  edge %d %d\n", sizeof(TGCnode         ), sizeof(TGCnode         ) * maxLabel, 
		                                                                sizeof(pair<int,double>), sizeof(pair<int,double>) * m_n_edges.size()  );
}


void TGCgraph2D::calcMinCut(const set<int> &forePixIds, 
							const set<int> &backPixIds, 
							   vector<int> &isForeBack)
{
	//fore/back pixels--> fore/back nodes
	set< int > foreNs, backNs;
	set<int>::const_iterator it;
	for( it = forePixIds.begin(); it != forePixIds.end(); ++it ) foreNs.insert( m_pixNodeID[ *it ] );
	for( it = backPixIds.begin(); it != backPixIds.end(); ++it ) backNs.insert( m_pixNodeID[ *it ] );

	//fore/back node --> core/back color
	vector<TGCColor> foreCs, backCs; 
	foreCs.reserve(foreNs.size());
	backCs.reserve(backNs.size());
	for( it = foreNs.begin(); it != foreNs.end(); ++it ) foreCs.push_back( m_n_nodes[ *it].m_c);
	for( it = backNs.begin(); it != backNs.end(); ++it ) backCs.push_back( m_n_nodes[ *it].m_c);

	//construct graph
	const int nSize = (int)m_n_nodes.size();
	Graph *g = new Graph();
	Graph::node_id *nodes = new Graph::node_id[ nSize ];

	// insert node and compute E1/////////////////////////////////////////////////////////////////////////
	for(int i = 0; i < nSize; ++i )
	{
		nodes[i] = g->add_node();
		
		if( m_n_nodes[i].m_isEnable ){
			double e1_f, e1_b;
			calcE1( foreCs, backCs, m_n_nodes[i].m_c, e1_f, e1_b);
			g->set_tweights( nodes[i], e1_f, e1_b);
		}
		else
			g->set_tweights( nodes[i], 0, FOREBACK_MAX);
	}

	// insert inf into fore/back constrained nodes
	for( it=foreNs.begin(); it != foreNs.end(); ++it) if(m_n_nodes[*it].m_isEnable){ g->set_tweights( nodes[*it], FOREBACK_MAX, 0); }
	for( it=backNs.begin(); it != backNs.end(); ++it) if(m_n_nodes[*it].m_isEnable){ g->set_tweights( nodes[*it], 0, FOREBACK_MAX); }

	// insert edge and compute E2/////////////////////////////////////////////////////////////////////////
	for(int i = 0; i < nSize; ++i) if(m_n_nodes[i].m_isEnable)
	{
		for( set<int>::const_iterator eIt = m_n_edges[i].begin(); eIt != m_n_edges[i].end(); ++eIt)
		if( i < *eIt && m_n_nodes[ *eIt ].m_isEnable )
		{
			double e2 = calcE2( m_n_nodes[ i ].m_c, m_n_nodes[ *eIt ].m_c, m_lambda);
			g->add_edge( nodes[ i ], nodes[ *eIt ], e2,e2); 
		}
	}

	//graph cut//
	Graph::flowtype flow = g->maxflow();  fprintf( stderr, "flow = %f\n", flow );
	
	//update fore/back nodes --> fore/back pixels
	for(int i = 0; i < nSize; ++i) m_n_nodes[i].m_isFore = ( g->what_segment( nodes[i] ) == Graph::SOURCE);

	isForeBack.resize( m_pixNodeID.size() );
	for( int i = 0, s = (int)m_pixNodeID.size(); i < s ; ++i) 
		isForeBack[ i ] = m_n_nodes[ m_pixNodeID[i] ].m_isFore ? 1 : 0;

	delete g;
	delete [] nodes;

}



/*
境界からbound width分の帯状の部分についてのみ計算する
*/
void TGCgraph2D::calcMinCut_pixLv_bound(int W, int H, byte* rgba,
										const set<int> &forePixIdxs, 
										const set<int> &backPixIdxs, vector<int> &isForeBack, int boundWidth)
{
	const int WH = W * H;
	
	//対象pixelを抽出//////////////////////////////////////

	//node以外はfore=-1 back=-2, nodeにはnodeIdxが入る
	int   *pixId2nodeId = new int  [ WH ]; 
	for( int i = 0; i < WH; ++i) pixId2nodeId[i] = m_n_nodes[ m_pixNodeID[i] ].m_isFore ? -1 : -2;

	//0: is not in the bound,  1: is in the bound;
	short *pixelStates  = new short[ WH ]; 
	memset( pixelStates, 0, sizeof( short )* WH );

	for( int y = 1; y < H; ++y)
	for( int x = 1; x < W ; ++x) //手前との差をみるから1からスタート
	{
		int idx = x + y * W;
		if( pixId2nodeId[ idx ] != pixId2nodeId[ idx-1 ] || 
			pixId2nodeId[ idx ] != pixId2nodeId[ idx-W ])
		{
			//flagを立てる
			for( int y_off = -boundWidth; y_off <= boundWidth; ++y_off)
			for( int x_off = -boundWidth; x_off <= boundWidth; ++x_off)
			if( 0 <= x+x_off && x+x_off < W && 0 <= y+y_off && y+y_off < H) pixelStates[idx + x_off + y_off * W] = 1;
		}
	}
	////////////計算すべきpixelのindexを計算//////////////////////////////////////
	vector< TGCpixelNode* > pixNodes;  
	pixNodes.reserve( WH / 4 );

	for( int i = 0; i < WH; ++i) if( pixelStates[i] == 1)
	{
		pixId2nodeId[ i ] = (int) pixNodes.size();
		pixNodes.push_back( new TGCpixelNode( i, rgba[4*i], rgba[4*i+1], rgba[4*i+2] ));
	}
	////////////////////////graph構築//////////////////////////////////////
	Graph::node_id *nodes = new Graph::node_id[ pixNodes.size() ];
	Graph *g = new Graph();
	for(int i = 0; i < (int)pixNodes.size(); ++i ) nodes[i] = g->add_node();	

	//pixel --> color 抽出
	vector< TGCColor > foreColors(forePixIdxs.size()),
		               backColors(backPixIdxs.size());
	{
		int i; set<int>::const_iterator it;
		for(i=0, it = forePixIdxs.begin(); it != forePixIdxs.end(); ++it,++i ) foreColors[i].Set( rgba[4*(*it)], rgba[4*(*it)+1], rgba[4*(*it)+2] );
		for(i=0, it = backPixIdxs.begin(); it != backPixIdxs.end(); ++it,++i ) backColors[i].Set( rgba[4*(*it)], rgba[4*(*it)+1], rgba[4*(*it)+2] );
	}

	for( int nodeI = 0; nodeI < (int) pixNodes.size(); ++nodeI)
	{
		int pixIdx = pixNodes[nodeI]->m_pixID;
		int pixX   = pixIdx  %  W;
		int pixY   = pixIdx  /  W;
		bool edgeToRight = false;
		bool edgeToBelow = false;
		
		int pixCondition = 0;// -1 backpixelが隣接  1 forepixelが隣接
		//上下左右をみてfore/backに隣接していたら、E1f/E1b = 無限
		if( pixX != 0   ) { if     ( pixId2nodeId[ pixIdx-1 ] == -1 ) pixCondition =  1;
			                else if( pixId2nodeId[ pixIdx-1 ] == -2 ) pixCondition = -1;}
		if( pixY != 0   ){  if     ( pixId2nodeId[ pixIdx-W ] == -1 ) pixCondition =  1;
			                else if( pixId2nodeId[ pixIdx-W ] == -2 ) pixCondition = -1;}

		if( pixX != W-1 ) { if     ( pixId2nodeId[ pixIdx+1 ] == -1 ) pixCondition =  1;
							else if( pixId2nodeId[ pixIdx+1 ] == -2 ) pixCondition = -1; 
							else edgeToRight = true; }
		if( pixY != H-1 ){  if     ( pixId2nodeId[ pixIdx+W ] == -1 ) pixCondition =  1;
							else if( pixId2nodeId[ pixIdx+W ] == -2 ) pixCondition = -1; 
							else edgeToBelow = true; }

		double E1_f,E1_b;
		if(     pixCondition == 1 ){ E1_f = FOREBACK_MAX; E1_b = 0; }
		else if(pixCondition ==-1 ){ E1_b = FOREBACK_MAX; E1_f = 0; }
		else {                      E1_f = E1_b = 0.5; }//calcE1( foreColors, backColors, pixNodes[nodeI]->m_c, E1_f, E1_b); }
		g->set_tweights( nodes[nodeI], E1_f, E1_b);
		
		//edgeを追加
		if( edgeToRight ){
			int nodeJ = pixId2nodeId[ pixIdx+1 ];
			double E2 = calcE2( pixNodes[nodeI]->m_c, pixNodes[nodeJ]->m_c, m_lambda); 
			g->add_edge( nodes[ nodeI ], nodes[ nodeJ ], E2, E2);
		}
		if( edgeToBelow ){
			int nodeJ = pixId2nodeId[ pixIdx+W ];
			double E2 = calcE2( pixNodes[nodeI]->m_c, pixNodes[nodeJ]->m_c, m_lambda ); 
			g->add_edge( nodes[ nodeI ], nodes[ nodeJ ], E2, E2);
		}	
	}

	//graph cut//
	fprintf( stderr, "node size == %d  ", (int)pixNodes.size() );
	Graph::flowtype flow = g->maxflow();  fprintf( stderr, "flow = %f\n", flow );

	isForeBack.resize( m_pixNodeID.size() );
	for( int i = 0; i < WH; ++i)
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




