#include "StdAfx.h"
#include "TWatershed.h"
#include "math.h"

using namespace std ;



inline static double ti_pow(double d,int n){
	if(n < 0){
		d = 1/d;
		n = -n;
	}
	
	double result = 1.0;
	for(int i = 0;i < n;++i) result *= d;
	return result;
}









/*
gradient magnitudeを利用してwatershedを行う
8近傍のうちで一番差の大きい物をgradient magnitudeとして利用.

rgbaは, byteの配列で、RGBA値が順に入る。
i番目の色は
R = rgba[4 * i+ 0]
G = rgba[4 * i+ 1]
B = rgba[4 * i+ 2]
A = rgba[4 * i+ 3]
*/
void TWaterShed_gradOfRgbaImage2D(int sizeW, int sizeH, const byte* rgba, vector<int> &labels)
{
	byte *img = new byte[ sizeW*sizeH ];
	
	for( int y = 0; y < sizeH; ++y)
	for( int x = 0; x < sizeW; ++x)
	{
		double a0 = 0, a1 = 0, a2 = 0, a3 = 0, 
			   a4 = 0, a5 = 0, a6 = 0, a7 = 0;
		
		int idx = 4 * ( y * sizeW + x );
		if( x > 0                     ) a0 = n_colorDist2( &rgba[idx], &rgba[idx - 4]            );
		if( x < sizeW-1               ) a1 = n_colorDist2( &rgba[idx], &rgba[idx + 4]            );
		if( y > 0                     ) a2 = n_colorDist2( &rgba[idx], &rgba[idx     - 4 * sizeW]);
		if( y < sizeH-1               ) a3 = n_colorDist2( &rgba[idx], &rgba[idx     + 4 * sizeW]);
		if( x > 0       && y > 0      ) a4 = n_colorDist2( &rgba[idx], &rgba[idx - 4 - 4 * sizeW]);
		if( x < sizeW-1 && y > 0      ) a5 = n_colorDist2( &rgba[idx], &rgba[idx + 4 - 4 * sizeW]);
		if( x > 0       && y < sizeH-1) a6 = n_colorDist2( &rgba[idx], &rgba[idx - 4 + 4 * sizeW]);
		if( x < sizeW-1 && y < sizeH-1) a7 = n_colorDist2( &rgba[idx], &rgba[idx + 4 + 4 * sizeW]);
		
		double a = max( max( max(a0, a1), max(a2, a3) ), max( max( a4, a5), max( a6, a7) ) );

		double dirMagnitude = sqrt((double)a ) / sqrt( 3.0 ) * 2;
		dirMagnitude *= 255;
		if( dirMagnitude > 255 ) dirMagnitude = 255;
		img[y * sizeW + x] = (byte) dirMagnitude;
	}

	TWatershed2D( sizeW, sizeH, img, labels);
	delete[] img;
}



/*　一ピクセルしか含まないregionを、watershed 領域として埋める
	labelを振り直す, 最小=1, 連続する正整数
*/
void postRemoveSinglePixelRegion( std::vector<int> &pixLabels)
{
	int maxLabel = 0;
	for( vector<int>::iterator it = pixLabels.begin(); it != pixLabels.end(); ++it ) maxLabel = max( maxLabel, *it);
	
	//各領域に含まれる要素数を計算
	vector<int> regionNum( maxLabel + 1);//0含む
	for( vector<int>::iterator it = pixLabels.begin(); it != pixLabels.end(); ++it ) ++regionNum[*it];
	
	//1 pixelしか入っていない領域をつぶす
	for( vector<int>::iterator it = pixLabels.begin(); it != pixLabels.end(); ++it )
	{
		if( regionNum[*it] == 1 )
		{
			regionNum[*it] = -1;
			*it = 0;//つぶす
		}
	}

	//labelを連続に振り直す
	int offset = 0;
	vector<int> newLabel( regionNum.size() );
	for( int i = 0; i < (int) regionNum.size(); ++i)
	{
		if( regionNum[i] == -1 ) offset++;
		else newLabel[i] = i-offset;
	}

	for( int i = 0; i < (int) pixLabels.size(); ++i)
		pixLabels[i] = newLabel[ pixLabels[i] ];


#if 0
	fprintf( stderr, "debug of posRemove...\n");
	//debug
	maxLabel = 0;
	for( vector<int>::iterator it = pixLabels.begin(); it != pixLabels.end(); ++it ) maxLabel = max( maxLabel, *it);
	for( int i = 0; i < (int) maxLabel; ++i)
	{
		bool canFind = false;
		for( vector<int>::iterator it = pixLabels.begin(); it != pixLabels.end(); ++it ) if( *it == i ) { canFind = true; break;}
			

		if( !canFind )
			fprintf( stderr, "this method contains bug, here !!!!\n");
	}
	
#endif
}







//watershed pixelを色の近い近傍領域に入れる
void collapseWsdPixels2D(int sizeW, int sizeH, const byte *rgba, std::vector<int> &labels)
{
	vector<int>  visFlg(sizeW * sizeH);

#if 0
	while( true )
	{
		int c = 0;
		for( vector<int>::iterator it = visFlg.begin(); it != visFlg.end(); ++it) *it = 0;

		int idx = 0;
		for( int y = 0; y < sizeH; ++y        )	
		for( int x = 0; x < sizeW; ++x, ++idx ) if( labels[ idx ] <= 0 )
		{	
			int a0 = -1, a1 = -1, a2 = -1, a3 = -1, lab = -1;
			if     ( x > 0       && labels[ idx-1    ] > 0 && visFlg[ idx-1     ] == 0) a0 = labels[ idx-1    ];//左のpixel
			else if( x < sizeW-1 && labels[ idx+1    ] > 0 && visFlg[ idx+1     ] == 0) a1 = labels[ idx+1    ];//右のpixel
			else if( y > 0       && labels[ idx-sizeW] > 0 && visFlg[ idx-sizeW ] == 0) a2 = labels[ idx-sizeW];//上のpixel
			else if( y < sizeH-1 && labels[ idx+sizeW] > 0 && visFlg[ idx+sizeW ] == 0) a3 = labels[ idx+sizeW];//下のpixel
			
			if( a0 == -1 && a1 == -1 && a2 == -1 && a3 == -1 ) continue;
			
			//同じラベルが2個あったらそれを入れる
			if( a0 != -1 && ( a0 == a1 || a0 == a2 || a0 == a3 ) ) lab = a0;
			if( a1 != -1 && ( a1 == a2 || a1 == a3             ) ) lab = a1;
			if( a2 != -1 && ( a2 == a3                         ) ) lab = a2;
			if( lab == -1 ){
				if( a0 != -1 ) lab = a0;
				if( a1 != -1 ) lab = a1;
				if( a2 != -1 ) lab = a2;
				if( a3 != -1 ) lab = a3;
			}

			labels[idx] = lab;
			visFlg[idx] = 1;
			++c;
		}
		if( c == 0 ) break;
	}
#else 
	while( true )
	{
		for( vector<int>::iterator it = visFlg.begin(); it != visFlg.end(); ++it) *it = 0;

		bool changed = false;

		int idx = 0;
		for( int y = 0; y < sizeH; ++y        )	
		for( int x = 0; x < sizeW; ++x, ++idx ) if( labels[ idx ] <= 0 )
		{	
			//左、右、上、下
			double colDiff = 1000000000000;
			for( int c = 0; c<4; ++c)
			{
				int neiIdx = 0;
				if     ( c == 0 && x >   0				      ) neiIdx = idx - 1    ;//左
				else if( c == 1 && x < sizeW-1			      ) neiIdx = idx + 1    ;//右
				else if( c == 2 && y >   0				      ) neiIdx = idx - sizeW;//下
				else if( c == 3 && y < sizeH-1                ) neiIdx = idx + sizeW;//上
				else if( c == 4 && x >   0     && y > 0       ) neiIdx = idx - sizeW - 1;
				else if( c == 5 && x < sizeW-1 && y > 0       ) neiIdx = idx - sizeW + 1;
				else if( c == 6 && x >   0     && y < sizeH-1 ) neiIdx = idx + sizeW - 1;
				else if( c == 7 && x < sizeW-1 && y < sizeH-1 ) neiIdx = idx + sizeW + 1;
				else continue;

				if( labels[ neiIdx ] > 0 && visFlg[ neiIdx ] == 0)
				{
					double d = colorDist2( &rgba[ idx*4  ], &rgba[ neiIdx*4 ] );
					if( d < colDiff )
					{ 
						colDiff = d; 
						labels[idx] = labels[ neiIdx ];
						visFlg[idx] = 1;
						changed = true;
					}
				}
			}
		}

		if( !changed ) break;
	}
#endif

}







////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
//3D.........///////////////////////////////////////////////////////////////////////////////////////



static int powerOfTwo_val(unsigned int val)
{
	unsigned int k = 1;
	for( int i = 0; i < 1024; ++i)
	{
		if( val <= k ) return i;
		k *= 2;
	}
	return 0;
}

void TWaterShed3D_forLargeVol(int W, int H, int D, const byte* gMagImg, std::vector<int> &labels, void (*func_setProgressPos)(const double) )
{
	if( func_setProgressPos != 0 ) func_setProgressPos( 0 );

	//Devide volume in the Z dir and compute Watershad one by one 
	int stepDepth = 0;
	int pW = powerOfTwo_val( W );
	int pH = powerOfTwo_val( H );
	int pD = powerOfTwo_val( D );
	int maxDepthStep_p = 23 - pW - pH;
	if( maxDepthStep_p < 0 ) {fprintf( stderr, "this volume is too big to compute water shed !\n"); return;}

	if( pD <= maxDepthStep_p ){ 
		TWatershed3D( W, H, D, gMagImg, labels );
		return;
	}

	int maxDepthStep = (int)ti_pow(2, maxDepthStep_p);

	int numOfIter = D / maxDepthStep;
	if( D % maxDepthStep > 0 ) numOfIter++;

	vector< vector<int> > segmentLavels( numOfIter, vector<int>() );
	
	for( int i = 0; i < numOfIter; ++i)
	{
		int pivDepth     =   i   * maxDepthStep;
		int pivDepthNext = (i+1) * maxDepthStep;
		if( pivDepthNext-1 <= D-1) TWatershed3D( W, H, maxDepthStep, &gMagImg[ W * H* pivDepth ], segmentLavels[i] );
		else                       TWatershed3D( W, H, D - pivDepth, &gMagImg[ W * H* pivDepth ], segmentLavels[i] );
		
		fprintf( stderr, "%f done\n", (i+1)/ (double)numOfIter * 100);
		if( func_setProgressPos != 0 ) func_setProgressPos( (i+1)/ (double)numOfIter );
	}

	//assemble  multiple label ID lists  into one
	labels.resize( W * H * D, 0);
	int labelOffset = 0;
	int idxOffset   = 0;
	for( int kk = 0; kk < numOfIter; ++kk)
	{
		int maxLabelVal = 0;
		for( int i = 0; i < (int) segmentLavels[kk].size(); ++i) 
		{
			if( segmentLavels[kk][i] == 0 ) labels[ i + idxOffset ] = 0;
			else                            labels[ i + idxOffset ] = segmentLavels[kk][i] + labelOffset;
			maxLabelVal = max( maxLabelVal, segmentLavels[kk][i] );
		}
		labelOffset += maxLabelVal;
		idxOffset   += (int) segmentLavels[kk].size();
	}
}


// insert  "watershed pixels" into  neighboring regions with nearest color
// watershed pixel' label = "0"
//*このインプリだと近傍pixelの色しか見ていない(本来はLabel領域全体を見るべき)
void collapseWsdPixels3D     (int W, int H, int D, const byte* r, const byte* g, const byte* b, std::vector<int> &labels)
{
	vector<byte>  visFlg(W * H * D);
	while( true )
	{
		for( vector<byte>::iterator it = visFlg.begin(); it != visFlg.end(); ++it) *it = 0;

		bool changed = false;

		int idx = 0;
		for( int z = 0; z < D; ++z        )	
		for( int y = 0; y < H; ++y        )	
		for( int x = 0; x < W; ++x, ++idx ) if( labels[ idx ] <= 0 )
		{	
			//左、右、上、下
			double colDiff = 1000000000000;
			for( int c = 0; c<6; ++c)
			{
				int neiIdx = 0;
				if     ( c == 0 && x >   0 ) neiIdx = idx - 1  ;//左
				else if( c == 1 && x < W-1 ) neiIdx = idx + 1  ;//右
				else if( c == 2 && y >   0 ) neiIdx = idx - W  ;//下
				else if( c == 3 && y < H-1 ) neiIdx = idx + W  ;//上
				else if( c == 4 && z >   0 ) neiIdx = idx - W*H;
				else if( c == 5 && z < D-1 ) neiIdx = idx + W*H;
				else continue;

				if( labels[ neiIdx ] <= 0 || visFlg[ neiIdx ] != 0) continue;
				
				double d = (r == 0) ?  0 : (r[idx] - r[neiIdx])*(r[idx] - r[neiIdx])  + 
					       (g == 0) ?  0 : (g[idx] - g[neiIdx])*(g[idx] - g[neiIdx])  + 
						   (b == 0) ?  0 : (b[idx] - b[neiIdx])*(b[idx] - b[neiIdx])  ;

				if( d < colDiff )
				{ 
					colDiff = d; 
					labels[idx] = labels[ neiIdx ];
					visFlg[idx] = 1;
					changed = true;
				}
			}
		}
		if( !changed ) break;
	}
}



void collapseWsdPixels3D(int sizeW, int sizeH, int sizeD, const byte *rgba, std::vector<int> &labels)
{
	vector<int>  visFlg(sizeW * sizeH * sizeD);
	while( true )
	{
		for( vector<int>::iterator it = visFlg.begin(); it != visFlg.end(); ++it) *it = 0;

		bool changed = false;

		int idx = 0;
		for( int z = 0; z < sizeD; ++z        )	
		for( int y = 0; y < sizeH; ++y        )	
		for( int x = 0; x < sizeW; ++x, ++idx ) if( labels[ idx ] <= 0 )
		{	
			//左、右、上、下
			double colDiff = 1000000000000;
			for( int c = 0; c<6; ++c)
			{
				int neiIdx = 0;
				if     ( c == 0 && x >   0	   ) neiIdx = idx - 1    ;//左
				else if( c == 1 && x < sizeW-1 ) neiIdx = idx + 1    ;//右
				else if( c == 2 && y >   0	   ) neiIdx = idx - sizeW;//下
				else if( c == 3 && y < sizeH-1 ) neiIdx = idx + sizeW;//上
				else if( c == 4 && z >   0     ) neiIdx = idx - sizeW*sizeH;
				else if( c == 5 && z < sizeD-1 ) neiIdx = idx + sizeW*sizeH;
				else continue;

				if( labels[ neiIdx ] > 0 && visFlg[ neiIdx ] == 0)
				{
					double d = colorDist2( &rgba[ idx*4  ], &rgba[ neiIdx*4 ] );
					if( d < colDiff )
					{ 
						colDiff = d; 
						labels[idx] = labels[ neiIdx ];
						visFlg[idx] = 1;
						changed = true;
					}
				}
			}
		}
		if( !changed ) break;
	}
}



















void TWaterShed_gradOfRgbaImage3D_forLargeVol(	int sizeW, 
												int sizeH, 
												int sizeD, 
												const byte* rgba , std::vector<int> &labels,
												void (*func_setProgressPos)(const double) )
{
	if( func_setProgressPos != 0 ) func_setProgressPos( 0 );

	//メモリの関係上 2 ^ 24 より大きいvolumeは一発でwater shadeできない.
	//depth方向に分割して、繰り返す
	int stepDepth = 0;
	int pW = powerOfTwo_val( sizeW );
	int pH = powerOfTwo_val( sizeH );
	int pD = powerOfTwo_val( sizeD );

	int maxDepthStep_p = 22 - pW - pH;//メモリがあれば23のほうが良い
	if( maxDepthStep_p < 0 ) {fprintf( stderr, "this volume is too big to compute water shed !\n"); return;}

	int maxDepthStep   = (int)ti_pow(2, maxDepthStep_p);

	if( pD <= maxDepthStep_p )
	{
		TWaterShed_gradOfRgbaImage3D( sizeW, sizeH, sizeD, rgba, labels );
		return;
	}

	int numOfIter = sizeD / maxDepthStep;
	if( sizeD % maxDepthStep > 0 ) numOfIter++;

	vector< vector<int> > segmentLavels( numOfIter, vector<int>() );
	
	for( int i = 0; i < numOfIter; ++i)
	{
		int pivDepth     =   i   * maxDepthStep;
		int pivDepthNext = (i+1) * maxDepthStep;
		if( pivDepthNext-1 <= sizeD-1)//indexにするには両方とも-1する必要がある
		{
			//踏み越えなし
			TWaterShed_gradOfRgbaImage3D( sizeW, sizeH, maxDepthStep, &rgba[ sizeW * sizeH * pivDepth * 4], segmentLavels[i] );
		}else{
			//踏み越え
			int lastSizeD = sizeD - pivDepth;
			TWaterShed_gradOfRgbaImage3D( sizeW, sizeH, lastSizeD, &rgba[ sizeW * sizeH * pivDepth * 4], segmentLavels[i] );
		}
		if( func_setProgressPos != 0 ) func_setProgressPos( (i+1)/ (double)numOfIter );
	}

	labels.resize( sizeW * sizeH * sizeD, 0);
	int labelOffset = 0;
	int idxOffset   = 0;
	for( int kk = 0; kk < numOfIter; ++kk)
	{
		int maxLabelVal = 0;
		//代入
		for( int i = 0; i < (int) segmentLavels[kk].size(); ++i) 
		{
			if( segmentLavels[kk][i] == 0 ) labels[ i + idxOffset ] = 0;
			else                            labels[ i + idxOffset ] = segmentLavels[kk][i] + labelOffset;
			maxLabelVal = max( maxLabelVal, segmentLavels[kk][i] );
		}
		labelOffset += maxLabelVal;
		idxOffset   += (int) segmentLavels[kk].size();
	}
}



void TWaterShed_gradOfRgbaImage3D(int sizeW, 
								  int sizeH, 
								  int sizeD, const byte* rgba, std::vector<int> &labels)
{
	//construct pixels and sort it///////////////////////////////////////////////
	const int size = sizeW * sizeH * sizeD;
	TWsPixel *pixels        = new TWsPixel [size];//
	TWsPixel **sortedPixPtr = new TWsPixel*[size];//sortしたもの
	for(int i = 0; i < size; ++i) sortedPixPtr[i] = &pixels[i];

	double a0, a1, a2, a3, a4, a5;
	int  zOffs = sizeW * sizeH;

	int _idx_ = 0; 
	for( int z = 0; z < sizeD; ++z)
	for( int y = 0; y < sizeH; ++y)
	for( int x = 0; x < sizeW; ++x, ++_idx_)
	{
		int idx = 4 * _idx_; //idx = 4 * ( z * sizeW * sizeH + y * sizeW + x );
		a0 = a1 = a2 = a3 = a4 = a5 = 0;
		
		if( x > 0       ) a0 = n_colorDist2( &rgba[idx], &rgba[idx - 4		 ]);
		if( y > 0       ) a2 = n_colorDist2( &rgba[idx], &rgba[idx - 4*sizeW ]);
		if( z > 0       ) a4 = n_colorDist2( &rgba[idx], &rgba[idx - 4*zOffs ]);

		if( x < sizeW-1 ) a1 = n_colorDist2( &rgba[idx], &rgba[idx + 4		 ]);
		if( y < sizeH-1 ) a3 = n_colorDist2( &rgba[idx], &rgba[idx + 4*sizeW ]);
		if( z < sizeD-1 ) a5 = n_colorDist2( &rgba[idx], &rgba[idx + 4*zOffs ]);


		double a = max( max(a0, a1), max(max(a2, a3),  max( a4, a5)));

		double gradMagnitude = sqrt( a ) / sqrt( 3.0 ) * 2;
		gradMagnitude *= 255; if( gradMagnitude > 255 ) gradMagnitude = 255;
		

		pixels[ _idx_ ].Set( (byte) gradMagnitude );
		for( int dz = -1; dz <=1 ; ++dz )
		for( int dy = -1; dy <=1 ; ++dy )
		for( int dx = -1; dx <=1 ; ++dx )if( dz != 0 || dy != 0 || dz != 0 )
		{
			if( z + dz >= 0 && z + dz < sizeD && 
				y + dy >= 0 && y + dy < sizeH && 
				x + dx >= 0 && x + dx < sizeW )
				
				pixels[_idx_].addNeighbour( &pixels[_idx_ + dx + dy * sizeW + dz * zOffs ] );
		}
	}

	qsort( sortedPixPtr, size, sizeof(TWsPixel*), cmpindex );

	//run watershed algorithm////////////////////////////////////////////////////
	runWatershedAlgorithm( size, pixels, sortedPixPtr, labels );

	delete[] sortedPixPtr;
	delete[] pixels;
}

