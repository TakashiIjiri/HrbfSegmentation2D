#include "StdAfx.h"
#include "TSparseMatrix.h"

#include <algorithm>
#ifdef _DEBUG
#define new DEBUG_NEW
#endif

TSparseMatrix::~TSparseMatrix(void)
{
	if( m_Ap_c != 0) delete[] m_Ap_c;
	if( m_Ai_c != 0) delete[] m_Ai_c;
	if( m_Ax_c != 0) delete[] m_Ax_c;

	if( m_umfPackSymbolic != 0 ){
		umfpack_di_free_symbolic (&m_umfPackSymbolic);
		umfpack_di_free_numeric  (&m_umfPackNumeric );
	}
}

//createNatrux (n*m) rowSize * colSize
TSparseMatrix::TSparseMatrix (int rowSize, int colSize): m_sizeRow(rowSize), m_sizeCol(colSize)
{ 
	m_M.resize( m_sizeRow );  
	fprintf( stderr, "mat created %d %d\n", rowSize,colSize );
	m_Ap_c = 0;
	m_Ai_c = 0;
	m_Ax_c = 0;

	m_umfPackSymbolic    = 0;
	m_umfPackNumeric     = 0;
}


TSparseMatrix::TSparseMatrix( const TSparseMatrix &src ): m_sizeRow(src.m_sizeRow), m_sizeCol(src.m_sizeCol)
{
	fprintf( stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
	fprintf( stderr, "!!!!!!!!TO IMPLEMENT! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
	fprintf( stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
	fprintf( stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
	fprintf( stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
	fprintf( stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
	fprintf( stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
	fprintf( stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
	fprintf( stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
}

int TSparseMatrix::calcMmultMt(TSparseMatrix &trgtMat) const
{
	//check
	if( trgtMat.m_sizeRow != m_sizeRow ){fprintf( stderr, "error 111123\n"); return 0;}  
	trgtMat.clear();
	Matrix& trgt = trgtMat.m_M;
	int zeroEntryNum = 0;
	int nonZeroNum = 0;
	for( int i = 0; i < m_sizeRow ; ++i )
	for( int j = i; j < m_sizeRow ; ++j )
	{
		double d = 0;
		bool   bExistEntry = false;
		MatRow::const_iterator itI = m_M[i].begin();
		MatRow::const_iterator itJ = m_M[j].begin();
		while( itI != m_M[i].end() && itJ != m_M[j].end() ) 
		{
			if(      itI->first == itJ->first )
			{
				d += itI->second * itJ->second;
				bExistEntry = true;
				++itI; ++itJ;
			}
			else if( itI->first < itJ->first ) ++itI;
			else                               ++itJ;
		}
		
		if( !bExistEntry ) { zeroEntryNum++; continue;}
		trgt[i].push_back( pair< int, double >( j, d ) ); ++nonZeroNum;
		if( i == j) continue;
		trgt[j].push_back( pair< int, double >( i, d ) ); ++nonZeroNum;
	}
	return nonZeroNum;
}



void TSparseMatrix::createFieldForSolveLinearSystem( const int n_nonZeroEntry)
{
	if( m_Ap_c != 0 ) delete[] m_Ap_c; 
	if( m_Ai_c != 0 ) delete[] m_Ai_c; 
	if( m_Ax_c != 0 ) delete[] m_Ax_c; 
	m_Ap_c = new int   [ m_sizeRow + 1  ];//offset 1 が必要
	m_Ai_c = new int   [ n_nonZeroEntry ];
	m_Ax_c = new double[ n_nonZeroEntry ];

	int index = 0;
	for( int i = 0; i < m_sizeRow; ++i)
	{
		m_Ap_c[i] = index;
		for( MatRow::iterator it = m_M[i].begin(); it != m_M[i].end(); ++it )
		{
			m_Ai_c[index] = it->first;
			m_Ax_c[index] = it->second;
			++index;
		}
	}
	m_Ap_c[m_sizeRow] = index;
}



void TSparseMatrix::multVec( const int v_size, const double* v, const int trgt_size, double* trgt) const
{
	if( v_size != m_sizeCol || trgt_size != m_sizeRow ){ fprintf( stderr, "error 11123\n"); return ; }
	
	if( m_Ai_c == 0){ 
		for( int i = 0; i < (int) m_M.size(); ++i)
		{
			trgt[i] = 0;
			for( MatRow::const_iterator it = m_M[i].begin(); it != m_M[i].end(); ++it)  
				trgt[i] += v[ it->first ] * it->second;
		}
	}
	else
	{
		for( int i = 0; i < m_sizeRow; ++i)
		{
			trgt[i] = 0;
			for( int idx = m_Ap_c[i]; idx < m_Ap_c[i+1]; ++idx)
				trgt[i] += v[ m_Ai_c[idx] ] * m_Ax_c[idx];
		}
	}
}




void TSparseMatrix::multVecToAdd( const int v_size, const double* v, const int trgt_size, double* trgt) const
{
	if( v_size != m_sizeCol || trgt_size != m_sizeRow ){ fprintf( stderr, "error 11123\n"); return ; }
	
	if( m_Ai_c == 0){ 
		for( int i = 0; i < (int) m_M.size(); ++i)
		{
			for( MatRow::const_iterator it = m_M[i].begin(); it != m_M[i].end(); ++it)  
				trgt[i] += v[ it->first ] * it->second;
		}
	}
	else
	{
		for( int i = 0; i < m_sizeRow; ++i)
		{
			for( int idx = m_Ap_c[i]; idx < m_Ap_c[i+1]; ++idx)
				trgt[i] += v[ m_Ai_c[idx] ] * m_Ax_c[idx];
		}
	}
}



/*------------------------------------------
boundary腸点では, laplacianは計算せず xi = xi0を計算する
------------------------------------------*/
void TSparseMatrix::calcLaplacianTransMat_fixBoudnary(const TTriangleMesh       &mesh, 
													  const vector<short      > &bBoundary,
													  const vector<vector<int>> &vert_1RingVs)
{
	const int         vSize = mesh.getVnum();
	const TVector3 *verts = mesh.m_verts;
	if( m_sizeRow != vSize|| m_sizeCol != vSize) {fprintf( stderr, "error1211\n");}
	
	clear();
	for( int i = 0; i < vSize; ++i)
	{
		push_back( i, i, 1 );
		if( !bBoundary[i] )
		{
			double sizeInv = - 1.0 / (double) vert_1RingVs[i].size();
			for( int k=0; k<(int) vert_1RingVs[i].size(); ++k) push_back( vert_1RingVs[i][k], i, sizeInv);
		}
	}
}


void TSparseMatrix::calcLaplacianTransMat(const vector<TVector3> &verts, const vector< vector< int > > &vert_1RingVs)
{
	if( m_sizeRow != (int) verts.size() || m_sizeCol != (int) verts.size() ) {fprintf( stderr, "error1211\n");}
	
	clear();
	for( int i = 0; i < (int) verts.size(); ++i)
	{
		push_back( i, i, 1 );
		bool oneRingContainSelf = false; // one ringが自分自身を含んでいる場合がある
		for(int j=0; j<(int) vert_1RingVs[i].size(); ++j )if( vert_1RingVs[i][j] == i){oneRingContainSelf = true; break;}
	
		double sizeInv = oneRingContainSelf ? - 1.0 / (double)( vert_1RingVs[i].size() -1):
			                                  - 1.0 / (double)  vert_1RingVs[i].size();	
		
		for(int j=0; j<(int) vert_1RingVs[i].size(); ++j )if( vert_1RingVs[i][j] != i) push_back( vert_1RingVs[i][j], i, sizeInv);
	}
}

void TSparseMatrix::calcLaplacianTransMat(const int vSize, TVector3* verts, const vector< vector< int > > &vert_1RingVs)
{
	if( m_sizeRow != vSize || m_sizeCol != vSize ) {fprintf( stderr, "error1211\n");}
	
	clear();
	for( int i = 0; i < vSize; ++i)
	{
		push_back( i, i, 1 );
		bool oneRingContainSelf = false; // one ringが自分自身を含んでいる場合がある
		for(int j=0; j<(int) vert_1RingVs[i].size(); ++j )if( vert_1RingVs[i][j] == i){oneRingContainSelf = true; break;}
	
		double sizeInv = oneRingContainSelf ? - 1.0 / (double)( vert_1RingVs[i].size() -1):
			                                  - 1.0 / (double)  vert_1RingVs[i].size();	
		
		for(int j=0; j<(int) vert_1RingVs[i].size(); ++j )if( vert_1RingVs[i][j] != i) push_back( vert_1RingVs[i][j], i, sizeInv);
	}
}

//Ltを作る
void TSparseMatrix::calcLaplacianTransMat(const TTriangleMesh &mesh)
{
	calcLaplacianTransMat( mesh.getVnum(), mesh.m_verts, mesh.m_v_1ring);
}


//Matrix should be 対称 
void TSparseMatrix::solveLinearEquationCG     (double* b, double* result )
{
	if( m_sizeRow != m_sizeCol) {fprintf( stderr, "error row size is not equal to col size\n"); return;}
	if( m_Ai_c == 0           ) {fprintf( stderr, "error mat field is not ready!\n          "); return;}

	t_solveLinearEquationByCG( m_sizeRow, m_Ax_c, m_Ap_c,  m_Ai_c, b, result, 0);
}


void TSparseMatrix::solveLinearEquationUmfPack(double* b, double* result )
{
	t_SolveLinearSystemUsingUmfpack( m_sizeRow, m_Ap_c, m_Ai_c,  m_Ax_c, b, result, 0);
}
//
////そのうちちゃんとした実装をする!
//void TSparseMatrix::calcInvertMatrix( TDenseMatrix &trgt )
//{
//	if( m_sizeRow      !=      m_sizeCol ||
//		trgt.m_colSize != trgt.m_rowSize || 
//		trgt.m_colSize !=      m_sizeRow  ){ fprintf( stderr, "cant calc inversion!!\n"); return ;}
//	if( m_Ai_c == 0           ) {fprintf( stderr, "error mat field is not ready!\n          "); return;}
//
//	const int N = m_sizeCol;
//	double *b   = new double[ N ];
//	double *res = new double[ N ];
//
//
//	for( int i = 0; i < N; ++i)
//	{
//		memset( res, 0, sizeof( double )* N );
//		memset( b  , 0, sizeof( double )* N );
//		b[i] = 1;
//
//		t_SolveLinearSystemUsingUmfpack( N, m_Ap_c, m_Ai_c,  m_Ax_c, b, res, 0);
//
//		int numNonZero = 0;
//		for( int j = 0; j < N; ++j ) trgt.m_data[ j * N + i ] = res[j];
//	}
//	delete[] b;
//	delete[] res;
//}
//



//Matrix for Implicit integration of simple spring network 
/*
Initialize 3N * 3N matrix
and its entries.

If edge vi-vj exist, this function set Kij(3*3) = {0,0,0  0,0,0  0,0,0}
This function also sets Kii = {0,0,0  0,0,0  0,0,0}
*/
void TSparseMatrix::springModel_prepareEntries_3N_3N(const int vSize, const vector<TWingEdge> &edges)
{
	if( m_sizeCol != 3 * vSize ) return;
	if( m_sizeRow != 3 * vSize ) return;

	//ii
	for( int i=0; i<vSize; ++i){
		const int piv = 3*i;
		m_M[piv  ].push_back( make_pair(piv, 0)); m_M[piv  ].push_back( make_pair(piv+1, 0)); m_M[piv  ].push_back( make_pair(piv+2, 0));
		m_M[piv+1].push_back( make_pair(piv, 0)); m_M[piv+1].push_back( make_pair(piv+1, 0)); m_M[piv+1].push_back( make_pair(piv+2, 0));
		m_M[piv+2].push_back( make_pair(piv, 0)); m_M[piv+2].push_back( make_pair(piv+1, 0)); m_M[piv+2].push_back( make_pair(piv+2, 0));
	}

	//ij & ji
	for( int e=0; e<(int)edges.size(); ++e){
		const int pivI = 3 * edges[e].v[0];
		const int pivJ = 3 * edges[e].v[1];
		//Kij
		m_M[pivI  ].push_back( make_pair(pivJ, 0)); m_M[pivI  ].push_back( make_pair(pivJ+1, 0)); m_M[pivI  ].push_back( make_pair(pivJ+2, 0));
		m_M[pivI+1].push_back( make_pair(pivJ, 0)); m_M[pivI+1].push_back( make_pair(pivJ+1, 0)); m_M[pivI+1].push_back( make_pair(pivJ+2, 0));
		m_M[pivI+2].push_back( make_pair(pivJ, 0)); m_M[pivI+2].push_back( make_pair(pivJ+1, 0)); m_M[pivI+2].push_back( make_pair(pivJ+2, 0));
		//Kji
		m_M[pivJ  ].push_back( make_pair(pivI, 0)); m_M[pivJ  ].push_back( make_pair(pivI+1, 0)); m_M[pivJ  ].push_back( make_pair(pivI+2, 0));
		m_M[pivJ+1].push_back( make_pair(pivI, 0)); m_M[pivJ+1].push_back( make_pair(pivI+1, 0)); m_M[pivJ+1].push_back( make_pair(pivI+2, 0));
		m_M[pivJ+2].push_back( make_pair(pivI, 0)); m_M[pivJ+2].push_back( make_pair(pivI+1, 0)); m_M[pivJ+2].push_back( make_pair(pivI+2, 0));
	}


	class pairIntDoubleComp {
	public:
		bool operator()(const pair<int,double>& riLeft, const pair<int,double>& riRight) const {
			return riLeft.first < riRight.first; // ここでは、元の < の意味に合わせてあるが、異なる判断も可能。
		}
	};
	//convert map --> vector
	for( int i=0; i<m_sizeRow; ++i) m_M[i].sort( pairIntDoubleComp() );

	//debug用//
	//for( int i=0; i<m_sizeRow; ++i){
	//	for( MatRow::iterator it = m_M[i].begin(); it != m_M[i].end(); ++it) fprintf( stderr, "%d ", it->first);
	//	fprintf( stderr, "\n" );
	//}

	int nonZeroEntries = 9 * vSize + 9 * 2 * (int) edges.size();
	createFieldForSolveLinearSystem( nonZeroEntries );

}




void TSparseMatrix::springModel_addSimpleSpring_dfdx   (const double stiffK, const double *verts, const vector<TWingEdge> &edges, const double *e_restLen)
{
	TMatrix9 Kij, Kii;
	TVector3 vij;
	for( int i=0; i<(int)edges.size(); ++i)
	{
		const TWingEdge &e = edges[i];
		const int idxI = e.v[0]*3;
		const int idxJ = e.v[1]*3;

		vij.data[0] = verts[idxJ  ] - verts[idxI  ]; 
		vij.data[1] = verts[idxJ+1] - verts[idxI+1]; 
		vij.data[2] = verts[idxJ+2] - verts[idxI+2]; 
		double L = vij.Length();

		double c =stiffK*( 1 - e_restLen[i]/L);
		Kij.data[0] = c; Kij.data[3] = 0; Kij.data[6] = 0;
		Kij.data[1] = 0; Kij.data[4] = c; Kij.data[7] = 0;
		Kij.data[2] = 0; Kij.data[5] = 0; Kij.data[8] = c;
		
		Kij.addMultVectors_withCoef( stiffK *e_restLen[i]/(L*L*L), vij, vij );
	
		Kii.Set( Kij );
		Kii.multScholar( -1 );

		add33BlockToExistEngry_onlyArrayField( idxI, idxI, Kii);
		add33BlockToExistEngry_onlyArrayField( idxJ, idxJ, Kii);
		add33BlockToExistEngry_onlyArrayField( idxI, idxJ, Kij);
		add33BlockToExistEngry_onlyArrayField( idxJ, idxI, Kij);//Kij_t=Kij
	}
}



void TSparseMatrix::springModel_addSimpleLaplacian_dfdx(const double stiffK, const vector<vector<int>>  &vert_1RingVs)
{

	TMatrix9 Kij, Kii;
	TVector3 vij;
	for( int i=0; i<(int)vert_1RingVs.size(); ++i)
	{
		const vector< int > &nei = vert_1RingVs[i];
		const int I = 3*i;
		addToExistEngry( I  ,I  , -stiffK * nei.size() );
		addToExistEngry( I+1,I+1, -stiffK * nei.size() );
		addToExistEngry( I+2,I+2, -stiffK * nei.size() );
		for( int k=0; k<(int)nei.size(); ++k) if( nei[k] != i)
		{
			const int J = 3*nei[k];
			addToExistEngry( I  ,J  , stiffK );
			addToExistEngry( I+1,J+1, stiffK );
			addToExistEngry( I+2,J+2, stiffK );
		}
	}
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Dense matrix//////////////////////////////////////////////////////////////////////////////////////////////////////////

//QR factorization by household method//
//Q: 直交行列   (allocated)
//R: 上三角行列 (allocated)
//note it's not necessary to compute Q when sloving Ax = b but compute Q^t b (こっちの方が明らかに軽い)
void TDenseMatrix::QR_factorization( TDenseMatrix &Q, TDenseMatrix &R)const 
{
	if( m_rowSize != Q.m_rowSize || m_rowSize != Q.m_colSize) return;
	if( m_rowSize != R.m_rowSize || m_colSize != R.m_colSize) return;
	if( m_rowSize < m_colSize ) return;

	R.copy( *this );
	Q.setIdentity();

	double *v   = new double [m_rowSize];
	double *vtR = new double [m_colSize];
	double *vtQ = new double [m_rowSize];

	for( int k = 0; k < m_colSize; ++k)
	{
		// H = |Ix   0    |    R' = HR   Q' = HQ  
		//     |0   I-vv^t|  
		const int n    = m_rowSize - k;
		double lenV_sq = 0;
		for( int i=0; i<n; ++i) { v[i] = R.m_data[k+i][k]; lenV_sq += v[i]*v[i];}
		double lenV = sqrt( lenV_sq );

		v[0] -= R.m_data[k][k]>0 ? -lenV : +lenV;

		//vを正規化
		lenV_sq += -R.m_data[k][k]*R.m_data[k][k] + v[0]*v[0];
		double coef = 2 / lenV_sq;

		//R' = H * R/////////////////////////////////////////////////////////////////
		for( int x = 0; x < m_colSize; ++x){
			vtR[x] = 0;
			for( int y=0; y<n; ++y) vtR[x] += v[y] * R.m_data[y+k][x];
			vtR[x] *= coef;
		}
		for( int y = 0; y < n        ; ++y) 
		for( int x = 0; x < m_colSize; ++x) R.m_data[y+k][x] -=  v[y] * vtR[x];

		//Q' = H * Q/////////////////////////////////////////////////////////////////
		for( int x = 0; x < m_rowSize; ++x){
			vtQ[x] = 0;
			for( int y=0; y<n; ++y) vtQ[x] += v[y] * Q.m_data[y+k][x];
			vtQ[x] *= coef;
		}
		for( int y = 0; y < n        ; ++y) 
		for( int x = 0; x < m_rowSize; ++x) Q.m_data[y+k][x] -=  v[y] * vtQ[x];
	}

	//今求めたQは  Hn Hn-1 ... H0 A = R --> QA = R というものなので A = Q^t R
	Q.transpose();

	delete[] v  ;
	delete[] vtR;
	delete[] vtQ;
}


void TDenseMatrix::solveLinearSystem_QR( double *b  ,//size:m_rowSize
	                                     double *res //size:m_colSize
											  )const 
{	
	clock_t aa = clock();
	double  *v    = new double [m_rowSize];
	double **R    = new double*[m_colSize];//転置した形で入れる

	for( int x=0; x<m_colSize; ++x){
		R[x] = new double[m_rowSize];
		for( int y=0; y<m_rowSize; ++y) R[x][y] = m_data[y][x];
	}
	//QR factorization/////////////////
	for( int k = 0; k < m_colSize; ++k)
	{
		// H = |Ix   0    |    R' = HR   Q' = HQ  
		//     |0   I-vv^t|  
		double lenV_sq = 0;
		for( int y=k; y<m_rowSize; ++y) { v[y] = R[k][y]; lenV_sq += v[y]*v[y];}//縦アクセス
		double lenV = sqrt( lenV_sq );
		v[k] -= R[k][k]>0 ? -lenV : +lenV; //v = x-y

		//vを正規化
		lenV_sq += -R[k][k] * R[k][k] + v[k]*v[k];  
		double coef = 2 / lenV_sq;

		//b' = H * b/////////////////////////////////////////////////////////////////
		double vt_b = 0;
		for( int y=k; y<m_rowSize; ++y) vt_b += v[y] * b[y];
		vt_b *= coef;	
		for( int y=k; y<m_rowSize; ++y) b[y] -= vt_b * v[y];

		//R' = H * R/////////////////////////////////////////////////////////////////
	    if( m_colSize - k < 100 ){
			for(int x = k; x < m_colSize; ++x){
				double vtR = 0;
				for(int y = k; y < m_rowSize; ++y) vtR += v[y] *R[x][y];//縦アクセス
				vtR *= coef;
				for(int y = k; y < m_rowSize; ++y) R[x][y] -= v[y] * vtR;//縦アクセス
			}
		}
		else{
			#pragma omp parallel shared( R, v)
			{
			#pragma omp for 
				for(int x = k; x < m_colSize; ++x){
					double vtR = 0;
					for(int y = k; y < m_rowSize; ++y) vtR += v[y] *R[x][y];//縦アクセス
					vtR *= coef;
					for(int y = k; y < m_rowSize; ++y) R[x][y] -= v[y] * vtR;//縦アクセス
				}
			}
		}
	}

	//solve R x = Qtb 後退代入
	for( int y = m_colSize-1; y>=0; --y)
	{
		double val = b[y];
		for( int x = y+1; x < m_colSize; ++x) val -= R[x][y] * res[x];//横アクセスになっちゃった・・・
		res[y] = val / R[y][y];
	}

	for( int i=0; i<m_rowSize; ++i)delete[] R[i];
	delete[] R   ;
	delete[] v   ;
	clock_t bb = clock();

	fprintf( stderr, "\nQR factorization%.10f\n\n", (bb-aa)/(double)CLOCKS_PER_SEC);
}



/*
	
	//行列Rを横にしかアクセスしないようにしたもの
	clock_t aa = clock();
	double  *v    = new double [m_rowSize];
	double *vtR   = new double [m_colSize];
	double **R    = new double*[m_rowSize];

	for( int i=0; i<m_rowSize; ++i){
		R[i] = new double[m_colSize];
		memcpy( R[i], m_data[i], sizeof( double ) * m_colSize );
	}
	clock_t sum1=0, sum2=0;
	//QR factorization/////////////////
	for( int k = 0; k < m_colSize; ++k)
	{
		// H = |Ix   0    |    R' = HR   Q' = HQ  
		//     |0   I-vv^t|  
		double lenV_sq = 0;
		for( int y=k; y<m_rowSize; ++y) { v[y] = R[y][k]; lenV_sq += v[y]*v[y];}//縦アクセス
		double lenV = sqrt( lenV_sq );
		v[k] -= R[k][k]>0 ? -lenV : +lenV; //v = x-y

		//vを正規化
		lenV_sq += -R[k][k] * R[k][k] + v[k]*v[k];  
		double coef = 2 / lenV_sq;

	clock_t t0 = clock();
		//R' = H * R/////////////////////////////////////////////////////////////////
		memset( &vtR[k], 0, sizeof( double ) * (m_colSize-k));

		for( int y = k; y < m_rowSize; ++y) 
		for( int x = k; x < m_colSize; ++x) vtR[x] += v[y] * R[y][x];//横アクセス

		for( int x = k; x < m_colSize; ++x) vtR[x] *= coef;

	clock_t t1 = clock();

		for( int y = k; y < m_rowSize; ++y) 
		for( int x = k; x < m_colSize; ++x) R[y][x] -= v[y] * vtR[x];//横アクセス
	clock_t t2 = clock();

		//b' = H * b/////////////////////////////////////////////////////////////////
		double vt_b = 0;
		for( int y=k; y<m_rowSize; ++y) vt_b += v[y] * b[y];
		vt_b *= coef;	
		for( int y=k; y<m_rowSize; ++y) b[y] -= vt_b * v[y];
	sum1 += t1-t0;
	sum2 += t2-t1;		
	}

	//solve R x = Qtb 後退代入
	for( int y = m_colSize-1; y>=0; --y)
	{
		double val = b[y];
		for( int x = y+1; x < m_colSize; ++x) val -= R[y][x] * res[x];
		res[y] = val / R[y][y];
	}

	for( int i=0; i<m_rowSize; ++i)delete[] R[i];
	delete[] R   ;
	delete[] vtR ;
	delete[] v   ;
	clock_t bb = clock();

	fprintf( stderr, "QR factorization%.10f %.10f %.10f\n\n", sum1/(double)CLOCKS_PER_SEC, 
		                                                      sum2/(double)CLOCKS_PER_SEC,
															  (bb-aa)/(double)CLOCKS_PER_SEC);

	//もとの実装 遅い
	double  *v    = new double [m_rowSize];
	double *vtR   = new double [m_colSize];
	double **R    = new double*[m_rowSize];

	clock_t t0 = clock();
	for( int i=0; i<m_rowSize; ++i){
		R[i] = new double[m_colSize];
		memcpy( R[i], m_data[i], sizeof( double ) * m_colSize );
	}

	//QR factorization/////////////////
	for( int k = 0; k < m_colSize; ++k)
	{
		// H = |Ix   0    |    R' = HR   Q' = HQ  
		//     |0   I-vv^t|  
		const int n    = m_rowSize - k;
		double lenV_sq = 0;

		for( int i=0; i<n; ++i) { v[i] = R[k+i][k]; lenV_sq += v[i]*v[i];}
		double lenV = sqrt( lenV_sq );

		v[0] -= R[k][k]>0 ? -lenV : +lenV;

		//vを正規化
		lenV_sq += -R[k][k] * R[k][k] + v[0]*v[0];  
		double coef = 2 / lenV_sq;

		//R' = H * R/////////////////////////////////////////////////////////////////
		for( int x = 0; x < m_colSize; ++x){
			vtR[x] = 0;
			for( int y=0; y<n; ++y) vtR[x] += v[y] * R[y+k][x];
			vtR[x] *= coef;
		}
		for( int y = 0; y < n        ; ++y) 
		for( int x = 0; x < m_colSize; ++x) R[y+k][x] -=  v[y] * vtR[x];

		//b' = H * b/////////////////////////////////////////////////////////////////
		double vt_b = 0;
		for( int y=0; y<n; ++y) vt_b   += v[y] * b[y+k];
		vt_b *= coef;	
		for( int y=0; y<n; ++y) b[y+k] -= vt_b * v[y];
	}

	clock_t t1 = clock();
	//solve R x = Qtb 後退代入
	for( int y = m_colSize-1; y>=0; --y)
	{
		double val = b[y];
		for( int x = y+1; x < m_colSize; ++x) val -= R[y][x] * res[x];
		res[y] = val / R[y][y];
	}

	for( int i=0; i<m_rowSize; ++i)delete[] R[i];
	delete[] R   ;
	delete[] vtR ;
	delete[] v   ;

	clock_t t2 = clock();
	fprintf( stderr, "QR factorization%.10f %.10f %.10f\n\n", (t1-t0)/(double)CLOCKS_PER_SEC, 
		                                                      (t2-t1)/(double)CLOCKS_PER_SEC);
*/



void TDenseMatrix::setMatMultMat( TDenseMatrix &L, TDenseMatrix &R)
{
	if( L.m_colSize != R.m_rowSize ) return;
	if( L.m_rowSize != m_rowSize || R.m_colSize != m_colSize ) return;

	for(int y = 0; y<m_rowSize; ++y)
	for(int x = 0; x<m_colSize; ++x)
	{
		double v =0;
		for( int k=0; k < L.m_colSize; ++k) v += L.m_data[y][k] * R.m_data[k][x];
		m_data[y][x] = v;
	}
}


void TDenseMatrix::solveLinearSystem_umfpack( const double *_b, double *res ) const 
{
	fprintf( stderr, "solve by umfpack %d %d\n", m_colSize, m_rowSize);
	if( m_colSize != m_rowSize ) return;
	if( m_colSize < 1 ) return;
	
	const int N = m_colSize;

	double *Ax_c = new double[ N * N ];
	int    *Ai_c = new int   [ N * N ];
	int    *Ap_c = new int   [ N +1  ];

	for( int y=0; y<N; ++y) Ai_c[y] = y;

	for( int y=0; y<N; ++y){
		Ap_c[y] = y * N;
		memcpy( &Ax_c[y*N], m_data[y], sizeof( double ) * N);
		memcpy( &Ai_c[y*N], &Ai_c [0], sizeof( int    ) * N);
	}
	Ap_c[N] = N*N;

	double *null ;
	void *Symbolic, *Numeric ;

	null = (double *) NULL ;
	int m = N , n = N ;
	(void) umfpack_di_symbolic( m, n,   Ap_c, Ai_c, Ax_c, &Symbolic,           null, null) ;
	(void) umfpack_di_numeric (         Ap_c, Ai_c, Ax_c,  Symbolic, &Numeric, null, null) ;
	(void) umfpack_di_solve (UMFPACK_A, Ap_c, Ai_c, Ax_c, res, _b,  Numeric, null, null) ;//UMFPACK_At　そもそも対称なので、必要ない
	umfpack_di_free_symbolic (&Symbolic);
	umfpack_di_free_numeric (&Numeric) ;

}


/*--------------------------------------------------------
LU_factorization_Full( const int *blockI, int *rowFlip)

行列をLU分解する
blockI  : block[i]の範囲でのみrow flipを行う
          もし0なら全体がflipの対象になる

rowFlip : 選択された軸が入る

--------------------------------------------------------*/
void TDenseMatrix::LU_factorization_Full( const int *blockI, int *rowFlip, void (*progressFunc)(double) )
{
	if( m_rowSize != m_colSize ) return;
	clock_t t0 = clock();
	
	const int N = m_rowSize;
	for( int i=0;i<N; ++i) rowFlip[i] = i;

	double *Myi = new double[N];
	double *m;
	double big = 0, tmp; 
	int    piv = 0     ;

	for( int I = 0; I < N; ++I)
	{
#pragma omp parallel for
		for( int y=0; y < N; ++y) Myi[y] = m_data[y][I];//縦アクセスはこれ以上しない

		double v = 0;
		for( int y=0; y <=I; ++y) {//U and L
			v = Myi[y];
			for( int k = 0; k < y; ++k) v -= m_data[y][k] * Myi[k];
			Myi[y] = v;
		}
		big = fabs( Myi[I] );
		piv = I;

		for( int y=I+1; y < N; ++y) {//U and L
			v = Myi[y];
			m = m_data[y];
			for( int k = 0; k < I; ++k, ++m) v -= *m * Myi[k];
			Myi[y] = v;
			if( ( blockI != 0 && y <= blockI[I]) && (tmp = fabs( Myi[y] )) > big ){ big = tmp;  piv = y;  }
		}
		if( piv != I ){
			swap(  m_data[I],  m_data[piv] );
			swap( rowFlip[I], rowFlip[piv] );
			swap(   Myi  [I],     Myi[piv] );
		}
		double coef = 1.0 / Myi[I];

#pragma omp parallel for
		for( int y = I+1; y < N; ++y) Myi[y] *= coef;
#pragma omp parallel for
		for( int y = 0  ; y < N; ++y) m_data[y][I] = Myi[y];

		if( I%100 == 0 && progressFunc != 0) progressFunc( I/(double)N );
	}
	delete[] Myi;

	fprintf( stderr, "LU factorization %.10f  \n", (clock()-t0)/(double)CLOCKS_PER_SEC);
}



/*--------------------------------------------------------
LU_factorization_RightBottom( const int *blockI, int *rowFlip)

行列をLU分解する
ただし、n0 * n0の領域は分解済みとする

blockI  : block[i]の範囲でのみrow flipを行う
          もし0なら全体がflipの対象になる

rowFlip : 選択された軸が入る

--------------------------------------------------------*/
void TDenseMatrix::LU_factorization_RightBottom( const int n0, const int *blockI, int *rowFlip, void (*progressFunc)(double))
{
	clock_t t0 = clock();

	if( m_rowSize != m_colSize ) return;
	const int N = m_rowSize;

	//M21について (fN, 0)  (fN, fn)    M12について ( 0, 0)  (fN, fn)
	//M21について (fN, N)  ( N, fn)    M12について (fN, N)  ( N, fn)

	for( int I = 0; I < n0; ++I)
	{
		double coef = 1.0 / m_data[I][I];
#pragma omp parallel for
		for( int xy = n0; xy < N; ++xy)
		{
			double a = 0, b = 0;
			for( int k=0; k < I; ++k){ 
				a += m_data[xy][k] * m_data[k][ I ];
				b -= m_data[ I][k] * m_data[k][ xy];
			}
			m_data[xy][ I ] -= a; m_data[xy][ I ] *= coef;
			m_data[I ][xy ] += b;
		}
	}

	double big = 0, tmp; 
	int    piv = 0     ;
	for( int I = n0; I < N; ++I)
	{
		double diagU = 0;
		for( int k=0; k < I; ++k) diagU += m_data[I][k] * m_data[k][I];
		m_data[I][I] -= diagU;

		//縦方向の計算と最大値検索を同時に
		big = fabs( m_data[I][I] ); piv = I;
#pragma omp parallel for
		for( int y = I+1; y < N; ++y)
		{
			double v = 0;
			for( int k=0; k < I; ++k) v += m_data[y][k] * m_data[k][I];
			m_data[y][I] -= v;

			if( ( blockI != 0 && y <= blockI[I]) && (tmp = fabs( m_data[y][I] )) > big ){ big = tmp;  piv = y;  }
		}

		if( piv != I ){ //縦方向(行) swap!
			swap(  m_data[I],  m_data[piv] );
			swap( rowFlip[I], rowFlip[piv] );
		}

		double coef = 1.0 / m_data[I][I];

		//横方向計算 (U)
#pragma omp parallel for
		for( int x = I+1; x < N; ++x){
			m_data[x][I] *= coef;
			double v = 0;
			for( int k = 0; k < I; ++k) v -= m_data[I][k] * m_data[k][x];
			m_data[I][x] += v;
		}
		if( I%100 == 0 && progressFunc != 0) progressFunc( I/(double)N );
	}

	fprintf( stderr, "LU factorization %.10f todo todo todo 縦アクセス入れたほうが多分早い！ \n", (clock()-t0)/(double)CLOCKS_PER_SEC);
}



/*--------------------------------------------------------
LU_SolveLinearSystem( const int *rowFlip, const int *colFlip, const double *b, double res)

LU分解がなされた下で、A res = bを解く
--------------------------------------------------------*/
void TDenseMatrix::LU_SolveLinearSystem( const int *rowFlip, const int *colFlip, const double *b, double *res)
{
	clock_t t0 = clock();

	const int N = m_rowSize	;
	double *f_B = new double[N];//fliped B
	double *f_Y = new double[N];//fliped Y
	double *f_X = new double[N];//fliped Y
	for( int i=0; i<N;++i) f_B[i] = b[ rowFlip[i] ];

	//前進代入 L*flipY = flipB --> a = L-1 b
	for( int y = 0; y < N; ++y){
		double v = f_B[y];
		for( int k=0; k<y; ++k) v -= m_data[y][k] * f_Y[k];
		f_Y[y] = v;
	}
	//solve U flipX = flipY (L-1 * flipB) 後退代入
	for( int y = N-1; y>=0; --y){
		double val = f_Y[y];
		for( int x = y+1; x < N; ++x) val -= m_data[y][x] * f_X[x];
		f_X[y] = val / m_data[y][y];
	}
	for( int i=0; i<N;++i) res[ colFlip[i] ] = f_X[ i ];
	delete[] f_B;
	delete[] f_Y;
	delete[] f_X;
	fprintf( stderr, "forward/backward substitution %.10f sec\n", (clock()-t0)/(double)CLOCKS_PER_SEC);
}



bool TDenseMatrix::LU_flipColumn( const int idx, const bool allowRowFlip, bool &doRowFlip)
{
	// L = 1 0  U = a c 
	//     z 1      0 d
	const int N = m_colSize;
	if( idx >= N-1 ) { fprintf( stderr, "strange input at t_flipRow\n"); return false; }
	
	double a = m_data[idx  ][idx], c = m_data[idx  ][idx+1];
	double z = m_data[idx+1][idx], d = m_data[idx+1][idx+1];
	double v = d + z*c ;

	if( !allowRowFlip && fabs( c ) < 0.0000001 ) {fprintf( stderr, "error strange imput!!\n"); return false;}

	double entrySum1 = fabs(d/c) + 2;
	double entrySum2 = fabs( c/v) + 1 + fabs( d/v) + fabs( z);
	if( !allowRowFlip || ( fabs( c ) >= 0.00000001 && entrySum1 <= entrySum2) )
	{
		doRowFlip = false;
		for( int i = 0; i < idx; ++i) swap( m_data[ i ][idx], m_data[ i ][idx+1] ); //x > y つまりUswap
		
		double m21 = d/c,  k21 = -d/c;

		m_data[idx  ][idx] = c  ;   m_data[idx  ][idx+1] =   a   ;
		m_data[idx+1][idx] = v/c;   m_data[idx+1][idx+1] = -a*d/c;

#pragma omp parallel for
		for( int i = idx+2; i < N; ++i){
			//L32 left / right  U23up / bottom
			double ll =       m_data[i][idx] + m_data[i][idx+1] * m21;
			double lr =                        m_data[i][idx+1]      ;
			double uu =       m_data[idx][i]               ;
			double ub = k21 * m_data[idx][i] + m_data[idx+1][i];

			m_data[idx  ][i] = uu;//U23up
			m_data[idx+1][i] = ub;//U23Bottom
			m_data[i][idx  ] = ll;//L32left 
			m_data[i][idx+1] = lr;//U32right
		}
	}else{
		doRowFlip = true;
		for( int i = 0; i < idx; ++i) {
			swap( m_data[idx][ i ], m_data[idx+1][ i   ] ); //x < y つまりLswap
			swap( m_data[ i ][idx], m_data[ i   ][idx+1] ); //x > y つまりUswap
		}
		double m11 = c/v           ,   k11 = z              ;  
		double m21 = d/v,  m22 = -z,   k21 = d/v, k22 = -c/v;
		m_data[idx  ][idx] = v  ;   m_data[idx  ][idx+1] =   a*z ;
		m_data[idx+1][idx] = c/v;   m_data[idx+1][idx+1] = a*d/v;

#pragma omp parallel for
		for( int i = idx+2; i < N; ++i)
		{
			//L32 left / right  U23up / bottom
			double ll = m_data[i][idx] * m11 + m_data[i][idx+1] * m21;
			double lr = m_data[i][idx]       + m_data[i][idx+1] * m22;
			double uu = k11 * m_data[idx][i] +       m_data[idx+1][i];
			double ub = k21 * m_data[idx][i] + k22 * m_data[idx+1][i];
			m_data[idx  ][i] = uu;//U23up
			m_data[idx+1][i] = ub;//U23Bottom
			m_data[i][idx  ] = ll;//L32left 
			m_data[i][idx+1] = lr;//U32right
		}
	}
	return true;
}


//row[ idx ] <--> row[ idx+1 ] フリップしつつ L U 性を保つ
//LとUの形によっては，col[idx] <--> col[idx+1]というフリップが不可避
bool TDenseMatrix::LU_flipRow( const int idx, const bool allowColFlip, bool &doColFlip )
{
	// L = 1 0  U = a c 
	//     z 1      0 d
	const int N = m_rowSize;
	if( idx >= N-1 ) { fprintf( stderr, "strange input at t_flipRow\n"); return false; }
	
	double a = m_data[idx  ][idx], c = m_data[idx  ][idx+1];
	double z = m_data[idx+1][idx], d = m_data[idx+1][idx+1];
	double v = d + z*c ;

	if( !allowColFlip && fabs( z ) < 0.00000001 ) return false;
	
	double entrySum1 = fabs( z  ) + fabs( 1/z) + 1;
	double entrySum2 = fabs( c/v) + 1 + fabs( d/v) + fabs( z);

	if( !allowColFlip || (fabs( z ) >= 0.000001 && entrySum1 <= entrySum2 ) )
	{
		doColFlip = false;
		for( int i = 0; i < idx; ++i) swap( m_data[idx][i], m_data[idx+1][i] );
		
		double m11 = 1/z           ,  k11 = z             ;
		double             m22 = -z,            k22 = -1/z;

		m_data[idx  ][idx] = a*z;   m_data[idx  ][idx+1] =   v ;
		m_data[idx+1][idx] = 1/z;   m_data[idx+1][idx+1] = -d/z;

#pragma omp parallel for
		for( int i = idx+2; i < N; ++i)
		{
			//L32 left / right  U23up / bottom
			double ll = m_data[i][idx] * m11                     ;
			double lr = m_data[i][idx]       + m_data[i][idx+1] * m22;
			double uu = k11 * m_data[idx][i] +       m_data[idx+1][i];
			double ub =                        k22 * m_data[idx+1][i];
			m_data[idx  ][i] = uu;//U23up
			m_data[idx+1][i] = ub;//U23Bottom
			m_data[i][idx  ] = ll;//L32left 
			m_data[i][idx+1] = lr;//U32right
		}
	}else{
		doColFlip = true;
		for( int i = 0; i < idx; ++i) {
			swap( m_data[idx][ i ], m_data[idx+1][ i   ] );
			swap( m_data[ i ][idx], m_data[ i   ][idx+1] );
		}
		double m11 = c/v           ,   k11 = z              ;  
		double m21 = d/v,  m22 = -z,   k21 = d/v, k22 = -c/v;
		m_data[idx  ][idx] = v  ;   m_data[idx  ][idx+1] =  a*z ;
		m_data[idx+1][idx] = c/v;   m_data[idx+1][idx+1] = a*d/v;

#pragma omp parallel for
		for( int i = idx+2; i < N; ++i)
		{
			//L32 left / right  U23up / bottom
			double ll = m_data[i][idx] * m11 + m_data[i][idx+1] * m21;
			double lr = m_data[i][idx]       + m_data[i][idx+1] * m22;
			double uu = k11 * m_data[idx][i] +       m_data[idx+1][i];
			double ub = k21 * m_data[idx][i] + k22 * m_data[idx+1][i];
			m_data[idx  ][i] = uu;//U23up
			m_data[idx+1][i] = ub;//U23Bottom
			m_data[i][idx  ] = ll;//L32left 
			m_data[i][idx+1] = lr;//U32right
		}
	}
	return true;
}




void TDenseMatrix::LU_matMult_dbg  ( const int *flipI, TDenseMatrix &trgt )
{
	const int N = m_colSize;

	for( int y = 0; y < N; ++y)
	for( int x = 0; x < N; ++x)
	{
		double v = 0;
		if( x == y ) {
			for( int k = 0; k < x; ++k) v += m_data[y][k] * 
				                             m_data[k][x];
			v += m_data[y][x];
		}else if( x > y ) {
			for( int k = 0; k < y; ++k) v += m_data[y][k] * 
				                             m_data[k][x];
			v += m_data[y][x];
		}else{
			for( int k = 0; k <=x; ++k) v += m_data[y][k] * 
				                             m_data[k][x];
		}

		trgt[y][x] = v;
	}
	double **tmp = new double*[N];
	for( int i=0; i< N; ++i) tmp[ flipI[i] ] = trgt.m_data[ i ];
	swap( trgt.m_data, tmp );
	delete[] tmp;
}
void TDenseMatrix::LU_matMult_dbg( const int *rowFlip, const int *colFlip, TDenseMatrix &trgt )
{
	const int N = m_colSize;
	double **tmp = new double*[N];
	for( int i=0; i<N; ++i){
		tmp[i] = new double[ N ];
	}

	for( int y = 0; y < N; ++y)
	for( int x = 0; x < N; ++x)
	{
		double v = 0;
		if( x == y ) {
			for( int k = 0; k < x; ++k) v += m_data[y][k] * 
				                             m_data[k][x];
			v += m_data[y][x];
		}else if( x > y ) {
			for( int k = 0; k < y; ++k) v += m_data[y][k] * 
				                             m_data[k][x];
			v += m_data[y][x];
		}else{
			for( int k = 0; k <=x; ++k) v += m_data[y][k] * 
				                             m_data[k][x];
		}

		tmp[y][x] = v;
	}
	for( int y = 0; y < N; ++y)
	for( int x = 0; x < N; ++x){
		trgt.m_data[ rowFlip[y] ][ colFlip[x] ] = tmp[y][x];
	}
	for( int i=0; i< N; ++i) delete[] tmp[i];
	delete[] tmp;
}





void TDenseMatrix::LDU_matMult_dbg  ( TDenseMatrix &trgt )
{
	const int N = m_colSize;

	for( int y = 0; y < N; ++y)
	for( int x = 0; x < N; ++x)
	{
		double v = 0;
		for( int k = 0; k < N; ++k)
		{
			double l = (k == y) ? 1 : (k < y) ? m_data[y][k] : 0 ;
			double r = (k == x) ? 1 : (k < x) ? m_data[k][x] : 0 ;
			v += m_data[k][k] * l*r;
		}
		trgt[y][x] = v;
	}

}



void TDenseMatrix::solveLinearSystem_LU( double *b, double *res ) 
{
	if( m_colSize != m_rowSize) return;
	const int N = m_colSize ;
	int *flipI = new int[ N ];
	for( int i=0;i<N; ++i) flipI[i] = i;

	clock_t t0 = clock();

	double **A = m_data; 
	double big = 0, tmp; 
	int piv = 0;
	for( int k=0; k<N; ++k) if( (tmp = fabs( A[k][0] ) ) > big ) { big = tmp; piv = k; } 

	if( piv != 0 ){
		swap(     A[0],     A[piv] );
		swap( flipI[0], flipI[piv] );
	}
	for( int k = 1; k< N; ++k) A[k][0] /= A[0][0];

	for( int I = 1; I < N; ++I )
	{
		double diagU = 0;
		for( int k=0; k < I; ++k) diagU += A[I][k] * A[k][I];
		A[I][I] -= diagU;

		//縦方向の計算と最大値検索を同時に
		big = fabs( A[I][I] ); piv = I;
#pragma omp parallel for
		for( int y = I+1; y < N; ++y)
		{
			double v = 0;
			for( int k=0; k < I; ++k) v += A[y][k] * A[k][I];
			A[y][I] -= v;

			if( (tmp = fabs( A[y][I] )) > big ){ big = tmp;  piv = y;  }
		}

		//縦方向(行) swap!
		if( piv != I ){
			swap(     A[I],      A[piv] );
			swap( flipI[I],  flipI[piv] );
		}

		double coef = 1.0 / A[I][I];
#pragma omp parallel for
		for( int y = I+1; y < N; ++y) A[y][I] *= coef;

		//横方向計算 (U)
#pragma omp parallel for
		for( int x = I+1; x < N; ++x){
			double v = 0;
			for( int k = 0; k < I; ++k) v -= A[I][k] * A[k][x];
			A[I][x] += v;
		}
	}

	clock_t t1 = clock();

	//L U x = b

	double *f_B = new double[N];//fliped B
	double *f_Y = new double[N];//fliped Y
	for( int i=0; i<N;++i) f_B[i] = b[ flipI[i] ];

	//前進代入 L*flipY = flipB --> a = L-1 b

	for( int y = 0; y < N; ++y){
		double v = f_B[y];
		for( int k=0; k<y; ++k) v -= A[y][k] * f_Y[k];
		f_Y[y] = v;
	}

	//solve U flipX = flipY (L-1 * flipB) 後退代入
	for( int y = N-1; y>=0; --y){
		double val = f_Y[y];
		for( int x = y+1; x < N; ++x) val -= A[y][x] * res[x];
		res[y] = val / A[y][y];
	}

	clock_t t2 = clock();

	fprintf( stderr, "LU factorization %.10f  代入%.10f\n\n", (t1-t0)/(double)CLOCKS_PER_SEC,
		                                                      (t2-t1)/(double)CLOCKS_PER_SEC);
	delete[] flipI;
	delete[] f_B;
	delete[] f_Y;
}



//the matrix will be factorized into LU form (diagonal entlies of L are 1)  
void TDenseMatCM::LU_fact( long* rowIdces )
{
	clock_t t0 = clock();

	if( m_Rsize != m_Csize ){ fprintf( stderr, "LU_fact_lapac is only for cubic matrix\n"); return; } 
	const int N = m_Rsize;
	for( int i=0;i<N; ++i) rowIdces[i] = i;

	double big  = 0, tmp; 
	int    piv  = 0     ;

	for( int I = 0; I < N; ++I)
	{
		const int I_N = I*N;
		//M(y, I) --> m_data[ I*N+y ]にデータが入る

		double v = 0;
		for( int y=0; y <=I; ++y) {//U and L
			v = 0;
			for( int k = 0; k < y; ++k) v -= m_data[k*N + y] * m_data[I_N + k];
			m_data[I_N + y] += v;
		}

		big = fabs( m_data[I*N + I] );
		piv = I;

		for( int y=I+1; y < N; ++y) {//U and L
			v = 0;
			for( int k = 0; k < I; ++k) v -= m_data[k*N+y] * m_data[I_N + k];
			m_data[I_N+y] += v;

			if( (tmp = fabs( m_data[I_N+y] )) > big ){ big = tmp;  piv = y;  }
		}

		if( piv != I ){
			for( int x=0; x < N; ++x) swap( m_data[x*N+I], m_data[x*N+piv] );
			swap( rowIdces[I], rowIdces[piv] );
		}

		double coef = 1.0 / m_data[I_N + I];
		for( int y = I+1; y < N; ++y) m_data[I_N + y] *= coef;

		if( I%150 == 0 ) fprintf( stderr, "%d/%d\n", I, N);
	}

	clock_t t1 = clock();
	fprintf( stderr, "LU factorization %.10f  \n", (t1-t0)/(double)CLOCKS_PER_SEC);
}



void TDenseMatCM::solveLinearSystem_LU(long* rowIdces, const double *b, double *res)
{
	if( m_Rsize != m_Csize ){ fprintf( stderr, "solveLinearSystem_LU is only for cubic matrix\n"); return; } 
	
	const int N = m_Csize;
	double *Y = new double[ N ];
	for( int i=0; i<N;++i) Y[i] = b[ rowIdces[i] ];

	//前進代入
	for( int ite = 1; ite < N; ++ite){
		const int ite_1 = (ite-1)*N;
		for( int k=ite; k < N; ++k) Y[k] -= Y[ite-1] * m_data[ ite_1 + k ];//L[k][ite-1]
	}

	memcpy( res, Y, sizeof( double ) * N );
	res[N-1] /= m_data[ (N-1)*N + N-1 ];

	//後進代入
	for( int ite = N-2; ite >= 0; --ite){
		const int ite_1 = (ite+1)*N;
		for( int k=ite; k >=0 ; --k) res[k] -= res[ ite+1 ] * m_data[ ite_1 + k ];//U[k][ite+1]
		res[ite] /= m_data[ ite*N + ite ];
	}
}



/*********************************************************************************/
/*  int dgetrf_(integer *m, integer *n, doublereal *a, integer *lda,             */
/*              integer *ipiv, integer *info);                                   */
/*  一般行列Aを三角分解(LU分解)する関数，分解は置換Pを含んで行なわれる           */
/*  すなわち P*A=L*U (Lの対角要素は1)                                            */
/*-------------------------------------------------------------------------------*/
/*  引数   : *m    行列Aの行数                                                   */
/*         : *n    行列Aの列数                                                   */
/*         : *a    行列Aへのポインタ                                             */
/*         : *lda  行列Aの行数                                                   */
/*         : *ipiv ピボットの添字，詳しくは後述                                  */
/*         : *info 例外情報                                                      */
/*********************************************************************************/
//the matrix will be factorized into LU form (diagonal entlies of L are 1)  
void TDenseMatCM::LU_fact_lapack( long* rowIdces, long *rowExchange){
	clock_t t0 = clock();

	if( m_Rsize != m_Csize ){ fprintf( stderr, "LU_fact_lapac is only for cubic matrix\n"); return; } 
	integer N = m_Rsize, INFO;
	dgetrf_( &N, &N, m_data, &N, rowExchange, &INFO );
	for( int i=0; i< N; ++i) rowIdces[i] = i;
	for( int i=0; i< N; ++i) swap( rowIdces[i], rowIdces[ rowExchange[i] ] );
	clock_t t1 = clock();
	fprintf( stderr, "LU factorization %.10f  \n", (t1-t0)/(double)CLOCKS_PER_SEC);
}

/*********************************************************************************/
/* int dgetrs_(char *trans, integer *n, integer *nrhs, doublereal *a,            */
/*             integer *lda, integer *ipiv, doublereal *b, integer *ldb,         */
/*             integer *info);                                                   */
/*  LU分解済みのAについて線型の方程式系 AX = B の求解を行なう                    */
/*-------------------------------------------------------------------------------*/
/*  引数   : *trans 行列Aの転置の方程式系を求解する．                            */
/*                  転置なし: 'N'，転置: 'T'，共役転置: 'C'                      */
/*         : *n     行列Aの次元                                                  */
/*         : *nrhs  右辺の行列Bの列数                                            */
/*         : *a     行列Aへのポインタ                                            */
/*         : *lda   行列Aの行数                                                  */
/*         : *ipiv  ピボットの添字，分解のときに出力されたもの                   */
/*         : *b     行列Bへのポインタ                                            */
/*         : *ldb   行列Bの行数                                                  */
/*         : *info  例外情報                                                     */
/*********************************************************************************/
void TDenseMatCM::solveLinearSystem_LU_lapack( long* rowExchange, const double *b, double *res)
{
	clock_t t0 = clock();
	if( m_Rsize != m_Csize ){ fprintf( stderr, "LU_fact_lapac is only for cubic matrix\n"); return; } 
	integer N = m_Rsize, N1 = 1, INFO;
	memcpy( res, b, sizeof( double ) * N );
	dgetrs_( "N", &N, &N1, m_data, &N, rowExchange, res, &N, &INFO );
	clock_t t1 = clock();
	fprintf( stderr, "前進行進代入 %.10f  \n", (t1-t0)/(double)CLOCKS_PER_SEC);

}


void TDenseMatCM::testLapack()
{
	integer N = 4;
	integer INFO ;    
	long   iPiv[4 ];
	double a   [16] = { 1.0 , 5.0 ,  9.0,  4.0, 
		                2.0 , 11.0,  7.0, 14.0, 
						3.0 , 10.0,  6.0, 15.0, 
						13.0, 2.0 , 12.0,  1.0};
  
	dgetrf_( &N, &N, a, &N, iPiv, &INFO );
        
	fprintf(stderr, "\n ");

    for (int y= 0; y< 4; ++y ){ 
		for (int x= 0; x< 4; ++x ) fprintf(stderr, "a[%d][%d]=%f ",y,x, a[ 4 * x + y] );
        fprintf(stderr, "\n ");
	}

    for (int y= 0; y< 4; ++y ){ 
		fprintf( stderr, "iPiv[%d]=%d \n", y, iPiv[y]);
	}
}







/*--------------------------------------------------------
*LU_factorization_Full( const int *blockI, int *rowFlip)
*
*行列をLU分解する
*blockI       : block[i]の範囲でのみrow flipを行う
*               blockI = 0 なら全体がflipの対象になる
*progressFunc : 進捗を表示する関数ポインタ
*               null出なければ進捗に応じて[0,1]の値を入れる
*
--------------------------------------------------------*/
void TDenseMatEx::LU_factorization_Full( const int *blockI, void (*progressFunc)(double) )
{
	clock_t t0 = clock();

	const int N = m_size;
	for( int i=0;i<N; ++i) m_RowFlip[i] = i;

	double *Myi = new double[N];
	double *m;
	double big = 0, tmp; 
	int    piv = 0     ;

	for( int I = 0; I < N; ++I)
	{
#pragma omp parallel for
		for( int y=0; y < N; ++y) Myi[y] = m_data[y][I];//縦アクセスはこれ以上しない

		double v = 0;
		for( int y=0; y <=I; ++y) {//U and L
			v = Myi[y];
			for( int k = 0; k < y; ++k) v -= m_data[y][k] * Myi[k];
			Myi[y] = v;
		}
		big = fabs( Myi[I] );
		piv = I;

		for( int y=I+1; y < N; ++y) {//U and L
			v = Myi[y];
			m = m_data[y];
			for( int k = 0; k < I; ++k, ++m) v -= *m * Myi[k];
			Myi[y] = v;
			if( ( ( blockI != 0 && y <= blockI[I]) || blockI == 0) && (tmp = fabs( Myi[y] )) > big ){ big = tmp;  piv = y;  }
		}
		if( piv != I ){
			swap( m_data   [I], m_data   [piv] );
			swap( m_RowFlip[I], m_RowFlip[piv] );
			swap(     Myi  [I],       Myi[piv] );
		}
		double coef = 1.0 / Myi[I];

#pragma omp parallel for
		for( int y = I+1; y < N; ++y) Myi[y] *= coef;
#pragma omp parallel for
		for( int y = 0  ; y < N; ++y) m_data[y][I] = Myi[y];

		if( I%100 == 0 && progressFunc != 0) progressFunc( I/(double)N );
	}
	delete[] Myi;

	fprintf( stderr, "LU factorization %.10f  \n", (clock()-t0)/(double)CLOCKS_PER_SEC);
}



/*-------------------------------------------------------------------------
**LU_factorization_RightBottom( const int *blockI, int *rowFlip)
**
**ただし、n0 * n0の領域は分解済みの下、行列をLU分解する
**
**blockI       : block[i]の範囲でのみrow flipを行う
**                0なら全体がflipの対象になる
**progressFunc : 進捗を表示する関数ポインタ
**               null出なければ進捗に応じて[0,1]の値を入れる
**          n0  n1
**  M = n0 M11 M12  = 済 未
**      n1 M21 M22    未 未
--------------------------------------------------------*/
void TDenseMatEx::LU_factorization_RightBottom( const int n0, const int *blockI, void (*progressFunc)(double))
{
	clock_t t0 = clock();

	const int N = m_size;

	double *Myi = new double[N];
	for( int I = 0; I < n0; ++I)
	{
		double coef = 1.0 / m_data[I][I];

#pragma omp parallel for
		for( int y=0; y < N; ++y) Myi[y] = m_data[y][I];//縦アクセスはこれ以上しない

#pragma omp parallel for
		for( int xy = n0; xy < N; ++xy)
		{
			double a = 0;
			for( int k=0; k < I; ++k) a += m_data[xy][k] * Myi[k];
			m_data[xy][ I ] -= a   ; 
			m_data[xy][ I ] *= coef;
		}
	}

	double *m;
	double big = 0, tmp; 
	int    piv = 0     ;

	for( int I = n0; I < N; ++I)
	{
#pragma omp parallel for
		for( int y=0; y < N; ++y) Myi[y] = m_data[y][I];//縦アクセスはこれ以上しない

		double v = 0;
		for( int y=0; y <=I; ++y) {//U and L
			v = Myi[y];
			for( int k = 0; k < y; ++k) v -= m_data[y][k] * Myi[k];
			Myi[y] = v;
		}
		big = fabs( Myi[I] );
		piv = I;

		for( int y=I+1; y < N; ++y) {//U and L
			v = Myi[y];
			m = m_data[y];
			for( int k = 0; k < I; ++k, ++m) v -= *m * Myi[k];
			Myi[y] = v;
			if( ( blockI != 0 && y <= blockI[I]) && (tmp = fabs( Myi[y] )) > big ){ big = tmp;  piv = y;  }
		}
		if( piv != I ){
			swap( m_data   [I], m_data   [piv] );
			swap( m_RowFlip[I], m_RowFlip[piv] );
			swap(     Myi  [I],       Myi[piv] );
		}
		double coef = 1.0 / Myi[I];

#pragma omp parallel for
		for( int y = I+1; y < N; ++y) Myi[y] *= coef;
#pragma omp parallel for
		for( int y = 0  ; y < N; ++y) m_data[y][I] = Myi[y];

		if( I%100 == 0 && progressFunc != 0) progressFunc( I/(double)N );
	}
	delete[] Myi;

	fprintf( stderr, "LU factorization %.10f  \n", (clock()-t0)/(double)CLOCKS_PER_SEC);
}



/*--------------------------------------------------------
**TDenseMatEx::LU_SolveLinearSystem( const double *b, double *res)
**LU分解がなされた下で、A res = bを解く
**--------------------------------------------------------*/
void TDenseMatEx::LU_SolveLinearSystem( const double *b, double *res)
{
	clock_t t0 = clock();

	const int N = m_size	;
	double *f_B = new double[N];//fliped B
	double *f_Y = new double[N];//fliped Y
	double *f_X = new double[N];//fliped Y
	for( int i=0; i<N;++i) f_B[i] = b[ m_RowFlip[i] ];
	//前進代入 L*flipY = flipB --> a = L-1 b
	for( int y = 0; y < N; ++y){
		double v = f_B[y];
		for( int k=0; k<y; ++k) v -= m_data[y][k] * f_Y[k];
		f_Y[y] = v;
	}
	//solve U flipX = flipY (L-1 * flipB) 後退代入
	for( int y = N-1; y>=0; --y){
		double val = f_Y[y];
		for( int x = y+1; x < N; ++x) val -= m_data[y][x] * f_X[x];
		f_X[y] = val / m_data[y][y];
	}
	for( int i=0; i<N;++i) res[ m_ColFlip[i] ] = f_X[ i ];
	delete[] f_B;
	delete[] f_Y;
	delete[] f_X;
	fprintf( stderr, "forward/backward substitution %.10f sec\n", (clock()-t0)/(double)CLOCKS_PER_SEC);
}


/*--------------------------------------------------------
**TDenseMatEx::LU_flipColumn( const int idx, const bool allowRowFlip, bool &doRowFlip)
**
** M'= Pr M Pc = LU 分解がされた下で、LU性を損なわないように、idx列目とidx+1列目をflipする。
**　その時に idx行目とidx+1行目をflipした方が数値的に安定であればflipする．
**
**    allowRowFlip : trueでかつ行flipを行った方が安定であれば行flipする
**    allowRowFlip = falseの場合、ゼロ割により計算できない事がある
　　　　　　　　　　 その場合falsを返す
**--------------------------------------------------------*/
bool TDenseMatEx::LU_flipColumn( const int idx, const bool allowRowFlip )
{
	// L = 1 0  U = a c 
	//     z 1      0 d
	const int N = m_size;
	if( idx >= N-1 ) { fprintf( stderr, "strange input at t_flipRow\n"); return false; }
	
	double a = m_data[idx  ][idx], c = m_data[idx  ][idx+1];
	double z = m_data[idx+1][idx], d = m_data[idx+1][idx+1];
	double v = d + z*c ;

	if( !allowRowFlip && fabs( c ) < 0.0000001 ) {fprintf( stderr, "error strange imput!!\n"); return false;}

	double entrySum1 = fabs(d/c) + 2;
	double entrySum2 = fabs( c/v) + 1 + fabs( d/v) + fabs( z);
	if( !allowRowFlip || ( fabs( c ) >= 0.00000001 && entrySum1 <= entrySum2) )
	{
		for( int i = 0; i < idx; ++i) swap( m_data[ i ][idx], m_data[ i ][idx+1] ); //x > y つまりUswap
		
		double m21 = d/c,  k21 = -d/c;

		m_data[idx  ][idx] = c  ;   m_data[idx  ][idx+1] =   a   ;
		m_data[idx+1][idx] = v/c;   m_data[idx+1][idx+1] = -a*d/c;

#pragma omp parallel for
		for( int i = idx+2; i < N; ++i){
			//L32 left / right  U23up / bottom
			double ll =       m_data[i][idx] + m_data[i][idx+1] * m21;
			double lr =                        m_data[i][idx+1]      ;
			double uu =       m_data[idx][i]               ;
			double ub = k21 * m_data[idx][i] + m_data[idx+1][i];

			m_data[idx  ][i] = uu;//U23up
			m_data[idx+1][i] = ub;//U23Bottom
			m_data[i][idx  ] = ll;//L32left 
			m_data[i][idx+1] = lr;//U32right
		}
	}else{
		for( int i = 0; i < idx; ++i) {
			swap( m_data[idx][ i ], m_data[idx+1][ i   ] ); //x < y つまりLswap
			swap( m_data[ i ][idx], m_data[ i   ][idx+1] ); //x > y つまりUswap
		}
		double m11 = c/v           ,   k11 = z              ;  
		double m21 = d/v,  m22 = -z,   k21 = d/v, k22 = -c/v;
		m_data[idx  ][idx] = v  ;   m_data[idx  ][idx+1] =   a*z ;
		m_data[idx+1][idx] = c/v;   m_data[idx+1][idx+1] = a*d/v;

#pragma omp parallel for
		for( int i = idx+2; i < N; ++i)
		{
			//L32 left / right  U23up / bottom
			double ll = m_data[i][idx] * m11 + m_data[i][idx+1] * m21;
			double lr = m_data[i][idx]       + m_data[i][idx+1] * m22;
			double uu = k11 * m_data[idx][i] +       m_data[idx+1][i];
			double ub = k21 * m_data[idx][i] + k22 * m_data[idx+1][i];
			m_data[idx  ][i] = uu;//U23up
			m_data[idx+1][i] = ub;//U23Bottom
			m_data[i][idx  ] = ll;//L32left 
			m_data[i][idx+1] = lr;//U32right
		}
		swap( m_RowFlip[ idx ], m_RowFlip[ idx+1 ] ); 
	}
	swap( m_ColFlip[ idx ], m_ColFlip[ idx+1 ] );
	return true;
}

//上の行flip版
bool TDenseMatEx::LU_flipRow( const int idx, const bool allowColFlip )
{
	// L = 1 0  U = a c 
	//     z 1      0 d
	const int N = m_size;
	if( idx >= N-1 ) { fprintf( stderr, "strange input at t_flipRow\n"); return false; }
	
	double a = m_data[idx  ][idx], c = m_data[idx  ][idx+1];
	double z = m_data[idx+1][idx], d = m_data[idx+1][idx+1];
	double v = d + z*c ;

	if( !allowColFlip && fabs( z ) < 0.00000001 ) return false;
	
	double entrySum1 = fabs( z  ) + fabs( 1/z) + 1;
	double entrySum2 = fabs( c/v) + 1 + fabs( d/v) + fabs( z);

	if( !allowColFlip || (fabs( z ) >= 0.000001 && entrySum1 <= entrySum2 ) )
	{
		for( int i = 0; i < idx; ++i) swap( m_data[idx][i], m_data[idx+1][i] );
		
		double m11 = 1/z           ,  k11 = z             ;
		double             m22 = -z,            k22 = -1/z;

		m_data[idx  ][idx] = a*z;   m_data[idx  ][idx+1] =   v ;
		m_data[idx+1][idx] = 1/z;   m_data[idx+1][idx+1] = -d/z;

#pragma omp parallel for
		for( int i = idx+2; i < N; ++i)
		{
			//L32 left / right  U23up / bottom
			double ll = m_data[i][idx] * m11                     ;
			double lr = m_data[i][idx]       + m_data[i][idx+1] * m22;
			double uu = k11 * m_data[idx][i] +       m_data[idx+1][i];
			double ub =                        k22 * m_data[idx+1][i];
			m_data[idx  ][i] = uu;//U23up
			m_data[idx+1][i] = ub;//U23Bottom
			m_data[i][idx  ] = ll;//L32left 
			m_data[i][idx+1] = lr;//U32right
		}
	}else{
		for( int i = 0; i < idx; ++i) {
			swap( m_data[idx][ i ], m_data[idx+1][ i   ] );
			swap( m_data[ i ][idx], m_data[ i   ][idx+1] );
		}
		double m11 = c/v           ,   k11 = z              ;  
		double m21 = d/v,  m22 = -z,   k21 = d/v, k22 = -c/v;
		m_data[idx  ][idx] = v  ;   m_data[idx  ][idx+1] =  a*z ;
		m_data[idx+1][idx] = c/v;   m_data[idx+1][idx+1] = a*d/v;

#pragma omp parallel for
		for( int i = idx+2; i < N; ++i)
		{
			//L32 left / right  U23up / bottom
			double ll = m_data[i][idx] * m11 + m_data[i][idx+1] * m21;
			double lr = m_data[i][idx]       + m_data[i][idx+1] * m22;
			double uu = k11 * m_data[idx][i] +       m_data[idx+1][i];
			double ub = k21 * m_data[idx][i] + k22 * m_data[idx+1][i];
			m_data[idx  ][i] = uu;//U23up
			m_data[idx+1][i] = ub;//U23Bottom
			m_data[i][idx  ] = ll;//L32left 
			m_data[i][idx+1] = lr;//U32right
		}
		swap( m_ColFlip[ idx ], m_ColFlip[ idx+1 ] );
	}
	swap( m_RowFlip[ idx ], m_RowFlip[ idx+1 ] ); 
	return true;
}

