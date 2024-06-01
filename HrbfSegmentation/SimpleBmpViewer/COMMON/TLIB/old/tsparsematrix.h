#pragma once

#include "tmath.h"
#include <map>
#include <vector>
#include <list>
#include "TTriangleMesh.h"

#include "matrixlib/UmfpackHead.h"
#include "matrixlib/clapackHead.h"

//column 縦の列
//row    横の行
class TSparseMatrix
{
public:
	typedef list< pair<int, double> > MatRow;
	typedef vector< MatRow> Matrix;

	Matrix m_M;
	const int m_sizeRow;
	const int m_sizeCol;

	int*        m_Ap_c; // note that m_Ap_c[ RosSize ] is the non zero nomber
	int*        m_Ai_c;
	double*     m_Ax_c;

	~TSparseMatrix(void);
	
	//createNatrux (n*m) rowSize * colSize
	TSparseMatrix(int rowSize, int colSize);
	TSparseMatrix( const TSparseMatrix &src );
	void solveLinearEquationCG     (double* b, double* result );//Matrix should be symmertric
	void solveLinearEquationUmfPack(double* b, double* result );

	void *m_umfPackSymbolic, *m_umfPackNumeric;


	inline void push_back( int i, int j, double val)
	{ 
		//if( m_M[i].size() != 0 && m_M[i].back().first >= j ) fprintf(stderr, "strangeInput!!!!!!!!!!! %d %d\n", m_M[i].back().first, j);
		m_M[i].push_back( pair<int,double>( j, val)); 
	}
	

	int  calcMmultMt(TSparseMatrix &trgtMat) const;
	void createFieldForSolveLinearSystem( const int n_nonZeroEntry);

	void multVec     ( const int v_size, const double* v, const int trgt_size, double* trgt) const;
	void multVecToAdd( const int v_size, const double* v, const int trgt_size, double* trgt) const;


	void calcLaplacianTransMat            (const TTriangleMesh &mesh);
	void calcLaplacianTransMat_fixBoudnary(const TTriangleMesh &mesh, const vector<short> &bBoundary, const vector<vector<int>> &vert_1RingVs);
	void calcLaplacianTransMat(const   vector<TVector3> &vertices, const vector< vector< int > > &vert_1RingVs);
	void calcLaplacianTransMat(const int vSize, TVector3*vertices, const vector< vector< int > > &vert_1RingVs);
	
	//void calcInvertMatrix( TDenseMatrix &trgt );

	//Matrix for Implicit integration of simple spring network 
	void springModel_prepareEntries_3N_3N(const int vSize, const vector<TWingEdge> &edges);
	void springModel_addSimpleSpring_dfdx   (const double stiffK, const double* verts, const vector<TWingEdge> &edges, const double *e_restLen);
	void springModel_addSimpleLaplacian_dfdx(const double stiffK, const vector<vector<int>>  &vert_1RingVs);

	void umfPack_Prepare(){

		if( m_umfPackSymbolic != 0 ){
			umfpack_di_free_symbolic (&m_umfPackSymbolic);
			umfpack_di_free_numeric  (&m_umfPackNumeric );
		}
		if( m_sizeRow != m_sizeCol) {fprintf( stderr, "error row size is not equal to col size\n"); return;}
		if( m_Ai_c == 0           ) {fprintf( stderr, "error mat field is not ready!\n          "); return;}
		
		double *null = (double *) NULL ;
		int m = m_sizeCol , n = m_sizeCol ;
		(void) umfpack_di_symbolic( m, n, m_Ap_c, m_Ai_c, m_Ax_c, &m_umfPackSymbolic,                    null, null) ;
		(void) umfpack_di_numeric (       m_Ap_c, m_Ai_c, m_Ax_c,  m_umfPackSymbolic, &m_umfPackNumeric, null, null) ;
	}
	void umfPack_Solve_forPrecomp(double* b, double* result )
	{
		double *null = (double *) NULL ;
		if( m_sizeRow != m_sizeCol) {fprintf( stderr, "error row size is not equal to col size\n"); return;}
		if( m_Ai_c == 0           ) {fprintf( stderr, "error mat field is not ready!\n          "); return;}
		if( m_umfPackSymbolic == 0) umfPack_Prepare();
		(void) umfpack_di_solve(UMFPACK_A, m_Ap_c, m_Ai_c, m_Ax_c, result, b, m_umfPackNumeric,null, null) ;//UMFPACK_At　そもそも対称なので、必要ない
	}
	void umfPack_Release(){
		if( m_umfPackSymbolic != 0 ){
			umfpack_di_free_symbolic (&m_umfPackSymbolic);
			umfpack_di_free_numeric  (&m_umfPackNumeric );
		}
	}

	inline void clear()
	{ 
		m_M.clear(); 
		m_M.resize( m_sizeRow );
		if( m_Ap_c != 0 ) delete[] m_Ap_c;   m_Ap_c = 0;
		if( m_Ai_c != 0 ) delete[] m_Ai_c;   m_Ai_c = 0;
		if( m_Ax_c != 0 ) delete[] m_Ax_c;   m_Ax_c = 0;
		m_umfPackSymbolic = 0;
		m_umfPackNumeric  = 0;
	}


	inline void createCopy( TSparseMatrix &trgt )
	{
		if( m_sizeRow != trgt.m_sizeRow || m_sizeCol != trgt.m_sizeCol )
		{
			fprintf( stderr, "cant creatCopy!!\n" );
			return;
		}
		trgt.clear(); 
		for( int i = 0; i < m_sizeRow; ++i) trgt.m_M[i] = m_M[i];
		if( m_Ap_c != 0 )
		{
			int nonZeroEntry = m_Ap_c[m_sizeRow];
			trgt.m_Ap_c = new int   [ m_sizeRow + 1  ];//offset 1 が必要
			trgt.m_Ai_c = new int   [ nonZeroEntry   ];
			trgt.m_Ax_c = new double[ nonZeroEntry   ];
			memcpy( trgt.m_Ap_c, m_Ap_c, sizeof(  int ) * ( m_sizeRow + 1 ) );
			memcpy( trgt.m_Ai_c, m_Ai_c, sizeof(  int ) * nonZeroEntry      );
			memcpy( trgt.m_Ax_c, m_Ax_c, sizeof(double) * nonZeroEntry      );
		}
	}

	inline bool addToExistEngry( const int i, const int j, double val)
	{
		bool thereIsEntryIJ = false;
		for( MatRow::iterator it = m_M[i].begin(); it != m_M[i].end(); ++it) if( it->first == j)
		{
			it->second += val;
			thereIsEntryIJ = true;
			break;
		}
		if( m_Ap_c == 0 ) return thereIsEntryIJ;
		
		if( !thereIsEntryIJ ) return false;
		for( int idx = m_Ap_c[i]; idx < m_Ap_c[i+1]; ++idx) if( j == m_Ai_c[ idx ] )
		{
			m_Ax_c[ idx ] += val;
			return true;
		}
		fprintf( stderr, "cant add to existing entry!!\n" );
		return false;
	}
	
	inline void add33BlockToExistEngry_onlyArrayField( const int i, const int j, TMatrix9 &K){
		for( int idx = m_Ap_c[i  ]; idx < m_Ap_c[i+1]; ++idx) if( j== m_Ai_c[idx]){  m_Ax_c[ idx ] += K.data[0]; m_Ax_c[ idx+1] += K.data[3]; m_Ax_c[ idx+2] += K.data[6];}
		for( int idx = m_Ap_c[i+1]; idx < m_Ap_c[i+2]; ++idx) if( j== m_Ai_c[idx]){  m_Ax_c[ idx ] += K.data[1]; m_Ax_c[ idx+1] += K.data[4]; m_Ax_c[ idx+2] += K.data[7];}
		for( int idx = m_Ap_c[i+2]; idx < m_Ap_c[i+3]; ++idx) if( j== m_Ai_c[idx]){  m_Ax_c[ idx ] += K.data[2]; m_Ax_c[ idx+1] += K.data[5]; m_Ax_c[ idx+2] += K.data[8];}
	}
	inline void setAllEntryZero()
	{
		for( int i=0; i<m_sizeRow; ++i) 
		for( MatRow::iterator it = m_M[i].begin(); it != m_M[i].end(); ++it) it->second = 0;

		if( m_Ap_c == 0 ) return ;
		
		const int nonZeroSize = m_Ap_c[m_sizeRow];//
		memset( m_Ax_c, 0, sizeof( double ) * nonZeroSize );
	}

	inline void multScalar(double c)
	{
		if( m_Ax_c == 0 ) {
			for( int i=0; i<m_sizeRow; ++i) for( MatRow::iterator it = m_M[i].begin(); it != m_M[i].end(); ++it) it->second *= c;
		}
		else{
			for( int i=0; i<m_Ap_c[m_sizeRow]; ++i) m_Ax_c[i] *= c;
		}
	}


	inline void trace_vectorField()
	{
		fprintf( stderr, "\n\n!trace matrix row*col = %d %d\n", m_sizeRow, m_sizeCol );
		for( int i=0; i<m_sizeRow; ++i)
		{
			fprintf( stderr, "row[%d] ", i);
			for( MatRow::iterator it= m_M[i].begin(); it != m_M[i].end(); ++it)
				fprintf( stderr, "[%d]%.1f ", it->first, it->second);
			fprintf( stderr, "\n");
		}

	}
	inline void trace_pointerField()
	{
		fprintf( stderr, "\n\n!trace matrix row*col = %d %d\n", m_sizeRow, m_sizeCol );
		if( m_Ax_c == 0) return;
		for( int i=0; i<m_sizeRow; ++i)
		{
			fprintf( stderr, "row[%d] ", i);
			for( int k= m_Ap_c[i]; k < m_Ap_c[i+1]; ++k)
				fprintf( stderr, "[%d]%.1f, ", m_Ai_c[k], m_Ax_c[k]);
			fprintf( stderr, "\n");
		}
	}
};

inline void t_SolveLinearSystemUsingUmfpack( int N, int *Ap_c, int* Ai_c, double* Ax_c, const double* d, double* result,double threshold )
{
	double *null ;
	void *Symbolic, *Numeric ;

	null = (double *) NULL ;
	int m = N , n = N ;
	(void) umfpack_di_symbolic( m, n,   Ap_c, Ai_c, Ax_c, &Symbolic,           null, null) ;
	(void) umfpack_di_numeric (         Ap_c, Ai_c, Ax_c,  Symbolic, &Numeric, null, null) ;
	(void) umfpack_di_solve (UMFPACK_A, Ap_c, Ai_c, Ax_c, result, d,  Numeric, null, null) ;//UMFPACK_At　そもそも対称なので、必要ない
	umfpack_di_free_symbolic (&Symbolic) ;
	umfpack_di_free_numeric (&Numeric) ;
}




class TDenseMatrix
{
	int m_rowSize;//行数 (縦長)  transposeのときに変更される可能性あり
	int m_colSize;//列数 (横長)

public:

	double **m_data;

	TDenseMatrix(int rowSize, int colSize): m_rowSize( rowSize ), m_colSize( colSize ){
		fprintf(stderr, "TDenseMatrix construct %d %d\n", m_rowSize,m_colSize );
		m_data = new double*[m_rowSize];
		for( int y = 0; y<m_rowSize; ++y){
			m_data[y] = new double[ m_colSize ];
			memset( m_data[y], 0, sizeof( double ) * m_colSize );
		}
	}
	TDenseMatrix(const TDenseMatrix &src ): m_rowSize( src.m_rowSize), m_colSize( src.m_colSize ){
		fprintf(stderr, "TDenseMatrix construct %d %d\n", m_rowSize,m_colSize );
		m_data = new double*[m_rowSize];
		for( int y = 0; y<m_rowSize; ++y){
			m_data[y] = new double[ m_colSize ];
			memcpy( m_data[y], src.m_data[y], sizeof( double ) * m_colSize );
		}
	}

	~TDenseMatrix(){
		fprintf(stderr, "TDenseMatrix destract %d %d\n", m_rowSize,m_colSize );
		for( int y = 0; y<m_rowSize; ++y) delete[] m_data[y];
		delete[] m_data;
	}

	inline void setAllEntryZero(){
		for( int y = 0; y<m_rowSize; ++y)memset( m_data[y], 0, sizeof( double ) * m_colSize );
	}
	inline void copy( const TDenseMatrix &src){
		if(m_colSize != src.m_colSize || m_rowSize != src.m_rowSize ) return;
		for( int y = 0; y<m_rowSize; ++y) memcpy( m_data[y], src.m_data[y], sizeof( double ) * m_colSize );
	}
	inline void setIdentity(){
		for( int y = 0; y<m_rowSize; ++y) {
			memset( m_data[y], 0, sizeof( double ) * m_colSize );
			if( y < m_colSize )m_data[y][y] = 1;
		}
	}
	inline void transpose(){
		if( m_rowSize ==  m_colSize ){
			for( int y = 0  ; y<m_rowSize; ++y)
			for( int x = y+1; x<m_colSize; ++x) swap( m_data[x][y], m_data[y][x]); 
		}else{
			double **tmp = new double*[m_colSize]; 
			for( int i=0; i<m_colSize; ++i) tmp[i] = new double[ m_rowSize ];

			for( int y = 0  ; y<m_rowSize; ++y)
			for( int x = y+1; x<m_colSize; ++x) tmp[x][y] = m_data[y][x];
			swap( m_data, tmp );

			for( int i=0; i<m_rowSize; ++i) delete[] tmp[i];
			delete[] tmp;
			swap( m_rowSize, m_colSize);
		}
	}

	//r = Mv
	inline void multVector( const double *v, double *r )
	{
		memset( r, 0, sizeof( double ) * m_rowSize );
		int idx = 0;
		for( int y = 0; y < m_rowSize; ++y )
		for( int x = 0; x < m_colSize; ++x ) r[y] += m_data[y][x] * v[x]; 
	}

	inline void Trace(int offset = 1){
		fprintf( stderr, "matSize %d %d\n", m_rowSize, m_colSize);
		for( int y = 0; y < m_rowSize; y += offset ){
			for( int x = 0; x < m_colSize; x+= offset ) fprintf(stderr, "%.2f ",m_data[y][x] );
			fprintf( stderr,"\n");
		}
	}
	inline void Trace4(int offset = 1){
		fprintf( stderr, "matSize %d %d\n", m_rowSize, m_colSize);
		for( int y = 0; y < m_rowSize; y += offset ){
			for( int x = 0; x < m_colSize; x+= offset ) fprintf(stderr, "%.4f ",m_data[y][x] );
			fprintf( stderr,"\n");
		}
	}

	inline void TraceSymmetry(int ofset=1){
		fprintf( stderr, "matSize %d %d\n", m_rowSize, m_colSize);
		if( m_rowSize != m_colSize ) return;

		for( int y = 0; y < m_rowSize; y += ofset ){
			for( int x = 0; x < m_colSize; x += ofset ) 
				if     ( m_data[y][x] == 0            ) fprintf( stderr, " 0");
				else if( m_data[y][x] == m_data[x][y] ) fprintf( stderr, " s");
				else                                    fprintf( stderr, " w");
			fprintf( stderr,"\n");
		}
	}
	inline void TraceSign(int ofst=1){
		fprintf( stderr, "\n\nmatSize %d %d-------------------------\n", m_rowSize, m_colSize);
		if( m_rowSize != m_colSize ) return;

		for( int y = 0; y < m_rowSize; y+=ofst ){
			for( int x = 0; x < m_colSize; x+=ofst ) 
				if     ( fabs( m_data[y][x])  <= 0.000001 ) fprintf( stderr, "0");
				else if( m_data[y][x] <0   ) fprintf( stderr, "-");
				else                         fprintf( stderr, "+");
			fprintf( stderr,"\n");
		}
	}
	inline double*  operator[](int id) const { return m_data[id]; }

	void setMatMultMat( TDenseMatrix &L, TDenseMatrix &R);

	void QR_factorization( TDenseMatrix &Q, TDenseMatrix &R)const ;//Q(Orthogonal) R(UpTriangle) are already allocated 


	void LU_factorization_Full       (               const int *blockI, int *rowFlip, void (*progressFunc)(double) = 0 ); 
	void LU_factorization_RightBottom( const int n0, const int *blockI, int *rowFlip, void (*progressFunc)(double) = 0); 
	void LU_SolveLinearSystem        ( const int *rowFlip, const int *colFlip, const double *b, double *res);
	bool LU_flipColumn( const int idx, const bool allowRowFlip, bool &doRowFlip);
	bool LU_flipRow   ( const int idx, const bool allowColFlip, bool &doColFlip );



	void LU_matMult_dbg      ( const int *flipI, TDenseMatrix &trgt ); 
	void LU_matMult_dbg      ( const int *rowFlip, const int *colFlip, TDenseMatrix &trgt ); 
	void LDU_matMult_dbg     ( TDenseMatrix &trgt ); 
	void solveLinearSystem_LU( double *b, double *res ) ;


	void solveLinearSystem_QR     (       double *_b, double *res )const ;
	void solveLinearSystem_umfpack( const double *_b, double *res )const ;

	inline int getRowSize()const{return m_rowSize;}

};


/* -----------class TDenseMatLU------------------------------------------------------------------
// LU分解のための密行列クラス (Row majour)
//
// Pr:RowFlip / Pc:ColFlipは、行/列入れ替えを表し
// TDenseMatEx M を初期化後、LU分解すると
// Pr M Pc = LU となる
// M : 入力行列
// Pr:RowFlipMat(m_RowFlip) {0, 2, 1, 5, 4. 3}なら　0,2,1,5,4,3という順でrowが入るという意味
// Pc:ColFlipMat(m_ColFlip)
// LU: 下三角/上三角行列 m_dataに格納される
//
// !!!!!!!!!!!!イメージ!!!!!!!!! 
元行列 M を PrとPcでフリップした行列 M' = Pr M Pc が LU分解されて m_dataに入っている 
LU分解後の i行目は 元行列 M の m_RowFlip[i]行目 
//
//
// 1 : 単純なLU分解
// 2 : 前進代入・後進代入による線形方程式を解く
// 3 : 左上の n0 * n0 成分が分解済みの行列に対するLU分解
// 4 : Pr M Pc = LU と分解済みの下、
//     Pr( i, i+1 ) Pr M Pc' = L'U'というように、行i, i+1を入れ替えつつLU性を保つ行入れ替え
// 5 : 列入れ替え
// をサポートする
-----------class TDenseMatLU------------------------------------------------------------------*/
class TDenseMatEx
{
public:
	const int m_size;
	double **m_data;
	int     *m_RowFlip ;//Row flip Index
	int     *m_ColFlip ;//Col flip Index 

	TDenseMatEx(int size): m_size( size )
	{
		fprintf(stderr, "TDenseMatEx construct %d \n", size );
		m_RowFlip = new int    [m_size]; 
		m_ColFlip = new int    [m_size]; 
		m_data    = new double*[m_size];
		for( int i = 0; i<m_size; ++i){
			m_data[i] = new double[ m_size ];
			memset( m_data[i], 0, sizeof( double ) * m_size);
			m_RowFlip[i] = i; 
			m_ColFlip[i] = i;
		}
	}
	TDenseMatEx(const TDenseMatEx &src ): m_size( src.m_size)
	{
		fprintf(stderr, "TDenseMatEx construct %d\n", m_size, m_size  );
		m_RowFlip = new int    [m_size]; memcpy( m_RowFlip, src.m_RowFlip, sizeof( int  )*m_size);
		m_ColFlip = new int    [m_size]; memcpy( m_ColFlip, src.m_ColFlip, sizeof( int  )*m_size);
		m_data    = new double*[m_size];
		for( int i = 0; i<m_size; ++i) { m_data[i] = new double[ m_size ];
			                             memcpy( m_data[i], src.m_data[i], sizeof(double)*m_size );}
	}
	~TDenseMatEx(){
		fprintf(stderr, "TDenseMatEx destract %d\n", m_size);
		for( int y = 0; y<m_size; ++y) delete[] m_data[y];
		delete[] m_data   ;
		delete[] m_RowFlip;
		delete[] m_ColFlip;
	}

	inline void setSym(const int &i, const int &j, const double &v){ m_data[i][j]=m_data[j][i]=v; }
	inline void setAllEntryZero(){
		for( int i = 0; i<m_size; ++i) { m_RowFlip[i] = i; m_ColFlip[i] = i; memset( m_data[i], 0, sizeof(double)*m_size );}
	}
	inline void copy( const TDenseMatEx &src){
		if( m_size != src.m_size ) return;
		for( int y = 0; y<m_size; ++y) memcpy( m_data[y], src.m_data[y], sizeof(double)*m_size);
			                           memcpy( m_RowFlip, src.m_RowFlip, sizeof( int  )*m_size);
			                           memcpy( m_ColFlip, src.m_ColFlip, sizeof( int  )*m_size);
	}
	inline void copy( const TDenseMatEx &src, int N){ //N*N分だけコピー
		if( N > m_size || N > src.m_size ) return;
		memcpy( m_RowFlip, src.m_RowFlip, sizeof( int  )*N);
		memcpy( m_ColFlip, src.m_ColFlip, sizeof( int  )*N);
		for( int y = 0; y<N; ++y) memcpy( m_data[y], src.m_data[y], sizeof(double)*N);
	}

	inline void setIdentity(){
		for( int i = 0; i<m_size; ++i){
			m_RowFlip[i] = i;
			m_ColFlip[i] = i;
			memset( m_data[i], 0, sizeof( double ) * m_size );
			m_data[i][i] = 1;
		}
	}
	inline void Trace(int offset = 1){
		fprintf( stderr, "-----TDenseMatEx Trace %d %d\n", m_size, m_size);
		for( int y = 0; y < m_size; y += offset ){
			for( int x = 0; x < m_size; x+= offset ) fprintf(stderr, "%.2f ",m_data[y][x] );
			fprintf( stderr,"\n");
		}
	}
	inline void Trace4(int offset = 1){
		fprintf( stderr, "-----TDenseMatEx Trace %d %d\n", m_size, m_size);
		for( int y = 0; y < m_size; y += offset ){
			for( int x = 0; x < m_size; x+= offset ) fprintf(stderr, "%.4f ",m_data[y][x] );
			fprintf( stderr,"\n");
		}
	}
	inline void TraceSymmetry(int ofset=1){
		fprintf( stderr, "-----TDenseMatEx Trace %d %d\n", m_size, m_size);
		for( int y = 0; y < m_size; y += ofset ){
			for( int x = 0; x < m_size; x += ofset ) 
				if     ( m_data[y][x] == 0            ) fprintf( stderr, " 0");
				else if( m_data[y][x] == m_data[x][y] ) fprintf( stderr, " s");
				else                                    fprintf( stderr, " w");
			fprintf( stderr,"\n");
		}
	}
	inline void TraceSign(int ofst=1){
		fprintf( stderr, "-----TDenseMatEx Trace %d %d\n", m_size, m_size);
		for( int y = 0; y < m_size; y+=ofst ){
			for( int x = 0; x < m_size; x+=ofst ) 
				if     ( fabs( m_data[y][x])  <= 0.000001 ) fprintf( stderr, "0");
				else if( m_data[y][x] <0   ) fprintf( stderr, "-");
				else                         fprintf( stderr, "+");
			fprintf( stderr,"\n");
		}
	}
	inline int getSize()const{return m_size;}
	inline double*  operator[](int id) const { return m_data[id]; }

	void LU_factorization_Full       (               const int *blockI, void (*progressFunc)(double) = 0 ); 
	void LU_factorization_RightBottom( const int n0, const int *blockI, void (*progressFunc)(double) = 0 ); 
	void LU_SolveLinearSystem        ( const double *b, double *res                                      );
	bool LU_flipColumn               ( const int idx, const bool allowRowFlip);
	bool LU_flipRow                  ( const int idx, const bool allowColFlip);

	inline void LU_FactAndSolve( const int *blockI, const double *b, double *res, void (*progressFunc)(double) = 0 ){
		LU_factorization_Full( blockI, progressFunc);
		LU_SolveLinearSystem ( b, res );
	}


	/*
	void setMatMultMat( TDenseMatrix &L, TDenseMatrix &R);
	void LU_matMult_dbg      ( const int *flipI, TDenseMatrix &trgt ); 
	void LU_matMult_dbg      ( const int *rowFlip, const int *colFlip, TDenseMatrix &trgt ); 
	void LDU_matMult_dbg     ( TDenseMatrix &trgt ); 
	void solveLinearSystem_LU( double *b, double *res ) ;
	*/

	/* todo 
	inline void transpose(){
		if( m_rowSize ==  m_colSize ){
			for( int y = 0  ; y<m_rowSize; ++y)
			for( int x = y+1; x<m_colSize; ++x) swap( m_data[x][y], m_data[y][x]); 
		}else{
			double **tmp = new double*[m_colSize]; 
			for( int i=0; i<m_colSize; ++i) tmp[i] = new double[ m_rowSize ];

			for( int y = 0  ; y<m_rowSize; ++y)
			for( int x = y+1; x<m_colSize; ++x) tmp[x][y] = m_data[y][x];
			swap( m_data, tmp );

			for( int i=0; i<m_rowSize; ++i) delete[] tmp[i];
			delete[] tmp;
			swap( m_rowSize, m_colSize);
		}
	}
	inline void multVector( const double *v, double *r )
	{
		memset( r, 0, sizeof( double ) * m_rowSize );
		int idx = 0;
		for( int y = 0; y < m_rowSize; ++y )
		for( int x = 0; x < m_colSize; ++x ) r[y] += m_data[y][x] * v[x]; 
	}
	*/
};












//matrix column majour type
//support lu factoriszation with CLAPACK
//
//   1  2  3  4
//M= 5  6  7  8  ⇒  m_data = {1,5,8 2,6,10, 3,7,11, 4,8,12} でが入る
//   9 10 11 12
//
// M[y][x] = m_data[ x * m_Rsize + y ] 
class TDenseMatCM
{
	int     m_Csize;//列数 (横長)
	int     m_Rsize;//行数 (縦長)  transposeのときに変更される可能性あり
public:
	double *m_data ;
	
	inline void setAllEntryZero(){ memset(  m_data, 0, sizeof( double ) * m_Csize * m_Rsize); }
	inline void copy( const TDenseMatCM &src){
		if(m_Csize != src.m_Csize || m_Rsize != src.m_Rsize ) return;
		memcpy( m_data,  src.m_data,  sizeof(double)*m_Csize*m_Rsize );
	}

	TDenseMatCM(int rowSize, int colSize): m_Rsize( rowSize ), m_Csize( colSize ){
		fprintf(stderr, "TDenseMatCM construct %d %d\n", m_Rsize, m_Csize );
		m_data = new double[ m_Rsize * m_Csize ];
		setAllEntryZero();
	}

	TDenseMatCM(const TDenseMatCM &src ): m_Rsize( src.m_Rsize ), m_Csize( src.m_Csize ) {
		fprintf(stderr, "TDenseMatrix construct %d %d\n", m_Rsize, m_Csize );
		m_data = new double[ m_Csize * m_Rsize ];
		copy( src );
	}
	~TDenseMatCM(){
		fprintf(stderr, "TDenseMatrix destract %d %d\n", m_Rsize, m_Csize );
		delete[] m_data;
	}

	inline int getRowSize(){return m_Rsize;}
	inline double& operator[](int idx) const { return m_data[idx]; }
	inline double get   ( const int &y, const int &x                 ){ return m_data[ x*m_Rsize + y];}
	inline void   set   ( const int &y, const int &x, const double &v){ m_data[ x * m_Rsize + y ]  = v; }
	inline void   add   ( const int &y, const int &x, const double &v){ m_data[ x * m_Rsize + y ] += v; }
	inline void   setSym( const int &y, const int &x, const double &v){ m_data[ x * m_Rsize + y ]  =    m_data[ y * m_Rsize + x ]  = v; }
	inline void   addSym( const int &y, const int &x, const double &v){ m_data[ x * m_Rsize + y ] += v; m_data[ y * m_Rsize + x ] += v;}
	inline void setIdentity(){
		setAllEntryZero();
		for( int i = 0, s = min( m_Rsize, m_Csize ); i < s; ++i) set( i,i,1 );
	}
	inline void transpose(){
		if( m_Rsize ==  m_Csize ){
			for( int y = 0  ; y<m_Rsize; ++y)
			for( int x = y+1; x<m_Csize; ++x) swap( m_data[x*m_Rsize + y], m_data[ y*m_Rsize + x] ); 
		}else{
			double *tmp = new double[m_Csize*m_Rsize]; 
			for( int y = 0  ; y<m_Rsize; ++y)
			for( int x = y+1; x<m_Csize; ++x) tmp[y*m_Csize + x] = m_data[x * m_Rsize + y];
			swap( m_data , tmp    );
			swap( m_Rsize, m_Csize);
			delete[] tmp;
		}
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////
	//Trace matrix////////////////////////////////////////////////////////////////////////////////////
	inline void Trace2(int offset = 1){
		fprintf( stderr, "matSize %d %d\n", m_Rsize, m_Csize);
		for( int y = 0; y < m_Rsize; y += offset ){
			for( int x = 0; x < m_Csize; x+= offset ) fprintf(stderr, "%.2f ", m_data[x*m_Rsize + y] );
			fprintf( stderr,"\n");
		}
	}
	inline void Trace4(int offset = 1){
		fprintf( stderr, "matSize %d %d\n", m_Rsize, m_Csize);
		for( int y = 0; y < m_Rsize; y += offset ){
			for( int x = 0; x < m_Csize; x+= offset ) fprintf(stderr, "%.4f ", m_data[x*m_Rsize + y] );
			fprintf( stderr,"\n");
		}
	}
	inline void TraceSymmetry(int ofset=1){
		fprintf( stderr, "matSize %d %d\n", m_Rsize, m_Csize); 
		if( m_Rsize != m_Csize ) return;

		for( int y = 0; y < m_Rsize; y += ofset ){
			for( int x = 0; x < m_Csize; x += ofset ) 
				if     ( m_data[x*m_Rsize + y] == 0                     ) fprintf( stderr, " 0");
				else if( m_data[x*m_Rsize + y] == m_data[y*m_Rsize + x] ) fprintf( stderr, " s");
				else                                                      fprintf( stderr, " w");
			fprintf( stderr,"\n");
		}
	}
	inline void TraceSign(int ofst=1){
		fprintf( stderr, "matSize %d %d\n", m_Rsize, m_Csize);
		for( int y = 0; y < m_Rsize; y+=ofst ){
			for( int x = 0; x < m_Csize; x+=ofst ) 
				if     ( fabs( m_data[x*m_Rsize + y])  <= 0.000001 ) fprintf( stderr, "0");
				else if(       m_data[x*m_Rsize + y]   <0          ) fprintf( stderr, "-");
				else                                                 fprintf( stderr, "+");
			fprintf( stderr,"\n");
		}
	}


	//rowExchange (rEx[i])は,row[i]とrow[rEx[r]]の入れ替えを表す
	//rowIndecesは、変換後のrow[i]はもとの行列のrowIndces[i]を指していることを表す
	void LU_fact       ( long* rowIdces                   );
	void LU_fact_lapack( long* rowIdces, long *rowExchange);//the matrix will be factorized into LU form (diagonal entlies of L are 1)  
	void solveLinearSystem_LU       ( long* rowIdces   , const double *b, double *res);
	void solveLinearSystem_LU_lapack( long* rowExchange, const double *b, double *res);
	
	static void testLapack();
	/*
	//r = Mv
	inline void multVector( const double *v, double *r );//todo Blas implementation

	void setMatMultMat( TDenseMatrix &L, TDenseMatrix &R);

	void QR_factorization( TDenseMatrix &Q, TDenseMatrix &R)const ;//Q(Orthogonal) R(UpTriangle) are already allocated 

	void LU_matMult_dbg      ( const int *flipI, TDenseMatrix &trgt ); 
	void LU_matMult_dbg      ( const int *rowFlip, const int *colFlip, TDenseMatrix &trgt ); 
	void LDU_matMult_dbg     ( TDenseMatrix &trgt ); 
	void solveLinearSystem_LU( double *b, double *res ) ;

	void solveLinearSystem_QR     (       double *_b, double *res )const ;
	void solveLinearSystem_umfpack( const double *_b, double *res )const ;

	inline int getRowSize()const{return m_rowSize;}
	*/
};