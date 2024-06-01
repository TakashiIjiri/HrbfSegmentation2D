#pragma once

#include "tmath.h"
#include <map>
#include <vector>
#include <list>
#include "TTriangleMesh.h"


#include "matrixlib/UmfpackHead.h"
#include "matrixlib/clapackHead.h"



inline void t_SolveLinearSystemUsingUmfpack( int N, int *Ap_c, int* Ai_c, double* Ax_c, const double* d, double* result,double threshold )
{
	double *null ;
	void *Symbolic, *Numeric ;

	null = (double *) NULL ;
	int m = N , n = N ;
	(void) umfpack_di_symbolic( m, n,   Ap_c, Ai_c, Ax_c, &Symbolic,           null, null) ;
	(void) umfpack_di_numeric (         Ap_c, Ai_c, Ax_c,  Symbolic, &Numeric, null, null) ;
	(void) umfpack_di_solve (UMFPACK_A, Ap_c, Ai_c, Ax_c, result, d,  Numeric, null, null) ;//UMFPACK_At�@���������Ώ̂Ȃ̂ŁA�K�v�Ȃ�
	umfpack_di_free_symbolic( &Symbolic );
	umfpack_di_free_numeric ( &Numeric  );
}


//column �c�̗�
//row    ���̍s
/*------------------------------------------------------------
TSparseMatrix �a�s��N���X

�e�s���Ƃ� 
list< int(column Idx), double (value)> m_M
�������� (0�v�f�����ɂ͉����z�u���Ȃ�)

�a�s��`����쐬���A
createFieldForSolveLinearSystem
���ĂԂƁA���`���������������߂̔z��(pointer field)
	int      *m_Ap_c; // note that m_Ap_c[ RosSize ] is the non zero nomber
	int      *m_Ai_c;
	double   *m_Ax_c;
���\�z�����B���ꂪ�\�z����Ĉȍ~�́A
�@�@addToExistEngry
�Ȃǂ̗v�f�ɑ΂��鑀��́Apoinder field�ɑ΂��čs����B
--------------------------------------------------------------*/
class TSparseMatrix
{
public:
	typedef list< pair<int, double> > MatRow;
	typedef vector< MatRow>           Matrix;

	const int m_sizeRow;
	const int m_sizeCol;
	Matrix    m_M;
	int      *m_Ap_c, *m_Ai_c; // note that m_Ap_c[ RosSize ] is the non zero nomber
	double   *m_Ax_c;

	void *m_umfPackSymbolic, *m_umfPackNumeric;
	//contructors------------------------------------------------------------
	TSparseMatrix(int rowSize, int colSize );
	TSparseMatrix(const TSparseMatrix &src );
	~TSparseMatrix(void){ clear();}
	inline void clear(){
		m_M.clear(); 
		m_M.resize( m_sizeRow );
		if( m_Ap_c != 0 ) delete[] m_Ap_c;   m_Ap_c = 0;
		if( m_Ai_c != 0 ) delete[] m_Ai_c;   m_Ai_c = 0;
		if( m_Ax_c != 0 ) delete[] m_Ax_c;   m_Ax_c = 0;

		if( m_umfPackSymbolic != 0 ){
			umfpack_di_free_symbolic (&m_umfPackSymbolic);   m_umfPackSymbolic = 0;
			umfpack_di_free_numeric  (&m_umfPackNumeric );   m_umfPackNumeric  = 0;
		}
	}
	void createFieldForSolveLinearSystem( const int n_nonZeroEntry);

	//Laplacian Matrix creation------------------------------------------------------------
	void calcLaplacianTransMat            (const TTriangleMesh &mesh);
	void calcLaplacianTransMat_fixBoudnary(const TTriangleMesh &mesh, const vector<short> &bBoundary, const vector<vector<int>> &vert_1RingVs);
	void calcLaplacianTransMat(const   vector<TVector3> &vertices, const vector< vector< int > > &vert_1RingVs);
	void calcLaplacianTransMat(const int vSize, TVector3*vertices, const vector< vector< int > > &vert_1RingVs);
	//Matrix for Implicit integration of simple MASS-SPRING network------------------------
	void springModel_prepareEntries_3N_3N   (const int vSize, const vector<TWingEdge> &edges);
	void springModel_addSimpleSpring_dfdx   (const double stiffK, const double* verts, const vector<TWingEdge> &edges, const double *e_restLen);
	void springModel_addSimpleLaplacian_dfdx(const double stiffK, const vector<vector<int>>  &vert_1RingVs);
	//Basic Interface----------------------------------------------------------------------
	int  calcMmultMt(TSparseMatrix &trgtMat) const;
	void multVec     ( const int v_size, const double* v, const int trgt_size, double* trgt) const;
	void multVecToAdd( const int v_size, const double* v, const int trgt_size, double* trgt) const;
	void multScalar  ( double c);
	void createCopy( TSparseMatrix &trgt );
	//CG solver(Matrix should be symmertric)-----------------------------------------------
	void solveLinearEquationCG     (double* b, double* result ){
		if( m_sizeRow != m_sizeCol) {fprintf( stderr, "error row size  is not equal to col size\n"); return;}
		if( m_Ai_c == 0           ) {fprintf( stderr, "error mat field is not ready!\n          "); return;}
		t_solveLinearEquationByCG( m_sizeRow, m_Ax_c, m_Ap_c,  m_Ai_c, b, result, 0);
	}

	//Umfpack solver------------------------------------------------------------------------
	void solveLinearEquationUmfPack(double* b, double* result ){
		t_SolveLinearSystemUsingUmfpack( m_sizeRow, m_Ap_c, m_Ai_c,  m_Ax_c, b, result, 0);
	}
	void umfPack_Prepare(){
		if( m_sizeRow != m_sizeCol) {fprintf( stderr, "error row size is not equal to col size\n"); return;}
		if( m_Ai_c == 0           ) {fprintf( stderr, "error mat field is not ready!\n          "); return;}

		if( m_umfPackSymbolic != 0 ){
			umfpack_di_free_symbolic (&m_umfPackSymbolic);
			umfpack_di_free_numeric  (&m_umfPackNumeric );
		}
		m_umfPackSymbolic = 0;
		m_umfPackNumeric  = 0;
		double *null = (double *) NULL ;
		int m = m_sizeCol , n = m_sizeCol ;
		(void) umfpack_di_symbolic( m, n, m_Ap_c, m_Ai_c, m_Ax_c, &m_umfPackSymbolic,                    null, null) ;
		(void) umfpack_di_numeric (       m_Ap_c, m_Ai_c, m_Ax_c,  m_umfPackSymbolic, &m_umfPackNumeric, null, null) ;
	}
	void umfPack_Solve_forPrecomp(double* b, double* result )
	{
		if( m_sizeRow != m_sizeCol) {fprintf( stderr, "error row size is not equal to col size\n"); return;}
		if( m_Ai_c == 0           ) {fprintf( stderr, "error mat field is not ready!\n          "); return;}
		double *null = (double *) NULL ;
		if( m_umfPackSymbolic == 0) umfPack_Prepare();
		(void) umfpack_di_solve(UMFPACK_At, m_Ap_c, m_Ai_c, m_Ax_c, result, b, m_umfPackNumeric,null, null) ;//UMFPACK_At
	}
	void umfPack_Release(){
		if( m_umfPackSymbolic != 0 ){
			umfpack_di_free_symbolic (&m_umfPackSymbolic);
			umfpack_di_free_numeric  (&m_umfPackNumeric );
		}
	}

	//interface for matrix entry------------------------------------------------------------------
	inline void push_back( int i, int j, double val)
	{ 
		assert( m_M[i].size() == 0 || m_M[i].back().first < j );
		m_M[i].push_back( pair<int,double>( j, val)); 
	}

	inline bool addToExistEngry( const int i, const int j, double val)
	{
		bool bIJexist = false;
		for( MatRow::iterator it = m_M[i].begin(); it != m_M[i].end(); ++it) if( it->first == j){ it->second += val; bIJexist = true; break; }

		if( m_Ap_c == 0 ) return bIJexist;
		if( !bIJexist   ) return false;

		for( int idx = m_Ap_c[i]; idx < m_Ap_c[i+1]; ++idx) if( j == m_Ai_c[ idx ] ){ m_Ax_c[ idx ] += val; return true; }
		fprintf( stderr, "cant add to existing entry!!\n" );
		return false;
	}

	inline void add33BlockToExistEngry_onlyPointerField( const int i, const int j, TMatrix9 &K){
		for( int idx = m_Ap_c[i  ]; idx < m_Ap_c[i+1]; ++idx) if( j== m_Ai_c[idx]){  m_Ax_c[ idx ] += K.data[0]; m_Ax_c[ idx+1] += K.data[3]; m_Ax_c[ idx+2] += K.data[6];}
		for( int idx = m_Ap_c[i+1]; idx < m_Ap_c[i+2]; ++idx) if( j== m_Ai_c[idx]){  m_Ax_c[ idx ] += K.data[1]; m_Ax_c[ idx+1] += K.data[4]; m_Ax_c[ idx+2] += K.data[7];}
		for( int idx = m_Ap_c[i+2]; idx < m_Ap_c[i+3]; ++idx) if( j== m_Ai_c[idx]){  m_Ax_c[ idx ] += K.data[2]; m_Ax_c[ idx+1] += K.data[5]; m_Ax_c[ idx+2] += K.data[8];}
	}
	inline void setAllEntryZero(){
		for( int i=0; i<m_sizeRow; ++i) 
		for( MatRow::iterator it = m_M[i].begin(); it != m_M[i].end(); ++it) it->second = 0;

		if( m_Ap_c == 0 ) return ;
		
		const int nonZeroSize = m_Ap_c[m_sizeRow];//
		memset( m_Ax_c, 0, sizeof( double ) * nonZeroSize );
	}

	//Trace matrix--------------------------------------------------------------------------------------------------
	inline void TraceMat_vector ();
	inline void TraceMat_pointer();
	static void TEST();
};





/*------------------------------------------------------------
TDenseMatrix ���s����N���X
m_rowSize �~ m_colSize�̖��s���\��

���`������solver������ Ax = b
LU����/QR����(�[���t�s��ɗ��p��)/umfpack/clapack
--------------------------------------------------------------*/
class TDenseMatrix
{
	int m_rowSize;//�s�� (�c��)  transpose�̂Ƃ��ɕύX�����\������
	int m_colSize;//�� (����)
public:
	double **m_data;
	inline int getRowSize()const{return m_rowSize;}
	inline double*  operator[](int id) const { return m_data[id]; }

	//Constructors-------------------------------------------------------------------------------------------//
	TDenseMatrix(int rowSize, int colSize): m_rowSize( rowSize ), m_colSize( colSize ){
		fprintf(stderr, "TDenseMatrix construct %d %d\n", m_rowSize,m_colSize );
		m_data = new double*[m_rowSize];
		for( int y = 0; y<m_rowSize; ++y){ m_data[y] = new double[ m_colSize ]; memset( m_data[y], 0, sizeof(double) * m_colSize );}
	}
	TDenseMatrix(const TDenseMatrix &src ): m_rowSize( src.m_rowSize), m_colSize( src.m_colSize ){
		fprintf(stderr, "TDenseMatrix construct %d %d\n", m_rowSize,m_colSize );
		m_data = new double*[m_rowSize];
		for( int y = 0; y<m_rowSize; ++y){ m_data[y] = new double[ m_colSize ]; memcpy( m_data[y], src.m_data[y], sizeof( double ) * m_colSize ); }
	}
	~TDenseMatrix(){
		fprintf(stderr, "TDenseMatrix destract %d %d\n", m_rowSize,m_colSize );
		for( int y = 0; y<m_rowSize; ++y) delete[] m_data[y];
		delete[] m_data;
	}

	//Set values---------------------------------------------------------------------------------------------//
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
	inline void multVector( const double *v, double *r ){ //r = Mv
		memset( r, 0, sizeof( double ) * m_rowSize );
		int idx = 0;
		for( int y = 0; y < m_rowSize; ++y )
		for( int x = 0; x < m_colSize; ++x ) r[y] += m_data[y][x] * v[x]; 
	}
	void setMatMultMat( TDenseMatrix &L, TDenseMatrix &R){
		if( L.m_colSize != R.m_rowSize || L.m_rowSize != m_rowSize || R.m_colSize != m_colSize) return;
		for(int y = 0; y<m_rowSize; ++y)
		for(int x = 0; x<m_colSize; ++x){
			double v =0;
			for( int k=0; k < L.m_colSize; ++k) v += L.m_data[y][k] * R.m_data[k][x];
			m_data[y][x] = v;
		}
	}
	
	//Trace Matrix------------------------------------------------------------------------------------//
	inline void Trace        (int offset = 1 );
	inline void Trace4       (int offset = 1 );
	inline void TraceSymmetry(int ofset  = 1 );
	inline void TraceSign    (int ofst   = 1 );

	//LU-factorization--------------------------------------
	void LU_factorization_Full       (               const int *blockI, int *rowFlip, void (*progressFunc)(double) = 0 ); 
	void LU_factorization_RightBottom( const int n0, const int *blockI, int *rowFlip, void (*progressFunc)(double) = 0 ); 
	void LU_SolveLinearSystem        ( const int *rowFlip, const double *b, double *res);
	void LU_matMult_dbg              ( const int *rowFlip, TDenseMatrix &trgt ); 

	//QR-factorization-------------------------------------
	void QR_factorization( TDenseMatrix &Q, TDenseMatrix &R)const ;//Q(Orthogonal) R(UpTriangle) are already allocated 

	void solveLinearSystem_QR     (       double *b, double *res )const ; //solve Ax = b by QR-factorization  (b-���������N����)
	void solveLinearSystem_LU     ( const double *b, double *res )const ; //solve Ax = b by LU-factorization

	void solveLinearSystem_umfpack( const double *b, double *res )const ; //solve Ax = b by umfpack
	void clapack_SolveLinearSystem( long* rowExchange, const double *b, double *res);
	void clapack_LU_factorization ( long* rowFlip  , long *rowExchange);

	static void TEST();

};

/* -----------class TDenseMatLU------------------------------------------------------------------
// LU�����̂��߂̖��s��N���X (Row Majour) (Row/Column flipping�@�\�t��)
//
// Pr:RowFlip / Pc:ColFlip�́A�s/�����ւ���\��
// TDenseMatEx M ����������ALU���������
// Pr M Pc = LU �ƂȂ�
// M : ���͍s��
// Pr:RowFlipMat(m_RowFlip) {0, 2, 1, 5, 4. 3}�Ȃ�@0,2,1,5,4,3�Ƃ�������row������Ƃ����Ӗ�
// Pc:ColFlipMat(m_ColFlip)
// LU: ���O�p/��O�p�s�� m_data�Ɋi�[�����
//
// !!!!!!!!!!!!�C���[�W!!!!!!!!! 
���s�� M �� Pr��Pc�Ńt���b�v�����s�� M' = Pr M Pc �� LU��������� m_data�ɓ����Ă��� 
LU������� i�s�ڂ� ���s�� M �� m_RowFlip[i]�s�� 
//
//
// 1 : �P����LU����
// 2 : �O�i����E��i����ɂ����`������������
// 3 : ����� n0 * n0 �����������ς݂̍s��ɑ΂���LU����
// 4 : Pr M Pc = LU �ƕ����ς݂̉��A
//     Pr( i, i+1 ) Pr M Pc' = L'U'�Ƃ����悤�ɁA�si, i+1�����ւ���LU����ۂs����ւ�
// 5 : �����ւ�
// ���T�|�[�g����
-----------class TDenseMatLU------------------------------------------------------------------*/
class TDenseMatEx
{
public:
	const int m_size;
	double  **m_data;
	int      *m_RowFlip ;//Row flip Index
	int      *m_ColFlip ;//Col flip Index 

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
	inline int getSize()const{return m_size;}
	inline double*  operator[](int id) const { return m_data[id]; }

	inline void copy( const TDenseMatEx &src){
		if( m_size != src.m_size ) return;
		for( int y = 0; y<m_size; ++y) memcpy( m_data[y], src.m_data[y], sizeof(double)*m_size);
			                           memcpy( m_RowFlip, src.m_RowFlip, sizeof( int  )*m_size);
			                           memcpy( m_ColFlip, src.m_ColFlip, sizeof( int  )*m_size);
	}
	
	inline void setSym(const int &i, const int &j, const double &v){ m_data[i][j]=m_data[j][i]=v; }
	inline void setAllEntryZero(){
		for( int i = 0; i<m_size; ++i) { m_RowFlip[i] = i; m_ColFlip[i] = i; memset( m_data[i], 0, sizeof(double)*m_size );}
	}

	inline void copy( const TDenseMatEx &src, int N){ //N*N�������R�s�[
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

	inline void Trace        (int offset = 1);
	inline void Trace4       (int offset = 1);
	inline void TraceSymmetry(int offset = 1);
	inline void TraceSign    (int offset = 1);

	void LU_factorization_Full       (               const int *blockI, void (*progressFunc)(double) = 0 ); 
	void LU_factorization_RightBottom( const int n0, const int *blockI, void (*progressFunc)(double) = 0 ); 
	void LU_SolveLinearSystem        ( const double *b, double *res                                      );
	bool LU_flipColumn               ( const int idx, const bool allowRowFlip);
	bool LU_flipRow                  ( const int idx, const bool allowColFlip);

	inline void LU_FactAndSolve( const int *blockI, const double *b, double *res, void (*progressFunc)(double) = 0 ){
		LU_factorization_Full( blockI, progressFunc);
		LU_SolveLinearSystem ( b, res );
	}

	//TODO �|���Z���͖�����
};









