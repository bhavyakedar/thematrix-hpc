#include<bits/stdc++.h> 
#include<chrono>
#include<omp.h>
using namespace std::chrono;
using namespace std;

#define CLK CLOCK_MONOTONIC

struct timespec diff(struct timespec start, struct timespec end){
	struct timespec temp;
	if((end.tv_nsec-start.tv_nsec)<0){
		temp.tv_sec = end.tv_sec-start.tv_sec-1;
		temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
	}
	else{
		temp.tv_sec = end.tv_sec-start.tv_sec;
		temp.tv_nsec = end.tv_nsec-start.tv_nsec;
	}
	return temp;
}

vector<vector<int>> normalMatMul(vector<vector<int>> &A, vector<vector<int>> &B)
{
    vector<vector<int>> C(A.size(),vector<int>(B[0].size(),0));
    for(int i = 0; i < A.size(); i++)
        for(int j = 0; j < B[0].size(); j++)
            for(int k = 0; k < A[0].size(); k++)
                C[i][j] += A[i][k]*B[k][j];
    return C;
}

vector<vector<int>> addSubMat(vector<vector<int>> &A, vector<vector<int>> &B, int X)
{
    vector<vector<int>> C(A.size(),vector<int>(A[0].size(),0));
    if(X == 1)
    {
        for(int i = 0; i < A.size(); i++)
            for(int j = 0; j < A[0].size(); j++)
                C[i][j] = A[i][j]+B[i][j];
    }
    else if(X == 0)
    {
        for(int i = 0; i < A.size(); i++)
            for(int j = 0; j < A[0].size(); j++)
                C[i][j] = A[i][j]-B[i][j];    
    }
    return C;
}

vector<vector<int>> matBlock(vector<vector<int>> &M, int s1, int s2)
{
    int r = M.size();
    int c = M[0].size();
    vector<vector<int>> MB(r/2, vector<int>(c/2,0));
    for(int i = 0, i1 = s1; i < r/2, i1 < r/2+s1; i++, i1++)
    {
        for(int j = 0, j1 = s2; j < c/2, j1 < c/2+s2; j++, j1++)
        {
            MB[i][j] = M[i1][j1];
        }
    }
    return MB;
}

vector<vector<int>> concMat(vector<vector<int>> &M11, vector<vector<int>> &M12, vector<vector<int>> &M21, vector<vector<int>> &M22)
{
    vector<vector<int>> M(M11.size()+M21.size(), vector<int>(M11[0].size()+M12[0].size(), 0));
    int r = M.size();
    int c = M[0].size();
    for(int i = 0; i < M.size(); i++)
        for(int j = 0; j < M[0].size(); j++)
        {
            M[i][j] = (i < r/2) ? ((j < c/2) ? M11[i][j] : M12[i][j-c/2]) : ((j < c/2) ? M21[i-r/2][j] : M22[i-r/2][j-c/2]);
        }
    return M;
}

void truncateRow(vector<vector<int>> &M)
{
    M.pop_back();
}

void truncateCol(vector<vector<int>> &M)
{
    for(int i = 0; i < M.size(); i++)
        M[i].pop_back();
}

vector<vector<int>> strassen(vector<vector<int>> &A, vector<vector<int>> &B, int l)
{
    int rA = A.size();
    int cA = A[0].size();
    int rB = B.size();
    int cB = B[0].size();

    if(cA != rB)
    {
        cout << "Invalid" << endl;
        vector<vector<int>> C = {{0}};
        return  C;
    }

    if(rA <= 32 || cA <= 32 || cB <= 32) return normalMatMul(A,B);

    bool flag1 = false, flag2 = false;

    if(rA%2 == 1)
    {
        vector<int> rA_zero(cA, 0);
        A.push_back(rA_zero);
        rA++;
        flag1 = true;
    }

    if(cA%2 == 1)
    {
        for(int i = 0; i < rA; i++)
            A[i].push_back(0);
        cA++;
        vector<int> rB_zero(cB, 0);
        B.push_back(rB_zero);
        rB++;
    }

    if(cB%2 == 1)
    {
        for(int i = 0; i < rB; i++)
            B[i].push_back(0);
        cB++;
        flag2 = true;
    }

    vector<vector<int>> A11 = matBlock(A,0,0);
    vector<vector<int>> A12 = matBlock(A,0,cA/2);
    vector<vector<int>> A21 = matBlock(A,rA/2,0);
    vector<vector<int>> A22 = matBlock(A,rA/2,cA/2);

    vector<vector<int>> B11 = matBlock(B,0,0);
    vector<vector<int>> B12 = matBlock(B,0,cB/2);
    vector<vector<int>> B21 = matBlock(B,rB/2,0);
    vector<vector<int>> B22 = matBlock(B,rB/2,cB/2);

    vector<vector<int>> T1 = addSubMat(B12,B11,0);
    vector<vector<int>> T2 = addSubMat(B22,T1,0);
    vector<vector<int>> T3 = addSubMat(B22,B12,0);
    vector<vector<int>> T4 = addSubMat(B21,T2,0);

    vector<vector<int>> S1 = addSubMat(A21,A22,1);
    vector<vector<int>> S2 = addSubMat(S1,A11,0);
    vector<vector<int>> S3 = addSubMat(A11,A21,0);
    vector<vector<int>> S4 = addSubMat(A12,S2,0);

	vector<vector<int>> P1;
	vector<vector<int>> P2;
	vector<vector<int>> P3;
	vector<vector<int>> P4;
	vector<vector<int>> P5;
	vector<vector<int>> P6;
	vector<vector<int>> P7;

	if(l < 3)
	{
		#pragma omp parallel
		{
			#pragma omp single
			{
				#pragma omp task
					P1 = strassen(A11,B11,l+1);
				#pragma omp task
					P2 = strassen(A12,B21,l+1);
				#pragma omp task
					P3 = strassen(S1,T1,l+1);
				#pragma omp task
					P4 = strassen(S2,T2,l+1);
				#pragma omp task
					P5 = strassen(S3,T3,l+1);
				#pragma omp task
					P6 = strassen(S4,B22,l+1);
				#pragma omp task
					P7 = strassen(A22,T4,l+1);

			}
		}
	}
	
	else
	{
		P1 = strassen(A11,B11,l+1);
		P2 = strassen(A12,B21,l+1);
		P3 = strassen(S1,T1,l+1);
		P4 = strassen(S2,T2,l+1);
		P5 = strassen(S3,T3,l+1);
		P6 = strassen(S4,B22,l+1);
		P7 = strassen(A22,T4,l+1);
	}
	
	

    vector<vector<int>> U1 = addSubMat(P1,P2,1);
    vector<vector<int>> U2 = addSubMat(P1,P4,1);
    vector<vector<int>> U3 = addSubMat(U2,P5,1);
    vector<vector<int>> U4 = addSubMat(U3,P7,1);
    vector<vector<int>> U5 = addSubMat(U3,P3,1);
    vector<vector<int>> U6 = addSubMat(U2,P3,1);
    vector<vector<int>> U7 = addSubMat(U6,P6,1);

    vector<vector<int>> C11 = U1;
    vector<vector<int>> C12 = U7;
    vector<vector<int>> C21 = U4;
    vector<vector<int>> C22 = U5;

    vector<vector<int>> C = concMat(C11, C12, C21, C22);

    if(flag1)
        truncateRow(C);
    
    if(flag2)
        truncateCol(C);

    return C;
}

int main(int argc, char* argv[])
{
    struct timespec start_e2e, end_e2e, start_alg, end_alg, e2e, alg;
	clock_gettime(CLK, &start_e2e);
	if(argc < 3){
		printf( "Usage: %s n p \n", argv[0] );
		return -1;
	}
	int N=atoi(argv[1]);	/* size of input matrices */
	int P=atoi(argv[2]);	/* number of processors*/
	omp_set_num_threads(P);

    char *problem_name = "matrix_multiplication";
	char *approach_name = "brute_force";
//	char buffer[10];
//	FILE* inputFile;
	FILE* outputFile;
	//	inputFile = fopen(argv[3],"r");

	char outputFileName[50];		
	sprintf(outputFileName,"output/%s_%s_%s_%s_output.txt",problem_name,approach_name,argv[1],argv[2]);

    int rA = N, cA = N, rB = N, cB = N;
    vector<vector<int>> A(rA, vector<int>(cA, 0));
    vector<vector<int>> B(rB, vector<int>(cB, 0));

    for(int i = 0; i < rA; i++)
        for(int j = 0; j < cA; j++)
            A[i][j] = 1+ (rand() % 100) - 200;
        
    for(int i = 0; i < rB; i++)
        for(int j = 0; j < cB; j++)
            B[i][j] = 1+ (rand() % 100) - 200;

    clock_gettime(CLK, &start_alg);	/* Start the algo timer */
    vector<vector<int>> C = strassen(A,B,0);
    clock_gettime(CLK, &end_alg);	/* End the algo timer */
	/* Ensure that only the algorithm is present between these two
	   timers. Further, the whole algorithm should be present. */


	/* Should end before anything else (printing comes later) */
	clock_gettime(CLK, &end_e2e);
	e2e = diff(start_e2e, end_e2e);
	alg = diff(start_alg, end_alg);

	/* problem_name,approach_name,n,p,e2e_sec,e2e_nsec,alg_sec,alg_nsec
	   Change problem_name to whatever problem you've been assigned
	   Change approach_name to whatever approach has been assigned
	   p should be 0 for serial codes!! 
	 */
	printf("%s,%s,%d,%d,%d,%ld,%d,%ld\n", problem_name, approach_name, N, P, e2e.tv_sec, e2e.tv_nsec, alg.tv_sec, alg.tv_nsec);

}