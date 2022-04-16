#include<bits/stdc++.h> 
#include<chrono>
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
    vector<vector<int>> C = normalMatMul(A,B);
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