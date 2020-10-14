#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

//define scores
#define matchScore 4
#define mismatchScore -1
#define gapScore -5
//Define direction constants
#define PATH -1
#define NONE 0
#define UP 1
#define LEFT 2
#define DIAG 3

void readFiles(char* queryFile, char* subjectFile);
void similarityScore(int i, int j, int* scoreMatrix, int* tbMatrix);
void backtrack(int* tbMatrix, int* scoreMatrix, long int* finalScore, char* queryResultReverse, char* subjectResultReverse);
void printMatrix(int* matrix);
void printTracebackMatrix(int* matrix);
void printResults(long int finalScore, double time, int numThreads, char* qrr, char* srr);
int matchMismatchScore(int i, int j);
int max(int x, int y);
int min(int x, int y);
void initialize(int *scoreMatrix);
long long int nElement(int i);
void calcFirstDiagElement(int *i, int *si, int *sj);

int querySize = 0;
int subjectSize = 0;
char* query, * subject;

int main(int argc, char* argv[]) {
	if (argc != 4) {
		printf("Please enter in this format: needleW <query_file_name> <subject_file_name> <num_threads>\n");
		return 1;
	}
	char* queryFile = argv[1];
	char* subjectFile = argv[2];
	int thread_count = atoi(argv[3]);
	readFiles(queryFile, subjectFile);


	//increment to add in 1 row and column
	querySize++;
	subjectSize++;

	//allocate flattened score matrix
	int *scoreMatrix = malloc(querySize * subjectSize * sizeof(int));
	int *tbMatrix = malloc(querySize * subjectSize * sizeof(int));

	//initialize variables
	long int finalScore = 0;
	int numThreads = 0;
	int si, sj, ai, aj, nEle;
	int nDiag = querySize + subjectSize - 3;
	//temporary allocation of string
	char* queryResultReverse = malloc(querySize*2);
	char* subjectResultReverse = malloc(subjectSize*2);
	//initialize matrix first row and column
	initialize(scoreMatrix);

	double initialTime = omp_get_wtime();

	#pragma omp parallel num_threads(thread_count) default(none) shared(scoreMatrix, tbMatrix, subjectSize, querySize, numThreads, nDiag) private(nEle, si, sj, ai, aj)
	{
		numThreads = omp_get_num_threads();
		for (int i=1; i <= nDiag; ++i) {
			nEle = nElement(i);
			calcFirstDiagElement(&i, &si, &sj);
			#pragma omp for
			for (int j=1; j<= nEle; ++j) {
				ai = si- j + 1;
				aj = sj + j -1;
				similarityScore(ai, aj, scoreMatrix, tbMatrix);
			}
		}
	}
	backtrack(tbMatrix, scoreMatrix, &finalScore, queryResultReverse, subjectResultReverse);


	double finalTime = omp_get_wtime();
	double timeElapsed = finalTime-initialTime;
	printResults(finalScore, timeElapsed, numThreads, queryResultReverse, subjectResultReverse);
	//printMatrix(scoreMatrix);
	//printTracebackMatrix(tbMatrix);

}

long long int nElement(int i) {
    if (i < querySize && i < subjectSize) {
        //Number of elements in the diagonal is increasing
        return i;
    }
    else if (i < max(querySize, subjectSize)) {
        //Number of elements in the diagonal is stable
        int val = min(querySize, subjectSize);
        return val - 1;
    }
    else {
        //Number of elements in the diagonal is decreasing
        int val = min(querySize, subjectSize);
        return 2 * val - i + abs(querySize - subjectSize) - 2;
    }
}

void calcFirstDiagElement(int *i, int *si, int *sj) {
    // Calculate the first element of diagonal
    if (*i < subjectSize) {
        *si = *i;
        *sj = 1;
    } else {
        *si = subjectSize - 1;
        *sj = *i - subjectSize + 2;
    }
}

void similarityScore(int i, int j, int* scoreMatrix, int* tbMatrix) {
    int up, left, diag;

    int index = querySize * i + j;

    //Get element above
    up = scoreMatrix[index-querySize] + gapScore;

    //Get element on the left
    left = scoreMatrix[index-1] + gapScore;

    //Get element on the diagonal
    diag = scoreMatrix[index-querySize-1] + matchMismatchScore(i, j);

    //Calculates the maximum
    int max;
    int pred;
    if (diag > left) {
    	max = diag;
    	pred = DIAG;
    }
    else {
    	max = left;
    	pred = LEFT;
    }
    if (up > max) {
    	max = up;
    	pred = UP;
    }
    //Inserts the value in the similarity and traceback matrixes
    scoreMatrix[index] = max;
    tbMatrix[index] = pred;
}

void backtrack(int* tbMatrix, int* scoreMatrix, long int* finalScore, char* queryResultReverse, char* subjectResultReverse) {
    int predPos;
	int resultSize = 0;
	//start from bottom right corner
	int currPos = querySize*subjectSize-1;
    *finalScore = scoreMatrix[currPos];
    //backtrack from btm right corner to top left corner
    do {
        if (tbMatrix[currPos] == DIAG) { //diagonal
            predPos = currPos - querySize - 1;
            //record character
            queryResultReverse[resultSize] = query[(currPos%querySize)-1];
            subjectResultReverse[resultSize++] = subject[((currPos-1)/querySize)-1];
        }
        else if (tbMatrix[currPos] == UP) { //up
            predPos = currPos - querySize;
            //insert - at subject string
            queryResultReverse[resultSize] = '-';
            subjectResultReverse[resultSize++] = subject[((currPos-1)/querySize)-1];
        }
        else if (tbMatrix[currPos] == LEFT) { //left
            predPos = currPos - 1;
            //insert - at query string
            queryResultReverse[resultSize] = query[(currPos%querySize)-1];
            subjectResultReverse[resultSize++] = '-';
        }
        tbMatrix[currPos] *= PATH;
        currPos = predPos;

    } while (currPos > 0);
    //null terminating the strings
    queryResultReverse[resultSize] = '\0';
    subjectResultReverse[resultSize] = '\0';
}

void initialize(int* scoreMatrix) {
	int i,j;

	for (i=0; i<querySize; i++) {
		for (j=0; j<subjectSize; j++) {
			int index = querySize * i + j;
			if (i==0 && j==0) {
				scoreMatrix[index] = 0;
			}
			else if (i==0) {
				scoreMatrix[index] = scoreMatrix[index-1] + gapScore;
			}
			else if (j==0) {
				scoreMatrix[index] = scoreMatrix[index-querySize] + gapScore;
			}
			else {
				scoreMatrix[index] = 0;
			}
		}
	}
}

void printMatrix(int* matrix) {
    int i, j;
	printf("\nSimilarity Matrix:\n");
    for (i = 0; i < querySize; i++) { //Lines
        for (j = 0; j < subjectSize; j++) {
            printf("%d\t", matrix[subjectSize * i +j ]);
        }
        printf("\n");
    }
}

void printTracebackMatrix(int* matrix) {
    int i, j, index;
    for (i = 0; i < querySize; i++) { //Lines
        for (j = 0; j < subjectSize; j++) {
            index = querySize * i + j;
            if(matrix[index] < 0) {
                if (matrix[index] == -UP)
                    printf("U ");
                else if (matrix[index] == -LEFT)
                    printf("L ");
                else if (matrix[index] == -DIAG)
                    printf("D ");
                else
                	printf("- ");
            }
            else {
            	printf("- ");
            }
        }
        printf("\n");
    }
}

void printResults(long int finalScore, double time, int numThreads, char* qrr, char* srr) {
	//reverse both strings
	strrev(qrr);
	strrev(srr);
	printf("\n======================================\n");
	printf("PROGRAM FINISHED\n");
	printf("Analyzed query and subject string of %d\n", querySize-1);
	printf("1) FINAL SCORE: %ld\n", finalScore);
	printf("2) ALIGNMENT STRING SIZE: %d\n", strlen(qrr));
	printf("3) ALIGNMENT STRING:\n");
	printf("\t%s\n", qrr);
	printf("\t");
	for (int i=0; i<strlen(qrr); i++) {
		if ((qrr[i] == '-') | (srr[i] == '-')) {
			printf(" ");
		}
		else if (qrr[i] == srr[i]) {
			printf("|");
		}
		else {
			printf("*");
		}
	}
	printf("\n\t%s\n", srr);
	printf("4) TIME ELAPSED: %fs\n", time);
	printf("5) NUMBER OF THREADS USED: %d\n", numThreads);
	printf("======================================\n");
}

int matchMismatchScore(int i, int j) {
	if (subject[i-1] == query[j-1])
        return matchScore;
    else
        return mismatchScore;
}

int max(int x, int y) {
	if (x > y)
		return x;
	else
		return y;
}

int min(int x, int y) {
	if (x > y)
		return y;
	else
		return x;
}

void readFiles(char* queryFile, char* subjectFile) {
	FILE* qfp = fopen(queryFile, "r");
	FILE* sfp = fopen(subjectFile, "r");

	if (qfp) {
		fseek(qfp, 0, SEEK_END);
		querySize = ftell(qfp);
		fseek(qfp, 0, SEEK_SET);
		query = malloc(querySize);
		if (query) {
			fread(query, 1, querySize, qfp);
		}
		fclose(qfp);
	}

	if (sfp) {
		fseek(sfp, 0, SEEK_END);
		subjectSize = ftell(sfp);
		fseek(sfp, 0, SEEK_SET);
		subject = malloc(subjectSize);
		if (subject) {
			fread(subject, 1, subjectSize, sfp);
		}
		fclose(sfp);
	}
}

