#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

//define scores
#define matchScore 5
#define mismatchScore -1
#define gapScore -5

void readFiles(char* queryFile, char* subjectFile);
void similarityScore(int i, int j, int* scoreMatrix);
void backtrack(int* scoreMatrix, long int* finalScore, char* queryResultReverse, char* subjectResultReverse);
void printMatrix(int* matrix);
void printResults(long int finalScore, double time, char* qrr, char* srr);
int matchMismatchScore(int i, int j);
int max(int x, int y);
void initialize(int *scoreMatrix);


int querySize = 0;
int subjectSize = 0;

char* query, * subject;

int main(int argc, char* argv[]) {
	if (argc != 3) {
		printf("Please enter in this format: needleW <query_file_name> <subject_file_name>\n");
		return 1;
	}
	char* queryFile = argv[1];
	char* subjectFile = argv[2];
	readFiles(queryFile, subjectFile);

	//increment to add in 1 row and column
	querySize++;
	subjectSize++;

	//allocate flattened score matrix
	int *scoreMatrix = malloc(querySize * subjectSize * sizeof(int));

	//initialize variables
	long int finalScore = 0;
	//temporary allocation of string
	char* queryResultReverse = malloc(querySize*2);
	char* subjectResultReverse = malloc(subjectSize*2);
	//initialize matrix first row and column
	initialize(scoreMatrix);

	double start = omp_get_wtime();

	for (int i=1; i<querySize; i++) {
		for (int j=1; j<subjectSize; j++) {
			similarityScore(i, j, scoreMatrix);
		}
	}
	printMatrix(scoreMatrix);
	backtrack(scoreMatrix, &finalScore, queryResultReverse, subjectResultReverse);
	double stop = omp_get_wtime();
	double timeElapsed = stop-start;
	printResults(finalScore, timeElapsed, queryResultReverse, subjectResultReverse);
}

void similarityScore(int i, int j, int* scoreMatrix) {
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
    if (diag > left) {
    	max = diag;
    }
    else {
    	max = left;
    }
    if (up > max) {
    	max = up;
    }
    //Inserts the value in the similarity and traceback matrixes
    scoreMatrix[index] = max;
}

void backtrack(int* scoreMatrix, long int* finalScore, char* queryResultReverse, char* subjectResultReverse) {
    int predPos;
	int resultSize = 0;
	//start from bottom right corner
	int currPos = querySize*subjectSize-1;
    *finalScore = scoreMatrix[currPos];
    //backtrack from btm right corner to top left corner
    do {
        if (((scoreMatrix[currPos] - matchScore) == scoreMatrix[currPos-querySize-1]) | ((scoreMatrix[currPos] - mismatchScore == scoreMatrix[currPos-querySize-1]))) { //diagonal
            predPos = currPos - querySize - 1;
            //record character
            queryResultReverse[resultSize] = query[(currPos%querySize)-1];
            subjectResultReverse[resultSize++] = subject[((currPos-1)/querySize)-1];
        }
        else if (scoreMatrix[currPos] - gapScore == scoreMatrix[currPos-querySize]) { //up
            predPos = currPos - querySize;
            //insert - at subject string
            queryResultReverse[resultSize] = '-';
            subjectResultReverse[resultSize++] = subject[((currPos-1)/querySize)-1];
        }
        else if (scoreMatrix[currPos] -gapScore == scoreMatrix[currPos-1]) { //left
            predPos = currPos - 1;
            //insert - at query string
            queryResultReverse[resultSize] = query[(currPos%querySize)-1];
            subjectResultReverse[resultSize++] = '-';
        }

        currPos = (predPos < 0) ? 0 : predPos;

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

void printResults(long int finalScore, double time, char* qrr, char* srr) {
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

