#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>

#define PATH -1
#define NONE 0
#define UP 1
#define LEFT 2
#define DIAGONAL 3
//Define scores
#define matchScore 5
#define mismatchScore -3
#define gapScore -4

void readFiles(char* queryFile, char* subjectFile);
void similarityScore(int i, int j, int* scoreMatrix, int* tbMatrix, int* maxPos);
int matchMismatchScore(int i, int j);
void backtrack(int* tbMatrix, int* scoreMatrix, int maxPos, long int* finalScore, char* queryResultReverse, char* subjectResultReverse);
void printMatrix(int* matrix);
void printTracebackMatrix(int* tbMatrix);
void printResults(long int finalScore, double time, char* qrr, char* srr);

int querySize = 0;
int subjectSize = 0;

char* query, * subject;

int main(int argc, char* argv[]) {
	if (argc != 3) {
		printf("Please enter in this format: SmithW <query_file_name> <subject_file_name>\n");
		return 1;
	}
	char* queryFile = argv[1];
	char* subjectFile = argv[2];
	readFiles(queryFile, subjectFile);

	//increment to include 0s in the first row and column
	querySize++;
	subjectSize++;

	//allocate flattened score matrix and traceback matrix
	int *scoreMatrix = calloc(querySize * subjectSize, sizeof(int));
	int *tbMatrix = calloc(querySize * subjectSize, sizeof(int));

	//initialize variables
	long int finalScore = 0;
	//temporary allocation of string
	char* queryResultReverse = malloc(querySize);
	char* subjectResultReverse = malloc(subjectSize);
	int maxPosition = 0;

	//start clock
	clock_t tic = clock();

	for (int i=1; i<querySize; i++) {
		for (int j=1; j<subjectSize; j++) {
			similarityScore(i, j, scoreMatrix, tbMatrix, &maxPosition);
		}
	}
	backtrack(tbMatrix, scoreMatrix, maxPosition, &finalScore, queryResultReverse, subjectResultReverse);

	//stop clock
	clock_t toc = clock();
	double timeElapsed = (double)(toc-tic)/CLOCKS_PER_SEC;

	#ifdef DEBUG
	printResults(finalScore, timeElapsed, queryResultReverse, subjectResultReverse);
	printMatrix(scoreMatrix);
	printTracebackMatrix(tbMatrix);
	#endif

    return 0;
}

void similarityScore(int i, int j, int* scoreMatrix, int* tbMatrix, int* maxPos) {
    int up, left, diag;

    int index = querySize * i + j;

    //Get element above
    up = scoreMatrix[index-querySize] + gapScore;

    //Get element on the left
    left = scoreMatrix[index-1] + gapScore;

    //Get element on the diagonal
    diag = scoreMatrix[index-querySize-1] + matchMismatchScore(i, j);

    //Calculates the maximum
    int max = NONE;
    int pred = NONE;

    if (diag > max) {
        max = diag;
        pred = DIAGONAL;
    }

    if (up > max) {
        max = up;
        pred = UP;
    }

    if (left > max) {
        max = left;
        pred = LEFT;
    }
    //Inserts the value in the similarity and traceback matrixes
    scoreMatrix[index] = max;
    tbMatrix[index] = pred;

    //Updates location of max value
    if (max > scoreMatrix[*maxPos]) {
        *maxPos = index;
    }

}

int matchMismatchScore(int i, int j) {
	if (subject[i-1] == query[j-1])
        return matchScore;
    else
        return mismatchScore;
}

void backtrack(int* tbMatrix, int* scoreMatrix, int maxPos, long int* finalScore, char* queryResultReverse, char* subjectResultReverse) {
    int predPos;
    int resultSize = 0;
    //backtrack from maxPos until reaches 0
    do {

        if (tbMatrix[maxPos] == DIAGONAL) {
            predPos = maxPos - querySize - 1;
        }
        else if (tbMatrix[maxPos] == UP) {
            predPos = maxPos - querySize;
        }
        else if (tbMatrix[maxPos] == LEFT) {
            predPos = maxPos - 1;
        }
        tbMatrix[maxPos] *= PATH;

        //collate results
//        *finalScore += scoreMatrix[maxPos];
//        queryResultReverse[resultSize] = query[((maxPos)%querySize)-1];
//        subjectResultReverse[resultSize++] = subject[(maxPos-1)/querySize];

        maxPos = predPos;

    } while (tbMatrix[maxPos] != NONE);
    //null terminating the strings
//    queryResultReverse[resultSize] = '\0';
//    subjectResultReverse[resultSize] = '\0';
}

void printResults(long int finalScore, double time, char* qrr, char* srr) {
	//reverse both strings
	strrev(qrr);
	strrev(srr);
	printf("\n======================================\n");
	printf("PROGRAM FINISHED\n");
//	printf("1) FINAL SCORE: %ld\n", finalScore);
//	printf("2) ALIGNMENT STRING SIZE: %d\n", strlen(qrr));
//	printf("3) ALIGNMENT STRING:\n");
//	printf("\t%s\n", qrr);
//	printf("\t");
//	for (int i=0; i<strlen(qrr); i++) {
//		if (qrr[i] == srr[i]) {
//			printf("*");
//		}
//		else {
//			printf("|");
//		}
//	}
//	printf("\n\t%s\n", srr);
	printf("4) TIME ELAPSED: %f\n", time);
  	printf("======================================\n");
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

void printTracebackMatrix(int* tbMatrix) {
    int i, j, index;
	printf("\nTraceback Matrix:\n");
    for (i = 0; i < querySize; i++) { //Lines
        for (j = 0; j < subjectSize; j++) {
        	index = subjectSize * i + j;
            if (tbMatrix[index] < 0) {
                if (tbMatrix[index] == -UP)
                    printf("N  ");
                else if (tbMatrix[index] == -LEFT)
                    printf("W  ");
                else if (tbMatrix[index] == -DIAGONAL)
                    printf("NW ");
                else
                    printf("-  ");
            }
            else {
                printf("-  ");
            }
        }
        printf("\n");
    }
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

