/***********************************************************************
 * Smith–Waterman algorithm
 * Purpose:     Local alignment of nucleotide or protein sequences
 * Authors:     Daniel Holanda, Hanoch Griner, Taynara Pinheiro
 ***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#define PATH -1
#define NONE 0
#define UP 1
#define LEFT 2
#define DIAGONAL 3
//Define scores
#define matchScore 5
#define mismatchScore -3
#define gapScore -4

void similarityScore(long long int i, long long int j, int** scoreMatrix, int** tbMatrix, long long int* maxPosRow, long long int* maxPosCol);
int matchMismatchScore(long long int i, long long int j);
void backtrack(int** tbMatrix, long long int maxPosRow, long long int maxPosCol);
void printMatrix(int** matrix);
void printPredecessorMatrix(int** matrix);
void readFiles(char* queryFile, char* subjectFile);

 //Defines size of strings to be compared
int querySize = 0; //Columns - Size of string a
int subjectSize = 0;  //Lines - Size of string b


//Strings over the Alphabet Sigma
char* query, * subject;

/*--------------------------------------------------------------------
 * Function:    main
 */
int main(int argc, char* argv[]) {
	if (argc != 3) {
		printf("Please enter in this format: SmithW <query_file_name> <subject_file_name>\n");
		return 1;
	}
	char* queryFile = argv[1];
	char* subjectFile = argv[2];
	readFiles(queryFile, subjectFile);

	//zeroes in the first row and column
	querySize++;
	subjectSize++;

	//create 2d score matrix and 2d traceback matrix
	int **scoreMatrix = (int**) calloc(querySize, sizeof(int));
	int **tbMatrix = (int **) calloc(querySize, sizeof(int));

	for (int row = 0; row < querySize; row++) {
		scoreMatrix[row] = (int *) calloc(subjectSize, sizeof(int));
		tbMatrix[row] = (int *) calloc(subjectSize, sizeof(int));
	}

	long long int maxPositionRow =0;
	long long int maxPositionCol =0;

	for (int i=1; i<subjectSize; i++) {
		for (int j=1; j<querySize; j++) {
			similarityScore(i, j, scoreMatrix, tbMatrix, &maxPositionRow, &maxPositionCol);
		}
	}
	backtrack(tbMatrix, maxPositionRow, maxPositionCol);

	#ifdef DEBUG
	printf("\nSimilarity Matrix:\n");
	printMatrix(scoreMatrix);

	printf("\nPredecessor Matrix:\n");
	printPredecessorMatrix(tbMatrix);
	#endif

	//free strings
	free(query);
	free(subject);

	//free matrix
	free(scoreMatrix);
	free(tbMatrix);

    return 0;
}



/*--------------------------------------------------------------------
 * Function:    SimilarityScore
 * Purpose:     Calculate  the maximum Similarity-Score H(i,j)
 */
void similarityScore(long long int i, long long int j, int** scoreMatrix, int** tbMatrix, long long int* maxPosRow, long long int* maxPosCol) {

    int up, left, diag;


    //Get element above
    up = scoreMatrix[i-1][j] + gapScore;

    //Get element on the left
    left = scoreMatrix[i][j-1] + gapScore;

    //Get element on the diagonal
    diag = scoreMatrix[i-1][j-1] + matchMismatchScore(i, j);

    //Calculates the maximum
    int max = NONE;
    int pred = NONE;
    /* === Matrix ===
     *      a[0] ... a[n]
     * b[0]
     * ...
     * b[n]
     *
     * generate 'a' from 'b', if 'W' insert e 'N' remove
     * a=GAATTCA
     * b=GACTT-A
     *
     * generate 'b' from 'a', if 'W' insert e 'N' remove
     * b=GACTT-A
     * a=GAATTCA
    */

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
    //Inserts the value in the similarity and predecessor matrixes
    scoreMatrix[i][j] = max;
    tbMatrix[i][j] = pred;

    //Updates maximum score to be used as seed on backtrack
    if (max > scoreMatrix[*maxPosRow][*maxPosCol]) {
        *maxPosRow = i;
        *maxPosCol = j;
    }

}

int matchMismatchScore(long long int i, long long int j) {
    if (subject[i] == query[j])
        return matchScore;
    else
        return mismatchScore;
}

/*--------------------------------------------------------------------
 * Function:    backtrack
 * Purpose:     Modify matrix to print, path change from value to PATH
 */
void backtrack(int** tbMatrix, long long int maxPosRow, long long int maxPosCol) {
    //hold maxPos value
    long long int predPosRow;
    long long int predPosCol;

    //backtrack from maxPos to startPos = 0
    do {
        if (tbMatrix[maxPosRow][maxPosCol] == DIAGONAL) {
            predPosRow = maxPosRow - 1;
            predPosCol = maxPosCol - 1;
        }
        else if (tbMatrix[maxPosRow][maxPosCol] == UP) {
            predPosRow = maxPosRow - 1;
        }
        else if (tbMatrix[maxPosRow][maxPosCol] == LEFT) {
            predPosCol = maxPosCol - 1;
        }
        tbMatrix[maxPosRow][maxPosCol] *= PATH;
        maxPosRow = predPosRow;
        maxPosCol = predPosCol;
    } while (tbMatrix[maxPosRow][maxPosCol] != NONE);
}

/*--------------------------------------------------------------------
 * Function:    printMatrix
 * Purpose:     Print Matrix
 */
void printMatrix(int** matrix) {
    long long int i, j;
    for (i = 0; i < subjectSize; i++) { //Lines
        for (j = 0; j < querySize; j++) {
            printf("%d\t", matrix[i][j]);
        }
        printf("\n");
    }

}

void printPredecessorMatrix(int** matrix) {
    long long int i, j;
    for (i = 0; i < subjectSize; i++) { //Lines
        for (j = 0; j < querySize; j++) {
            if (matrix[i][j] < 0) {
                if (matrix[i][j] == -UP)
                    printf("N  ");
                else if (matrix[i][j] == -LEFT)
                    printf("W  ");
                else if (matrix[i][j] == -DIAGONAL)
                    printf("NW ");
                else
                    printf("-  ");
            }
            else {
                if (matrix[i][j] == UP)
                    printf("N  ");
                else if (matrix[i][j] == LEFT)
                    printf("W  ");
                else if (matrix[i][j] == DIAGONAL)
                    printf("NW ");
                else
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
		printf("%d\n", querySize);
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
		printf("%d\n", subjectSize);
		fseek(sfp, 0, SEEK_SET);
		subject = malloc(subjectSize);
		if (subject) {
			fread(subject, 1, subjectSize, sfp);
		}
		fclose(sfp);
	}
}

