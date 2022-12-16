#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define NUM_NBR 4
#define BUFF 10

void rainFall(int, float, int, char*);
void readInputFile(char*, int ***, int);
void debugPrintElevations(int**, int);
void debugPrintRainDrops(float**, int);
void printAbsorptionMatrix(float**, int);
double calc_time(struct timespec, struct timespec);

int main(int argc, char** argv) {
    int rainPeriod, gridDim;
    float absorptionRate;
    char* fname;

    if (argc != 6) {
        printf("Usage: ./rainfall <P> <M> <A> <N> <elevation file>\n");
        exit(1);
    }

    // initialize variables from command line arguments
    rainPeriod = atoi(argv[2]);
    absorptionRate = atof(argv[3]);
    gridDim = atoi(argv[4]);
    fname = argv[5];

    rainFall(rainPeriod, absorptionRate, gridDim, fname);
}

void rainFall(int rainPeriod, float absorptionRate, int gridDim, char* inputFile) {
    struct timespec start_time, end_time;
    double elapsed_time;
    int row, col, nbr;
    int **elevations;
    float **rainDrops;
    float **trickleDrops;
    float **cumulativeAbsorption;
    int nbrLocs[8]; // row, col pairs for each N/S/E/W nbr stored in array of structures format
    int minElevation;
    int elevationNbrs[4];
    int minNbrs[4];
    int numWaysTie;
    float trickleDrop;
    float trickleFrac;
    float absorbAmt;
    int timestep = 0;
    int stillDraining = 1;
    
    elevations = (int**)malloc(gridDim * sizeof(int *));
    rainDrops = (float**)malloc(gridDim * sizeof(float *));
    cumulativeAbsorption = (float**)malloc(gridDim * sizeof(float *));

    for (row = 0; row <gridDim; row++) {
        elevations[row] = (int*)calloc(gridDim, sizeof(int));
        rainDrops[row] = (float*)calloc(gridDim, sizeof(float));
        cumulativeAbsorption[row] = (float*)calloc(gridDim, sizeof(float));
    }

    // initialize elevations from input file
    readInputFile(inputFile, &elevations, gridDim);
    // DEBUG
    //debugPrintElevations(elevations, gridDim);

    clock_gettime(CLOCK_MONOTONIC, &start_time);
    while (stillDraining) {
        trickleDrops = (float**)malloc(gridDim * sizeof(float *));
        for (row = 0; row < gridDim; row++) {
            trickleDrops[row] = (float*)calloc(gridDim, sizeof(float));
        }
        for (row = 0; row < gridDim; row++) {
            for (col = 0; col < gridDim; col++) {
                // receive a new rain drop if it is still raining
                if (timestep < rainPeriod) {
                    rainDrops[row][col] += 1.0;
                }
                // absorb water into the ground according to the absorption rate if there are any rain drops
                if (rainDrops[row][col] > 0.0) {
                    absorbAmt = absorptionRate > rainDrops[row][col] ? rainDrops[row][col] : absorptionRate;
                    rainDrops[row][col] -= absorbAmt;
                    cumulativeAbsorption[row][col] += absorbAmt;
                }

                // compute neighbors
                // N
                nbrLocs[0] = row == 0 ? -1 : row - 1;
                nbrLocs[1] = col;
                // S
                nbrLocs[2] = row == gridDim - 1 ? -1 : row + 1;
                nbrLocs[3] = col;
                // E
                nbrLocs[4] = row;
                nbrLocs[5] = col == gridDim - 1 ? -1 : col + 1;
                // W
                nbrLocs[6] = row;
                nbrLocs[7] = col == 0 ? -1 : col - 1;

                if (rainDrops[row][col] > 0.0) { // check to see if any remaining rain drops at this point
                    // loop over neighbors and store elevations
                    for (nbr = 0; nbr < NUM_NBR; nbr++) {
                        elevationNbrs[nbr] = -1;
                        if (nbrLocs[nbr*2] == -1 || nbrLocs[nbr*2 + 1] == -1) {
                            continue;
                        }
                        else {
                            elevationNbrs[nbr] = elevations[nbrLocs[nbr*2]][nbrLocs[nbr*2+1]];
                        }
                    }

                    // determine minimum elevation
                    minElevation = 99999999;
                    for (nbr = 0; nbr < NUM_NBR; nbr++) {
                        if (elevationNbrs[nbr] >= 0) {
                            if (minElevation > elevationNbrs[nbr]) {
                                minElevation = elevationNbrs[nbr];
                            }
                        }
                    }

                    // make sure we don't trickle if we are at the lowest elevation
                    if (elevations[row][col] > minElevation) {
                        // determine if there is a tie
                        numWaysTie = 0;
                        for (nbr = 0; nbr < NUM_NBR; nbr++) {
                            minNbrs[nbr] = 0;
                            if (elevationNbrs[nbr] < 0) {
                                continue;
                            }
                            else if (elevationNbrs[nbr] == minElevation) { // should never have a situation where numWaysTie = 0
                                minNbrs[nbr] = 1;
                                numWaysTie++;
                            }
                        }
                        // stream trickle updates to another array
                        trickleFrac = (rainDrops[row][col] < 1.0 ? rainDrops[row][col] : 1.0) / numWaysTie;
                        for (nbr = 0; nbr < NUM_NBR; nbr++) {
                            if (minNbrs[nbr] == 1) {
                                trickleDrops[nbrLocs[nbr*2]][nbrLocs[nbr*2+1]] += trickleFrac;
                                rainDrops[row][col] -= trickleFrac;
                            }
                        }

                    }

                }

            }
        }
        // make second traversal over all landscape points to update
        for (row = 0; row < gridDim; row++) {
            for (col = 0; col < gridDim; col++) {
                rainDrops[row][col] += trickleDrops[row][col];
            }
        }

        // check to see if still draining
        stillDraining = 0;
        for (row = 0; row < gridDim; row++) {
            for (col=0; col < gridDim; col++) {
                if (rainDrops[row][col] > 0.0) {
                    stillDraining = 1;
                }
            }
        }
        timestep++;
        for (row = 0; row < gridDim; row++) {
            free(trickleDrops[row]);
        }
        free(trickleDrops);
    }
    clock_gettime(CLOCK_MONOTONIC, &end_time);
    elapsed_time = calc_time(start_time, end_time) / 1000000000;

    // output simulation results
    printf("Rainfall simulation took %d time steps to complete.\n", timestep);
    printf("Runtime = %f seconds\n\n", elapsed_time);
    printf("The following grid shows the number of raindrops absorbed at each point:\n");
    printAbsorptionMatrix(cumulativeAbsorption, gridDim);
    for (row = 0; row < gridDim; row++) {
        free(elevations[row]);
        free(rainDrops[row]);
        free(cumulativeAbsorption[row]);
    }
    free(elevations);
    free(rainDrops);
    free(cumulativeAbsorption);
}

void readInputFile(char* inFile, int*** elevationMat, int gridDim) {
    FILE* myFile;
    int row, col;
    char c;
    char * digits;
    int numDigits = 0;

    myFile = fopen(inFile, "r");
    if (!myFile) {
        printf("ERROR: Unable to open file!\n");
        exit(1);
    }

    for (row = 0; row < gridDim; row++) {
        for (col = 0; col < gridDim; col++) {
            digits = (char*)calloc(BUFF, sizeof(char));
            c = fgetc(myFile);
            if (c == ' ') {
                do{
                    c = fgetc(myFile);
                } while (c == ' ');
            }
            while (c!= ' ' && c!= '\n'){
                digits[numDigits] = c;
                numDigits++;
                c = fgetc(myFile);

            }
            (*elevationMat)[row][col] = atoi(digits);
            numDigits = 0;
            free(digits);
        }
    }
    fclose(myFile);

}

void debugPrintElevations(int** elevationMat, int gridDim) {
    int row, col;
    for (row = 0; row < gridDim; row++) {
        for (col = 0; col < gridDim; col++) {
            printf("%d ", elevationMat[row][col]);
        }
        printf("\n");
    }
}

void debugPrintRainDrops(float** rainDrops, int gridDim) {
    int row, col;
    for (row = 0; row < gridDim; row++) {
        for (col = 0; col < gridDim; col++) {
            printf("%g ", rainDrops[row][col]);
        }
        printf("\n");
    }
}

void printAbsorptionMatrix(float** cumulativeAbs, int gridDim) {
    int row, col;
    for (row = 0; row < gridDim; row++) {
        for (col = 0; col < gridDim; col++) {
            printf("%8g", cumulativeAbs[row][col]);
        }
        printf("\n");
    }
}

double calc_time(struct timespec start, struct timespec end) {
  double start_sec = (double)start.tv_sec*1000000000.0 + (double)start.tv_nsec;
  double end_sec = (double)end.tv_sec*1000000000.0 + (double)end.tv_nsec;
  if (end_sec < start_sec) {
    return 0;
  } else {
    return end_sec - start_sec;
  }
}