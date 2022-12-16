#include <stdio.h>
#include <pthread.h>
#include <stdlib.h>
#include <time.h>

#define NUM_NBR 4
#define BUFF 10
#define HALO_DEPTH 1

pthread_barrier_t barrier;
int *draining;
double *elapsed_times;
float **rainDrops;
float **trickleDrops;
float **cumulativeAbsorption;
int** elevations;

typedef struct inputs {
    int rainTime;
    int gridSize;
    int haloSize;
    int tid;
    int numThread;
    float absorbRate;
} inputs;

void* rainFall(void* ptr);
void readInputFile(char*, int ***, int);
//DEBUG
//void debugPrintElevations(int**, int);
//void debugPrintRainDrops(float**, int, int, int, int);
void printAbsorptionMatrix(float**, int);
double calc_time(struct timespec, struct timespec);

int main(int argc, char** argv) {
    pthread_t *threads;

    int numThreads, rainPeriod, gridDim, haloDim;
    float absorptionRate;
    char* fname;
    int i;

    if (argc != 6) {
        printf("Usage: ./rainfall <P> <M> <A> <N> <elevation file>\n");
        exit(1);
    }

    // initialize variables from command line arguments
    numThreads = atoi(argv[1]);
    rainPeriod = atoi(argv[2]);
    absorptionRate = atof(argv[3]);
    gridDim = atoi(argv[4]);
    fname = argv[5];

    int tids[numThreads];

    haloDim = gridDim + 2 * numThreads - 2;

    threads = (pthread_t*)malloc(numThreads * sizeof(pthread_t));
    elevations = (int**)malloc(gridDim * sizeof(int *));    
    rainDrops = (float**)malloc(haloDim * sizeof(float *));
    trickleDrops = (float**)malloc(haloDim * sizeof(float *));
    cumulativeAbsorption = (float**)malloc(gridDim * sizeof(float *));
    elapsed_times = (double*)calloc(numThreads, sizeof(double));
    draining = (int*)calloc(numThreads, sizeof(int));
    for (i=0; i <numThreads; i++) {
        draining[i] = 1;
    }

    pthread_barrier_init(&barrier, NULL, numThreads);

    for (i = 0; i < haloDim; i++) {
        rainDrops[i] = (float*)calloc(gridDim, sizeof(float));
        trickleDrops[i] = (float*)calloc(gridDim, sizeof(float));
    }

    for (i = 0; i <gridDim; i++) {
        elevations[i] = (int*)calloc(gridDim, sizeof(int));
        cumulativeAbsorption[i] = (float*)calloc(gridDim, sizeof(float));
    }

    readInputFile(fname, &elevations, gridDim);
    inputs* threadArgs[numThreads];
    for (i=0; i<numThreads; i++) {
        threadArgs[i] = malloc(sizeof(inputs));
        threadArgs[i]->rainTime = rainPeriod;
        threadArgs[i]->absorbRate = absorptionRate;
        threadArgs[i]->gridSize = gridDim;
        threadArgs[i]->haloSize = haloDim;
        threadArgs[i]->numThread = numThreads;
        threadArgs[i]->tid = i;
        pthread_create(&threads[i], NULL, rainFall, (void*)threadArgs[i]);
    }
    for (i=0; i<numThreads; i++) {
        pthread_join(threads[i], NULL);
    }
    pthread_barrier_destroy(&barrier);
    for (i=0; i<gridDim; i++) {
        free(elevations[i]);
        free(cumulativeAbsorption[i]);
    }
    for (i=0; i<haloDim; i++) {
        free(rainDrops[i]);
        free(trickleDrops[i]);
    }
    for (i=0; i<numThreads; i++) {
        free(threadArgs[i]);
    }
    free(threads);
    free(elevations);
    free(cumulativeAbsorption);
    free(rainDrops);
    free(trickleDrops);
}

void* rainFall(void* ptr) {
    inputs* myArgs = (inputs*)ptr;

    int rainPeriod = myArgs->rainTime;
    float absorptionRate = myArgs->absorbRate;
    int gridDim = myArgs->gridSize;
    int haloDim = myArgs->haloSize;
    int tidx = myArgs->tid;
    int numThreads = myArgs->numThread;
    
    struct timespec start_time, end_time;
    double elapsed_time;
    int row, col, nbr, td;
    int startrow, stoprow;
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
    int offset = 2*HALO_DEPTH*tidx;

    double chunk_size = haloDim/(double)numThreads;
    startrow = tidx == 0 ? chunk_size * tidx : chunk_size * tidx + HALO_DEPTH;
    stoprow = tidx == numThreads - 1 ? chunk_size * (tidx + 1) : chunk_size * (tidx + 1) - HALO_DEPTH;

    if (tidx == numThreads - 1) {
        if (stoprow != haloDim) stoprow = haloDim;
    }

    clock_gettime(CLOCK_MONOTONIC, &start_time);
    while (stillDraining){
        int end_idx = (tidx == numThreads - 1 ? stoprow : stoprow + 2 * HALO_DEPTH);
        // RESET THE TRICKLE DROPS
        for (row = startrow; row < end_idx; row++) {
            for (col = 0; col < gridDim; col++) {
                trickleDrops[row][col] = 0.0;
            }
        }
        pthread_barrier_wait(&barrier);
        for (row = startrow; row < stoprow; row++) {
            for (col = 0; col < gridDim; col++) {
                // receive a new rain drop if it is still raining
                if (timestep < rainPeriod) {
                    rainDrops[row][col] += 1.0;
                }
                // absorb water into the ground according to the absorption rate if there are any rain drops
                if (rainDrops[row][col] > 0.0) {
                    absorbAmt = absorptionRate > rainDrops[row][col] ? rainDrops[row][col] : absorptionRate;
		            rainDrops[row][col] -= absorbAmt;
                    cumulativeAbsorption[row - offset][col] += absorbAmt;
                }

                // compute neighbors
                // N
                nbrLocs[0] = row == 0 ? -1 : row - 1;
                nbrLocs[1] = col;
                // S
                nbrLocs[2] = row == haloDim - 1 ? -1 : row + 1;
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
                            elevationNbrs[nbr] = elevations[nbrLocs[nbr*2] - offset][nbrLocs[nbr*2+1]];
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
                    if (elevations[row - offset][col] > minElevation) {
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
        pthread_barrier_wait(&barrier);
        // HALO EXCHANGE OF TRICKLE DROPS ACROSS NEIGHBORING THREADS
        if (tidx != numThreads - 1) {
            for (col = 0; col < gridDim; col++) {
                trickleDrops[stoprow-HALO_DEPTH][col] += trickleDrops[stoprow + HALO_DEPTH][col];
            }
        }
        if (tidx != 0) {
            for (col = 0; col < gridDim; col++) {
                trickleDrops[startrow][col] += trickleDrops[startrow - 2*HALO_DEPTH][col];
            }
        }
        pthread_barrier_wait(&barrier);

        // make second traversal over all landscape points to update
        for (row = startrow; row < stoprow; row++) {
            for (col = 0; col < gridDim; col++) {
                rainDrops[row][col] += trickleDrops[row][col];
            }
        }

        // check to see if still draining
        stillDraining = 0;
        draining[tidx] = 0;
        for (row = startrow; row < stoprow; row++) {
            for (col=0; col < gridDim; col++) {
                if (rainDrops[row][col] > 0.0) {
                    draining[tidx] = 1;
                }
            }
        }
        pthread_barrier_wait(&barrier);
        for (td = 0; td < numThreads; td++) {
            if (draining[td] == 1) {
                stillDraining = 1;
                break;
            }
        }
        timestep++;
    }
    pthread_barrier_wait(&barrier);
    clock_gettime(CLOCK_MONOTONIC, &end_time);
    elapsed_time = calc_time(start_time, end_time) / 1000000000;
    elapsed_times[tidx] = elapsed_time;
    pthread_barrier_wait(&barrier);
    if (tidx == 0) {
        // get minimum elapsed time
        double minTime = 99999999;
        for (td = 0; td < numThreads; td++) {
            if (minTime > elapsed_times[td]) {
                minTime = elapsed_times[td];
            }
        }
        // output simulation results
        printf("Rainfall simulation took %d time steps to complete.\n", timestep);
        printf("Runtime = %f seconds\n\n", minTime);
        printf("The following grid shows the number of raindrops absorbed at each point:\n");
        printAbsorptionMatrix(cumulativeAbsorption, gridDim);
    }
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

//DEBUG
/*void debugPrintElevations(int** elevationMat, int gridDim) {
    int row, col;
    for (row = 0; row < gridDim; row++) {
        for (col = 0; col < gridDim; col++) {
            printf("%d ", elevationMat[row][col]);
        }
        printf("\n");
    }
}

void debugPrintRainDrops(float** rainDrops, int start, int stop, int gridDim, int numCols) {
    int row, col;
    for (row = 0; row < gridDim; row++) {
        for (col = 0; col < numCols; col++) {
            printf("%g ", rainDrops[row][col]);
        }
        printf("\n");
    }
}*/

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
