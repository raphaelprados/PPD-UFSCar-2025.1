/*
    Student: Raphael Alexsander Prado dos Santos
    Class: PPD - 2025.1 - UFSCar - Msc. Computer Science
    Professor: PhD. Hermes Senger
    
    Application description: A Parallel Laplace Solver Using Simple Jacobi Iteration
    Parallel Library: PThreads (C++)

    Some of the functions used here are recycled from HPC-101/examples/laplace/, as it was
    the reference for this implementation. 
    https://github.com/HPCSys-Lab/HPC-101/tree/main/examples/laplace
*/

#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<sys/time.h>

#include<pthread.h>

#define ITER_MAX 3000
#define CONV_THRESHOLD 1.0e-5f

// Custom datatypes

typedef struct chunk{
    int id;
    int size;
    // Identifies the 1D position of the first stencil
    long first_pos;
    // 'u' -> unselected, 'p' -> pending, 'c' -> computed
    char status;
    // Stablishes the maximum error for this chunk
    double local_err;
} chunk_info;

typedef struct thread_data{
    chunk_info *chunks;
    int n;
}thread_data;

// ------------------------------ Global variables ------------------------------

int size;
int iterations;
int global_chunk_size;

double *grid;
double *new_grid;
double total_error = 1.0e-1f;

pthread_mutex_t chunk_mutex;
pthread_cond_t chunks_available,
               chunks_unavailable;

bool reset_chunks = true;

// ------------------------------ Utilitary functions ------------------------------

double max(double a, double b) { return a > b ? a : b; }
double absolute(double a) { return a < 0 ? a * -1.0f : a ; }
int ijTo1D(int i, int j) { return size*i + j; }

// ------------------------------ Operation functions ------------------------------

void init_grid() {
    // seed for random generator
    srand(10);

    int linf = size / 2;
    int lsup = linf + size / 10;
    int pos_1D;
    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
            
            pos_1D = size*i + j;
            
            // inicializa região de calor no centro do grid
            if ( i>=linf && i < lsup && j>=linf && j<lsup)
                grid[pos_1D] = 100;
            else
                grid[pos_1D] = 0;
            new_grid[pos_1D] = 0.0;
        }
    }
}

// allocate memory for the grid
void allocate_memory(){
    grid = (double *) malloc(size * size * sizeof(double));
    new_grid = (double *) malloc(size * size * sizeof(double));
}

void deallocate_memory(){
    free(grid);
    free(new_grid);
}

// save the grid in a file
void save_grid(){

    char file_name[30];
    sprintf(file_name, "grid_laplace.txt");

    // save the result
    FILE *file;
    file = fopen(file_name, "w");

    for(int i = 0; i < size; i++){
        for(int j = 0; j < size; j++){
            fprintf(file, "%lf ", grid[size*i + j]);
        }
        fprintf(file, "\n");
    }

    fclose(file);
}

void print_grid(int size) {
    int i, j;

    for(i = 0; i < size; i++) {
        for(j = 0; j < size; j++) {
            printf("%lf", grid[size*i + j]);
        }
        printf("\n");
    }
}

int set_first_pos(chunk_info &ck) {

    int stencil_pos_int, stencil_pos, stencil_row, stencil_col, first_pos;

    stencil_pos_int = (ck.size*ck.id);
    stencil_pos = size + 1 + stencil_pos_int/(size-2)*2 + stencil_pos_int;

    ck.first_pos = stencil_pos;

    return stencil_pos;
}

int get_n_pos(chunk_info ck, int n) {
    int stencil_pos_int, stencil_pos, stencil_row, stencil_col, first_pos;

    stencil_pos_int = (ck.size*ck.id);
    stencil_pos = size + 1 + stencil_pos_int/(size-2)*2 + stencil_pos_int;

    return stencil_pos;
}

void chunk_traverse(chunk_info ck) {
    
    int stencil_pos_int, stencil_pos, stencil_row, stencil_col, first_pos;
    int temp;

    for(int i = 0; i < ck.size; i++) {
        stencil_pos_int = (global_chunk_size*ck.id + i);
        stencil_pos = size + 1 + stencil_pos_int/(size-2)*2 + stencil_pos_int;

        if(i == 0)
            first_pos = stencil_pos;

        stencil_row = stencil_pos / size;
        stencil_col = stencil_pos % size;
        
        if(i == 0 || i == ck.size-1)
            printf("%s: %d (%d, %d) %c ", (i == 0 ? "start" : "end"), stencil_pos, stencil_row, 
                    stencil_col, (i == 0 ? '|' : '\n'));
    }
}

void printChunk(chunk_info ck) {
    
    printf("Chunk %d : {f_pos: %d, size: %d, state: %c}\n", ck.id, ck.first_pos, ck.size, ck.status);

    // int stencil_pos_int, stencil_pos, stencil_row, stencil_col;

    // for(int i = 0; i < ck.first_pos + ck.size; i++) {
        
    //     stencil_pos_int = (ck.first_pos + i);
    //     stencil_pos = size + 1 + stencil_pos_int/(size-2)*2 + stencil_pos_int;

    //     stencil_row = stencil_pos / size;
    //     stencil_col = stencil_pos % size;
        
    //     if(i == 0 || i == ck.size-1)
    //         printf("%s: %d (%d, %d) %c ", (i == 0 ? "start" : "end"), stencil_pos, stencil_row, 
    //                 stencil_col, (i == 0 ? '|' : '\n'));
    // }
}

chunk_info* create_chunks(int n, int size) {
    chunk_info *chunks = (chunk_info*) malloc(n*sizeof(chunk_info));
    // Removing the halo region's points from the chunk creation
    int total_points = size*(size-2)-2*(size-2);
    // Setting the chunk size 
    int chunk_size = total_points / n;
    global_chunk_size = chunk_size;

    // Input error checking
    if(n <= 0 || size <= 0 || total_points/n == 0)
        return nullptr;

    for(int i = 0; i < n; i++) {
        chunks[i] = {i, chunk_size, 0, 'u'};
        set_first_pos(chunks[i]);
    }
    chunks[n-1].size += (total_points % chunk_size == 0 ? 0 : total_points % chunk_size);

    return chunks;
}

void* producer(void *arg) {
    
    thread_data *ptd = (thread_data*)arg;
    int computed_chunks;

    while(CONV_THRESHOLD < total_error) {

        computed_chunks = 0;    
        iterations++;

        printf("Producer Started\n");

        pthread_mutex_lock(&chunk_mutex);

            // Uses a thread_cond to check whether it's time to reset the chunks
            while(!reset_chunks) {
                pthread_cond_wait(&chunks_unavailable, &chunk_mutex);
            }
            
            // Resetando chunk para próxima iteração
            for(int i = 0; i < ptd->n; i++) {
                ptd->chunks[i].status = 'u';
                total_error = max(total_error, ptd->chunks[i].local_err);
                ptd->chunks[i].local_err = 1.0e-1f;
            }
            
            reset_chunks = false;
            iterations++;

            for(int i = 0; i < ptd->n; i++) {
                printChunk(ptd->chunks[i]);
                chunk_traverse(ptd->chunks[i]);
            }

            pthread_cond_signal(&chunks_available);
        pthread_mutex_unlock(&chunk_mutex);

        printf("Producer Finished\n");
    }

    return NULL;    
}

double compute_stencil(chunk_info thread_chunk) {
    double err = 0.0;           // Correctude variables

    double *result = (double*)malloc(thread_chunk.size * sizeof(double));
    int stencil_row, stencil_col, stencil_pos_int, stencil_pos; 

    for(int i = 0; i < thread_chunk.size; i++) {
        
        stencil_pos = get_n_pos(thread_chunk, i);

        stencil_row = stencil_pos / size;
        stencil_col = stencil_pos % size;
        
        result[i] = 0.25 * (grid[ijTo1D(stencil_row - 1, stencil_col)] +
                            grid[ijTo1D(stencil_row + 1, stencil_col)] + 
                            grid[ijTo1D(stencil_row, stencil_col - 1)] + 
                            grid[ijTo1D(stencil_row, stencil_col + 1)]);

        err = max(err, absolute(result[stencil_pos] - grid[stencil_pos]));
    }

    // Escreve os resultados da computação no grid
    for(int i = 0; i < thread_chunk.size; i++) {
        stencil_pos = get_n_pos(thread_chunk, i);
        new_grid[stencil_pos] = result[stencil_pos]; 
    }

    return err;
}

void* consumer(void * arg) {

    thread_data *t_data = (thread_data*) arg;
    chunk_info *t_chunk;
    int done_chunks;
    bool compute;
    double local_err;

    printf("Threads daher\n");

    while(CONV_THRESHOLD < total_error) {
        
        done_chunks = 0;
        compute = true;
        
        printf("Threads there\n");

        if(!reset_chunks) {
            printf("Threads here\n");
            pthread_mutex_lock(&chunk_mutex);
                
                for(int i = 0; i < t_data->n; i++) {
                    // Selects an available chunk
                    if(t_data->chunks[i].status == 'u') {
                        t_chunk = &t_data->chunks[i];
                        t_chunk->status = 'c';
                        printf("%dth chunk for computing\n", t_chunk->id);
                    }
                    // Checks if all chunks were already processed. If so, signals to the
                    //      producer. 
                    if(done_chunks == t_data->n) {
                        pthread_cond_signal(&chunks_unavailable);
                        reset_chunks = true;
                        compute = false;
                    }
                }
            pthread_mutex_unlock(&chunk_mutex);

            // Calcula o chunk e anota o erro
            if(compute) {
                t_chunk->local_err = compute_stencil(*t_chunk);
            }
        }   
    }
    return NULL;
}

int main(int argc, char const *argv[]) {
    
    if(argc != 3) {
        printf("Usage: ./laplace_pthreads S N T\n");
        printf("       S: The size of each side of the domain (grid)\n"
               "       N: The number of chunks that divide the domain\n"
               "       T: The number of threads\n");
        exit(-1);
    }

    // Gets the size of the grid as the first input during the program call
    size = atoi(argv[1]);

    // Sets the number of chunks 
    int n_chunks = atoi(argv[2]);

    // Gets the number of threads to be executed as the second program parameter
    int n_threads;
    n_threads = atoi(argv[3]);

    // Chunks storage variable used in the producer-consumer  
    chunk_info* chunks = create_chunks(n_chunks, size);

    // Thread variables
    pthread_t consumer_Threads[n_threads];
    pthread_t producer_Thread;

    // Thread information variable (shared between all threads)
    thread_data *thread_info = new thread_data; 
    thread_info->chunks = chunks; 
    thread_info->n = n_chunks;


    // -------------------------------- Pre-Processing Area --------------------------------


    // Allocating the grids
    allocate_memory();

    // Inserting the grid data
    init_grid();

    // Mutex initiation
    pthread_mutex_init(&chunk_mutex, NULL);
    // Condition mutex initiation
    pthread_cond_init(&chunks_available, NULL);
    pthread_cond_init(&chunks_unavailable, NULL);


    // -------------------------------- Processing Area --------------------------------


    printf("Laplace execution starting with grid[%d][%d] and %d threads!\n", size, size, n_threads);

    pthread_create(&producer_Thread, NULL, producer, (void*)thread_info);

    for(int i = 0; i < n_threads; i++) {
        pthread_create(&consumer_Threads[i], NULL, &consumer, (void*)chunks);
    }

    pthread_join(producer_Thread, NULL);

    for(int i = 0; i < n_threads; i++) {
        printf("Thread %d instantiated\n", i);
        pthread_join(consumer_Threads[i], NULL);
        
    }


    // -------------------------------- Post-Processing Area --------------------------------

    printf("After %d iterations, error is %d. Finished\n", iterations, total_error);

    // Deallocate the data 
    deallocate_memory();

    return 0;
}

/* 

Code Sources

[1] https://stackoverflow.com/questions/36890624/malloc-a-2d-array-in-c
[2] https://stackoverflow.com/questions/15937309/errorpe513-a-value-of-type-void-cannot-be-assigned-to-an-entity-of-type

*/