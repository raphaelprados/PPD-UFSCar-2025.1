/*
    This program solves Laplace's equation on a regular 2D grid using simple Jacobi iteration.

    The stencil calculation stops when  ( err >= CONV_THRESHOLD OR  iter > ITER_MAX )
*/

// Must be used to allow barrier to work. Source: https://stackoverflow.com/questions/61647896/unknown-type-name-pthread-barrier-t
#define _XOPEN_SOURCE 600

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <pthread.h>
#include <stdbool.h>
#include <sys/time.h>

#define ITER_MAX 3000 // number of maximum iterations
#define CONV_THRESHOLD 1.0e-5f // threshold of convergence

// --------------------------------- Custom data types ---------------------------------

typedef struct thread_data{
    int begin, end;
    int id;
}t_data;

// --------------------------------- Global variables ---------------------------------

// matrix to be solved
double *grid;

// auxiliary matrix
double *new_grid;

// thread error buffer
double *error_buffer;

// size of each side of the grid
int size;

// number of threads
int n_threads;

// number of current iterations performed
int iter = 0;

// control variable for threads to select a row
int next_row = 1;

// the global error for all rows
double err;

// Auxiliar variable for printing debug information 
bool debug;

// Barrier variable that ensure synchronyzation
pthread_barrier_t barrier;

// Mutex used for row selection
pthread_mutex_t row_selec;

// return the maximum value
double max(double a, double b){
    if(a > b)
        return a;
    return b;
}

// return the absolute value of a number
double absolute(double num){
    if(num < 0)
        return -1.0 * num;
    return num;
}

// allocate memory for the grid
void allocate_memory(){
    grid = (double *) malloc(size * size * sizeof(double));
    new_grid = (double *) malloc(size * size * sizeof(double));
}

// initialize the grid
void initialize_grid(){
    // seed for random generator
    srand(10);

    int linf = size / 2;
    int lsup = linf + size / 10;
    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
            // inicializa região de calor no centro do grid
            if ( i>=linf && i < lsup && j>=linf && j<lsup)
                grid[i*size+j] = 100;
            else
               grid[i*size+j] = 0;
            new_grid[i*size+j] = 0.0;
        }
    }
}

// Initialize the id, begin and end of a thread computing
void initialize_thread_data(t_data *thread_info) {
    
    int base_chunk_size = (size - 2) / n_threads;
    int i;

    for(i = 0; i < n_threads; i++) {
        thread_info[i].begin = i*base_chunk_size+1;
        thread_info[i].end = thread_info[i].begin + base_chunk_size + ((i == n_threads - 1) ? ((size - 2) % n_threads) : 0);
        thread_info[i].id = i;
        // printf("Thread %d: {%d, %d}\n", thread_info[i].id, thread_info[i].begin, thread_info[i].end);
    }
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
            fprintf(file, "%lf ", grid[i*size+j]);
        }
        fprintf(file, "\n");
    }

    fclose(file);
}

void print_grid() {
    int i, j;

    for(i = 0; i < size; i++) {
        for(j = 0; j < size; j++) {
            printf("| %.1lf |", grid[i*size+j]);
        }
        printf("\n");
    }
}

void free_grid() {
    free(grid);
    free(new_grid);
}

double stencil(int i) {
    double local_error = 0.0;

    // Jacobi iteration
    // This loop will end if either the maximum change reaches below a set threshold (convergence)
    // or a fixed number of maximum iterations have completed
    for(int j = 1; j < size-1; j++) {

        new_grid[i*size+j] = 0.25 * (grid[i*size + j+1] + grid[i*size + j-1] +
                                    grid[(i-1)*size+j] + grid[(i+1)*size+j]);

        local_error = max(local_error, absolute(new_grid[i*size+j] - grid[i*size+j]));
    }

    return local_error;
}

void update_grid() {
    
    int i, j;
    
    // Swaping the grids -> O(n²)
    for( int i = 1; i < size-1; i++) {
        for( int j = 1; j < size-1; j++) {
            grid[i*size+j] = new_grid[i*size+j];
        }
    }

    // Increasing the iteration count
    iter++;

    // next
    next_row = 1;

    // Resets the global error for the current iteration calculation
    err = 0.0;

    // Calculates the global error
    for(i = 0; i < n_threads; i++) {
        err = max(err, error_buffer[i]);
    }

    // if(iter % 100 == 0)
        // printf("Iteration %d with %.10lf error\n", iter, err);
    
    if(debug) {
        print_grid();
        getchar();
    }
}

void* thread_manager(void *arg) {
    
    t_data *thread_info = (t_data*)arg;
    double local_error = 0.0;
    int i;

    do {
        // Resets the error for the current iteration
        local_error = 0.0;

        // Takes continuous rows for processing. Predefined row assignment
        i = 1;
        while(true) {
            pthread_mutex_lock(&row_selec);
                i = next_row;
                next_row++;
                // printf("Thread %d picked the %dth row\n", thread_info->id, i);
            pthread_mutex_unlock(&row_selec);

            if(i < size - 1)
                local_error = max(local_error, stencil(i));
            else
                break;

            if(debug)
                printf("[Iter: %d, thread: %d]: error: %.6lf\n", iter, thread_info->id, local_error);
        }

        // Prevents false sharing by only updating the buffer once
        error_buffer[thread_info->id] = local_error;

        if(debug)
            printf("Thread %d, err: %.6lf\n", thread_info->id, local_error);

        // First barrier (synchronize for update)
        pthread_barrier_wait(&barrier);
        
        // Thread 0 updates control variables, the grid and the global error
        if(thread_info->id == 0)
            update_grid();

        // Second barrier (synchronize for next iteration)
        pthread_barrier_wait(&barrier);
    } while ( err > CONV_THRESHOLD && iter <= ITER_MAX );

}

int main(int argc, char *argv[], char *envp[]){

    if(argc != 3){
        printf("Usage: ./laplace_seq N\n");
        printf("N: The size of each side of the domain (grid)\n"
               "T: The number of threads used\n");
        exit(-1);
    }

    size = atoi(argv[1]);

    n_threads = atoi(argv[2]);

    debug = false;

    // --------------------------------- Variable declaration ---------------------------------

    // variables to measure execution time
    struct timeval time_start;
    struct timeval time_end;

    int i;

    pthread_t *threads = malloc(n_threads * sizeof(pthread_t));

    // --------------------------------- Pre-Processing ---------------------------------

    // allocate memory to the grid (matrix)
    allocate_memory();

    // allocates the thread error buffer
    error_buffer = malloc(n_threads * sizeof(double)); 

    // set grid initial conditions
    initialize_grid();

    // Thread data initialization
    t_data *thread_info = malloc(n_threads*sizeof(t_data)); 
    initialize_thread_data(thread_info);

    // Mutex initialization
    pthread_mutex_init(&row_selec, NULL);

    // --------------------------------- Processing (Threads) ---------------------------------

    pthread_barrier_init(&barrier, NULL, n_threads);

    printf("Jacobi relaxation calculation: %d x %d grid\n", size, size);

    // get the start time
    gettimeofday(&time_start, NULL);

    // Thread initialization
    for(i = 0; i < n_threads; i++) {
        // printf("Thread %d created\n", thread_info[i].id);
        pthread_create(&threads[i], NULL, &thread_manager, (void*)&(thread_info[i]));
    }

    // Thread gathering
    for(i = 1; i <= n_threads; i++) 
        pthread_join(threads[i], NULL);

    // get the end time
    gettimeofday(&time_end, NULL);

    // --------------------------------- Post-Processing ---------------------------------

    double exec_time = (double) (time_end.tv_sec - time_start.tv_sec) +
                       (double) (time_end.tv_usec - time_start.tv_usec) / 1000000.0;

    //save the final grid in file
    save_grid();

    free_grid();
    free(error_buffer);
    free(threads);

    printf("\nKernel executed in %lf seconds with %d iterations and error of %0.10lf\n", exec_time, iter, err);

    return 0;
}
