#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <windows.h>

#include "H1.h"
#include "b1.h"

// #include "H80p13.h"
// #include "H140p7.h"
// #include "H200p5.h"
// #include "H260p2.h"
// #include "H320p2.h"
// #include "H380p2.h"
// #include "b80.h"
// #include "b140.h"
// #include "b200.h"
// #include "b260.h"
// #include "b320.h"
// #include "b380.h"

#include "H440p2.h"
#include "H440p5.h"
#include "H440p15.h"
#include "H440p50.h"
#include "H440p95.h"
#include "b440.h"

// #include "H512p2.h"
// #include "b512.h"

// #include "H1024p2.h"
// #include "H1024p15.h"
// #include "H1024p50.h"
// #include "b1024.h"

#include "H2048p2.h"
//#include "H2048p15.h"
// #include "H2048p50.h"
#include "b2048.h"

//#include "H4096p2.h"
//#include "H4096p15.h"
// #include "H4096p50.h"
// #include "b4096.h"

#define H H2048p2
#define b b2048
#define N 2048
#define IDX(i,j,n) ((i)*(n)+(j))

typedef struct {
    int n;
    int nnz;
    float *values;
    int *col_idx;
    int *row_ptr;
} CSRMatrix;

CSRMatrix* dense_to_csr(float *dense, int n) {
    // count non-zero number
    int nnz = 0;
    for (int i = 0; i < n * n; i++) {
        if (fabsf(dense[i]) > 1e-8f) {
            nnz++;
        }
    }

    // allocate CSR matrix size
    CSRMatrix *csr = malloc(sizeof(CSRMatrix));
    if (!csr) return NULL;

    csr->n = n;
    csr->nnz = nnz;
    csr->values  = malloc(sizeof(float) * nnz);
    csr->col_idx = malloc(sizeof(int)   * nnz);
    csr->row_ptr = malloc(sizeof(int)   * (n + 1));

    if (!csr->values || !csr->col_idx || !csr->row_ptr) {
        free(csr->values);
        free(csr->col_idx);
        free(csr->row_ptr);
        free(csr);
        return NULL;
    }

    // fill in CSR with non-zero
    int pos = 0;
    csr->row_ptr[0] = 0;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            float val = dense[IDX(i, j, n)];
            if (fabsf(val) > 1e-8f) {
                csr->values[pos]  = val;
                csr->col_idx[pos] = j;
                pos++;
            }
        }
        csr->row_ptr[i + 1] = pos;
    }

    if (pos != nnz) {
        fprintf(stderr, "dense_to_csr: nnz mismatch (%d != %d)\n", pos, nnz);
    }

    return csr;
}
void free_csr(CSRMatrix *csr){
    free(csr->values);
    free(csr->col_idx);
    free(csr->row_ptr);
    free(csr);
}

// ------- RCM -------
typedef struct {
    int *queue;
    int front, rear;
} Queue;

void queue_init(Queue *q, int n){
    q->queue = malloc(sizeof(int)*n);
    q->front = 0; q->rear = 0;
}
void queue_push(Queue *q, int v){
    q->queue[q->rear++] = v;
}
int queue_pop(Queue *q){
    return q->queue[q->front++];
}
int queue_empty(Queue *q){
    return q->front == q->rear;
}
void queue_free(Queue *q){
    free(q->queue);
}

void compute_degrees(int *degrees, int *adj, int n){
    for(int i=0; i<n; i++){
        int deg=0;
        for(int j=0; j<n; j++){
            if(adj[IDX(i,j,n)]) deg++;
        }
        degrees[i] = deg;
    }
}

int find_min_degree_node(int *degrees, int n){
    int min_deg = degrees[0];
    int min_node = 0;
    for(int i=1; i<n; i++){
        if(degrees[i] < min_deg){
            min_deg = degrees[i];
            min_node = i;
        }
    }
    return min_node;
}

void rcm(int *adj, int n, int *order){
    int *visited = calloc(n, sizeof(int));
    int *degrees = malloc(sizeof(int)*n);
    compute_degrees(degrees, adj, n);

    int order_index = 0;
    Queue q;
    queue_init(&q, n);

    int *neighbors = malloc(sizeof(int)*n);

    while(order_index < n) {
        // find the start node with minimum degree
        int start = find_min_degree_node(degrees, n);;
        int min_deg = INT_MAX;
        for(int i=0; i<n; i++){
            if(!visited[i] && degrees[i] < min_deg){
                min_deg = degrees[i];
                start = i;
            }
        }
        if(start == -1) break;  // all visited

        queue_push(&q, start);
        visited[start] = 1;

        while(!queue_empty(&q)){
            int u = queue_pop(&q);
            order[order_index++] = u;

            int nb_count = 0;
            for(int v=0; v<n; v++){
                if(adj[IDX(u,v,n)] && !visited[v]){
                    neighbors[nb_count++] = v;
                }
            }
            // order neighobrs with increase degree
            for(int i=0; i<nb_count-1; i++){
                for(int j=i+1; j<nb_count; j++){
                    if(degrees[neighbors[i]] > degrees[neighbors[j]]){
                        int tmp = neighbors[i];
                        neighbors[i] = neighbors[j];
                        neighbors[j] = tmp;
                    }
                }
            }
            for(int i=0; i<nb_count; i++){
                visited[neighbors[i]] = 1;
                queue_push(&q, neighbors[i]);
            }
        }
    }

    free(neighbors);
    queue_free(&q);
    free(visited);
    free(degrees);

    // printf("CM order is:\n");
    // for(int i=0; i<n; i++){
    //     printf("%d ", order[i]);
    // }
    // printf("\n");

    // inverse CM to get RCM
    for(int i=0; i<n/2; i++){
        int tmp = order[i];
        order[i] = order[n-1 - i];
        order[n-1 - i] = tmp;
    }

    // printf("RCM order is:\n");
    // for(int i=0; i<n; i++){
    //     printf("%d ", order[i]);
    // }
    // printf("\n");
}


void reorder_matrix_vector(float *H_in, float *b_in, float *H_out, float *b_out, int *order, int n){
    for(int i=0; i<n; i++){
        b_out[i] = b_in[order[i]];
        for(int j=0; j<n; j++){
            H_out[IDX(i,j,n)] = H_in[IDX(order[i], order[j], n)];
        }
    }
}
void inverse_reorder_vector(float *x_in, float *x_out, int *order, int n){
    for(int i=0; i<n; i++){
        x_out[order[i]] = x_in[i];
    }
}

void reorder_adj_matrix(int *adj_in, int *adj_out, int *order, int n){
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            adj_out[IDX(i,j,n)] = adj_in[IDX(order[i], order[j], n)];
        }
    }
}

// calculate bandwidth
int compute_bandwidth(int *order, int n, CSRMatrix *csr){
    int bandwidth = 0;
    for(int i=0; i<n; i++){
        int row_start = csr->row_ptr[i];
        int row_end = csr->row_ptr[i+1];
        for(int idx=row_start; idx<row_end; idx++){
            int col = csr->col_idx[idx];
            int dist = abs(i - col);
            if(dist > bandwidth) bandwidth = dist;
        }
    }
    return bandwidth;
}

// navice way Cholesky O^3 for dense matrix, less for sparse matrix
int cholesky_csr_plain(CSRMatrix *H_csr, float *L_dense) {
    int n = H_csr->n;
    memset(L_dense, 0, sizeof(float)*n*n);

    for(int i=0; i<n; i++){  // iteration in H matrix - row
        int row_start = H_csr->row_ptr[i];
        int row_end = H_csr->row_ptr[i+1];
        //printf("row start:%d, row end:%d\n",row_start,row_end);
        for(int idx_j = row_start; idx_j < row_end; idx_j++){  // iteration in H matrix - col
            int j = H_csr->col_idx[idx_j];
            if(j > i) continue; // compute lower L matrix only

            float sum = 0.f;
            for(int k=0; k<j; k++){ // calculate sum
                sum += L_dense[IDX(i,k,n)] * L_dense[IDX(j,k,n)];
            }

            if(i == j){
                float val = H_csr->values[idx_j] - sum;
                if(val <= 0) return -1;
                L_dense[IDX(i,j,n)] = sqrtf(val);
            } else {
                if(fabs(L_dense[IDX(j,j,n)]) < 1e-12) return -1;
                L_dense[IDX(i,j,n)] = (H_csr->values[idx_j] - sum) / L_dense[IDX(j,j,n)];
            }
        }
    }
    return 0;
}



// Cholesky with bandwidth constrain
int cholesky_csr_bandwidth(CSRMatrix *H_csr, float *L_dense, int bandwidth) {
    int n = H_csr->n;
    memset(L_dense, 0, sizeof(float) * n * n);

    for (int i = 0; i < n; i++) { // iteration in H matrix - row
        int row_start_i = H_csr->row_ptr[i];    // RCM order
        int row_end_i   = H_csr->row_ptr[i + 1];// RCM order

        // bandwidth constrain, which col to start
        int k_lower_bound = (i > bandwidth) ? (i - bandwidth) : 0;

        for (int idx_j = row_start_i; idx_j < row_end_i; idx_j++) { // iteration in H matrix - col
            int j = H_csr->col_idx[idx_j];      // RCM order
            if (j > i) continue; // lower L only

            int diff = i - j;
            if (diff > bandwidth) continue; // skip if the start location is outside the bandwidth

            float sum = 0.f;

            // sparse row from start to end 
            int row_start_j = H_csr->row_ptr[j]; 
            int row_end_j   = H_csr->row_ptr[j + 1];

            int pi = row_start_i;
            int pj = row_start_j;

            while (pi < row_end_i && pj < row_end_j) { // calculate sum
                int col_i = H_csr->col_idx[pi];
                int col_j = H_csr->col_idx[pj];

                if (col_i >= j || col_j >= j) break; // skip if col is outside bandwidth
                if (col_i < k_lower_bound || col_j < k_lower_bound) {
                    // jump to k_lower_bound
                    if (col_i < col_j) { pi++; continue; }
                    if (col_j < col_i) { pj++; continue; }
                }

                if (col_i == col_j) {
                    sum += L_dense[i * n + col_i] * L_dense[j * n + col_j];
                    pi++; pj++;
                } else if (col_i < col_j) {
                    pi++;
                } else {
                    pj++;
                }
            }

            if (i == j) {
                float val = H_csr->values[idx_j] - sum;
                if (val < -1e-12) return -1; // avoid zero division
                if (val < 0) val = 0;
                L_dense[i * n + j] = sqrtf(val);
            } else {
                float diag_j = L_dense[j * n + j];
                if (fabsf(diag_j) < 1e-12) return -1; // avoid zero division
                L_dense[i * n + j] = (H_csr->values[idx_j] - sum) / diag_j;
            }
        }
    }
    return 0;
}

// slove linear system
void forward_substitution(float *L, float *b, float *y, int n){
    for(int i=0; i<n; i++){
        float sum = 0.f;
        for(int j=0; j<i; j++){
            sum += L[IDX(i,j,n)] * y[j];
        }
        y[i] = (b[i] - sum) / L[IDX(i,i,n)];
    }
}

void backward_substitution(float *L, float *y, float *x, int n){
    for(int i=n-1; i>=0; i--){
        float sum = 0.f;
        for(int j=i+1; j<n; j++){
            sum += L[IDX(j,i,n)] * x[j];
        }
        x[i] = (y[i] - sum) / L[IDX(i,i,n)];
    }
}



int main(){
    printf("Solving H * delta_x = -b with matrix size %d\n", N);

    //float *H_new = malloc(N * N * sizeof(float));
    //float *b_new = malloc(N * sizeof(float));
    float *L = malloc(N * N * sizeof(float));
    float *y = malloc(N * sizeof(float));
    float *delta_x1 = malloc(N * sizeof(float));
    float *delta_x2 = malloc(N * sizeof(float));
    float *H2 = malloc(N * N * sizeof(float));
    float *b2 = malloc(N * sizeof(float));
    int *order = malloc(N * sizeof(int));

    if( !L || !y || !delta_x1 || !delta_x2 || !H2 || !b2 || !order){
        printf("Memory allocation failed!\n");
        return -1;
    }

    int ret;

    // generate adjacent matrix
    int *adj_orig = malloc(N*N*sizeof(int));
    if(!adj_orig){
        printf("Memory alloc failed\n");
        return -1;
    }
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            adj_orig[i*N+j] = (fabs(H[i][j]) > 1e-4) ? 1 : 0;
        }
    }

    // generate CSR matrix from H, compress
    CSRMatrix *csr_H = dense_to_csr(H, N);
    if(!csr_H){
        printf("CSR malloc failed\n");
        free(adj_orig);
        return -1;
    }
        //bandwith check
    int original_bandwidth = compute_bandwidth(NULL, N, csr_H);
    printf("original bandwidth of H: %d\n",original_bandwidth);

    // way 1: directly cholesky
        float *neg_b = malloc(N * sizeof(float));
        if(!neg_b){
            printf("Memory alloc failed\n");
            free_csr(csr_H);
            free(adj_orig);
            return -1;
        }
        for(int i=0; i<N; i++) neg_b[i] = -b[i];

        Sleep(1000);

        LARGE_INTEGER freq1, start1, end1;
        QueryPerformanceFrequency(&freq1);
        QueryPerformanceCounter(&start1);

        ret = cholesky_csr_plain(csr_H, L);

        QueryPerformanceCounter(&end1);
        double time_ms1 = (double)(end1.QuadPart - start1.QuadPart) * 1000.0 / freq1.QuadPart;

        if(ret != 0){
            printf("Cholesky CSR plain failed!\n");
            free(neg_b);
            free_csr(csr_H);
            free(adj_orig);
            goto cleanup;
        }
        LARGE_INTEGER freq2, start2, end2;
        QueryPerformanceFrequency(&freq2);
        QueryPerformanceCounter(&start2);

        forward_substitution(L, neg_b, y, N);
        backward_substitution(L, y, delta_x1, N);
        
        QueryPerformanceCounter(&end2);
        double time_ms2 = (double)(end2.QuadPart - start2.QuadPart) * 1000.0 / freq2.QuadPart;

        free(neg_b);
    
    Sleep(2000);  
        
    // way 2: RCM -> cholesky
        int *adj_reordered = malloc(N*N*sizeof(int));
        for(int i=0; i<N; i++) b2[i] = -b2[i];


        LARGE_INTEGER freq3, start3, end3;
        QueryPerformanceFrequency(&freq3);
        QueryPerformanceCounter(&start3);

        rcm(adj_orig, N, order);
        reorder_matrix_vector(H, b, H2, b2, order, N);

        QueryPerformanceCounter(&end3);
        double time_ms3 = (double)(end3.QuadPart - start3.QuadPart) * 1000.0 / freq3.QuadPart;
        
        if(!adj_reordered){
            printf("Memory alloc failed\n");
            free_csr(csr_H);
            free(adj_orig);
            goto cleanup;
        }
        reorder_adj_matrix(adj_orig, adj_reordered, order, N);
        CSRMatrix *csr_H2 = dense_to_csr(H2, N);
        if(!csr_H2){
            printf("CSR malloc failed\n");
            free(adj_reordered);
            free_csr(csr_H);
            free(adj_orig);
            goto cleanup;
        }

        int bandwidth = compute_bandwidth(order, N, csr_H2);
        printf("Computed bandwidth after RCM: %d\n", bandwidth);

        LARGE_INTEGER freq4, start4, end4;
        QueryPerformanceFrequency(&freq4);
        QueryPerformanceCounter(&start4);

        ret = cholesky_csr_bandwidth(csr_H2, L, bandwidth);

        QueryPerformanceCounter(&end4);
        double time_ms4 = (double)(end4.QuadPart - start4.QuadPart) * 1000.0 / freq4.QuadPart;

        if(ret != 0){
            printf("Cholesky CSR with bandwidth failed!\n");
            free_csr(csr_H2);
            free(adj_reordered);
            free_csr(csr_H);
            free(adj_orig);
            goto cleanup;
        }

        LARGE_INTEGER freq5, start5, end5;
        QueryPerformanceFrequency(&freq5);
        QueryPerformanceCounter(&start5);

        forward_substitution(L, b2, y, N);
        backward_substitution(L, y, delta_x2, N);

        QueryPerformanceCounter(&end5);
        double time_ms5 = (double)(end5.QuadPart - start5.QuadPart) * 1000.0 / freq5.QuadPart;

        float *delta_x2_inv = malloc(N * sizeof(float));
        if(!delta_x2_inv){
            printf("Memory alloc failed\n");
            free_csr(csr_H2);
            free(adj_reordered);
            free_csr(csr_H);
            free(adj_orig);
            goto cleanup;
        }
        inverse_reorder_vector(delta_x2, delta_x2_inv, order, N);

        free(delta_x2_inv);
        free_csr(csr_H2);
        free(adj_reordered);

        printf("\nMethon 1: directly CSR Cholesky decomposition -> solve linear system\n");
        printf("        1A - Cholesky decomposition ................................  %.3f ms\n", time_ms1);
        printf("        1B - Solve linear system ...................................  %.3f ms\n", time_ms2);
        printf("        Method 1 take total time ............................................ %.3f ms\n", time_ms1 + time_ms2);
        printf("\nMethon 2: RCM -> Cholesky decomposition -> solve linear system\n");
        printf("        2A - time (RCM):............................................  %.3f ms\n", time_ms3);
        printf("        2B - time (CSR Cholesky decomposition by RCM constrain):....  %.3f ms\n", time_ms4);
        printf("        2C - solve linear system ...................................  %.3f ms\n", time_ms5);
        printf("        Method 2 take total time ...........................................  %.3f ms\n", time_ms3 + time_ms4 + time_ms5);
    

cleanup:
    //free(H);
    //free(b);
    free(L);
    free(y);
    free(delta_x1);
    free(delta_x2);
    free(H2);
    free(b2);
    free(order);
    free_csr(csr_H);
    free(adj_orig);
    return 0;
}
