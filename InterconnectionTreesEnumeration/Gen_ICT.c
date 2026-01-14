#define _POSIX_C_SOURCE 200809L

#include <limits.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>

static volatile sig_atomic_t time_limit_hit = 0;
long long int compteurSol = 0;

#ifdef WEIGHTED
struct InterconnectionTree {
    int *edges;
    double total_weight;
};

static struct InterconnectionTree *g_interconnection_trees = NULL;
static long long g_tree_capacity = 0;
static struct timespec start_stock, end_stock;
static double elapsed_stock_ms = 0.0;

static void ensure_tree_capacity(long long required) {
    if (required <= g_tree_capacity) {
        return;
    }

    long long new_capacity = g_tree_capacity > 0 ? g_tree_capacity : 1024;
    while (new_capacity < required) {
        if (new_capacity > LLONG_MAX / 2) {
            fprintf(stderr, "Tree capacity overflow\n");
            exit(EXIT_FAILURE);
        }
        new_capacity *= 2;
    }

    struct InterconnectionTree *tmp = realloc(g_interconnection_trees, (size_t)new_capacity * sizeof(*tmp));
    if (tmp == NULL) {
        fprintf(stderr, "Failed to allocate memory for %lld trees\n", new_capacity);
        exit(EXIT_FAILURE);
    }

    g_interconnection_trees = tmp;
    g_tree_capacity = new_capacity;
}

static void store_interconnection_tree(const int *edges, int edge_slot_count, double total_weight) {
    ensure_tree_capacity(compteurSol + 1);
    struct InterconnectionTree *tree = &g_interconnection_trees[compteurSol];

    clock_gettime(CLOCK_MONOTONIC, &start_stock);
    tree->edges = malloc((size_t)edge_slot_count * sizeof(int));
    if (tree->edges == NULL) {
        fprintf(stderr, "Failed to allocate memory for tree edges\n");
        exit(EXIT_FAILURE);
    }
    memcpy(tree->edges, edges, (size_t)edge_slot_count * sizeof(int));
    tree->total_weight = total_weight;
    clock_gettime(CLOCK_MONOTONIC, &end_stock);
    elapsed_stock_ms += (end_stock.tv_sec - start_stock.tv_sec) * 1000.0 +
                        (end_stock.tv_nsec - start_stock.tv_nsec) / 1e6;
}

static void free_weighted_resources(long long tree_count) {
    if (g_interconnection_trees == NULL) {
        return;
    }

    for (long long i = 0; i < tree_count; ++i) {
        free(g_interconnection_trees[i].edges);
        g_interconnection_trees[i].edges = NULL;
    }

    free(g_interconnection_trees);
    g_interconnection_trees = NULL;
    g_tree_capacity = 0;
}

static void merge_interconnection_trees(struct InterconnectionTree *arr, struct InterconnectionTree *tmp, long long left, long long mid, long long right) {
    long long i = left;
    long long j = mid;
    long long k = left;

    while (i < mid && j < right) {
        if (arr[i].total_weight <= arr[j].total_weight) {
            tmp[k++] = arr[i++];
        } else {
            tmp[k++] = arr[j++];
        }
    }

    while (i < mid) {
        tmp[k++] = arr[i++];
    }

    while (j < right) {
        tmp[k++] = arr[j++];
    }

    for (long long idx = left; idx < right; ++idx) {
        arr[idx] = tmp[idx];
    }
}

static void merge_sort_interconnection_trees(struct InterconnectionTree *arr, struct InterconnectionTree *tmp, long long left, long long right) {
    if (right - left <= 1) {
        return;
    }

    long long mid = left + (right - left) / 2;
    merge_sort_interconnection_trees(arr, tmp, left, mid);
    merge_sort_interconnection_trees(arr, tmp, mid, right);
    merge_interconnection_trees(arr, tmp, left, mid, right);
}

#endif

static void handle_time_limit(int sig) {
    (void)sig;
    time_limit_hit = 1;
}

static void load_partition(int** S, int** C, int* num_vertex, int* num_components, const char* input_file) {
    FILE *file = fopen(input_file, "r");
    if (file == NULL) {
        perror("Failed to open input file");
        exit(EXIT_FAILURE);
    }

    if (fscanf(file, "%d %d", num_vertex, num_components) != 2) {
        fprintf(stderr, "Malformed header in %s\n", input_file);
        fclose(file);
        exit(EXIT_FAILURE);
    }

    if (*num_vertex <= 0 || *num_components <= 0) {
        fprintf(stderr, "Invalid sizes in %s: vertices=%d components=%d\n", input_file, *num_vertex, *num_components);
        fclose(file);
        exit(EXIT_FAILURE);
    }

    *S = malloc((size_t)(*num_vertex) * sizeof(int));
    if (*S == NULL) {
        fprintf(stderr, "Failed to allocate vertex-to-component array\n");
        fclose(file);
        exit(EXIT_FAILURE);
    }

    *C = calloc((size_t)(*num_components), sizeof(int));
    if (*C == NULL) {
        fprintf(stderr, "Failed to allocate component sizes\n");
        fclose(file);
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < *num_vertex; ++i) {
        int component = -1;
        if (fscanf(file, "%d", &component) != 1) {
            fprintf(stderr, "Malformed component assignment at vertex %d in %s\n", i, input_file);
            fclose(file);
            exit(EXIT_FAILURE);
        }

        if (component < 0 || component >= *num_components) {
            fprintf(stderr, "Component index %d out of bounds for %s\n", component, input_file);
            fclose(file);
            exit(EXIT_FAILURE);
        }

        (*S)[i] = component;
        (*C)[component]++;
    }

    fclose(file);
}

// Update vertex in place with a backup
static int UpdateVertex(int* S, int S1, int S2, int num_vertex, int* backup_S1, int* backup_S2, int* affected_indices, int* affected_count) {
    int min, max;

    if (S[S1] < S[S2]) {
        min = S[S1];
        max = S[S2];
    } else {
        min = S[S2];
        max = S[S1];
    }

    *backup_S1 = S[S1];
    *backup_S2 = S[S2];

    // Store the indices of vertices that are modified
    *affected_count = 0;
    for (int i = 0; i < num_vertex; i++) {
        if (S[i] == max) {
            affected_indices[*affected_count] = i;
            (*affected_count)++;
            S[i] = min;  // Replace max with min
        }
    }

    // Mark the selected vertices as -1
    S[S1] = -1;
    S[S2] = -1;

    return max; // Return the original max value for restoring
}

// Restore the vertex array after modification
static void RestoreVertex(int* S, int S1, int S2, int backup_S1, int backup_S2, int* affected_indices, int affected_count, int backup_max) {
    // Restore the vertices S1 and S2
    S[S1] = backup_S1;
    S[S2] = backup_S2;

    // Restore the affected vertices that were set to min back to max
    for (int i = 0; i < affected_count; i++) {
        S[affected_indices[i]] = backup_max;
    }
}

// Update component array in place with a backup
static void UpdateSubset(int* C, int P1, int P2, int* backup_C_max, int* backup_C_min) {
    int min, max;

    if (P1 < P2) {
        min = P1;
        max = P2;
    } else {
        min = P2;
        max = P1;
    }

    *backup_C_min = C[min];
    *backup_C_max = C[max];

    C[min] = C[min] + C[max] - 2;
    C[max] = 0;
}

// Restore the component array after modification
static void RestoreSubset(int* C, int P1, int P2, int backup_C_max, int backup_C_min) {
    int min, max;

    if (P1 < P2) {
        min = P1;
        max = P2;
    } else {
        min = P2;
        max = P1;
    }

    C[min] = backup_C_min;
    C[max] = backup_C_max;
}

static void Ban_s_Vertex(int* S, int s, int* backup_S) {
    *backup_S = S[s];
    S[s] = -1;
}

static void RestoreBan_s_Vertex(int* S, int s, int backup_S) {
    S[s] = backup_S;
}

static void Ban_s_Subset(int* C, int p, int* backup_C) {
    *backup_C = C[p];
    C[p]--;
}

static void RestoreBan_s_Subset(int* C, int p, int backup_C) {
    C[p] = backup_C;
}

#ifdef WEIGHTED
static void EnumArbresInterconnexion(int* S, int* C, int* t, int num_vertex, int k, int K, int s, int marge, double (*edge_weight)[num_vertex], double total_weight) {
    if (time_limit_hit) {
        return;
    }

    if (K >= (k - 1)) {
        store_interconnection_tree(t, 2 * (k - 1), total_weight);
        compteurSol++;
        return;
    }

    int next_s;
    for (int v = s + 1; v < num_vertex; v++) {
        int comp_sum = C[S[s]] + C[S[v]] - 2;
        if ((S[s] != -1) && (S[v] != -1) && (S[s] != S[v]) && ((K == (k - 2)) || (comp_sum > 0))) {
            t[K * 2] = s;
            t[K * 2 + 1] = v;

            int backup_S1, backup_S2, backup_max;
            int affected_indices[num_vertex];
            int affected_count = 0;
            int backup_C_max, backup_C_min;

            UpdateSubset(C, S[s], S[v], &backup_C_max, &backup_C_min);
            backup_max = UpdateVertex(S, s, v, num_vertex, &backup_S1, &backup_S2, affected_indices, &affected_count);

            next_s = s;
            do {
                next_s = next_s + 1;
            } while (next_s < num_vertex && S[next_s] == -1);

            EnumArbresInterconnexion(S, C, t, num_vertex, k, K + 1, next_s, marge, edge_weight, total_weight + edge_weight[s][v]);

            RestoreVertex(S, s, v, backup_S1, backup_S2, affected_indices, affected_count, backup_max);
            RestoreSubset(C, S[s], S[v], backup_C_max, backup_C_min);
        }
    }

    next_s = s;
    if ((marge > 0) && (C[S[s]] > 0)) {
        int backup_S, backup_C;
        Ban_s_Subset(C, S[s], &backup_C);
        Ban_s_Vertex(S, s, &backup_S);
        do {
            next_s = next_s + 1;
        } while (next_s < num_vertex && S[next_s] == -1);
        EnumArbresInterconnexion(S, C, t, num_vertex, k, K, next_s, marge - 1, edge_weight, total_weight);
        RestoreBan_s_Vertex(S, s, backup_S);
        RestoreBan_s_Subset(C, S[s], backup_C);
    }
}
#else
static void EnumArbresInterconnexion(int* S, int* C, int* t, int num_vertex, int k, int K, int s, int marge) {
    if (time_limit_hit) {
        return;
    }

    if (K >= (k - 1)) {
        compteurSol++;
        return;
    }

    int next_s;
    for (int v = s + 1; v < num_vertex; v++) {
        int comp_sum = C[S[s]] + C[S[v]] - 2;
        if ((S[s] != -1) && (S[v] != -1) && (S[s] != S[v]) && ((K == (k - 2)) || (comp_sum > 0))) {
            t[K * 2] = s;
            t[K * 2 + 1] = v;

            int backup_S1, backup_S2, backup_max;
            int affected_indices[num_vertex];
            int affected_count = 0;
            int backup_C_max, backup_C_min;

            UpdateSubset(C, S[s], S[v], &backup_C_max, &backup_C_min);
            backup_max = UpdateVertex(S, s, v, num_vertex, &backup_S1, &backup_S2, affected_indices, &affected_count);

            next_s = s;
            do {
                next_s = next_s + 1;
            } while (next_s < num_vertex && S[next_s] == -1);

            EnumArbresInterconnexion(S, C, t, num_vertex, k, K + 1, next_s, marge);

            RestoreVertex(S, s, v, backup_S1, backup_S2, affected_indices, affected_count, backup_max);
            RestoreSubset(C, S[s], S[v], backup_C_max, backup_C_min);
        }
    }

    next_s = s;
    if ((marge > 0) && (C[S[s]] > 0)) {
        int backup_S, backup_C;
        Ban_s_Subset(C, S[s], &backup_C);
        Ban_s_Vertex(S, s, &backup_S);
        do {
            next_s = next_s + 1;
        } while (next_s < num_vertex && S[next_s] == -1);
        EnumArbresInterconnexion(S, C, t, num_vertex, k, K, next_s, marge - 1);
        RestoreBan_s_Vertex(S, s, backup_S);
        RestoreBan_s_Subset(C, S[s], backup_C);
    }
}
#endif


int main(int argc, char** argv) {
    if (argc < 2 || argc > 3) {
        fprintf(stderr, "Usage: %s <input_file> [time_limit_seconds]\n", argv[0]);
        return EXIT_FAILURE;
    }

    double time_limit_sec = 0.0;
    if (argc == 3) {
        time_limit_sec = atof(argv[2]);
        if (time_limit_sec <= 0.0) {
            fprintf(stderr, "time limit must be > 0\n");
            return EXIT_FAILURE;
        }

        struct sigaction sa;
        memset(&sa, 0, sizeof(sa));
        sa.sa_handler = handle_time_limit;
        sigaction(SIGALRM, &sa, NULL);

        struct itimerval timer;
        memset(&timer, 0, sizeof(timer));
        timer.it_value.tv_sec = (time_t)time_limit_sec;
        timer.it_value.tv_usec = (suseconds_t)((time_limit_sec - timer.it_value.tv_sec) * 1e6);
        setitimer(ITIMER_REAL, &timer, NULL);
    }

    time_limit_hit = 0;
    compteurSol = 0;
    // Example usage:
    int num_vertex;
    int num_components;   
    
    int *S = NULL;        // Array for storing vertex to component mapping
    int *C = NULL;    // Array for storing the count of vertices per component

    char* input_file = argv[1];

    // Load the partition description
    load_partition(&S, &C, &num_vertex, &num_components, input_file);

    //find min and max in C
    int min_size = C[0];
    int max_size = C[0];
    for(int i=1;i<num_components;i++){
        if(C[i]<min_size){
            min_size = C[i];
        }
        if(C[i]>max_size){
            max_size = C[i];
        }
    }

   // printf("vertex %d, components %d, min %d, max %d\n",num_vertex,num_components, min_size, max_size);
    
    // Print S and C arrays to verify the results
    // printf("S array (vertex to component mapping):\n");
    // for (int i = 0; i < num_vertex; i++) {
    //     printf("Vertex %d: Component %d\n", i, S[i]);
    // }
    // printf("\nC array (component vertex count):\n");
    // for (int i = 0; i < num_components; i++) {
    //     printf("Component %d: %d vertices\n", i, C[i]);
    // }

    #ifdef WEIGHTED
    double edge_weight[num_vertex][num_vertex];
    srand((unsigned)time(NULL));
    for (int i = 0; i < num_vertex; i++) {
        for (int j = 0; j < num_vertex; j++) {
            edge_weight[i][j] = ((double)(rand() % 1000) / 1000.0) * 9.9 + 0.1;
        }
    }
    struct timespec start_sort = {0}, end_sort = {0};
    double elapsed_sort_ms = 0.0;
    elapsed_stock_ms = 0.0;
    #endif

    size_t edge_slot_count = (size_t)(2 * (num_components - 1));
    int* t = edge_slot_count ? malloc(edge_slot_count * sizeof(int)) : NULL;
    if (edge_slot_count > 0 && t == NULL) {
        fprintf(stderr, "Failed to allocate edge buffer\n");
        free(S);
        free(C);
        return EXIT_FAILURE;
    }

 //  printf("\n Start Algo\n");
 //   printf("marge %d \n", num_vertex - (2 * (num_components - 1)));

    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);
    #ifdef WEIGHTED
        EnumArbresInterconnexion(S, C, t, num_vertex, num_components, 0, 0, num_vertex - (2 * (num_components - 1)), edge_weight, 0.0);
        if (compteurSol > 1) {
            struct InterconnectionTree* tmp = malloc((size_t)compteurSol * sizeof(*tmp));
            if (tmp == NULL) {
                fprintf(stderr, "Failed to allocate buffer for sorting\n");
                free(t);
                free(S);
                free(C);
                free_weighted_resources(compteurSol);
                return EXIT_FAILURE;
            }
            clock_gettime(CLOCK_MONOTONIC, &start_sort);
            merge_sort_interconnection_trees(g_interconnection_trees, tmp, 0, compteurSol);
            clock_gettime(CLOCK_MONOTONIC, &end_sort);
            elapsed_sort_ms = (end_sort.tv_sec - start_sort.tv_sec) * 1000.0 +
                              (end_sort.tv_nsec - start_sort.tv_nsec) / 1e6;
            free(tmp);
        }
    #else
        EnumArbresInterconnexion(S, C, t, num_vertex, num_components, 0, 0, num_vertex - (2 * (num_components - 1)));
    #endif
    clock_gettime(CLOCK_MONOTONIC, &end);
    double elapsed_ms = (end.tv_sec - start.tv_sec) * 1000.0 + (end.tv_nsec - start.tv_nsec) / 1e6;

    //printf("Nombre de solutions : %d\n", compteurSol);
    int timed_out = time_limit_hit;
    if (time_limit_sec > 0.0 && !timed_out) {
        struct itimerval timer = {0};
        setitimer(ITIMER_REAL, &timer, NULL);
    }

    double avg_delay_ms = (compteurSol > 0) ? (elapsed_ms / compteurSol) : elapsed_ms;

    #ifdef WEIGHTED
    printf("%s,%d,%d,%d,%d,%lld,%.3f,%.6f,%s,%.6f,%.6f\n",
           input_file, num_vertex, num_components,
           min_size, max_size, compteurSol,
           elapsed_ms, avg_delay_ms, timed_out ? "timeout" : "ok", elapsed_sort_ms, elapsed_stock_ms);
    #else
    printf("%s,%d,%d,%d,%d,%lld,%.3f,%.6f,%s\n",
           input_file, num_vertex, num_components,
           min_size, max_size, compteurSol,
           elapsed_ms, avg_delay_ms, timed_out ? "timeout" : "ok");
    #endif    
        

    // Free allocated memory
    free(S);
    free(C);
    free(t);
#ifdef WEIGHTED
    free_weighted_resources(compteurSol);
#endif
    
    return 0;
}
