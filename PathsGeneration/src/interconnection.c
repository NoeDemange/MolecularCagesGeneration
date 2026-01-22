#include <float.h>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "assembly.h"
#include "constant.h"
#include "discretization.h"
#include "distance.h"
#include "interconnection.h"
#include "util.h"

/**
 * @file interconnection.c
 * @brief Functions for enumerating interconnection trees.
 */

// FIND CONNECTED COMPONENTS
/**
 * @brief Performs a depth-first search (DFS) to identify connected components in a molecular cage.
 *
 * This function traverses the molecular cage starting from a given atom index, marking all atoms
 * in the same connected component with the corresponding component number. It explores neighbors
 * recursively using DFS.
 *
 * @param s Pointer to the `Cage_t` structure representing the molecular cage.
 * @param index Starting atom index for the DFS traversal.
 * @param V Pointer to an array where each element corresponds to an atom and stores the component number.
 *          Unvisited atoms are marked with `-1`.
 * @param num_components The current component number used to label the connected atoms.
 *
 * @details
 * - Marks the current atom with the `num_components` label in the `V` array.
 * - Recursively visits all unvisited neighbors of the current atom.
 * - Assumes that the cage is represented as a graph where atoms are nodes and bonds are edges.
 *
 * ### Usage Notes:
 * - This function is used internally by `connectedComponents` and assumes that `V` is pre-initialized
 *   with `-1` for all atoms.
 *
 * @see connectedComponents
 */
void dfSearch(Cage_t *s, const int INDEX, int *V, const int NUMCOMPONENTS) {
  AtomCage_t *a = atom(s, INDEX);
  V[INDEX] = NUMCOMPONENTS;
  for (int i = 0; forEachNeighbor(a, i); i++) {
    if (V[neighbor(a, i)] == -1) {
      dfSearch(s, neighbor(a, i), V, NUMCOMPONENTS);
    }
  }
}

/**
 * @brief Counts the number of connected components in a molecular cage.
 *
 * This function identifies and counts the connected components of a molecular cage
 * using a depth-first search. Each atom is assigned a component number in the `V` array,
 * where atoms in the same connected component share the same number.
 *
 * @param s Pointer to the `Cage_t` structure representing the molecular cage.
 * @param V Pointer to an array where each element corresponds to an atom and stores the
 *          component number. Unvisited atoms are initialized to `-1`.
 * @param num_components Pointer to an integer that stores the total number of connected components
 *                      in the cage after execution.
 *
 * @details
 * - Iterates through all atoms in the cage and starts a new DFS traversal for each unvisited atom
 *   with a valid flag (i.e., not `NOT_DEF_F`).
 * - Each DFS traversal marks all atoms in the same connected component with the current component number.
 * - Increments the `num_components` counter for each new connected component found.
 *
 * ### Usage Notes:
 * - The `V` array must be pre-allocated with a size equal to the number of atoms in the cage and
 *   initialized to `-1` for all elements.
 * - The function updates `num_components` to include the total number of connected components.
 *
 * @see dfSearch
 */
void connectedComponents(Cage_t *s, int *V, int *num_components) {
  ////Initialisation
  (*num_components) = -1;
  ////Run
  for (int i = 0; i < size(s); i++) {
    if (V[i] == -1 && (flag(atom(s, i)) != NOT_DEF_F)) {
      (*num_components)++;
      dfSearch(s, i, V, (*num_components));
    }
  }
  ////Adapt for table
  (*num_components)++;
}

// ENUMERATION

int cpt_inter_tree_egv = 0; // Extern global variable to count the number of interconnection trees

typedef struct {
  int *edges;        // Copied edge list (pairs of vertices)
  double total_size; // Sum of distances between the start and end of each path
} StoredInterconnectionTree;

typedef struct {
  StoredInterconnectionTree *items;
  size_t count;
  size_t capacity;
  int edge_slot_count;
} InterconnectionTreeStore;

/**
 * @brief Initialize the bookkeeping structure that caches candidate trees.
 *
 * @param store Store to initialize; fields are reset to empty values.
 * @param edge_slot_count Number of integers required to encode one tree (2 * edges).
 */
static void initTreeStore(InterconnectionTreeStore *store, int edge_slot_count) {
  store->items = NULL;
  store->count = 0;
  store->capacity = 0;
  store->edge_slot_count = edge_slot_count;
}

/**
 * @brief Release all memory held by the tree store.
 *
 * @param store Store to destroy; safe to pass NULL.
 */
static void freeTreeStore(InterconnectionTreeStore *store) {
  if (!store) {
    return;
  }
  for (size_t i = 0; i < store->count; i++) {
    free(store->items[i].edges);
    store->items[i].edges = NULL;
  }
  free(store->items);
  store->items = NULL;
  store->count = 0;
  store->capacity = 0;
  store->edge_slot_count = 0;
}

/**
 * @brief Grow the storage so it can hold at least @p required entries.
 *
 * @param store Target store.
 * @param required Desired number of slots.
 */
static void ensureTreeStoreCapacity(InterconnectionTreeStore *store, size_t required) {
  if (required <= store->capacity) {
    return;
  }

  size_t new_capacity = store->capacity ? store->capacity : 64;
  while (new_capacity < required) {
    if (new_capacity > (SIZE_MAX / 2)) {
      fprintf(stderr, "Tree storage capacity overflow\n");
      exit(EXIT_FAILURE);
    }
    new_capacity *= 2;
  }

  StoredInterconnectionTree *tmp = realloc(store->items, new_capacity * sizeof(*tmp));
  if (!tmp) {
    fprintf(stderr, "Unable to allocate memory for %zu interconnection trees\n", new_capacity);
    exit(EXIT_FAILURE);
  }
  store->items = tmp;
  store->capacity = new_capacity;
}

/**
 * @brief Append a new tree snapshot to the store.
 *
 * @param store Destination store (must be initialized).
 * @param edges Array describing the tree edges (length edge_slot_count).
 * @return Pointer to the stored entry for further mutation.
 */
static StoredInterconnectionTree *appendTree(InterconnectionTreeStore *store, const int *edges) {
  if (!store || store->edge_slot_count <= 0) {
    return NULL;
  }
  ensureTreeStoreCapacity(store, store->count + 1);
  StoredInterconnectionTree *entry = &store->items[store->count++];
  size_t bytes = (size_t)store->edge_slot_count * sizeof(int);
  entry->edges = malloc(bytes);
  if (!entry->edges) {
    fprintf(stderr, "Unable to allocate memory for interconnection tree edges\n");
    exit(EXIT_FAILURE);
  }
  memcpy(entry->edges, edges, bytes);
  entry->total_size = 0.0;
  return entry;
}

/**
 * @brief Compare two stored trees by their aggregated length.
 */
static int compareStoredTrees(const void *lhs, const void *rhs) {
  const StoredInterconnectionTree *a = (const StoredInterconnectionTree *)lhs;
  const StoredInterconnectionTree *b = (const StoredInterconnectionTree *)rhs;

  if (a->total_size < b->total_size) {
    return -1;
  }
  if (a->total_size > b->total_size) {
    return 1;
  }
  return 0;
}

/**
 * @brief Compute the distance between two candidate vertices using the active metric.
 *
 * @param start_id Cage vertex index acting as start.
 * @param end_id Cage vertex index acting as end.
 * @param cage Cage definition for coordinate queries.
 * @param paths Paths workspace (contains grids/heaps for A* variants).
 * @param grid_sub Grid used for collision checks in hybrid modes.
 * @param substrat_t Substrat atom positions.
 * @return Distance cost; DBL_MAX when no valid path exists.
 */
static double computePairDistance(int start_id, int end_id, Cage_t *cage, Paths_t *paths, GridSubstrat *grid_sub,
                                  double ***substrat_t) {
  Point_t start = coords(atom(cage, start_id));
  Point_t end = coords(atom(cage, end_id));
  DistanceType type = get_current_distance_type();
  if (type == DISTANCE_HYBRID) {
    type = distanceHybridLineSpheresCollisionTest(start, end, cage, paths, grid_sub, substrat_t);
  }

  if (type == DISTANCE_EUCLIDEAN) {
    return dist(start, end);
  }

  if ((type == DISTANCE_A_STAR || type == DISTANCE_SSMTA_STAR) && paths && paths->grids && paths->minHeaps &&
      paths->grids[0] && paths->minHeaps[0]) {
    double cost = aStarDistance(start, end, paths->grids[0], paths->minHeaps[0]);
    if (cost >= 0.0) {
      return cost;
    }
    return DBL_MAX;
  }

  return dist(start, end);
}

/**
 * @brief Sum the pair distances for every edge within a stored tree.
 *
 * @param tree Tree snapshot to evaluate.
 * @param num_paths Number of edges in the tree.
 * @param cage Cage definition for coordinates.
 * @param paths Paths workspace.
 * @param grid_sub Grid for hybrid collision checks.
 * @param substrat_t Substrat atom positions.
 * @return Aggregated size or DBL_MAX when any segment is invalid.
 */
static double computeTreeTotalSize(StoredInterconnectionTree *tree, int num_paths, Cage_t *cage, Paths_t *paths,
                                   GridSubstrat *grid_sub, double ***substrat_t) {
  if (!tree || !tree->edges || num_paths <= 0) {
    return 0.0;
  }

  double total = 0.0;
  for (int i = 0; i < num_paths; i++) {
    double edge_size = computePairDistance(tree->edges[i * 2], tree->edges[i * 2 + 1], cage, paths, grid_sub,
                                           substrat_t);
    if (edge_size == DBL_MAX) {
      return DBL_MAX;
    }
    total += edge_size;
  }
  return total;
}

/**
 * @brief Check whether a stored tree reuses any edge marked as banned.
 *
 * @param edges Edge array (pairs of vertex indices).
 * @param num_paths Number of edges in the tree.
 * @param banned_edges Flattened list of banned start/end vertex pairs.
 * @param banned_count Number of banned edges.
 * @return 1 when the tree contains a banned edge, 0 otherwise.
 */
static int treeContainsBannedEdge(const int *edges, int num_paths, const int *banned_edges, int banned_count) {
  if (!edges || !banned_edges || num_paths <= 0 || banned_count <= 0) {
    return 0;
  }
  for (int path_idx = 0; path_idx < num_paths; path_idx++) {
    int start = edges[path_idx * 2];
    int end = edges[path_idx * 2 + 1];
    for (int b = 0; b < banned_count; b++) {
      if (banned_edges[b * 2] == start && banned_edges[b * 2 + 1] == end) {
        return 1;
      }
    }
  }
  return 0;
}

/**
 * @brief Updates the vertex sets by merging the sets of two vertices and propagating changes.
 *
 * This function merges the two sets of vertices (S1 and S2) by replacing their associated
 * values with the minimum of the two. It also modifies other vertices connected to the
 * vertex with the higher value by replacing their values accordingly. The original values
 * of S1 and S2 are backed up, and the indices of affected vertices are stored.
 *
 * @param tab_s The array representing the vertex set, where each vertex is associated with two values.
 * @param S1 The index of the first vertex to merge.
 * @param S2 The index of the second vertex to merge.
 * @param num_vertex The total number of vertices in the set.
 * @param backup_s1 A pointer to store the original value of vertex S1.
 * @param backup_s2 A pointer to store the original value of vertex S2.
 * @param affected_indices An array to store the indices of vertices affected by the update.
 * @param affected_count A pointer to an integer that will store the number of affected vertices.
 *
 * @return The original maximum value of the two merged vertices (before the update).
 */
int updateVertex(int *tab_s, int S1, int S2, int num_vertex, int *backup_s1, int *backup_s2, int *affected_indices,
                 int *affected_count) {
  int min, max;

  if (tab_s[S1 * 2 + 1] < tab_s[S2 * 2 + 1]) {
    min = tab_s[S1 * 2 + 1];
    max = tab_s[S2 * 2 + 1];
  } else {
    min = tab_s[S2 * 2 + 1];
    max = tab_s[S1 * 2 + 1];
  }

  *backup_s1 = tab_s[S1 * 2 + 1];
  *backup_s2 = tab_s[S2 * 2 + 1];

  // Store the indices of vertices that are modified
  *affected_count = 0;
  for (int i = 0; i < num_vertex; i++) {
    if (tab_s[i * 2 + 1] == max) {
      affected_indices[*affected_count] = i;
      (*affected_count)++;
      tab_s[i * 2 + 1] = min; // Replace max with min
    }
  }

  // Mark the selected vertices as -1
  tab_s[S1 * 2 + 1] = -1;
  tab_s[S2 * 2 + 1] = -1;

  return max; // Return the original max value for restoring
}

/**
 * @brief Restores the vertex array to its original state after modification.
 *
 * This function restores the values of two vertices (S1 and S2) and the affected vertices
 * that were modified during an update operation. It uses the backed-up values to restore the original state.
 *
 * @param tab_s The array representing the vertex set, where each vertex is associated with two values.
 * @param S1 The index of the first vertex that was modified.
 * @param S2 The index of the second vertex that was modified.
 * @param backup_s1 The original value of vertex S1 before modification.
 * @param backup_s2 The original value of vertex S2 before modification.
 * @param affected_indices An array containing the indices of the vertices that were affected by the update.
 * @param affected_count The number of vertices in the affected_indices array.
 * @param backup_max The original maximum value of the affected vertices before modification.
 */
void restoreVertex(int *tab_s, int S1, int S2, int backup_s1, int backup_s2, int *affected_indices, int affected_count,
                   int backup_max) {
  // Restore the vertices S1 and S2
  tab_s[S1 * 2 + 1] = backup_s1;
  tab_s[S2 * 2 + 1] = backup_s2;

  // Restore the affected vertices that were set to min back to max
  for (int i = 0; i < affected_count; i++) {
    tab_s[affected_indices[i] * 2 + 1] = backup_max;
  }
}

/**
 * @brief Updates the component array in place and backs up the original values.
 *
 * This function updates the component array by merging two components (P1 and P2) and
 * modifies them in place. The original values of these components are backed up for restoration.
 *
 * @param tab_c The array representing the component set.
 * @param P1 The index of the first component to merge.
 * @param P2 The index of the second component to merge.
 * @param backup_c_max A pointer to store the original value of the larger component.
 * @param backup_c_min A pointer to store the original value of the smaller component.
 */
void updateSubset(int *tab_c, int P1, int P2, int *backup_c_max, int *backup_c_min) {
  int min, max;

  if (P1 < P2) {
    min = P1;
    max = P2;
  } else {
    min = P2;
    max = P1;
  }

  *backup_c_min = tab_c[min];
  *backup_c_max = tab_c[max];

  tab_c[min] = tab_c[min] + tab_c[max] - 2;
  tab_c[max] = 0;
}

/**
 * @brief Restores the component array to its original state after modification.
 *
 * This function restores the component array after an update operation using backed-up values.
 * It reverts the changes made during the update and restores the components at P1 and P2.
 *
 * @param tab_c The array representing the component set.
 * @param P1 The index of the first component that was modified.
 * @param P2 The index of the second component that was modified.
 * @param backup_c_max The original value of the larger component before modification.
 * @param backup_c_min The original value of the smaller component before modification.
 */
void restoreSubset(int *tab_c, int P1, int P2, int backup_c_max, int backup_c_min) {
  int min, max;

  if (P1 < P2) {
    min = P1;
    max = P2;
  } else {
    min = P2;
    max = P1;
  }

  tab_c[min] = backup_c_min;
  tab_c[max] = backup_c_max;
}

/**
 * @brief Recursively enumerates all possible interconnection trees between vertices in a molecular graph.
 *
 * This function generates all valid interconnection trees by linking vertices across components
 * while ensuring specific constraints are satisfied. It explores all combinations of edges
 * between components using recursion, applies conditions for valid connections, and backtracks
 * after exploring each possibility. Depending on the configured execution mode, a completed tree
 * is either pushed directly to `generatePaths` or stored for later sorting and processing.
 *
 * @param vertices_link_id_comp An array where each pair of elements represents a vertex ID and its corresponding
 * component ID. Format: `[vertex1, componentID1, vertex2, componentID2, ...]`.
 * @param tab_c An array representing the number of vertices in each component. Each index corresponds to a component
 * ID.
 * @param inter_tree An array for storing the current tree solution. Each pair of elements corresponds to two vertices
 *                  connected by an edge in the interconnection tree.
 * @param num_vertex The total number of vertices in the molecular graph.
 * @param k The total number of components in the graph.
 * @param K The current depth of the recursion (number of edges already added to the tree).
 * @param s The index of the current vertex being processed in the `vertices_link_id_comp` array.
 * @param marge The maximum number of vertices that can be temporarily excluded ("banned") from the tree.
 * @param cage Pointer to the `Cage_t` structure representing the molecular cage.
 * @param paths Pointer to the `Paths_t` structure for managing paths generated from the interconnection tree.
 * @param grid_sub Pointer to the `GridSubstrat` structure representing the molecular grid substrate.
 * @param substrat_t A pointer on table of positions of substrat's atoms.
 * @param options Structure containing configuration options for the tree generation process.
 * @param list_banned_edges An array to keep track of edges that are temporarily banned during the enumeration.
 * @param size_list_banned_edges Pointer to an integer that stores the size of the `list_banned_edges` array.
 *
 * @details
 * - The function operates recursively, adding edges between vertices from different components
 *   to form valid interconnection trees.
 * - Edges are added only if:
 *   1. The two vertices belong to different components.
 *   2. The sum of the sizes of the two components after the connection is valid.
 *   3. The distance between the two vertices does not exceed a predefined threshold.
 * - If all edges for a valid interconnection tree are found (`K == k - 1`), the behavior depends on
 *   the provided storage pointer: a non-null store collects the tree for deferred sorting/execution,
 *   while a null store triggers an immediate call to `generatePaths` with the current tree.
 * - Backtracking restores the state of the `vertices_link_id_comp` and `tab_c` arrays to explore other combinations.
 *
 * ### Recursion Logic:
 * - At each recursive step:
 *   1. Try adding an edge between the current vertex `s` and any other vertex `v`.
 *   2. If the edge satisfies the constraints, update the `vertices_link_id_comp` and `tab_c` arrays
 *      to reflect the new connection, and proceed with the next vertex.
 *   3. If no valid edge can be added, the function attempts to "ban" the current vertex and
 *      continues with the next vertex, decrementing the `marge` counter.
 * - After exploring all possibilities, the function backtracks by restoring the previous state
 *   of the `vertices_link_id_comp` and `tab_c` arrays.
 *
 * ### Usage Notes:
 * - The function assumes that the `vertices_link_id_comp` and `tab_c` arrays are pre-initialized with
 *   valid vertex-component mappings and component sizes, respectively.
 * - The `inter_tree` array must have sufficient size to store up to `(k - 1) * 2` elements.
 * - The `marge` parameter allows for flexibility in tree generation by permitting some vertices
 *   to be excluded temporarily. Setting `marge` to 0 disables this feature.
 *
 * ### Output:
 * - Each valid interconnection tree is stored and later processed by the `generatePaths` function,
 *   which generates paths based on the tree structure.
 * - The function also prints debug information, such as the tree number (`cpt_inter_tree_egv`).
 *
 * @see generatePaths
 * @see updateSubset
 * @see restoreSubset
 * @see updateVertex
 * @see restoreVertex
 */
void enumInterconnectionTrees(int *vertices_link_id_comp, int *tab_c, int *inter_tree, const int NUM_VERTEX, int k,
                              int K, int s, int marge, Cage_t *cage, Paths_t *paths, GridSubstrat *grid_sub,
                              double ***substrat_t, Options_t options, int *list_banned_edges,
                              int *size_list_banned_edges, InterconnectionTreeStore *tree_store) {
  // Check if we have a banned edge in the current interconnection tree
  if (options.isBannedEdges == 1 && *size_list_banned_edges > 0 && K > 0) {
    // Check if any edge in the inter_tree is banned
    for (int i = 0; i < (K - 1); i++) {
      for (int b = 0; b < *size_list_banned_edges;
           b++) { // TO DO verify if we need to check all banned edges, I think only the last one is enough
        if ((list_banned_edges[b * 2] == inter_tree[i * 2]) &&
            (list_banned_edges[b * 2 + 1] == inter_tree[i * 2 + 1])) {
          // Skip this edge if it is banned
          return;
        }
      }
    }
  }

  // Base case: if we have enough edges,
  if (K >= (k - 1)) {
    cpt_inter_tree_egv++;
    printf("New Interconnection Tree %d: \n", cpt_inter_tree_egv);
    if (tree_store) {
      appendTree(tree_store, inter_tree);
    } else {
      generatePaths(cage, inter_tree, paths, grid_sub, substrat_t, options, list_banned_edges,
                    size_list_banned_edges);
    }
    return;
  }
  int next_s;
  for (int v = s + 1; v < NUM_VERTEX; v++) {
    // Early exit optimizations
    if (UNLIKELY(vertices_link_id_comp[s * 2 + 1] == -1) || UNLIKELY(vertices_link_id_comp[v * 2 + 1] == -1) ||
        UNLIKELY(vertices_link_id_comp[s * 2 + 1] == vertices_link_id_comp[v * 2 + 1])) {
      continue;
    }

    // Pre-calculate component sum for efficiency
    int comp_sum = tab_c[vertices_link_id_comp[s * 2 + 1]] + tab_c[vertices_link_id_comp[v * 2 + 1]] - 2;

    double pair_distance = dist(coords(atom(cage, vertices_link_id_comp[s * 2])),
                  coords(atom(cage, vertices_link_id_comp[v * 2])));
    if (LIKELY((K == (k - 2)) || (comp_sum > 0)) &&
      LIKELY(!is_path_boundary_filter_enabled() ||
           pair_distance <= DIST_PATH_BOUNDARY(options.sizeMaxPath))) {
      // DIST_SIMPLE_PATTERN * options.sizeMaxPath + DIST_SIMPLE + DIST_ERROR)) { // Test if distance is valid
      // Add edge to the inter_tree
      inter_tree[K * 2] = vertices_link_id_comp[s * 2];
      inter_tree[K * 2 + 1] = vertices_link_id_comp[v * 2];

      // Backup and update
      int backup_s1, backup_s2, backup_max;
      int affected_indices[NUM_VERTEX]; // Array to store affected indices
      int affected_count = 0;           // Number of affected indices
      int backup_c_max, backup_c_min;

      updateSubset(tab_c, vertices_link_id_comp[s * 2 + 1], vertices_link_id_comp[v * 2 + 1], &backup_c_max,
                   &backup_c_min);
      backup_max = updateVertex(vertices_link_id_comp, s, v, NUM_VERTEX, &backup_s1, &backup_s2, affected_indices,
                                &affected_count);

      next_s = s;
      do {
        next_s = next_s + 1;
      } while (next_s < NUM_VERTEX && vertices_link_id_comp[next_s * 2 + 1] == -1);
      enumInterconnectionTrees(vertices_link_id_comp, tab_c, inter_tree, NUM_VERTEX, k, K + 1, next_s, marge, cage,
               paths, grid_sub, substrat_t, options, list_banned_edges, size_list_banned_edges,
               tree_store);
      // Restore changes
      restoreVertex(vertices_link_id_comp, s, v, backup_s1, backup_s2, affected_indices, affected_count, backup_max);
      restoreSubset(tab_c, vertices_link_id_comp[s * 2 + 1], vertices_link_id_comp[v * 2 + 1], backup_c_max,
                    backup_c_min);
    }
  }
  next_s = s;
  if ((marge > 0) && (tab_c[vertices_link_id_comp[s * 2 + 1]] > 0)) {
    int backup_s, backup_c;
    // Ban_s_Subset
    backup_c = tab_c[vertices_link_id_comp[s * 2 + 1]];
    tab_c[vertices_link_id_comp[s * 2 + 1]]--;
    // Ban_s_Vertex
    backup_s = vertices_link_id_comp[s * 2 + 1];
    vertices_link_id_comp[s * 2 + 1] = -1;
    do {
      next_s = next_s + 1;
    } while (next_s < NUM_VERTEX && vertices_link_id_comp[next_s * 2 + 1] == -1);
    enumInterconnectionTrees(vertices_link_id_comp, tab_c, inter_tree, NUM_VERTEX, k, K, next_s, marge - 1, cage, paths,
                 grid_sub, substrat_t, options, list_banned_edges, size_list_banned_edges, tree_store);
    // RestoreBan_s_Vertex
    vertices_link_id_comp[s * 2 + 1] = backup_s;
    // RestoreBan_s_Subset
    tab_c[vertices_link_id_comp[s * 2 + 1]] = backup_c;
  }
}

/**
 * @brief Finds and generates interconnection trees among the linkable vertices in a molecular cage structure.
 *
 * This function identifies linkable vertices in a cage structure, determines their connected components,
 * and generates interconnection trees using the `enumInterconnectionTrees` function.
 * The trees connect vertices across components, representing all possible valid interconnection patterns.
 *
 * @param cage A pointer to the `Cage_t` structure representing the molecular system.
 * @param substrat_t A pointer on table of positions of substrat's atoms.
 * @param grid_sub A pointer to the `GridSubstrat` structure representing the grid substrate of the molecular system.
 * @param options An `Options_t` structure containing configuration parameters for the computation.
 *
 * @details
 * ### Function Workflow:
 * 1. **Initialization and Component Detection:**
 *    - The function begins by determining the number of vertices in the cage and their connected components.
 *    - It initializes an array to map all vertices to their respective components using the `connectedComponents`
 * function.
 *
 * 2. **Filtering Linkable Vertices:**
 *    - It filters out the vertices flagged as `LINKABLE_F` (vertices eligible for interconnection) and maps them to
 * their components.
 *    - This mapping is stored in the `vertices_link_id_comp` array, where:
 *      - `vertices_link_id_comp[2 * i]` stores the vertex ID.
 *      - `vertices_link_id_comp[2 * i + 1]` stores the component ID of the vertex.
 *    - Additionally, it counts the number of vertices in each component using the `tab_c` array.
 *
 * 3. **Interconnection Tree Generation:**
 *    - The function initializes an array `inter_tree` to store the edges of the interconnection tree and, when
 *      requested, a dynamic store to retain every valid tree.
 *    - A `Paths_t` structure is created to manage paths generated by the interconnection trees.
 *    - The `enumInterconnectionTrees` function is called to recursively explore all valid interconnection trees.
 *      Depending on user configuration, trees are either processed on-the-fly or stored, sorted by their total
 *      length, and then passed to `generatePaths`.
 *
 * 4. **Debug Information:**
 *    - Debug output is printed, including:
 *      - The total number of components.
 *      - The mapping of linkable vertices to their components (`vertices_link_id_comp`).
 *      - The count of vertices in each component (`tab_c`).
 *
 * 5. **Finalization:**
 *    - After all trees are generated, memory allocated for arrays and data structures is freed.
 *
 * ### Memory Management:
 * - Dynamically allocated memory for `all_vertices`, `vertices_link_id_comp`, `tab_c`, `inter_tree`, and `paths` is
 * freed before the function exits to prevent memory leaks.
 *
 * ### Usage Notes:
 * - Ensure the `cage` structure is correctly initialized with vertex and component information before calling this
 * function.
 * - The `options` parameter must be properly configured, particularly the `sizeMaxPath` field.
 * - The `grid_sub` parameter provides additional structural context for path generation within the molecular system.
 *
 * @see connectedComponents
 * @see enumInterconnectionTrees
 */
void findInterconnection(Cage_t *cage, GridSubstrat *grid_sub, double ***substrat_t, Options_t options) {
  // Example usage:
  int num_vertex = size(cage);
  int num_components;

  int *all_vertices = (int *)malloc(num_vertex * sizeof(int)); // Array for storing vertex to component mapping
  for (int i = 0; i < num_vertex; i++) {
    all_vertices[i] = -1;
  }
  connectedComponents(cage, all_vertices, &num_components);


  int num_vertex_linkable = 0;
  for (int i = 0; i < num_vertex; i++) {
    if (flag(atom(cage, i)) == LINKABLE_F)
      num_vertex_linkable++; // To Do add here a test to know if this vertex can be use (no collision) for start or
                             // end a path
  }
  int *vertices_link_id_comp =
      (int *)malloc(2 * num_vertex_linkable *
                    sizeof(int)); // Array for storing linkable vertex to component mapping 2D first ID, second Part
  int *tab_c = (int *)calloc(num_components, sizeof(int)); // Array for storing the count of vertices per component
  int pos_link = 0;
  int part;
  for (int i = 0; i < num_vertex; i++) {
    if (flag(atom(cage, i)) == LINKABLE_F) {
      part = all_vertices[i];
      vertices_link_id_comp[pos_link * 2] = i;
      vertices_link_id_comp[pos_link * 2 + 1] = part;
      tab_c[part]++;
      pos_link++;
    }
  }

  int size_list_banned_edges = 0;
  int *list_banned_edges = (int *)malloc((((num_vertex_linkable * (num_vertex_linkable - 1)) / 2) * 2) *
                                         sizeof(int)); // Array for storing banned edges, size is maximum possible edges

  int min_vertices_component = INT_MAX;
  int max_vertices_component = 0;
  for (int i = 0; i < num_components; i++) {
    int count = tab_c[i];
    if (count < min_vertices_component) {
      min_vertices_component = count;
    }
    if (count > max_vertices_component) {
      max_vertices_component = count;
    }
  }
  if (min_vertices_component == INT_MAX) {
    min_vertices_component = 0;
  }

  printf("COMPONENT_SUMMARY components=%d min_vertices=%d max_vertices=%d\n", num_components,
         min_vertices_component, max_vertices_component);
  fflush(stdout);

  free(all_vertices);
  // printf("num comp %d \n", num_components);

  // printf("tab_s array (vertex to component mapping):\n");
  // for (int i = 0; i < num_vertex_linkable; i++) {
  //   printf("Vertex %d: ID %d: Component %d\n", i, vertices_link_id_comp[2 * i], vertices_link_id_comp[2 * i + 1]);
  // }
  // printf("\nC array (component vertex count):\n");
  // for (int i = 0; i < num_components; i++) {
  //   printf("Component %d: %d vertices\n", i, tab_c[i]);
  // }

  int edge_slot_count = 2 * (num_components - 1);
  int *inter_tree = malloc(edge_slot_count * sizeof(int)); // Interconnection tree working buffer
  Paths_t *paths = pthCreate(options.sizeMaxPath, num_components);
  if (get_current_distance_type() != DISTANCE_EUCLIDEAN) {
    createGrid(paths->grids[0], cage, paths, substrat_t, grid_sub);
    initMinHeap(paths->minHeaps[0], paths->grids[0]->depth * paths->grids[0]->width * paths->grids[0]->height);
    // writeGridToMol2(paths->grids[0], "grid.mol2", 0);
    // writeGridToMol2(paths->grids[0], "gridFull.mol2", 1);
  }

  InterconnectionTreeStore tree_store;
  InterconnectionTreeStore *store_ptr = NULL;
  if (options.sortInterTreesBeforePaths) {
    initTreeStore(&tree_store, edge_slot_count);
    store_ptr = &tree_store;
  }

  enumInterconnectionTrees(vertices_link_id_comp, tab_c, inter_tree, num_vertex_linkable, num_components, 0, 0,
                           num_vertex_linkable - (2 * (num_components - 1)), cage, paths, grid_sub, substrat_t,
                           options, list_banned_edges, &size_list_banned_edges, store_ptr);

  if (store_ptr) {
    int num_paths = (edge_slot_count > 0) ? edge_slot_count / 2 : 0;
    for (size_t i = 0; i < tree_store.count; i++) {
      if (num_paths > 0) {
        pthInit(paths, tree_store.items[i].edges, cage);
      }
      tree_store.items[i].total_size =
          computeTreeTotalSize(&tree_store.items[i], num_paths, cage, paths, grid_sub, substrat_t);
    }

    if (tree_store.count > 1) {
      qsort(tree_store.items, tree_store.count, sizeof(StoredInterconnectionTree), compareStoredTrees);
    }

    if (num_paths > 0) {
      for (size_t i = 0; i < tree_store.count; i++) {
        if (options.isBannedEdges == 1 && size_list_banned_edges > 0 &&
            treeContainsBannedEdge(tree_store.items[i].edges, num_paths, list_banned_edges, size_list_banned_edges)) {
          continue;
        }
        generatePaths(cage, tree_store.items[i].edges, paths, grid_sub, substrat_t, options, list_banned_edges,
                      &size_list_banned_edges);
      }
    }
    freeTreeStore(&tree_store);
  }

  // Free allocated memory
  free(list_banned_edges);
  free(vertices_link_id_comp);
  free(tab_c);
  free(inter_tree);
  pthDelete(paths);
}
