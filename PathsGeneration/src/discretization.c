#include "discretization.h"
#include "distance.h"

#include <float.h>
#include <math.h>

#define MAX(a, b) (((a) > (b)) ? (a) : (b))

void initGrid(Grid_t *grid, int size_x, int size_y, int size_z) {
  // Free existing grid->nodes if already allocated
  if (grid->nodes) {
    for (int x = 0; x < grid->width; x++) {
      for (int y = 0; y < grid->height; y++) {
        for (int z = 0; z < grid->depth; z++) {
          free(grid->nodes[x][y][z].indexStartsCandidates);
        }
        free(grid->nodes[x][y]);
      }
      free(grid->nodes[x]);
    }
    free(grid->nodes);
    grid->nodes = NULL; // Set to NULL to avoid dangling pointers
  }

  // Update grid dimensions
  grid->width = 2 * size_x + 1;
  grid->height = 2 * size_y + 1;
  grid->depth = 2 * size_z + 1;

  grid->offset_x = size_x;
  grid->offset_y = size_y;
  grid->offset_z = size_z;

  // Allocate 3D array
  grid->nodes = malloc(grid->width * sizeof(Node **));
  for (int x = 0; x < grid->width; x++) {
    grid->nodes[x] = malloc(grid->height * sizeof(Node *));
    for (int y = 0; y < grid->height; y++) {
      grid->nodes[x][y] = malloc(grid->depth * sizeof(Node));
      for (int z = 0; z < grid->depth; z++) {
        grid->nodes[x][y][z].x = x; // - grid->offset_x;
        grid->nodes[x][y][z].y = y; // - grid->offset_y;
        grid->nodes[x][y][z].z = z; // - grid->offset_z;
        grid->nodes[x][y][z].walkable = 1;
        grid->nodes[x][y][z].closed = 0; // 1 = closed, 0 = open for A* algorithm
        grid->nodes[x][y][z].opened = 0; // 1 = opened, 0 = closed for A* algorithm
        grid->nodes[x][y][z].parent = NULL;
        grid->nodes[x][y][z].heapIndex = -1;
        grid->nodes[x][y][z].nbTimesCandidate = 0;
        grid->nodes[x][y][z].indexStartsCandidates = (int *)malloc(MAX_POSITION_KEEP_360 * sizeof(int));
        if (grid->nodes[x][y][z].indexStartsCandidates == NULL) {
          fprintf(stderr, "Memory allocation failed for indexStartsCandidates\n");
          exit(EXIT_FAILURE);
        }
      }
    }
  }

  // Initialize visited nodes tracking
  if (grid->visitedNodes) {
    free(grid->visitedNodes);
    grid->visitedNodes = NULL;
  }
  initVisitedNodes(grid, 1024); // Start with capacity for 1024 visited nodes
}

/**
 * @brief reset values used for A* algorithm.
 *
 * This function resets the values used for the A* algorithm in the grid nodes.
 * It sets the gCost to DBL_MAX, fCost to 0, hCost to 0, parent to NULL,
 * closed to 0, and opened to 0 for all nodes in the grid.
 *
 * @param grid Pointer to the Grid_t structure representing the grid.
 */
void resetGridState(Grid_t *grid) {

  // print the grid before reset
  // printf("Resetting grid state...\n");
  // for (int x = 0; x < grid->width; x++) {
  //   for (int y = 0; y < grid->height; y++) {
  //     for (int z = 0; z < grid->depth; z++)
  //       if (grid->nodes[x][y][z].closed == 1 || grid->nodes[x][y][z].opened == 1)
  //         printf("Node (%d, %d, %d): gCost = %f, fCost = %f, hCost = %f, closed = %d, opened = %d\n", x, y, z,
  //                grid->nodes[x][y][z].gCost, grid->nodes[x][y][z].fCost, grid->nodes[x][y][z].hCost,
  //                grid->nodes[x][y][z].closed, grid->nodes[x][y][z].opened);
  //   }
  // }

  // Reset only the nodes that were actually visited during the last pathfinding
  for (int i = 0; i < grid->numVisited; i++) {
    Node *node = grid->visitedNodes[i];
    node->gCost = DBL_MAX;
    node->fCost = 0;
    node->hCost = 0;
    node->parent = NULL;
    node->closed = 0;
    node->opened = 0;
    node->nbTimesCandidate = 0;
    // printf("Resetting node (%d, %d, %d): gCost = %f, fCost = %f, hCost = %f, closed = %d, opened = %d\n", node->x,
    //        node->y, node->z, node->gCost, node->fCost, node->hCost, node->closed, node->opened);
  }
  // Reset the visited count for the next pathfinding
  grid->numVisited = 0;
}

/**
 * @brief Find the minimum and maximum coordinates for a cage, paths, and substrate.
 *
 * This function calculates the minimum and maximum coordinates in the x, y, and z dimensions
 * for a given cage, paths, and substrate. It updates the provided pointers with the calculated
 * values.
 *
 * @param cage Pointer to the Cage_t structure representing the cage.
 * @param paths Pointer to the Paths_t structure representing the paths.
 * @param substrat_t Pointer to a 3D array representing the substrate.
 * @param substrat_size Size of the substrate array.
 * @param min_x Pointer to store the minimum x-coordinate.
 * @param max_x Pointer to store the maximum x-coordinate.
 * @param min_y Pointer to store the minimum y-coordinate.
 * @param max_y Pointer to store the maximum y-coordinate.
 * @param min_z Pointer to store the minimum z-coordinate.
 * @param max_z Pointer to store the maximum z-coordinate.
 */
void findBounds(Cage_t *cage, Paths_t *paths, GridSubstrat *grid_sub, double *min_x, double *max_x, double *min_y,
                double *max_y, double *min_z, double *max_z) {
  if (grid_sub->substratSize != 0) {
    *min_x = grid_sub->xMinWOGap;
    *max_x = grid_sub->xMaxWOGap;
    *min_y = grid_sub->yMinWOGap;
    *max_y = grid_sub->yMaxWOGap;
    *min_z = grid_sub->zMinWOGap;
    *max_z = grid_sub->zMaxWOGap;

    if (cage->xMinWOGap < *min_x)
      *min_x = cage->xMinWOGap;
    if (cage->xMaxWOGap > *max_x)
      *max_x = cage->xMaxWOGap;
    if (cage->yMinWOGap < *min_y)
      *min_y = cage->yMinWOGap;
    if (cage->yMaxWOGap > *max_y)
      *max_y = cage->yMaxWOGap;
    if (cage->zMinWOGap < *min_z)
      *min_z = cage->zMinWOGap;
    if (cage->zMaxWOGap > *max_z)
      *max_z = cage->zMaxWOGap;
  } else {
    *min_x = cage->xMinWOGap;
    *max_x = cage->xMaxWOGap;
    *min_y = cage->yMinWOGap;
    *max_y = cage->yMaxWOGap;
    *min_z = cage->zMinWOGap;
    *max_z = cage->zMaxWOGap;
  }

  // cage to do improve don't calculate every time the max and min for cage but only one time
  for (int i = 0; i < size(cage); i++) {
    if (coords(atom(cage, i)).x < *min_x)
      *min_x = coords(atom(cage, i)).x;
    if (coords(atom(cage, i)).x > *max_x)
      *max_x = coords(atom(cage, i)).x;
    if (coords(atom(cage, i)).y < *min_y)
      *min_y = coords(atom(cage, i)).y;
    if (coords(atom(cage, i)).y > *max_y)
      *max_y = coords(atom(cage, i)).y;
    if (coords(atom(cage, i)).z < *min_z)
      *min_z = coords(atom(cage, i)).z;
    if (coords(atom(cage, i)).z > *max_z)
      *max_z = coords(atom(cage, i)).z;
  }
  // paths
  int size;
  for (int k = 0; k < paths->currentPath; k++) {
    size = paths->curPthPos[k];
    Point_t p_tmp;
    for (int i = 2; i <= size; i++) {
      // int size = (path->patternNum[i] == CYCLE_PATTERN) ? MAX_NB_ATOMS_PATTERN : 3;
      for (int j = 0; j < MAX_NB_ATOMS_PATTERN; j++) {
        p_tmp = paths->patterns[indexPointPaths(k, i, paths->positionCurNum[(indexPathPosition(k, i, paths->sizeMax))],
                                                j, paths->sizeMax)];
        if (p_tmp.x < *min_x)
          *min_x = p_tmp.x;
        if (p_tmp.x > *max_x)
          *max_x = p_tmp.x;
        if (p_tmp.y < *min_y)
          *min_y = p_tmp.y;
        if (p_tmp.y > *max_y)
          *max_y = p_tmp.y;
        if (p_tmp.z < *min_z)
          *min_z = p_tmp.z;
        if (p_tmp.z > *max_z)
          *max_z = p_tmp.z;
      }
    }
    for (int i = size + 1; i <= size + 2; i++) {
      p_tmp = paths->patterns[indexPointPaths(k, i, paths->positionCurNum[(indexPathPosition(k, i, paths->sizeMax))], 1,
                                              paths->sizeMax)];
      if (p_tmp.x < *min_x)
        *min_x = p_tmp.x;
      if (p_tmp.x > *max_x)
        *max_x = p_tmp.x;
      if (p_tmp.y < *min_y)
        *min_y = p_tmp.y;
      if (p_tmp.y > *max_y)
        *max_y = p_tmp.y;
      if (p_tmp.z < *min_z)
        *min_z = p_tmp.z;
      if (p_tmp.z > *max_z)
        *max_z = p_tmp.z;
      p_tmp = paths->patterns[indexPointPaths(k, i, paths->positionCurNum[(indexPathPosition(k, i, paths->sizeMax))], 2,
                                              paths->sizeMax)];
      if (p_tmp.x < *min_x)
        *min_x = p_tmp.x;
      if (p_tmp.x > *max_x)
        *max_x = p_tmp.x;
      if (p_tmp.y < *min_y)
        *min_y = p_tmp.y;
      if (p_tmp.y > *max_y)
        *max_y = p_tmp.y;
      if (p_tmp.z < *min_z)
        *min_z = p_tmp.z;
      if (p_tmp.z > *max_z)
        *max_z = p_tmp.z;
    }
  }
  // printf("min_x %f max_x %f min_y %f max_y %f min_z %f max_z %f\n",*min_x,*max_x,*min_y,*max_y,*min_z,*max_z);
}

// Helper function to mark a sphere of voxels as non-walkable
void markSphereWorld(Grid_t *grid, double cx_w, double cy_w, double cz_w, double radius) {
  // Convert world center to voxel indices
  int cx = WORLD_TO_VOXEL(cx_w, grid->offset_x, STEP_GRID_VOXEL);
  int cy = WORLD_TO_VOXEL(cy_w, grid->offset_y, STEP_GRID_VOXEL);
  int cz = WORLD_TO_VOXEL(cz_w, grid->offset_z, STEP_GRID_VOXEL);

  // Radius in voxel units
  int r_voxels = (int)ceil(radius / STEP_GRID_VOXEL);

  for (int x = cx - r_voxels; x <= cx + r_voxels; x++) {
    for (int y = cy - r_voxels; y <= cy + r_voxels; y++) {
      for (int z = cz - r_voxels; z <= cz + r_voxels; z++) {
        // Bounds check
        if (x >= 0 && x < grid->width && y >= 0 && y < grid->height && z >= 0 && z < grid->depth) {
          // Convert voxel index back to world-space position
          double dx = VOXEL_TO_WORLD(x, grid->offset_x, STEP_GRID_VOXEL) - cx_w;
          double dy = VOXEL_TO_WORLD(y, grid->offset_y, STEP_GRID_VOXEL) - cy_w;
          double dz = VOXEL_TO_WORLD(z, grid->offset_z, STEP_GRID_VOXEL) - cz_w;

          if ((dx * dx + dy * dy + dz * dz) <= radius * radius) {
            grid->nodes[x][y][z].walkable = 0;
          }
        }
      }
    }
  }
}

/**
 * @brief Updates the grid to mark non-walkable voxels based on the substrate, cage, and paths.
 *
 * This function processes the substrate, cage, and paths to determine which voxels in the grid
 * should be marked as non-walkable. It considers the spatial arrangement and the radius of atoms
 * to ensure accurate marking of the grid.
 *
 * @param grid Pointer to the Grid_t structure that represents the voxel grid.
 * @param cage Pointer to the Cage_t structure containing information about the molecular cage.
 * @param paths Pointer to the Paths_t structure containing the paths to be considered.
 * @param substrat_t Pointer to a 3D array representing the substrate's spatial data.
 * @param substrat_size The size of the substrate array in each dimension.
 */
void markVoxels(Grid_t *grid, Cage_t *cage, Paths_t *paths, double ***substrat_t, int substrat_size) {

  // Mark substrate voxels
  for (int i = 0; i < substrat_size; i++) {
    markSphereWorld(grid, (*substrat_t)[i][0], (*substrat_t)[i][1], (*substrat_t)[i][2], DIST_GAP_SUBSTRATE);
  }

  // Mark cage voxels
  for (int i = 0; i < size(cage); i++) {
    AtomCage_t *a = atom(cage, i);
    if (flag(a) != LINKABLE_F) {
      markSphereWorld(grid, coords(a).x, coords(a).y, coords(a).z, DIST_GAP_CAGE);
    }
  }

  // Mark path voxels
  for (int k = 0; k < paths->currentPath; k++) {
    int size = paths->curPthPos[k];

    for (int i = 2; i <= size; i++) {
      for (int j = 0; j < MAX_NB_ATOMS_PATTERN; j++) {
        Point_t p = paths->patterns[indexPointPaths(
            k, i, paths->positionCurNum[indexPathPosition(k, i, paths->sizeMax)], j, paths->sizeMax)];
        markSphereWorld(grid, p.x, p.y, p.z, DIST_GAP_CAGE);
      }
    }

    for (int i = size + 1; i <= size + 2; i++) {
      for (int j = 1; j <= 2; j++) {
        Point_t p = paths->patterns[indexPointPaths(
            k, i, paths->positionCurNum[indexPathPosition(k, i, paths->sizeMax)], j, paths->sizeMax)];
        markSphereWorld(grid, p.x, p.y, p.z, DIST_GAP_CAGE);
      }
    }
  }
}

/**
 * @brief Creates a grid based on the cage, paths, and substrate.
 *
 * This function initializes a grid structure and marks the voxels as walkable or non-walkable
 * based on the provided cage, paths, and substrate data. It calculates the grid size based on
 * the bounding box of the cage and paths.
 *
 * @param grid Pointer to the Grid_t structure representing the grid.
 * @param cage Pointer to the Cage_t structure representing the cage.
 * @param paths Pointer to the Paths_t structure representing the paths.
 * @param substrat_t Pointer to a 3D array representing the substrate's spatial data.
 * @param substrat_size The size of the substrate array in each dimension.
 */
void createGrid(Grid_t *grid, Cage_t *cage, Paths_t *paths, double ***substrat_t, GridSubstrat *grid_sub) {
  double min_x, max_x, min_y, max_y, min_z, max_z;
  findBounds(cage, paths, grid_sub, &min_x, &max_x, &min_y, &max_y, &min_z, &max_z);

  // Calculate the total extent of the grid (min-max coordinates)
  double grid_min_x = min_x - DIST_GAP_MAX;
  double grid_max_x = max_x + DIST_GAP_MAX;
  double grid_min_y = min_y - DIST_GAP_MAX;
  double grid_max_y = max_y + DIST_GAP_MAX;
  double grid_min_z = min_z - DIST_GAP_MAX;
  double grid_max_z = max_z + DIST_GAP_MAX;

  // Calculate grid sizes (half-extents)
  int size_x = (int)ceil(MAX(fabs(grid_max_x), fabs(grid_min_x)) / STEP_GRID_VOXEL);
  int size_y = (int)ceil(MAX(fabs(grid_max_y), fabs(grid_min_y)) / STEP_GRID_VOXEL);
  int size_z = (int)ceil(MAX(fabs(grid_max_z), fabs(grid_min_z)) / STEP_GRID_VOXEL);

  // Create the grid with these calculated sizes
  initGrid(grid, size_x, size_y, size_z);
  markVoxels(grid, cage, paths, substrat_t, grid_sub->substratSize);
}

/**
 * @brief Prints the details of a node in the grid.
 *
 * This function prints the coordinates, real-world coordinates, walkable status,
 * gCost, fCost, hCost, closed status, opened status, parent pointer, and heap index
 * of a given node.
 *
 * @param node Pointer to the Node structure to be printed.
 * @param grid Pointer to the Grid_t structure for coordinate conversion.
 */
void printNode(Node *node, Grid_t *grid) {
  if (node == NULL) {
    printf("Node is NULL\n");
    return;
  }
  printf("Node: (%d, %d, %d), realWorld (%lf,%lf,%lf) Walkable: %d\n", node->x, node->y, node->z,
         VOXEL_TO_WORLD(node->x, grid->offset_x, STEP_GRID_VOXEL),
         VOXEL_TO_WORLD(node->y, grid->offset_y, STEP_GRID_VOXEL),
         VOXEL_TO_WORLD(node->z, grid->offset_z, STEP_GRID_VOXEL), node->walkable);
  printf("gCost: %f, fCost: %f, hCost: %f\n", node->gCost, node->fCost, node->hCost);
  printf("closed: %d, opened: %d\n", node->closed, node->opened);
  printf("parent: %p, heapIndex: %d\n", node->parent, node->heapIndex);
}

/**
 * @brief Writes the grid to a .mol2 file, either walkable or non-walkable nodes.
 *
 * This function iterates through the grid and writes the world coordinates of the
 * selected nodes (walkable or non-walkable) to a .mol2 file.
 *
 * @param grid Pointer to the Grid_t structure.
 * @param filename Name of the .moc2 file to write to.
 * @param writeWalkable If 1, writes walkable nodes; if 0, writes non-walkable nodes.
 */
void writeGridToMol2(Grid_t *grid, const char *filename, int writeWalkable) {
  FILE *file = fopen(filename, "w");
  if (!file) {
    perror("Error opening file");
    return;
  }

  // Count matching nodes
  int atom_count = 0;
  for (int x = 0; x < grid->width; x++) {
    for (int y = 0; y < grid->height; y++) {
      for (int z = 0; z < grid->depth; z++) {
        if (grid->nodes[x][y][z].walkable == writeWalkable) {
          atom_count++;
        }
      }
    }
  }

  // Write MOL2 header
  fprintf(file, "@<TRIPOS>MOLECULE\n");
  fprintf(file, "grid_export\n");
  fprintf(file, "%d 0 0 0 0\n", atom_count); // atoms, bonds, etc.
  fprintf(file, "SMALL\n");
  fprintf(file, "NO_CHARGES\n\n");

  // Atom block
  fprintf(file, "@<TRIPOS>ATOM\n");
  int id = 1;
  for (int x = 0; x < grid->width; x++) {
    for (int y = 0; y < grid->height; y++) {
      for (int z = 0; z < grid->depth; z++) {
        if (grid->nodes[x][y][z].walkable == writeWalkable) {
          double wx = VOXEL_TO_WORLD(x, grid->offset_x, STEP_GRID_VOXEL);
          double wy = VOXEL_TO_WORLD(y, grid->offset_y, STEP_GRID_VOXEL);
          double wz = VOXEL_TO_WORLD(z, grid->offset_z, STEP_GRID_VOXEL);

          fprintf(file, "%7d G%d       %10.4f %10.4f %10.4f Gr      1 <0>        0.0000\n", id, id, wx, wy, wz);
          id++;
        }
      }
    }
  }

  fclose(file);
}

/**
 * @brief Frees the memory allocated for a Grid_t structure.
 *
 * This function deallocates the memory used by the Grid_t structure, including its nodes.
 *
 * @param grid Pointer to the Grid_t structure to be freed.
 */
void freeGrid(Grid_t *grid) {
  if (grid == NULL) {
    return; // Nothing to free
  }

  // Free visited nodes tracking array
  if (grid->visitedNodes != NULL) {
    free(grid->visitedNodes);
    grid->visitedNodes = NULL;
  }

  // Free the 3D nodes array if it exists
  if (grid->nodes != NULL) {
    for (int x = 0; x < grid->width; x++) {
      if (grid->nodes[x] != NULL) {
        for (int y = 0; y < grid->height; y++) {
          if (grid->nodes[x][y] != NULL) {
            for (int z = 0; z < grid->depth; z++) {
              if (grid->nodes[x][y][z].indexStartsCandidates != NULL) {
                free(grid->nodes[x][y][z].indexStartsCandidates);
                grid->nodes[x][y][z].indexStartsCandidates = NULL;
              }
            }
            free(grid->nodes[x][y]);
            grid->nodes[x][y] = NULL;
          }
        }
        free(grid->nodes[x]);
        grid->nodes[x] = NULL;
      }
    }
    free(grid->nodes);
    grid->nodes = NULL;
  }

  free(grid);
}

/**
 * @brief Add a node to the visited nodes list for efficient reset.
 *
 * This function adds a node to the visitedNodes array for tracking during pathfinding.
 * It will automatically resize the array if needed.
 *
 * @param grid Pointer to the Grid_t structure.
 * @param node Pointer to the node to mark as visited.
 */
void addVisitedNode(Grid_t *grid, Node *node) {
  // Check if we need to resize the visitedNodes array
  if (grid->numVisited >= grid->maxVisited) {
    // Double the capacity if needed
    int new_max_visited = grid->maxVisited == 0 ? 1024 : grid->maxVisited * 2;
    Node **new_visited_nodes = realloc(grid->visitedNodes, new_max_visited * sizeof(Node *));
    if (new_visited_nodes == NULL) {
      fprintf(stderr, "Memory allocation failed for visitedNodes expansion\n");
      exit(EXIT_FAILURE);
    }
    grid->visitedNodes = new_visited_nodes;
    grid->maxVisited = new_max_visited;
  }

  // Add the node to the visited list
  grid->visitedNodes[grid->numVisited] = node;
  grid->numVisited++;
}

/**
 * @brief Initialize the visited nodes tracking for a grid.
 *
 * This function initializes the visitedNodes array with an initial capacity.
 *
 * @param grid Pointer to the Grid_t structure.
 * @param initialCapacity Initial capacity for the visitedNodes array.
 */
void initVisitedNodes(Grid_t *grid, int initialCapacity) {
  grid->visitedNodes = malloc(initialCapacity * sizeof(Node *));
  if (grid->visitedNodes == NULL) {
    fprintf(stderr, "Memory allocation failed for visitedNodes\n");
    exit(EXIT_FAILURE);
  }
  grid->maxVisited = initialCapacity;
  grid->numVisited = 0;
}