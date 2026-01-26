#include "distance.h"
#include "discretization.h"
#include "structure.h"

#include <ctype.h>
#include "float.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef ENABLE_STATS
extern int nb_path_boundary_allowed;
extern int nb_path_boundary_blocked;
#endif

#define SQRT2 1.4142135623730951
#define SQRT3 1.7320508075688772

// Default: Euclidean distance
DistanceFunc distance_func = dist;
DistanceType current_distance_type = DISTANCE_EUCLIDEAN;
static int path_boundary_filter_enabled = 1;

/**
 * @brief Squared Euclidean Distance between two points.
 *
 * @param point1 Double table of the 3 coords.
 * @param point2 Double table of the 3 coords.
 */
double squaredEuclideanDistance(double point1[3], double point2[3]) {
  return (point1[0] - point2[0]) * (point1[0] - point2[0]) + (point1[1] - point2[1]) * (point1[1] - point2[1]) +
         (point1[2] - point2[2]) * (point1[2] - point2[2]);
}

/**
 * @brief Squared Euclidean Distance between three double for point 1 and a table of double for point 2.
 *
 * @param p1x Double x coord of point 1.
 * @param p1y Double y coord of point 1.
 * @param p1z Double z coord of point 1.
 * @param point2 Double table of the 3 coords.
 */
double squaredEuclideanDistanceCoordsPoint(double p1x, double p1y, double p1z, double point2[3]) {
  return (p1x - point2[0]) * (p1x - point2[0]) + (p1y - point2[1]) * (p1y - point2[1]) +
         (p1z - point2[2]) * (p1z - point2[2]);
}

/**
 * @brief Squared Euclidean Distance between three double for point 1 and a table of double for point 2.
 *
 * @param p1 Point_t structure
 * @param point2 Double table of the 3 coords.
 */
double squaredEuclideanDistancefloatPointT(Point_t p1, double point2[3]) {
  return (p1.x - point2[0]) * (p1.x - point2[0]) + (p1.y - point2[1]) * (p1.y - point2[1]) +
         (p1.z - point2[2]) * (p1.z - point2[2]);
}

/**
 * @brief Compute the Square Euclidean distance between two points.
 *
 * This function calculates the Square Euclidean distance between two 3D points A and B.
 *
 * @param A The first point.
 * @param B The second point.
 * @return The Euclidean distance between A and B.
 */
double squaredEuclideanDistancePointTPointT(Point_t A, Point_t B) {
  return (A.x - B.x) * (A.x - B.x) + (A.y - B.y) * (A.y - B.y) + (A.z - B.z) * (A.z - B.z); // euclidian
}

/**
 * @brief Compute the Euclidean distance between two points.
 *
 * This function calculates the Euclidean distance between two 3D points A and B.
 *
 * @param A The first point.
 * @param B The second point.
 * @return The Euclidean distance between A and B.
 */
double dist(Point_t A, Point_t B) {
  double dx = A.x - B.x;
  double dy = A.y - B.y;
  double dz = A.z - B.z;
  return sqrt(dx * dx + dy * dy + dz * dz); // euclidian
}

/**
 * @brief Checks if the squared distance between two 3D points is less than a specified threshold.
 *
 * This function calculates the squared Euclidean distance between two points `A` and `B` in 3D space
 * and checks if the squared distance is less than a specified threshold `dist`. The comparison is made
 * without taking the square root, which can be computationally expensive, and instead uses the squared
 * distance directly.
 *
 * @param A The first `Point_t` structure (3D point).
 * @param B The second `Point_t` structure (3D point).
 * @param dist The threshold distance. The squared distance between `A` and `B` is compared with `dist^2`.
 *
 * @return `1` (true) if the squared distance between `A` and `B` is less than `dist^2`, `0` (false) otherwise.
 *
 * @note This function performs a comparison based on squared distances and avoids computing the square root, which
 * is often used in distance comparisons. The threshold `dist` should be considered as the *maximum distance* allowed
 * between the two points.
 *
 * @see Point_t
 */
int distInf(Point_t A, Point_t B, double dist) {
  return (A.x - B.x) * (A.x - B.x) + (A.y - B.y) * (A.y - B.y) + (A.z - B.z) * (A.z - B.z) < dist * dist;
}

double inline heuristic(Node *a, Node *b) {
  // Manhattan distance heuristic
  // return (fabs(a->x - b->x) + fabs(a->y - b->y) + fabs(a->z - b->z))*STEP_GRID_VOXEL;

  // Euclidean distance heuristic - optimized to avoid redundant calculations
  double dx = a->x - b->x;
  double dy = a->y - b->y;
  double dz = a->z - b->z;
  return sqrt(dx * dx + dy * dy + dz * dz) * STEP_GRID_VOXEL;

  /*The voxel distance is composed of distances ∆x, ∆y, and ∆z in each respective axis.
   * It is calculated as: h(n, n′) = (√3−√2) * dmin + (√2−1) * dmid + dmax
   * where dmax = max(∆x, ∆y, ∆z), dmin = min(∆x, ∆y, ∆z), and dmid = {∆x, ∆y, ∆z} \ {dmax, dmin}.
   * src: The Jump Point Search Pathfinding System in 3D.
   */
  // int dx = abs(a->x - b->x);
  // int dy = abs(a->y - b->y);
  // int dz = abs(a->z - b->z);

  // int dmin = dx < dy ? (dx < dz ? dx : dz) : (dy < dz ? dy : dz);
  // int dmax = dx > dy ? (dx > dz ? dx : dz) : (dy > dz ? dy : dz);
  // int dmid = dx + dy + dz - dmin - dmax;

  // return ((SQRT3 - SQRT2) * dmin + (SQRT2 - 1) * dmid + dmax) * STEP_GRID_VOXEL;
}

double moveCost(Node *a, Node *b) {
  // printf("moveCost: a (%d, %d, %d) b (%d, %d, %d)\n", a->x, a->y, a->z, b->x, b->y, b->z);
  int dx = abs(a->x - b->x);
  int dy = abs(a->y - b->y);
  int dz = abs(a->z - b->z);

  int axes_moved = dx + dy + dz;
  double move_cost = DBL_MAX;
  if (axes_moved == 1)
    move_cost = 1.0;
  else if (axes_moved == 2)
    move_cost = SQRT2;
  else if (axes_moved == 3)
    move_cost = SQRT3;
  else if (axes_moved == 0)
    move_cost = 0.0;
  else
    fprintf(stderr, "Error: moveCost: axes_moved = %d\n", axes_moved);

  return move_cost * STEP_GRID_VOXEL;
}

/**
 * @brief Reconstruct the path from the end node to the start node.
 *
 * This function reconstructs the path from the end node to the start node by following
 * the parent pointers of each node. It prints the coordinates of each voxel in the path.
 *
 * @param end Pointer to the end node.
 */
void reconstructPath(Node *end, Grid_t *grid) {
  Node *current = end;
  while (current) {
    // printf("Path voxel: (%d, %d, %d)\n", current->x, current->y, current->z);
    printNode(current, grid);
    current = current->parent;
  }
}

static inline int clampIndex(int value, int limit) {
  if (limit <= 0) {
    return 0;
  }
  if (value < 0) {
    return 0;
  }
  if (value >= limit) {
    return limit - 1;
  }
  return value;
}

static Node *findClosestWalkableVertex(Point_t point, Grid_t *grid, double *distance_out) {
  double scaled_x = (point.x / STEP_GRID_VOXEL) + grid->offset_x;
  double scaled_y = (point.y / STEP_GRID_VOXEL) + grid->offset_y;
  double scaled_z = (point.z / STEP_GRID_VOXEL) + grid->offset_z;

  int base_x = (int)floor(scaled_x);
  int base_y = (int)floor(scaled_y);
  int base_z = (int)floor(scaled_z);

  if (base_x < 0 || base_y < 0 || base_z < 0 || base_x + 1 >= grid->width || base_y + 1 >= grid->height ||
      base_z + 1 >= grid->depth) {
    return NULL;
  }

  double best_dist = DBL_MAX;
  Node *best_node = NULL;

  for (int dx = 0; dx <= 1; dx++) {
    for (int dy = 0; dy <= 1; dy++) {
      for (int dz = 0; dz <= 1; dz++) {
        int vx = base_x + dx;
        int vy = base_y + dy;
        int vz = base_z + dz;

        Node *vertex = &grid->nodes[vx][vy][vz];
        if (!vertex->walkable) {
          continue;
        }

        Point_t vertex_world = {VOXEL_TO_WORLD(vx, grid->offset_x, STEP_GRID_VOXEL),
                                VOXEL_TO_WORLD(vy, grid->offset_y, STEP_GRID_VOXEL),
                                VOXEL_TO_WORLD(vz, grid->offset_z, STEP_GRID_VOXEL)};
        double tmp_dist = dist(point, vertex_world);
        if (tmp_dist < best_dist) {
          best_node = vertex;
          best_dist = tmp_dist;
        }
      }
    }
  }

  if (best_node && distance_out) {
    *distance_out = best_dist;
  }
  return best_node;
}

static Node *findClosestWalkableNeighbor(Point_t point, Grid_t *grid, double *distance_out) {
  Node *vertex = findClosestWalkableVertex(point, grid, distance_out);
  if (vertex) {
    return vertex;
  }

  int center_x = clampIndex(WORLD_TO_VOXEL(point.x, grid->offset_x, STEP_GRID_VOXEL), grid->width);
  int center_y = clampIndex(WORLD_TO_VOXEL(point.y, grid->offset_y, STEP_GRID_VOXEL), grid->height);
  int center_z = clampIndex(WORLD_TO_VOXEL(point.z, grid->offset_z, STEP_GRID_VOXEL), grid->depth);

  double best_dist = DBL_MAX;
  Node *best_node = NULL;

  for (int dx = -1; dx <= 1; dx++) {
    int neighbor_x = center_x + dx;
    if (neighbor_x < 0 || neighbor_x >= grid->width) {
      continue;
    }
    for (int dy = -1; dy <= 1; dy++) {
      int neighbor_y = center_y + dy;
      if (neighbor_y < 0 || neighbor_y >= grid->height) {
        continue;
      }
      for (int dz = -1; dz <= 1; dz++) {
        int neighbor_z = center_z + dz;
        if (neighbor_z < 0 || neighbor_z >= grid->depth) {
          continue;
        }

        Node *neighbor = &grid->nodes[neighbor_x][neighbor_y][neighbor_z];
        if (!neighbor->walkable) {
          continue;
        }

        Point_t neighbor_world = {VOXEL_TO_WORLD(neighbor_x, grid->offset_x, STEP_GRID_VOXEL),
                                  VOXEL_TO_WORLD(neighbor_y, grid->offset_y, STEP_GRID_VOXEL),
                                  VOXEL_TO_WORLD(neighbor_z, grid->offset_z, STEP_GRID_VOXEL)};
        double tmp_dist = dist(point, neighbor_world);
        if (tmp_dist < best_dist) {
          best_node = neighbor;
          best_dist = tmp_dist;
        }
      }
    }
  }

  if (best_node && distance_out) {
    *distance_out = best_dist;
  }
  return best_node;
}

/**
 * @brief Compute the A* distance between a start point and an end point.
 *
 * This function calculates the shortest path distance between a start point (in world coordinates)
 * and an end point (in voxel coordinates) using the A* algorithm. The start point is converted to
 * the nearest voxel in the grid if it is outside the grid bounds.
 *
 * @param start The start point in world coordinates.
 * @param end The end point in voxel coordinates.
 * @param grid Pointer to the 3D grid structure.
 * @param heap Pointer to the MinHeap used for the A* algorithm.
 * @return The shortest path distance between the start and end points, or -1 if no path exists.
 */
double aStarDistance(Point_t start, Point_t end, Grid_t *grid, MinHeap_t *heap) {
  // Initialize the MinHeap and reset the grid state
  clearMinHeap(heap);
  if (grid->numVisited != 0) {
    resetGridState(grid);
  }

  double end_dist = DBL_MAX;
  Node *end_node = findClosestWalkableNeighbor(end, grid, &end_dist);
  if (end_node == NULL) {
    return -1; // No path exists, end point is not walkable
  }

  double start_dist = DBL_MAX;
  Node *start_node = findClosestWalkableNeighbor(start, grid, &start_dist);
  if (start_node == NULL) {
    printf("Error: No walkable node found for start point\n");
    printf("Start point: (%f, %f, %f)\n", start.x, start.y, start.z);
    return -1;
  }

  start_node->gCost = start_dist;
  start_node->hCost = heuristic(start_node, end_node);
  start_node->fCost = start_node->gCost + start_node->hCost;
  start_node->opened = 1;
  addVisitedNode(grid, start_node);
  insertMinHeap(heap, start_node);
  

  // A* algorithm loop
  while (heap->size > 0) {
    // Extract the node with the lowest fCost
    Node *current = extractMin(heap);
    current->closed = 1;
    // Track closed nodes for efficient reset (if not already tracked)
    if (current->opened == 0) {
      addVisitedNode(grid, current);
    }

    // Check if we reached the end point
    if (current->x == end_node->x && current->y == end_node->y && current->z == end_node->z) {
      return current->gCost + end_dist; // Return the shortest path distance + end distance from voxel to world
    }

    // Explore neighbors
    for (int dx = -1; dx <= 1; dx++) {
      for (int dy = -1; dy <= 1; dy++) {
        for (int dz = -1; dz <= 1; dz++) {
          if (dx == 0 && dy == 0 && dz == 0) {
            continue; // Skip the current node
          }

          int neighbor_x = current->x + dx;
          int neighbor_y = current->y + dy;
          int neighbor_z = current->z + dz;

          // Check if the neighbor is within bounds
          if (neighbor_x < 0 || neighbor_x >= grid->width || neighbor_y < 0 || neighbor_y >= grid->height ||
              neighbor_z < 0 || neighbor_z >= grid->depth) {
            continue;
          }

          Node *neighbor = &grid->nodes[neighbor_x][neighbor_y][neighbor_z];

          // Skip non-walkable or already closed nodes
          if (!neighbor->walkable || neighbor->closed) {
            continue;
          }

          // Calculate the tentative gCost
          double tentative_g_cost = current->gCost + moveCost(current, neighbor);

          // If the neighbor is not opened or the tentative gCost is better
          if (!neighbor->opened || tentative_g_cost < neighbor->gCost) {
            neighbor->parent = current;

            if (!neighbor->opened) {
              neighbor->gCost = tentative_g_cost;
              neighbor->hCost = heuristic(neighbor, end_node);
              neighbor->fCost = neighbor->gCost + neighbor->hCost;
              neighbor->opened = 1;
              addVisitedNode(grid, neighbor); // Track this node for efficient reset
              insertMinHeap(heap, neighbor);
            } else {
              decreaseKeyMinHeap(heap, neighbor, tentative_g_cost);
            }
          }
        }
      }
    }
  }

  // If we exit the loop, no path was found
  return -1;
}

/**
 * @brief Check if a candidate set contains a specific node.
 *
 * This function checks if a given node is present in the candidate set.
 *
 * @param candidates Array of candidate nodes.
 * @param nb_candidates Number of candidates in the array.
 * @param node The node to check for presence in the candidate set.
 * @return 1 if the node is found in the candidate set, 0 otherwise.
 */
int candidateSetContains(Node **candidates, int nb_candidates, Node *node) {
  for (int i = 0; i < nb_candidates; i++) {
    if (node->x == candidates[i]->x && node->y == candidates[i]->y && node->z == candidates[i]->z) {
      return 1;
    }
  }
  return 0;
}

/**
 * @brief Compute the A* distance between multiple start points and an end point.
 *
 * This function calculates the shortest path distance between multiple start points (in world coordinates)
 * and an end point (in world coordinates) using the SSMTA* algorithm.
 *
 * @param end The end point in world coordinates.
 * @param paths Pointer to the Paths_t structure containing the start points, candidates, results, grid, and min heap.
 * @param num_starts Number of start points.
 * @return The number of results found.
 */
int dSSMTAstar(Point_t end, Paths_t *paths, int num_starts, int growth_limit) {
  int nb_results = 0;
  int nb_candidates = 0;
  int nb_candidates_results = 0;

  Point_t *starts = paths->starts;
  Node **candidates = paths->Candidates;
  Point_t *results = paths->results_pos;
  Grid_t *grid = paths->grids[paths->currentPath];
  MinHeap_t *heap = paths->minHeaps[paths->currentPath];
  double *distances_multi_candidates = paths->distancesMultiCandidates;

  // print starts
  // for (int i = 0; i < num_starts; i++) {
  //   printf("Start %d: (%f, %f, %f)\n", i, starts[i].x, starts[i].y, starts[i].z);
  // }
  // printf("End: (%f, %f, %f)\n", end.x, end.y, end.z);

  // Initialize the MinHeap and reset the grid state if needed
  clearMinHeap(heap);
  if (grid->numVisited != 0) {
    resetGridState(grid);
  }

  // Gestion of the candidates from the start points

  // To do adapt for the SSMTA* algorithm
  //  Convert the start point from world coordinates to voxel coordinates
  int boundary_enabled = is_path_boundary_filter_enabled();
  double boundary_limit = 0.0;
  if (boundary_enabled) {
    int remaining_slots = growth_limit - paths->curPthPos[paths->currentPath]+1; //include slot at growth_limit
    if (remaining_slots < 0) {
      remaining_slots = 0;
    }
    boundary_limit = DIST_PATH_BOUNDARY(remaining_slots);
  }

  for (int i = 0; i < num_starts; i++) {

    // Test to cut branch if the distance is too large, not use this candidate
    double computed_dist = dist(starts[i], end);
    if (!boundary_enabled || computed_dist <= boundary_limit) {
#ifdef ENABLE_STATS
      if (boundary_enabled) {
        nb_path_boundary_allowed++;
      }
#endif

      Node *start_node = findClosestWalkableNeighbor(starts[i], grid, NULL);
      if (start_node == NULL) {
        printf("Error: No walkable node found for start point %d\n", i);
        printf("Start point: (%f, %f, %f)\n", starts[i].x, starts[i].y, starts[i].z);
        continue;
      }
      if (start_node->nbTimesCandidate == 0) {
        addVisitedNode(grid, start_node);
      }
      if (start_node->nbTimesCandidate >= MAX_POSITION_KEEP_360) {
        fprintf(stderr, "Warning: candidate node hit capacity (%d) while mapping start %d; skipping this start.\n",
                MAX_POSITION_KEEP_360, i);
        continue;
      }
      start_node->indexStartsCandidates[start_node->nbTimesCandidate] = i;
      start_node->nbTimesCandidate++;
      if (start_node->nbTimesCandidate == 1) {
        candidates[nb_candidates] = start_node;
        nb_candidates++;
      }
    } else {
  #ifdef ENABLE_STATS
      nb_path_boundary_blocked++;
  #endif
    }
  }

  if (nb_candidates == 0) {
    return 0; // No path exists, end point is not walkable, results is empty.
  }

  // End point in voxel coordinates
  double end_dist = DBL_MAX;
  Node *end_node = findClosestWalkableNeighbor(end, grid, &end_dist);

  if (end_node == NULL) {
    return 0; // No path exists, end point is not walkable, results is empty.
  }

  // Compute the minimum Euclidean distance from end to any candidate
  double dmin = DBL_MAX;
  for (int i = 0; i < nb_candidates; i++) {
    double dist_to_candidate = heuristic(end_node, candidates[i]);
    if (dist_to_candidate < dmin) {
      dmin = dist_to_candidate;
    }
  }

  // Put cost of end and add to the PQ
  end_node->gCost = end_dist;
  end_node->hCost = dmin;
  end_node->fCost = end_node->gCost + end_node->hCost;
  end_node->opened = 1;
  addVisitedNode(grid, end_node); // Track this node for efficient reset
  // printNode(end_node, grid);
  insertMinHeap(heap, end_node);

  // SSMTA* algorithm loop
  while (heap->size > 0) {
    // Extract the node with the lowest fCost
    Node *current = extractMin(heap);
    if (current->closed == 1) {
      continue; // Skip already closed nodes
    }
    current->closed = 1;
    // Track closed nodes for efficient reset (if not already tracked)
    if (current->opened == 0) {
      addVisitedNode(grid, current);
    }

    // Check if we reached a candidate point
    if (candidateSetContains(candidates, nb_candidates, current)) {
      nb_candidates_results++;
      if (current->nbTimesCandidate > 1) {
        // I need to find the points in starts corresponding to the candidate and calculating their distance to sort
        // it in ascending order.
        for (int i = 0; i < current->nbTimesCandidate; i++) {
          int index_start = current->indexStartsCandidates[i];
          double dist_to_start =
              current->gCost +
              dist(starts[index_start], (Point_t){VOXEL_TO_WORLD(current->x, grid->offset_x, STEP_GRID_VOXEL),
                                                  VOXEL_TO_WORLD(current->y, grid->offset_y, STEP_GRID_VOXEL),
                                                  VOXEL_TO_WORLD(current->z, grid->offset_z, STEP_GRID_VOXEL)});
          distances_multi_candidates[index_start] = dist_to_start;
          results[nb_results] = starts[index_start];
          nb_results++;
        }

      } else { // printf("Candidate found: (%d, %d, %d)\n", current->x, current->y, current->z);
        // I need to find the world coordinates of the candidate, so the point in starts corresponding to the
        // candidate.
        int index_start = current->indexStartsCandidates[0];
        double dist_to_start =
            current->gCost +
            dist(starts[index_start], (Point_t){VOXEL_TO_WORLD(current->x, grid->offset_x, STEP_GRID_VOXEL),
                                                VOXEL_TO_WORLD(current->y, grid->offset_y, STEP_GRID_VOXEL),
                                                VOXEL_TO_WORLD(current->z, grid->offset_z, STEP_GRID_VOXEL)});
        distances_multi_candidates[index_start] = dist_to_start;
        results[nb_results] = starts[index_start];
        nb_results++;
      }
      // If all candidates are processed, return the result set
      if (nb_candidates_results == nb_candidates) {
        // printf("All candidates reached.\n");
        // Sort the distances in ascending order
  //      printf("Sorting %d results based on distances.\n", num_starts);

        // Total candidates 
        int total_candidates = 0;
        for (int i = 0; i < nb_candidates; i++) {
          total_candidates += candidates[i]->nbTimesCandidate;
        }

        int pos_res[total_candidates];
        int pos = 0;
        for (int i = 0; i < nb_candidates; i++) {
          for (int j = 0; j < candidates[i]->nbTimesCandidate; j++) {
            pos_res[pos] = candidates[i]->indexStartsCandidates[j];
            pos++;
          }
        }

        // for (int i = 0; i < num_starts; i++) {
        //   printf("Distance to candidate %d: %f\n", pos_res[i], distances_multi_candidates[pos_res[i]]);
        // }

        // Sort the results based on distances
        for (int i = 0; i < total_candidates - 1; i++) {
          for (int j = i + 1; j < total_candidates; j++) {
            if (distances_multi_candidates[pos_res[i]] > distances_multi_candidates[pos_res[j]]) {
              // Swap the indices
              int temp = pos_res[i];
              pos_res[i] = pos_res[j];
              pos_res[j] = temp;
            }
          }
        }

        // Recopy starts in right order
        for (int i = 0; i < total_candidates; i++) {
          // printf("Distance to candidate %d: %f\n", pos_res[i], distances_multi_candidates[pos_res[i]]);
          results[i] = starts[pos_res[i]];
        }
        // exit(0); // For debug purpose, remove this line in production
        return nb_results;
      }
      // Renew the heap with the optimized fast version
      renewMinHeapFast(heap, candidates, nb_candidates);
    }
    // Explore neighbors
    for (int dx = -1; dx <= 1; dx++) {
      for (int dy = -1; dy <= 1; dy++) {
        for (int dz = -1; dz <= 1; dz++) {
          if (dx == 0 && dy == 0 && dz == 0) {
            continue; // Skip the current node
          }

          int neighbor_x = current->x + dx;
          int neighbor_y = current->y + dy;
          int neighbor_z = current->z + dz;

          // Check if the neighbor is within bounds
          if (neighbor_x < 0 || neighbor_x >= grid->width || neighbor_y < 0 || neighbor_y >= grid->height ||
              neighbor_z < 0 || neighbor_z >= grid->depth) {
            continue;
          }

          Node *neighbor = &grid->nodes[neighbor_x][neighbor_y][neighbor_z];

          // Skip non-walkable or already closed nodes
          if (!neighbor->walkable || neighbor->closed) {
            continue;
          }

          // Calculate the tentative gCost
          double tentative_g_cost = current->gCost + moveCost(current, neighbor);

          // If the neighbor is not opened or the tentative gCost is better
          if (!neighbor->opened || tentative_g_cost < neighbor->gCost) {
            neighbor->parent = current;

            if (!neighbor->opened) {
              neighbor->gCost = tentative_g_cost;
              // find distance min to candidate
              double dmin_nn = DBL_MAX;
              for (int i = 0; i < nb_candidates; i++) {
                // verify candidate is not closed
                if (candidates[i]->closed) {
                  continue;
                }
                double dist_nn_to_candidate = heuristic(neighbor, candidates[i]);
                if (dist_nn_to_candidate < dmin_nn) {
                  dmin_nn = dist_nn_to_candidate;
                }
              }
              neighbor->hCost = dmin_nn;
              neighbor->fCost = neighbor->gCost + neighbor->hCost;
              neighbor->opened = 1;
              addVisitedNode(grid, neighbor); // Track this node for efficient reset
              insertMinHeap(heap, neighbor);
            } else {
              decreaseKeyMinHeap(heap, neighbor, tentative_g_cost);
            }
          }
        }
      }
    }
  }
  // If we exit the loop, not all candidates were reached
  return nb_results; // The results_pos return is not necessary ordered, normally this return is don't use
}

/**
 * @brief Check if a line segment intersects with a sphere in 3D space.
 * This function checks if the line segment defined by two points (x1, y1, z1) and (x2, y2, z2)
 * intersects with a sphere defined by its center (cx, cy, cz) and radius r.
 * @param x1, y1, z1 Coordinates of the first point of the line segment.
 * @param x2, y2, z2 Coordinates of the second point of the line segment.
 * @param cx, cy, cz Coordinates of the center of the sphere.
 * @param r Radius of the sphere, distance of collision.
 * @return 1 if the line segment intersects with the sphere, 0 otherwise.
 */
static inline int segmentIntersectsSphereFast(double x1, double y1, double z1, double x2, double y2, double z2,
                                              double cx, double cy, double cz, double r) {
  double dx = x2 - x1;
  double dy = y2 - y1;
  double dz = z2 - z1;

  double ox = x1 - cx;
  double oy = y1 - cy;
  double oz = z1 - cz;

  double a = dx * dx + dy * dy + dz * dz;
  double b = 2.0 * (dx * ox + dy * oy + dz * oz);
  double c = ox * ox + oy * oy + oz * oz - r * r;

  double discriminant = b * b - 4.0 * a * c;

  if (discriminant < 0.0)
    return 0; // no real roots → no intersection

  double sqrt_discriminant = sqrt(discriminant);
  double inv_2a = 0.5 / a;

  double t1 = (-b - sqrt_discriminant) * inv_2a;
  double t2 = (-b + sqrt_discriminant) * inv_2a;

  // Check if either intersection point lies within segment [0,1]
  return (t1 >= 0.0 && t1 <= 1.0) || (t2 >= 0.0 && t2 <= 1.0);
}

/**
 * @brief Check for collisions between a line segment and spheres in a cage, paths, and substrate.
 *
 * This function checks if the line segment defined by two points intersects with any spheres
 * in the cage, paths, or substrate. If a collision is detected, it returns a specific distance type.
 *
 * @param start The starting point of the line segment.
 * @param end The ending point of the line segment.
 * @param cage Pointer to the cage structure containing spheres.
 * @param paths Pointer to the paths structure containing patterns.
 * @param grid_sub Pointer to the substrate grid structure.
 * @param substrat_t Pointer to the substrate data.
 * @return DistanceType indicating the type of distance based on collisions detected.
 */
DistanceType distanceHybridLineSpheresCollisionTest(Point_t start, Point_t end, Cage_t *cage, Paths_t *paths,
                                                    GridSubstrat *grid_sub, double ***substrat_t) {
  // Check for substrat collision
  for (int i = 0; i < grid_sub->substratSize; i++) {
    if (segmentIntersectsSphereFast(start.x, start.y, start.z, end.x, end.y, end.z, (*substrat_t)[i][0],
                                    (*substrat_t)[i][1], (*substrat_t)[i][2], DIST_GAP_SUBSTRATE)) {
      return DISTANCE_SSMTA_STAR; // Use SSMTA* distance if collision with substrate
    }
  }
  Point_t sphere_fix;
  // Check for cage collision
  for (int i = 0; i < size(cage); i++) {
    sphere_fix = coords(atom(cage, i));
    // check sphere_fix is not the same as start and end
    if (sphere_fix.x == start.x && sphere_fix.y == start.y && sphere_fix.z == start.z) {
      continue; // Skip if the sphere is the same as start
    }
    if (sphere_fix.x == end.x && sphere_fix.y == end.y && sphere_fix.z == end.z) {
      continue; // Skip if the sphere is the same as end
    }
    if (segmentIntersectsSphereFast(start.x, start.y, start.z, end.x, end.y, end.z, sphere_fix.x, sphere_fix.y,
                                    sphere_fix.z, DIST_GAP_CAGE)) {
      return DISTANCE_SSMTA_STAR; // Use SSMTA* distance if collision with cage
    }
  }
  // Check for paths collision

  // paths processing
  int size_p;
  for (int k = 0; k <= paths->currentPath; k++) {
    size_p = paths->curPthPos[k]; // IDK if we check the end atom of path for his hydrogens.
    for (int i = 2; i < size_p; i++) {
      for (int j = 0; j < MAX_NB_ATOMS_PATTERN; j++) {
        sphere_fix = paths->patterns[indexPointPaths(
            k, i, paths->positionCurNum[(indexPathPosition(k, i, paths->sizeMax))], j, paths->sizeMax)];
        // check sphere_fix is not the same as start and end
        if (sphere_fix.x == start.x && sphere_fix.y == start.y && sphere_fix.z == start.z) {
          continue; // Skip if the sphere is the same as start
        }
        if (sphere_fix.x == end.x && sphere_fix.y == end.y && sphere_fix.z == end.z) {
          continue; // Skip if the sphere is the same as end
        }
        if (segmentIntersectsSphereFast(start.x, start.y, start.z, end.x, end.y, end.z, sphere_fix.x, sphere_fix.y,
                                        sphere_fix.z, DIST_GAP_CAGE)) {
          return DISTANCE_SSMTA_STAR; // Use SSMTA* distance if collision with cage
        }
      }
    }
  }
  for (int k = 0; k < paths->currentPath; k++) {
    size_p = paths->curPthPos[k];
    for (int i = size_p; i <= size_p + 2; i++) {
      if (i == size_p) {
        sphere_fix = paths->patterns[indexPointPaths(
            k, i, paths->positionCurNum[(indexPathPosition(k, i, paths->sizeMax))], 0, paths->sizeMax)];
        // check sphere_fix is not the same as start and end
        if (sphere_fix.x == start.x && sphere_fix.y == start.y && sphere_fix.z == start.z) {
          continue; // Skip if the sphere is the same as start
        }
        if (sphere_fix.x == end.x && sphere_fix.y == end.y && sphere_fix.z == end.z) {
          continue; // Skip if the sphere is the same as end
        }
        if (segmentIntersectsSphereFast(start.x, start.y, start.z, end.x, end.y, end.z, sphere_fix.x, sphere_fix.y,
                                        sphere_fix.z, DIST_GAP_CAGE)) {
          return DISTANCE_SSMTA_STAR; // Use SSMTA* distance if collision with cage
        }
      }
      sphere_fix = paths->patterns[indexPointPaths(
          k, i, paths->positionCurNum[(indexPathPosition(k, i, paths->sizeMax))], 1, paths->sizeMax)];
      // check sphere_fix is not the same as start and end
      if (sphere_fix.x == start.x && sphere_fix.y == start.y && sphere_fix.z == start.z) {
        continue; // Skip if the sphere is the same as start
      }
      if (sphere_fix.x == end.x && sphere_fix.y == end.y && sphere_fix.z == end.z) {
        continue; // Skip if the sphere is the same as end
      }
      if (segmentIntersectsSphereFast(start.x, start.y, start.z, end.x, end.y, end.z, sphere_fix.x, sphere_fix.y,
                                      sphere_fix.z, DIST_GAP_CAGE)) {
        return DISTANCE_SSMTA_STAR; // Use SSMTA* distance if collision with cage
      }

      sphere_fix = paths->patterns[indexPointPaths(
          k, i, paths->positionCurNum[(indexPathPosition(k, i, paths->sizeMax))], 2, paths->sizeMax)];
      // check sphere_fix is not the same as start and end
      if (sphere_fix.x == start.x && sphere_fix.y == start.y && sphere_fix.z == start.z) {
        continue; // Skip if the sphere is the same as start
      }
      if (sphere_fix.x == end.x && sphere_fix.y == end.y && sphere_fix.z == end.z) {
        continue; // Skip if the sphere is the same as end
      }
      if (segmentIntersectsSphereFast(start.x, start.y, start.z, end.x, end.y, end.z, sphere_fix.x, sphere_fix.y,
                                      sphere_fix.z, DIST_GAP_CAGE)) {
        return DISTANCE_SSMTA_STAR; // Use SSMTA* distance if collision with cage
      }
    }
  }

  return DISTANCE_EUCLIDEAN; // Default to Euclidean distance
}

/**
 * @brief Set the distance function based on the environment variable DISTANCE_TYPE.
 *
 * This function checks the environment variable DISTANCE_TYPE and sets the distance function
 * accordingly. If the variable is set to "A*", it uses the A* distance function; otherwise, it defaults
 * to the Euclidean distance function.
 */
void set_distance_function_from_env() {
  char *env = getenv("DISTANCE_TYPE");
  if (env) {
    if (strcmp(env, "A*") == 0) {
      // distance_func = aStarDistance;
      current_distance_type = DISTANCE_A_STAR;
      printf("Using A* distance.\n");
    } else if (strcmp(env, "SSMTA*") == 0) {
      // distance_func = aStarDistance;
      current_distance_type = DISTANCE_SSMTA_STAR;
      printf("Using SSMTA* distance.\n");
    } else if (strcmp(env, "HYBRID") == 0) {
      // distance_func = aStarDistance;
      current_distance_type = DISTANCE_HYBRID;
      printf("Using Hybrid distance.\n");
    } else {
      distance_func = dist;
      current_distance_type = DISTANCE_EUCLIDEAN;
      printf("Using Euclidean distance.\n");
    }
  } else {
    printf("Using default Euclidean distance.\n");
  }
}

/**
 * @brief Get the current distance type.
 *
 * This function returns the current distance type being used in the program.
 *
 * @return The current DistanceType enum value.
 */
DistanceType get_current_distance_type() { return current_distance_type; }

/**
 * @brief Convert a DistanceType enum to a string representation.
 *
 * This function takes a DistanceType enum value and returns its string representation.
 *
 * @param type The DistanceType enum value to convert.
 * @return A string representing the distance type.
 */
const char *distance_type_to_string(DistanceType type) {
  switch (type) {
  case DISTANCE_EUCLIDEAN:
    return "Euclidean";
  case DISTANCE_A_STAR:
    return "A*";
  case DISTANCE_SSMTA_STAR:
    return "SSMTA*";
  case DISTANCE_HYBRID:
    return "HYBRID";
  default:
    return "Unknown";
  }
}

static int parse_bool_flag(const char *value, int default_value) {
  if (!value || value[0] == '\0') {
    return default_value;
  }
  char c = (char)tolower((unsigned char)value[0]);
  if (c == '0') {
    return 0;
  }
  if (c == '1') {
    return 1;
  }
  if (c == 'f' || c == 'n') {
    return 0;
  }
  if (c == 't' || c == 'y') {
    return 1;
  }
  return default_value;
}

void configure_path_boundary_filter_from_env() {
  const char *explicit_flag = getenv("CAGE_PATH_BOUNDARY");
  if (explicit_flag) {
    set_path_boundary_filter_enabled(parse_bool_flag(explicit_flag, path_boundary_filter_enabled));
  }
  const char *legacy_disable = getenv("CAGE_DISABLE_PATH_BOUNDARY");
  if (legacy_disable && parse_bool_flag(legacy_disable, 0)) {
    set_path_boundary_filter_enabled(0);
  }
}

void set_path_boundary_filter_enabled(int enabled) { path_boundary_filter_enabled = enabled ? 1 : 0; }

int is_path_boundary_filter_enabled() { return path_boundary_filter_enabled; }
