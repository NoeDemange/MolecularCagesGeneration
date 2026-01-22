#include "assembly.h"
#include "constant.h"
#include "discretization.h"
#include "distance.h"
#include "intervalHandler.h"
#include "output.h"
#include "thetaSelection.h"
#include "util.h"

#include <float.h>
#include <math.h>

/**
 * @file assembly.c
 * @brief This file contains functions for generating connected cages.
 */

#ifdef ENABLE_STATS
// variable global statistique
double cumul_nb_interval = 0.0;
int max_nb_interval = 0;
int nb_use_inter = 0;
int nb_no_inter = 0;
double max_covered = 0.0;
double cumul_covered = 0.0;
int cpt_collision = 0;
int cpt_no_next_point = 0;
int nb_branches = 0;
int nb_path_boundary_allowed = 0;
int nb_path_boundary_blocked = 0;
#endif

static inline int base_growth_limit(const Paths_t *paths) {
  return paths ? (paths->sizeMax - 3) : 0;
}

int effective_growth_limit(Paths_t *paths, int path_index, Options_t options) {
  int fallback = base_growth_limit(paths);
  if (!paths || path_index < 0)
    return fallback;
  if (!options.enableDynamicPathLimit || !paths->maxGrowthLimit)
    return fallback;
  int stored = paths->maxGrowthLimit[path_index];
  if (stored <= 0 || stored > fallback)
    return fallback;
  return stored;
}

static void record_best_path_length(Paths_t *paths, int path_index) {
  if (!paths || !paths->bestPathLength || !paths->maxGrowthLimit || path_index < 0)
    return;
  int curPos = paths->curPthPos[path_index];
  int realLen = curPos - 1; // number of patterns in the path, excluding the start and start's neighbors but it starts count from 0 so 0 and 1 is not count as pattern and the value is the last case where we have a pattern (example for a path of size 3 we have patterns in case 2, 3, 4 so curPthPos =4 and number of patterns =3)
  if (realLen < 1)
    return;
  int prev = paths->bestPathLength[path_index];
  if (prev == -1 || realLen < prev) {
    paths->bestPathLength[path_index] = realLen;
    int fallback = base_growth_limit(paths);
    int new_limit = curPos;
    if (new_limit > fallback)
      new_limit = fallback;
    if (new_limit <= 0)
      new_limit = fallback;
    paths->maxGrowthLimit[path_index] = new_limit;
  }
}

/**************************************/
/******** Collision & Intervals ********/
/**************************************/

/**
 * @brief Determines if a given point (atom) is too close to the substrate,
 * the molecular cage, or any part of the path under construction.
 *
 * This function ensures spatial constraints are respected by checking if the
 * point overlaps or is within a defined minimum distance (`DIST_GAP_CAGE`)
 * from:
 * 1. The substrate atoms.
 * 2. Existing atoms in the molecular cage.
 * 3. The current path under construction.
 *
 * @param cage Pointer to the molecular cage being generated. Contains atoms and their positions.
 * @param grid_sub Pointer to the grid representation of the substrate. Used for spatial collision checks.
 * @param paths Pointer to the paths under construction.
 * @param p The point (atom) being tested for spatial constraints.
 * @return int Returns:
 * - `1` if the point is too close (hindered) to the substrate, cage, or path.
 * - `0` if the point respects the spatial constraints and is not hindered.
 *
 * @details
 * - **Substrate check**: Calls `CheckGridCollisionSubstrat_Point_t` to verify if the point is within
 *   `STEP_GRID` of any substrate atom.
 * - **Cage check**: Loops through all atoms in the cage to ensure the point is outside the
 *   minimum distance (`DIST_GAP_CAGE`).
 * - **Path check**: Iterates through all paths and their patterns to ensure the point does not
 *   overlap with existing path components.
 *
 * This function ensures that atoms maintain proper separation during the generation process,
 * avoiding overlaps or improper placements.
 */
int isHindered(Cage_t *cage, GridSubstrat *grid_sub, Paths_t *paths, Point_t p) {
  // check for substrate
  if (checkGridCollisionSubstratPointT(p, grid_sub, STEP_GRID))
    return 1;
  // check for actual cage
  for (int i = 0; i < size(cage); i++) {
    Point_t a = coords(atom(cage, i));
    if (distInf(a, p, DIST_GAP_CAGE)) {
      return 1;
    }
  }
  // check for paths under construction
  int size;
  for (int k = 0; k <= paths->currentPath; k++) {
    size = paths->curPthPos[k];
    for (int i = 2; i < size; i++) {
      // int size = (path->patternNum[i] == CYCLE_PATTERN) ? MAX_NB_ATOMS_PATTERN : 3;
      for (int j = 0; j < MAX_NB_ATOMS_PATTERN; j++) {
        if (distInf(paths->patterns[indexPointPaths(
                        k, i, paths->positionCurNum[(indexPathPosition(k, i, paths->sizeMax))], j, paths->sizeMax)],
                    p, DIST_GAP_CAGE)) // paths->patterns[i][path->positionNum[i]][j]
          return 1;
      }
    }
  }
  // For paths fully completed, check the last two atoms
  for (int k = 0; k < paths->currentPath; k++) { // Warning adapt for simple pattern
    size = paths->curPthPos[k];
    for (int i = size; i <= size + 2; i++) {
      if (i == size) {
        if (distInf(paths->patterns[indexPointPaths(
                        k, i, paths->positionCurNum[(indexPathPosition(k, i, paths->sizeMax))], 0, paths->sizeMax)],
                    p, DIST_GAP_CAGE))
          return 1;
      }
      if (distInf(paths->patterns[indexPointPaths(
                      k, i, paths->positionCurNum[(indexPathPosition(k, i, paths->sizeMax))], 1, paths->sizeMax)],
                  p, DIST_GAP_CAGE))
        return 1;
      if (distInf(paths->patterns[indexPointPaths(
                      k, i, paths->positionCurNum[(indexPathPosition(k, i, paths->sizeMax))], 2, paths->sizeMax)],
                  p, DIST_GAP_CAGE))
        return 1;
    }
  }
  return 0;
}

/* Intervals */

/**
 * @brief Safely adds an interval to the array, handling wraparound cases
 */
void addIntervalSafe(Interval *intervals, int *count, int max_size, double theta_min, double theta_max) {
  if (*count >= max_size - 1) {
    return; // Not enough space
  }

  // Handle wraparound intervals (crossing 0/2π boundary)
  if (theta_min > theta_max) {
    intervals[*count] = (Interval){theta_min, 2 * M_PI}; // Wraparound case, add two intervals
    (*count)++;
    if (*count >= max_size) {
      return;
    }
    intervals[*count] = (Interval){0, theta_max};
  } else {
    intervals[*count] = (Interval){theta_min, theta_max};
  }
  (*count)++;
}

/**
 * @brief Processes a single atom for collision detection
 */
void processAtomCollision(Point_t center, Point_t v1, Point_t v2, Point_t atom_position, double gap_distance,
                          double circle_radius_squared, double max_collision_dist_sq, Interval *intervals,
                          int *interval_count, int max_intervals) {
  // Quick distance check to avoid expensive calculations for distant atoms
  if (squaredEuclideanDistancePointTPointT(center, atom_position) >= max_collision_dist_sq) {
    return; // No collision possible, continue processing
  }

  Interval interval;
  if (findValidThetaForSphere(center, circle_radius_squared, v1, v2, atom_position, &interval, gap_distance)) {
    addIntervalSafe(intervals, interval_count, max_intervals, interval.theta_min, interval.theta_max);
  }
  return; // Success, continue processing
}

/**
 * @brief Collects all collision intervals from substrate, cage, and path atoms
 * @return Number of intervals found, or -1 on error
 */
int collectCollisionIntervals(Point_t center, Point_t v1, Point_t v2, Cage_t *cage, Paths_t *paths,
                              GridSubstrat *grid_sub, double ***substrat_t, Interval **intervals_out,
                              char *type_of_circle) {

  // Estimate the maximum number of intervals needed
  int total_size_path = 0;
  for (int k = 0; k <= paths->currentPath; k++) {
    total_size_path += paths->curPthPos[k];
  }
  int max_intervals =
      (total_size_path + grid_sub->substratSize + size(cage)) * 2;    // Factor of 2 for potential wraparound splits
  Interval *intervals_tmp = malloc(max_intervals * sizeof(Interval)); // Allocate memory for intervals
  if (!intervals_tmp) {
    printf("Error: Failed to allocate memory for collision intervals\n");
    return -1;
  }

  int interval_count = 0; // Count of intervals found

  double circle_radius = 0;
  double circle_radius_squared = 0;
  if (strcmp(type_of_circle, "carbon") == 0) {
    // Handle carbon circle case
    circle_radius = CIRCLE_RADIUS_C;
    circle_radius_squared = CIRCLE_RADIUS_SQUARED_C;

  } else if (strcmp(type_of_circle, "hydrogen") == 0) {
    // Handle hydrogen circle case
    circle_radius = CIRCLE_RADIUS_H;
    circle_radius_squared = CIRCLE_RADIUS_SQUARED_H;
  } else {
    printf("Error: Invalid type of circle '%s'\n", type_of_circle);
    free(intervals_tmp);
    *intervals_out = NULL;
    return -1; // Invalid type of circle
  }

  double max_collision_dist_sq = 0; // Maximum distance squared for collision checks

  // Process substrate atoms
  max_collision_dist_sq = (circle_radius + DIST_GAP_SUBSTRATE) * (circle_radius + DIST_GAP_SUBSTRATE);
  for (int i = 0; i < grid_sub->substratSize && interval_count < max_intervals - 1; i++) {
    Point_t substrate_position = {(*substrat_t)[i][0], (*substrat_t)[i][1], (*substrat_t)[i][2]};
    processAtomCollision(center, v1, v2, substrate_position, DIST_GAP_SUBSTRATE, circle_radius_squared,
                         max_collision_dist_sq, intervals_tmp, &interval_count, max_intervals);
  }

  max_collision_dist_sq = (circle_radius + DIST_GAP_CAGE) * (circle_radius + DIST_GAP_CAGE);
  // Process cage atoms
  for (int i = 0; i < size(cage) && interval_count < max_intervals - 1; i++) {
    Point_t cage_position = coords(atom(cage, i));
    processAtomCollision(center, v1, v2, cage_position, DIST_GAP_CAGE, circle_radius_squared, max_collision_dist_sq,
                         intervals_tmp, &interval_count, max_intervals);
  }

  // Process atoms in current path construction
  for (int k = 0; k <= paths->currentPath && interval_count < max_intervals - 1; k++) {
    int path_size = paths->curPthPos[k];

    // Process atoms in the main path body (skip start and start_neighbor at indices 0,1)
    for (int i = 2; i < path_size && interval_count < max_intervals - 1; i++) {
      for (int j = 0; j < MAX_NB_ATOMS_PATTERN && interval_count < max_intervals - 1; j++) {
        Point_t path_position = paths->patterns[indexPointPaths(
            k, i, paths->positionCurNum[indexPathPosition(k, i, paths->sizeMax)], j, paths->sizeMax)];

        processAtomCollision(center, v1, v2, path_position, DIST_GAP_CAGE, circle_radius_squared, max_collision_dist_sq,
                             intervals_tmp, &interval_count, max_intervals);
      }
    }
  }

  // Process completed path terminal atoms
  for (int k = 0; k < paths->currentPath && interval_count < max_intervals - 1; k++) {
    int path_size = paths->curPthPos[k];

    for (int i = path_size; i <= path_size + 2 && interval_count < max_intervals - 1; i++) {
      int atom_start = (i == path_size) ? 0 : 1;

      for (int j = atom_start; j < 3 && interval_count < max_intervals - 1; j++) {
        Point_t terminal_position = paths->patterns[indexPointPaths(
            k, i, paths->positionCurNum[indexPathPosition(k, i, paths->sizeMax)], j, paths->sizeMax)];

        processAtomCollision(center, v1, v2, terminal_position, DIST_GAP_CAGE, circle_radius_squared,
                             max_collision_dist_sq, intervals_tmp, &interval_count, max_intervals);
      }
    }
  }

  *intervals_out = intervals_tmp;
  return interval_count;
}

/**
 * @brief Computes the optimal angle (theta) for a point in a circular projection
 * based on the center, destination, and two vectors defining the plane of the circle.
 *
 * This function calculates the angle that minimizes the distance between the
 * projection of a point onto a circle defined by its center and two orthogonal vectors.
 *
 * @param center The center of the circle.
 * @param destination The point to project onto the circle.
 * @param v1 First vector defining the plane of the circle.
 * @param v2 Second vector defining the plane of the circle.
 * @return double The optimal angle (theta) in radians.
 *
 * @details
 * - The function computes the direction vector from the center to the destination.
 * - It then calculates the angle using the dot product of this direction with the two vectors defining the circle's
 * plane.
 * - The angle is adjusted to ensure it is within [0, 2π].
 */
double optimalTheta(Point_t center, Point_t destination, Point_t v1, Point_t v2) {
  // Calculate the vector from center to destination
  Point_t direction = vector(destination, center);
  double alpha = direction.x * v1.x + direction.y * v1.y + direction.z * v1.z;
  double beta = direction.x * v2.x + direction.y * v2.y + direction.z * v2.z;

  double theta_optimal;
  if (alpha == 0.0) { // IMPORTANT pour éviter division par 0 et nan
    theta_optimal = atan(beta);
  } else
    theta_optimal = atan(beta / alpha);

  if (alpha > 0) {
    theta_optimal = theta_optimal + M_PI;
  }
  if (theta_optimal < 0) {
    theta_optimal = theta_optimal + 2 * (M_PI);
  }

  return theta_optimal;
}

/**************************************/
/******** Projections location ********/
/**************************************/

/**
 * @brief Check if a theta angle is within any forbidden interval.
 *
 * @param theta The angle to check
 * @param intervals Array of forbidden intervals
 * @param nb_intervals Number of intervals
 * @return 1 if theta is in a forbidden interval, 0 otherwise
 */
// static int isInForbiddenInterval(double theta, Interval *intervals, int nb_intervals) {
//   for (int i = 0; i < nb_intervals; i++) {
//     if (theta > intervals[i].theta_min && theta < intervals[i].theta_max) {
//       return 1; // theta is in a forbidden interval
//     }
//   }
//   return 0; // theta is not in any forbidden interval
// }

/**
 * @brief Generate a position.
 *
 * @param theta The angle to generate position for
 * @param carbon_center Center of the carbon circle
 * @param v1 First vector defining the plane
 * @param v2 Second vector defining the plane
 * @param start Starting atom position
 * @param start_neighbor Neighbor of the starting atom
 * @param cage Molecular cage structure
 * @param grid_sub Grid substrate structure
 * @param paths Paths structure
 * @param growth_limit Maximum allowed growth for the current path
 * @param positions List to store valid positions
 * @param point_dest Destination point
 * @param distance_type Current distance calculation type
 * @param nb_starts Pointer to counter for SSMTA*
 * @param count_results Pointer to result counter
 * @return 1 if position is valid and added, 0 otherwise
 */
static int generatePosition(double theta, Point_t carbon_center, Point_t v1, Point_t v2, Point_t start,
                                    Point_t start_neighbor, Cage_t *cage, GridSubstrat *grid_sub, Paths_t *paths,
                                    int growth_limit, List_s *positions, Point_t point_dest, DistanceType distance_type,
                                    int *nb_starts, int *count_results) {
  Point_t new_carbon_position;

  // Calculate carbon position using circle equation
  new_carbon_position.x = carbon_center.x + CIRCLE_RADIUS_C * cos(theta) * v1.x + CIRCLE_RADIUS_C * sin(theta) * v2.x;
  new_carbon_position.y = carbon_center.y + CIRCLE_RADIUS_C * cos(theta) * v1.y + CIRCLE_RADIUS_C * sin(theta) * v2.y;
  new_carbon_position.z = carbon_center.z + CIRCLE_RADIUS_C * cos(theta) * v1.z + CIRCLE_RADIUS_C * sin(theta) * v2.z;

  (*count_results)++;

#ifdef ENABLE_STATS
  // Count every branch exploration attempt regardless of distance backend
  nb_branches++;
#endif

  if (distance_type != DISTANCE_SSMTA_STAR) {
    lstsAddElementInOrder(positions, new_carbon_position, point_dest, paths, growth_limit,
          distance_type_to_string(distance_type));
  } else {
    if(*nb_starts > MAX_POSITION_KEEP_360){
      printf("SSMTA*: Adding position %d at (%.3f, %.3f, %.3f)\n", *nb_starts, new_carbon_position.x,
           new_carbon_position.y, new_carbon_position.z);
    }

    paths->starts[*nb_starts] = new_carbon_position;
    (*nb_starts)++;
  }
  return 1; // Position added successfully
}

/**
 * @brief Generate candidate positions for the next carbon atom based on available intervals.
 *
 * This function explores positions around the optimal theta, testing both the optimal
 * position and nearby positions within allowed intervals. It generates positions for
 * carbon atoms and their associated hydrogens while checking collision constraints.
 *
 * @param carbon_center Center of the carbon circle
 * @param v1 First vector defining the plane of the circle
 * @param v2 Second vector defining the plane of the circle
 * @param intervals Array of valid intervals
 * @param nb_intervals Number of intervals in the array
 * @param start Starting atom position
 * @param start_neighbor Neighbor of the starting atom
 * @param cage Molecular cage structure
 * @param grid_sub Grid substrate structure
 * @param paths Paths structure
 * @param growth_limit Maximum allowed growth (dynamic boundary) for the current path
 * @param positions List to store valid positions (for non-SSMTA* algorithms)
 * @param distance_type Current distance calculation type
 * @param nb_starts Pointer to counter for SSMTA* algorithm
 * @param end Destination point for the path
 */
static void generateCandidatePositions(Point_t carbon_center, Point_t v1, Point_t v2, Interval *intervals, int nb_intervals, Point_t start, Point_t start_neighbor,
                                       Cage_t *cage, GridSubstrat *grid_sub, Paths_t *paths, int growth_limit,
                                       List_s *positions, DistanceType distance_type, int *nb_starts, Point_t end) {
  if (nb_intervals <= 0) {
    return; // No valid intervals, exit early
  }
  int count_results = 0;

  if(distance_type == DISTANCE_EUCLIDEAN || DISCRITIZATION_TYPE == 0){
    //Test if we allow enough positions to have one per interval and print a warning if not
  if (MAX_POSITION_KEEP_360 < nb_intervals) {
    printf("Warning: MAX_POSITION_KEEP_360 (%d) is less than the number of valid intervals (%d). Some intervals may not be sampled.\n",
           MAX_POSITION_KEEP_360, nb_intervals);
  }
    //Calculate the size of all intervals
    double total_interval_size = 0.0;
    int max_positions;
    for (int i = 0; i < nb_intervals; i++) {
      total_interval_size += intervals[i].theta_max - intervals[i].theta_min;
    }

    //Find shift angle to place MAX_POSITION_KEEP_360 positions evenly across total intervals
    double shift_angle = total_interval_size / MAX_POSITION_KEEP_360;
    //Check if shift_angle is less than THRESHOLD_ANGLE, if yes set it to THRESHOLD_ANGLE
    if (shift_angle < THRESHOLD_ANGLE) {
      shift_angle = THRESHOLD_ANGLE;
      max_positions = (int)(total_interval_size / shift_angle)+1; //+1 to avoid losing positions due to truncation
    }else{
      max_positions = MAX_POSITION_KEEP_360;
    }
    
    if(distance_type == DISTANCE_EUCLIDEAN){// For Euclidean distances, prioritize optimal theta
      // Find the optimal theta for the next atom
      double theta_optimal = optimalTheta(carbon_center, end, v1, v2);

      // Discretize valid intervals around theta_optimal with the adaptive shift_angle,
      // enforcing at least shift_angle between positions.
      double thetas[NUMBER_POSITION_AX1E3];
      int n_sel = discretizeAroundThetaFromValid(
          intervals,                // valid intervals (sorted, merged)
          nb_intervals,
          theta_optimal,
          shift_angle,            // minimum angle between positions
          (max_positions < NUMBER_POSITION_AX1E3 ? max_positions : NUMBER_POSITION_AX1E3),    // limit: number of positions to keep
          (int)(max_positions/2),     // search depth (k), in practice half of max_positions because we search on both sides of optimal theta
          thetas);

      for (int i = 0; i < n_sel; ++i) {
        generatePosition(thetas[i], carbon_center, v1, v2, start, start_neighbor, cage, grid_sub, paths,
              growth_limit, positions, end, distance_type, nb_starts, &count_results);
      }

    }else{
      // For non-Euclidean distances, with global discretization,
      // discretize valid intervals using shift_angle
      // and cap the total number of generated positions to max_positions.
      for (int i = 0; i < nb_intervals && count_results < max_positions; i++) {
        double theta_start = intervals[i].theta_min;
        double theta_end = intervals[i].theta_max;

        // Generate positions at regular intervals within the valid range
        for (double theta = theta_start; theta <= theta_end && count_results < max_positions; theta += shift_angle) {
          generatePosition(theta, carbon_center, v1, v2, start, start_neighbor, cage, grid_sub, paths, growth_limit,
                           positions, end, distance_type, nb_starts, &count_results);
        }
      }
    }
  }else{ 
    // For non-Euclidean distances, with local discretization
    /* We have k intervals, of size t1,t2,...,tk and n positions to place (n=MAX_POSITION_KEEP_360)
     * We suppose k <= MAX_POSITION_KEEP_360 (else we warn before)
     * We place at least one position per interval, then we distribute the remaining positions
     * proportionally to the size of each interval.
     * We have a function f : [k] -> [n] function that gives the number of positions to place in each interval.
     1. f(i) >= 1 for all i in [k]
     2. f(1) + f(2) + ... + f(k) = n
     And we search to minimize on the f(i) the following cost function :
     C= sum over i=1 to k of (ti/(2+f(i)))
      * With ti always superior to threshold_angle between two positions in the same interval
     */
      const double two_pi = 2.0 * M_PI;
      Interval normalized_intervals[nb_intervals > 0 ? nb_intervals : 1];
      int interval_count = 0;
      int has_wraparound = (nb_intervals > 1) && (fabs(intervals[0].theta_min) < EPSILON_ANGLE) &&
                           (fabs(two_pi - intervals[nb_intervals - 1].theta_max) < EPSILON_ANGLE);

      if (has_wraparound) {
        // Preserve middle intervals in order, then append merged wrap interval
        for (int i = 1; i < nb_intervals - 1; i++) {
          normalized_intervals[interval_count++] = intervals[i];
        }
        normalized_intervals[interval_count++] =
            (Interval){intervals[nb_intervals - 1].theta_min, intervals[0].theta_max + two_pi};
      } else {
        for (int i = 0; i < nb_intervals; i++) {
          normalized_intervals[interval_count++] = intervals[i];
        }
      }

      if (MAX_POSITION_KEEP_360 < interval_count) {
        printf("Warning: MAX_POSITION_KEEP_360 (%d) is less than the number of valid intervals (%d). Some intervals may not be sampled.\n",
              MAX_POSITION_KEEP_360, interval_count);
      }

      int active_intervals = (interval_count < MAX_POSITION_KEEP_360) ? interval_count : MAX_POSITION_KEEP_360;
      if (active_intervals <= 0) {
        return;
      }

      // Step 1: Calculate size of each interval with one position already allocated
      double interval_sizes[active_intervals];
      int f[active_intervals];
      for (int i = 0; i < active_intervals; i++) {
        f[i] = 1; // At least one position per interval
        interval_sizes[i] = (normalized_intervals[i].theta_max - normalized_intervals[i].theta_min) / 2.0; // Initial size
      }
      int remaining_positions = MAX_POSITION_KEEP_360 - active_intervals;
      if (remaining_positions < 0) {
        remaining_positions = 0;
      }
    // Step 2: Distribute remaining positions to minimize cost function
      while (remaining_positions > 0) {
      // Find interval with maximum size
      int max_index = 0;
        for (int i = 1; i < active_intervals; i++) {
        if (interval_sizes[i] > interval_sizes[max_index]) {
          max_index = i;
        }
      }
      //Test if the maximum size is less than threshold angle, if yes break the loop
      if(interval_sizes[max_index] < THRESHOLD_ANGLE){
        break;
      }

      // Allocate one more position to the interval with the maximum size
      f[max_index]++;
      double span = normalized_intervals[max_index].theta_max - normalized_intervals[max_index].theta_min;
      interval_sizes[max_index] = span / (1 + f[max_index]);
      remaining_positions--;
    }
    // Step 3: Generate positions based on the distribution f
    for (int i = 0; i < active_intervals; i++) {
      double theta_start = normalized_intervals[i].theta_min;
      double step = interval_sizes[i]; // Step size based on allocated positions
      for (int j = 1; j <= f[i]; j++) {
        double theta = theta_start + j * step;
        double wrapped_theta = fmod(theta, two_pi);
        if (wrapped_theta < 0) {
          wrapped_theta += two_pi;
        }
        generatePosition(wrapped_theta, carbon_center, v1, v2, start, start_neighbor, cage, grid_sub, paths,
              growth_limit, positions, end, distance_type, nb_starts, &count_results);
      }
    }


  }
}

/**
 * @brief Process the candidate results and store them in the appropriate data structures.
 *
 * This function handles the results from candidate position generation, storing them
 * in the paths structure for further processing. It handles both regular distance
 * algorithms and the SSMTA* algorithm differently.
 *
 * @param positions List of candidate positions (for non-SSMTA* algorithms)
 * @param paths Paths structure to store results
 * @param cage Molecular cage structure
 * @param end Destination point for the path
 * @param distance_type Current distance calculation type
 * @param nb_starts Number of starting positions for SSMTA*
 * @param growth_limit Maximum allowed growth for the current path (dynamic boundary)
 * @param start Starting atom position
 * @param start_neighbor Neighbor of the starting atom
 */
static void processCandidateResults(List_s *positions, Paths_t *paths, Cage_t *cage, Point_t end,
                                    DistanceType distance_type, int nb_starts, int growth_limit, Point_t start,
                                    Point_t start_neighbor) {
  Point_t new_carbon_position, hydrogen1, hydrogen2;
  int i;

  if (distance_type != DISTANCE_SSMTA_STAR) {
    // Process results for regular distance algorithms
    for (i = 0; i < NUMBER_POSITION_AX1E3 && positions->first; i++) {
      new_carbon_position = positions->first->position;

      // Store carbon position
      paths->patterns[indexPointPaths(paths->currentPath, paths->curPthPos[paths->currentPath], i, 0, paths->sizeMax)] =
          new_carbon_position;

      // Calculate and store hydrogen positions
      hydrogen1 = funAX2E2(start, start_neighbor, new_carbon_position, DIST_ATOM_H);
      hydrogen2 = funAX3E1(start, start_neighbor, new_carbon_position, hydrogen1, DIST_ATOM_H);

      paths->patterns[indexPointPaths(paths->currentPath, paths->curPthPos[paths->currentPath], i, 1, paths->sizeMax)] =
          hydrogen1;
      paths->patterns[indexPointPaths(paths->currentPath, paths->curPthPos[paths->currentPath], i, 2, paths->sizeMax)] =
          hydrogen2;

      lstsRemoveFirst(positions);
    }

    // Set the maximum number of positions found
    paths->maxPositions[indexPathPosition(paths->currentPath, paths->curPthPos[paths->currentPath], paths->sizeMax)] =
        (i > 0) ? i - 1 : -1;

    lstsDelete(positions);
  } else {
    // Process results for SSMTA* algorithm
    int nb_results = dSSMTAstar(end, paths, nb_starts, growth_limit);
    // if (nb_results == -5000) {
    //   printf("Error: SSMTA* failed to find a valid path\n");
    //   writeGridToMol2(paths->grids[paths->currentPath], "grid_NoWalkable.mol2", 0);
    //   Options_t options = {"error", "0", 9, 1};
    //   writeCageOutput(cage, paths, tmpinterTree, options);
    //   exit(EXIT_FAILURE); // Handle error appropriately
    // }

    // Store the best results from SSMTA*
    for (i = 0; i < NUMBER_POSITION_AX1E3 && i < nb_results; i++) {
      new_carbon_position = paths->results_pos[i];

      // Store carbon position
      paths->patterns[indexPointPaths(paths->currentPath, paths->curPthPos[paths->currentPath], i, 0, paths->sizeMax)] =
          new_carbon_position;

      // Calculate and store hydrogen positions
      hydrogen1 = funAX2E2(start, start_neighbor, new_carbon_position, DIST_ATOM_H);
      hydrogen2 = funAX3E1(start, start_neighbor, new_carbon_position, hydrogen1, DIST_ATOM_H);

      paths->patterns[indexPointPaths(paths->currentPath, paths->curPthPos[paths->currentPath], i, 1, paths->sizeMax)] =
          hydrogen1;
      paths->patterns[indexPointPaths(paths->currentPath, paths->curPthPos[paths->currentPath], i, 2, paths->sizeMax)] =
          hydrogen2;
    }

    // Set the maximum number of positions found
    paths->maxPositions[indexPathPosition(paths->currentPath, paths->curPthPos[paths->currentPath], paths->sizeMax)] =
        (i > 0) ? i - 1 : -1;

    lstsDelete(positions); // Clean up (initialized but not used for SSMTA*)
  }
}

/**
 * @brief Projection function for AX1E3 pattern.
 *
 * This function generates candidate positions for the next carbon atom in the AX1E3 pattern,
 * considering the current path, cage, and substrate constraints. It handles both regular
 * distance algorithms and the SSMTA* algorithm.
 *
 * @param cage Pointer to the molecular cage structure.
 * @param interTree Array of indices for the current path.
 * @param paths Pointer to the paths structure containing patterns and positions.
 * @param grid_sub Pointer to the grid representation of the substrate.
 * @param substrat_t 3D array representing substrate atom coordinates.
 * @param growth_limit Maximum allowed growth for the current path, used for dynamic boundary checks.
 */
void projectionAX1E3(Cage_t *cage, int *interTree, Paths_t *paths, GridSubstrat *grid_sub, double ***substrat_t,
                     int growth_limit) {

  // get the start atom and its neighbor = origine of the projection
  Point_t start =
      paths
          ->patterns[indexPointPaths(paths->currentPath, paths->curPthPos[paths->currentPath] - 1,
                                     paths->positionCurNum[indexPathPosition(
                                         paths->currentPath, paths->curPthPos[paths->currentPath] - 1, paths->sizeMax)],
                                     0, paths->sizeMax)];
  Point_t start_neighbor =
      paths
          ->patterns[indexPointPaths(paths->currentPath, paths->curPthPos[paths->currentPath] - 2,
                                     paths->positionCurNum[indexPathPosition(
                                         paths->currentPath, paths->curPthPos[paths->currentPath] - 2, paths->sizeMax)],
                                     0, paths->sizeMax)];

  // get the end atom = destination of the path
  Point_t end = coords(atom(cage, interTree[(2 * paths->currentPath) + 1]));

  // Plane vectors for circles
  Point_t v1, v2;
  planeVectors(start_neighbor, start, &v1, &v2);
  Point_t direction = vector(start_neighbor, start);

  // Center of the carbon circle
  Point_t carbon_direction = normalization(direction, DIST_START_CENTERCIRCLE_C);
  Point_t carbon_center = {start.x + carbon_direction.x, start.y + carbon_direction.y, start.z + carbon_direction.z};

  // Center of the hydrogen circle
  Point_t hydrogen_direction = normalization(direction, DIST_START_CENTERCIRCLE_H);
  Point_t hydrogen_center = {start.x + hydrogen_direction.x, start.y + hydrogen_direction.y,
                             start.z + hydrogen_direction.z};

  // Collect intervals for carbon circle
  Interval *intervals_tmp_carbon = NULL;
  int nb_intervals_carbon = collectCollisionIntervals(carbon_center, v1, v2, cage, paths, grid_sub, substrat_t,
                                                      &intervals_tmp_carbon, "carbon");
  if (nb_intervals_carbon < 0) {
    printf("Error collecting carbon intervals\n");
    free(intervals_tmp_carbon); // Free previously allocated intervals
    exit(EXIT_FAILURE);         // Handle error appropriately
  }

  // Collect intervals for hydrogen circle
  Interval *intervals_tmp_hydrogen = NULL;
  int nb_intervals_hydrogen = collectCollisionIntervals(hydrogen_center, v1, v2, cage, paths, grid_sub, substrat_t,
                                                        &intervals_tmp_hydrogen, "hydrogen");
  if (nb_intervals_hydrogen < 0) {
    printf("Error collecting hydrogen intervals\n");
    free(intervals_tmp_carbon);   // Free previously allocated intervals
    free(intervals_tmp_hydrogen); // Free previously allocated intervals
    exit(EXIT_FAILURE);           // Handle error appropriately
  }

  // Convert hydrogen intervals to carbon circle coordinate system
  // Use a reasonable size to avoid zero-sized arrays
  int max_converted_hydrogen_size = (nb_intervals_hydrogen > 0) ? nb_intervals_hydrogen * 4 : 1;
  Interval converted_hydrogen_intervals[max_converted_hydrogen_size];
  int nb_converted_hydrogen = 0;

  if (nb_intervals_hydrogen > 0) {
    // First merge and sort hydrogen intervals
    sortIntervals(intervals_tmp_hydrogen, nb_intervals_hydrogen);
    Interval merged_hydrogen_intervals[nb_intervals_hydrogen > 0 ? nb_intervals_hydrogen : 1];
    int nb_merged_hydrogen =
        mergeOverlappingIntervals(intervals_tmp_hydrogen, nb_intervals_hydrogen, merged_hydrogen_intervals);

    // Convert hydrogen intervals to carbon circle using ANGLE_SHIFT
    nb_converted_hydrogen = convertHydrogenIntervalsToCarbon(merged_hydrogen_intervals, nb_merged_hydrogen,
                                                             converted_hydrogen_intervals, max_converted_hydrogen_size);
    if (nb_converted_hydrogen < 0) {
      printf("Error converting hydrogen intervals to carbon\n");
      free(intervals_tmp_carbon);
      free(intervals_tmp_hydrogen);
      exit(EXIT_FAILURE);
    }
  }

  // Merge carbon intervals with converted hydrogen intervals
  int max_final_intervals_size = nb_intervals_carbon + nb_converted_hydrogen;
  if (max_final_intervals_size <= 0)
    max_final_intervals_size = 1; // Avoid zero-sized array
  Interval final_intervals[max_final_intervals_size];
  int nb_final_intervals = 0;

  if (nb_intervals_carbon > 0 && nb_converted_hydrogen > 0) { // Both carbon and converted hydrogen intervals exist
    // printf("Merging carbon and converted hydrogen intervals\n");
    // Sort carbon intervals first
    sortIntervals(intervals_tmp_carbon, nb_intervals_carbon);
    Interval merged_carbon_intervals[nb_intervals_carbon > 0 ? nb_intervals_carbon : 1];
    int nb_merged_carbon =
        mergeOverlappingIntervals(intervals_tmp_carbon, nb_intervals_carbon, merged_carbon_intervals);

    // Merge both carbon and converted hydrogen intervals
    nb_final_intervals = mergeIntervalArrays(merged_carbon_intervals, nb_merged_carbon, converted_hydrogen_intervals,
                                             nb_converted_hydrogen, final_intervals, max_final_intervals_size);
    if (nb_final_intervals < 0) {
      printf("Error merging interval arrays\n");
      free(intervals_tmp_carbon);
      free(intervals_tmp_hydrogen);
      exit(EXIT_FAILURE);
    }
  } else if (nb_intervals_carbon > 0) {
    // Only carbon intervals
    // printf("Only carbon intervals available\n");
    sortIntervals(intervals_tmp_carbon, nb_intervals_carbon);
    nb_final_intervals = mergeOverlappingIntervals(intervals_tmp_carbon, nb_intervals_carbon, final_intervals);
  } else if (nb_converted_hydrogen > 0) {
    // Only converted hydrogen intervals
    // printf("Only converted hydrogen intervals available\n");
    sortIntervals(converted_hydrogen_intervals, nb_converted_hydrogen);
    nb_final_intervals =
        mergeOverlappingIntervals(converted_hydrogen_intervals, nb_converted_hydrogen, final_intervals);
  }

  // Free temporary interval arrays
  free(intervals_tmp_carbon);
  free(intervals_tmp_hydrogen);

  // Check if all positions are blocked
  if (nb_final_intervals > 0 && final_intervals[0].theta_min == 0 && final_intervals[0].theta_max == 2 * M_PI) {
#ifdef ENABLE_STATS
   // printf("Pas d'emplacement pour un prochain sommet ! \n");
    cpt_no_next_point++;
#endif
    return;
  }

#ifdef ENABLE_STATS
  // statistiques
  if (nb_final_intervals > 0) {
    if (final_intervals[0].theta_min == 0 && final_intervals[nb_final_intervals - 1].theta_max ==
                                                 2 * M_PI) { // ATTENTION fonctionne à condition d'utiliser des double
      cumul_nb_interval += (nb_final_intervals - 1);
      nb_use_inter++;
      if (max_nb_interval < (nb_final_intervals - 1)) {
        max_nb_interval = (nb_final_intervals - 1);
      }
    } else {
      // printf("2PI %lf\n",2*M_PI);
      cumul_nb_interval += nb_final_intervals;
      nb_use_inter++;
      if (max_nb_interval < nb_final_intervals) {
        max_nb_interval = nb_final_intervals;
      }
    }
    double total_covered = 0.0;
    // Sum up the lengths of all merged intervals
    for (int i = 0; i < nb_final_intervals; i++) {
      total_covered += (final_intervals[i].theta_max - final_intervals[i].theta_min);
    }
    // Compute percentage
    total_covered = (total_covered / (2 * M_PI)) * 100.0;
    cumul_covered += total_covered;
    if (max_covered < total_covered) {
      max_covered = total_covered;
    }
  }
#endif
  // Compute valid intervals as the complement of final forbidden intervals
  int max_valid_intervals_size = nb_final_intervals + 1;
  Interval valid_intervals[max_valid_intervals_size];
  int nb_valid_intervals = complementFromMerged(final_intervals, nb_final_intervals, valid_intervals, max_valid_intervals_size);
  if (nb_valid_intervals < 0) {
      printf("Error computing valid intervals\n");
      exit(EXIT_FAILURE);
    }
  if( nb_valid_intervals == 0){
    //All positions are blocked
    printf("All positions are blocked after computing valid intervals! (Not normal)\n");
    return;
  }

  // Start working to find the next position
  // Choose the distance
  DistanceType tmp_distance_type;
  if (get_current_distance_type() == DISTANCE_HYBRID) {
    // Determine the distance type based on the collision test
    // If there is no collision, Euclidean distance type is used. If there is a collision, SSMTA* distance type is used
    // This allows to use the SSMTA* algorithm only when necessary
    tmp_distance_type = distanceHybridLineSpheresCollisionTest(start, end, cage, paths, grid_sub, substrat_t);
  } else {
    tmp_distance_type = get_current_distance_type();
  }

  // Initialize variables for position selection
  // Initialize the number of candidates for SSMTA* distance type
  int nb_starts = 0;

  // Initialize list for storing candidate positions (for non-SSMTA* algorithms)
  List_s *positions = lstsInit();

  // Generate candidate positions based on available intervals
  generateCandidatePositions(carbon_center, v1, v2, valid_intervals, nb_valid_intervals, start,
                             start_neighbor, cage, grid_sub, paths, growth_limit, positions, tmp_distance_type,
                             &nb_starts, end);

  // Process the results based on the distance type
  processCandidateResults(positions, paths, cage, end, tmp_distance_type, nb_starts, growth_limit, start,
                          start_neighbor);
}

/**************************************/
/********** Generate paths ************/
/**************************************/

/**
 * @brief Computes the positions of the next atoms in the path based on the geometry of the last added atom.
 *
 * This function determines the geometry of the next positions in a molecular path. Currently,
 * it uses the AX1E3 projection exclusively, which calculates the positions for a single atom
 * with one neighbor and its associated hydrogens. The function is designed to accommodate
 * additional projection geometries (e.g., AX2E2, AX3E1) based on the number of neighbors of
 * the last atom, but these are not yet implemented.
 *
 * @param cage Pointer to the molecular cage being generated. Used for spatial constraints.
 * @param interTree Pointer to the interconnection tree array.
 * @param paths Pointer to the `Paths_t` structure that stores path-related data, including
 *              patterns and positions of atoms.
 * @param grid_sub Pointer to the grid substrate structure, used to ensure new positions
 *                do not overlap with the substrate.
 * @param substrat_t A pointer on table of positions of substrat's atoms.
 * @param growth_limit Maximum allowed growth for the active path.
 *
 * @details
 * - **AX1E3 Projection**: The function currently calls `projectionAX1E3`, which calculates the
 *   positions of the next atom and its associated hydrogens based on a single bonded atom
 *   with one neighbor.
 * - **Planned Extensions**: The function includes placeholders for handling other geometries
 *   (AX2E2, AX3E1) based on the number of neighbors of the last added atom. These would allow
 *   for more versatile path extensions.
 *
 * @note Future implementations should extend this function to dynamically choose
 * the projection type based on the geometry of the last added atom. The commented-out code
 * indicates potential directions for such extensions.
 *
 * @see projectionAX1E3
 * @see projectionAX2E2
 * @see projectionAX3E1
 *
 * @todo Implement projections for other geometries (e.g., AX2E2 and AX3E1) and update this function
 * to dynamically choose the appropriate projection.
 */
void choosePositions(Cage_t *cage, int *interTree, Paths_t *paths, GridSubstrat *grid_sub, double ***substrat_t,
                     int growth_limit) {

  // if (path->size > 2) {
  projectionAX1E3(cage, interTree, paths, grid_sub, substrat_t, growth_limit);
  // }
  // else {
  // 	int startNbNeighbors = LST_nbElements(neighborhood(atom(processedMoc, path->idStart)));
  // 	if (startNbNeighbors == 1) {
  // 		projectionAX1E3(processedMoc, path, substrate, voxelGrid, vMap, nodeHeap);
  // 	}
  // 	else if (startNbNeighbors == 2) {
  // 		projectionAX2E2(processedMoc, path, substrate); // TODO check if we need to know the atom flag.
  // 	}
  // 	else {
  // 		projectionAX3E1(processedMoc, path, substrate);
  // 	}
  // }
}

/**
 * @brief Adds the current edge to the list of banned edges.
 *
 * This function adds the edge defined by the start and end atom IDs of the current path
 * to the list of banned edges. This is used when an edge can be construct to prevent reusing edges in subsequent path
 * generations.
 *
 * @param cage Pointer to the `Cage_t` structure representing the molecular cage.
 * @param interTree Pointer to the interconnection tree array, which maps relationships
 *                  between atoms in the cage.
 * @param num_path_not_end The index of the path not end in the interconnection tree.
 * @param list_banned_edges Array of edges that are banned from being used in paths.
 * @param size_list_banned_edges Pointer to the size of the `list_banned_edges` array, which will be updated.
 */
void addbannedEdges(Cage_t *cage, int *interTree, int num_path_not_end, int *list_banned_edges,
                    int *size_list_banned_edges) {
  int pos = *size_list_banned_edges;
  int id_start = interTree[(2 * num_path_not_end)];
  int id_end = interTree[(2 * num_path_not_end) + 1];
  // Check if the edge already exists in the banned edges list
  for (int i = 0; i < pos; i++) {
    if (list_banned_edges[i * 2] == id_start && list_banned_edges[i * 2 + 1] == id_end) {
      return; // Edge already exists, do not add again
    }         // only check this orientation, normally the edge is created only once
  }
  list_banned_edges[pos * 2] = id_start;
  list_banned_edges[pos * 2 + 1] = id_end;
  // printf("Adding banned edge: %d - %d\n", id_start + 1, id_end + 1);
  (*size_list_banned_edges)++;
}

/**
 * @brief Generates all possible paths for connecting atoms in a molecular cage.
 *
 * This function generates paths within a molecular cage based on the interconnection tree
 * and ensures that the paths adhere to spatial constraints imposed by the cage and substrate.
 * It systematically explores possible paths, handles backtracking when constraints are violated,
 * and writes valid paths to the output once the generation process is complete.
 *
 * @param cage Pointer to the `Cage_t` structure representing the molecular cage being generated.
 * @param interTree Pointer to the interconnection tree array, which maps the relationships
 *                  between atoms in the cage.
 * @param paths Pointer to the `Paths_t` structure containing the patterns, positions,
 *              and metadata for path generation.
 * @param grid_sub Pointer to the `GridSubstrat` structure representing the substrate molecule,
 *                ensuring spatial constraints are respected.
 * @param substrat_t A pointer on table of positions of substrat's atoms.
 * @param options Options for path generation and output, provided as an `Options_t` structure.
 * @param list_banned_edges Array of edges that are banned from being used in the paths.
 * @param size_list_banned_edges Size of the `list_banned_edges` array.
 *
 * @details
 * ### Key Steps in Path Generation:
 * 1. **Initialization**: The function initializes paths using `PTH_init`.
 * 2. **Path Expansion**:
 *    - Attempts to calculate the next positions for atoms in the current path using
 *      the `choosePositions` function.
 *    - Checks for spatial constraints and ensures new positions are not hindered.
 *    - Handles the addition of terminal atoms and their associated hydrogens when
 *      reaching the end of a path.
 * 3. **Backtracking**:
 *    - If a path cannot be extended further or violates constraints, the function
 *      backtracks to explore alternative paths.
 *    - This involves resetting the current path's state and revisiting earlier steps
 *      in the interconnection tree.
 * 4. **Path Completion**:
 *    - Once a path is successfully generated, it is written to the output using
 *      `writeCageOutput`.
 *    - The function ensures all possible paths in the interconnection tree are explored
 *      before termination.
 *
 * ### Constraints:
 * - Ensures that atoms in the path respect the spatial geometry and do not overlap with
 *   the cage or the substrate.
 * - Validates angles and distances when adding terminal atoms and their hydrogens.
 *
 * ### Notes:
 * The function includes commented-out code for potential extensions, such as handling
 * aromatic rings or different path patterns. These are not yet implemented.
 *
 * @todo
 * Implement support for more complex path patterns (e.g., aromatic rings).
 *
 * @see PTH_init
 * @see choosePositions
 * @see writeCageOutput
 */
void generatePaths(Cage_t *cage, int *interTree, Paths_t *paths, GridSubstrat *grid_sub, double ***substrat_t,
                   Options_t options, int *list_banned_edges, int *size_list_banned_edges) {
  pthInit(paths, interTree, cage);
  int max_path_reached = 0;
  int write_output = 0;

  while (paths->currentPath >= 0) {
    int exists_new_start = (paths->maxPositions[indexPathPosition(
                                paths->currentPath, paths->curPthPos[paths->currentPath], paths->sizeMax)] >= 0);
    int growth_limit = effective_growth_limit(paths, paths->currentPath, options);

    // if (paths->curPthPos >= 2 &&
    // paths->patternNum[indexPathPosition(paths->currentPath,paths->curPthPos,paths->sizeMax)] == CYCLE_PATTERN) {
    // 	exists_new_start = addAromaticRing(processedMoc, path, substrate);
    // }
    // else
    if (paths->curPthPos[paths->currentPath] <= growth_limit /*it's start counting from 0*/ && !exists_new_start) {
      // writeCageOutput(cage, paths, interTree, options);
      choosePositions(cage, interTree, paths, grid_sub, substrat_t, growth_limit);
      exists_new_start = (paths->maxPositions[indexPathPosition(
                              paths->currentPath, paths->curPthPos[paths->currentPath], paths->sizeMax)] >= 0);
    }
    int close_to_the_end = 0;
    if (exists_new_start) {

      Point_t new_start =
          paths
              ->patterns[indexPointPaths(paths->currentPath, paths->curPthPos[paths->currentPath],
                                         paths->positionCurNum[indexPathPosition(
                                             paths->currentPath, paths->curPthPos[paths->currentPath], paths->sizeMax)],
                                         0, paths->sizeMax)]; // path->size][path->positionNum[path->size]][0];
      int id_end = interTree[(2 * paths->currentPath) + 1];
      Point_t end = coords(atom(cage, id_end));

      // writeCageOutput(cage,paths,interTree,options);

      //calculate end distance to test if the path can be finished
      double distance_to_end = dist(new_start, end);
      if (distance_to_end < DIST_SIMPLE + DIST_ERROR) {
        close_to_the_end = 1;
        Point_t neighbor_new_start = paths->patterns[indexPointPaths(
            paths->currentPath, paths->curPthPos[paths->currentPath] - 1,
            paths->positionCurNum[indexPathPosition(paths->currentPath, paths->curPthPos[paths->currentPath] - 1,
                                                    paths->sizeMax)],
            0, paths->sizeMax)]; // path->size - 1][path->positionNum[path->size - 1]][0];
        double before_last_angle = angle(new_start, end, neighbor_new_start);
        double last_angle = angle(end, new_start, coordsNeighbor(cage, id_end, 0));
        Point_t hydrogen1 = funAX2E2(new_start, neighbor_new_start, end, DIST_ATOM_H);
        Point_t hydrogen2 = funAX3E1(new_start, neighbor_new_start, end, hydrogen1, DIST_ATOM_H);
        Point_t hydro_end1 = funAX2E2(end, coordsNeighbor(cage, id_end, 0), new_start, DIST_ATOM_H);
        Point_t hydro_end2 = funAX3E1(end, coordsNeighbor(cage, id_end, 0), new_start, hydro_end1, DIST_ATOM_H);

        if (!isHindered(cage, grid_sub, paths, hydrogen1) && !isHindered(cage, grid_sub, paths, hydrogen2) &&
            !isHindered(cage, grid_sub, paths, hydro_end1) && !isHindered(cage, grid_sub, paths, hydro_end2)) {
          if (before_last_angle >= (END_ANGLE - ANGLE_ERROR) && before_last_angle <= (END_ANGLE + ANGLE_ERROR) &&
              last_angle <= (END_ANGLE + ANGLE_ERROR) && last_angle >= (END_ANGLE - ANGLE_ERROR)) {
            paths->patterns[indexPointPaths(
                paths->currentPath, paths->curPthPos[paths->currentPath] + 1,
                paths->positionCurNum[indexPathPosition(paths->currentPath, paths->curPthPos[paths->currentPath] + 1,
                                                        paths->sizeMax)],
                1, paths->sizeMax)] = hydrogen1; // add hydrogen1 to path
            paths->patterns[indexPointPaths(
                paths->currentPath, paths->curPthPos[paths->currentPath] + 1,
                paths->positionCurNum[indexPathPosition(paths->currentPath, paths->curPthPos[paths->currentPath] + 1,
                                                        paths->sizeMax)],
                2, paths->sizeMax)] = hydrogen2; // add hydrogen2 to path
            paths->patterns[indexPointPaths(
                paths->currentPath, paths->curPthPos[paths->currentPath] + 1,
                paths->positionCurNum[indexPathPosition(paths->currentPath, paths->curPthPos[paths->currentPath] + 1,
                                                        paths->sizeMax)],
                0, paths->sizeMax)] = end; // add END to path
            paths->patterns[indexPointPaths(
                paths->currentPath, paths->curPthPos[paths->currentPath] + 2,
                paths->positionCurNum[indexPathPosition(paths->currentPath, paths->curPthPos[paths->currentPath] + 2,
                                                        paths->sizeMax)],
                1, paths->sizeMax)] = hydro_end1; // add hydrogenEND1 to path
            paths->patterns[indexPointPaths(
                paths->currentPath, paths->curPthPos[paths->currentPath] + 2,
                paths->positionCurNum[indexPathPosition(paths->currentPath, paths->curPthPos[paths->currentPath] + 2,
                                                        paths->sizeMax)],
                2, paths->sizeMax)] = hydro_end2; // add hydrogenEND2 to path
            // end is put in paths->curPthPos[k]+1 but paths->curPthPos[k] is not increase by 1. End's hydrogen is in
            // paths->curPthPos[k]+2

            //update the path MSD 
            paths->pathMSD[paths->currentPath] = ((DIST_SIMPLE - distance_to_end) * (DIST_SIMPLE - distance_to_end) +
                                                  (END_ANGLE - before_last_angle) * (END_ANGLE - before_last_angle) +
                                                  (END_ANGLE - last_angle) * (END_ANGLE - last_angle)
                                                );

            record_best_path_length(paths, paths->currentPath);

            if (paths->currentPath < paths->numPaths - 1) {
              // Change current path to progress.
              close_to_the_end = 2;
            }
            if (paths->currentPath == paths->numPaths - 1) {
              // output
              writeCageOutput(cage, paths, interTree, options);
              write_output = 1;
              if (options.oneCageByInterconnectionTree == 1) {
                // If we want to stop after the first path found.
                pthReboot(paths);
                return;
              }

              // PTH_printOrWriteAll(paths, "test_path");   //Comment when not use to debug this line and for-loop after
              // for (int i = 0; i < paths->numPaths; i++) {
              // 	    for (int j = 0; j < 2; j++) {
              // 	        printf("%d ", interTree[i*2+j]+1);
              // 	    }
              // 	    printf("\n");
              // 	}
              // exit(EXIT_SUCCESS);
            }
          }
        }
      }
    }
        if (paths->curPthPos[paths->currentPath] ==
          growth_limit /*don't count end, end's hydroge and start's neighbor, it's start counting from 0*/
        || close_to_the_end || !exists_new_start) {
      if (close_to_the_end == 2) {
        (paths->currentPath)++;
        max_path_reached = paths->currentPath; // Update max_path_reached to the last path index
        if (get_current_distance_type() != DISTANCE_EUCLIDEAN && paths->currentPath > 0) {
          createGrid(paths->grids[paths->currentPath], cage, paths, substrat_t, grid_sub);
          initMinHeap(paths->minHeaps[paths->currentPath], paths->grids[paths->currentPath]->depth *
                                                               paths->grids[paths->currentPath]->width *
                                                               paths->grids[paths->currentPath]->height);
          // char grid_name[25];
          // sprintf(grid_name, "grid_%d.mol2", paths->currentPath);
          // writeGridToMol2(paths->grids[paths->currentPath], grid_name, 0);
          // writeCageOutput(cage, paths, interTree, options);
        }
      } else {
        while ((paths->currentPath >= 0 /*&& existAPathInProgress*/ &&
                ((/*paths->patternCurNum[path->size] == CYCLE_PATTERN &&*/ paths->positionCurNum[indexPathPosition(
                      paths->currentPath, paths->curPthPos[paths->currentPath], paths->sizeMax)] ==
                  paths->maxPositions[indexPathPosition(paths->currentPath, paths->curPthPos[paths->currentPath],
                                                        paths->sizeMax)]) ||
                 paths->maxPositions[indexPathPosition(paths->currentPath, paths->curPthPos[paths->currentPath],
                                                       paths->sizeMax)] < 0))) {

          paths->patternCurNum[indexPathPosition(paths->currentPath, paths->curPthPos[paths->currentPath],
                                                 paths->sizeMax)] = 0;
          paths->positionCurNum[indexPathPosition(paths->currentPath, paths->curPthPos[paths->currentPath],
                                                  paths->sizeMax)] = 0;
          paths->maxPositions[indexPathPosition(paths->currentPath, paths->curPthPos[paths->currentPath],
                                                paths->sizeMax)] = -1;
          (paths->curPthPos[paths->currentPath])--;

          if (paths->curPthPos[paths->currentPath] <= 1) { // not count start and its neighbor
            (paths->currentPath)--;                        // backtracking, changing path
          }
        }
        if (paths->currentPath >= 0) {
          // if
          // (paths->patternCurNum[indexPathPosition(paths->currentPath,paths->curPthPos[paths->currentPath],paths->sizeMax)]
          // == SIMPLE_PATTERN) {
          // 	paths->patternCurNum[indexPathPosition(paths->currentPath,paths->curPthPos[paths->currentPath],paths->sizeMax)]
          // = CYCLE_PATTERN;
          // }
          // else {
          (paths->positionCurNum[indexPathPosition(paths->currentPath, paths->curPthPos[paths->currentPath],
                                                   paths->sizeMax)])++;
          // 	path->patternNum[path->size] = SIMPLE_PATTERN;
          // }
        }
      }
    } else {
      (paths->curPthPos[paths->currentPath])++;
    }
  }
  // PTH_printOrWritePaths(paths, NULL);
  // PTH_printOrWritePathTables(paths, NULL);
  // PTH_printOrWriteAll(paths, "test_path");
  // exit(EXIT_SUCCESS);
  if (!write_output) {
    // printf("max_path_reached: %d\n", max_path_reached);
    addbannedEdges(cage, interTree, max_path_reached, list_banned_edges, size_list_banned_edges);
  }
}