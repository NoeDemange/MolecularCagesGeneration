#ifndef __DISTANCE_H
#define __DISTANCE_H

#include "structure.h"
#include "substrat.h"

typedef double (*DistanceFunc)(Point_t, Point_t);

// Enum for distance types
typedef enum { DISTANCE_EUCLIDEAN, DISTANCE_A_STAR, DISTANCE_SSMTA_STAR, DISTANCE_HYBRID } DistanceType;

// Function declarations
/**
 * @brief Squared Euclidean Distance between two points.
 *
 * @param point1 Double table of the 3 coords.
 * @param point2 Double table of the 3 coords.
 */
double squaredEuclideanDistance(double point1[3], double point2[3]);

/**
 * @brief Squared Euclidean Distance between three double for point 1 and a table of double for point 2.
 *
 * @param p1x Double x coord of point 1.
 * @param p1y Double y coord of point 1.
 * @param p1z Double z coord of point 1.
 * @param point2 Double table of the 3 coords.
 */
double squaredEuclideanDistanceCoordsPoint(double p1x, double p1y, double p1z, double point2[3]);

/**
 * @brief Squared Euclidean Distance between three double for point 1 and a table of double for point 2.
 *
 * @param p1 Point_t structure
 * @param point2 Double table of the 3 coords.
 */
double squaredEuclideanDistancefloatPointT(Point_t p1, double point2[3]);

/**
 * @brief Compute the Square Euclidean distance between two points.
 *
 * This function calculates the Square Euclidean distance between two 3D points A and B.
 *
 * @param A The first point.
 * @param B The second point.
 * @return The Euclidean distance between A and B.
 */
double squaredEuclideanDistancePointTPointT(Point_t A, Point_t B);

/**
 * @brief Compute the Euclidean distance between two points.
 *
 * This function calculates the Euclidean distance between two 3D points A and B.
 *
 * @param A The first point.
 * @param B The second point.
 * @return The Euclidean distance between A and B.
 */
double dist(Point_t, Point_t);

int distInf(Point_t A, Point_t B, double dist);

double heuristic(Node *a, Node *b);

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
double aStarDistance(Point_t start, Point_t end, Grid_t *grid, MinHeap_t *heap);

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
int dSSMTAstar(Point_t end, Paths_t *paths, int num_starts, int growth_limit);

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
                                                    GridSubstrat *grid_sub, double ***substrat_t);

// Global function pointer and type variable
extern DistanceFunc distance_func;
extern DistanceType current_distance_type;

// Functions to set and get distance function

/**
 * @brief Set the distance function based on the environment variable DISTANCE_TYPE.
 *
 * This function checks the environment variable DISTANCE_TYPE and sets the distance function
 * accordingly. If the variable is set to "A*", it uses the A* distance function; otherwise, it defaults
 * to the Euclidean distance function.
 */
void set_distance_function_from_env();

/**
 * @brief Get the current distance type.
 *
 * This function returns the current distance type being used in the program.
 *
 * @return The current DistanceType enum value.
 */
DistanceType get_current_distance_type();

/**
 * @brief Convert a DistanceType enum to a string representation.
 *
 * This function takes a DistanceType enum value and returns its string representation.
 *
 * @param type The DistanceType enum value to convert.
 * @return A string representing the distance type.
 */
const char *distance_type_to_string(DistanceType type);

/**
 * @brief Configure the DIST_PATH_BOUNDARY filter from environment variables.
 */
void configure_path_boundary_filter_from_env();

/**
 * @brief Enable or disable the DIST_PATH_BOUNDARY pruning logic.
 */
void set_path_boundary_filter_enabled(int enabled);

/**
 * @brief Check if the DIST_PATH_BOUNDARY pruning is currently enabled.
 */
int is_path_boundary_filter_enabled();

#endif