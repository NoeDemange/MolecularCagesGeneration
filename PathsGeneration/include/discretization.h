#ifndef __DISCRETIZATION_H
#define __DISCRETIZATION_H

#include "structure.h"
#include "substrat.h"

// Convert from world coordinate to voxel index (grid index)
#define WORLD_TO_VOXEL(coord, offset, step) ((int)round((coord) / (step)) + (offset))

// Convert from voxel index to world coordinate (center of voxel)
#define VOXEL_TO_WORLD(index, offset, step) (((index) - (offset)) * (step))

void findBounds(Cage_t *cage, Paths_t *paths, GridSubstrat *grid_sub, double *min_x, double *max_x, double *min_y,
                double *max_y, double *min_z, double *max_z);

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
void createGrid(Grid_t *grid, Cage_t *cage, Paths_t *paths, double ***substrat_t, GridSubstrat *grid_sub);

/**
 * @brief reset values used for A* algorithm.
 *
 * This function resets the values used for the A* algorithm in the grid nodes.
 * It sets the gCost to DBL_MAX, fCost to 0, hCost to 0, parent to NULL,
 * closed to 0, and opened to 0 for all nodes in the grid.
 *
 * @param grid Pointer to the Grid_t structure representing the grid.
 */
void resetGridState(Grid_t *grid);

/**
 * @brief Add a node to the visited nodes list for efficient reset.
 *
 * This function adds a node to the visitedNodes array for tracking during pathfinding.
 * It will automatically resize the array if needed.
 *
 * @param grid Pointer to the Grid_t structure.
 * @param node Pointer to the node to mark as visited.
 */
void addVisitedNode(Grid_t *grid, Node *node);

/**
 * @brief Initialize the visited nodes tracking for a grid.
 *
 * This function initializes the visitedNodes array with an initial capacity.
 *
 * @param grid Pointer to the Grid_t structure.
 * @param initialCapacity Initial capacity for the visitedNodes array.
 */
void initVisitedNodes(Grid_t *grid, int initialCapacity);

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
void printNode(Node *node, Grid_t *grid);

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
void writeGridToMol2(Grid_t *grid, const char *filename, int writeWalkable);

/**
 * @brief Frees the memory allocated for a Grid_t structure.
 *
 * This function deallocates the memory used by the Grid_t structure, including its nodes.
 *
 * @param grid Pointer to the Grid_t structure to be freed.
 */
void freeGrid(Grid_t *grid);

#endif