#ifndef _SUBSTRAT_H
#define _SUBSTRAT_H

#include "structure.h"

/**
 * @file substrat.h
 * @brief This file contains functions header for substrat
 */

//#define GET_NUM_POINTS_F(grid, x, y, z) ((grid).gridS[(x)][(y)][(z)][0][0]) //a double
/**
 * @brief Macro to calculate the index of a coordinate in the grid.
 *
 * This macro computes the linear index for a given 3D coordinate (i, j, k) along with
 * additional parameters (l, m) for accessing a specific coordinate in the gridS structure.
 *
 * @param grid The GridSubstrat structure containing the grid dimensions.
 * @param i The x-index in the grid.
 * @param j The y-index in the grid.
 * @param k The z-index in the grid.
 * @param l The index for the specific point in the grid.
 * @param m The coordinate index (0, 1, or 2 for x, y, z respectively).
 * @return The computed index for the specified coordinate.
 */
#define indexCoord(grid, i, j, k, l, m)                                                                                \
  (((((i) * (grid).ySize + (j)) * (grid).zSize + (k)) * ((grid).caseSize + 1) + (l)) * 3 + (m))

/**
 * @brief Macro to calculate the index of a point in the grid.
 *
 * This macro computes the linear index for a given 3D coordinate (i, j, k) along with
 * an additional parameter (l) for accessing a specific point in the gridS structure.
 *
 * @param grid The GridSubstrat structure containing the grid dimensions.
 * @param i The x-index in the grid.
 * @param j The y-index in the grid.
 * @param k The z-index in the grid.
 * @param l The index for the specific point in the grid.
 * @return The computed index for the specified point.
 */
#define indexPoint(grid, i, j, k, l)                                                                                   \
  (((((i) * (grid).ySize + (j)) * (grid).zSize + (k)) * ((grid).caseSize + 1) + (l)) * 3)

/**
 * @brief Macro to retrieve the number of points for a specific grid cell.
 *
 * This macro returns the number of points stored in the gridS for the specified 3D grid
 * coordinates (i, j, k). It retrieves this value from the first coordinate of the
 * corresponding cell in the gridS structure.
 *
 * @param grid The GridSubstrat structure containing the grid data.
 * @param i The x-index in the grid.
 * @param j The y-index in the grid.
 * @param k The z-index in the grid.
 * @return The number of points stored in the grid for the specified cell.
 */
#define getNumPointsF(grid, i, j, k) ((grid).gridS[indexCoord((grid), (i), (j), (k), 0, 0)])

/**
 * @brief Macro to obtain the address of the first point in a specific grid cell.
 *
 * This macro provides a pointer to the beginning of the points data for a specific
 * grid cell defined by the coordinates (i, j, k).
 *
 * @param grid The GridSubstrat structure containing the grid data.
 * @param i The x-index in the grid.
 * @param j The y-index in the grid.
 * @param k The z-index in the grid.
 * @return A pointer to the first point's data in the specified grid cell.
 */
#define getNumPointsAddre(grid, i, j, k) (&(grid).gridS[indexCoord((grid), (i), (j), (k), 0, 0)])

/**
 * @typedef GridS
 * @brief A type definition for a pointer to an array of floats representing grid data.
 */
typedef double *GridS;

/**
 * @struct GridSubstrat
 * @brief Represents a 3D grid containing vertex data for a substrate.
 *
 * This structure holds the dimensions and boundaries of the grid along with
 * the actual data representing points in the substrate.
 */
typedef struct GridSubstrat GridSubstrat;
struct GridSubstrat {
  int substratSize; /**< size of the substrat */
  int caseSize;     /**< Max Points by case */

  double xMin;      /**< The x min of the grid (min x vertex + dist_gap_susbtrat). */
  double yMin;      /**< The y min of the grid (min y vertex + dist_gap_susbtrat). */
  double zMin;      /**< The z min of the grid (min z vertex + dist_gap_susbtrat). */
  double xMax;      /**< The x max of the grid (max x vertex + dist_gap_susbtrat). */
  double yMax;      /**< The y max of the grid (max y vertex + dist_gap_susbtrat). */
  double zMax;      /**< The z max of the grid (max z vertex + dist_gap_susbtrat). */
  double xMinWOGap; /**< The x min of the grid (min x vertex). */
  double yMinWOGap; /**< The y min of the grid (min y vertex). */
  double zMinWOGap; /**< The z min of the grid (min z vertex). */
  double xMaxWOGap; /**< The x max of the grid (max x vertex). */
  double yMaxWOGap; /**< The y max of the grid (max y vertex). */
  double zMaxWOGap; /**< The z max of the grid (max z vertex). */
  int xSize;        /**< The x size. */
  int ySize;        /**< The y size. */
  int zSize;        /**< The z size. */

  GridS gridS; /**< The grid, [x][y][z][points][coords], points is size of substrat + 1 to stock number of vertex in
                  gridS[][][][0][0]. */
};

/**
 * Retrieves the file containing the molecule data.
 * Must have an .mol2 extension.
 *
 * @param inputname Name of the file containing the molecule data.
 * @param m Address of the molecule.
 */
double **readSubstratMol2(char *inputname, int *atom_count);

/**
 * @brief Initialize a grid for a substrate based on provided coordinates.
 *
 * This function calculates the minimum and maximum coordinates (x, y, z)
 * from a given set of substrate points, creates a 3D grid based on those
 * coordinates with a specified step size, and populates the grid with
 * the substrate points. The grid is adjusted with a defined gap to ensure
 * that the points are adequately represented in the grid.
 *
 * @param substrat_t A pointer to a 2D array containing the substrate points,
 *                   where each point is represented as a double array of size 3 (x, y, z).
 * @param substrat_size The number of substrate points contained in substrat_t.
 * @param step The size of each grid cell.
 *
 * @return GridSubstrat A structure containing the initialized grid data,
 *                      including sizes, min/max coordinates, and point data.
 *
 * @note The function allocates memory for the grid and returns a GridSubstrat structure.
 *       It is the caller's responsibility to free any allocated memory once it is no longer needed.
 *
 * @warning Ensure that the number of substrate points does not exceed the memory limits
 *          as this function uses dynamic memory allocation.
 */
GridSubstrat initGridSubstrat(double **substrat_t, int substrat_size, double step); // meilleur allocation mémoire

/**
 * @brief Checks for collision of a point with a substrate grid.
 *
 * This function determines if a given point is within a defined grid substrate.
 * It calculates the grid cell corresponding to the point and checks for collisions
 * with any points in that cell.
 *
 * @param point A 3-element array representing the coordinates of the point to check.
 *              It should be in the form [x, y, z].
 * @param gridSubstrat A pointer to the GridSubstrat structure that contains information
 *                     about the grid and its properties, such as dimensions and point data.
 * @param step The size of each grid cell. It determines how the 3D point maps to the grid.
 * @return Returns 1 if a collision is detected (i.e., the point is close to any points
 *         in the grid cell), or 0 if no collision is detected (i.e., the point is
 *         outside the grid or there are no nearby points).
 */
int checkGridCollisionSubstratFloatTable(double point[3], GridSubstrat *gridSubstrat, double step);

/**
 * @brief Checks for collision of a point with a substrate grid.
 *
 * This function determines if a given point is within a defined grid substrate.
 * It calculates the grid cell corresponding to the point and checks for collisions
 * with any points in that cell.
 *
 * @param point A Point_t structure.
 * @param gridSubstrat A pointer to the GridSubstrat structure that contains information
 *                     about the grid and its properties, such as dimensions and point data.
 * @param step The size of each grid cell. It determines how the 3D point maps to the grid.
 * @return Returns 1 if a collision is detected (i.e., the point is close to any points
 *         in the grid cell), or 0 if no collision is detected (i.e., the point is
 *         outside the grid or there are no nearby points).
 */
int checkGridCollisionSubstratPointT(Point_t p, GridSubstrat *gridSubstrat, double step);

/**
 * @brief Generates statistics about the points in the substrate grid.
 *
 * This function computes several statistics for a given GridSubstrat. It evaluates
 * the maximum, minimum, and average number of points per grid cell, as well as the
 * number of cells that are adjacent to empty cells and those that are surrounded by
 * non-empty cells.
 *
 * @param gridSubstrat The GridSubstrat structure that contains information about the grid
 *                     and its point data.
 * @return A dynamically allocated string containing the statistics of the grid. The
 *         string format includes the size of the grid, memory allocation size,
 *         maximum and minimum points per cell, and average points per cell for both
 *         empty and non-empty neighboring cells. The caller is responsible for freeing
 *         this string. Returns NULL in case of memory allocation failure.
 */
char *statGrid(GridSubstrat gridSubstrat);

/**
 * @brief Frees the memory allocated for a GridSubstrat.
 *
 * This function deallocates the memory used by the grid within the GridSubstrat structure.
 * It ensures that there are no memory leaks in the program when a GridSubstrat is no longer needed.
 *
 * @param gridSubstrat The GridSubstrat structure whose grid memory will be freed.
 */
void freeGridSubstrat(GridSubstrat gridSubstrat);

/**
 * @brief Writes the contents of a GridSubstrat to a specified file.
 *
 * This function creates a directory structure for saving results and writes the properties
 * and data of the GridSubstrat to a file in a human-readable format. The output file contains
 * details about the substrate dimensions, coordinates, and atom positions in each grid cell.
 *
 * @param filename The name of the output file (without extension) where the grid data will be saved.
 * @param gridSubstrat A pointer to the GridSubstrat structure that contains the data to be written.
 *                      This structure should have been properly initialized and populated before
 *                      calling this function.
 * @return void This function does not return a value. It will terminate the program if there
 *         is an error opening the output file.
 */
void writeGridSubstratToFile(char *filename, GridSubstrat *gridSubstrat);

/**
 * @brief Imports substrate data from a specified file and initializes a GridSubstrat.
 *
 * This function reads substrate data from a specified input file (e.g., in MOL2 format) and
 * initializes a new GridSubstrat structure based on the read data. The function assumes the
 * input file contains valid substrate data and uses a predefined grid step size.
 *
 * @param inputname The name of the input file containing substrate data.
 * @param substrat_t A pointer on table of positions of substrat's atoms
 * @return A GridSubstrat structure that has been initialized with the imported data.
 *         The returned structure must be properly managed (freed) by the user to avoid memory leaks.
 */
GridSubstrat importSubstratToGrid(char *inputname, double ***substrat_t);

#endif