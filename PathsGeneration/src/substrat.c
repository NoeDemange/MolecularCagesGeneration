/**
 * @file substrat.c
 * @brief Contains functions related to molecule substrate handling.
 *
 * This file implements functionality for reading and processing
 * molecule data from .mol2 files.
 * Various helper functions are also defined to perform substrate-related
 * operations.
 */

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include "constant.h"
#include "distance.h"
#include "substrat.h"
#include "util.h"

/**************************************/
/* INITIALIZATION OF MOLECULE *********/
/**************************************/

/**
 * @brief Retrieves the molecule data from a file with .mol2 extension.
 *
 * This function reads the molecule file and parses the data to create
 * an array of molecules. The file path is constructed by concatenating
 * the input name with a prefix and suffix. The number of atoms is also
 * retrieved.
 *
 * @param inputname The name of the file containing the molecule data.
 * @param atom_count Pointer to store the number of atoms.
 * @return A pointer to the array of molecules.
 */
double **readSubstratMol2(char *inputname, int *atom_count) {
  int ret;
  // Calculate the length of the fixed parts and the inputname
  const char *suffix = ".mol2";
  int suffix_len = strlen(suffix);
  int inputname_len = strlen(inputname);
  char *name = getBasename(inputname);
  // Calculate the total length needed for the filepath
  int filepath_len;
  if (inputname[inputname_len - 1] == '/') {
    filepath_len = inputname_len + strlen(name) + suffix_len + 1;
  } else {
    filepath_len = inputname_len + 1 + strlen(name) + suffix_len + 1;
  }

  // Allocate memory for the filepath
  char *filepath = (char *)malloc(filepath_len);
  if (filepath == NULL) {
    printf("Error allocating memory");
    exit(EXIT_FAILURE);
  }

  // Construct the file path
  if (inputname[inputname_len - 1] == '/') {
    sprintf(filepath, "%s%s%s", inputname, name, suffix);
  } else {
    sprintf(filepath, "%s/%s%s", inputname, name, suffix);
  }

  // Open the file
  FILE *filestream = fopen(filepath, "r");
  if (!filestream) {
    fprintf(stderr, "The file %s doesn't exist.\n", filepath);
    exit(EXIT_FAILURE);
  }

  free(filepath);
  free(name);

  // Define variable use for reading
  char line[256];

  // Read to number of Atom (jump 2 first lines and read on third)
  if (fgets(line, sizeof(line), filestream) == NULL) {
    printf("Error reading file");
    fclose(filestream);
    exit(EXIT_FAILURE);
  }
  if (fgets(line, sizeof(line), filestream) == NULL) {
    printf("Error reading file");
    fclose(filestream);
    exit(EXIT_FAILURE);
  }

  ret = fscanf(filestream, "%d", atom_count);

  // Jump next 4 lines
  for (int i = 0; i < 5; i++) {
    if (fgets(line, sizeof(line), filestream) == NULL) {
      printf("Error reading file");
      fclose(filestream);
      exit(EXIT_FAILURE);
    }
  }

  // Init double** coords array and variable to read
  double **substrat = (double **)malloc((*atom_count) * sizeof(double *));
  if (substrat == NULL) {
    printf("Error allocating memory");
    exit(EXIT_FAILURE);
  }
  for (int i = 0; i < (*atom_count); i++) {
    substrat[i] = (double *)malloc(3 * sizeof(double));
    if (substrat[i] == NULL) {
      printf("Error allocating memory");
      exit(EXIT_FAILURE);
    }
  }
  int index;
  char atom_name[10];

  for (int i = 0; i < (*atom_count); i++) {
    ret =
        fscanf(filestream, "%d %s %lf %lf %lf %s", &index, atom_name, &substrat[i][0], &substrat[i][1], &substrat[i][2],
               atom_name); // x=0,y=1,z=2
    // Consume the newline character left in the input buffer
    int c;
    while ((c = fgetc(filestream)) != '\n' && c != EOF)
      ;
  }

  fclose(filestream);

  if (ret < 0) {
    fprintf(stderr, "An error occured while reading %s.\n", inputname);
    exit(EXIT_FAILURE);
  }

  return substrat;
}

/**************************************/
/* Discretization             *********/
/**************************************/

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
GridSubstrat initGridSubstrat(double **substrat_t, int substrat_size, double step) {

  double x_min = INFINITY;
  double y_min = INFINITY;
  double z_min = INFINITY;
  double x_max = -INFINITY;
  double y_max = -INFINITY;
  double z_max = -INFINITY;

  for (int i = 0; i < substrat_size; i++) {
    if (substrat_t[i][0] < x_min)
      x_min = substrat_t[i][0];
    if (substrat_t[i][1] < y_min)
      y_min = substrat_t[i][1];
    if (substrat_t[i][2] < z_min)
      z_min = substrat_t[i][2];
    if (substrat_t[i][0] > x_max)
      x_max = substrat_t[i][0];
    if (substrat_t[i][1] > y_max)
      y_max = substrat_t[i][1];
    if (substrat_t[i][2] > z_max)
      z_max = substrat_t[i][2];
  }
  // printf("max/min, x:%lf/%lf y:%lf/%lf z:%lf/%lf\n",x_max,x_min,y_max,yMin,z_max,z_min);

  GridSubstrat grid_substrat;
  grid_substrat.substratSize = substrat_size;
  grid_substrat.xMax = x_max + DIST_GAP_SUBSTRATE;
  grid_substrat.xMin = x_min - DIST_GAP_SUBSTRATE;
  grid_substrat.yMax = y_max + DIST_GAP_SUBSTRATE;
  grid_substrat.yMin = y_min - DIST_GAP_SUBSTRATE;
  grid_substrat.zMax = z_max + DIST_GAP_SUBSTRATE;
  grid_substrat.zMin = z_min - DIST_GAP_SUBSTRATE;
  grid_substrat.xMaxWOGap = x_max;
  grid_substrat.xMinWOGap = x_min;
  grid_substrat.yMaxWOGap = y_max;
  grid_substrat.yMinWOGap = y_min;
  grid_substrat.zMaxWOGap = z_max;
  grid_substrat.zMinWOGap = z_min;

  grid_substrat.xSize = (int)ceil((grid_substrat.xMax - grid_substrat.xMin) / step);
  grid_substrat.ySize = (int)ceil((grid_substrat.yMax - grid_substrat.yMin) / step);
  grid_substrat.zSize = (int)ceil((grid_substrat.zMax - grid_substrat.zMin) / step);

  // printf("max/min, x:%lf/%lf y:%lf/%lf
  // z:%lf/%lf\n",grid_substrat.x_max,grid_substrat.x_min,grid_substrat.y_max,grid_substrat.yMin,grid_substrat.z_max,grid_substrat.z_min);
  // printf("ecart, x:%lf y:%lf
  // z:%lf\n",grid_substrat.x_max-grid_substrat.x_min,grid_substrat.y_max-grid_substrat.yMin,grid_substrat.z_max-grid_substrat.z_min);
  // printf("size, x:%d y:%d z:%d\n",grid_substrat.xSize,grid_substrat.ySize,grid_substrat.zSize);

  // Nb points max par case
  int *max_points = (int *)calloc(grid_substrat.xSize * grid_substrat.ySize * grid_substrat.zSize, sizeof(int));
  if (max_points == NULL) {
    printf("Error allocating memory\n");
    exit(EXIT_FAILURE);
  }
  int x_pos, y_pos, z_pos;
  double x_temp, y_temp, z_temp;
  for (int i = 0; i < substrat_size; i++) {
    x_pos = (int)((substrat_t[i][0] - grid_substrat.xMin) / step);
    y_pos = (int)((substrat_t[i][1] - grid_substrat.yMin) / step);
    z_pos = (int)((substrat_t[i][2] - grid_substrat.zMin) / step);
    // printf("Atom %d, x:%d y:%d z:%d\n",i,x_pos,y_pos,z_pos);
    max_points[(((x_pos) * (grid_substrat).ySize + (y_pos)) * (grid_substrat).zSize + (z_pos))]++;
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        for (int l = 0; l < 3; l++) {
          if ((x_pos - 1 + j) < 0 || (y_pos - 1 + k) < 0 || (z_pos - 1 + l) < 0 ||
              (x_pos - 1 + j) >= grid_substrat.xSize || (y_pos - 1 + k) >= grid_substrat.ySize ||
              (z_pos - 1 + l) >= grid_substrat.zSize || (j == 1 && k == 1 && l == 1))
            continue;
          if (j == 0) {
            x_temp = grid_substrat.xMin + (step * x_pos);
          } else if (j == 1) {
            x_temp = substrat_t[i][0];
          } else {
            x_temp = grid_substrat.xMin + (step * (x_pos + 1));
          }
          if (k == 0) {
            y_temp = grid_substrat.yMin + (step * y_pos);
          } else if (k == 1) {
            y_temp = substrat_t[i][1];
          } else {
            y_temp = grid_substrat.yMin + (step * (y_pos + 1));
          }
          if (l == 0) {
            z_temp = grid_substrat.zMin + (step * z_pos);
          } else if (l == 1) {
            z_temp = substrat_t[i][2];
          } else {
            z_temp = grid_substrat.zMin + (step * (z_pos + 1));
          }
          if (squaredEuclideanDistanceCoordsPoint(x_temp, y_temp, z_temp, substrat_t[i]) <
              DIST_GAP_SUBSTRATE * DIST_GAP_SUBSTRATE) {
            max_points[(((x_pos - 1 + j) * (grid_substrat).ySize + (y_pos - 1 + k)) * (grid_substrat).zSize +
                        (z_pos - 1 + l))]++;
          }
        }
      }
    }
  }

  int max_points_val = 0;
  for (int j = 0; j < grid_substrat.xSize; j++) {
    for (int k = 0; k < grid_substrat.ySize; k++) {
      for (int l = 0; l < grid_substrat.zSize; l++) {
        if (max_points_val < max_points[(((j) * (grid_substrat).ySize + (k)) * (grid_substrat).zSize + (l))]) {
          max_points_val = max_points[(((j) * (grid_substrat).ySize + (k)) * (grid_substrat).zSize + (l))];
        }
      }
    }
  }
  grid_substrat.caseSize = max_points_val;
  // printf("max %d",grid_substrat.max_points);
  free(max_points);
  // Allocation
  grid_substrat.gridS = (double *)malloc(grid_substrat.xSize * grid_substrat.ySize * grid_substrat.zSize *
                                         (grid_substrat.caseSize + 1) * 3 * sizeof(double));
  if (grid_substrat.gridS == NULL) {
    printf("Error allocating memory\n");
    exit(EXIT_FAILURE);
  }
  // Initialize the number of points to 0
  for (int i = 0; i < grid_substrat.xSize; i++) {
    for (int j = 0; j < grid_substrat.ySize; j++) {
      for (int k = 0; k < grid_substrat.zSize; k++) {
        grid_substrat.gridS[indexCoord(grid_substrat, i, j, k, 0, 0)] = 0;
      }
    }
  }

  for (int i = 0; i < substrat_size; i++) {
    x_pos = (int)((substrat_t[i][0] - grid_substrat.xMin) / step);
    y_pos = (int)((substrat_t[i][1] - grid_substrat.yMin) / step);
    z_pos = (int)((substrat_t[i][2] - grid_substrat.zMin) / step);
    // printf("Atom %d, x:%d y:%d z:%d\n",i,x_pos,y_pos,z_pos);

    getNumPointsF(grid_substrat, x_pos, y_pos, z_pos)++;
    memcpy(&grid_substrat.gridS[indexPoint(grid_substrat, x_pos, y_pos, z_pos,
                                           (int)getNumPointsF(grid_substrat, x_pos, y_pos, z_pos))],
           substrat_t[i], 3 * sizeof(double));
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        for (int l = 0; l < 3; l++) {
          if ((x_pos - 1 + j) < 0 || (y_pos - 1 + k) < 0 || (z_pos - 1 + l) < 0 ||
              (x_pos - 1 + j) >= grid_substrat.xSize || (y_pos - 1 + k) >= grid_substrat.ySize ||
              (z_pos - 1 + l) >= grid_substrat.zSize || (j == 1 && k == 1 && l == 1))
            continue;
          if (j == 0) {
            x_temp = grid_substrat.xMin + (step * x_pos);
          } else if (j == 1) {
            x_temp = substrat_t[i][0];
          } else {
            x_temp = grid_substrat.xMin + (step * (x_pos + 1));
          }
          if (k == 0) {
            y_temp = grid_substrat.yMin + (step * y_pos);
          } else if (k == 1) {
            y_temp = substrat_t[i][1];
          } else {
            y_temp = grid_substrat.yMin + (step * (y_pos + 1));
          }
          if (l == 0) {
            z_temp = grid_substrat.zMin + (step * z_pos);
          } else if (l == 1) {
            z_temp = substrat_t[i][2];
          } else {
            z_temp = grid_substrat.zMin + (step * (z_pos + 1));
          }

          if (squaredEuclideanDistanceCoordsPoint(x_temp, y_temp, z_temp, substrat_t[i]) <
              DIST_GAP_SUBSTRATE * DIST_GAP_SUBSTRATE) {
            getNumPointsF(grid_substrat, (x_pos - 1 + j), (y_pos - 1 + k), (z_pos - 1 + l))++;
            memcpy(&grid_substrat.gridS[indexPoint(
                       grid_substrat, (x_pos - 1 + j), (y_pos - 1 + k), (z_pos - 1 + l),
                       (int)getNumPointsF(grid_substrat, (x_pos - 1 + j), (y_pos - 1 + k), (z_pos - 1 + l)))],
                   substrat_t[i], 3 * sizeof(double));
          }
        }
      }
    }
  }

  return grid_substrat;
}

/**
 * @brief Checks for collision of a point with a substrate grid.
 *
 * This function determines if a given point is within a defined grid substrate.
 * It calculates the grid cell corresponding to the point and checks for collisions
 * with any points in that cell.
 *
 * @param point A 3-element array representing the coordinates of the point to check.
 *              It should be in the form [x, y, z].
 * @param grid_substrat A pointer to the GridSubstrat structure that contains information
 *                     about the grid and its properties, such as dimensions and point data.
 * @param step The size of each grid cell. It determines how the 3D point maps to the grid.
 * @return Returns 1 if a collision is detected (i.e., the point is close to any points
 *         in the grid cell), or 0 if no collision is detected (i.e., the point is
 *         outside the grid or there are no nearby points).
 */
int checkGridCollisionSubstratFloatTable(double point[3], GridSubstrat *grid_substrat, double step) {

  int x_pos = (int)((point[0] - grid_substrat->xMin) / step);
  int y_pos = (int)((point[1] - grid_substrat->yMin) / step);
  int z_pos = (int)((point[2] - grid_substrat->zMin) / step);
  if (x_pos < 0 || x_pos >= grid_substrat->xSize || y_pos < 0 || y_pos >= grid_substrat->ySize || z_pos < 0 ||
      z_pos >= grid_substrat->zSize)
    return 0;
  double *address = getNumPointsAddre(*grid_substrat, x_pos, y_pos, z_pos);
  double *max_address = address + (3 * ((int)*address));
  for (double *i = address + 3; i <= max_address; i += 3) { // itération sur les adresses
    if (squaredEuclideanDistance(point, i) < (DIST_GAP_SUBSTRATE * DIST_GAP_SUBSTRATE))
      return 1;
  }
  return 0;
}

/**
 * @brief Checks for collision of a point with a substrate grid.
 *
 * This function determines if a given point is within a defined grid substrate.
 * It calculates the grid cell corresponding to the point and checks for collisions
 * with any points in that cell.
 *
 * @param point A Point_t structure.
 * @param grid_substrat A pointer to the GridSubstrat structure that contains information
 *                     about the grid and its properties, such as dimensions and point data.
 * @param step The size of each grid cell. It determines how the 3D point maps to the grid.
 * @return Returns 1 if a collision is detected (i.e., the point is close to any points
 *         in the grid cell), or 0 if no collision is detected (i.e., the point is
 *         outside the grid or there are no nearby points).
 */
int checkGridCollisionSubstratPointT(Point_t p, GridSubstrat *grid_substrat, double step) {

  int x_pos = (int)((p.x - grid_substrat->xMin) / step);
  int y_pos = (int)((p.y - grid_substrat->yMin) / step);
  int z_pos = (int)((p.z - grid_substrat->zMin) / step);
  if (x_pos < 0 || x_pos >= grid_substrat->xSize || y_pos < 0 || y_pos >= grid_substrat->ySize || z_pos < 0 ||
      z_pos >= grid_substrat->zSize)
    return 0;
  double *address = getNumPointsAddre(*grid_substrat, x_pos, y_pos, z_pos);
  double *max_address = address + (3 * ((int)*address));
  for (double *i = address + 3; i <= max_address; i += 3) { // itération sur les adresses
    if (squaredEuclideanDistancefloatPointT(p, i) < (DIST_GAP_SUBSTRATE * DIST_GAP_SUBSTRATE))
      return 1;
  }
  return 0;
}

/**
 * @brief Generates statistics about the points in the substrate grid.
 *
 * This function computes several statistics for a given GridSubstrat. It evaluates
 * the maximum, minimum, and average number of points per grid cell, as well as the
 * number of cells that are adjacent to empty cells and those that are surrounded by
 * non-empty cells.
 *
 * @param grid_substrat The GridSubstrat structure that contains information about the grid
 *                     and its point data.
 * @return A dynamically allocated string containing the statistics of the grid. The
 *         string format includes the size of the grid, memory allocation size,
 *         maximum and minimum points per cell, and average points per cell for both
 *         empty and non-empty neighboring cells. The caller is responsible for freeing
 *         this string. Returns NULL in case of memory allocation failure.
 */
char *statGrid(GridSubstrat grid_substrat) {
  int max_nb_points = INT_MIN;
  int min_nb_points = INT_MAX;
  int sum_nb_points = 0;
  int sum_nb_points_void_box_around = 0;
  int case_nb_points_void_box_around = 0;
  int sum_nb_points_no_void_box_around = 0;
  int case_nb_points_no_void_box_around = 0;

  for (int x_pos = 0; x_pos < grid_substrat.xSize; x_pos++) {
    for (int y_pos = 0; y_pos < grid_substrat.ySize; y_pos++) {
      for (int z_pos = 0; z_pos < grid_substrat.zSize; z_pos++) {
        int taille = (int)getNumPointsF(grid_substrat, x_pos, y_pos, z_pos);
        if (taille > max_nb_points)
          max_nb_points = taille;
        if (taille < min_nb_points)
          min_nb_points = taille;
        sum_nb_points += taille;

        // void
        if ((x_pos - 1) < 0 || (y_pos - 1) < 0 || (z_pos - 1) < 0 || (x_pos + 1) >= grid_substrat.xSize ||
            (y_pos + 1) >= grid_substrat.ySize || (z_pos + 1) >= grid_substrat.zSize ||
            (int)getNumPointsF(grid_substrat, x_pos - 1, y_pos - 1, z_pos - 1) == 0 ||
            (int)getNumPointsF(grid_substrat, x_pos - 1, y_pos - 1, z_pos) == 0 ||
            (int)getNumPointsF(grid_substrat, x_pos - 1, y_pos - 1, z_pos + 1) == 0 ||
            (int)getNumPointsF(grid_substrat, x_pos - 1, y_pos, z_pos - 1) == 0 ||
            (int)getNumPointsF(grid_substrat, x_pos - 1, y_pos + 1, z_pos - 1) == 0 ||
            (int)getNumPointsF(grid_substrat, x_pos, y_pos - 1, z_pos - 1) == 0 ||
            (int)getNumPointsF(grid_substrat, x_pos + 1, y_pos - 1, z_pos - 1) == 0 ||
            (int)getNumPointsF(grid_substrat, x_pos, y_pos - 1, z_pos) == 0 ||
            (int)getNumPointsF(grid_substrat, x_pos, y_pos - 1, z_pos + 1) == 0 ||
            (int)getNumPointsF(grid_substrat, x_pos, y_pos, z_pos - 1) == 0 ||
            (int)getNumPointsF(grid_substrat, x_pos, y_pos + 1, z_pos - 1) == 0 ||
            (int)getNumPointsF(grid_substrat, x_pos + 1, y_pos - 1, z_pos) == 0 ||
            (int)getNumPointsF(grid_substrat, x_pos + 1, y_pos - 1, z_pos + 1) == 0 ||
            (int)getNumPointsF(grid_substrat, x_pos + 1, y_pos, z_pos - 1) == 0 ||
            (int)getNumPointsF(grid_substrat, x_pos + 1, y_pos + 1, z_pos - 1) == 0 ||
            (int)getNumPointsF(grid_substrat, x_pos, y_pos, z_pos + 1) == 0 ||
            (int)getNumPointsF(grid_substrat, x_pos, y_pos + 1, z_pos) == 0 ||
            (int)getNumPointsF(grid_substrat, x_pos, y_pos + 1, z_pos + 1) == 0 ||
            (int)getNumPointsF(grid_substrat, x_pos - 1, y_pos + 1, z_pos) == 0 ||
            (int)getNumPointsF(grid_substrat, x_pos - 1, y_pos + 1, z_pos + 1) == 0 ||
            (int)getNumPointsF(grid_substrat, x_pos - 1, y_pos, z_pos) == 0 ||
            (int)getNumPointsF(grid_substrat, x_pos - 1, y_pos, z_pos + 1) == 0 ||
            (int)getNumPointsF(grid_substrat, x_pos + 1, y_pos, z_pos) == 0 ||
            (int)getNumPointsF(grid_substrat, x_pos + 1, y_pos, z_pos + 1) == 0 ||
            (int)getNumPointsF(grid_substrat, x_pos + 1, y_pos + 1, z_pos) == 0 ||
            (int)getNumPointsF(grid_substrat, x_pos + 1, y_pos + 1, z_pos + 1) == 0) {
          sum_nb_points_void_box_around += (int)getNumPointsF(grid_substrat, x_pos, y_pos, z_pos);
          case_nb_points_void_box_around++;
        }

        // Novoid
        if ((x_pos - 1) >= 0 && (y_pos - 1) >= 0 && (z_pos - 1) >= 0 && (x_pos + 1) < grid_substrat.xSize &&
            (y_pos + 1) < grid_substrat.ySize && (z_pos + 1) < grid_substrat.zSize &&
            (int)getNumPointsF(grid_substrat, x_pos - 1, y_pos - 1, z_pos - 1) != 0 &&
            (int)getNumPointsF(grid_substrat, x_pos - 1, y_pos - 1, z_pos) != 0 &&
            (int)getNumPointsF(grid_substrat, x_pos - 1, y_pos - 1, z_pos + 1) != 0 &&
            (int)getNumPointsF(grid_substrat, x_pos - 1, y_pos, z_pos - 1) != 0 &&
            (int)getNumPointsF(grid_substrat, x_pos - 1, y_pos + 1, z_pos - 1) != 0 &&
            (int)getNumPointsF(grid_substrat, x_pos, y_pos - 1, z_pos - 1) != 0 &&
            (int)getNumPointsF(grid_substrat, x_pos + 1, y_pos - 1, z_pos - 1) != 0 &&
            (int)getNumPointsF(grid_substrat, x_pos, y_pos - 1, z_pos) != 0 &&
            (int)getNumPointsF(grid_substrat, x_pos, y_pos - 1, z_pos + 1) != 0 &&
            (int)getNumPointsF(grid_substrat, x_pos, y_pos, z_pos - 1) != 0 &&
            (int)getNumPointsF(grid_substrat, x_pos, y_pos + 1, z_pos - 1) != 0 &&
            (int)getNumPointsF(grid_substrat, x_pos + 1, y_pos - 1, z_pos) != 0 &&
            (int)getNumPointsF(grid_substrat, x_pos + 1, y_pos - 1, z_pos + 1) != 0 &&
            (int)getNumPointsF(grid_substrat, x_pos + 1, y_pos, z_pos - 1) != 0 &&
            (int)getNumPointsF(grid_substrat, x_pos + 1, y_pos + 1, z_pos - 1) != 0 &&
            (int)getNumPointsF(grid_substrat, x_pos, y_pos, z_pos + 1) != 0 &&
            (int)getNumPointsF(grid_substrat, x_pos, y_pos + 1, z_pos) != 0 &&
            (int)getNumPointsF(grid_substrat, x_pos, y_pos + 1, z_pos + 1) != 0 &&
            (int)getNumPointsF(grid_substrat, x_pos - 1, y_pos + 1, z_pos) != 0 &&
            (int)getNumPointsF(grid_substrat, x_pos - 1, y_pos + 1, z_pos + 1) != 0 &&
            (int)getNumPointsF(grid_substrat, x_pos - 1, y_pos, z_pos) != 0 &&
            (int)getNumPointsF(grid_substrat, x_pos - 1, y_pos, z_pos + 1) != 0 &&
            (int)getNumPointsF(grid_substrat, x_pos + 1, y_pos, z_pos) != 0 &&
            (int)getNumPointsF(grid_substrat, x_pos + 1, y_pos, z_pos + 1) != 0 &&
            (int)getNumPointsF(grid_substrat, x_pos + 1, y_pos + 1, z_pos) != 0 &&
            (int)getNumPointsF(grid_substrat, x_pos + 1, y_pos + 1, z_pos + 1) != 0) {
          sum_nb_points_no_void_box_around += (int)getNumPointsF(grid_substrat, x_pos, y_pos, z_pos);
          case_nb_points_no_void_box_around++;
        }
      }
    }
  }
  size_t taille_mem = grid_substrat.xSize * grid_substrat.ySize * grid_substrat.zSize * (grid_substrat.caseSize + 1) *
                      3 * sizeof(double);

  // Allocate memory for the result string
  char *result = (char *)malloc(512 * sizeof(char));
  if (result == NULL) {
    // Handle memory allocation error
    return NULL;
  }

  // Format the statistics into the result string
  snprintf(
      result, 512,
      "STATS :\n GRID : xSize %d, ySize %d, zSize %d, memoire %ld octets\n Points/Case : Max %d, Min %d, Moyenne "
      "%lf\n Case vide voisine? : nb %d, MoyPoints/case %lf\n NON Case vide voisine? : nb %d, MoyPoints/case %lf\n",
      grid_substrat.xSize, grid_substrat.ySize, grid_substrat.zSize, taille_mem, max_nb_points, min_nb_points,
      sum_nb_points / (double)(grid_substrat.xSize * grid_substrat.ySize * grid_substrat.zSize),
      case_nb_points_no_void_box_around, sum_nb_points_void_box_around / (double)case_nb_points_void_box_around,
      case_nb_points_no_void_box_around, sum_nb_points_no_void_box_around / (double)case_nb_points_no_void_box_around);
  return result;
  // printf("STATS :\n GRID : xSize %d, ySize %d, zSize %d, mémoire %ld octets\n Points/Case : Max %d, Min %d, Moyenne
  // %lf\n Case vide voisine? : nb %d, MoyPoints/case %lf\n NON Case vide voisine? : nb %d, MoyPoints/case
  // %lf\n",grid_substrat.xSize,grid_substrat.ySize,grid_substrat.zSize,taille_mem,max_nb_points,min_nb_points,sum_nb_points/(double)(grid_substrat.xSize*grid_substrat.ySize*grid_substrat.zSize),case_nb_points_no_void_box_around,um_nb_points_void_box_around/(double)case_nb_points_no_void_box_around,case_nb_points_no_void_box_around,um_nb_points_no_void_box_around/(double)case_nb_points_no_void_box_around);
}

/**
 * @brief Frees the memory allocated for a GridSubstrat.
 *
 * This function deallocates the memory used by the grid within the GridSubstrat structure.
 * It ensures that there are no memory leaks in the program when a GridSubstrat is no longer needed.
 *
 * @param grid_substrat The GridSubstrat structure whose grid memory will be freed.
 */
void freeGridSubstrat(GridSubstrat grid_substrat) { free(grid_substrat.gridS); }

/**
 * @brief Writes the contents of a GridSubstrat to a specified file.
 *
 * This function creates a directory structure for saving results and writes the properties
 * and data of the GridSubstrat to a file in a human-readable format. The output file contains
 * details about the substrate dimensions, coordinates, and atom positions in each grid cell.
 *
 * @param filename The name of the output file (without extension) where the grid data will be saved.
 * @param grid_substrat A pointer to the GridSubstrat structure that contains the data to be written.
 *                      This structure should have been properly initialized and populated before
 *                      calling this function.
 * @return void This function does not return a value. It will terminate the program if there
 *         is an error opening the output file.
 */
void writeGridSubstratToFile(char *filename, GridSubstrat *grid_substrat) {
  mkdir("./results", 0755);
  mkdir("./results/Write", 0755);
  char outputname[512];
  sprintf(outputname, "./results/Write/%s.txt", filename);
  FILE *file = fopen(outputname, "w");
  if (file == NULL) {
    printf("Error opening file %s for writing\n", outputname);
    exit(EXIT_FAILURE);
  }

  // Write the scalar fields with commentary
  fprintf(file, "substratSize: %d\n", grid_substrat->substratSize);
  fprintf(file, "x_min: %lf    x_max: %lf\n", grid_substrat->xMin, grid_substrat->xMax);
  fprintf(file, "yMin: %lf    y_max: %lf\n", grid_substrat->yMin, grid_substrat->yMax);
  fprintf(file, "z_min: %lf    z_max: %lf\n", grid_substrat->zMin, grid_substrat->zMax);
  fprintf(file, "xSize: %d\n", grid_substrat->xSize);
  fprintf(file, "ySize: %d\n", grid_substrat->ySize);
  fprintf(file, "zSize: %d\n\n", grid_substrat->zSize);

  // Write the gridS data with commentary
  for (int i = 0; i < grid_substrat->xSize; i++) {
    for (int j = 0; j < grid_substrat->ySize; j++) {
      for (int k = 0; k < grid_substrat->zSize; k++) {
        fprintf(file, "Case [%d][%d][%d], Nombre d'atomes : %d\n", i, j, k,
                (int)getNumPointsF(*grid_substrat, i, j, k));
        if ((int)getNumPointsF(*grid_substrat, i, j, k) == 0)
          continue;
        fprintf(file, "  Les atomes sont : \n");
        for (int points = 1; points <= (int)getNumPointsF(*grid_substrat, i, j, k); points++) {
          fprintf(file, "    [%d] : ", points);
          for (int coords = 0; coords < 3; coords++) { // Assuming 3 coordinates
            fprintf(file, "%lf ", grid_substrat->gridS[indexCoord(*grid_substrat, i, j, k, points, coords)]);
          }
          fprintf(file, "\n");
        }
      }
    }
  }

  fclose(file);
}

/**
 * @brief Imports substrate data from a specified file and initializes a GridSubstrat.
 *
 * This function reads substrate data from a specified input file (e.g., in MOL2 format) and
 * initializes a new GridSubstrat structure based on the read data. The function assumes the
 * input file contains valid substrate data and uses a predefined grid step size.
 *
 * @param inputname The name of the input file containing substrate data.
 * @param substrat_t A pointer on table of positions of substrat's atoms.
 * @return A GridSubstrat structure that has been initialized with the imported data.
 *         The returned structure must be properly managed (freed) by the user to avoid memory leaks.
 */
GridSubstrat importSubstratToGrid(char *inputname, double ***substrat_t) {
  int substrat_size = 0;
  *substrat_t = readSubstratMol2(inputname, &substrat_size);
  return initGridSubstrat(*substrat_t, substrat_size, STEP_GRID);
}