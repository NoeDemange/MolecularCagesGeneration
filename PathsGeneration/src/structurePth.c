#include "discretization.h"
#include "distance.h"
#include "structure.h"
//#include "util.h"

/**
 * @file structurePth.c
 * @brief Grouping of functions used for manipulating Path_t.
 *
 * This file contains functions for creating, manipulating, and printing Path_t element.
 */

/**
 * @brief create a new Paths_t.
 *
 * @param size Maximum number of patterns in a path.
 * @param numComponents Number of components.
 * @return (Path_t*) Pointer to the new paths.
 */
Paths_t *pthCreate(int size, int numComponents) {
  Paths_t *paths = malloc(sizeof(Paths_t));
  // paths->positionsBuffer = malloc(MAX_POSITION_KEEP_360 * sizeof(Point_t)); // Unused for the moment. TODO use
  // it to store positions computed in projectionAX1E3.
  paths->sizeMax = size + 4; // 2 additional columns for starting atom, its neighbor and end and hydrogen of end. think
                             // to decrease by two for iteration because 0 is start and n is end.
  paths->numPaths = numComponents - 1; // Number of part in the partition -1 = number of paths to create.
  // int sizePosition = NUMBER_POSITION_AX1E3;//(NUMBER_POSITION_AX1E3 < 2 ? 2 : NUMBER_POSITION_AX1E3);
  paths->patterns = (Point_t *)malloc(paths->numPaths * paths->sizeMax * NUMBER_POSITION_PATHS * MAX_NB_ATOMS_PATTERN *
                                      sizeof(Point_t));
  paths->positionCurNum = (int *)calloc(paths->numPaths * paths->sizeMax, sizeof(int));
  paths->patternCurNum = (int *)calloc(paths->numPaths * paths->sizeMax, sizeof(int));
  paths->maxPositions = (int *)malloc(paths->numPaths * paths->sizeMax * sizeof(int));
  paths->curPthPos = (int *)malloc(paths->numPaths * sizeof(int)); // one int for each path.
  paths->pathRMSD = (double *)malloc(paths->numPaths * sizeof(double)); // one double for each path.
  paths->bestPathLength = (int *)malloc(paths->numPaths * sizeof(int));
  paths->maxGrowthLimit = (int *)malloc(paths->numPaths * sizeof(int));

  for (int k = 0; k < paths->numPaths; k++) {
    paths->bestPathLength[k] = -1;
    paths->maxGrowthLimit[k] = paths->sizeMax - 3;
    for (int s = 0; s < paths->sizeMax; s++) {
      for (int p = 0; p < NUMBER_POSITION_PATHS; p++) {
        for (int a = 0; a < MAX_NB_ATOMS_PATTERN; a++) {
          paths->patterns[indexPointPaths(k, s, p, a, paths->sizeMax)] = ptInit(0);
        }
      }
      paths->maxPositions[indexPathPosition(k, s, paths->sizeMax)] = -1;
    }
  }
  // 	for (int i = 0; i < sizeMax; i++) {
  // 		paths->patterns[i] = malloc(sizePosition * sizeof(Point_t*));
  //     for (int j = 0; j < sizePosition; j++) {
  //       paths->patterns[i][j] = malloc(MAX_NB_ATOMS_PATTERN * sizeof(Point_t));
  //       for (int k = 0; k <  MAX_NB_ATOMS_PATTERN; k++) {
  //         paths->patterns[i][j][k] = PT_init(0);
  //       }
  //     }
  //     paths->positionNum[i] = 0;
  //     paths->patternNum[i] = 0;
  //     paths->maxPositions[i] = -1;
  // 	}
  //   paths->idStart = pair->start;
  //   paths->idEnd = pair->end;
  //  paths->sizeMax = sizeMax - 1; // Because indexing starts at 0. Pas sûr de besoin juste faire attention indice

  if (get_current_distance_type() != DISTANCE_EUCLIDEAN) {
    paths->grids = (Grid_t **)malloc((paths->numPaths) *
                                     sizeof(Grid_t *)); // numpaths-1+1 for the initial grid and not the last path
    for (int i = 0; i < paths->numPaths; i++) {
      paths->grids[i] = malloc(sizeof(Grid_t)); // Alloc one time and after just change the values.
      paths->grids[i]->nodes = NULL;
      paths->grids[i]->visitedNodes = NULL;
    }

    paths->minHeaps = (MinHeap_t **)malloc((paths->numPaths) * sizeof(MinHeap_t *));
    for (int i = 0; i < paths->numPaths; i++) {
      paths->minHeaps[i] = malloc(sizeof(MinHeap_t));
      paths->minHeaps[i]->nodes = NULL;
    }

    if (get_current_distance_type() == DISTANCE_SSMTA_STAR || get_current_distance_type() == DISTANCE_HYBRID) {
      paths->starts = (Point_t *)malloc(MAX_POSITION_KEEP_360 * sizeof(Point_t));
      paths->results_pos = (Point_t *)malloc(MAX_POSITION_KEEP_360 * sizeof(Point_t));
      paths->Candidates = (Node **)malloc(MAX_POSITION_KEEP_360 * sizeof(Node *));
      paths->distancesMultiCandidates = (double *)malloc(MAX_POSITION_KEEP_360 * sizeof(double));
    }
  }

  paths->currentPath = 0;

  return paths;
}

/**
 * @brief Reboots the paths structure by resetting the current path position and other related arrays.
 *
 * This function resets the current path position for each path in the Paths_t structure, allowing for a fresh start
 * in path processing. It also resets the position and pattern numbers for each path position.
 *
 * @param paths Pointer to the Paths_t structure to be rebooted.
 */
void pthReboot(Paths_t *paths) {
  // Reset the current path position for each path
  for (int k = 0; k < paths->numPaths; k++) {
    paths->bestPathLength[k] = -1;
    paths->maxGrowthLimit[k] = paths->sizeMax - 3;
    for (int s = 0; s < paths->sizeMax; s++) {
      paths->maxPositions[indexPathPosition(k, s, paths->sizeMax)] = -1;
      paths->positionCurNum[indexPathPosition(k, s, paths->sizeMax)] = 0;
      paths->patternCurNum[indexPathPosition(k, s, paths->sizeMax)] = 0;
    }
  }
}

/**
 * @brief Deletes the paths and frees the associated memory.
 *
 * @param paths Pointer to the paths to be deleted.
 */
void pthDelete(Paths_t *paths) {

  if (paths) {

    free(paths->patterns);
    // free(path->positionsBuffer);
    free(paths->positionCurNum);
    free(paths->patternCurNum);
    free(paths->maxPositions);
    free(paths->curPthPos);
    free(paths->pathRMSD);
    free(paths->bestPathLength);
    free(paths->maxGrowthLimit);

    if (get_current_distance_type() != DISTANCE_EUCLIDEAN) {
      // Free the grids
      for (int i = 0; i < paths->numPaths; i++) {
        if (paths->grids[i]) {
          freeGrid(paths->grids[i]);
        }
      }
      free(paths->grids);

      // Free the minHeaps
      for (int i = 0; i < paths->numPaths; i++) {

        if (paths->minHeaps[i]->nodes) {
          free(paths->minHeaps[i]->nodes);
        }
        free(paths->minHeaps[i]);
      }
      free(paths->minHeaps);

      if (get_current_distance_type() == DISTANCE_SSMTA_STAR || get_current_distance_type() == DISTANCE_HYBRID) {
        free(paths->starts);
        free(paths->results_pos);
        free(paths->Candidates);
        free(paths->distancesMultiCandidates);
      }
    }
    free(paths);
  }
}

/**
 * @brief Initializes a `Paths_t` structure with an interconnection tree (starting points and neighbors).
 *
 * This function initializes a `Paths_t` structure by filling it with an interconnection tree. The tree is
 * constructed based on the `interTree` array.
 *
 * @param paths A pointer to the `Paths_t` structure that will be initialized with the interconnection tree.
 * @param interTree An array representing the interconnection tree. Each pair of integers in this array
 *                  corresponds to a start point and end point.
 * @param cage A pointer to the `Cage_t` structure which contains the atoms and coordinates. This structure
 *             is used to retrieve atom coordinates for the interconnection tree.
 *
 * @details
 * The function processes the `interTree` array and uses it to populate the `Paths_t` structure:
 * - The start points and their corresponding neighbor coordinates are retrieved and stored in the `patterns` array.
 * - The start points' current position in the path is initialized, and the maximum positions are set.
 * - The function assumes that for each path, there is at least one starting point and one neighbor to begin the path.
 *
 * The paths are initialized such that the `patterns` array holds the coordinates of the start point's
 * neighbor and the coordinates of the start point itself. The `maxPositions` array is also updated to
 * track the maximum positions of the current paths.
 *
 * @note The function sets the initial path position (`curPthPos`) to `2` for each path to indicate that the
 *       first two positions are filled (start point and its neighbor). This allows the path creation to
 *       continue from this point in subsequent operations.
 *
 */
void pthInit(Paths_t *paths, int *interTree, Cage_t *cage) {

  for (int k = 0; k < paths->numPaths; k++) {
    int start_id = interTree[2 * k];
    // printf("k %d , Start ID: %d, x %.2f, y %.2f, z %.2f\n",k, start_id, coords(atom(cage, start_id)).x
    // ,coords(atom(cage, start_id)).y, coords(atom(cage, start_id)).z );
    paths->patterns[indexPointPaths(k, 0, 0, 0, paths->sizeMax)] = coordsNeighbor(cage, start_id, 0);
    paths->patterns[indexPointPaths(k, 1, 0, 0, paths->sizeMax)] = coords(atom(cage, start_id));
    paths->maxPositions[indexPathPosition(k, 1, paths->sizeMax)] = 0;
    paths->curPthPos[k] = 2; // Because 1 is start and 0 start's neighbor init with PTH_init
    paths->bestPathLength[k] = -1;
    paths->maxGrowthLimit[k] = paths->sizeMax - 3;
    // for( int s = 0; s<paths->sizeMax; s++){
    //     paths->positionCurNum[indexPathPosition(k,s,paths->sizeMax)] = 0;
    //     paths->patternCurNum[indexPathPosition(k,s,paths->sizeMax)] = 0;
    // }
  }
  paths->currentPath = 0;
  // return paths;
}

/**
 * @brief Prints or writes the paths structure to a file or console with detailed alignment.
 *
 * This function prints the `Paths_t` structure to the console or writes it to a specified file.
 * The output is a neatly formatted table with paths organized by their components, positions, and atoms.
 *
 * @param paths Pointer to the `Paths_t` structure containing the paths data.
 * @param file Pointer to the open file for writing. If `NULL`, the data is printed to console.
 *
 * The format of the output table includes:
 * - A header row with indices of `sizeMax` components.
 * - Rows for each position (`Pos`) with corresponding atoms (`Atom`).
 * - Each cell contains the coordinates (`x, y, z`) of a point in the pattern.
 *
 * @note The function dynamically switches between printing to the console and writing to a file.
 *       If the file cannot be opened, an error message is displayed.
 *
 * @warning Ensure the `paths` structure is properly initialized before calling this function.
 *          Any invalid indices or uninitialized memory in the `patterns` array may cause undefined behavior.
 */
void pthPrintOrWritePaths(Paths_t *paths, FILE *file) {
#define PRINT(fmt, ...)                                                                                                \
  do {                                                                                                                 \
    if (file)                                                                                                          \
      fprintf(file, fmt, ##__VA_ARGS__);                                                                               \
    else                                                                                                               \
      printf(fmt, ##__VA_ARGS__);                                                                                      \
  } while (0)

  // Print the paths structure
  for (int k = 0; k < paths->numPaths; k++) {
    PRINT("\nPath %2d:", k);

    // Print the sizeMax indices above the table
    PRINT("         +"); // Offset for left labels
    for (int m = 0; m < paths->sizeMax; m++) {
      PRINT("         %4d         +", m);
    }
    PRINT("\n");

    // Top horizontal border
    PRINT("-----------------+");
    for (int m = 0; m < paths->sizeMax; m++) {
      PRINT("----------------------+");
    }
    PRINT("\n");

    for (int p = 0; p < NUMBER_POSITION_PATHS; p++) {
      for (int a = 0; a < MAX_NB_ATOMS_PATTERN; a++) {
        // Print Pos label only once, on the same row as Atom 0
        if (a == 0) {
          PRINT("Pos %2d | Atom %2d |", p, a);
        } else {
          PRINT("       | Atom %2d |", a);
        }

        for (int m = 0; m < paths->sizeMax; m++) {
          int index = indexPointPaths(k, m, p, a, paths->sizeMax);
          Point_t point = paths->patterns[index];

          // Print point values with one decimal precision, aligned
          PRINT(" %6.1f,%6.1f,%6.1f |", point.x, point.y, point.z);
        }
        PRINT("\n");
      }

      // Divider between rows
      PRINT("-----------------+");
      for (int m = 0; m < paths->sizeMax; m++) {
        PRINT("----------------------+");
      }
      PRINT("\n");
    }
  }

#undef PRINT
}

/**
 * @brief Prints or writes the path tables section of data to a file or console.
 *
 * This function prints or writes the path tables data. If writing to a file,
 * it appends the data after the existing content.
 *
 * @param paths Pointer to the `Paths_t` structure containing the data.
 * @param file Pointer to the open file for writing. If `NULL`, the data is printed to console.
 */
void pthPrintOrWritePathTables(Paths_t *paths, FILE *file) {

#define PRINT(fmt, ...)                                                                                                \
  do {                                                                                                                 \
    if (file)                                                                                                          \
      fprintf(file, fmt, ##__VA_ARGS__);                                                                               \
    else                                                                                                               \
      printf(fmt, ##__VA_ARGS__);                                                                                      \
  } while (0)

  // Print table headers
  PRINT("\nPath Tables:\n");

  for (int k = 0; k < paths->numPaths; k++) {
    // Print the path index
    PRINT("\nPath %2d:\n", k);

    // Print sizeMax indices above the table
    PRINT("         +");
    for (int s = 0; s < paths->sizeMax; s++) {
      PRINT("     %4d     +", s);
    }
    PRINT("\n");

    // Print horizontal border
    PRINT("---------+");
    for (int s = 0; s < paths->sizeMax; s++) {
      PRINT("--------------+");
    }
    PRINT("\n");

    // Print positionCurNum
    PRINT("CurPos   |");
    for (int s = 0; s < paths->sizeMax; s++) {
      int index = indexPathPosition(k, s, paths->sizeMax);
      PRINT("     %4d     |", paths->positionCurNum[index]);
    }
    PRINT("\n");

    // Print patternCurNum
    PRINT("CurPat   |");
    for (int s = 0; s < paths->sizeMax; s++) {
      int index = indexPathPosition(k, s, paths->sizeMax);
      PRINT("     %4d     |", paths->patternCurNum[index]);
    }
    PRINT("\n");

    // Print maxPositions
    PRINT("MaxPos   |");
    for (int s = 0; s < paths->sizeMax; s++) {
      int index = indexPathPosition(k, s, paths->sizeMax);
      PRINT("     %4d     |", paths->maxPositions[index]);
    }
    PRINT("\n");

    // Print horizontal border below each path
    PRINT("---------+");
    for (int s = 0; s < paths->sizeMax; s++) {
      PRINT("--------------+");
    }
    PRINT("\n");
  }

#undef PRINT
}

/**
 * @brief Prints or writes both the paths and path tables to a file or console.
 *
 * This function calls `PTH_printOrWritePaths` and `PTH_printOrWritePathTables`
 * to output the data stored in the `Paths_t` structure. The output can be directed
 * to the console or a specified file.
 *
 * @param paths Pointer to the `Paths_t` structure containing the data.
 * @param filename Pointer to a string specifying the output file name.
 *                 If `NULL`, the output is printed to the console.
 */
void pthPrintOrWriteAll(Paths_t *paths, const char *filename) {
  FILE *file = NULL;
  int to_file = 0;

  // Open the file if a filename is provided
  if (filename) {
    file = fopen(filename, "w");
    if (!file) {
      printf("Error opening file %s for writing.\n", filename);
      return;
    }
    to_file = 1;
  }

// Helper macro for printing to file or console
#define PRINT(fmt, ...)                                                                                                \
  do {                                                                                                                 \
    if (to_file)                                                                                                       \
      fprintf(file, fmt, ##__VA_ARGS__);                                                                               \
    else                                                                                                               \
      printf(fmt, ##__VA_ARGS__);                                                                                      \
  } while (0)

  // Print the Paths data section
  PRINT("=== Paths Data ===\n");
  pthPrintOrWritePaths(paths, file); // Pass the file pointer to print paths

  // Add extra newline after paths section
  PRINT("\n");

  // Print the Path Tables data section
  PRINT("=== Path Tables ===\n");
  pthPrintOrWritePathTables(paths, file); // Pass the file pointer to print path tables

  // Close the file if writing to a file
  if (to_file)
    fclose(file);

#undef PRINT
}
