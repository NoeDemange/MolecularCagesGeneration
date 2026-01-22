#include "output.h"
#include "assembly.h"
#include "constant.h"
#include "distance.h"
#include "interconnection.h"
#include "util.h"

#include <errno.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>

typedef struct {
  int initialized;
  int statsOnly;
  int resultsCount;
  double bestRmsd;
  double worstRmsd;
  int minPathLength;
  int maxPathLength;
  double sumPathLength;
  int pathSamples;
  double startMonotonicMs;
  double firstResultMs;
  int hasFirstResult;
  int minPathCount;
  double bestRmsdAtMinPath;
  int minAtomCount;
  int maxAtomCount;
  double sumAtomCount;
  int atomSamples;
  int minAtomCountHits;
  double bestRmsdAtMinAtom;
} ResultSummary;

static ResultSummary g_result_summary = {0};

/**
 * @brief Retrieve the current monotonic time in milliseconds.
 *
 * @return Current monotonic time relative to an arbitrary origin.
 */
static double monotonicMilliseconds(void) { return monotonic_now_ms(); }

#ifdef ENABLE_STATS
static void emitStatsSummary(FILE *stream);
#endif

/**
 * @brief Lazily initialize the in-memory statistics accumulator.
 *
 * This prepares the global ResultSummary structure the first time it is needed
 * and detects whether the process runs in stats-only mode (CAGE_DISABLE_WRITE).
 */
static void ensureResultSummaryInitialized(void) {
  if (g_result_summary.initialized)
    return;
  const char *env = getenv("CAGE_DISABLE_WRITE");
  g_result_summary.statsOnly = (env != NULL && env[0] != '\0' && env[0] != '0');
  g_result_summary.resultsCount = 0;
  g_result_summary.bestRmsd = DBL_MAX;
  g_result_summary.worstRmsd = -DBL_MAX;
  g_result_summary.minPathLength = INT_MAX;
  g_result_summary.maxPathLength = INT_MIN;
  g_result_summary.sumPathLength = 0.0;
  g_result_summary.pathSamples = 0;
  g_result_summary.startMonotonicMs = monotonicMilliseconds();
  g_result_summary.firstResultMs = 0.0;
  g_result_summary.hasFirstResult = 0;
  g_result_summary.minPathCount = 0;
  g_result_summary.bestRmsdAtMinPath = DBL_MAX;
  g_result_summary.minAtomCount = INT_MAX;
  g_result_summary.maxAtomCount = INT_MIN;
  g_result_summary.sumAtomCount = 0.0;
  g_result_summary.atomSamples = 0;
  g_result_summary.minAtomCountHits = 0;
  g_result_summary.bestRmsdAtMinAtom = DBL_MAX;
  g_result_summary.initialized = 1;
}

/**
 * @brief Tell whether the solver runs in stats-only mode.
 *
 * @return 1 when file generation is disabled, 0 otherwise.
 */
static int statsOnlyMode(void) {
  ensureResultSummaryInitialized();
  return g_result_summary.statsOnly;
}

/**
 * @brief Update aggregated metrics for one generated cage result.
 *
 * @param rmsd Root-mean-square deviation of the current result.
 * @param min_len Minimum path length observed within this result.
 * @param max_len Maximum path length observed within this result.
 * @param total_len Sum of all path lengths, used to compute the average.
 * @param samples Number of paths contributing to the length statistics.
 * @param atom_count Number of atoms added by all paths for this result.
 */
static void recordResultMetrics(double rmsd, int min_len, int max_len, int total_len, int samples, int atom_count) {
  if (!statsOnlyMode())
    return;

  ResultSummary *summary = &g_result_summary;
  if (samples == 0) {
    min_len = 0;
    max_len = 0;
  }
  if (atom_count < 0)
    atom_count = 0;

  if (summary->resultsCount == 0) {
    summary->bestRmsd = rmsd;
    summary->worstRmsd = rmsd;
  } else {
    if (rmsd < summary->bestRmsd)
      summary->bestRmsd = rmsd;
    if (rmsd > summary->worstRmsd)
      summary->worstRmsd = rmsd;
  }

  if (samples > 0) {
    if (summary->minPathLength == INT_MAX || min_len < summary->minPathLength) {
      summary->minPathLength = min_len;
      summary->minPathCount = 1;
      summary->bestRmsdAtMinPath = rmsd;
    } else if (min_len == summary->minPathLength) {
      summary->minPathCount++;
      if (rmsd < summary->bestRmsdAtMinPath)
        summary->bestRmsdAtMinPath = rmsd;
    }
    if (summary->maxPathLength == INT_MIN || max_len > summary->maxPathLength)
      summary->maxPathLength = max_len;
    summary->sumPathLength += total_len;
    summary->pathSamples += samples;
  }

  summary->resultsCount++;
  if (!summary->hasFirstResult) {
    summary->firstResultMs = monotonicMilliseconds() - summary->startMonotonicMs;
    summary->hasFirstResult = 1;
  }

  if (summary->minAtomCount == INT_MAX || atom_count < summary->minAtomCount) {
    summary->minAtomCount = atom_count;
    summary->minAtomCountHits = 1;
    summary->bestRmsdAtMinAtom = rmsd;
  } else if (atom_count == summary->minAtomCount) {
    summary->minAtomCountHits++;
    if (rmsd < summary->bestRmsdAtMinAtom)
      summary->bestRmsdAtMinAtom = rmsd;
  }

  if (summary->maxAtomCount == INT_MIN || atom_count > summary->maxAtomCount)
    summary->maxAtomCount = atom_count;
  summary->sumAtomCount += atom_count;
  summary->atomSamples++;
}

/**
 * @brief Print the consolidated RESULT_SUMMARY line for stats-only runs.
 *
 * @param total_wall_ms Total elapsed monotonic time in milliseconds.
 * @param total_cpu_ms Total CPU time consumed in milliseconds.
 */
static void printResultSummary(double total_wall_ms, double total_cpu_ms) {
  if (!statsOnlyMode())
    return;

  ResultSummary *summary = &g_result_summary;
  char best_buf[32];
  char worst_buf[32];
  char min_buf[16];
  char max_buf[16];
  char avg_buf[32];
  char first_buf[32];
  char min_count_buf[16];
  char min_best_buf[32];
  char min_atoms_buf[16];
  char max_atoms_buf[16];
  char avg_atoms_buf[32];
  char min_atoms_count_buf[16];
  char min_atoms_best_buf[32];

  if (summary->resultsCount > 0) {
    snprintf(best_buf, sizeof(best_buf), "%.6f", summary->bestRmsd);
    snprintf(worst_buf, sizeof(worst_buf), "%.6f", summary->worstRmsd);
  } else {
    strcpy(best_buf, "NA");
    strcpy(worst_buf, "NA");
  }

  if (summary->pathSamples > 0) {
    snprintf(min_buf, sizeof(min_buf), "%d", summary->minPathLength);
    snprintf(max_buf, sizeof(max_buf), "%d", summary->maxPathLength);
    double avg = summary->sumPathLength / summary->pathSamples;
    snprintf(avg_buf, sizeof(avg_buf), "%.2f", avg);
  } else {
    strcpy(min_buf, "NA");
    strcpy(max_buf, "NA");
    strcpy(avg_buf, "NA");
  }

  if (summary->hasFirstResult) {
    snprintf(first_buf, sizeof(first_buf), "%.3f", summary->firstResultMs);
  } else {
    strcpy(first_buf, "NA");
  }

  if (summary->minPathCount > 0) {
    snprintf(min_count_buf, sizeof(min_count_buf), "%d", summary->minPathCount);
    snprintf(min_best_buf, sizeof(min_best_buf), "%.6f", summary->bestRmsdAtMinPath);
  } else {
    strcpy(min_count_buf, "NA");
    strcpy(min_best_buf, "NA");
  }

  if (summary->atomSamples > 0) {
    snprintf(min_atoms_buf, sizeof(min_atoms_buf), "%d", summary->minAtomCount);
    snprintf(max_atoms_buf, sizeof(max_atoms_buf), "%d", summary->maxAtomCount);
    double avg_atoms = summary->sumAtomCount / summary->atomSamples;
    snprintf(avg_atoms_buf, sizeof(avg_atoms_buf), "%.2f", avg_atoms);
  } else {
    strcpy(min_atoms_buf, "NA");
    strcpy(max_atoms_buf, "NA");
    strcpy(avg_atoms_buf, "NA");
  }

  if (summary->minAtomCountHits > 0) {
    snprintf(min_atoms_count_buf, sizeof(min_atoms_count_buf), "%d", summary->minAtomCountHits);
    snprintf(min_atoms_best_buf, sizeof(min_atoms_best_buf), "%.6f", summary->bestRmsdAtMinAtom);
  } else {
    strcpy(min_atoms_count_buf, "NA");
    strcpy(min_atoms_best_buf, "NA");
  }

  printf("RESULT_SUMMARY results=%d best_rmsd=%s worst_rmsd=%s min_path=%s min_path_count=%s min_path_best_rmsd=%s "
         "max_path=%s avg_path=%s min_atoms=%s min_atoms_count=%s min_atoms_best_rmsd=%s max_atoms=%s avg_atoms=%s "
         "first_result_ms=%s total_time_ms=%.3f cpu_time_ms=%.3f\n",
         summary->resultsCount, best_buf, worst_buf, min_buf, min_count_buf, min_best_buf, max_buf, avg_buf,
         min_atoms_buf, min_atoms_count_buf, min_atoms_best_buf, max_atoms_buf, avg_atoms_buf, first_buf,
         total_wall_ms, total_cpu_ms);
}

/**
 * @brief Finalize a stats-only execution by emitting timings and summaries.
 *
 * @param print_timings When non-zero, also print human-readable timing lines.
 */
static void finalizeStatsOnlyRun(int print_timings) {
  if (!statsOnlyMode())
    return;

  struct timespec end_ts;
  clock_gettime(CLOCK_MONOTONIC, &end_ts);
  clock_t end_clock = clock();

  double wall_ms = timespec_diff_ms(start_time_ts, end_ts);
  double cpu_ms = ((double)(end_clock - start_clock)) * 1000.0 / CLOCKS_PER_SEC;

  if (print_timings) {
    printf("\nExecution time (monotonic): %.3f ms\n", wall_ms);
    printf("Execution time(clock): %.3f ms\n", cpu_ms);
  }
  printResultSummary(wall_ms, cpu_ms);
#ifdef ENABLE_STATS
  emitStatsSummary(stdout);
#endif
}

/**
 * @file output.c
 * @brief Functions for generating and writing output files.
 *
 * This file contains functions for creating directories, writing output files in various formats,
 * and printing molecule and graph information to the console.
 * The functions here are used for generating cages and storing the results.
 *
 */

/**
 * @brief Ensure that a directory exists for the provided absolute or relative path.
 *
 * @param path Directory path that must exist after the call.
 */
static void ensureDirectoryExists(const char *path) {
  if (mkdir(path, 0755) != 0 && errno != EEXIST) {
    fprintf(stderr, "Failed to create directory %s: %s\n", path, strerror(errno));
    exit(EXIT_FAILURE);
  }
}

/**
 * @brief Create (or reuse) the top-level results directory for the current run.
 *
 * Depending on configuration it either reuses CAGE_RESULTS_DIR, builds a dated
 * folder inside resultsTest (stats builds) or creates results/<input>.
 *
 * @param input Base name of the processed input file.
 * @return Newly allocated path to the directory; caller must free it.
 */
char *createDir(char *input) {
  const char *override_dir = getenv("CAGE_RESULTS_DIR");
  if (override_dir != NULL && override_dir[0] != '\0') {
    char *dir_name = malloc((strlen(override_dir) + 1) * sizeof(char));
    if (dir_name == NULL) {
      fprintf(stderr, "Unable to allocate memory for directory name\n");
      exit(EXIT_FAILURE);
    }
    strcpy(dir_name, override_dir);
    ensureDirectoryExists(dir_name);
    return dir_name;
  }

  char *dir_name = malloc(256 * sizeof(char));

#ifdef ENABLE_STATS
  char timestamp[32];
  char dist_type[32];
  if (get_current_distance_type() == DISTANCE_A_STAR) {
    strcpy(dist_type, "AStar");
  } else if (get_current_distance_type() == DISTANCE_SSMTA_STAR) {
    strcpy(dist_type, "SSMTAStar");
  } else if (get_current_distance_type() == DISTANCE_HYBRID) {
    strcpy(dist_type, "HYBRID");
  } else {
    strcpy(dist_type, "Eucli");
  }
  // Get current date-time
  time_t now = time(NULL);
  struct tm *t = localtime(&now);
  strftime(timestamp, sizeof(timestamp), "%Y%m%d_%Hh%Mm%Ss", t);
  ensureDirectoryExists("./resultsTest");
  // Create the input_DATE_TIME directory inside resultsTest
  sprintf(dir_name, "./resultsTest/%s_%s_%s", input, dist_type, timestamp);
  ensureDirectoryExists(dir_name);
#else
  ensureDirectoryExists("./results");
  sprintf(dir_name, "./results/%s", input);
#endif
  ensureDirectoryExists(dir_name);

  return dir_name;
}

/**
 * @brief Create a subdirectory that groups outputs by motif/atom count.
 *
 * @param input Parent directory created by createDir.
 * @param nbmotif The motif number.
 * @return A dynamically allocated string representing the created subdirectory name.
 */
char *createUnderDir(char *input, int nbmotif) {
  char *dir_name = malloc(256 * sizeof(char));
  sprintf(dir_name, "%s/%d", input, nbmotif);
  ensureDirectoryExists(dir_name);

  return dir_name;
}

// /**
//  * @brief Copy the input file to the specified directory with the given name.
//  *
//  * @param inputFile The input file to be copied.
//  * @param dir_name The destination directory.
//  * @param name The name of the file to be copied.
//  */
// void copytoDir(char* inputFile, char* dir_name, char* name) {
//   char cmd[512];

//   sprintf(cmd, "cp %s %s/%s.xyz", inputFile, dir_name, name);
//   if (system(cmd) == -1) {
//     printf("%s",cmd);
//     exit(0);
//   }
// }

// /**
//  * @brief Write the contents of a list to the console.
//  *
//  * @param l The list to be written.
//  */
// void LST_write(List_t* l) {

//   int i;

//   printf("Size : %d, NbElement : %d\n", size(l), LST_nbElements(l));

//   for (i=0; i<size(l); i++)
//     printf(" %d,", elts(l,i));
//   printf("\n");
// }

/**
 * @brief Write the information of an atom in a cage to the console.
 *
 * @param a The atom to be written.
 */
void cageWriteAtom(AtomCage_t *a) {
  int i;

  printf("(%d) (%8.4f, %8.4f, %8.4f) :", flag(a), atomX(a), atomY(a), atomZ(a));

  for (i = 0; i < neighborhoodSize(a); i++)
    // if (neighbor(a,i) != -1)
    printf(" %d,", neighbor(a, i));

  printf("\n");
}

/**
 * @brief Write the information of a cage to the console.
 *
 * @param s The cage to be written.
 */
void cageWrite(Cage_t *s) {
  int i;

  printf("Size : %d, SizeMax : %d\n", cageNbAtom(s), size(s));
  for (i = 0; i < size(s); i++) {
    // if (flag(atom(s,i)) != NOT_DEF_F) {
    printf("%d ", i);
    cageWriteAtom(atom(s, i));
    //}
  }
  printf("\n");
}

/**
 * @brief Writes the contents of a cage structure to a .mol2 file format.
 *
 * This function generates a .mol2 file representing the molecular structure of the cage, including atoms,
 * paths, and bonds. The .mol2 format is commonly used for molecular modeling and visualization.
 *
 * @param output The path to the output .mol2 file.
 * @param s A pointer to the `Cage_t` structure representing the molecular system.
 * @param paths A pointer to the `Paths_t` structure containing path information (e.g., atom paths and their positions).
 * @param interTree An array representing the interconnection tree, where each pair of integers defines the start
 *                  and end points of a connection.
 * @param nb_atoms_paths The number of additional atoms in paths to be included in the file.
 * @param cumul_MSD The cumulative Mean Squared Deviation (MSD) for all paths, used for analysis.
 *
 * @details
 * ### Function Workflow:
 * 1. **File Initialization**:
 *    - The function opens the specified file for writing. If the file cannot be opened, an error message is printed,
 *      and the function exits.
 *
 * 2. **Molecule Section**:
 *    - Writes the general header of the .mol2 file, including the number of atoms and bonds.
 *    - Assumes the "SMALL" and "GASTEIGER" molecule types in the .mol2 format.
 *
 * 3. **Atom Section**:
 *    - Iterates over all vertices in the cage and writes atom information such as ID, type, and coordinates.
 *    - Determines the type of each atom based on its flag (`CYCLE_F`, `HYDRO_PATTERN_F`, etc.).
 *    - For atoms in paths, their positions are derived from the `Paths_t` structure, and additional atoms (e.g.,
 * hydrogens) are handled specifically.
 *
 * 4. **Bond Section**:
 *    - Iterates over the neighborhood of each atom in the cage to define bonds between connected atoms.
 *    - Includes bonds for paths, connecting the atoms defined in the `interTree` array and additional atoms in
 * `Paths_t`.
 *
 * 5. **Error Handling**:
 *    - If any writing operation fails, an error message is printed, and the program exits with an error code.
 *
 * 6. **Memory Management**:
 *    - Allocates memory for indexing cage vertices to atom IDs.
 *    - Ensures all dynamically allocated memory (e.g., `index`) is freed before the function exits.
 *
 * ### Notes:
 * - The function assumes that all vertices and paths in the input structures are properly initialized.
 * - The atom types (e.g., `H`, `C`, `O`, etc.) and flags (`CYCLE_F`, `HYDRO_PATTERN_F`) must conform to the expected
 * format.
 * - The generated file is in the .mol2 format, compatible with most molecular modeling tools.
 *
 * @see CAGE_nbAtom
 * @see CAGE_nbEdges
 * @see Paths_t
 */
void cageWriteMol2(char *output, Cage_t *s, Paths_t *paths, int *interTree, int nb_atoms_paths, double cumul_MSD) {
  FILE *filestream = NULL;
  int ret, i, j, l, ind_last_id_cage;
  int *index = malloc(size(s) * sizeof(int));

  filestream = fopen(output, "w");
  if (filestream == NULL) {
    printf("The file %s could not be open for writting.\n", output);
    if (index)
      free(index);
    return;
  }

  //Ecriture de l'entete
  ret = fprintf(filestream, "# Molecule generated by CagePathGen\n");
  ret = fprintf(filestream, "# Number of paths: %d\n", paths->numPaths);
  ret = fprintf(filestream, "# RMSD: %f\n", cumul_MSD);

  const char *path_length_line = "NA";
  char *allocated_lengths = NULL;
  if (paths->numPaths > 0) {
    size_t capacity = (size_t)paths->numPaths * 12 + 32; // generous per-path buffer
    allocated_lengths = (char *)malloc(capacity);
    if (allocated_lengths != NULL) {
      allocated_lengths[0] = '\0';
      size_t offset = 0;
      for (int path_index = 0; path_index < paths->numPaths; path_index++) {
        if (path_index > 0) {
          int added = snprintf(allocated_lengths + offset, capacity - offset, ", ");
          if (added < 0 || (size_t)added >= capacity - offset) {
            strcpy(allocated_lengths, "TRUNCATED");
            break;
          }
          offset += (size_t)added;
        }

        int pattern_count = paths->curPthPos[path_index] - 2;
        if (pattern_count < 0) {
          pattern_count = 0;
        }
        int written = snprintf(allocated_lengths + offset, capacity - offset, "%d", pattern_count);
        if (written < 0 || (size_t)written >= capacity - offset) {
          strcpy(allocated_lengths, "TRUNCATED");
          break;
        }
        offset += (size_t)written;
      }
      path_length_line = allocated_lengths;
    }
  }

  ret = fprintf(filestream, "# Path_Lengths: %s\n", path_length_line);
  ret = fprintf(filestream, "\n");
  free(allocated_lengths);

  // Ecriture de la molecule
  ret = fprintf(filestream, "@<TRIPOS>MOLECULE\n*****\n");
  ret = fprintf(filestream, " %d %d 0 0 0\n", cageNbAtom(s) + nb_atoms_paths,
                cageNbEdges(s) + (nb_atoms_paths + paths->numPaths));
  ret = fprintf(filestream, "SMALL\nGASTEIGER\n\n");

  // Ecriture des sommets
  ret = fprintf(filestream, "@<TRIPOS>ATOM\n");
  for (i = 0, j = 1; i < size(s); i++) {
    if (flag(atom(s, i)) != NOT_DEF_F) {
      index[i] = j;
      if (flag(atom(s, i)) == CYCLE_F)
        ret = fprintf(filestream, " %3d S", j);
      else if (flag(atom(s, i)) == HYDRO_PATTERN_F)
        if (lstNbElements(neighborhood(atom(s, i))) == 1) {
          if (flag(atom(s, neighbor(atom(s, i), 0))) == HYDRO_PATTERN_F) {
            ret = fprintf(filestream, " %3d H", j);
          } else {
            ret = fprintf(filestream, " %3d U", j);
          }
        } else
          ret = fprintf(filestream, " %3d U", j);
      else if (flag(atom(s, i)) == LINKABLE_F)
        /*if(LST_nbElements(neighborhood(atom(s,i)))>1){
          ret = fprintf(filestream, " %3d C", j);
        } else {*/
        ret = fprintf(filestream, " %3d P", j);
      //}
      else if (flag(atom(s, i)) == OXYGEN_F)
        ret = fprintf(filestream, " %3d O", j);
      else if (flag(atom(s, i)) == NITROGEN_F)
        ret = fprintf(filestream, " %3d N", j);
      else if (flag(atom(s, i)) == CARBON_F)
        ret = fprintf(filestream, " %3d C", j);
      else if (flag(atom(s, i)) == HYDROGEN_F)
        ret = fprintf(filestream, " %3d H", j);
      else
        ret = fprintf(filestream, " %3d Al", j);
      ret = fprintf(filestream, "   %3.4f", atomX(atom(s, i)));
      ret = fprintf(filestream, "   %3.4f", atomY(atom(s, i)));
      ret = fprintf(filestream, "   %3.4f", atomZ(atom(s, i)));
      if (flag(atom(s, i)) == CYCLE_F)
        ret = fprintf(filestream, "   S\n");
      else if (flag(atom(s, i)) == HYDRO_PATTERN_F)
        if (lstNbElements(neighborhood(atom(s, i))) == 1) {
          if (flag(atom(s, neighbor(atom(s, i), 0))) == HYDRO_PATTERN_F) {
            ret = fprintf(filestream, "   H\n");
          } else {
            ret = fprintf(filestream, "   U\n");
          }
        } else
          ret = fprintf(filestream, "   U\n");
      else if (flag(atom(s, i)) == LINKABLE_F)
        /*if(LST_nbElements(neighborhood(atom(s,i)))>1){
          ret = fprintf(filestream, "   C\n");
        }else */
        ret = fprintf(filestream, "   P\n");
      else if (flag(atom(s, i)) == OXYGEN_F)
        ret = fprintf(filestream, "   O\n");
      else if (flag(atom(s, i)) == NITROGEN_F)
        ret = fprintf(filestream, "   N\n");
      else if (flag(atom(s, i)) == CARBON_F)
        ret = fprintf(filestream, "   C\n");
      else if (flag(atom(s, i)) == HYDROGEN_F)
        ret = fprintf(filestream, "   H\n");
      else
        ret = fprintf(filestream, "   Al\n");
      j++;
    } else
      index[i] = -1;
  }
  ind_last_id_cage = j;
  Point_t tmp_p;
  for (int k = 0; k < paths->numPaths; k++) {
    // tmp_p = paths->patterns[indexPointPaths(k,1,0,1,paths->sizeMax)]; //first hydrogen start
    // ret = fprintf(filestream, " %3d H   %3.4f   %3.4f   %3.4f   H\n",j,tmp_p.x,tmp_p.y,tmp_p.z);
    // j++;
    // tmp_p = paths->patterns[indexPointPaths(k,1,0,2,paths->sizeMax)]; //second hydrogen start
    // ret = fprintf(filestream, " %3d H   %3.4f   %3.4f   %3.4f   H\n",j,tmp_p.x,tmp_p.y,tmp_p.z);
    // j++;

    for (int s = 2; s <= paths->curPthPos[k]; s++) {
      for (int a = 0; a < MAX_NB_ATOMS_PATTERN; a++) { // working for simple pattern
        if (a == 0) {
          tmp_p = paths->patterns[indexPointPaths(k, s, paths->positionCurNum[indexPathPosition(k, s, paths->sizeMax)],
                                                  a, paths->sizeMax)];
          ret = fprintf(filestream, " %3d C   %3.4f   %3.4f   %3.4f   C\n", j, tmp_p.x, tmp_p.y, tmp_p.z);
          j++;
        } else {
          tmp_p = paths->patterns[indexPointPaths(k, s, paths->positionCurNum[indexPathPosition(k, s, paths->sizeMax)],
                                                  a, paths->sizeMax)];
          ret = fprintf(filestream, " %3d H   %3.4f   %3.4f   %3.4f   H\n", j, tmp_p.x, tmp_p.y, tmp_p.z);
          j++;
        }
      }
    }
    tmp_p = paths->patterns[indexPointPaths(k, paths->curPthPos[k] + 1, 0, 1,
                                            paths->sizeMax)]; // first hydrogen neighbor end
    ret = fprintf(filestream, " %3d H   %3.4f   %3.4f   %3.4f   H\n", j, tmp_p.x, tmp_p.y, tmp_p.z);
    j++;
    tmp_p = paths->patterns[indexPointPaths(k, paths->curPthPos[k] + 1, 0, 2,
                                            paths->sizeMax)]; // second hydrogen neighbor end
    ret = fprintf(filestream, " %3d H   %3.4f   %3.4f   %3.4f   H\n", j, tmp_p.x, tmp_p.y, tmp_p.z);
    j++;
    tmp_p = paths->patterns[indexPointPaths(k, paths->curPthPos[k] + 2, 0, 1, paths->sizeMax)]; // first hydrogen end
    ret = fprintf(filestream, " %3d H   %3.4f   %3.4f   %3.4f   H\n", j, tmp_p.x, tmp_p.y, tmp_p.z);
    j++;
    tmp_p = paths->patterns[indexPointPaths(k, paths->curPthPos[k] + 2, 0, 2, paths->sizeMax)]; // second hydrogen end
    ret = fprintf(filestream, " %3d H   %3.4f   %3.4f   %3.4f   H\n", j, tmp_p.x, tmp_p.y, tmp_p.z);
    j++;
  }

  // Ecriture des liens
  ret = fprintf(filestream, "\n@<TRIPOS>BOND\n");
  for (i = 0, l = 1; i < size(s); i++) {
    for (j = 0; j < neighborhoodSize(atom(s, i)); j++) {
      if (i < neighbor(atom(s, i), j)) {
        ret = fprintf(filestream, " %3d %3d %3d %3d\n", l, index[i], index[neighbor(atom(s, i), j)], 1);
        l++;
      }
    }
  }
  j = ind_last_id_cage;
  int prec_a, end, cur_a;
  for (int k = 0; k < paths->numPaths; k++) {
    prec_a = interTree[2 * k] + 1;
    // s
    for (int s = 2; s <= paths->curPthPos[k]; s++) {
      for (int a = 0; a < MAX_NB_ATOMS_PATTERN; a++) { // working for simple pattern
        if (a == 0) {
          ret = fprintf(filestream, " %3d %3d %3d %3d\n", l, j, prec_a, 1);
          cur_a = j;
          j++;
          l++;
        } else {
          ret = fprintf(filestream, " %3d %3d %3d %3d\n", l, j, prec_a, 1);
          j++;
          l++;
        }
      }
      prec_a = cur_a;
    }
    end = interTree[2 * k + 1] + 1;
    ret = fprintf(filestream, " %3d %3d %3d %3d\n", l, j, prec_a, 1);
    j++;
    l++;
    ret = fprintf(filestream, " %3d %3d %3d %3d\n", l, j, prec_a, 1);
    j++;
    l++;
    ret = fprintf(filestream, " %3d %3d %3d %3d\n", l, j, end, 1);
    j++;
    l++;
    ret = fprintf(filestream, " %3d %3d %3d %3d\n", l, j, end, 1);
    j++;
    l++;
    ret = fprintf(filestream, " %3d %3d %3d %3d\n", l, prec_a, end, 1);
    l++;
  }
  if (ret < 0) {
    printf("Writting of file %s did not go well.\n", output);
    exit(2);
  }

  free(index);
  fclose(filestream);
}

int cpt_result_cage_egv = 0; // Extern global variable to count the number of cages generated

/**
 * @brief Write the execution time and statistics to a file.
 *
 * This function calculates the execution time of the program and writes it to a specified output file.
 * It also includes the number of interconnection trees generated during the execution.
 *
 * @param outputname The name of the output file where the execution time and statistics will be written.
 *
 * @details
 * ### Function Workflow:
 * 1. **Open File**:
 *    - Opens the specified file in write mode. If the file cannot be opened, an error message is printed.
 *
 * 2. **Calculate Execution Time**:
 *    - Calculates the total execution time in seconds using `difftime`.
 *    - Calculates CPU time used using `clock()`.
 *
 * 3. **Write Execution Time**:
 *    - Writes the formatted execution time to the file, including hours, minutes, and seconds.
 *
 * 4. **Write Interconnection Trees**:
 *    - Writes the number of interconnection trees generated during execution to the file.
 *
 * 5. **Close File**:
 *    - Closes the file after writing all information.
 *
 */
void writeTime(const char *outputname) {
  ensureResultSummaryInitialized();
  if (statsOnlyMode()) {
    finalizeStatsOnlyRun(0);
    return;
  }

  FILE *file = fopen(outputname, "w"); // Open file in write mode
  if (file == NULL) {
    perror("Error opening file");
    return;
  }
  fprintf(file, "Temps\n");
  struct timespec end_ts;
  clock_gettime(CLOCK_MONOTONIC, &end_ts);
  clock_t end_clock = clock();

  double wall_ms = timespec_diff_ms(start_time_ts, end_ts);
  double cpu_ms = ((double)(end_clock - start_clock)) * 1000.0 / CLOCKS_PER_SEC;

  long long total_seconds = (long long)(wall_ms / 1000.0);
  long long hours = total_seconds / 3600;
  total_seconds -= hours * 3600;
  long long minutes = total_seconds / 60;
  long long seconds = total_seconds - minutes * 60;
  fprintf(file, "Execution time : %lld hour(s) %lld minute(s) %lld second(s)\n", hours, minutes, seconds);
  fprintf(file, "Execution time (monotonic): %.3f ms\n", wall_ms);
  fprintf(file, "Execution time(clock): %.3f ms\n", cpu_ms);

  fprintf(file, "\nInterconnection Trees: %d\n", cpt_inter_tree_egv);

  fclose(file); // Close the file
}

/**
 * @brief Writes output files for a given cage structure and its paths.
 *
 * This function generates and writes a `.mol2` file for the given cage structure,
 * including all atoms and paths, and ensures proper directory organization.
 * It updates atom flags as needed for writing and restores them afterward.
 *
 * @param cage A pointer to the `Cage_t` structure representing the molecular system.
 * @param paths A pointer to the `Paths_t` structure containing path information.
 * @param interTree An array representing the interconnection tree, where each pair of integers
 *                  defines the start and end points of a connection.
 * @param options The `Options_t` structure containing input file name, output directory,
 *                and the maximum number of results to be written.
 *
 * @details
 * ### Function Workflow:
 * 1. **Prepare Output Directory**:
 *    - Creates a subdirectory based on the input file name and the number of atoms in the paths.
 *
 * 2. **Count Additional Atoms**:
 *    - Calculates the number of atoms and patterns added by paths (`nb_atoms_paths` and `nb_pattern`).
 *    - Updates atom flags in the cage for writing (sets `CARBON_F` for start and end points of paths).
 *
 * 3. **Write `.mol2` File**:
 *    - Constructs the output file name using the input file name and result index.
 *    - Calls `cageWriteMol2` to write the cage and its paths to a `.mol2` file.
 *
 * 4. **Restore Atom Flags**:
 *    - Resets the flags of start and end points in the cage to `LINKABLE_F` after writing.
 *
 * 5. **Track Results**:
 *    - Tracks the number of results written using a static counter.
 *    - If the maximum number of results is reached, the program exits successfully.
 *
 * ### Notes:
 * - The function ensures that directory creation and file naming are consistent with the input file name.
 * - Atom flags are temporarily modified for writing and restored afterward.
 * - If the maximum number of results (`options.maxResults`) is reached, the program exits with a success message.
 *
 * @see createDir
 * @see createUnderDir
 * @see cageWriteMol2
 *
 */
void writeCageOutput(Cage_t *cage, Paths_t *paths, int *interTree, Options_t options) {
  ensureResultSummaryInitialized();
  const int stats_only = statsOnlyMode();

  char outputname[512];
  char *name = getBasename(options.input);
  char *dir_name = NULL;
  char *under_dir_name = NULL;

  if (!stats_only) {
    dir_name = createDir(name);
#ifdef ENABLE_STATS
    sprintf(outputname, "%s/%s_stats.txt", dir_name, name);
    writeStats(outputname);
#endif
  }

  int nb_atoms_paths = 0;
  double cumul_MSD = 0.0;
  int min_path_length = INT_MAX;
  int max_path_length = INT_MIN;
  int total_path_length = 0;
  int path_length_samples = 0;

  for (int k = 0; k < paths->numPaths; k++) {
    nb_atoms_paths += 4; // for the two hydrogen of start and of end.
    int path_len = paths->curPthPos[k] - 1; // number of patterns in the path, excluding the start and start's neighbors but it starts count from 0 so 0 and 1 is not count as pattern and the value is the last case where we have a pattern (example for a path of size 3 we have patterns in case 2, 3, 4 so curPthPos =4 and number of patterns =3)
    if (path_len < 0)
      path_len = 0;
    if (path_len < min_path_length)
      min_path_length = path_len;
    if (path_len > max_path_length)
      max_path_length = path_len;
    total_path_length += path_len;
    path_length_samples++;

    for (int s = 2; s <= paths->curPthPos[k]; s++) {
      nb_atoms_paths += MAX_NB_ATOMS_PATTERN;
    }
    cumul_MSD += paths->pathMSD[k];
    flag(atom(cage, interTree[2 * k])) = CARBON_F;
    flag(atom(cage, interTree[2 * k + 1])) = CARBON_F;
  }

  if (path_length_samples == 0) {
    min_path_length = 0;
    max_path_length = 0;
  }

  if (!stats_only) {
    under_dir_name = createUnderDir(dir_name, nb_atoms_paths);
  }

  if (paths->numPaths == 0) {
    fprintf(stderr, "Warning: writeCageOutput called without any paths.\n");
    if (dir_name)
      free(dir_name);
    if (under_dir_name)
      free(under_dir_name);
    free(name);
    return;
  }

  cumul_MSD /= paths->numPaths; // divide by the number of paths to have the mean SD
  cumul_MSD = sqrt(cumul_MSD);  // square root to have the RMSD

  if (!stats_only) {
    sprintf(outputname, "%s/%s_mot%d.mol2", under_dir_name, name, cpt_result_cage_egv);
    cageWriteMol2(outputname, cage, paths, interTree, nb_atoms_paths, cumul_MSD);
  } else {
    recordResultMetrics(cumul_MSD, min_path_length, max_path_length, total_path_length, path_length_samples,
                       nb_atoms_paths);
  }

  for (int k = 0; k < paths->numPaths; k++) {
    flag(atom(cage, interTree[2 * k])) = LINKABLE_F;
    flag(atom(cage, interTree[2 * k + 1])) = LINKABLE_F;
  }
  printf("### Result %d generated ###\n", cpt_result_cage_egv);

  cpt_result_cage_egv++;
  if (cpt_result_cage_egv == options.maxResults) {
#ifdef ENABLE_STATS
    if (stats_only) {
      finalizeStatsOnlyRun(1);
    } else {
      printStats();
    }
#else
    if (stats_only) {
      finalizeStatsOnlyRun(1);
    } else {
      sprintf(outputname, "%s/%s_time.txt", dir_name, name);
      writeTime(outputname);
    }
#endif
    if (dir_name)
      free(dir_name);
    if (under_dir_name)
      free(under_dir_name);
    free(name);
    printf("\nCages Results: %d\n", cpt_result_cage_egv);
    printf("Interconnection Trees: %d\n", cpt_inter_tree_egv);
    exit(EXIT_SUCCESS);
  }

  if (dir_name)
    free(dir_name);
  if (under_dir_name)
    free(under_dir_name);
  free(name);
}

/**
 * @brief Persist or print the solver parameters for traceability.
 *
 * @param options Parsed CLI options describing the current run.
 */
void writeParameters(Options_t options) {
  char outputname[512];
  char *name = getBasename(options.input);
  ensureResultSummaryInitialized();
  DistanceType type = get_current_distance_type();
  const char *distance_str = distance_type_to_string(type);
  const char *tree_mode = options.sortInterTreesBeforePaths ? "store-sort" : "on-the-fly";

  if (statsOnlyMode()) {
    printf("=== Parameters (stats-only) ===\n");
    printf("Input file        : %s\n", name);
    printf("Mocap Number (n)  : %s\n", options.numMoc);
    printf("Max Path Size (s) : %d\n", options.sizeMaxPath);
    printf("Max Results (r)   : %d\n", options.maxResults);
    printf("Distance Type     : %s\n", distance_str);
    printf("Path boundary filt: %s\n", options.enablePathBoundary ? "enabled" : "disabled");
    printf("Best path cutoff  : %s\n", options.enableDynamicPathLimit ? "enabled" : "disabled");
    printf("Inter-tree mode   : %s\n", tree_mode);
    free(name);
    return;
  }

  char *dir_name = createDir(name);

  sprintf(outputname, "%s/%s_parametres.txt", dir_name, name);
  FILE *filestream = NULL;
  filestream = fopen(outputname, "w");

  // Get current date/time
  time_t now = time(NULL);
  struct tm *t = localtime(&now);
  char datetime[64];
  strftime(datetime, sizeof(datetime), "%Y-%m-%d %H:%M:%S", t);

  char git_hash[64] = "unknown";
  FILE *git_pipe = popen("git rev-parse --short HEAD 2>/dev/null", "r");
  if (git_pipe) {
    if (fgets(git_hash, sizeof(git_hash), git_pipe)) {
      // Remove newline if present
      size_t len = strlen(git_hash);
      if (git_hash[len - 1] == '\n')
        git_hash[len - 1] = '\0';
    }
    pclose(git_pipe);
  }

  int ret = fprintf(filestream,
                    "=== Parameters Used ===\n"
                    "Input file        : %s\n"
                    "Mocap Number (n)  : %s\n"
                    "Max Path Size (s) : %d\n"
                    "Max Results (r)   : %d\n"
                    "Distance Type     : %s\n"
                    "Path boundary filt: %s\n"
                    "Best path cutoff  : %s\n"
                    "Inter-tree mode   : %s\n"
                    "Date & Time       : %s\n"
                    "Git Commit        : %s\n",
                    name, options.numMoc, options.sizeMaxPath, options.maxResults, distance_str,
                    options.enablePathBoundary ? "enabled" : "disabled",
                    options.enableDynamicPathLimit ? "enabled" : "disabled", tree_mode, datetime, git_hash);
  if (ret < 0) {
    printf("Writting of file %s did not go well.\n", outputname);
    exit(2);
  }
  fclose(filestream);

  free(name);
  free(dir_name);
}

/**
 * @brief Flush runtime statistics when an external signal requests termination.
 *
 * When running in stats-only mode this emits RESULT_SUMMARY and STATS_SUMMARY
 * immediately. Otherwise it falls back to printStats (if enabled) so at least
 * the interval metrics are printed before exiting.
 */
void flushStatsOnSignal(void) {
  ensureResultSummaryInitialized();
  if (statsOnlyMode()) {
    finalizeStatsOnlyRun(1);
#ifdef ENABLE_STATS
  } else {
    printStats();
#endif
  }
}

#ifdef ENABLE_STATS
/**
 * @brief Compute mean interval statistics derived from global counters.
 *
 * @param mean_with_intersection Optional output for branches containing intersections.
 * @param mean_overall Optional output averaged over all branches.
 * @param mean_coverage Optional output for average coverage percentage.
 */
static void computeIntervalStats(double *mean_with_intersection, double *mean_overall, double *mean_coverage) {
  if (mean_with_intersection)
    *mean_with_intersection = (nb_use_inter > 0) ? (cumul_nb_interval / nb_use_inter) : 0.0;
  if (mean_overall)
    *mean_overall = (nb_use_inter + nb_no_inter) > 0 ? (cumul_nb_interval / (nb_use_inter + nb_no_inter)) : 0.0;
  if (mean_coverage)
    *mean_coverage = (nb_use_inter > 0) ? (cumul_covered / nb_use_inter) : 0.0;
}

/**
 * @brief Emit the STATS_SUMMARY line to the requested stream.
 *
 * @param stream Destination stream (stdout or a file) receiving the summary.
 */
static void emitStatsSummary(FILE *stream) {
  double mean_with = 0.0;
  double mean_all = 0.0;
  double mean_cov = 0.0;
  computeIntervalStats(&mean_with, &mean_all, &mean_cov);
  fprintf(stream,
          "STATS_SUMMARY branches=%d without_intersections=%d with_intersections=%d mean_intervals_with=%.2f "
          "mean_intervals_overall=%.2f max_intervals=%d mean_coverage_pct=%.2f max_coverage_pct=%.2f"
      " total_collisions=%d circle_blocked=%d boundary_allowed=%d boundary_blocked=%d\n",
          nb_branches, nb_no_inter, nb_use_inter, mean_with, mean_all, max_nb_interval, mean_cov, max_covered,
      cpt_collision, cpt_no_next_point, nb_path_boundary_allowed, nb_path_boundary_blocked);
}

/**
 * @brief Print detailed interval statistics to stdout.
 */
void printStats() {
  double mean_with = 0.0;
  double mean_all = 0.0;
  double mean_cov = 0.0;
  computeIntervalStats(&mean_with, &mean_all, &mean_cov);

  printf("\nInterval statistics:\n");
  printf("Branches explored: %d\n", nb_branches);
  printf("Iterations without intersections: %d, with intersections: %d\n", nb_no_inter, nb_use_inter);
  printf("Intervals per branch (with intersections): %.2f\n", mean_with);
  printf("Intervals per branch (overall): %.2f (max %d)\n", mean_all, max_nb_interval);
  printf("Coverage: mean %.2f%%, max %.2f%%\n", mean_cov, max_covered);
  printf("Total collisions: %d\n", cpt_collision);
  printf("Circle blocked (no next position): %d\n", cpt_no_next_point);
    printf("DIST_PATH_BOUNDARY filter: allowed %d, blocked %d\n", nb_path_boundary_allowed,
      nb_path_boundary_blocked);
  emitStatsSummary(stdout);
}

/**
 * @brief Write interval and timing statistics to a dedicated file.
 *
 * @param outputname Path of the destination file.
 */
void writeStats(const char *outputname) {
  FILE *file = fopen(outputname, "w"); // Open file in write mode
  if (file == NULL) {
    perror("Error opening file");
    return;
  }

  double mean_with = 0.0;
  double mean_all = 0.0;
  double mean_cov = 0.0;
  computeIntervalStats(&mean_with, &mean_all, &mean_cov);

  // Print statistics to the file
  fprintf(file, "\nInterval statistics:\n");
  fprintf(file, "Branches explored %d\n", nb_branches);
  fprintf(file, "Iterations without intersections %d, with intersections %d\n", nb_no_inter, nb_use_inter);
  fprintf(file, "Intervals per branch (with intersections) %.2f\n", mean_with);
  fprintf(file, "Intervals per branch (overall) %.2f (max %d)\n", mean_all, max_nb_interval);
  fprintf(file, "Coverage: mean %.2f%%%%, max %.2f%%%%\n", mean_cov, max_covered);
  fprintf(file, "Total collisions %d\n", cpt_collision);
  fprintf(file, "Circle blocked (no next position) %d\n", cpt_no_next_point);
    fprintf(file, "DIST_PATH_BOUNDARY filter allowed %d blocked %d\n", nb_path_boundary_allowed,
      nb_path_boundary_blocked);
  emitStatsSummary(file);
  fprintf(file, "\nTemps\n");
  struct timespec end_ts;
  clock_gettime(CLOCK_MONOTONIC, &end_ts);
  clock_t end_clock = clock();

  double wall_ms = timespec_diff_ms(start_time_ts, end_ts);
  double cpu_ms = ((double)(end_clock - start_clock)) * 1000.0 / CLOCKS_PER_SEC;

  long long total_seconds = (long long)(wall_ms / 1000.0);
  long long hours = total_seconds / 3600;
  total_seconds -= hours * 3600;
  long long minutes = total_seconds / 60;
  long long seconds = total_seconds - minutes * 60;
  fprintf(file, "Execution time : %lld hour(s) %lld minute(s) %lld second(s)\n", hours, minutes, seconds);
  fprintf(file, "Execution time (monotonic): %.3f ms\n", wall_ms);
  fprintf(file, "Execution time(clock): %.3f ms\n", cpu_ms);

  fclose(file);                                      // Close the file
  printf("Statistics written to: %s\n", outputname); // Notify user
}
#endif