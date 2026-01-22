/**
 * @file main.c
 * @brief Main program for generating cages for a specific substrate.
 *
 * This program reads a substrate molecule from an input .moc2 file, reads a
 * partial cage from an input .moc2 and generates whole cages.
 */

#include <getopt.h>
#include <math.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "constant.h"
#include "distance.h"
#include "interconnection.h"
#include "main.h"
#include "output.h"
#include "structure.h"
#include "substrat.h"
#include "util.h"

// for stats
#include "assembly.h"

struct timespec start_time_ts;
clock_t start_clock;

/**
 * Main function to run the cage generation process.
 *
 * @param argc Number of command-line arguments.
 * @param argv Array of command-line arguments.
 * @return Exit status of the program.
 */
int main(int argc, char **argv) {

  // // Set the default distance function based on environment variable
  set_distance_function_from_env();
  configure_path_boundary_filter_from_env();
  // // Print the distance type being used
  // if (get_current_distance_type() == DISTANCE_A_STAR) {
  //   printf("Running A* specific code...\n");
  //   // Execute A*-specific logic
  // } else {
  //   printf("Running Euclidean-specific code...\n");
  //   // Execute Euclidean-specific logic
  // }

  // Set up the signal handler for graceful termination
  signal(SIGTERM, timeoutHandler);
  signal(SIGINT, timeoutHandler);

  clock_gettime(CLOCK_MONOTONIC, &start_time_ts);
  start_clock = clock();

  /******************************** Options *****/
  int opt;
  Options_t options = {NULL,
                       NULL,
                       DEFLT_SIZEMAX,
                       DEFLT_MAX_RESULTS,
                       DEFLT_BANNED_EDGES,
                       DEFLT_ONE_CAGE_BY_INTERCONNECTION_TREE,
                       DEFLT_PATH_BOUNDARY,
                       DEFLT_DYNAMIC_PATH_LIMIT,
                       DEFLT_SORT_INTERCONNECTION_TREES};
  options.enablePathBoundary = is_path_boundary_filter_enabled();
  int boundary_cli_override = 0;

  while ((opt = getopt(argc, argv, OPTSTR)) != EOF) {
    switch (opt) {
    case 'i':
      options.input = optarg;
      break;
    case 'n':
      options.numMoc = optarg;
      break;
    case 's':
      options.sizeMaxPath = atoi(optarg);
      break;
    case 'r':
      options.maxResults = atoi(optarg);
      break;
    case 'b':
      options.isBannedEdges = atoi(optarg);
      break;
    case 't':
      options.oneCageByInterconnectionTree = atoi(optarg);
      break;
    case 'p':
      options.enablePathBoundary = atoi(optarg);
      set_path_boundary_filter_enabled(options.enablePathBoundary != 0);
      boundary_cli_override = 1;
      break;
    case 'l':
      options.enableDynamicPathLimit = atoi(optarg);
      break;
    case 'g':
      options.sortInterTreesBeforePaths = atoi(optarg);
      break;

    case 'h':
    default:
      usage();
      break;
    }
  }

  if (!boundary_cli_override) {
    options.enablePathBoundary = is_path_boundary_filter_enabled();
  }

  if (options.input == NULL) {
    fprintf(stderr, "Input name of the substrate is missing. Give a name like "
                    "the folders in demos\n");
    usage();
    exit(EXIT_FAILURE);
  }

  if (options.numMoc == NULL) {
    int result = fprintf(stderr, "Moc number is missing. Give a moc number in the substrate folder\n");
    if (result < 0) {
      perror("Error writing\n");
    }
    usage();
    exit(EXIT_FAILURE);
  }

  printf("num position %d\n", NUMBER_POSITION_AX1E3);

  /******************************* Execution *****/
  writeParameters(options);
  double **substrat_t = NULL;
  GridSubstrat grid_sub = importSubstratToGrid(options.input, &substrat_t);

  Cage_t *import_cage = cageImport(options.input, options.numMoc);
  findInterconnection(import_cage, &grid_sub, &substrat_t, options);
  freeGridSubstrat(grid_sub);
  free2DDouble(substrat_t, grid_sub.substratSize);
  cageDelete(import_cage);

#ifdef ENABLE_STATS
  printStats();
#endif

  printf("\nCages Results: %d\n", cpt_result_cage_egv);
  printf("Interconnection Trees: %d\n", cpt_inter_tree_egv);
  /************************************ Time *****/

  struct timespec end_ts;
  clock_gettime(CLOCK_MONOTONIC, &end_ts);
  clock_t end_clock = clock();

  double elapsed_ms = timespec_diff_ms(start_time_ts, end_ts);
  double cpu_time_ms = ((double)(end_clock - start_clock)) * 1000.0 / CLOCKS_PER_SEC;
  printf("\nExecution time (monotonic): %.3f ms\n", elapsed_ms);
  printf("Execution time(clock): %.3f ms\n", cpu_time_ms);

  const char *disable_write = getenv("CAGE_DISABLE_WRITE");
  int stats_only = (disable_write != NULL && disable_write[0] != '\0' && disable_write[0] != '0');
  if (stats_only) {
    writeTime(NULL);
  } else {
    char outputname[512];
    char *name = getBasename(options.input);
    char *dir_name = createDir(name);
    sprintf(outputname, "%s/%s_time.txt", dir_name, name);
    writeTime(outputname);
    free(name);
    free(dir_name);
  }

  /************************************ Exit *****/
  return EXIT_SUCCESS;
}

/**
 * Print the usage of the program.
 */
void usage() {
        fprintf(stderr, USAGE_FMT, DEFLT_SIZEMAX, DEFLT_MAX_RESULTS, DEFLT_BANNED_EDGES,
          DEFLT_ONE_CAGE_BY_INTERCONNECTION_TREE, DEFLT_PATH_BOUNDARY, DEFLT_DYNAMIC_PATH_LIMIT,
          DEFLT_SORT_INTERCONNECTION_TREES);
  exit(EXIT_FAILURE);
}

/**
 * @brief Handles termination signals (SIGTERM/SIGINT).
 *
 * This function is invoked when a timeout signal is received. It prints
 * the values of global variables `cptResultCage_EGV` and `cptInterTree_EGV`
 * (representing the number of results and interconnection trees, respectively)
 * to `stderr`, and then gracefully terminates the program.
 *
 * @param signum The signal number that triggered the handler (e.g., SIGTERM).
 *
 * @note This function uses the `exit(1)` call to terminate the program.
 *       Ensure that any required cleanup has been completed before invoking
 * this handler.
 */
void timeoutHandler(int signum) {
  fprintf(stderr, "Termination signal (%d) received!\n", signum);
  flushStatsOnSignal();
  fprintf(stderr, "Global Variable Values:\n");
  fprintf(stderr, "  Cages Results: %d\n", cpt_result_cage_egv);
  fprintf(stderr, "  Interconnection Trees: %d\n", cpt_inter_tree_egv);
  exit(EXIT_FAILURE); // Exit the program
}