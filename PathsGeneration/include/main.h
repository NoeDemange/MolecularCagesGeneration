#ifndef _MAIN_H
#define _MAIN_H

/**
 * @file main.h
 * @brief Main Header File
 *
 * This file contains the definition of the Options_t structure and function
 * related to program usage.
 */

/**
 * @brief Structure to store command-line options.
 *
 * This structure holds command-line options that can be passed to the program.
 * It includes fields for input filename, alpha value, maximum size, and maximum results.
 */
typedef struct {
  char *input;                      /**< Input filename */
  char *numMoc;                     /**< Moc Number */
  int sizeMaxPath;                  /**< Maximum size of path*/
  int maxResults;                   /**< Maximum results */
  int isBannedEdges;                /**< Banned edges 0 = false, 1 = true*/
  int oneCageByInterconnectionTree; /**< One cage by interconnection tree 0 = false, 1 = true */
  int enablePathBoundary;           /**< Enable DIST_PATH_BOUNDARY pruning 0 = false, 1 = true */
  int enableDynamicPathLimit;       /**< Enable adaptive cut-off once a shorter path is discovered */
  int sortInterTreesBeforePaths;    /**< 1 = store & sort interconnection trees, 0 = generate on-the-fly */
} Options_t;

/**
 * @brief Print the usage of the program.
 *
 * This function prints information about how to use the program, including the available
 * command-line options and their descriptions.
 */
void usage();

void timeoutHandler(int signum);

#endif