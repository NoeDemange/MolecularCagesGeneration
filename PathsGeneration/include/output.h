#ifndef __OUTPUT_H
#define __OUTPUT_H

#include "main.h"
#include "structure.h"

/**
 * @file output.h
 * @brief Output Header File
 *
 * This file contains function for handling various output-related operations.
 */

extern int cpt_result_cage_egv; // Extern global variable to count the number of cages generated

/**
 * @brief Create a directory with the given input name in the "../results" folder.
 *
 * @param input The input filename.
 * @return A dynamically allocated string representing the created directory name.
 */
char *createDir(char *);

/**
 * @brief Write the contents of a list to the console.
 *
 * @param l The list to be written.
 */
void lstWrite(List_t *);

/**
 * @brief Write the information of a cage to the console.
 *
 * @param s The cage to be written.
 */
void cageWrite(Cage_t *);

void cageWriteMol2(char *output, Cage_t *s, Paths_t *paths, int *interTree, int nbAtomPaths, double cumul_MSD);

/**
 * @brief Write the output files for a given cage.
 *
 * @param options The options with input file name and number of results.
 * @param s The cage to be written.
 * @param paths The paths of the cage.
 */
void writeCageOutput(Cage_t *cage, Paths_t *paths, int *interTree, Options_t options);

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
void writeTime(const char *outputname);

/**
 * @brief Write the parameters of execution.
 *
 * @param options Options_t type.
 */
void writeParameters(Options_t options);

/**
 * @brief Emit partial statistics when the solver receives a termination signal.
 */
void flushStatsOnSignal(void);

#ifdef ENABLE_STATS
/**
 * @brief print the stats of execution.
 */
void printStats();

/**
 * @brief write the stats of execution.
 *
 * @param outpuname the path of the file to print.
 */
void writeStats(const char *outputname);
#endif

#endif