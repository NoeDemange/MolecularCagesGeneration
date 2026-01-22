#ifndef __ASSEMBLY_H
#define __ASSEMBLY_H

#include "main.h"
#include "structure.h"
#include "substrat.h"

/** @file assembly.h
 *  @brief Functions related to generating whole cages.
 */

#ifdef ENABLE_STATS
// variable global statistique
extern double cumul_nb_interval;
extern int max_nb_interval;
extern int nb_use_inter;
extern int nb_no_inter;
extern double max_covered;
extern double cumul_covered;
extern int cpt_collision;
extern int cpt_no_next_point;
extern int nb_branches;
extern int nb_path_boundary_allowed;
extern int nb_path_boundary_blocked;
#endif

int effective_growth_limit(Paths_t *paths, int path_index, Options_t options);

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
                   Options_t options, int *list_banned_edges, int *size_list_banned_edges);

#endif