#ifndef _INTERCONNECTION_H
#define _INTERCONNECTION_H

#include "main.h"
#include "structure.h"
#include "substrat.h"

extern int cpt_inter_tree_egv; // Extern global variable to count the number of interconnection trees

/**
 * @brief Finds interconnection trees among the linkable vertices in a cage structure.
 *
 * This function identifies linkable vertices in a cage structure, determines their connected components,
 * and then calls the EnumInterconnectionTrees function to generate all interconnection trees.
 *
 * @param cage A pointer to the cage structure representing the molecular system.
 * @param gridSub A pointer to the `GridSubstrat` structure representing the grid substrate of the molecular system.
 * @param substrat_t A pointer on table of positions of substrat's atoms.
 * @param options An `Options_t` structure containing configuration parameters for the computation.
 *
 * @details
 * The function first computes the connected components of the vertices in the cage and then
 * filters out the linkable vertices (those flagged as `LINKABLE_F`). It sets up the component arrays
 * and initializes a recursive search for interconnection trees, which are trees that span all components
 * by linking vertices between them. It frees allocated memory after the computation.
 */
void findInterconnection(Cage_t *cage, GridSubstrat *gridSub, double ***substrat_t, Options_t options);

#endif