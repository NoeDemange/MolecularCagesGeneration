#ifndef _CONSTANT_H
#define _CONSTANT_H

/** @file constant.h
 *  @brief Constants and options for the program.
 */

// Main : Options for execution

/** @def OPTSTR
 *  @brief Command-line options string for getopt.
 */
#define OPTSTR "i:n:s:r:b:t:p:l:g:h"

/** @def USAGE_FMT
 *  @brief Usage format string for displaying program options.
 */
#define USAGE_FMT                                                                                                      \
  "usage : [-i input directory] [-n moc number] [-s sizemaxpath (default : %d)] [-r maxresults (default : %d)] [-b "   \
  "isBannedEdges (default : %d)] [-t oneCageByInterconnectionTree (default : %d)] [-p enablePathBoundary (default : %d)] [-l " \
  "bestPathCutoff (default : %d)] [-g interTreeMode (default : %d, 0=on-the-fly, 1=store-sort)] " \
  "[-h]\n"

/** @def DEFLT_SIZEMAX
 *  @brief Default maximum size of a path created in number of patterns.
 */
#define DEFLT_SIZEMAX 5

/** @def DEFLT_MAX_RESULTS
 *  @brief Default number of results generated.
 */
#define DEFLT_MAX_RESULTS 10

/** @def DEFLT_BANNED_EDGES
 *  @brief Default value for banned edges option.
 *  0 = false, 1 = true
 */
#define DEFLT_BANNED_EDGES 1 // Default value for banned edges option

/** @def DEFLT_ONE_CAGE_BY_INTERCONNECTION_TREE
 *  @brief Default value for one cage by interconnection tree option.
 *  0 = false, 1 = true
 */
#define DEFLT_ONE_CAGE_BY_INTERCONNECTION_TREE 1 // Default value for one cage by interconnection tree option

/** @def DEFLT_PATH_BOUNDARY
 *  @brief Default toggle for the DIST_PATH_BOUNDARY pruning (1 = enabled).
 */
#define DEFLT_PATH_BOUNDARY 1

/** @def DEFLT_DYNAMIC_PATH_LIMIT
 *  @brief Default toggle for the adaptive best-path cut-off (0 = disabled).
 */
#define DEFLT_DYNAMIC_PATH_LIMIT 1

/** @def DEFLT_SORT_INTERCONNECTION_TREES
 *  @brief Default strategy for interconnection trees (0 = on-the-fly, 1 = store & sort).
 *  The default value keeps the legacy on-the-fly behavior; switch to 1 to buffer and sort.
 */
#define DEFLT_SORT_INTERCONNECTION_TREES 0

// Structure

/** @def REALLOCSIZE
 *  @brief Size of reallocation for dynamic memory allocation.
 */
#define REALLOCSIZE 4 // TODO could it be decreased?

// Distance

/** @def DIST_SIMPLE
 *  @brief Simple covalent bond size.
 */
#define DIST_SIMPLE 1.5

/** @def DIST_ERROR
 *  @brief Acceptable error for the covalent bond size.
 */
#define DIST_ERROR 0.5

/** @def MINDIS
 *  @brief Minimal distance between two atoms (otherwise they are merged).
 */
#define MINDIS 0.75

/** @def DIST_ATOM_H
 *  @brief Distance between one atom and an atom of hydrogen.
 */
#define DIST_ATOM_H (DIST_SIMPLE / 2) + (MINDIS / 2)

/** @def DIST_GAP_CAGE
 *  @brief Distance between two cage's atoms (defined by trial). distance H-C.
 */
#define DIST_GAP_CAGE ((DIST_SIMPLE / 2) + (MINDIS / 2) - 0.0001)

/** @def DIST_GAP_SUBSTRATE
 *  @brief Distance between an atom of the cage and an atom of the substrate (at least a hydrogen bond size by trial).
 */
#define DIST_GAP_SUBSTRATE 1.7 // We said 1.8 but with cycle linking pattern is doesn't work.

/** @def DIST_GAP
 *  @brief Maximum distance between two atoms (defined by trial).
 */
#define DIST_GAP_MAX ((DIST_GAP_CAGE > DIST_GAP_SUBSTRATE) ? DIST_GAP_CAGE : DIST_GAP_SUBSTRATE)

/** @def DIST_SIMPLE_PATTERN
 *  @brief AC/2, with ABC a triangle where each of its vertex is an atom of the path involved in a simple pattern,
 *         AB = BC = 1.5 and angle ABC = 109°C.
 *            B
 *         \ / \ /
 *          A   C
 */
#define DIST_SIMPLE_PATTERN 1.22

/** @def DIST_PATH_BOUNDARY
 * @brief Distance maximum (euclidean) between two atoms to be connected by a path with nb_paterns_remaining.
 *
 */
#define DIST_PATH_BOUNDARY(nb_patterns_remaining) DIST_SIMPLE_PATTERN *nb_patterns_remaining + DIST_SIMPLE_PATTERN

// GridSubstrat

#define STEP_GRID DIST_GAP_SUBSTRATE // minimum DIST_GAP_SUBSTRATE

// Angle

/** @def END_ANGLE
 *  @brief The end angle for path choice.
 */
#define END_ANGLE 109.47

/** @def ANGLE_ERROR
 *  @brief Acceptable error for angle calculations.
 */
#define ANGLE_ERROR 10

// Path choice

// Constants for next positions
/** @def DIST_START_CENTERCIRCLE_C
 * @brief Distance from the start atom to the center of the circle used for generating positions.
 * This distance is used to ensure that the next atom is positioned correctly in relation to the start atom.
 * This value is derived from the geometry of the molecular structure being modeled.
 */
#define DIST_START_CENTERCIRCLE_C 0.499969871352356 // distance from start to center of circle

/** @def CIRCLE_RADIUS_C
 * @brief Defines the radius of the circle used in calculations.
 *
 * This macro represents the fixed radius of the circle, used for generating
 * positions and distance calculations.
 */
#define CIRCLE_RADIUS_C 1.414224214097577

/** @def CIRCLE_RADIUS_SQUARED_C
 * @brief Defines the squared radius of the circle.
 *
 * This macro precomputes the squared value of the circle's radius
 * to optimize distance calculations and avoid redundant multiplications.
 */
#define CIRCLE_RADIUS_SQUARED_C (CIRCLE_RADIUS_C * CIRCLE_RADIUS_C)

/** @def DIST_START_CENTERCIRCLE_H
 * @brief Distance from the start atom to the center of the circle used for generating hydrogen positions.
 *
 * This value is used to ensure that the hydrogen atoms are positioned correctly in relation to the start atom.
 * It is derived from the geometry of the molecular structure being modeled.
 */
#define DIST_START_CENTERCIRCLE_H 0.366644572325061 // distance from start to center of circle Hydrogen

/** @def CIRCLE_RADIUS_H
 * @brief Defines the radius of the circle used for hydrogen positions.
 *
 * This macro represents the fixed radius of the circle used for generating
 * hydrogen positions in relation to the start atom.
 */
#define CIRCLE_RADIUS_H 1.037097757004890

/** @def CIRCLE_RADIUS_SQUARED_H
 * @brief Defines the squared radius of the circle used for hydrogen positions.
 *
 * This macro precomputes the squared value of the circle's radius for hydrogen positions
 * to optimize distance calculations and avoid redundant multiplications.
 */
#define CIRCLE_RADIUS_SQUARED_H (CIRCLE_RADIUS_H * CIRCLE_RADIUS_H)

/** @def ANGLE_SHIFT
 * @brief Defines the angle shift applied to theta when converting from carbon to hydrogen positions.
 *
 * This shift is necessary to account for the geometric arrangement of hydrogen atoms
 * in relation to the carbon atom in the molecular structure.
 * It is equivalent to 119.997010192191553 degrees.
 * This value is used to adjust the angle when calculating hydrogen positions
 * based on the carbon atom's position.
 * It ensures that the hydrogen atoms are positioned correctly in relation to the carbon atom.
 */
#define ANGLE_SHIFT 2.094342920402936

// Floating-point precision
/** @def EPSILON_ANGLE
 *  @brief Small tolerance for floating-point angle comparisons to handle precision errors.
 */
#define EPSILON_ANGLE 1e-12

/** @def POSITIONNING_CORRECTION
 *  @brief Defines the shift apply to theta when the optimal is in a interval and we need to pass from 0 to 2*Pi or from
 * 2Pi to 0.
 */
#define POSITIONNING_CORRECTION 0.001;

//

/** @def NUMBER_POSITION_AX1E3
 *  @brief Number of positions to keep for AX1E3 patterns.
 */
#ifndef NUMBER_POSITION_AX1E3
#define NUMBER_POSITION_AX1E3 3
#endif

/** @def NUMBER_POSITION_PATHS
 *  @brief Number of positions to keep to make a path.
 */
#ifndef NUMBER_POSITION_PATHS
#define NUMBER_POSITION_PATHS NUMBER_POSITION_AX1E3 // MAke a max between the differents patterns.
#endif

/** @def THRESHOLD_ANGLE
 *  @brief Threshold angle in radians for position generation. The minimum space between two positions should be at least this value.
 */
#ifndef THRESHOLD_ANGLE
#define THRESHOLD_ANGLE (15 *(M_PI / 180)) // in radian
#endif

/** @def MAX_POSITION_KEEP_360
 *  @brief Maximum number of positions to keep for 360° rotation.
 */
#ifndef MAX_POSITION_KEEP_360
#define MAX_POSITION_KEEP_360 12
#endif

/** @def DISCRITIZATION_TYPE
 *  @brief Type of discretization to use for generating positions when distance is not euclidean.
 */
#ifndef DISCRITIZATION_TYPE
#define DISCRITIZATION_TYPE 1 // 0 : global, 1 : local
#endif

/********* not to be modified (the incremental order must be preserved) */

// Flags atoms in the envelope and cage

/** @def NOT_DEF_F
 *  @brief Atom not used.
 */
#define NOT_DEF_F -1

/** @def SHELL_F
 *  @brief Atom of the shell.
 */
#define SHELL_F 0

/** @def LINKABLE_F
 *  @brief Atom at the edge of a pattern that can still make another bond (unless it's a hydrogen).
 */
#define LINKABLE_F 1

/** @def CYCLE_F
 *  @brief Atom in an aromatic ring that can't make another bond.
 */
#define CYCLE_F 2

/** @def HYDRO_PATTERN_F
 *  @brief Atom involved in a hydrogen pattern that can't make another bond.
 */
#define HYDRO_PATTERN_F 3

/***********************************************************************/

// Path (in the cage)

/** @def NB_PATTERNS
 *  @brief Number of patterns available to make a path.
 */
#define NB_PATTERNS 1 // 2 //with cycle

/**
 * @def MAX_NB_ATOMS_PATTERN
 * @brief Maximum number of atoms in a pattern (number of atoms in cycle : 7 C + 6 H).
 */
#define MAX_NB_ATOMS_PATTERN 3 // 13

/**
 * @def SIMPLE_PATTERN
 * @brief Number for the simple pattern, a carbon.
 */
#define SIMPLE_PATTERN 0

/** @def CYCLE_PATTERN
 *  @brief Number for the cycle pattern.
 */
#define CYCLE_PATTERN 1

// Flags atoms in paths, ! must be different of flags in the envelope

/** @def CARBON_F
 *  @brief Flag for carbon atom in paths.
 */
#define CARBON_F 6

/** @def NITROGEN_F
 *  @brief Flag for nitrogen atom in paths.
 */
#define NITROGEN_F 5

/** @def OXYGEN_F
 *  @brief Flag for oxygen atom in paths.
 */
#define OXYGEN_F 4

/** @def HYDROGEN_F
 *  @brief Flag for hydrogen atom in paths.
 */
#define HYDROGEN_F 7

// Discretization
/** @def STEP_GRID_VOXEL
 *  @brief Step size for the grid discretization.
 */
#ifndef STEP_GRID_VOXEL
#define STEP_GRID_VOXEL 0.5
#endif
// Need to be smaller than DIST_GAP_CAGE, because we need to be able to access linkable atoms

#endif