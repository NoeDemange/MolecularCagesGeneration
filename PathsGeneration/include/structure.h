#ifndef __STRUCTURE_H
#define __STRUCTURE_H

#include "constant.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * @file structure.h
 * @brief Structure Header File
 *
 * This file contains data structures and macros used throughout the application.
 */

// MACROS FOR ACCESSING STRUCTURE ELEMENTS

/**
 * @brief Access the address of an atom in the cage structure.
 *
 * This macro returns the address of the atom at the specified index in the cage structure.
 *
 * @param o The cage structure.
 * @param i The index of the atom.
 * @return The address of the atom.
 */
#define atom(o, i) ((o)->atoms + (i))

/**
 * @brief Access the neighborhood of an atom.
 *
 * This macro returns the address of the neighborhood list of an atom in the cage structure.
 *
 * @param a The atom in the cage structure.
 * @return The address of the neighborhood list of the atom.
 */
#define neighborhood(a) (a)->neighborhood

/**
 * @brief Access the size of the cage structure.
 *
 * This macro returns the size of the cage structure.
 *
 * @param o The cage structure.
 * @return The size of the cage structure.
 */
#define size(o) (o)->size

/**
 * @def elts(l, i)
 * @brief Access the i-th element of the list l.
 */
#define elts(l, i) (l)->elts[(i)]

/**
 * @def forEachElement(l, i)
 * @brief Macro to iterate through a list until it finds an unused value (-1).
 * Only used to scroll through lists where unused values are at the end.
 */
#define forEachElement(l, i) (i) < size((l)) && elts((l), (i)) != -1

/**
 * @def coords(a)
 * @brief Access the coordinates of an atom (a).
 */
#define coords(a) (a)->coords

/**
 * @def atomX(a)
 * @brief Access the coordinate x of an atom (a).
 */
#define atomX(a) (a)->coords.x

/**
 * @def atom(a)
 * @brief Access the coordinate y of an atom (a).
 */
#define atomY(a) (a)->coords.y

/**
 * @def atomZ(a)
 * @brief Access the coordinate z of an atom (a).
 */
#define atomZ(a) (a)->coords.z

/**
 * @def neighbor(a, i)
 * @brief Access the i-th neighbor of an atom (a).
 */
#define neighbor(a, i) (a)->neighborhood->elts[(i)]

/**
 * @def neighborhoodSize(a)
 * @brief Get the size of the neighborhood of an atom (a).
 */
#define neighborhoodSize(a) (a)->neighborhood->size

/**
 * @def forEachNeighbor(a, i)
 * @brief Macro to iterate through the neighbors of an atom until it finds an unused value (-1).
 * The unused values of neighbors must be at the end.
 */
#define forEachNeighbor(a, i) (i) < neighborhoodSize((a)) && neighbor((a), (i)) != -1

#define coordsNeighbor(s, a, i) coords(atom((s), neighbor(atom((s), (a)), (i))))

/**
 * @def flag(a)
 * @brief Access the flag of a cage atom (a).
 */
#define flag(a) (a)->flag

/**
 * @def parentAtom(a)
 * @brief Access the parent atom index of a cage atom (a).
 */
#define parentAtom(a) (a)->parentAtom

/**
 * @brief Macro to access the index of a Point_t in a flattened 4D array.
 *
 * This macro computes the index of a `Point_t` element stored in the
 * 4D array `path->patterns` within the structure.
 *
 * @param k Index in the first dimension.
 * @param s Index in the second dimension.
 * @param p Index in the third dimension.
 * @param a Index in the fourth dimension.
 * @param sizeMax Size of the second dimension.
 * @return Index to the specified `Point_t` in the array.
 */
#define indexPointPaths(k, s, p, a, sizeMax)                                                                           \
  (((k) * (sizeMax) * (NUMBER_POSITION_PATHS) * (MAX_NB_ATOMS_PATTERN)) +                                              \
   ((s) * (NUMBER_POSITION_PATHS) * (MAX_NB_ATOMS_PATTERN)) + ((p) * (MAX_NB_ATOMS_PATTERN)) + (a))

/**
 * @brief Macro to access a path positions in the Paths.
 *
 * This macro computes the index to an `int` element in a 2D array,
 *  which is stored as a flattened 1D array.
 *
 * @param k Index in the first dimension, number of paths.
 * @param s Index in the second dimension, number of path position.
 * @param sizeMax Size of the second dimension (`paths->sizeMax`).
 * @return Pointer to the specified `int` in the array.
 */
#define indexPathPosition(k, s, sizeMax) ((k) * (sizeMax) + (s))
/**************************************/
/* POINT ******************************/
/**************************************/

/**
 * @struct Point_t
 * @brief Represents a 3D point with floating-point coordinates.
 */
typedef struct {
  double x; /**< The x-coordinate of the point. */
  double y; /**< The y-coordinate of the point. */
  double z; /**< The z-coordinate of the point. */
} Point_t;

/**************************************/
/* LISTE ******************************/
/**************************************/
/**
 * @struct List_t
 * @brief Represents a list of integers.
 */
typedef struct {
  int *elts;     /**< An array of integers representing the elements in the list. */
  unsigned size; /**< The size of the list. */
} List_t;

/**
 * @struct Elem_s
 * @brief Represents an element in the list of intermediate vertices.
 */
typedef struct Elem_s Elem_s;
struct Elem_s {
  Point_t position; /**< The position of the intermediate vertex. */
  double distance;  /**< The distance of the intermediate vertex from the starting point. */
  Elem_s *next;     /**< Pointer to the next element in the list. */
};

/**
 * @struct List_s
 * @brief Represents a list of intermediate vertices.
 */
typedef struct {
  Elem_s *first; /**< Pointer to the first element in the list. */
} List_s;

/**************************************/
/* CAGE ******************************/
/**************************************/
/**
 * @struct AtomCage_t
 * @brief Represents an atom in a cage.
 *
 * This structure represents an atom within a cage, including its flag, coordinates,
 * parent atom index, and neighborhood.
 */
typedef struct {
  int flag;             /**< Flag representing the state of the atom. */
  Point_t coords;       /**< Coordinates of the atom. */
  List_t *neighborhood; /**< List of neighbors of the atom in the cage. */
} AtomCage_t;

/**
 * @struct Cage_t
 * @brief Represents a cage.
 *
 * This structure holds information about a cage, including its atoms, cycles, bonds, and size.
 */
typedef struct {
  AtomCage_t *atoms; /**< Array of atoms in the cage. */
  unsigned size;     /**< Size of the cage (number of atoms). */
  double xMinWOGap;  /**< The x min of the cage (min x vertex). */
  double yMinWOGap;  /**< The y min of the cage (min y vertex). */
  double zMinWOGap;  /**< The z min of the cage (min z vertex). */
  double xMaxWOGap;  /**< The x max of the cage (max x vertex). */
  double yMaxWOGap;  /**< The y max of the cage (max y vertex). */
  double zMaxWOGap;  /**< The z max of the cage (max z vertex). */
} Cage_t;

/**************************************/
/* GRID *******************************/
/**************************************/

/**
 * @struct Node
 * @brief Represents a node in a 3D grid.
 *
 * This structure represents a node in a 3D grid used for voxelization, including its coordinates,
 * walkability status, and values for the A* algorithm.
 *
 * @note The `walkable` field indicates whether the node is walkable (1) or blocked (0).
 * The `closed` and `opened` fields are used for the A* algorithm to track the status of the node.
 * The `gCost`, `hCost`, and `fCost` fields are used to store the costs associated with the node
 * during the A* pathfinding process.
 */
typedef struct Node Node;
struct Node {
  int x, y, z;  /**< Coordinates in the grid */
  int walkable; /**< 1 = walkable, 0 = blocked */
  // Values for A* algorithm
  int closed;                 /**< 1 = closed, 0 = open for A* algorithm */
  int opened;                 /**< 1 = opened, 0 = closed for A* algorithm */
  double gCost;               /**< The g-value of the node. Not init for the moment.*/
  double hCost;               /**< The h-value of the node. Not init for the moment.*/
  double fCost;               /**< The f-value of the node. Not init for the moment.*/
  Node *parent;               /**< Pointer to the parent node in the path. */
  int heapIndex;              /**< Index in the min-heap for A* algorithm. */
  int nbTimesCandidate;       /**< Number of times this node is a candidate. */
  int *indexStartsCandidates; /**< Index of the candidate node. */
};

/**
 * @struct Grid_t
 * @brief Represents a 3D grid for voxelization.
 *
 * This structure represents a 3D grid used for voxelization, including its dimensions,
 * offset, and nodes.
 */
typedef struct {
  int width, height, depth;         /**< Total voxel count per axis */
  int offset_x, offset_y, offset_z; /**< Offset to center the grid */
  Node ***nodes;                    /**< 3D array of nodes in the grid */
  Node **visitedNodes;              /**< Array of pointers to visited nodes */
  int numVisited;                   /**< Number of currently visited nodes */
  int maxVisited;                   /**< Maximum capacity of visitedNodes array */
} Grid_t;

/**************************************/
/* MINHEAP ****************************/
/**************************************/

/**
 * @struct MinHeap_t
 * @brief Represents a min-heap.
 *
 * This structure represents a min-heap used for the A* algorithm,
 * including its size, maximum size, and an array of nodes.
 */
typedef struct {
  int size;     /**< Number of elements in the heap */
  int maxSize;  /**< Maximum size of the heap */
  Node **nodes; /**< Array of nodes in the heap */
} MinHeap_t;

/**************************************/
/* PATHS ******************************/
/**************************************/

/**
 * @struct Paths_t
 * @brief Represents a set of paths for cage creation.
 *
 * This structure holds information about the paths used for cage creation,
 * including the number of paths, maximum size, and an array of patterns.
 */
typedef struct {
  // int idStart;
  // int idEnd;
  // ADD interconnectionTree+numComponent?
  int currentPath;   /**< path being processed. */
  int numPaths;      /**< Number of paths to have a cage */
  int sizeMax;       /**< Number max of patterns in a path */
  Point_t *patterns; /**< Array of patterns. */
  // Point_t* positionsBuffer; Not use
  //  Point_t* solution; Not use
  // int* orientations; // Cycle orientation.
  int *curPthPos;                   /**< Curent path position in current path. Array size numPaths. */
  int *positionCurNum;              /**< Current pattern position for each path position for each path. */
  int *patternCurNum;               /**< Current pattern number for each path position for each path. */
  int *maxPositions;                /**< number of saved pattern positions for each path position for each path. */
  Grid_t **grids;                   /**< Array of grids, one for each step of path creation. */
  MinHeap_t **minHeaps;             /**< Array of min-heaps, one for each step of path creation. */
  Point_t *starts;                  /**< Array of starting points used each times we add atoms in the path. */
  Point_t *results_pos;             /**< Array of results positions. */
  Node **Candidates;                /**< Array of candidate node addresses for the SSMTA* algorithm. */
  double *distancesMultiCandidates; /**< Array of distances for the SSMTA* algorithm. */
  double *pathMSD;              /**< Array of MSD (Mean Squared Deviation) for each path. Not initialized, just allocated. */
  int *bestPathLength;          /**< Best number of real patterns (without start/start neighbor) found per path, -1 if none. */
  int *maxGrowthLimit;          /**< Current depth limit (curPthPos) allowed for each path. */
} Paths_t;

// DEFINITION
// Point
/**
 * @brief Initializes a new
Point_t with equal scalar values.
 *
 * This function initializes a new Point_t with the same scalar value for all coordinates (x, y, z).
 *
 * @param scal The scalar value to set for all coordinates of the Point_t.
 * @return The initialized Point_t.
 */
Point_t ptInit(double);

/**
 * @brief Adds two Point_t together and returns the result.
 *
 * This function adds two Point_t (A and B) together and returns the resulting Point_t.
 *
 * @param A The first Point_t operand.
 * @param B The second Point_t operand.
 * @return The result of the addition operation as a new Point_t.
 */
Point_t ptAdd(Point_t, Point_t);

/**
 * @brief Subtracts one Point_t from another and returns the result.
 *
 * This function subtracts Point_t B from Point_t A and returns the resulting Point_t.
 *
 * @param A The first Point_t operand.
 * @param B The second Point_t operand.
 * @return The result of the subtraction operation as a new Point_t.
 */
Point_t ptSub(Point_t, Point_t);

/**
 * @brief Multiplies a Point_t by a scalar value and returns the result.
 *
 * This function multiplies each coordinate of the Point_t A by the scalar value "scal" and returns the resulting
 * Point_t.
 *
 * @param A The Point_t operand.
 * @param scal The scalar value to multiply each coordinate of the Point_t A.
 * @return The result of the multiplication operation as a new Point_t.
 */
Point_t ptMul(Point_t, double);

/**
 * @brief Divides a Point_t by a scalar value and returns the result.
 *
 * This function divides each coordinate of the Point_t A by the scalar value "scal" and returns the resulting Point_t.
 * If the scalar value is 0, it returns a Point_t with all coordinates set to 0.
 *
 * @param A The Point_t operand.
 * @param scal The scalar value to divide each coordinate of the Point_t A.
 * @return The result of the division operation as a new Point_t.
 */
Point_t ptDiv(Point_t, double);
int ptCompare(Point_t A, Point_t B);

/**
 * @brief Merges two Point_t by taking their average and returns the result.
 *
 * This function takes the average of Point_t A and Point_t B (by adding them and dividing by 2) and returns the
 * resulting Point_t.
 *
 * @param A The first Point_t operand.
 * @param B The second Point_t operand.
 * @return The result of the merging operation as a new Point_t.
 */
Point_t ptMerge(Point_t A, Point_t B);

/**
 * @brief Checks if two Point_t are equal.
 *
 * This function checks if two Point_t (A and B) are equal by comparing their x, y, and z coordinates.
 * If they are equal, the function returns 1; otherwise, it returns 0.
 *
 * @param A The first Point_t operand.
 * @param B The second Point_t operand.
 * @return 1 if Point_t A and Point_t B are equal; otherwise, 0.
 */
int ptEqual(Point_t A, Point_t B);

// Liste
/**
 * @brief Initializes a new list.
 *
 * This function initializes a new list by setting its elements to NULL and size to 0.
 *
 * @param l Pointer to the list to initialize.
 */
void lstInit(List_t *);

/**
 * @brief Adds allocated memory to the list.
 *
 * This function adds allocated memory to the list to accommodate new elements.
 *
 * @param l Pointer to the list to which memory is added.
 */
void lstAddAlloc(List_t *);

/**
 * @brief Returns the number of elements in the list.
 *
 * This function returns the number of elements in the list.
 *
 * @param l Pointer to the list.
 * @return The number of elements in the list.
 */
unsigned lstNbElements(List_t *);

/**
 * @brief Returns the index of the first available free element in the list.
 *
 * This function returns the index of the first available free element in the list.
 * If no free element is available, it adds and allocates memory for more elements.
 *
 * @param l Pointer to the list.
 * @return The index of the first available free element in the list.
 */
unsigned lstGetIndiceFree(List_t *);

/**
 * @brief Returns the index of the element with the specified identifier in the list.
 *
 * This function returns the index of the element with the specified identifier in the list.
 *
 * @param l Pointer to the list.
 * @param id Identifier of the element to search for.
 * @return The index of the element with the specified identifier, or -1 if not found.
 */
unsigned lstGetIndice(List_t *, unsigned);

/**
 * @brief Checks if an element with the specified identifier exists in the list.
 *
 * This function checks if an element with the specified identifier exists in the list.
 *
 * @param l Pointer to the list.
 * @param id Identifier of the element to check for.
 * @return 1 if the element with the specified identifier exists, 0 otherwise.
 */
unsigned lstCheck(List_t *, unsigned);

/**
 * @brief Adds an element to the list if it does not exist.
 *
 * This function adds an element to the list if it does not already exist.
 *
 * @param l Pointer to the list.
 * @param id Identifier of the element to be added.
 */
void lstAddElement(List_t *, unsigned);

/**
 * @brief Removes an element with the specified identifier from the list.
 *
 * This function removes an element with the specified identifier from the list.
 *
 * @param l Pointer to the list.
 * @param id Identifier of the element to be removed.
 */
void lstRemoveElement(List_t *, unsigned);

/**
 * @brief Creates a new integer list.
 *
 * This function allocates memory for a new integer list and initializes its fields.
 *
 * @return Pointer to the newly allocated integer list.
 */
List_t *lstCreate();

/**
 * @brief Creates a copy of an integer list.
 *
 * This function creates a copy of an integer list by copying its elements.
 *
 * @param l Pointer to the integer list to be copied.
 * @return Pointer to the newly created copy of the integer list.
 */
List_t *lstCopy(List_t *);

/**
 * @brief Copies the list with an offset of the numbering of elements
 * according to the value of the array passed as an argument.
 *
 * @param l List to copy.
 * @param shifts Array of offsets.
 * @return (List_t*) Copied list with shifts.
 */
List_t *lstCopyWithShift(List_t *l, int *mod_pos_nei);

/**
 * @brief Merges two integer lists.
 *
 * This function merges two integer lists, adding elements from the second list to the first.
 * The function then deletes both input lists and returns the merged list.
 *
 * @param l1 First integer list.
 * @param l2 Second integer list.
 * @return Merged integer list.
 */
List_t *lstAddList(List_t *, List_t *);

/**
 * @brief Deletes an integer list and frees the associated memory.
 *
 * This function deletes an integer list and frees the memory associated with its elements.
 *
 * @param l Pointer to the integer list to be deleted.
 */
void lstDelete(List_t *);

/**
 * @brief Initializes a new List_s (Point_t) list.
 *
 * This function initializes a new List_s (used with Point_t elements) by setting its first element to NULL.
 *
 * @return Pointer to the newly allocated List_s.
 */
List_s *lstsInit();

/**
 * @brief Adds a Point_t element at the beginning of the List_s (Point_t) list.
 *
 * This function adds a Point_t element to the List_s (Point_t) list.
 *
 * @param list Pointer to the List_s (Point_t) list.
 * @param sommet Point_t element to add.
 */
void lstsAddElement(List_s *list, Point_t sommet);

/**
 * @brief Adds a Point_t element to the List_s (Point_t) list in ascending order of A* distance.
 *
 * This function adds a Point_t element to the List_s (Point_t) list in ascending order of A* distance.
 *
 * @param list Pointer to the List_s (Point_t) list.
 * @param startPos Position of the atom.
 * @param endPos Position of the objective atom.
 * @param voxelGrid Grid of voxelization.
 * @param paths Pointer to the Paths_t structure containing the grids and minHeaps.
 * @param growth_limit Current dynamic cap (curPthPos index) allowed for the active path.
 * @param string_distance_type Type of distance to use (e.g., Euclidean, A*).
 */
void lstsAddElementInOrder(List_s *list, Point_t startPos, Point_t endPos, Paths_t *paths, int growth_limit,
                           const char *string_distance_type);

/**
 * @brief Removes the first Point_t element from the List_s (Point_t) list.
 *
 * This function removes the first Point_t element from the List_s (Point_t) list and frees the associated memory.
 *
 * @param list Pointer to the List_s (Point_t) list.
 */
void lstsRemoveFirst(List_s *list);

/**
 * @brief Deletes the List_s (Point_t) list and frees the associated memory.
 *
 * This function deletes the List_s (Point_t) list and frees the memory associated with its elements.
 *
 * @param list Pointer to the List_s (Point_t) list to be deleted.
 */
void lstsDelete(List_s *list);

/**
 * @brief Removes a Point_t element from the List_s (Point_t) list.
 *
 * This function removes a Point_t element from the List_s (Point_t) list.
 *
 * @param list Pointer to the List_s (Point_t) list.
 * @param p Point_t element to be removed.
 */
void lstsRemoveElement(List_s *list, Point_t p);

/**
 * @brief Prints the elements of the List_s (Point_t) list.
 *
 * This function prints the elements of the List_s (Point_t) list to the standard output.
 *
 * @param list Pointer to the List_s (Point_t) list to be printed.
 */
void lstsPrint(List_s *list);

// Cage

/**
 * @brief Gets the number of neighbors in the neighborhood of an AtomCage_t.
 *
 * This function returns the number of neighbors in the neighborhood of the given AtomCage_t.
 *
 * @param a Pointer to the AtomCage_t.
 * @return The number of neighbors in the AtomCage_t's neighborhood.
 */
int cageNbNeighborhood(AtomCage_t *);

/**
 * @brief Get the number of atoms in the Cage structure.
 *
 * This function returns the number of atoms in the specified Cage_t structure.
 *
 * @param s Pointer to the Cage_t structure.
 * @return The number of atoms in the Cage structure.
 */
int cageNbAtom(Cage_t *);

/**
 * @brief Get the number of edges in the Cage structure.
 *
 * This function returns the number of edges (connections between atoms) in the specified Cage_t structure.
 *
 * @param s Pointer to the Cage_t structure.
 * @return The number of edges in the Cage structure.
 */
int cageNbEdges(Cage_t *);

/**
 * @brief Adds an edge (bond) between two vertices with the specified IDs.
 *
 * This function adds an edge (bond) between two vertices (atoms) with the given IDs to the Cage data structure.
 * It also adds the reverse edge to create an undirected graph representation.
 *
 * @param s Pointer to the Cage data structure.
 * @param id1 The ID of the first vertex (atom).
 * @param id2 The ID of the second vertex (atom).
 */
void cageAddEdge(Cage_t *, unsigned, unsigned);

/**
 * @brief Removes an edge (bond) between two vertices with the specified IDs.
 *
 * This function removes an edge (bond) between two vertices (atoms) with the given IDs from the Cage data structure.
 * It also removes the reverse edge to maintain the undirected graph representation.
 *
 * @param s Pointer to the Cage data structure.
 * @param id1 The ID of the first vertex (atom).
 * @param id2 The ID of the second vertex (atom).
 */
void cageRemoveEdge(Cage_t *, unsigned, unsigned);

/**
 * @brief Adds a new atom to the Cage with the specified coordinates and parent atom ID.
 *
 * This function adds a new atom to the Cage data structure with the specified coordinates and parent atom ID.
 * It also assigns a unique ID (indice) to the new atom.
 *
 * @param s Pointer to the Cage data structure.
 * @param coords The coordinates (Point_t) of the new atom.
 * @param parent The ID of the parent atom for the new atom.
 * @return The unique ID (indice) assigned to the newly added atom.
 */
unsigned cageAddAtom(Cage_t *, Point_t); //, unsigned);

/**
 * @brief Removes an atom from the Cage with the specified ID.
 *
 * This function removes an atom with the given ID from the Cage data structure.
 * It also updates the neighborhood and other properties accordingly.
 *
 * @param s Pointer to the Cage data structure.
 * @param id The ID of the atom to be removed.
 */
void cageRemoveAtom(Cage_t *, unsigned);

/**
 * @brief Adds a vertex to the Cage with the specified ID.
 *
 * This function adds a vertex with the given ID to the Cage data structure.
 *
 * @param s Pointer to the Cage data structure.
 * @param id The ID of the vertex to be added.
 * @return The ID of the added vertex.
 */
unsigned cageAddVertex(Cage_t *, unsigned);

/**
 * @brief Removes a vertex from the Cage with the specified ID.
 *
 * This function removes a vertex with the given ID from the Cage data structure.
 *
 * @param s Pointer to the Cage data structure.
 * @param id The ID of the vertex to be removed.
 */
void cageRemoveVertex(Cage_t *, unsigned);

/**
 * @brief Adds a bond between two vertices with the specified IDs.
 *
 * This function adds a bond between two vertices (atoms) with the given IDs to the Cage data structure.
 *
 * @param s Pointer to the Cage data structure.
 * @param id1 The ID of the first vertex (atom).
 * @param id2 The ID of the second vertex (atom).
 */
void cageAddBond(Cage_t *, unsigned, unsigned);

/**
 * @brief Removes a bond between two vertices with the specified IDs.
 *
 * This function removes a bond between two vertices (atoms) with the given IDs from the Cage data structure.
 *
 * @param s Pointer to the Cage data structure.
 * @param id1 The ID of the first vertex (atom).
 * @param id2 The ID of the second vertex (atom).
 */
void cageRemoveBond(Cage_t *, unsigned, unsigned);

/**
 * @brief Adds a cycle to the Cage with the specified atom ID.
 *
 * This function adds a cycle to the Cage data structure with the specified atom ID.
 *
 * @param s Pointer to the Cage data structure.
 * @param id The ID of the atom to be added to the cycle.
 */
void cageAddCycle(Cage_t *, unsigned);

/**
 * @brief Merges two atoms in the Cage.
 *
 * This function merges (combines) two atoms in the Cage data structure.
 * It moves all the edges and properties from the eaten atom to the eater atom,
 * and then removes the eaten atom from the Cage.
 *
 * @param s Pointer to the Cage data structure.
 * @param eater The ID of the atom that will "eat" (absorb) the other atom.
 * @param eaten The ID of the atom that will be "eaten" (absorbed) by the other atom.
 */
void cageMergeAtom(Cage_t *, unsigned, unsigned);
// void CAGE_testDis(Cage_t*);

/**
 * @brief Creates a new Cage data structure.
 *
 * This function creates and initializes a new Cage data structure.
 *
 * @return A pointer to the newly created Cage data structure.
 */
Cage_t *cageCreate();

/**
 * @brief Creates a deep copy of a Cage data structure.
 *
 * This function creates a deep copy of the given Cage data structure.
 *
 * @param s Pointer to the original Cage data structure to be copied.
 * @return A pointer to the newly created deep copy of the Cage data structure.
 */
Cage_t *cageCopy(Cage_t *);

/**
 * @brief Creates a trimmed copy of a Cage data structure.
 *
 * This function creates a trimmed copy of the given Cage data structure by removing unused atoms or atoms belonging to
 * the envelope.
 *
 * @param s Pointer to the original Cage data structure containing unused atoms.
 * @return A pointer to the newly created trimmed copy of the Cage data structure.
 */
Cage_t *cageCopyCageAtoms(Cage_t *s);

/**
 * @brief Deletes a Cage data structure.
 *
 * This function deletes the entire Cage data structure, including all atoms and the graph representation.
 *
 * @param s Pointer to the Cage data structure to be deleted.
 */
void cageDelete(Cage_t *);

/**
 * @brief Deletes an atom in the Cage.
 *
 * This function deletes an atom (node) in the Cage data structure.
 *
 * @param a Pointer to the atom to be deleted.
 */
void cageDeleteAtom(AtomCage_t *a);

Cage_t *cageImport(char *inputname, char *mocNum);

// Path

Paths_t *pthCreate(int size, int numComponents);
void pthDelete(Paths_t *paths);
void pthInit(Paths_t *paths, int *interTree, Cage_t *cage);
void pthPrintOrWritePaths(Paths_t *paths, FILE *file);
void pthPrintOrWritePathTables(Paths_t *paths, FILE *file);
void pthPrintOrWriteAll(Paths_t *paths, const char *filename);
/**
 * @brief Reboots the paths structure by resetting the current path position and other related arrays.
 *
 * This function resets the current path position for each path in the Paths_t structure, allowing for a fresh start
 * in path processing. It also resets the position and pattern numbers for each path position.
 *
 * @param paths Pointer to the Paths_t structure to be rebooted.
 */
void pthReboot(Paths_t *paths);
// int PTH_countAroRings(Path_t*);
// void PTH_addPath(Shell_t* moc, Path_t* path);
// int PTH_isHindered(Path_t* path, Point_t end, Point_t endNeighbor);

// minHeap

/**
 * @brief Initialize a new MinHeap.
 *
 * This function initializes a new MinHeap with the specified maximum size.
 *
 * @param heap Pointer to the MinHeap to be initialized.
 * @param maxSize The maximum size of the MinHeap.
 */
void initMinHeap(MinHeap_t *heap, int maxSize);

/**
 * @brief Clear the MinHeap.
 *
 * This function clears the MinHeap by resetting its size to 0.
 *
 * @param heap Pointer to the MinHeap to be cleared.
 */
void clearMinHeap(MinHeap_t *heap);

/**
 * @brief Insert a node into the MinHeap.
 *
 * This function inserts a new node into the MinHeap while maintaining the MinHeap property.
 *
 * @param heap Pointer to the MinHeap.
 * @param node Pointer to the node to be inserted.
 */
void insertMinHeap(MinHeap_t *heap, Node *node);

/**
 * @brief Extract the minimum node from the MinHeap.
 *
 * This function extracts the minimum node from the MinHeap and maintains the MinHeap property.
 *
 * @param heap Pointer to the MinHeap.
 * @return Pointer to the extracted minimum node.
 */
Node *extractMin(MinHeap_t *heap);

/**
 * @brief Decrease the priority of a node in the MinHeap.
 *
 * This function decreases the `fCost` of a node in the MinHeap and adjusts its position
 * to maintain the MinHeap property.
 *
 * @param heap Pointer to the MinHeap.
 * @param index The index of the node whose priority is to be decreased.
 * @param newFCost The new `fCost` value for the node.
 */
void decreaseKeyMinHeap(MinHeap_t *heap, Node *node, double newGCost);

/**
 * @brief Renew the MinHeap with new candidates.
 *
 * This function renews the MinHeap with new candidates by recalculating the hCost
 * for each node in the heap based on the minimum distance to a candidate that is not closed.
 *
 * @param heap Pointer to the MinHeap to be renewed.
 * @param candidates Array of candidate nodes.
 * @param nb_candidates Number of candidates in the array.
 */
void renewMinHeap(MinHeap_t *heap, Node **candidates, int nb_candidates);

/**
 * @brief Fast version of renewMinHeap for cases with many candidates.
 *
 * This function provides an even faster version when there are many candidates
 * by using spatial optimization and candidate filtering.
 *
 * @param heap Pointer to the MinHeap to be renewed.
 * @param candidates Array of candidate nodes.
 * @param nb_candidates Number of candidates in the array.
 */
void renewMinHeapFast(MinHeap_t *heap, Node **candidates, int nb_candidates);

void printMinHeap(MinHeap_t *heap);

/**
 * @brief Free the memory allocated for the MinHeap.
 *
 * This function frees the memory allocated for the MinHeap and its nodes.
 *
 * @param heap Pointer to the MinHeap to be freed.
 */
void freeMinHeap(MinHeap_t *heap);

#endif