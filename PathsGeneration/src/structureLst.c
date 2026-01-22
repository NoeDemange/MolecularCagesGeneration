#include "assembly.h"
#include "distance.h"
#include "structure.h"
#include "util.h"

#include <limits.h>

/**
 * @file structureLst.c
 * @brief List operations functions for manipulating linked lists.
 *
 * This file contains functions for creating, manipulating, and checking linked lists.
 */

/**************************************/
/* LIST *******************************/
/**************************************/
/**
 * @brief Initializes a new list.
 *
 * This function initializes a new list by setting its elements to NULL and size to 0.
 *
 * @param l Pointer to the list to initialize.
 */
void lstInit(List_t *l) {
  l->elts = NULL;
  l->size = 0;
}

/**
 * @brief Adds allocated memory to the list.
 *
 * This function adds allocated memory to the list to accommodate new elements.
 *
 * @param l Pointer to the list to which memory is added.
 */
void lstAddAlloc(List_t *l) {
  int i;

  l->elts = realloc(l->elts, (size(l) + REALLOCSIZE) * sizeof(int));

  for (i = 0; i < REALLOCSIZE; i++)
    l->elts[l->size + i] = -1;

  size(l) += REALLOCSIZE;
}

/**
 * @brief Returns the number of elements in the list.
 *
 * This function returns the number of elements in the list.
 *
 * @param l Pointer to the list.
 * @return The number of elements in the list.
 */
unsigned lstNbElements(List_t *l) {
  int cpt;

  for (cpt = size(l); cpt > 0 && elts(l, cpt - 1) == -1; cpt--)
    ;

  return cpt;
}

/**
 * @brief Returns the index of the first available free element in the list.
 *
 * This function returns the index of the first available free element in the list.
 * If no free element is available, it adds and allocates memory for more elements.
 *
 * @param l Pointer to the list.
 * @return The index of the first available free element in the list.
 */
unsigned lstGetIndiceFree(List_t *l) {
  int i;

  for (i = 0; i < size(l); i++)
    if (elts(l, i) == -1)
      return i;

  lstAddAlloc(l);
  return i;
}

/**
 * @brief Returns the index of the element with the specified identifier in the list.
 *
 * This function returns the index of the element with the specified identifier in the list.
 *
 * @param l Pointer to the list.
 * @param id Identifier of the element to search for.
 * @return The index of the element with the specified identifier, or -1 if not found.
 */
unsigned lstGetIndice(List_t *l, unsigned id) {
  int i;

  for (i = 0; i < size(l) && elts(l, i) != -1; i++)
    if (elts(l, i) == id)
      return i;

  return -1;
}

/**
 * @brief Checks if an element with the specified identifier exists in the list.
 *
 * This function checks if an element with the specified identifier exists in the list.
 *
 * @param l Pointer to the list.
 * @param id Identifier of the element to check for.
 * @return 1 if the element with the specified identifier exists, 0 otherwise.
 */
unsigned lstCheck(List_t *l, unsigned id) {

  if (lstGetIndice(l, id) == -1)
    return 0;
  return 1;
}

/**
 * @brief Adds an element to the list if it does not exist.
 *
 * This function adds an element to the list if it does not already exist.
 *
 * @param l Pointer to the list.
 * @param id Identifier of the element to be added.
 */
void lstAddElement(List_t *l, unsigned id) {

  int i;
  if (lstGetIndice(l, id) == -1) {
    i = lstGetIndiceFree(l);
    l->elts[i] = id;
  }
}

/**
 * @brief Removes an element with the specified identifier from the list.
 *
 * This function removes an element with the specified identifier from the list.
 *
 * @param l Pointer to the list.
 * @param id Identifier of the element to be removed.
 */
void lstRemoveElement(List_t *l, unsigned id) {

  int i = lstGetIndice(l, id);

  if (i != -1) {
    while (i < size(l) - 1 && elts(l, i) != -1) {
      elts(l, i) = elts(l, i + 1);
      i++;
    }
    elts(l, i) = -1;
  }
}

/**
 * @brief Creates a new integer list.
 *
 * This function allocates memory for a new integer list and initializes its fields.
 *
 * @return Pointer to the newly allocated integer list.
 */
List_t *lstCreate() {

  List_t *l = malloc(sizeof(List_t));
  lstInit(l);

  return l;
}

/**
 * @brief Creates a copy of an integer list.
 *
 * This function creates a copy of an integer list by copying its elements.
 *
 * @param l Pointer to the integer list to be copied.
 * @return Pointer to the newly created copy of the integer list.
 */
List_t *lstCopy(List_t *l) {

  int i;
  List_t *copy = lstCreate();

  size(copy) = lstNbElements(l);

  copy->elts = malloc(size(copy) * sizeof(int));

  for (i = 0; i < size(copy); i++)
    elts(copy, i) = elts(l, i);

  return copy;
}

/**
 * @brief Copies the list with an offset of the numbering of elements
 * according to the value of the array passed as an argument.
 *
 * @param l List to copy.
 * @param shifts Array of offsets.
 * @return (List_t*) Copied list with shifts.
 */
List_t *lstCopyWithShift(List_t *l, int *shifts) {

  int i;
  List_t *copy = lstCreate();

  size(copy) = lstNbElements(l);

  copy->elts = malloc(size(copy) * sizeof(int));

  for (i = 0; i < size(copy); i++)
    elts(copy, i) = elts(l, i) - shifts[elts(l, i)];

  return copy;
}

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
List_t *lstAddList(List_t *l1, List_t *l2) {

  int i;

  List_t *out = lstCopy(l1);

  for (i = 0; i < size(l2) && elts(l2, i) != -1; i++)
    lstAddElement(out, elts(l2, i));

  lstDelete(l2);
  lstDelete(l1);
  return out;
}

/**
 * @brief Deletes an integer list and frees the associated memory.
 *
 * This function deletes an integer list and frees the memory associated with its elemencurrentElemts.
 *
 * @param l Pointer to the integer list to be deleted.
 */
void lstDelete(List_t *l) {

  free(l->elts);
  free(l);
}

/******************************/

/**
 * @brief Initializes a new List_s, a list of Point_t.
 *
 * This function initializes a new List_s (used with Point_t elements) by setting its first element to NULL.
 *
 * @return Pointer to the newly allocated List_s.
 */
List_s *lstsInit() {

  List_s *list = malloc(sizeof(List_s));
  if (!list) {
    fprintf(stderr, "Memory allocation failed for List_s\n");
    exit(EXIT_FAILURE);
  }
  list->first = NULL;

  return list;
}

// Ajout au début
/**
 * @brief Adds a Point_t element at the beginning of the List_s (Point_t) list.
 *
 * This function adds a Point_t element to the List_s (Point_t) list.
 *
 * @param list Pointer to the List_s (Point_t) list.
 * @param sommet Point_t element to add.
 */
void lstsAddElement(List_s *list, Point_t sommet) {

  Elem_s *elem = malloc(sizeof(Elem_s));

  elem->position.x = sommet.x;
  elem->position.y = sommet.y;
  elem->position.z = sommet.z;
  elem->next = list->first;

  list->first = elem;
}

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
 * @param distance_type Type of distance to use (e.g., Euclidean, A*).
 */
void lstsAddElementInOrder(List_s *list, Point_t startPos, Point_t endPos, Paths_t *paths, int growth_limit,
                           const char *string_distance_type) {
  double computed_dist = dist(startPos, endPos);
  int boundary_enabled = is_path_boundary_filter_enabled();
  double boundary_limit = 0.0;
  if (boundary_enabled) {
    int current_path = paths->currentPath;
    int remaining_slots = growth_limit - paths->curPthPos[current_path]+1; //include slot at growth_limit
    if (remaining_slots < 0) {
      remaining_slots = 0;
    }
    boundary_limit = DIST_PATH_BOUNDARY(remaining_slots);
  }
  if (!boundary_enabled || computed_dist <= boundary_limit) {
    // choix entre distance euclidienne et distance A*, modify computed_dist if A*
    if (strcmp(string_distance_type, "Euclidean") != 0) {
      computed_dist =
          aStarDistance(startPos, endPos, paths->grids[paths->currentPath], paths->minHeaps[paths->currentPath]);
    }
    // printf("paths size max %d, current %d\n", paths->sizeMax, paths->curPthPos[paths->currentPath]);

    Elem_s *current_elem = list->first;
    Elem_s *previous_elem = NULL;
    while (current_elem) {
      if (current_elem->distance < computed_dist) {
        previous_elem = current_elem;
        current_elem = current_elem->next;
      } else {
        break;
      }
    }
    Elem_s *elem = malloc(sizeof(Elem_s));
    elem->position = startPos;
    elem->distance = computed_dist;
    elem->next = current_elem;
    if (previous_elem) {
      previous_elem->next = elem;
    } else {
      list->first = elem;
    }
    if (boundary_enabled) {
#ifdef ENABLE_STATS
      nb_path_boundary_allowed++;
#endif
    }
  } else {
#ifdef ENABLE_STATS
    nb_path_boundary_blocked++;
#endif
  }
}

/**
 * @brief Removes the first Point_t element from the List_s (Point_t) list.
 *
 * This function removes the first Point_t element from the List_s (Point_t) list and frees the associated memory.
 *
 * @param list Pointer to the List_s (Point_t) list.
 */
void lstsRemoveFirst(List_s *list) {

  Elem_s *suppr = list->first;
  list->first = list->first->next;
  free(suppr);
}

/**
 * @brief Deletes the List_s (Point_t) list and frees the associated memory.
 *
 * This function deletes the List_s (Point_t) list and frees the memory associated with its elements.
 *
 * @param list Pointer to the List_s (Point_t) list to be deleted.
 */
void lstsDelete(List_s *list) {

  while (list->first) {
    lstsRemoveFirst(list);
  }
  free(list);
}

/**
 * @brief Removes a Point_t element from the List_s (Point_t) list.
 *
 * This function removes a Point_t element from the List_s (Point_t) list.
 *
 * @param list Pointer to the List_s (Point_t) list.
 * @param p Point_t element to be removed.
 */
void lstsRemoveElement(List_s *list, Point_t p) {

  Elem_s *cursor = list->first;
  Elem_s *suppr = NULL;
  if (cursor) {
    if (cursor->position.x == p.x && cursor->position.y == p.y && cursor->position.z == p.z) {
      lstsRemoveFirst(list);
    } else {
      while (cursor->next && !suppr) {
        if (cursor->next->position.x == p.x && cursor->next->position.y == p.y && cursor->next->position.z == p.z) {
          suppr = cursor->next;
          cursor->next = cursor->next->next;
        } else
          cursor = cursor->next;
      }
    }
  }
  if (suppr) {
    free(suppr);
  }
}

/**
 * @brief Prints the elements of the List_s (Point_t) list.
 *
 * This function prints the elements of the List_s (Point_t) list to the standard output.
 *
 * @param list Pointer to the List_s (Point_t) list to be printed.
 */
void lstsPrint(List_s *list) {
  Elem_s *cursor = list->first;
  printf("List of elements:\n");
  while (cursor) {
    printf("x: %f, y: %f, z: %f, dist %lf\n", cursor->position.x, cursor->position.y, cursor->position.z,
           cursor->distance);
    cursor = cursor->next;
  }
}