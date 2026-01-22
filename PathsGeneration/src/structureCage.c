#include "output.h"
#include "structure.h"
#include "util.h"

#include <math.h>

/**
 * @file structureCage.c
 * @brief Cage Data Structure Implementation
 *
 * This file contains the implementation of the Cage data structure and related operations.
 * The Cage data structure represents a cage or a molecular cage with atoms and their connections.
 */

/**************************************/
/* SHELL ******************************/
/**************************************/

/**
 * @brief Initializes an AtomCage_t with default values.
 *
 * This function initializes an AtomCage_t with default values for its attributes.
 *
 * @param a Pointer to the AtomCage_t to be initialized.
 */
void cageInitAtom(AtomCage_t *a) {

  flag(a) = NOT_DEF_F;

  atomX(a) = 0;
  atomY(a) = 0;
  atomZ(a) = 0;

  // parentAtom(a) = -1;

  neighborhood(a) = lstCreate();
}

/**
 * @brief Gets the number of neighbors in the neighborhood of an AtomCage_t.
 *
 * This function returns the number of neighbors in the neighborhood of the given AtomCage_t.
 *
 * @param a Pointer to the AtomCage_t.
 * @return The number of neighbors in the AtomCage_t's neighborhood.
 */
int cageNbNeighborhood(AtomCage_t *a) { return lstNbElements(neighborhood(a)); }

/**
 * @brief Gets the index of a free neighbor in the neighborhood of an AtomCage_t.
 *
 * This function returns the index of a free neighbor in the neighborhood of the given AtomCage_t.
 *
 * @param a Pointer to the AtomCage_t.
 * @return The index of a free neighbor in the AtomCage_t's neighborhood.
 */
int cageGetIndiceFreeNeighbor(AtomCage_t *a) { return lstGetIndiceFree(a->neighborhood); }

/**
 * @brief Get the index of an atom in the neighborhood list.
 *
 * This function retrieves the index of an atom with the given ID in the neighborhood list
 * of the specified AtomCage_t structure.
 *
 * @param a Pointer to the AtomCage_t structure.
 * @param id The ID of the atom to search for in the neighborhood list.
 * @return The index of the atom in the neighborhood list, or -1 if not found.
 */
int cageGetIndice(AtomCage_t *a, unsigned id) { return lstGetIndice(a->neighborhood, id); }

/**
 * @brief Add a neighbor to the atom's neighborhood.
 *
 * This function adds a neighbor with the given ID to the neighborhood list of the specified AtomCage_t structure.
 *
 * @param a Pointer to the AtomCage_t structure.
 * @param id The ID of the neighbor atom to add to the neighborhood list.
 */
void cageAddNeighbor(AtomCage_t *a, unsigned id) { lstAddElement(a->neighborhood, id); }

/**
 * @brief Remove a neighbor from the atom's neighborhood.
 *
 * This function removes a neighbor with the given ID from the neighborhood list of the specified AtomCage_t structure.
 *
 * @param a Pointer to the AtomCage_t structure.
 * @param id The ID of the neighbor atom to remove from the neighborhood list.
 */
void cageRemoveNeighbor(AtomCage_t *a, unsigned id) { lstRemoveElement(a->neighborhood, id); }

/**
 * @brief Allocate and add new atoms to the Cage structure.
 *
 * This function allocates memory for new AtomCage_t structures and adds them to the Cage structure.
 * It dynamically reallocates memory for the atoms array in the Cage structure to accommodate the new atoms.
 * If memory reallocation fails, the function prints an error message and exits the program.
 *
 * @param s Pointer to the Cage_t structure.
 */
void cageAddAllocAtom(Cage_t *s) {

  AtomCage_t *tmp = realloc(s->atoms, (size(s) + REALLOCSIZE) * sizeof(AtomCage_t));
  if (tmp == NULL) {
    fprintf(stderr, "A problem occurred during the reallocation (structureCage.c:115).\n");
    exit(EXIT_FAILURE);
  } else {
    s->atoms = tmp;
  }

  for (int i = 0; i < REALLOCSIZE; i++) {
    cageInitAtom(atom(s, size(s) + i));
  }

  size(s) += REALLOCSIZE;
}

/**
 * @brief Get the number of atoms in the Cage structure.
 *
 * This function returns the number of atoms in the specified Cage_t structure.
 *
 * @param s Pointer to the Cage_t structure.
 * @return The number of atoms in the Cage structure.
 */
int cageNbAtom(Cage_t *s) {
  int i, cpt = 0;

  for (i = 0; i < size(s); i++)
    if (flag(atom(s, i)) != NOT_DEF_F)
      cpt++;

  return cpt;
}

/**
 * @brief Get the number of edges in the Cage structure.
 *
 * This function returns the number of edges (connections between atoms) in the specified Cage_t structure.
 *
 * @param s Pointer to the Cage_t structure.
 * @return The number of edges in the Cage structure.
 */
int cageNbEdges(Cage_t *s) {
  int i, cpt = 0;

  for (i = 0; i < size(s); i++)
    cpt += cageNbNeighborhood(atom(s, i));

  return cpt / 2;
}

/**
 * @brief Get the index of a free atom in the Cage structure.
 *
 * This function searches for a free (not defined) atom in the specified Cage_t structure and returns its index.
 * If no free atom is found, it allocates new atoms using CAGE_addAllocAtom and returns the index of the first newly
 * added atom.
 *
 * @param s Pointer to the Cage_t structure.
 * @return The index of the free atom in the Cage structure.
 */
int cageGetIndiceFreeAtom(Cage_t *s) {
  int i;

  for (i = 0; i < size(s); i++)
    if (flag(atom(s, i)) == NOT_DEF_F)
      return i;

  cageAddAllocAtom(s);
  return i;
}

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
void cageAddEdge(Cage_t *s, unsigned id1, unsigned id2) {

  if (id1 < size(s) && id2 < size(s) && id1 != id2) {

    cageAddNeighbor(atom(s, id1), id2);
    cageAddNeighbor(atom(s, id2), id1);
  }
}

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
void cageRemoveEdge(Cage_t *s, unsigned id1, unsigned id2) {

  if (id1 < size(s) && id2 < size(s)) {

    cageRemoveNeighbor(atom(s, id1), id2);
    cageRemoveNeighbor(atom(s, id2), id1);
  }
}

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
unsigned cageAddAtom(Cage_t *s, Point_t coords) { //, unsigned parent) {

  unsigned indice = cageGetIndiceFreeAtom(s);

  flag(atom(s, indice)) = SHELL_F;
  coords(atom(s, indice)) = coords;
  // parentAtom(atom(s,indice)) = parent;

  return indice;
}

/**
 * @brief Removes an atom from the Cage with the specified ID.
 *
 * This function removes an atom with the given ID from the Cage data structure.
 * It also updates the neighborhood and other properties accordingly.
 *
 * @param s Pointer to the Cage data structure.
 * @param id The ID of the atom to be removed.
 */
void cageRemoveAtom(Cage_t *s, unsigned id) {

  int i;

  if (id < size(s)) {
    AtomCage_t *a = atom(s, id);

    for (i = 0; i < neighborhoodSize(a); i++)
      if (neighbor(a, i) != -1)
        cageRemoveNeighbor(atom(s, neighbor(a, i)), id);

    lstDelete(neighborhood(a));

    cageInitAtom(a);
  }
}

/**
 * @brief Seeks the border atoms in the Cage for a given atom ID.
 *
 * This function seeks the border atoms in the Cage for a given atom ID and returns them as a list.
 *
 * @param s Pointer to the Cage data structure.
 * @param in Pointer to the input list.
 * @param id The ID of the atom to start seeking from.
 * @return A list of border atoms (border of the Cage) starting from the specified atom ID.
 */
List_t *cageSeekBorder(Cage_t *s, List_t *in, unsigned id) {

  int i;
  List_t *out = lstCreate();
  AtomCage_t *a = atom(s, id);

  if (flag(a) == SHELL_F || flag(a) == LINKABLE_F) {
    lstAddElement(out, id);
    return out;
  }

  lstAddElement(in, id);
  for (i = 0; i < neighborhoodSize(a) && neighbor(a, i) != -1; i++) {
    if (!lstCheck(in, neighbor(a, i)))
      out = lstAddList(out, cageSeekBorder(s, in, neighbor(a, i)));
  }

  return out;
}

/**
 * @brief Creates a new Cage data structure.
 *
 * This function creates and initializes a new Cage data structure.
 *
 * @return A pointer to the newly created Cage data structure.
 */
Cage_t *cageCreate() {

  Cage_t *a = malloc(sizeof(Cage_t));

  a->size = 0;
  a->atoms = NULL;

  return a;
}

/**
 * @brief Creates a deep copy of a Cage data structure.
 *
 * This function creates a deep copy of the given Cage data structure.
 *
 * @param s Pointer to the original Cage data structure to be copied.
 * @return A pointer to the newly created deep copy of the Cage data structure.
 */
Cage_t *cageCopy(Cage_t *s) {

  int i;
  Cage_t *copy = malloc(sizeof(Cage_t));

  size(copy) = size(s);
  copy->atoms = malloc(size(copy) * sizeof(AtomCage_t));

  for (i = 0; i < size(s); i++) {

    flag(atom(copy, i)) = flag(atom(s, i));
    coords(atom(copy, i)) = coords(atom(s, i));
    // parentAtom(atom(copy,i)) = parentAtom(atom(s,i));
    neighborhood(atom(copy, i)) = lstCopy(neighborhood(atom(s, i)));
  }

  return copy;
}

/**
 * @brief Creates a trimmed copy of a Cage data structure.
 *
 * This function creates a trimmed copy of the given Cage data structure by removing unused atoms or atoms belonging to
 * the envelope.
 *
 * @param s Pointer to the original Cage data structure containing unused atoms.
 * @return A pointer to the newly created trimmed copy of the Cage data structure.
 */
Cage_t *cageCopyCageAtoms(Cage_t *s) {

  int not_defs_counter = 0;
  int index = 0;
  Cage_t *copy = malloc(sizeof(Cage_t));
  int *relative_empty_positions = malloc(size(s) * sizeof(int));

  // Remove the envelope's atoms (change them to unused).
  for (int j = 0; j < size(s); j++) {
    if (flag(atom(s, j)) == SHELL_F) {
      cageRemoveAtom(s, j);
    }
  }

  size(copy) = cageNbAtom(s);
  copy->atoms = malloc(size(copy) * sizeof(AtomCage_t));

  // Count the offset in position for each kept atom (to copy the list of neighbors).
  for (int i = 0; i < size(s); i++) {
    if (flag(atom(s, i)) == NOT_DEF_F) {
      not_defs_counter++;
    } else {
      relative_empty_positions[i] = not_defs_counter;
      index++;
    }
  }

  // Copy only used atoms.
  index = 0;
  for (int i = 0; i < size(s); i++) {
    if ((flag(atom(s, i)) != NOT_DEF_F)) {
      flag(atom(copy, index)) = flag(atom(s, i));
      coords(atom(copy, index)) = coords(atom(s, i));
      // parentAtom(atom(copy,index)) = parentAtom(atom(s,i));
      neighborhood(atom(copy, index)) = lstCopyWithShift(neighborhood(atom(s, i)), relative_empty_positions);
      index++;
    }
  }
  free(relative_empty_positions);
  return copy;
}

/**
 * @brief Deletes an atom in the Cage.
 *
 * This function deletes an atom (node) in the Cage data structure.
 *
 * @param a Pointer to the atom to be deleted.
 */
void cageDeleteAtom(AtomCage_t *a) { lstDelete(neighborhood(a)); }

/**
 * @brief Deletes a Cage data structure.
 *
 * This function deletes the entire Cage data structure, including all atoms and the graph representation.
 *
 * @param s Pointer to the Cage data structure to be deleted.
 */
void cageDelete(Cage_t *s) {

  int i;

  if (s->atoms != NULL) {
    for (i = 0; i < size(s); i++)
      cageDeleteAtom(atom(s, i));
    free(s->atoms);
  }

  free(s);
}

/**************************************/
/* INITIALIZATION OF MOC **************/
/**************************************/

/**
 * @brief Imports a cage structure from a `.mol2` file and initializes its atoms and edges.
 *
 * This function reads a `.mol2` file for a given input name and modification number (e.g., _moc1),
 * processes the atoms and edges, and returns a pointer to a newly created `Cage_t` structure.
 * The cage structure is populated with atom data, their coordinates, and interconnections (edges).
 *
 * @param inputname The base name of the input file, which is used to construct the file path.
 * @param mocNum The modification number (e.g., "_moc1") that is appended to the input file name to form the complete
 * file path.
 *
 * @return A pointer to a `Cage_t` structure that contains the imported atoms and edges.
 *
 * @details
 * ### File Path Construction:
 * The function constructs a file path using the provided `inputname`, `mocNum`, and some predefined prefixes and
 * suffixes:
 * - Prefix: `./demos/`
 * - Suffix: `.mol2`
 * - The resulting file path follows the format: `./demos/{inputname}/{inputname}_moc{mocNum}.mol2`.
 *
 * ### File Reading:
 * - The function opens the file for reading. If the file is not found, an error message is printed and the program
 * exits.
 * - It skips a few header lines and reads the number of atoms and edges.
 * - Atom types and their coordinates are read, with specific atom flags being set based on the atom type:
 *   - `"S"` atoms are flagged as `CYCLE_F`.
 *   - `"H"` and `"U"` atoms are flagged as `HYDRO_PATTERN_F`.
 *   - `"P"` atoms are flagged as `LINKABLE_F`.
 * - The atoms are then initialized with their coordinates.
 * - After reading the atoms, the function reads the edges between atoms and establishes connections in the cage.
 *
 * ### Error Handling:
 * - If memory allocation fails during atom creation or reading of the file, the program will print an error message and
 * exit.
 * - If the file does not exist or cannot be read correctly, the program will exit with an error message.
 *
 * @see cageCreate
 * @see cageInitAtom
 * @see CAGE_addEdge
 * @see flag
 *
 */
Cage_t *cageImport(char *inputname, char *mocNum) {
  // Calculate the length of the fixed parts and the inputname
  const char *suffix = ".mol2";
  const char *access_moc = "_moc";
  int suffix_len = strlen(suffix);
  int access_moc_len = strlen(access_moc);
  int inputname_len = strlen(inputname);
  int moc_num_len = strlen(mocNum);
  char *name = getBasename(inputname);

  // Calculate the total length needed for the filepath
  int filepath_len;
  if (inputname[inputname_len - 1] == '/') {
    filepath_len = inputname_len + strlen(name) + access_moc_len + moc_num_len + suffix_len + 1;
  } else {
    filepath_len = inputname_len + 1 + strlen(name) + access_moc_len + moc_num_len + suffix_len + 1;
  }

  // Allocate memory for the filepath
  char *filepath = (char *)malloc(filepath_len);
  if (filepath == NULL) {
    printf("Error allocating memory");
    exit(EXIT_FAILURE);
  }

  // Construct the file path
  // Construct the file path
  if (inputname[inputname_len - 1] == '/') {
    sprintf(filepath, "%s%s%s%s%s", inputname, name, access_moc, mocNum, suffix);
  } else {
    sprintf(filepath, "%s/%s%s%s%s", inputname, name, access_moc, mocNum, suffix);
  }

  // Open the file
  FILE *filestream = fopen(filepath, "r");
  if (!filestream) {
    fprintf(stderr, "The file %s doesn't exist.\n", filepath);
    exit(EXIT_FAILURE);
  }

  free(filepath);
  free(name);
  int ret, nb_edges;
  Cage_t *import_cage = cageCreate();

  // Define variable use for reading
  char line[256];

  // Read to number of Atom (jump 2 first lines and read on third)
  if (fgets(line, sizeof(line), filestream) == NULL) {
    printf("Error reading file");
    fclose(filestream);
    exit(EXIT_FAILURE);
  }
  if (fgets(line, sizeof(line), filestream) == NULL) {
    printf("Error reading file");
    fclose(filestream);
    exit(EXIT_FAILURE);
  }

  ret = fscanf(filestream, "%d %d", &size(import_cage), &nb_edges);
  import_cage->atoms = (AtomCage_t *)malloc(size(import_cage) * sizeof(AtomCage_t));
  if (import_cage->atoms == NULL) {
    fprintf(stderr, "A problem occurred during the reallocation (structureCage.c:115).\n");
    exit(EXIT_FAILURE);
  }

  for (int i = 0; i < size(import_cage); i++) {
    cageInitAtom(atom(import_cage, i));
  }

  // Jump next 5 lines
  for (int i = 0; i < 5; i++) {
    if (fgets(line, sizeof(line), filestream) == NULL) {
      printf("Error reading file");
      fclose(filestream);
      exit(EXIT_FAILURE);
    }
  }
  int index;
  char atom_name[2];
  // int atom_flag;

  // init min/max
  import_cage->xMaxWOGap = -INFINITY;
  import_cage->xMinWOGap = INFINITY;
  import_cage->yMaxWOGap = -INFINITY;
  import_cage->yMinWOGap = INFINITY;
  import_cage->zMaxWOGap = -INFINITY;
  import_cage->zMinWOGap = INFINITY;

  for (int i = 0; i < size(import_cage); i++) {
    ret = fscanf(filestream, "%d %s", &index, atom_name);
    if (strcmp(atom_name, "S") == 0) {
      flag((atom(import_cage, index - 1))) = CYCLE_F;
    } else if ((strcmp(atom_name, "H") == 0) || (strcmp(atom_name, "U") == 0)) {
      flag((atom(import_cage, index - 1))) = HYDRO_PATTERN_F;
    } else if (strcmp(atom_name, "P") == 0) {
      flag((atom(import_cage, index - 1))) = LINKABLE_F;
    } else if (strcmp(atom_name, "C") == 0) {
      flag((atom(import_cage, index - 1))) = CARBON_F;
    } else {
      printf("unknown atom type: %s\n", atom_name);
      exit(EXIT_FAILURE);
    }
    ret = fscanf(filestream, "%lf %lf %lf %s", &atomX(atom(import_cage, index - 1)),
                 &atomY(atom(import_cage, index - 1)), &atomZ(atom(import_cage, index - 1)), atom_name);
    // Find min and max for each axis
    if (atomX(atom(import_cage, index - 1)) < import_cage->xMinWOGap)
      import_cage->xMinWOGap = atomX(atom(import_cage, index - 1));
    if (atomX(atom(import_cage, index - 1)) > import_cage->xMaxWOGap)
      import_cage->xMaxWOGap = atomX(atom(import_cage, index - 1));
    if (atomY(atom(import_cage, index - 1)) < import_cage->yMinWOGap)
      import_cage->yMinWOGap = atomY(atom(import_cage, index - 1));
    if (atomY(atom(import_cage, index - 1)) > import_cage->yMaxWOGap)
      import_cage->yMaxWOGap = atomY(atom(import_cage, index - 1));
    if (atomZ(atom(import_cage, index - 1)) < import_cage->zMinWOGap)
      import_cage->zMinWOGap = atomZ(atom(import_cage, index - 1));
    if (atomZ(atom(import_cage, index - 1)) > import_cage->zMaxWOGap)
      import_cage->zMaxWOGap = atomZ(atom(import_cage, index - 1));
    // Consume the newline character left in the input buffer
    int c;
    while ((c = fgetc(filestream)) != '\n' && c != EOF)
      ;
  }
  // Jump next 2 lines
  for (int i = 0; i < 2; i++) {
    if (fgets(line, sizeof(line), filestream) == NULL) {
      printf("Error reading file");
      fclose(filestream);
      exit(EXIT_FAILURE);
    }
  }

  int tmp, id1, id2;
  for (int i = 0; i < nb_edges; i++) {
    ret = fscanf(filestream, "%d %d %d", &tmp, &id1, &id2);
    cageAddEdge(import_cage, id1 - 1, id2 - 1);
    // Consume the newline character left in the input buffer
    int c;
    while ((c = fgetc(filestream)) != '\n' && c != EOF)
      ;
  }
  // CAGE_write(import_cage);

  fclose(filestream);

  if (ret < 0) {
    fprintf(stderr, "An error occured while reading %s.\n", inputname);
    exit(EXIT_FAILURE);
  }

  return import_cage;
}