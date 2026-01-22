#include "structure.h"

#include "structure.h"
#include <distance.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/**
 * @file structureMH.c
 * @brief MinHeap Data Structure Implementation
 *
 * This file contains the implementation of the MinHeap data structures
 * and related operations. The MinHeap represents a priority queue.
 */

/**
 * @brief Swap two nodes in the MinHeap.
 *
 * This function swaps two nodes in the MinHeap.
 *
 * @param a Pointer to the first node.
 * @param b Pointer to the second node.
 */
static void swapNodes(Node **a, Node **b) {
  Node *temp = *a;
  *a = *b;
  *b = temp;

  // Update heapIndex for the swapped nodes
  int temp_index = (*a)->heapIndex;
  (*a)->heapIndex = (*b)->heapIndex;
  (*b)->heapIndex = temp_index;
}

/**
 * @brief Get the index of the parent node.
 *
 * This function returns the index of the parent node for a given index.
 *
 * @param i The index of the child node.
 * @return The index of the parent node.
 */
static int parent(int i) { return (i - 1) / 2; }

/**
 * @brief Get the index of the left child node.
 *
 * This function returns the index of the left child node for a given index.
 *
 * @param i The index of the parent node.
 * @return The index of the left child node.
 */
static int leftChild(int i) { return 2 * i + 1; }

/**
 * @brief Get the index of the right child node.
 *
 * This function returns the index of the right child node for a given index.
 *
 * @param i The index of the parent node.
 * @return The index of the right child node.
 */
static int rightChild(int i) { return 2 * i + 2; }

/**
 * @brief Initialize a new MinHeap.
 *
 * This function initializes a new MinHeap with the specified maximum size.
 *
 * @param heap Pointer to the MinHeap to be initialized.
 * @param maxSize The maximum size of the MinHeap.
 */
void initMinHeap(MinHeap_t *heap, int maxSize) {
  // Free existing heap->nodes if already allocated
  if (heap->nodes) {
    free(heap->nodes);
    heap->nodes = NULL; // Set to NULL to avoid dangling pointers
  }
  heap->size = 0;
  heap->maxSize = maxSize;
  heap->nodes = (Node **)malloc(maxSize * sizeof(Node *));
}

/**
 * @brief Clear the MinHeap.
 *
 * This function clears the MinHeap by resetting its size to 0.
 *
 * @param heap Pointer to the MinHeap to be cleared.
 */
void clearMinHeap(MinHeap_t *heap) { heap->size = 0; }

/**
 * @brief Insert a node into the MinHeap.
 *
 * This function inserts a new node into the MinHeap while maintaining the MinHeap property.
 *
 * @param heap Pointer to the MinHeap.
 * @param node Pointer to the node to be inserted.
 */
void insertMinHeap(MinHeap_t *heap, Node *node) {
  if (heap->size == heap->maxSize) {
    printf("Heap overflow: Cannot insert node\n");
    return;
  }

  // Insert the new node at the end
  heap->nodes[heap->size] = node;
  node->heapIndex = heap->size; // Set the heapIndex
  int i = heap->size;
  heap->size++;

  // Fix the MinHeap property if violated
  while (i != 0 && heap->nodes[parent(i)]->fCost > heap->nodes[i]->fCost) {
    swapNodes(&heap->nodes[i], &heap->nodes[parent(i)]);
    i = parent(i);
  }
}

/**
 * @brief Min-Heapify a subtree rooted at index i.
 *
 * This function maintains the MinHeap property for a subtree rooted at index i.
 *
 * @param heap Pointer to the MinHeap.
 * @param i The index of the root of the subtree.
 */
static void minHeapify(MinHeap_t *heap, int i) {
  // Iterative version to avoid recursion overhead and stack overflow
  while (1) {
    int smallest = i;
    int left = leftChild(i);
    int right = rightChild(i);

    // Check left child
    if (left < heap->size && heap->nodes[left]->fCost < heap->nodes[smallest]->fCost) {
      smallest = left;
    }

    // Check right child
    if (right < heap->size && heap->nodes[right]->fCost < heap->nodes[smallest]->fCost) {
      smallest = right;
    }

    // If heap property is satisfied, we're done
    if (smallest == i) {
      break;
    }

    // Swap and continue with the child that was smallest
    swapNodes(&heap->nodes[i], &heap->nodes[smallest]);
    i = smallest;
  }
}

/**
 * @brief Extract the minimum node from the MinHeap.
 *
 * This function extracts the minimum node from the MinHeap and maintains the MinHeap property.
 *
 * @param heap Pointer to the MinHeap.
 * @return Pointer to the extracted minimum node.
 */
Node *extractMin(MinHeap_t *heap) {
  if (heap->size <= 0) {
    return NULL;
  }
  if (heap->size == 1) {
    heap->size--;
    heap->nodes[0]->heapIndex = -1; // Mark as removed
    return heap->nodes[0];
  }

  // Store the minimum value and remove it from the heap
  Node *root = heap->nodes[0];
  heap->nodes[0] = heap->nodes[heap->size - 1];
  heap->nodes[0]->heapIndex = 0; // Update heapIndex for the new root
  heap->size--;
  minHeapify(heap, 0);

  root->heapIndex = -1; // Mark as removed
  return root;
}

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
void decreaseKeyMinHeap(MinHeap_t *heap, Node *node, double newGCost) {
  int index = node->heapIndex;
  if (index < 0 || index >= heap->size) {
    printf("Error: Index out of bounds\n");
    return;
  }

  if (heap->nodes[index]->gCost <= newGCost) {
    printf("Error: New gCost is not smaller than the current gCost\n");
    return;
  }

  // Update the gCost and fCost of the node
  node->gCost = newGCost;
  node->fCost = node->gCost + node->hCost;

  // Fix the MinHeap property by moving the node up
  while (index != 0 && heap->nodes[parent(index)]->fCost > heap->nodes[index]->fCost) {
    swapNodes(&heap->nodes[index], &heap->nodes[parent(index)]);
    index = parent(index);
  }
}

/**
 * @brief Renew the MinHeap with new candidates.
 *
 * This function renews the MinHeap with new candidates by recalculating the hCost
 * for each node in the heap based on the minimum distance to a candidate that is not closed.
 * Optimized version with early termination and reduced heapification.
 *
 * @param heap Pointer to the MinHeap to be renewed.
 * @param candidates Array of candidate nodes.
 * @param nb_candidates Number of candidates in the array.
 */
void renewMinHeap(MinHeap_t *heap, Node **candidates, int nb_candidates) {
  // Early return if no candidates or empty heap
  if (nb_candidates == 0 || heap->size == 0) {
    return;
  }

  // Count non-closed candidates first for early optimization
  int active_candidates = 0;
  for (int j = 0; j < nb_candidates; j++) {
    if (candidates[j]->closed == 0) {
      active_candidates++;
    }
  }

  // If no active candidates, nothing to do
  if (active_candidates == 0) {
    return;
  }

  int heap_changed = 0;

  // Update hCost and fCost for each non-closed node in the heap
  for (int i = 0; i < heap->size; i++) {
    Node *node = heap->nodes[i];
    if (node->closed == 0) {
      double old_f_cost = node->fCost;
      double dmin = DBL_MAX;

      // Find minimum distance to active candidates
      for (int j = 0; j < nb_candidates; j++) {
        if (candidates[j]->closed == 0) {
          double dist_to_candidate = heuristic(node, candidates[j]);
          if (dist_to_candidate < dmin) {
            dmin = dist_to_candidate;
            // Early termination if we find a very close candidate
            if (dmin < 0.1) { // Small threshold for "very close"
              break;
            }
          }
        }
      }

      // Update costs only if they changed significantly
      if (fabs(node->hCost - dmin) > 1e-6) {
        node->hCost = dmin;
        node->fCost = node->gCost + node->hCost;

        // Mark that heap structure might need adjustment
        if (fabs(old_f_cost - node->fCost) > 1e-6) {
          heap_changed = 1;
        }
      }
    }
  }

  // Only rebuild heap if costs actually changed
  if (heap_changed) {
    // Use bottom-up heapification which is more efficient than full rebuild
    for (int i = (heap->size / 2) - 1; i >= 0; i--) {
      minHeapify(heap, i);
    }
  }

  // printf("Heap renewed with %d active candidates\n", active_candidates);
  // printMinHeap(heap);
}

void printMinHeap(MinHeap_t *heap) {
  printf("MinHeap: ");
  for (int i = 0; i < heap->size; i++) {
    printf("(%d, %d, %d) fCost %lf index %d", heap->nodes[i]->x, heap->nodes[i]->y, heap->nodes[i]->z,
           heap->nodes[i]->fCost, heap->nodes[i]->heapIndex);
    if (i % 4 == 3) {
      printf("\n");
    } else {
      printf(", ");
    }
  }
  printf("\n");
}

/**
 * @brief Free the memory allocated for the MinHeap.
 *
 * This function frees the memory allocated for the MinHeap and its nodes.
 *
 * @param heap Pointer to the MinHeap to be freed.
 */
void freeMinHeap(MinHeap_t *heap) {
  if (heap) {
    free(heap->nodes);
    free(heap);
  }
}

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
void renewMinHeapFast(MinHeap_t *heap, Node **candidates, int nb_candidates) {
  // Early return if no candidates or empty heap
  if (nb_candidates == 0 || heap->size == 0) {
    return;
  }

  // For small numbers of candidates, use the regular optimized version
  if (nb_candidates <= 5) {
    renewMinHeap(heap, candidates, nb_candidates);
    return;
  }

  // Pre-filter active candidates
  Node *active_candidates[nb_candidates];
  int active_count = 0;

  for (int j = 0; j < nb_candidates; j++) {
    if (candidates[j]->closed == 0) {
      active_candidates[active_count] = candidates[j];
      active_count++;
    }
  }

  if (active_count == 0) {
    return;
  }

  int heap_changed = 0;

  // Update hCost and fCost for each non-closed node in the heap
  for (int i = 0; i < heap->size; i++) {
    Node *node = heap->nodes[i];
    if (node->closed == 0) {
      double old_f_cost = node->fCost;
      double dmin = DBL_MAX;

      // Use the pre-filtered active candidates
      for (int j = 0; j < active_count; j++) {
        double dist_to_candidate = heuristic(node, active_candidates[j]);
        if (dist_to_candidate < dmin) {
          dmin = dist_to_candidate;
          // Early termination for very close candidates
          if (dmin < 0.05) {
            break;
          }
        }
      }

      // Update costs only if they changed significantly
      if (fabs(node->hCost - dmin) > 1e-6) {
        node->hCost = dmin;
        node->fCost = node->gCost + node->hCost;

        if (fabs(old_f_cost - node->fCost) > 1e-6) {
          heap_changed = 1;
        }
      }
    }
  }

  // Only rebuild heap if costs actually changed
  if (heap_changed) {
    // Use bottom-up heapification
    for (int i = (heap->size / 2) - 1; i >= 0; i--) {
      minHeapify(heap, i);
    }
  }
}