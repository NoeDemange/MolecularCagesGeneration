#include "structure.h"

/**
 * @file structurePT.c
 * @brief Point_t Data Structure Implementation and related operations
 *
 * This file contains the implementation of the Point_t data structure and related operations.
 * Point_t represents a 3D point with double coordinates (x, y, z).
 */

/**
 * @brief Initializes a new Point_t with equal scalar values.
 *
 * This function initializes a new Point_t with the same scalar value for all coordinates (x, y, z).
 *
 * @param scal The scalar value to set for all coordinates of the Point_t.
 * @return The initialized Point_t.
 */
Point_t ptInit(double scal) {
  Point_t new_point;

  new_point.x = scal;
  new_point.y = scal;
  new_point.z = scal;

  return new_point;
}

/**
 * @brief Adds two Point_t together and returns the result.
 *
 * This function adds two Point_t (A and B) together and returns the resulting Point_t.
 *
 * @param A The first Point_t operand.
 * @param B The second Point_t operand.
 * @return The result of the addition operation as a new Point_t.
 */
Point_t ptAdd(Point_t A, Point_t B) {
  Point_t new_point;

  new_point.x = A.x + B.x;
  new_point.y = A.y + B.y;
  new_point.z = A.z + B.z;

  return new_point;
}

/**
 * @brief Subtracts one Point_t from another and returns the result.
 *
 * This function subtracts Point_t B from Point_t A and returns the resulting Point_t.
 *
 * @param A The first Point_t operand.
 * @param B The second Point_t operand.
 * @return The result of the subtraction operation as a new Point_t.
 */
Point_t ptSub(Point_t A, Point_t B) {
  Point_t new_point;

  new_point.x = A.x - B.x;
  new_point.y = A.y - B.y;
  new_point.z = A.z - B.z;

  return new_point;
}

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
Point_t ptMul(Point_t A, double scal) {
  Point_t new_point;

  new_point.x = scal * A.x;
  new_point.y = scal * A.y;
  new_point.z = scal * A.z;

  return new_point;
}

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
Point_t ptDiv(Point_t A, double scal) {
  Point_t new_point;

  if (scal == 0)
    return ptInit(0);

  new_point.x = A.x / scal;
  new_point.y = A.y / scal;
  new_point.z = A.z / scal;

  return new_point;
}

/**
 * @brief Compares two 3D points for equality.
 *
 * This function compares the coordinates of two `Point_t` structures (`A` and `B`)
 * to check if they are equal. It compares the `x`, `y`, and `z` components of both points.
 *
 * @param A The first `Point_t` structure to compare.
 * @param B The second `Point_t` structure to compare.
 *
 * @return `1` (true) if the points are identical, `0` (false) if they are different.
 *
 * @note This comparison function performs an exact match, meaning that floating-point coordinates
 * are compared directly. For applications involving floating-point precision errors, consider using
 * a tolerance-based comparison.
 *
 * @see Point_t
 */
int ptCompare(Point_t A, Point_t B) { return A.x == B.x && A.y == B.y && A.z == B.z; }

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
Point_t ptMerge(Point_t A, Point_t B) {
  Point_t new_point;

  new_point.x = (A.x + B.x) / 2;
  new_point.y = (A.y + B.y) / 2;
  new_point.z = (A.z + B.z) / 2;

  return new_point;
}

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
int ptEqual(Point_t A, Point_t B) { // 1 if equal, if not 0
  if (A.x == B.x && A.y == B.y && A.z == B.z)
    return 1;
  return 0;
}