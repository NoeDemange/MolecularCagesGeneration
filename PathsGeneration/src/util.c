#include "util.h"
#include "distance.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

double monotonic_now_ms(void) {
#ifdef CLOCK_MONOTONIC
  struct timespec now;
  clock_gettime(CLOCK_MONOTONIC, &now);
  return (double)now.tv_sec * 1000.0 + (double)now.tv_nsec / 1e6;
#else
  return (double)clock() * 1000.0 / CLOCKS_PER_SEC;
#endif
}

double timespec_diff_ms(struct timespec start, struct timespec end) {
  time_t sec = end.tv_sec - start.tv_sec;
  long nsec = end.tv_nsec - start.tv_nsec;
  if (nsec < 0) {
    sec -= 1;
    nsec += 1000000000L;
  }
  return (double)sec * 1000.0 + (double)nsec / 1e6;
}
/**
 * @file util.c
 * @brief Utility Functions File
 */

/**
 * @brief Free a 2D array of double.
 *
 * @param arr Double pointer to the 2D Double.
 * @param rows Number of rows in the array.
 */
void free2DDouble(double **arr, int rows) {
  for (int i = 0; i < rows; i++) {
    free(arr[i]);
  }
  free(arr);
};

/**
 * @brief Convert radians to degrees.
 *
 * This function converts an angle from radians to degrees.
 *
 * @param a The angle in radians.
 * @return The angle in degrees.
 */
double radianToDegre(double a) { return a * 180 / M_PI; }

/**
 * @brief Convert degrees to radians.
 *
 * This function converts an angle from degrees to radians.
 *
 * @param a The angle in degrees.
 * @return The angle in radians.
 */
double degreToRadian(double a) { return a * M_PI / 180; }

/**
 * @brief Normalize a vector to a specified length.
 *
 * This function normalizes a given vector to the specified length. The resulting vector will have
 * the same direction as the original vector but with a magnitude equal to the specified length.
 *
 * @param normal The vector to be normalized (Point_t structure).
 * @param length The desired length of the normalized vector.
 * @return The normalized vector (Point_t structure).
 */
Point_t normalization(Point_t normal, double length) {
  Point_t a;
  double z;

  z = sqrt((length * length) / (normal.x * normal.x + normal.y * normal.y + normal.z * normal.z));
  a.x = z * normal.x;
  a.y = z * normal.y;
  a.z = z * normal.z;
  return a;
}

/**
 * @brief Compute the angle at A of a ABC triangle.
 *
 * * This function calculates the angle at vertex A of a triangle ABC using the law of cosines.
 *
 * @param A The first vertex of the triangle.
 * @param B The second vertex of the triangle.
 * @param C The third vertex of the triangle.
 * @return The angle at vertex A in degrees.
 */
double angle(Point_t A, Point_t B, Point_t C) {
  double ab = dist(A, B), ac = dist(A, C), bc = dist(B, C);

  return acos((ac * ac + ab * ab - bc * bc) / (2 * ac * ab)) * 180 / M_PI;
}

/**
 * @brief Compute the vector between two points.
 *
 * This function calculates the vector (direction) between two 3D points A and B.
 *
 * @param A The starting point.
 * @param B The ending point.
 * @return The vector from A to B.
 */
Point_t vector(Point_t A, Point_t B) {
  Point_t new_point;

  new_point.x = B.x - A.x;
  new_point.y = B.y - A.y;
  new_point.z = B.z - A.z;

  return new_point;
}

/**
 * @brief Computes the cross product of two 3D vectors.
 * @param u The first vector.
 * @param v The second vector.
 * @return The cross product vector.
 */
Point_t crossProduct(Point_t u, Point_t v) {
  Point_t w;
  w.x = u.y * v.z - u.z * v.y;
  w.y = u.z * v.x - u.x * v.z;
  w.z = u.x * v.y - u.y * v.x;
  return w;
}

/**
 * @brief Compute the normal vector of a plane defined by three points.
 *
 * This function calculates the normal vector of a plane defined by three 3D points A, B, and C.
 *
 * @param A The first point of the plane.
 * @param B The second point of the plane.
 * @param C The third point of the plane.
 * @return The normal vector of the plane.
 */
Point_t planNormal(Point_t A, Point_t B, Point_t C) {
  Point_t normal;

  normal.x = (B.y - A.y) * (C.z - A.z) - (B.z - A.z) * (C.y - A.y);
  normal.y = (B.z - A.z) * (C.x - A.x) - (B.x - A.x) * (C.z - A.z);
  normal.z = (B.x - A.x) * (C.y - A.y) - (B.y - A.y) * (C.x - A.x);

  return normalization(normal, 1);
}

/**
 * @brief Computes two orthogonal vectors defining a plane.
 *
 * This function calculates two vectors (`v1` and `v2`) that define a plane
 * orthogonal to the direction vector from `nei` to `start`.
 *
 * @param nei The reference point (neighbor).
 * @param start The starting point.
 * @param v1 Pointer to the first orthogonal vector of the plane.
 * @param v2 Pointer to the second orthogonal vector of the plane.
 *
 * @note The function assumes that `nei` and `start` are not the same point to avoid
 *       division by zero in the normal computation.
 *
 * @warning If `dir.x` is very close to zero, the normal vector computation may
 *          produce unstable results.
 */
void planeVectors(Point_t nei, Point_t start, Point_t *v1, Point_t *v2) {
  Point_t dir = normalization(vector(nei, start), 1);

  // Choix du vecteur de référence pour minimiser l'erreur numérique
  Point_t ref;
  if (fabsl(dir.x) <= fabsl(dir.y) && fabsl(dir.x) <= fabsl(dir.z)) {
    ref = (Point_t){0, -dir.z, dir.y}; // Favorise X si possible
  } else if (fabsl(dir.y) <= fabsl(dir.z)) {
    ref = (Point_t){-dir.z, 0, dir.x}; // Favorise Y
  } else {
    ref = (Point_t){dir.y, -dir.x, 0}; // Favorise Z
  }

  // Normalisation du premier vecteur orthogonal
  *v1 = normalization(ref, 1);

  // Calcul du second vecteur orthogonal via le produit vectoriel
  *v2 = crossProduct(dir, *v1);
}

/**
 * @brief Perform rotation around a point.
 *
 * This function rotates a point A around a given vector vec by a specified angle alpha.
 *
 * @param vec The vector used for rotation.
 * @param alpha The angle of rotation in degrees.
 * @param A The point to be rotated.
 * @return The rotated point.
 */
Point_t rotation(Point_t vec, double alpha, Point_t A) {
  Point_t rot;
  alpha = degreToRadian(alpha);
  vec = normalization(vec, 1);

  rot.x = (vec.x * vec.x + (1 - vec.x * vec.x) * cos(alpha)) * A.x +
          (vec.x * vec.y * (1 - cos(alpha)) - vec.z * sin(alpha)) * A.y +
          (vec.x * vec.z * (1 - cos(alpha)) + vec.y * sin(alpha)) * A.z;
  rot.y = (vec.y * vec.y + (1 - vec.y * vec.y) * cos(alpha)) * A.y +
          (vec.x * vec.y * (1 - cos(alpha)) + vec.z * sin(alpha)) * A.x +
          (vec.y * vec.z * (1 - cos(alpha)) - vec.x * sin(alpha)) * A.z;
  rot.z = (vec.z * vec.z + (1 - vec.z * vec.z) * cos(alpha)) * A.z +
          (vec.x * vec.z * (1 - cos(alpha)) - vec.y * sin(alpha)) * A.x +
          (vec.y * vec.z * (1 - cos(alpha)) + vec.x * sin(alpha)) * A.y;

  return rot;
}

/**
 * @brief Compute the point obtained by moving a distance length from point a in the direction of vector x1.
 *
 * This function computes the point obtained by moving a distance length from point a in the direction
 * of the normalized vector x1.
 *
 * @param a The starting point.
 * @param x1 The vector direction.
 * @param length The distance to move.
 * @return The new point after the move.
 */
Point_t funAX1E1(Point_t a, Point_t x1, double length) { return ptAdd(a, normalization(vector(x1, a), length)); }

/**
 * @brief Compute the point obtained by moving a distance length from point a in the direction of the sum of vectors x1
 * and x2.
 *
 * This function computes the point obtained by moving a distance length from point a in the direction
 * of the normalized sum of vectors x1 and x2.
 *
 * @param a The starting point.
 * @param x1 The first vector direction.
 * @param x2 The second vector direction.
 * @param length The distance to move.
 * @return The new point after the move.
 */
Point_t funAX2E1(Point_t a, Point_t x1, Point_t x2, double length) {

  Point_t v1, v2;

  v1 = normalization(vector(x1, a), 1);
  v2 = normalization(vector(x2, a), 1);

  return ptAdd(a, normalization(ptAdd(v1, v2), length));
}

/**
 * @brief Compute the point obtained by moving a distance length from point a in the direction of the rotated vector x1.
 *
 * This function computes the point obtained by moving a distance length from point a in the direction
 * of the normalized vector x1 rotated around the specified normal vector.
 *
 * @param a The starting point.
 * @param x1 The vector direction to rotate.
 * @param normal The normal vector around which to perform the rotation.
 * @param length The distance to move.
 * @return The new point after the move.
 */
Point_t funAX1E2(Point_t a, Point_t x1, Point_t normal, double length) {

  Point_t v1;

  v1 = normalization(vector(a, x1), 1);

  return ptAdd(a, normalization(rotation(normal, 120, v1), length));
}

/**
 * @brief Compute the point obtained by moving a distance length from point a in the direction of the sum of vectors x1,
 * x2, and x3.
 *
 * This function computes the point obtained by moving a distance length from point a in the direction
 * of the normalized sum of vectors x1, x2, and x3.
 *
 * @param a The starting point.
 * @param x1 The first vector direction.
 * @param x2 The second vector direction.
 * @param x3 The third vector direction.
 * @param length The distance to move.
 * @return The new point after the move.
 */
Point_t funAX3E1(Point_t a, Point_t x1, Point_t x2, Point_t x3, double length) {

  Point_t v1, v2, v3;

  v1 = normalization(vector(x1, a), 1);
  v2 = normalization(vector(x2, a), 1);
  v3 = normalization(vector(x3, a), 1);

  return ptAdd(a, normalization(ptAdd(v1, ptAdd(v2, v3)), length));
}

/**
 * @brief Compute the point obtained by moving a distance length from point a in the direction of the sum of vectors x1
 * and x2 and performing a rotation.
 *
 * This function computes the point obtained by moving a distance length from point a in the direction
 * of the normalized sum of vectors x1 and x2, and then performing a rotation around a calculated
 * normal vector.
 *
 * @param a The starting point.
 * @param x1 The first vector direction.
 * @param x2 The second vector direction.
 * @param length The distance to move.
 * @return The new point after the move.
 */
Point_t funAX2E2(Point_t a, Point_t x1, Point_t x2, double length) {

  Point_t v1, v2, other, normal, zero = {0, 0, 0};
  double angle = 180 - (109.47 / 2);

  v1 = normalization(vector(a, x1), 1);
  v2 = normalization(vector(a, x2), 1);

  other = normalization(ptAdd(v1, v2), 1);
  normal = normalization(planNormal(zero, planNormal(a, x1, x2), other), 1);

  return ptAdd(a, normalization(rotation(normal, angle, other), length));
}

/**
 * @brief Compute the point obtained by moving a distance length from point a in the direction of the rotated vector x1.
 *
 * This function computes the point obtained by moving a distance length from point a in the direction
 * of the normalized vector x1 rotated around the specified normal vector.
 *
 * @param a The starting point.
 * @param x1 The vector direction to rotate.
 * @param normal The normal vector around which to perform the rotation.
 * @param length The distance to move.
 * @return The new point after the move.
 */
Point_t funAX1E3(Point_t a, Point_t x1, Point_t normal, double length) {

  Point_t v1;

  v1 = normalization(vector(a, x1), 1);

  return ptAdd(a, normalization(rotation(normal, 109.47, v1), length));
}

/// STRING

/**
 * @brief Extract the basename from the given input string.
 *
 * @param in The input string.
 * @return A dynamically allocated string representing the extracted basename.
 */
char *getBasename(const char *in) {
  if (in == NULL || *in == '\0') {
    return NULL; // Handle empty string case
  }

  // Copy input to a temporary buffer to modify it
  char *temp = strdup(in);
  if (!temp) {
    return NULL; // Memory allocation failure
  }

  // Remove trailing slash if it exists
  size_t len = strlen(temp);
  if (len > 1 && temp[len - 1] == '/') {
    temp[len - 1] = '\0';
  }

  // Find the last occurrence of '/'
  char *start = strrchr(temp, '/');
  if (start) {
    start++; // Move past the '/'
  } else {
    start = temp; // No '/' found, use the whole string
  }

  // Allocate memory for the result
  char *result = strdup(start);
  free(temp); // Free temp buffer

  return result;
}