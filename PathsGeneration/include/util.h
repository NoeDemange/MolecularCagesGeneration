#ifndef __UTIL_H
#define __UTIL_H

#include "structure.h"
#include <time.h>

/**
 * @file util.h
 * @brief Utility Functions Header File
 *
 * This file contains function prototypes for various utility functions used in the project.
 */

/**
 * @brief Global variable to store the start time of the program.
 */
extern struct timespec start_time_ts;
/**
 * @brief Global variable to store the start clock of the program.
 */
extern clock_t start_clock;

double monotonic_now_ms(void);
double timespec_diff_ms(struct timespec start, struct timespec end);

/**
 * @brief Free a 2D array of double.
 *
 * @param arr Double pointer to the 2D Double.
 * @param cols Number of columns in the array.
 */
void free2DDouble(double **arr, int cols);

/**
 * @brief Convert radians to degrees.
 *
 * This function converts an angle from radians to degrees.
 *
 * @param a The angle in radians.
 * @return The angle in degrees.
 */
double radianToDegre(double);

/**
 * @brief Convert degrees to radians.
 *
 * This function converts an angle from degrees to radians.
 *
 * @param a The angle in degrees.
 * @return The angle in radians.
 */
double degreToRadian(double);

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
Point_t normalization(Point_t, double);

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
double angle(Point_t, Point_t, Point_t);

/**
 * @brief Compute the vector between two points.
 *
 * This function calculates the vector (direction) between two 3D points A and B.
 *
 * @param A The starting point.
 * @param B The ending point.
 * @return The vector from A to B.
 */
Point_t vector(Point_t, Point_t);

/**
 * @brief Computes the cross product of two 3D vectors.
 * @param u The first vector.
 * @param v The second vector.
 * @return The cross product vector.
 */
Point_t crossProduct(Point_t u, Point_t v);

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
Point_t planNormal(Point_t, Point_t, Point_t);

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
void planeVectors(Point_t nei, Point_t start, Point_t *v1, Point_t *v2);

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
Point_t rotation(Point_t, double, Point_t);

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
Point_t funAX1E1(Point_t, Point_t, double);

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
Point_t funAX2E1(Point_t, Point_t, Point_t, double);

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
Point_t funAX1E2(Point_t, Point_t, Point_t, double);

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
Point_t funAX3E1(Point_t, Point_t, Point_t, Point_t, double);

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
Point_t funAX2E2(Point_t, Point_t, Point_t, double);

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
Point_t funAX1E3(Point_t, Point_t, Point_t, double);

/**
 * @brief Extract the basename from the given input string.
 *
 * @param in The input string.
 * @return A dynamically allocated string representing the extracted basename.
 */
char *getBasename(const char *in);

/**
 * @brief Inline optimization macros.
 *
 * These macros are used to provide hints to the compiler about the expected likelihood of certain
 * conditions, allowing for potential optimization of branch predictions.
 */
#define LIKELY(x)   __builtin_expect(!!(x), 1) /**< Likely branch prediction hint. */
#define UNLIKELY(x) __builtin_expect(!!(x), 0) /**< Unlikely branch prediction hint. */

#endif