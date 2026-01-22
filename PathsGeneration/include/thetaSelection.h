#ifndef _THETASELECTION_H
#define _THETASELECTION_H

#include "intervalHandler.h"
#include "structure.h"

/**
 * @brief Generates discrete positions along a circular path.
 *
 * This function calculates and prints discrete positions along a circle
 * defined by the center (`circle_center`), and two
 * orthogonal vectors (`v1` and `v2`). The positions are computed at
 * 10-degree increments from 0 to 360 degrees.
 *
 * @param circle_center The center point of the circle.
 * @param circle_radius The radius of the circle.
 * @param v1 The first orthogonal vector defining the plane of the circle.
 * @param v2 The second orthogonal vector defining the plane of the circle.
 *
 * @note The function assumes `v1` and `v2` are already normalized and orthogonal.
 *
 * @warning The function directly prints the computed positions in a predefined format.
 *          If visualization or further processing is needed, consider modifying the function
 *          to return the computed positions instead of printing them.
 */
void cerclePositionDiscret(Point_t circle_center, double circle_radius, Point_t v1, Point_t v2);

/**
 * @brief Computes valid theta ranges for a given sphere.
 *
 * This function determines the valid angular range (theta) for a moving sphere
 * such that it does not intersect with a fixed sphere. The calculation is based
 * on the quadratic equation derived from the geometric constraints.
 *
 * @param circle_center Center of the moving sphere.
 * @param circle_squared_radius Squared radius of the moving sphere.
 * @param v1 First orthogonal vector defining the movement plane.
 * @param v2 Second orthogonal vector defining the movement plane.
 * @param sphere_fixed Center of the fixed sphere.
 * @param interval Pointer to store the resulting valid theta range.
 * @param dist_coli The distance of collision between the two centers of spheres.
 * @return 1 if a valid interval is found, 0 otherwise.
 */
int findValidThetaForSphere(Point_t circle_center, double circle_squared_radius, Point_t v1, Point_t v2,
                            Point_t sphere_fixed, Interval *interval, double dist_coli);

#endif