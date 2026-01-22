#include "thetaSelection.h"
#include "constant.h"
#include "util.h"

#include <math.h>

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
void cerclePositionDiscret(Point_t circle_center, double circle_radius, Point_t v1, Point_t v2) {
  Point_t circle;
  double rad = 0.0;
  for (int theta = 0; theta < 360; theta = theta + 10) {
    rad = theta * M_PI / 180;
    circle.x = circle_center.x + circle_radius * cos(rad) * v1.x + circle_radius * sin(rad) * v2.x;
    circle.y = circle_center.y + circle_radius * cos(rad) * v1.y + circle_radius * sin(rad) * v2.y;
    circle.z = circle_center.z + circle_radius * cos(rad) * v1.z + circle_radius * sin(rad) * v2.z;
    printf("pseudoatom theta%d, pos=[%lf,%lf,%lf]\n", theta, circle.x, circle.y, circle.z);
  }
}

/**
 * @brief Solves a quadratic inequality.
 *
 * This function finds the real roots of the quadratic equation:
 *
 * \f$ a \cdot x^2 + b \cdot x + c = 0 \f$
 *
 * If real solutions exist, they are stored in `tan_theta1` and `tan_theta2`.
 *
 * @param a Coefficient of the quadratic term.
 * @param b Coefficient of the linear term.
 * @param c Constant term.
 * @param tan_theta1 Pointer to store the first root.
 * @param tan_theta2 Pointer to store the second root.
 * @return 1 if real solutions exist, 0 otherwise.
 */
int solveQuadratic(double a, double b, double c, double *tan_theta1, double *tan_theta2) {
  double discriminant = (b * b) - (4 * a * c);

  if (discriminant < 0) {
    // No real solutions
    return 0;
  }

  // Compute the roots of the quadratic equation
  *tan_theta1 = (-b - sqrt(discriminant)) / (2 * a);
  *tan_theta2 = (-b + sqrt(discriminant)) / (2 * a);
  return 1;
}

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
                            Point_t sphere_fixed, Interval *interval, double dist_coli) {
  Point_t di = vector(sphere_fixed, circle_center);
  // printf("pseudoatom di, pos=[%lf,%lf,%lf]\n", di.x, di.y, di.z);

  double alpha = di.x * v1.x + di.y * v1.y + di.z * v1.z;
  double beta = di.x * v2.x + di.y * v2.y + di.z * v2.z;
  double gamma = (dist_coli * dist_coli) - (circle_squared_radius) - (di.x * di.x + di.y * di.y + di.z * di.z);

  // Quadratic equation coefficients
  double a = (4 * circle_squared_radius * beta * beta) - (gamma * gamma);
  double b = 8 * circle_squared_radius * alpha * beta;
  double c = (4 * circle_squared_radius * alpha * alpha) - (gamma * gamma);

  double tan_theta1, tan_theta2;
  if (!solveQuadratic(a, b, c, &tan_theta1, &tan_theta2)) {
    return 0; // No valid intervals
  }

  // Determine theta_min and theta_max based on conditions
  if (a > 0) {
    if (beta > 0) {
      interval->theta_min = atan(tan_theta2) + M_PI; // Shift by 180°
      interval->theta_max = atan(tan_theta1);
    } else {
      interval->theta_min = atan(tan_theta2);
      interval->theta_max = atan(tan_theta1) + M_PI; // Shift by 180°
    }
  } else {
    if (alpha > 0) {
      // Shift by 180° and reverse interval direction
      interval->theta_min = atan(tan_theta2) + M_PI;
      interval->theta_max = atan(tan_theta1) + M_PI;
    } else {
      // Reverse interval direction
      interval->theta_min = atan(tan_theta2);
      interval->theta_max = atan(tan_theta1);
    }
  }

  // Ensure theta is within [0, 2π]
  if (interval->theta_min < 0) {
    interval->theta_min += 2 * M_PI;
  }
  if (interval->theta_max < 0) {
    interval->theta_max += 2 * M_PI;
  }

  return 1;
}
