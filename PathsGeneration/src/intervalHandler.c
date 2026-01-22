#include "intervalHandler.h"
#include "constant.h"
#include <math.h>
#include <stdio.h>

/**
 * @brief Merges two sorted subarrays into a single sorted array.
 *
 * This function is part of the merge sort algorithm. It takes two
 * sorted halves of an array and merges them into a single sorted
 * sequence in-place.
 *
 * @param inter Array of intervals to be sorted.
 * @param left Starting index of the first subarray.
 * @param mid Middle index dividing the subarrays.
 * @param right Ending index of the second subarray.
 * @param temp Temporary array used for merging.
 */
void merge(Interval inter[], int left, int mid, int right, Interval temp[]) {
  int i = left;    // Index of first subarray
  int j = mid + 1; // Index of second subarray
  int k = left;    // Index for temporary array

  // Merge both halves while maintaining order
  while (i <= mid && j <= right) {
    if (inter[i].theta_min < inter[j].theta_min ||
        (inter[i].theta_min == inter[j].theta_min && inter[i].theta_max <= inter[j].theta_max)) {
      temp[k++] = inter[i++];
    } else {
      temp[k++] = inter[j++];
    }
  }

  // Add remaining elements from the first subarray
  while (i <= mid) {
    temp[k++] = inter[i++];
  }

  // Add remaining elements from the second subarray
  while (j <= right) {
    temp[k++] = inter[j++];
  }

  // Copy sorted elements back to the original array
  for (i = left; i <= right; i++) {
    inter[i] = temp[i];
  }
}

/**
 * @brief Recursive function implementing merge sort.
 *
 * This function recursively sorts an array of intervals using the
 * merge sort algorithm, ensuring efficient sorting in O(n log n) time.
 *
 * @param inter Array of intervals to be sorted.
 * @param left Starting index of the array.
 * @param right Ending index of the array.
 * @param temp Temporary array used for merging.
 */
void mergeSort(Interval inter[], int left, int right, Interval temp[]) {
  if (left < right) {
    int mid = left + (right - left) / 2;

    // Sort first half
    mergeSort(inter, left, mid, temp);
    // Sort second half
    mergeSort(inter, mid + 1, right, temp);
    // Merge both halves
    merge(inter, left, mid, right, temp);
  }
}

/**
 * @brief Sorts an array of intervals using merge sort.
 *
 * This function serves as a utility to sort an array of intervals in
 * ascending order based on theta_min.
 *
 * @param inter Array of intervals to be sorted.
 * @param size Number of elements in the array.
 */
void sortIntervals(Interval inter[], int size) {
  Interval temp[size]; // Temporary array to avoid dynamic allocation
  mergeSort(inter, 0, size - 1, temp);
}

/**
 * @brief Merges overlapping intervals in a sorted array.
 *
 * Given a sorted array of intervals, this function merges any overlapping
 * intervals and returns a new array containing only the non-overlapping
 * merged intervals.
 *
 * @param intervals Sorted array of intervals.
 * @param n Number of intervals in the input array.
 * @param merged Output array to store merged intervals.
 * @return The number of merged intervals.
 */
int mergeOverlappingIntervals(Interval intervals[], int n, Interval merged[]) {
  if (n == 0)
    return 0; // No intervals to process

  int index = 0;                // Index for merged intervals array
  merged[index] = intervals[0]; // Initialize with the first interval

  // Traverse all intervals
  for (int i = 1; i < n; i++) {
    // Check if current interval overlaps with the last merged interval
    // Use small epsilon to handle floating-point precision issues
    if (intervals[i].theta_min <= merged[index].theta_max + EPSILON_ANGLE) {
      // Extend the merged interval if necessary
      if (intervals[i].theta_max > merged[index].theta_max) {
        merged[index].theta_max = intervals[i].theta_max;
      }
    } else {
      // Add a new non-overlapping interval
      index++;
      merged[index] = intervals[i];
    }
  }

  return index + 1; // Total number of merged intervals
}

/**
 * @brief Prints an array of intervals.
 *
 * This function displays a list of intervals in the format [theta_min, theta_max].
 *
 * @param inter Array of intervals to print.
 * @param size Number of intervals in the array.
 */
void printIntervals(Interval *inter, int size) {
  for (int i = 0; i < size; i++) {
    printf("[%lf, %lf]\n", inter[i].theta_min, inter[i].theta_max);
  }
}

/**
 * @brief Normalizes an angle to the range [0, 2π].
 *
 * @param angle Input angle in radians.
 * @return Normalized angle in [0, 2π] range.
 */
double normalizeAngle(double angle) {
  while (angle < 0) {
    angle += 2 * M_PI;
  }
  while (angle >= 2 * M_PI) {
    angle -= 2 * M_PI;
  }
  return angle;
}

/**
 * @brief Adds an interval to an array, handling wraparound cases.
 *
 * @param intervals Target array.
 * @param count Current count of intervals.
 * @param max_count Maximum capacity.
 * @param theta_min Minimum angle.
 * @param theta_max Maximum angle.
 * @return 0 on success, -1 if capacity exceeded.
 */
int addInterval(Interval intervals[], int *count, int max_count, double theta_min, double theta_max) {
  if (*count >= max_count) {
    return -1; // Capacity exceeded
  }

  // Handle wraparound case
  if (theta_min > theta_max) {
    // Split into two intervals: [theta_min, 2π] and [0, theta_max]
    if (*count + 1 >= max_count) {
      return -1; // Not enough space for two intervals
    }
    intervals[*count].theta_min = theta_min;
    intervals[*count].theta_max = 2 * M_PI;
    (*count)++;
    intervals[*count].theta_min = 0;
    intervals[*count].theta_max = theta_max;
    (*count)++;
  } else {
    intervals[*count].theta_min = theta_min;
    intervals[*count].theta_max = theta_max;
    (*count)++;
  }

  return 0;
}

/**
 * @brief Converts hydrogen circle intervals to carbon circle intervals using angle shift.
 *
 * For tetrahedral carbon with 2 hydrogens, this function maps hydrogen circle intervals
 * back to the carbon circle by applying the ANGLE_SHIFT transformation. Each hydrogen
 * interval is mapped to two carbon intervals (+ and - ANGLE_SHIFT).
 *
 * @param hydrogen_intervals Array of hydrogen intervals.
 * @param nb_hydrogen_intervals Number of hydrogen intervals.
 * @param carbon_intervals Output array to store converted carbon intervals.
 * @param max_carbon_intervals Maximum capacity of carbon_intervals array.
 * @return Number of carbon intervals created, or -1 if capacity exceeded.
 */
int convertHydrogenIntervalsToCarbon(Interval hydrogen_intervals[], int nb_hydrogen_intervals,
                                     Interval carbon_intervals[], int max_carbon_intervals) {
  int count = 0;

  for (int i = 0; i < nb_hydrogen_intervals; i++) {
    // For each hydrogen interval, create two carbon intervals:
    // 1. hydrogen_theta + ANGLE_SHIFT corresponds to carbon_theta
    // 2. hydrogen_theta - ANGLE_SHIFT corresponds to carbon_theta

    // First mapping: carbon_theta = hydrogen_theta + ANGLE_SHIFT
    double carbon_min_1 = normalizeAngle(hydrogen_intervals[i].theta_min + ANGLE_SHIFT);
    double carbon_max_1 = normalizeAngle(hydrogen_intervals[i].theta_max + ANGLE_SHIFT);

    if (addInterval(carbon_intervals, &count, max_carbon_intervals, carbon_min_1, carbon_max_1) < 0) {
      return -1; // Capacity exceeded
    }

    // Second mapping: carbon_theta = hydrogen_theta - ANGLE_SHIFT
    double carbon_min_2 = normalizeAngle(hydrogen_intervals[i].theta_min - ANGLE_SHIFT);
    double carbon_max_2 = normalizeAngle(hydrogen_intervals[i].theta_max - ANGLE_SHIFT);

    if (addInterval(carbon_intervals, &count, max_carbon_intervals, carbon_min_2, carbon_max_2) < 0) {
      return -1; // Capacity exceeded
    }
  }

  return count;
}

/**
 * @brief Merges two sorted arrays of intervals into one sorted array.
 *
 * This function combines intervals from carbon circle and converted hydrogen intervals,
 * sorts them, and merges overlapping intervals.
 *
 * @param intervals1 First array of intervals.
 * @param size1 Number of intervals in first array.
 * @param intervals2 Second array of intervals.
 * @param size2 Number of intervals in second array.
 * @param merged_output Output array for merged intervals.
 * @param max_output_size Maximum capacity of output array.
 * @return Number of merged intervals, or -1 if capacity exceeded.
 */
int mergeIntervalArrays(Interval intervals1[], int size1, Interval intervals2[], int size2, Interval merged_output[],
                        int max_output_size) {

  // Check if we have enough space for all intervals
  if (size1 + size2 > max_output_size) {
    return -1; // Not enough space
  }

  // Combine both arrays into a temporary array
  Interval combined[size1 + size2];
  int combined_count = 0;

  // Copy intervals from first array
  for (int i = 0; i < size1; i++) {
    combined[combined_count++] = intervals1[i];
  }

  // Copy intervals from second array
  for (int i = 0; i < size2; i++) {
    combined[combined_count++] = intervals2[i];
  }

  // Sort the combined array
  sortIntervals(combined, combined_count);

  // Merge overlapping intervals
  int merged_count = mergeOverlappingIntervals(combined, combined_count, merged_output);

  return merged_count;
}

/**
 * @brief Compute complement assuming `merged` is sorted and non-overlapping.
 *
 * This function is intended for the case where the caller already has a
 * sorted/merged list of forbidden intervals (e.g. output of
 * `mergeOverlappingIntervals`). It simply returns the gaps between those
 * intervals on [0,2π).
 */
int complementFromMerged(Interval *merged, int n_merged, Interval *valid, int max_valid) {
  if (max_valid <= 0) {
    return -1;
  }

  // Empty forbidden list -> whole circle valid
  if (n_merged <= 0) {
    if (max_valid < 1)
      return -1;
    valid[0].theta_min = 0.0;
    valid[0].theta_max = 2 * M_PI;
    return 1;
  }

  // If the single merged interval covers the whole circle -> no valid intervals 
  // Normally case already handled by caller, but just in case
  if (n_merged == 1 && merged[0].theta_min <= 0.0 + EPSILON_ANGLE &&
      merged[0].theta_max >= 2 * M_PI - EPSILON_ANGLE) {
    return 0;
  }

  int valid_count = 0;
  double cursor = 0.0;

  for (int i = 0; i < n_merged; i++) {
    double a = merged[i].theta_min;
    double b = merged[i].theta_max;

    if (a > cursor + EPSILON_ANGLE) {
      if (valid_count >= max_valid) return -1;
      valid[valid_count].theta_min = cursor;
      valid[valid_count].theta_max = a;
      valid_count++;
    }

    if (b > cursor) cursor = b;
  }

  // Tail gap
  if (cursor < 2 * M_PI - EPSILON_ANGLE) {
    if (valid_count >= max_valid) return -1;
    valid[valid_count].theta_min = cursor;
    valid[valid_count].theta_max = 2 * M_PI;
    valid_count++;
  }

  return valid_count;
}

static inline double angular_distance(double a, double b) {
  double d = fabs(a - b);
  if (d > M_PI) d = (2 * M_PI) - d;
  return d;
}

static int is_in_any_valid(double theta, Interval *valid, int n_valid) {
  for (int i = 0; i < n_valid; ++i) {
    if (theta + EPSILON_ANGLE >= valid[i].theta_min && theta <= valid[i].theta_max + EPSILON_ANGLE) {
      return 1;
    }
  }
  return 0;
}

// Given a target angle, find the closest angle that lies inside the valid intervals
// by clamping target to each [a,b] and keeping the closest one.
static double closest_point_in_valid(double target, Interval *valid, int n_valid) {
  double best_theta = target;
  double best_dist = 1e300;
  for (int i = 0; i < n_valid; ++i) {
    double a = valid[i].theta_min;
    double b = valid[i].theta_max;
    double cand = target;
    if (cand < a) cand = a;
    if (cand > b) cand = b;
    double d = angular_distance(cand, target);
    if (d < best_dist) {
      best_dist = d;
      best_theta = cand;
    }
  }
  return best_theta;
}

static double next_valid_in_direction(double target, int sign, Interval *valid, int n_valid) {
  target = normalizeAngle(target);

  if (is_in_any_valid(target, valid, n_valid)) {
    return target;
  }

  if (sign >= 0) {
    for (int i = 0; i < n_valid; ++i) {
      if (valid[i].theta_min > target + EPSILON_ANGLE) {
        return valid[i].theta_min;
      }
    }
    return valid[0].theta_min; // wrap around
  } else {
    for (int i = n_valid - 1; i >= 0; --i) {
      if (valid[i].theta_max < target - EPSILON_ANGLE) {
        return valid[i].theta_max;
      }
    }
    return valid[n_valid - 1].theta_max; // wrap around backwards
  }
}

int discretizeAroundThetaFromValid(Interval *valid, int n_valid, double theta_opt, double variation, int out_limit,
                                  int max_steps, double *out_thetas) {
  if (out_limit <= 0) return 0;
  if (n_valid <= 0) return 0;

  theta_opt = normalizeAngle(theta_opt);

  int count = 0;

  // First angle: the closest point to theta_opt that is inside valid
  double first = is_in_any_valid(theta_opt, valid, n_valid)
                     ? theta_opt
                     : closest_point_in_valid(theta_opt, valid, n_valid);
  out_thetas[count++] = first;

  // Explore ±k * variation
  for (int k = 1; k <= max_steps && count < out_limit; ++k) {
    double target_plus = normalizeAngle(theta_opt + (+1) * k * variation);
    double target_minus = normalizeAngle(theta_opt + (-1) * k * variation);

    double cand_plus = next_valid_in_direction(target_plus, +1, valid, n_valid);
    double cand_minus = next_valid_in_direction(target_minus, -1, valid, n_valid);
    // double cand_plus = is_in_any_valid(target_plus, valid, n_valid) ? target_plus : closest_point_in_valid(target_plus, valid, n_valid);
    // double cand_minus = is_in_any_valid(target_minus, valid, n_valid) ? target_minus : closest_point_in_valid(target_minus, valid, n_valid);

    double d_plus = angular_distance(cand_plus, theta_opt);
    double d_minus = angular_distance(cand_minus, theta_opt);

    // Determine order: pick closer to theta_opt first
    int first_is_plus = (d_plus <= d_minus);
    for (int pass = 0; pass < 2 && count < out_limit; ++pass) {
      double candidate = first_is_plus ? (pass == 0 ? cand_plus : cand_minus)
                                       : (pass == 0 ? cand_minus : cand_plus);

      // Enforce minimal separation to all previously chosen thetas
      int ok = 1;
      for (int i = 0; i < count; ++i) {
        if (angular_distance(candidate, out_thetas[i]) < variation - EPSILON_ANGLE) {
          ok = 0;
          break;
        }
      }
      if (!ok) continue;

      // Avoid duplicates due to clamping to same boundary
      int dup = 0;
      for (int i = 0; i < count; ++i) {
        if (angular_distance(candidate, out_thetas[i]) < EPSILON_ANGLE) {
          dup = 1;
          break;
        }
      }
      if (dup) continue;

      out_thetas[count++] = candidate;
    }
  }

  return count;
}