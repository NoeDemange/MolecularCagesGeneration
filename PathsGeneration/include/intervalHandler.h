#ifndef _INTERVALHANDLER_H
#define _INTERVALHANDLER_H

#include "structure.h"

/**
 * @struct Interval
 * @brief Represents an angular range.
 */
typedef struct {
  double theta_min;
  double theta_max;
} Interval;

/**
 * @brief Sorts an array of intervals using merge sort.
 *
 * This function serves as a utility to sort an array of intervals in
 * ascending order based on theta_min.
 *
 * @param inter Array of intervals to be sorted.
 * @param size Number of elements in the array.
 */
void sortIntervals(Interval inter[], int size);

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
int mergeOverlappingIntervals(Interval intervals[], int n, Interval merged[]);

/**
 * @brief Prints an array of intervals.
 *
 * This function displays a list of intervals in the format [theta_min, theta_max].
 *
 * @param inter Array of intervals to print.
 * @param size Number of intervals in the array.
 */
void printIntervals(Interval *inter, int size);

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
                                     Interval carbon_intervals[], int max_carbon_intervals);

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
                        int max_output_size);

/**
 * @brief Compute the complement assuming `forbidden` is already sorted and merged.
 *
 * This variant expects `forbidden` to be non-overlapping, sorted by theta_min,
 * and with theta_min <= theta_max for every entry. It's cheaper than the
 * general `complementIntervals` because it skips merging/splitting logic.
 *
 * @param merged Array of sorted, merged forbidden intervals.
 * @param n_merged Number of intervals in `merged`.
 * @param valid Output array to store valid intervals (complement).
 * @param max_valid Capacity of the `valid` array.
 * @return Number of valid intervals written, 0 if none, -1 on error.
 */
int complementFromMerged(Interval *merged, int n_merged, Interval *valid, int max_valid);

/**
 * @brief Discretize valid intervals around a target angle.
 *
 * Generates up to `out_limit` angles starting from `theta_opt` and then
 * alternating by ±k·variation (k = 1, 2, ...), always selecting an angle that
 * lies within at least one valid interval and keeping at least `variation`
 * separation between selected angles. Search attempts are bounded by
 * `max_steps` (maximum k explored in each direction).
 *
 * Pre-conditions:
 * - `valid` must be sorted, non-overlapping, with 0 <= theta_min <= theta_max <= 2π
 *
 * @param valid Array of valid (merged) intervals.
 * @param n_valid Number of intervals in `valid`.
 * @param theta_opt Target angle to center the discretization.
 * @param variation Minimal angular separation between chosen thetas.
 * @param out_limit Maximum number of angles to generate (e.g., NUMBER_POSITION_AX1E3).
 * @param max_steps Maximum step index k to explore (e.g., MAX_POSITION_KEEP_360 / 2).
 * @param out_thetas Output array for angles (size >= out_limit).
 * @return Number of angles written to `out_thetas` (0..out_limit).
 */
int discretizeAroundThetaFromValid(Interval *valid, int n_valid, double theta_opt, double variation, int out_limit,
                                  int max_steps, double *out_thetas);

#endif