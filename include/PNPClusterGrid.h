//
// Created by suriya on 1/31/25.
//

#pragma once

#include <iostream>
#include <vector>
#include <cmath>

class PNPClusterGrid {
public:
    PNPClusterGrid(double x_start, double x_end, int N, double epsilon, int min_boundary_nodes = 9)
        : x_start(x_start), x_end(x_end), no_of_nodes(N), epsilon(epsilon), min_boundary_nodes(min_boundary_nodes) {}

    // Generate the cosine-spaced grid
    std::vector<double> generateGrid() const {

        //double Lx = x_end - x_start;
        //double h_ = Lx / (no_of_nodes - 1);

        //if (epsilon/h_ > 8.0) {
            //return computeLinearSpacing();
        //}

        //double optimal_amp = fitAmplitude(min_boundary_nodes);
        //return computeCosineSpacing(optimal_amp);
        double optimal_offset = fitOffset(min_boundary_nodes);
        return computeLogSpacing(optimal_offset);
    }

private:
    double x_start, x_end, epsilon;
    int no_of_nodes, min_boundary_nodes;

    // Function to iteratively adjust the amplitude to fit boundary nodes
    double fitAmplitude(int target_nodes, double initial_amp = 1.0, int max_iter = 1000) const {
        double amp = initial_amp;
        for (int iter = 0; iter < max_iter; iter++) {
            std::vector<double> grid = computeCosineSpacing(amp);

            // Count nodes in the epsilon region
            int num_left = 0, num_right = 0;
            for (double x : grid) {
                if (x <= x_start + epsilon) num_left++;
                if (x >= x_end - epsilon) num_right++;
            }

            // Check if the condition is met
            if (num_left >= target_nodes && num_right >= target_nodes) {
                return amp;
            }

            // Increase amplitude slightly to push more points to boundaries
            amp *= 1.1;
        }

        // Warning message if convergence not achieved
        std::cerr << "Warning: Max iterations reached in fitAmplitude(). Result may be suboptimal.\n";
        return amp;  // Return last value if exact match not found
    }

    // Function to compute the cosine spacing given an amplitude
    std::vector<double> computeCosineSpacing(double amplitude) const {
        std::vector<double> grid(no_of_nodes);
        for (int i = 0; i < no_of_nodes; i++) {
            double s = static_cast<double>(i) / (no_of_nodes - 1);
            grid[i] = x_start + (x_end - x_start) * amplitude * (1 - std::cos(M_PI * s) ) / 2;
        }
        return grid;
    }

    double fitOffset(int target_nodes, double initial_offset = 0.0, int max_iter = 1000) const {


        double offset = initial_offset;
        for (int iter = 0; iter < max_iter; iter++) {
            std::vector<double> grid = computeLogSpacing(offset);

            // Count nodes in the epsilon region
            int num_left = 0, num_right = 0;
            for (double x : grid) {
                if (x <= x_start + epsilon) num_left++;
                if (x >= x_end - epsilon) num_right++;
            }

            // Check if the condition is met
            if (num_left >= target_nodes && num_right >= target_nodes) {
                return offset;
            }

            // Increase amplitude slightly to push more points to boundaries
            offset += 0.01;
        }

        // Warning message if convergence not achieved
        std::cerr << "Warning: Max iterations reached in fitOffset(). Result may be suboptimal.\n";
        return offset;  // Return last value if exact match not found

    }

    std::vector<double> computeLinearSpacing() const {
        std::vector<double> grid(no_of_nodes);

        double h = (x_end - x_start) / (no_of_nodes - 1);
        for (int i = 0; i < no_of_nodes; i++) {
            grid[i] = x_start + i * h;
        }

        return grid;
    }


    std::vector<double> computeLogSpacing(double offset) const {

        if (offset == 0.0) {
            return computeLinearSpacing();
        }

        std::vector<double> grid(no_of_nodes);

        grid[0] = x_start;
        grid[no_of_nodes - 1] = x_end;

        int first_half = no_of_nodes / 2;
        int second_half = no_of_nodes - first_half;

        double half_length = (x_end - x_start) / 2;

        double nominal_grid_spacing = (x_end - x_start) / (no_of_nodes - 1);

        double smallest_grid_spacing = pow(10, log10(nominal_grid_spacing) - offset);
        double largest_grid_spacing = half_length;

        double log_grid_spacing_1 = (log10(largest_grid_spacing) - log10(smallest_grid_spacing))/(first_half - 2);
        double log_grid_spacing_2 = (log10(largest_grid_spacing) - log10(smallest_grid_spacing))/(second_half - 2);

        for (int i = 1; i < first_half; i++) {
            double exponent = log10(smallest_grid_spacing) + log_grid_spacing_1 * (i-1);
            double h = pow(10,exponent);
            grid[i] = x_start + h;
        }

        for (int i = 1; i < second_half; i++) {
            double exponent = log10(smallest_grid_spacing) + log_grid_spacing_2 * (i-1);
            double h = pow(10,exponent);
            grid[no_of_nodes - 1 - i] = x_end - h;
        }

        return grid;
    }
};
