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
        : x_start(x_start), x_end(x_end), N(N), epsilon(epsilon), min_boundary_nodes(min_boundary_nodes) {}

    // Generate the cosine-spaced grid
    std::vector<double> generateGrid() const {
        double optimal_amp = fitAmplitude(min_boundary_nodes);
        return computeCosineSpacing(optimal_amp);
    }

private:
    double x_start, x_end, epsilon;
    int N, min_boundary_nodes;

    // Function to iteratively adjust the amplitude to fit boundary nodes
    double fitAmplitude(int target_nodes, double initial_amp = 1.0, int max_iter = 100) const {
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
            amp *= 1.05;
        }

        // Warning message if convergence not achieved
        std::cerr << "Warning: Max iterations reached in fitAmplitude(). Result may be suboptimal.\n";
        return amp;  // Return last value if exact match not found
    }

    // Function to compute the cosine spacing given an amplitude
    std::vector<double> computeCosineSpacing(double amplitude) const {
        std::vector<double> grid(N);
        for (int i = 0; i < N; i++) {
            double s = static_cast<double>(i) / (N - 1);
            grid[i] = x_start + (x_end - x_start) * (1 - std::cos(M_PI * s * amplitude)) / 2;
        }
        return grid;
    }
};
