#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

// Define the objective function (Example: Sphere function)
double objective_function(double* x, int dim) {
    double sum = 0.0;
    for (int i = 0; i < dim; i++) {
        sum += x[i] * x[i];
    }
    return sum;
}

// Generate a random double between min and max
double rand_double(double min, double max) {
    return min + ((double)rand() / RAND_MAX) * (max - min);
}

// Initialize the population
void initialize_population(double** population, int pop_size, int dim, double bounds[2]) {
    for (int i = 0; i < pop_size; i++) {
        for (int j = 0; j < dim; j++) {
            population[i][j] = rand_double(bounds[0], bounds[1]);
        }
    }
}

// Update the position of grasshoppers
void update_position(double** population, double* best_solution, int pop_size, int dim, double c_max, double c_min, int max_iter, int current_iter) {
    double c = c_max - (c_max - c_min) * ((double)current_iter / max_iter);
    for (int i = 0; i < pop_size; i++) {
        double* s_i = (double*)calloc(dim, sizeof(double));
        for (int j = 0; j < pop_size; j++) {
            if (i != j) {
                double dist = 0.0;
                for (int k = 0; k < dim; k++) {
                    dist += (population[i][k] - population[j][k]) * (population[i][k] - population[j][k]);
                }
                dist = sqrt(dist);
                double* r_ij = (double*)malloc(dim * sizeof(double));
                for (int k = 0; k < dim; k++) {
                    r_ij[k] = (population[j][k] - population[i][k]) / (dist + 1e-10);
                }
                double x_ij = (2 * rand_double(0.0, 1.0) - 1) * exp(-dist / c);
                for (int k = 0; k < dim; k++) {
                    s_i[k] += x_ij * r_ij[k];
                }
                free(r_ij);
            }
        }
        for (int k = 0; k < dim; k++) {
            population[i][k] = best_solution[k] * exp(-c * sqrt((best_solution[k] - population[i][k]) * (best_solution[k] - population[i][k]))) + s_i[k];
        }
        free(s_i);
    }
}

// Grasshopper Optimization Algorithm
void grasshopper_optimization(double (*objective_function)(double*, int), double bounds[2], int pop_size, int dim, int max_iter, double c_max, double c_min) {
    double** population = (double**)malloc(pop_size * sizeof(double*));
    for (int i = 0; i < pop_size; i++) {
        population[i] = (double*)malloc(dim * sizeof(double));
    }

    initialize_population(population, pop_size, dim, bounds);

    double* fitness = (double*)malloc(pop_size * sizeof(double));
    for (int i = 0; i < pop_size; i++) {
        fitness[i] = objective_function(population[i], dim);
    }

    double* best_solution = (double*)malloc(dim * sizeof(double));
    double best_fitness = DBL_MAX;
    for (int i = 0; i < pop_size; i++) {
        if (fitness[i] < best_fitness) {
            best_fitness = fitness[i];
            for (int j = 0; j < dim; j++) {
                best_solution[j] = population[i][j];
            }
        }
    }

    for (int current_iter = 0; current_iter < max_iter; current_iter++) {
        update_position(population, best_solution, pop_size, dim, c_max, c_min, max_iter, current_iter);
        for (int i = 0; i < pop_size; i++) {
            fitness[i] = objective_function(population[i], dim);
        }
        for (int i = 0; i < pop_size; i++) {
            if (fitness[i] < best_fitness) {
                best_fitness = fitness[i];
                for (int j = 0; j < dim; j++) {
                    best_solution[j] = population[i][j];
                }
            }
        }
        printf("Iteration: %d, Best Fitness: %f\n", current_iter + 1, best_fitness);
    }

    printf("Best Solution: ");
    for (int i = 0; i < dim; i++) {
        printf("%f ", best_solution[i]);
    }
    printf("\nBest Fitness: %f\n", best_fitness);

    free(best_solution);
    free(fitness);
    for (int i = 0; i < pop_size; i++) {
        free(population[i]);
    }
    free(population);
}

int main() {
    double bounds[2] = {-10.0, 10.0};
    int pop_size = 30;
    int dim = 30;
    int max_iter = 500;
    double c_max = 1.0;
    double c_min = 0.00004;

    grasshopper_optimization(objective_function, bounds, pop_size, dim, max_iter, c_max, c_min);

    return 0;
}
