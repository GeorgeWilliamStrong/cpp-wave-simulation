// Filename: main.cpp
// Author: George Strong
// Created: June 3rd, 2018
// Description: Minimal 2D wave equation solver for variable velocity model using the 
// finite difference method and dynamic memory allocation
// Copyright © 2018 George Strong. All rights reserved.

#include <iostream>
#include <cmath>
#include <algorithm>
#include <fstream>


// Variable declarations
double lx, // Size of the model in x dimension (in m)
       ly, // Size of the model in y dimension (in m)
       dx, // Spatial grid spacing in x dimension (in m)
       dy, // Spatial grid spacing in y dimension (in m)
       t,  // Total time to run simulation (in s)
       dt, // Temporal grid spacing (in s)
       c,  // Courant number
       v1, // Velocity of upper layer (in m/s)
       v2, // Velocity of upper-mid layer (in m/s)
       v3, // Velocity of lower-mid layer (in m/s)
       v4, // Velocity of lower layer (in m/s)
       a,  // Source amplitude
       s,  // Source spread (in m)
       ox, // Source location in x dimension (given in cells)
       oy; // Source location in y dimension (given in cells)

int nt, // Number of time samples
    nx, // Number of grid cells in x dimension
    ny; // Number of grid cells in y dimension

std::string parameter_loading_method;


/**
 * @brief Prompts the user to manually enter parameters for a simulation model.
 *
 * This function interacts with the user through the console, prompting them to enter various parameters
 * required for the simulation. The entered parameters include dimensions, grid spacing, time,
 * Courant number, layer velocities and source properties.
 */
void enter_parameters(void)
{
    // Prompt user for model dimensions
    std::cout << "Enter model length in x dimension: ";
    std::cin >> lx;
    std::cout << "Enter model length in y dimension: ";
    std::cin >> ly;

    // Prompt user for grid spacing and assume isotropy
    std::cout << "Enter model grid spacing: ";
    std::cin >> dx;
    dy = dx; // Assume the model is isotropic

    // Prompt user for time and Courant number
    std::cout << "Enter time: ";
    std::cin >> t;
    std::cout << "Enter Courant number: ";
    std::cin >> c;

    // Prompt user for layer velocities
    std::cout << "Enter velocity of upper layer: ";
    std::cin >> v1;
    std::cout << "Enter velocity of upper-mid layer: ";
    std::cin >> v2;
    std::cout << "Enter velocity of lower-mid layer: ";
    std::cin >> v3;
    std::cout << "Enter velocity of lower layer: ";
    std::cin >> v4;

    // Square the velocities as they appear in the wave equation
    v1 = pow(v1, 2);
    v2 = pow(v2, 2);
    v3 = pow(v3, 2);
    v4 = pow(v4, 2);

    // Prompt user for source properties
    std::cout << "Enter source amplitude: ";
    std::cin >> a;
    std::cout << "Enter source spread: ";
    std::cin >> s;
    std::cout << "Enter source location in x dimension: ";
    std::cin >> ox;
    std::cout << "Enter source location in y dimension: ";
    std::cin >> oy;

    // Calculate time increment based on Courant number and maximum layer velocities
    dt = (c * dx) / std::max({v1, v2, v3, v4});

    // Calculate number of time steps, and grid dimensions
    nt = ceil(t / dt);
    nx = floor(lx / dx);
    ny = floor(ly / dy);
}


/**
 * @brief Extract parameters for a simulation model from the "runfile.txt" configuration file.
 *
 * This function reads simulation parameters from a configuration file named "runfile.txt".
 * The file should contain 13 numeric values representing various parameters required for the simulation model.
 * The extracted parameters include dimensions, grid spacing, time, Courant number, layer velocities,
 * source properties, and more.
 *
 * @note The file format should be as follows:
 * ```
 * <lx> <ly> <dx> <t> <c> <v1> <v2> <v3> <v4> <a> <s> <ox> <oy>
 * ```
 *
 * @details The function opens "runfile.txt" and reads the parameters. If the file opening fails,
 * it displays an error message and exits the program. The loaded parameters are then used to set the
 * corresponding variables for the simulation model, such as lx, ly, dx, dy, t, c, v1, v2, v3, v4, a, s, ox, oy.
 * Additionally, the time increment (dt) and grid dimensions (nt, nx, ny) are calculated based on the loaded values.
 *
 * @see You should ensure that the "runfile.txt" file exists and is properly formatted.
 */
void load_parameters(void)
{
    std::ifstream Runfile;
    Runfile.open("runfile.txt");
    if (Runfile.fail()) // Check whether "runfile.txt" has been opened successfully
    {
        std::cout << "Error: Could not open 'runfile.txt'" << std::endl;
        exit(0); // Abort program
    }

    double parameters[13];
    for(int i = 0; i < 13; i++)
        Runfile >> parameters[i];
    Runfile.close(); // Close "runfile.txt"

    // Assign parameters to variables for the simulation model
    lx = parameters[0];
    ly = parameters[1];
    dx = parameters[2];
    dy = dx; // Assume the model is isotropic
    t = parameters[3];
    c = parameters[4];
    v1 = pow(parameters[5], 2); // Squared velocity appears in the wave equation
    v2 = pow(parameters[6], 2); 
    v3 = pow(parameters[7], 2);
    v4 = pow(parameters[8], 2);
    a = parameters[9];
    s = parameters[10];
    ox = parameters[11];
    oy = parameters[12];

    // Calculate time increment based on Courant number and maximum layer velocities
    dt = (c * dx) / std::max({v1, v2, v3, v4});

    // Calculate number of time steps, and grid dimensions
    nt = ceil(t / dt);
    nx = floor(lx / dx);
    ny = floor(ly / dy);
}


/**
 * @brief Performs basic checks to ensure that simulation parameters are valid and sensible.
 *
 * This function checks various simulation parameters to ensure that they meet certain criteria. 
 * If any of the checks fail, an error message is displayed, and the program
 * is terminated using the exit() function.
 *
 * The following checks are performed:
 * - lx and ly must be greater than 0.
 * - dx and dy must be greater than 0 and less than lx and ly, respectively.
 * - c must be greater than 0 and less than 1.
 * - ox and oy must be greater than 0 and less than lx and ly, respectively.
 * - s must be greater than 2*dx.
 * - t must be greater than 0.
 *
 * @see It is assumed that the relevant parameters (lx, ly, dx, dy, c, ox, oy, s, t) are declared
 *      and defined in the program's scope before calling this function.
 */
void check_parameters(void)
{
    if (lx <= 0 || ly <= 0)
    {
        std::cout << "Error: lx and ly must be greater than 0" << std::endl;
        exit(0);
    }
    else if (dx <= 0 || dx >= lx || dy <= 0 || dy >= ly)
    {
        std::cout << "Error: dx and dy must be greater than 0 and less than lx and ly, respectively" << std::endl;
        exit(0);
    }
    else if (c <= 0 || c >= 1)
    {
        std::cout << "Error: c must be greater than 0 and less than 1" << std::endl;
        exit(0);
    }
    else if (ox <= 0 || ox >= lx || oy <= 0 || oy >= ly)
    {
        std::cout << "Error: ox and oy must be greater than 0 and less than lx and ly, respectively" << std::endl;
        exit(0);
    }
    else if (s <= 2 * dx)
    {
        std::cout << "Error: s must be greater than 2*dx" << std::endl;
        exit(0);
    }
    else if (t <= 0)
    {
        std::cout << "Error: t must be greater than 0" << std::endl;
        exit(0);
    }
}


/**
 * @brief Loads model and simulation parameters from "runfile.txt" or prompts the user for manual entry.
 *
 * This function provides a user-friendly interface for obtaining simulation parameters. It prompts the user
 * with the option to either manually enter parameters or load them from a configuration file named "runfile.txt".
 *
 * @note The function uses the enter_parameters() function if the user chooses to manually enter parameters ("y").
 *       It uses the load_parameters() function if the user chooses to load parameters from "runfile.txt" ("n").
 *       Additionally, it performs basic checks on the loaded or manually entered parameters using check_parameters().
 *
 * @see The functions enter_parameters(), load_parameters(), and check_parameters() are assumed to be
 *      declared and defined in the program's scope before calling this function.
 */
void get_parameters(void)
{
    std::cout << "Manually enter runfile parameters? (y/n): ";
    std::cin >> parameter_loading_method;

    if (parameter_loading_method == "y")
        enter_parameters();
    else if (parameter_loading_method == "n")
        load_parameters();
    else
    {
        std::cout << "Invalid input. Please enter 'y' or 'n'." << std::endl;
        exit(0);
    }

    // Perform basic checks on the loaded or manually entered parameters
    check_parameters();
}


/**
 * @brief Allocates memory for a 2D array of double values.
 *
 * This function dynamically allocates memory for a 2D array of doubles with dimensions nx times ny.
 * The memory is allocated as an array of pointers to arrays, creating a matrix structure.
 *
 * @param nx The number of rows in the 2D array.
 * @param ny The number of columns in the 2D array.
 * @return   A pointer to the dynamically allocated 2D array.
 *
 * @warning It is the responsibility of the caller to free the allocated memory.
 *
 * @see free_memory()
 */
double** allocate_2d_array(int nx, int ny)
{
    // Allocate memory for the array of pointers to arrays
    double** array = new double*[nx];

    // Allocate memory for each row in the 2D array
    for(int i = 0; i < nx; i++)
        array[i] = new double[ny];

    return array;
}


/**
 * @brief Fills a 2D dynamically allocated array with zeros.
 *
 * This function takes a dynamically allocated 2D array of doubles with dimensions nx times ny
 * and sets each element to zero, effectively initializing the entire array with zeros.
 *
 * @param array A pointer to the dynamically allocated 2D array to be filled with zeros.
 * @param nx    The number of rows in the 2D array.
 * @param ny    The number of columns in the 2D array.
 *
 * @note The function assumes that the memory for the 2D array has been properly allocated before calling.
 *       It does not perform any allocation and only initializes the values in the existing array.
 */
void zeros(double** array, int nx, int ny)
{
    // Fill the 2D array with zeros
    for(int i = 0; i < nx; i++)
        for(int j = 0; j < ny; j++)
            array[i][j] = 0.0;
}


/**
 * @brief Fills a 2D array with four equally spaced variable velocities.
 *
 * This function populates a 2D array with four layers of equally spaced variable velocities along the vertical direction.
 * The velocities vary from v1 (top layer) to v4 (bottom layer) with equal spacing between the layers.
 *
 * @param v_model A pointer to the 2D array to be filled with velocity values.
 * @param nx      The number of rows (horizontal dimension) in the 2D array.
 * @param ny      The number of columns (vertical dimension) in the 2D array.
 * @param v1      Velocity of the top layer.
 * @param v2      Velocity of the second layer.
 * @param v3      Velocity of the third layer.
 * @param v4      Velocity of the bottom layer.
 *
 * @note The function assumes that the memory for the 2D array has been properly allocated before calling.
 *       It does not perform any allocation and only initializes the velocity values in the existing array.
 */
void four_layer_vel(double** v_model, int nx, int ny, double v1, double v2, double v3, double v4)
{
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            if (j <= ny / 4.)
                v_model[i][j] = v1;
            else if (j <= (2 * ny) / 4.)
                v_model[i][j] = v2;
            else if (j <= (3 * ny) / 4.)
                v_model[i][j] = v3;
            else
                v_model[i][j] = v4;
        }
    }
}


/**
 * @brief Fills a 1D array with equally spaced numbers from 0 to lxy.
 *
 * This function populates a 1D array with equally spaced numbers ranging from 0 to lxy.
 * The array represents a set of coordinates along the x or y axis with uniform spacing.
 *
 * @param array A pointer to the 1D array to be filled with coordinates.
 * @param nxy   The length of the 1D array.
 * @param lxy   The maximum value for the coordinates (inclusive).
 *
 * @note The function assumes that the memory for the 1D array has been properly allocated before calling.
 *       It does not perform any allocation and only initializes the coordinate values in the existing array.
 */
void fill_xy(double* array, int nxy, double lxy)
{
    for (int i = 0; i < nxy; i++)
    {
        array[i] = i * (lxy / (nxy - 1.0));
    }
}


/**
 * @brief Calculates the Ricker wavelet value for a given model location.
 *
 * This function computes the Ricker wavelet value for a specified location (x, y) in a 2D model.
 * The Ricker wavelet is commonly used in seismic exploration to represent a wavelet with a central peak.
 *
 * @param x  The x-coordinate of the model location.
 * @param y  The y-coordinate of the model location.
 * @param ox The x-coordinate of the wavelet's origin.
 * @param oy The y-coordinate of the wavelet's origin.
 * @param a  The amplitude of the wavelet.
 * @param s  The spread or width of the wavelet.
 * @return   The Ricker wavelet value for the given location.
 */
double ricker(double x, double y, double ox, double oy, double a, double s)
{
    double distance_squared = pow((x - ox), 2) + pow((y - oy), 2);
    double exponent_term = exp(-(distance_squared) / (2.0 * pow(s, 2)));

    return (a / (M_PI * pow(s, 2))) * (1 - 0.5 * (distance_squared / pow(s, 2))) * exponent_term;
}


/**
 * @brief Sets initial conditions for a wave propagation simulation.
 *
 * This function initializes the simulation by setting the initial conditions:
 * 1) \( u(x, y) \) is set to a Ricker wavelet at \( t = 0 \).
 * 2) \( \frac{\partial u}{\partial t}(x, y) \) is set to zero at \( t = 0 \),
 *    and \( u_{\text{old}}(t = 0) \) is determined using this condition.
 *
 * @param u       A pointer to the current state of the wavefield.
 * @param u_old   A pointer to the previous state of the wavefield.
 * @param v_model A pointer to the 2D model of velocity values.
 * @param x       A pointer to the array of x-coordinates.
 * @param y       A pointer to the array of y-coordinates.
 * @param nx      The number of grid points in the x-direction.
 * @param ny      The number of grid points in the y-direction.
 * @param ox      The x-coordinate of the wavelet's origin.
 * @param oy      The y-coordinate of the wavelet's origin.
 * @param a       The amplitude of the Ricker wavelet.
 * @param s       The spread or width of the Ricker wavelet.
 * @param dx      The grid spacing in the x-direction.
 * @param dy      The grid spacing in the y-direction.
 * @param dt      The time increment for the simulation.
 *
 * @note The function assumes that the memory for u, u_old, v_model, x, and y has been properly allocated before calling.
 *       It does not perform any allocation and only initializes the values in the existing arrays.
 */
void initial_conditions(double** u, double** u_old, double** v_model, double* x, double* y, int nx, int ny,
                        double ox, double oy, double a, double s, double dx, double dy, double dt)
{
    // Set u(t = 0) equal to the Ricker wavelet function across the whole model domain
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            u[i][j] = ricker(x[i], y[j], ox, oy, a, s);
        }
    }

    // Set u_old(t = 0) expression across the model inner domain
    for (int i = 1; i < nx - 1; i++)
    {
        for (int j = 1; j < ny - 1; j++)
        {
            u_old[i][j] = u[i][j] +
                          0.5 * (pow((dt / dx), 2) * (((v_model[i + 1][j] + v_model[i][j]) / 2.) * (u[i + 1][j] - u[i][j]) -
                                                      ((v_model[i - 1][j] + v_model[i][j]) / 2.) * (u[i][j] - u[i - 1][j])) +
                                 pow((dt / dy), 2) * (((v_model[i][j + 1] + v_model[i][j]) / 2.) * (u[i][j + 1] - u[i][j]) -
                                                      ((v_model[i][j - 1] + v_model[i][j]) / 2.) * (u[i][j] - u[i][j - 1])));
        }
    }
}


/**
 * @brief Advances the simulation incrementally in time using the finite-difference method.
 *
 * This function performs a time loop to update the wavefield at each grid point in the model over a specified time period.
 * It uses the finite-difference method to numerically solve the wave equation.
 *
 * @param u_new   A pointer to the wavefield at the next time step.
 * @param u       A pointer to the current state of the wavefield.
 * @param u_old   A pointer to the previous state of the wavefield.
 * @param v_model A pointer to the 2D model of velocity values.
 * @param nx      The number of grid points in the x-direction.
 * @param ny      The number of grid points in the y-direction.
 * @param nt      The total number of time steps.
 * @param dx      The grid spacing in the x-direction.
 * @param dy      The grid spacing in the y-direction.
 * @param dt      The time increment for the simulation.
 *
 * @note The function assumes that the memory for u_new, u, u_old, and v_model has been properly allocated before calling.
 *       It does not perform any allocation and only updates the values in the existing arrays during the time loop.
 */
void time_loop(double** u_new, double** u, double** u_old, double** v_model, int nx, int ny, int nt, double dx, double dy, double dt)
{
    // Begin time loop
    for (int k = 0; k < nt + 1; k++)
    {
        std::cout << "time-step " << k << " of " << nt << "\n";

        // Update all inner points using the finite-difference method
        for (int i = 1; i < nx - 1; i++)
        {
            for (int j = 1; j < ny - 1; j++)
            {
                u_new[i][j] = 2 * u[i][j] - u_old[i][j] +
                               (pow((dt / dx), 2) * (((v_model[i + 1][j] + v_model[i][j]) / 2.) * (u[i + 1][j] - u[i][j]) -
                                                     ((v_model[i - 1][j] + v_model[i][j]) / 2.) * (u[i][j] - u[i - 1][j])) +
                                pow((dt / dy), 2) * (((v_model[i][j + 1] + v_model[i][j]) / 2.) * (u[i][j + 1] - u[i][j]) -
                                                     ((v_model[i][j - 1] + v_model[i][j]) / 2.) * (u[i][j] - u[i][j - 1])));
            }
        }

        // Initialize u_old, u, u_new for the next time step
        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                u_old[i][j] = u[i][j];
                u[i][j] = u_new[i][j];
            }
        }
    }
}


/**
 * @brief Writes the final array to a file named "FDWE.txt".
 *
 * This function takes the final 2D array 'u' and writes its values to a file named "FDWE.txt" in the working directory.
 *
 * @param nx The number of rows in the 2D array.
 * @param ny The number of columns in the 2D array.
 * @param u  A pointer to the 2D array to be written to the file.
 */
void write_file(int nx, int ny, double** u)
{
    std::ofstream myfile;
    myfile.open("FDWE.txt");

    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            myfile << u[i][j];
            myfile << ", "; // delimiter = ", "
        }
        myfile << "\n"; // insert a new line between each row of u
    }

    myfile.close();
}


/**
 * @brief Frees all dynamically allocated memory.
 *
 * This function frees the memory allocated for 2D arrays 'u_new', 'u', 'u_old', 'v_model', and 1D arrays 'x' and 'y'.
 *
 * @param u_new   A pointer to the 2D array 'u_new'.
 * @param u       A pointer to the 2D array 'u'.
 * @param u_old   A pointer to the 2D array 'u_old'.
 * @param v_model A pointer to the 2D array 'v_model'.
 * @param x       A pointer to the 1D array 'x'.
 * @param y       A pointer to the 1D array 'y'.
 * @param nx      The size of the outermost level of the arrays.
 */
void free_memory(double** u_new, double** u, double** u_old, double** v_model, double* x, double* y, int nx)
{
    for (int i = 0; i < nx; i++) // for 2D arrays
    {
        delete[] u_new[i];
        delete[] u[i];
        delete[] u_old[i];
        delete[] v_model[i];
    }

    delete[] u_new;
    delete[] u;
    delete[] u_old;
    delete[] v_model;
    delete[] x; // 1D arrays
    delete[] y;
}


int main(void)
{
    // Load model and source wavelet parameters from runfile.txt or manually enter them
    get_parameters();

    // Dynamically allocate memory for 2D arrays u_new, u, u_old, v_model, and 1D arrays x, y
    double** u_new = allocate_2d_array(nx, ny);
    double** u = allocate_2d_array(nx, ny);
    double** u_old = allocate_2d_array(nx, ny);
    double** v_model = allocate_2d_array(nx, ny);
    double* x = new double[nx];
    double* y = new double[ny];

    // Fill u_new, u, u_old arrays with zeros
    zeros(u_new, nx, ny);
    zeros(u, nx, ny);
    zeros(u_old, nx, ny);

    // Fill v_model with four equally spaced variable velocities from v1 (top) to v4 (bottom)
    four_layer_vel(v_model, nx, ny, v1, v2, v3, v4);

    // Fill x, y from 0 to lxy with nxy values
    fill_xy(x, nx, lx);
    fill_xy(y, ny, ly);

    // Set initial conditions: 1) u(x,y) = ricker(x,y) at t = 0
    //                         2) use du/dt (x,y) = 0 at t = 0 to find the expression for u_old(t = 0)
    initial_conditions(u, u_old, v_model, x, y, nx, ny, ox, oy, a, s, dx, dy, dt);

    // Enter time-stepping loop to propagate wave
    time_loop(u_new, u, u_old, v_model, nx, ny, nt, dx, dy, dt);

    // Write the final array 'u' to the file "FDWE.txt"
    write_file(nx, ny, u);

    // Free all dynamically allocated memory
    free_memory(u_new, u, u_old, v_model, x, y, nx);

    return 0;
}
