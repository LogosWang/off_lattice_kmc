#include "Site.h" // Always include the corresponding header file
#include <cmath>  // Include cmath if you might use math functions (e.g., for complex moves)

// --- Constructor Definition ---
// This defines how a new Site object is created and its members are initialized.
Site::Site(int site_id, double initial_x, double initial_y, double initial_z, int site_type)
    : id(site_id), x(initial_x), y(initial_y), z(initial_z), type(site_type)
{
    // Any additional setup logic for a new site can go here.
    // For simple initialization, the initializer list (after the colon) is efficient.
}

// --- Move Function Definition ---
// This defines how a Site object changes its position.
// 'direction_axis': Typically 0 for X, 1 for Y, 2 for Z.
// 'step_size': The magnitude of the displacement along the chosen axis.
void Site::move(std::vector<double> direction, double step_size) {
    x += step_size*direction[0];
    y += step_size*direction[1];
    z += step_size*direction[2];
}