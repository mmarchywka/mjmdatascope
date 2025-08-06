#include <iostream>     // For input/output operations (e.g., std::cout, std::cerr)
#include <vector>       // For dynamic arrays (used temporarily for dense coefficient representation)
#include <complex>      // For complex number arithmetic (std::complex<double>)
#include <map>          // For storing polynomial coefficients efficiently (degree -> coefficient)
#include <cmath>        // For mathematical functions like std::abs (for magnitude of complex numbers)
#include <algorithm>    // For algorithms like std::find_if (though not directly used in final version)
#include <cstdlib>      // For std::rand and std::srand (for random number generation)
#include <ctime>        // For std::time (to seed the random number generator)

namespace gemini_root_ns {
// Define a small epsilon value for floating-point comparisons.
// Numbers smaller than this are considered effectively zero.
double EPSILON = 1e-9;
void set_noise_floor(const double & eps) {EPSILON=eps; } 

/**
 * @class Polynomial
 * @brief Represents a polynomial with complex coefficients.
 *
 * This class stores polynomial coefficients in a sparse manner using a std::map,
 * where the key is the degree and the value is the complex coefficient.
 * It provides methods for adding terms, evaluating the polynomial, computing its derivative,
 * and performing polynomial deflation (division by a linear factor).
 */
class Polynomial {
public:
// MEMBERS
    // Stores coefficients: map from degree (int) to complex coefficient (std::complex<double>).
    std::map<int, std::complex<double>> coeffs;
    // The highest degree of the polynomial.
    int degree;
bool m_real_coefs;
    /**
     * @brief Constructor for the Polynomial class.
     * Initializes a zero polynomial with degree 0.
     */
    Polynomial() : degree(0),m_real_coefs(true) {}

    /**
     * @brief Adds a term (coefficient * x^degree) to the polynomial.
     * @param deg The degree of the term.
     * @param coeff The complex coefficient of the term.
     * If the coefficient is effectively zero, the term is removed if it exists.
     */
    void add_term(int deg, std::complex<double> coeff) {
// mjm want to accumulate, could zero out... 
	coeffs[deg]+=coeff;
       /* // Check if the coefficient is effectively zero
        if (std::abs(coeff.real()) < EPSILON && std::abs(coeff.imag()) < EPSILON) {
            // If it's zero and the term exists, remove it
            if (coeffs.count(deg)) {
                coeffs.erase(deg);
            }
        } else {
            // Otherwise, add or update the term
            coeffs[deg] = coeff;
        } */
        // Always update the polynomial's degree after modifying coefficients
        update_degree();
    }

    /**
     * @brief Updates the polynomial's degree based on the highest non-zero coefficient.
     * If all coefficients are zero, the degree is set to 0.
// really just check updated coef... 
     */
    void update_degree() {
        degree = 0;
m_real_coefs=true;
        if (coeffs.empty()) {
            return; // Empty map means zero polynomial, degree 0
        }

        int max_deg = 0;
        bool has_non_zero_coeff = false;
        // Iterate through all stored coefficients to find the highest non-zero one
        for (auto const& [deg, coeff] : coeffs) {
            const double  imabs= std::abs(coeff.imag());
            if (std::abs(coeff.real()) > EPSILON || std::abs(coeff.imag()) > EPSILON) {
                if (deg > max_deg) { max_deg = deg; }
				if (imabs!=0) { m_real_coefs=false; }
                has_non_zero_coeff = true;
            }
        }
        // Set the degree. If no non-zero coefficients, it's a zero polynomial (degree 0).
        degree = has_non_zero_coeff ? max_deg : 0;
        // If it's effectively a zero polynomial, clear the map for consistency.
        if (!has_non_zero_coeff) {
            coeffs.clear();
        }
    }

    /**
     * @brief Evaluates the polynomial at a given complex value x.
     * @param x The complex value at which to evaluate the polynomial.
     * @return The complex result of P(x).
     */
    std::complex<double> evaluate(std::complex<double> x) const {
        std::complex<double> result = 0.0;
        // Iterate through all terms and sum them up
        for (auto const& [deg, coeff] : coeffs) {
            result += coeff * std::pow(x, deg);
        }
        return result;
    }

    /**
     * @brief Computes the derivative of the polynomial.
     * @return A new Polynomial object representing the derivative P'(x).
     */
    Polynomial derivative() const {
        Polynomial p_prime;
        // For each term a_k * x^k, its derivative is (k * a_k) * x^(k-1)
        for (auto const& [deg, coeff] : coeffs) {
            if (deg > 0) { // Derivative of a constant term (deg 0) is 0
                p_prime.add_term(deg - 1, coeff * static_cast<double>(deg));
            }
        }
        return p_prime;
    }

    /**
     * @brief Performs polynomial deflation: divides P(x) by (x - root).
     * This method uses synthetic division to find the quotient Q(x) such that P(x) = (x - root)Q(x) + R.
     * If 'root' is an actual root, the remainder R should be zero.
     * @param root The complex root to divide by.
     * @return A new Polynomial object representing the quotient Q(x).
     */
    Polynomial deflate(std::complex<double> root) const {
        if (degree == 0) {
            return Polynomial(); // Cannot deflate a constant polynomial
        }

        Polynomial Q;
        // Create a temporary dense vector for coefficients of P(x) for synthetic division.
        // The size is degree + 1, where index k corresponds to coefficient of x^k.
        std::vector<std::complex<double>> a(degree + 1, 0.0);
        for (auto const& [deg, coeff] : coeffs) {
            a[deg] = coeff;
        }

        // Vector to store coefficients for Q(x). Q(x) will have degree - 1.
        std::vector<std::complex<double>> b(degree, 0.0);

        // Synthetic division:
        // The highest degree coefficient of Q(x) (b_{n-1}) is the same as P(x)'s highest (a_n).
        b[degree - 1] = a[degree];

        // Calculate remaining coefficients of Q(x) by iterating backwards.
        // b_k = a_{k+1} + root * b_{k+1}
        for (int k = degree - 2; k >= 0; --k) {
            b[k] = a[k + 1] + root * b[k + 1];
        }

        // Add the calculated coefficients to the new Polynomial object Q.
        for (int k = 0; k < degree; ++k) {
            Q.add_term(k, b[k]);
        }
    // Q.print();
        return Q;
    }

    /**
     * @brief Prints the polynomial in a human-readable format.
     * Displays terms from highest degree to lowest.
     */
    void print() const {
        // Handle the case of a zero polynomial
        if (coeffs.empty() || (degree == 0 && coeffs.count(0) && std::abs(coeffs.at(0).real()) < EPSILON && std::abs(coeffs.at(0).imag()) < EPSILON)) {
            std::cout << "0" << std::endl;
            return;
        }
        bool first_term = true;
        // Iterate from the highest degree down to 0 for printing
        for (int deg = degree; deg >= 0; --deg) {
            // Check if a term exists at this degree and its coefficient is non-zero
            if (coeffs.count(deg) && (std::abs(coeffs.at(deg).real()) > EPSILON || std::abs(coeffs.at(deg).imag()) > EPSILON)) {
                if (!first_term) {
                    // Add " + " or " " depending on the sign of the real part of the coefficient
                    if (coeffs.at(deg).real() >= 0 || (std::abs(coeffs.at(deg).real()) < EPSILON && coeffs.at(deg).imag() >= 0)) {
                        std::cout << " + ";
                    } else {
                        std::cout << " "; // Negative sign will be part of the complex number output
                    }
                }

                std::cout << coeffs.at(deg); // Print the complex coefficient
                if (deg == 1) {
                    std::cout << "x";
                } else if (deg > 1) {
                    std::cout << "x^" << deg;
                }
                first_term = false;
            }
        }
        std::cout << std::endl;
    }
};

/**
 * @brief Implements the Newton-Raphson method to find a single complex root of a polynomial.
 * @param p The polynomial for which to find a root.
 * @param initial_guess The starting complex guess for the root.
 * @param tolerance The maximum allowed magnitude of P(x) at the found root for convergence.
 * @param max_iterations The maximum number of iterations to perform.
 * @return The complex root found, or the best guess after max_iterations.
 */
std::complex<double> newton_raphson(
    const Polynomial& p,
    std::complex<double> initial_guess,
    double tolerance,
    int max_iterations
) {
    std::complex<double> x = initial_guess;
    Polynomial p_prime = p.derivative(); // Compute the derivative once

    for (int i = 0; i < max_iterations; ++i) {
        std::complex<double> fx = p.evaluate(x);
        std::complex<double> f_prime_x = p_prime.evaluate(x);

        // Check for convergence: if P(x) is close to zero
        if (std::abs(fx.real()) < tolerance && std::abs(fx.imag()) < tolerance) {
            return x;
        }

        // Avoid division by zero or very small derivative, which indicates a flat region
        // or a multiple root where Newton-Raphson struggles.
        if (std::abs(f_prime_x.real()) < EPSILON && std::abs(f_prime_x.imag()) < EPSILON) {
            // Derivative is zero or very small. Perturb the guess slightly to escape.
            x += std::complex<double>(0.001, 0.001); // Add a small complex perturbation
            continue; // Continue to the next iteration with the perturbed guess
        }

        // Newton-Raphson iteration formula: x_{n+1} = x_n - P(x_n) / P'(x_n)
        x = x - fx / f_prime_x;
    }
    return x; // Return the best guess found after max_iterations
}

/**
 * @brief Finds all complex roots of a polynomial using Newton-Raphson and polynomial deflation.
 * It iteratively finds one root, deflates the polynomial, and repeats until all roots are found
 * or the polynomial is reduced to a constant.
 * @param p The polynomial for which to find all roots. A copy is made as it will be modified.
 * @param tolerance The tolerance for root convergence in Newton-Raphson.
 * @param max_iterations The maximum iterations for each Newton-Raphson call.
 * @param num_initial_guesses_per_root The number of random initial guesses to try for each root.
 * @param guess_range The range for generating random initial guesses (e.g., -10.0 to 10.0 for real/imaginary parts).
 * @return A vector of complex numbers representing the found roots.
 */
std::vector<std::complex<double>> find_all_roots(
    Polynomial p, // Pass by value to allow modification (deflation)
    double tolerance = 1e-7,
    int max_iterations = 1000,
    int num_initial_guesses_per_root = 20,
    double guess_range = 10.0
) {
    std::vector<std::complex<double>> roots;
    unsigned int  original_degree = p.degree; // mjm unsigned

    // Handle constant polynomial (degree 0). A non-zero constant has no roots.
    // A zero polynomial (P(x) = 0) has all complex numbers as roots, which is not what
    // this function typically aims to find.
    if (original_degree == 0) {
        return roots;
    }

// mjm find conjugate if poly coefs are real ... doh should stay real that wy 

    // Seed the random number generator using the current time for varied guesses.
    std::srand(static_cast<unsigned int>(std::time(nullptr)));

    // Loop until the polynomial is reduced to a constant (degree 0) or we have found
    // as many roots as the original polynomial's degree.
    while (p.degree > 0 && roots.size() < original_degree) {
        std::complex<double> current_root;
        bool root_found_and_unique = false;

        // Try multiple random initial guesses for the current polynomial.
        // This heuristic helps in finding different roots and avoiding local minima.
        for (int i = 0; i < num_initial_guesses_per_root; ++i) {
            // Generate random real and imaginary parts within the specified guess_range.
            double rand_real = (static_cast<double>(rand()) / RAND_MAX) * (2 * guess_range) - guess_range;
            double rand_imag = (static_cast<double>(rand()) / RAND_MAX) * (2 * guess_range) - guess_range;
            std::complex<double> guess(rand_real, rand_imag);

            // Attempt to find a root using Newton-Raphson with the current guess.
            std::complex<double> found_root = newton_raphson(p, guess, tolerance, max_iterations);

            // Verify if the found_root is indeed a root of the current polynomial.
            if (std::abs(p.evaluate(found_root).real()) < tolerance && std::abs(p.evaluate(found_root).imag()) < tolerance) {
                // Check if this root is a duplicate of an already found root (due to convergence to the same root).
                bool is_duplicate = false;
                for (const auto& existing_root : roots) {
                    // Compare roots with a small tolerance
                    if (std::abs(found_root - existing_root) < tolerance) {
                        is_duplicate = true;
                        break;
                    }
                }

                if (!is_duplicate) {
                    current_root = found_root;
                    root_found_and_unique = true;
                    break; // Found a new, unique root, exit guess loop
                }
            }
        }

        if (!root_found_and_unique) {
            // If no new unique root could be found after all guesses, issue a warning and stop.
            // This can happen due to numerical instability, roots outside the guess range,
            // or limitations of the Newton-Raphson method for certain polynomial characteristics.
            std::cerr << "Warning: Could not find a new unique root for remaining polynomial of degree " << p.degree << ". Stopping." << std::endl;
            break;
        }

        // Add the newly found unique root to the list.
        roots.push_back(current_root);

        // Deflate the polynomial by dividing it by (x - current_root).
        // This reduces the polynomial's degree, allowing us to find the next root.
        p = p.deflate(current_root);
    }

    return roots;
}

// Main function to demonstrate the polynomial root finder
int main() {
    std::cout << "--- Polynomial Root Finder ---" << std::endl;

    // Example 1: A simple quadratic polynomial: x^2 - 4 = 0 (roots: 2, -2)
    std::cout << "\nExample 1: P(x) = x^2 - 4" << std::endl;
    Polynomial p1;
    p1.add_term(2, {1.0, 0.0}); // 1 * x^2
    p1.add_term(0, {-4.0, 0.0}); // -4 * x^0
    std::cout << "Polynomial: ";
    p1.print();
    std::vector<std::complex<double>> roots1 = find_all_roots(p1);
    std::cout << "Roots found (" << roots1.size() << "):" << std::endl;
    for (const auto& root : roots1) {
        std::cout << "  " << root << std::endl;
    }

    // Example 2: A cubic polynomial with complex roots: x^3 + 1 = 0 (roots: -1, 0.5+0.866i, 0.5-0.866i)
    std::cout << "\nExample 2: P(x) = x^3 + 1" << std::endl;
    Polynomial p2;
    p2.add_term(3, {1.0, 0.0}); // 1 * x^3
    p2.add_term(0, {1.0, 0.0}); // 1 * x^0
    std::cout << "Polynomial: ";
    p2.print();
    std::vector<std::complex<double>> roots2 = find_all_roots(p2);
    std::cout << "Roots found (" << roots2.size() << "):" << std::endl;
    for (const auto& root : roots2) {
        std::cout << "  " << root << std::endl;
    }

    // Example 3: A high-degree sparse polynomial: x^100 + x^5 - 2 = 0
    // This will be challenging for numerical methods due to high degree.
    std::cout << "\nExample 3: P(x) = x^100 + x^5 - 2" << std::endl;
    Polynomial p3;
    p3.add_term(100, {1.0, 0.0}); // 1 * x^100
    p3.add_term(5, {1.0, 0.0});   // 1 * x^5
    p3.add_term(0, {-2.0, 0.0});  // -2 * x^0
    std::cout << "Polynomial: ";
    p3.print();
    std::vector<std::complex<double>> roots3 = find_all_roots(p3, 1e-7, 2000, 50, 5.0); // Adjusted parameters
    std::cout << "Roots found (" << roots3.size() << "):" << std::endl;
    for (const auto& root : roots3) {
        std::cout << "  " << root << std::endl;
    }

    // Example 4: A polynomial with a multiple root: (x-1)^2 = x^2 - 2x + 1 = 0 (root: 1, 1)
    std::cout << "\nExample 4: P(x) = x^2 - 2x + 1" << std::endl;
    Polynomial p4;
    p4.add_term(2, {1.0, 0.0});
    p4.add_term(1, {-2.0, 0.0});
    p4.add_term(0, {1.0, 0.0});
    std::cout << "Polynomial: ";
    p4.print();
    std::vector<std::complex<double>> roots4 = find_all_roots(p4);
    std::cout << "Roots found (" << roots4.size() << "):" << std::endl;
    for (const auto& root : roots4) {
        std::cout << "  " << root << std::endl;
    }

    // Example 5: A polynomial with a high degree and few terms, testing max degree
    // P(x) = x^500 + 1 (roots are 500th roots of -1)
    std::cout << "\nExample 5: P(x) = x^500 + 1" << std::endl;
    Polynomial p5;
    p5.add_term(500, {1.0, 0.0});
    p5.add_term(0, {1.0, 0.0});
    std::cout << "Polynomial: ";
    p5.print();
    // For very high degrees, finding all roots can be very time-consuming and numerically challenging.
    // The random guess range and number of guesses might need to be adjusted significantly.
    // This example might not find all 500 roots reliably without more advanced techniques.
    std::vector<std::complex<double>> roots5 = find_all_roots(p5, 1e-7, 2000, 50, 2.0); // Roots of x^N + C are typically within a radius of |C|^(1/N)
    std::cout << "Roots found (" << roots5.size() << "):" << std::endl;
    for (const auto& root : roots5) {
        std::cout << "  " << root << std::endl;
    }

    return 0;
}
}; // gemin_root_ns
