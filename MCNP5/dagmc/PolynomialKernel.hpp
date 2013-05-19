// MCNP5/dagmc/PolynomialKernel.hpp

#ifndef DAGMC_POLYNOMIAL_KERNEL_H
#define DAGMC_POLYNOMIAL_KERNEL_H

#include <vector>

#include "KDEKernel.hpp"

/**
 * \class PolynomialKernel
 * \brief Defines a polynomial kernel function
 * 
 * PolynomialKernel is a Derived class that can create any polynomial kernel
 * function that is defined by a smoothness factor s and kernel order 2r.  The
 * domain of these kernel functions is limited to u = [-1, 1].
 *
 * Common polynomial kernels include uniform (s = 0), epanechnikov (s = 1),
 * biweight (s = 2), and triweight (s = 3).  All of these cases and more can be
 * constructed using the following general formula
 *
 *               (1/2)_s+1
 *      K_s(u) = --------- (1 - u^2)^s
 *                   s!
 *
 * where (1/2)_s+1 is known as the Pochhammer symbol
 *
 *      (x)_n = x(x + 1)(x + 2)...(x + n - 1)
 *
 * Any kernel function based on K_s(u) is said to be of 2nd-order.  Higher-
 * order cases (r > 1) can be derived from these 2nd-order kernels by simply
 * multiplying them with a second polynomial
 *
 *      K_2r,s(u) = B_r,s(u) * K_s(u)
 *
 * For more information on how to construct higher-order kernel functions,
 * please refer to the following reference
 * 
 *     B. E. Hansen, "Exact Mean Integrated Squared Error of Higher Order
 *     "Kernel Estimators", Econometric Theory, 21, pp. 1031-1057 (2005)
 */
class PolynomialKernel : public KDEKernel
{
  public:
    /**
     * \brief Constructor
     * \param s the level of smoothness of the kernel
     * \param r defines a kernel of 2rth-order
     */
    PolynomialKernel(unsigned int s, unsigned int r);

    /**
     * \brief Destructor
     */
    ~PolynomialKernel();

    // >>> DERIVED PUBLIC INTERFACE from KDEKernel.hpp

    /**
     * \brief evaluate the kernel function
     * \param u the value at which the kernel will be evaluated
     * \return K_2r,s(u)
     */
    virtual double evaluate(double u) const;

    /**
     * \brief get_kernel_name()
     * \return string representing kernel name
     */
    virtual std::string get_kernel_name() const;

    /**
     * \brief integrates the ith moment function for this polynomial kernel
     * \param a, b the lower and upper integration limits
     * \param i the index representing the ith moment function
     * \return integral of the ith moment function evaluated from a to b
     */
    virtual double integrate_moment(double a, double b, unsigned int i) const;
  
  private:
    /// Smoothness factor for this kernel
    unsigned int s;

    /// Related to the order of this kernel (order = 2r)
    unsigned int r;

    /// Represents constant multiplier of kernel function
    const double multiplier;

    /// Coefficients of the polynomial generated for kernels of order > 2
    std::vector<double> coefficients;

    /// Quadrature set for integrating moment functions
    Quadrature* quadrature;

    // >>> PRIVATE FUNCTIONS

    /**
     * \brief sets the common multiplier term for the kernel function
     */
    double set_multiplier();

    /**
     * \brief evaluates the Pochhammer symbol
     * \param x the value for which to evaluate the Pochhammer symbol
     * \param n any non-negative integer (n >= 0)
     * \return (x)_n
     *
     * Computes (x)_n = x(x + 1)(x + 2)...(x + n - 1) where (x)_0 = 1.
     */
    double pochhammer(double x, unsigned int n);

    /**
     * \brief evaluates the factorial function
     * \param n any non-negative integer (n >=0)
     * \return n!
     *
     * Computes n! = n(n-1)(n-2)...(2)(1) where 0! = 1.
     */
    long double factorial(unsigned int n);
};

#endif // DAGMC_POLYNOMIAL_KERNEL_H

// end of MCNP5/dagmc/PolynomialKernel.hpp
