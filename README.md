Ex1 – Polynomial Utilities (Ariel University, 2026)
Welcome to my GitHub repository for Ex1 – Arrays, Static Functions & JUnit, part of the Introduction to Computer Science course at Ariel University. This project implements a full set of utility methods for handling polynomials represented as arrays of doubles, where each index corresponds to a coefficient.

📚 Function Overview

1. f(double[] poly, double x) Evaluates the polynomial at a given value of x using Horner’s method.
2. root_rec(double[] p, double x1, double x2, double eps) Recursively finds a root of the polynomial in the interval [x1, x2] using the bisection method. Assumes the polynomial changes sign on the interval.
3. PolynomFromPoints(double[] xx, double[] yy) Builds a polynomial (degree ≤ 2) that passes through 2 or 3 given points. Returns null if there is no unique polynomial.
4. equals(double[] p1, double[] p2) Checks whether two polynomials are mathematically equal (within epsilon tolerance).
5. poly(double[] poly) Creates a readable string such as: -1.2x^2 + 3.1x + 2.0
6. sameValue(double[] p1, double[] p2, double x1, double x2, double eps) Finds an x where the two polynomials have (approximately) the same value. Uses the bisection method on f1(x) – f2(x).
7. length(double[] p, double x1, double x2, int n) Approximates the arc length of the polynomial curve between x1 and x2 using n line segments.
8. area(double[] p1, double[] p2, double x1, double x2, int n) Computes the approximate area between two polynomials over [x1, x2] using a trapezoid-based Riemann method. Handles crossings where the functions swap dominance.
9. getPolynomFromString(String s) Parses a polynomial string such as "-1x^2 + 3x + 2" back into its array representation.
10. add(double[] p1, double[] p2) Returns the coefficient-wise sum of two polynomials.
11. mul(double[] p1, double[] p2) Performs polynomial multiplication using convolution.
12. derivative(double[] p) Returns the derivative polynomial. Example: {2, 3, 4} → derivative {3, 8} (3 + 8x).
    
🤝 Contributions This repository represents my personal solution and improvements to the assignment. Suggestions, optimizations, or pull requests are welcome! 
⭐ Enjoy the Code! Thanks for visiting! Feel free to explore, learn, or reuse parts of the implementation. If you find this helpful, consider starring the repo ★

![Image_Alt](https://github.com/AnanstasSs777/I2CS_Ex1/blob/35a2a72e0712e9485efe74c108d84ec9ca8b33c8/WhatsApp%20Image%202025-11-29%20at%2018.10.33.jpeg) 
