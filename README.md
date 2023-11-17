
# Leverrier-Faddeev Method Implementation

## Overview

This repository contains a C++ program implementing the Leverrier-Faddeev method, a numerical algorithm used for finding the characteristic polynomial of a matrix. The program also includes the Faddeev method for solving polynomial equations and the Newton method for finding roots of polynomials.

## Usage

To use this program, you can follow these steps:

1. Clone the repository:

   ```bash
   git clone https://github.com/your-username/leverrier-faddeev-method.git
   ```

2. Compile the C++ program:

   ```bash
   g++ -o leverrier_faddeev leverrier_faddeev.cpp
   ```

3. Run the executable:

   ```bash
   ./leverrier_faddeev
   ```

4. Follow the on-screen prompts to enter the coefficients of the polynomial and observe the roots found using the Leverrier-Faddeev and Newton methods.

## Features

- **Matrix Operations**: The program provides functions for matrix operations, including matrix multiplication and matrix power.

- **Leverrier's Method**: The `leverrier` function calculates the coefficients of the characteristic polynomial using Leverrier's method.

- **Faddeev's Method**: The `fadeef` function computes the coefficients of the characteristic polynomial using Faddeev's method.

- **Newton's Method**: The `newton` function uses Newton's method to find roots of a given polynomial.

- **Sturm Sequence**: The program includes functions to compute the Sturm sequence of a polynomial.

- **Polynomial Division and Differentiation**: Functions for polynomial division and differentiation are provided to support various polynomial-related computations.

## Examples

The `main` function demonstrates the usage of the Leverrier-Faddeev method with a sample matrix. You can modify the matrix or integrate the methods into your own projects as needed.

## Performance Evaluation

The code includes a commented-out section for performance evaluation, comparing the execution times of Leverrier's and Faddeev's methods for random matrices of increasing sizes.

Feel free to uncomment and modify this section for further performance analysis.

