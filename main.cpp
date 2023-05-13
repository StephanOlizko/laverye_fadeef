#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <string>
#include <chrono>

using namespace std;


vector<vector<double>> identity_matrix(int n) {
    vector<vector<double>> I(n, vector<double>(n, 0.0));
    for (int i = 0; i < n; i++) {
        I[i][i] = 1.0;
    }
    return I;
}


vector<vector<double>>
add_scaled_matrices(vector<vector<double>> A, vector<vector<double>> B, double alpha_A, double alpha_B) {
    int n = A.size();
    int m = A[0].size();
    vector<vector<double>> C(n, vector<double>(m));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            C[i][j] = alpha_A * A[i][j] + alpha_B * B[i][j];
        }
    }
    return C;
}

vector<vector<double>> matrixMultiply(vector<vector<double>> &A, vector<vector<double>> &B) {
    int n = A.size();
    int m = A[0].size();
    int k = B[0].size();
    vector<vector<double>> res(n, vector<double>(k, 0));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < k; j++) {
            for (int p = 0; p < m; p++) {
                res[i][j] += A[i][p] * B[p][j];
            }
        }
    }
    return res;
}

vector<vector<double>> matrixPower(vector<vector<double>> &A, int t) {
    int n = A.size();
    vector<vector<double>> res(n, vector<double>(n, 0));
    // Заполняем диагональные элементы единицами
    for (int i = 0; i < n; i++) {
        res[i][i] = 1;
    }
    // Быстрое возведение матрицы в степень
    while (t > 0) {
        if (t % 2 == 1) {
            res = matrixMultiply(res, A);
        }
        A = matrixMultiply(A, A);
        t /= 2;
    }
    return res;
}

double trace(vector<vector<double>> A) {
    int n = A.size();
    double sum = 0;
    for (int i = 0; i < n; i++) {
        sum += A[i][i];
    }
    return sum;
}

void printMatrix(vector<vector<double>> matrix) {
    for (auto row: matrix) {
        for (auto elem: row) {
            cout << elem << " ";
        }
        cout << endl;
    }
}


vector<double> leverrier(vector<vector<double>> matrix) {
    int n = matrix.size();
    vector<double> res(n), s(n + 1);
    res[0] = pow(-1, n);
    vector<vector<double>> tmp_M = matrix;
    for (int i = 0; i <= n; i++) {
        s[i] = trace(tmp_M);
        tmp_M = matrixMultiply(tmp_M, matrix);
    }
    for (int i = 0; i < n; i++) {
        double sm = 0;
        for (int j = 0; j < i; j++) {
            sm = sm + res[j] * s[i - j - 1];
        }
        res[i] = (s[i] - sm) / (i + 1);
    }

    vector<double> q(n + 1);
    q[0] = pow(-1, n);
    for (int i = 1; i < q.size(); i++) {
        q[i] = res[i - 1] * pow(-1, n + 1);
    }
    return q;
}

vector<double> fadeef(vector<vector<double>> matrix) {
    int n = matrix.size();
    vector<double> res(n);
    vector<vector<vector<double>>> A(n + 1), B(n + 1);
    B[0] = identity_matrix(n);
    for (int i = 1; i <= n; i++) {
        A[i] = matrixMultiply(matrix, B[i - 1]);
        res[i - 1] = trace(A[i]) / i;
        B[i] = add_scaled_matrices(A[i], identity_matrix(n), 1, -res[i - 1]);
    }
    vector<double> q(n + 1);
    q[0] = pow(-1, n);
    for (int i = 1; i < q.size(); i++) {
        q[i] = res[i - 1] * pow(-1, n + 1);
    }
    return q;
}


vector<vector<double>> random_matrix(int n) {
    // Инициализируем генератор случайных чисел
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(-10.0, 10.0); // Равномерное распределение от -1.0 до 1.0

    vector<vector<double>> A(n, vector<double>(n));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A[i][j] = dis(gen); // Заполняем элемент случайным числом
        }
    }
    return A;
}


vector<double> roots;

// Define the polynomial to solve
double polynomial(double x, vector<double> coeffs) {
    double result = 0.0;
    int degree = coeffs.size() - 1;
    for (int i = 0; i < coeffs.size(); i++) {
        result += coeffs[i] * pow(x, (degree - i));
    }
    return result;
}

// Define the derivative of the polynomial
double derivative(double x, vector<double> coeffs) {
    double result = 0.0;
    int degree = coeffs.size() - 1;
    for (int i = 0; i < coeffs.size() - 1; i++) {
        result += (degree - i) * coeffs[i] * pow(x, (degree - i - 1));
    }
    return result;
}

// Get the coefficients of the polynomial from the user
vector<double> get_coefficients() {
    vector<double> coeffs;
    string input;

    cout << "Enter the coefficients of the polynomial (highest degree first), separated by spaces: ";
    getline(cin, input);

    int last_pos = 0;
    int pos = input.find(" ");
    while (pos != -1) {
        string str_coeff = input.substr(last_pos, pos - last_pos);
        double coeff = stod(str_coeff);
        coeffs.push_back(coeff);
        last_pos = pos + 1;
        pos = input.find(" ", last_pos);
    }

    string str_coeff = input.substr(last_pos);
    double coeff = stod(str_coeff);
    coeffs.push_back(coeff);

    return coeffs;
}

// Find the roots of the polynomial using Newton's method
void newton(vector<double> coeffs, double x0, double epsilon, int max_iterations) {
    double xi = x0;
    double xi_plus_1;
    double error = numeric_limits<double>::max(); // set error to a large value to begin with
    int num_iterations = 0;

    while (error > epsilon) {
        xi_plus_1 = xi - polynomial(xi, coeffs) / derivative(xi, coeffs); // calculate the next iteration
        error = abs(polynomial(xi_plus_1, coeffs)); // calculate the error
        xi = xi_plus_1; // update x_i to the next iteration
        num_iterations++;
        if (error < epsilon) roots.push_back(xi); // if the error is within the tolerance, add this root to the list
    }
}

vector<double> divide_polynomials(vector<double> dividend, vector<double> divisor) {
    vector<double> quotient;
    int dividend_degree = dividend.size() - 1;
    int divisor_degree = divisor.size() - 1;
    int quotient_size = dividend_degree - divisor_degree + 1;
    quotient.resize(quotient_size);
    vector<double> remainder = dividend;
    for (int i = dividend_degree; i >= divisor_degree; i--) {
        double factor = remainder[i] / divisor[divisor_degree];
        quotient[i - divisor_degree] = factor;
        for (int j = 0; j <= divisor_degree; j++) {
            remainder[i - j] -= factor * divisor[divisor_degree - j];
        }
    }
    return remainder;
}

vector<double> differentiation(vector<double> coeffs) {
    unsigned int n = coeffs.size() - 1;
    for (int i = 0; i < n; ++i) {
        coeffs[i] = (n - i) * coeffs[i];
    }
    coeffs.pop_back();
    return coeffs;
}

// Compute the Sturm sequence of a polynomial
vector<vector<double>> sturm_sequence(vector<double> coeffs) {
    vector<vector<double>> sturm_seq;
    sturm_seq.push_back(coeffs);
    sturm_seq.push_back(differentiation(coeffs));
    for (int i = 2; i < coeffs.size(); i++) {
        if (sturm_seq[sturm_seq.size() - 1].size() == 1) {
            break;
        }
        vector<double> sturm_i = divide_polynomials(sturm_seq[i - 2], sturm_seq[i - 1]);
        while ((sturm_i[sturm_i.size() - 1] == 0) && (sturm_i.size() > 1)) {
            sturm_i.pop_back();
        }
        for (int j = 0; j < sturm_i.size(); ++j) {
            sturm_i[j] = -sturm_i[j];
        }
        sturm_seq.push_back(sturm_i);
    }
    return sturm_seq;
}

double evaluate(double x, vector<double> coeffs) {
    double result = 0.0;
    for (int i = 0; i < coeffs.size(); i++) {
        result += coeffs[i] * pow(x, (coeffs.size() - i - 1));
    }
    return result;
}

int sign_changes(vector<double> sequence) {
    int changes = 0;
    int sign = 0; // the sign of the previous term
    for (int i = 0; i < sequence.size(); i++) {
        if (sign * sequence[i] < 0.0) { // if the sign has changed
            changes++;
            sign = sign * -1; // flip the sign
        } else if (sequence[i] != 0.0) { // if the sign hasn't changed, but the term is nonzero
            sign = sequence[i] < 0.0 ? -1 : 1;
        }
    }
    return changes;
}

int localisation(vector<vector<double>> &polinoms, double a, double b) {
    vector<double> value_1(polinoms.size());
    vector<double> value_2(polinoms.size());
    int num_roots = 0;
    int inverse_1 = 0, inverse_2 = 0;
    for (int i = 0; i < polinoms.size(); ++i) {
        value_2[i] = evaluate(b, polinoms[i]);
        value_1[i] = evaluate(a, polinoms[i]);
    }
    for (int i = 0; i < polinoms.size() - 1; ++i) {
        if (value_1[i] * value_1[i + 1] < 0) {
            inverse_1++;
        }
        if (value_2[i] * value_2[i + 1] < 0) {
            inverse_2++;
        }
    }
    num_roots = abs(inverse_1 - inverse_2);
    return num_roots;
}

void solve(vector<double> coeffs) {
    double epsilon = 0.01;
    const int max_iterations = 1000000000;

    vector<vector<double>> a;

    cout << "The polynomial you entered is: ";
    vector<double> location;
    cout << "The roots of the polynomial P(x) = ";
    for (int i = 0; i < coeffs.size(); i++) {
        if (coeffs[i] > 0.0) {
            cout << "+ " << coeffs[i] << "x^" << (coeffs.size() - i - 1)
                 << " ";
        } else if (coeffs[i] < 0.0) {
            cout << "- " << -coeffs[i] << "x^" << (coeffs.size() - i - 1) << " ";
        }
    }
    cout << endl;
    double left = -10, right = 10, mid, stop;
    mid = (right - left) / 2;
    vector<vector<double>> polinoms = sturm_sequence(coeffs);
    int num = localisation(polinoms, left, right);
    while (num == 0) {
        left *= 10;
        right *= 10;
        num = localisation(polinoms, left, right);
    }
    stop = right;
    location.push_back(left);
    while (num > 0) {
        int num_roots = localisation(polinoms, left, right);
        if (right > stop) {
            break;
        }
        if (num_roots == 1) {
            location.push_back(right);
            num--;
            mid = (right - left) / 2;
            left = right;
            right += mid;
        } else if (num_roots == 0) {
            right += mid + stop * epsilon;
            mid = (right - left) / 2;
        } else {
            right -= mid + stop * epsilon;
            mid = mid / 2;
        }
    }
    for (int i = 1; i < location.size(); ++i) {
        newton(coeffs, location[i], epsilon, max_iterations);
    }
    for (double r: roots) cout << "Root found: " << r << endl;
}


int main() {
    // Матрица для тестирования

    vector<vector<double>> matrix{{1, 1, 0, 0, 0, 0},
                                  {-2, 3, 0, 0, 0, 0},
                                  {0, 0, 0, 1, 0, 0},
                                  {0, 0, 0, 0, 1, 0},
                                  {0, 0, 0, 0, 0, 1},
                                  {0, 0, 0, 0, -5, 4}};
    printMatrix(matrixPower(matrix, 2));
    vector<double> a = leverrier(matrix);


    //for (int i = 0; i < a.size(); i++) {
    //    cout << a[i] << " ";
    //}

    //cout << endl;
    //solve(a);
    //roots.clear();
    /*
    matrix = {{2, 1,   1},
             {1, 2.5, 1},
             {1, 1,   3}};
    a = leverrier(matrix);

    for (int i = 0; i < a.size(); i++) {
        cout << a[i] << " ";
    }
    cout << endl;
    solve(a);
    */
    /*
    vector<vector<double>> matrix = random_matrix(100);
    //printMatrix(matrix); cout << endl;
    int mxn = 20;
    vector<double> ld(mxn);
    vector<double> fd(mxn);

    for (int j = 1; j < mxn; j++) {
        for (int i = 0; i < 1000; i++) {
            ld[j] = 0;
            fd[j] = 0;
            vector<vector<double>> matrix = random_matrix(j);

            auto start_time = std::chrono::high_resolution_clock::now();
            vector<double> a = leverrier(matrix);
            auto end_time = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
            ld[j] = ld[j] + duration;

            start_time = std::chrono::high_resolution_clock::now();
            vector<double> b = fadeef(matrix);
            end_time = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
            fd[j] = fd[j] + duration;


        }
        ld[j] = ld[j] / 100;
        fd[j] = fd[j] / 100;
    }
    for (int i = 0; i < ld.size(); i++) {
        cout << ld[i] << " ";
    }
    cout << endl;

    for (int i = 0; i < fd.size(); i++) {
        cout << fd[i] << " ";
    }
    cout << endl;

*/
    return 0;
}

