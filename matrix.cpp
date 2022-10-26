#include "matrix.h"
#include <algorithm>
#include <stdexcept>

Matrix::Matrix() {

}

Matrix::~Matrix() {

}

void Matrix::addRow(std::vector<double> newRow) {
    if (newRow.size() > cols) {
        for (auto& r: matrix) {
            while (r.size() < newRow.size()) {
                r.push_back(0);
            }
        }
        cols = newRow.size();
    }
    else if (newRow.size() < cols) {
        while (newRow.size() < cols) {
            newRow.push_back(0);
        }
    }
    matrix.push_back(newRow);
    ++rows;
}

std::vector<double> Matrix::gaussElim(const std::vector<double>& B) {
    std::vector<double> C = B;
    if (rows != cols) throw std::invalid_argument("eliminate called on invalid matrix");
    int n = rows;
    for (int i = 0; i < n; i++) {
        if (matrix[i][i] == 0) {
            int c = 1;
            while ((i+c) < n && matrix[i+c][i] == 0) {
                c++;
            }
            if ((i+c) == n) throw std::invalid_argument("special case encountered");
            int j = i;
            for (int k = 0; k < n; k++) {
                std::swap(matrix[j][k], matrix[j+c][k]);
            }
            std::swap(C[j], C[j+c]);
        }
        for (int j = 0; j < n; j++) {
            if (i != j) {
                double ratio = matrix[j][i]/matrix[i][i];
                for (int k = 0; k < n; k++) {
                    matrix[j][k] = matrix[j][k] - ratio*matrix[i][k];
                }
                C[j] = C[j] - ratio*C[i];
            }
        }
    }
    std::vector<double> solutions;
    for (int i = 0; i < n; i++) {
        solutions.push_back(C[i]/matrix[i][i]);
    }
    return solutions;
}

Matrix transpose(const Matrix& A) { // Only implemented for square matrices
    Matrix B(A);
    if (A.rows != A.cols) throw std::invalid_argument("transpose of non-square matrix");
    for (unsigned int i = 0; i < B.rows; i++) {
        for (unsigned int j = i + 1; j < B.cols; j++) {
            double temp = B.matrix[i][j];
            B.matrix[i][j] = B.matrix[j][i];
            B.matrix[j][i] = temp;
        }
    }
    return B;
}

std::vector<double> multiply(const Matrix& A, const std::vector<double> vec) {
    std::vector<double> solutions;
    if (A.cols != vec.size()) throw std::invalid_argument("invalid matrix multiplication");
    for (unsigned int i = 0; i < A.rows; i++) {
        double temp = 0;
        for (unsigned int j = 0; j < A.cols; j++) {
            temp += A.matrix[i][j] * vec[j];
        }
        solutions.push_back(temp);
    }
    return solutions;
}

Matrix multiply(const Matrix& A, const Matrix& B) {
    Matrix C;
    if (A.cols != B.rows) throw std::invalid_argument("invalid matrix multiplication");
    for (unsigned int i = 0; i < A.rows; i++) {
        std::vector<double> row;
        for (unsigned int j = 0; j < B.cols; j++) {
            double temp = 0;
            for (int k = 0; k < A.rows; k++) {
                temp += A.matrix[i][k] * B.matrix[k][j];
            }
            row.push_back(temp);
        }
        C.addRow(row);
    }
    C.rows = A.rows;
    C.cols = B.cols;
    return C;
}

std::vector<double>& Matrix::operator[](int index) {
    return matrix[index];
}

std::vector<double> operator*(const Matrix& A, const std::vector<double>& vec) {
    std::vector<double> solution = multiply(A, vec);
    return solution;
}

Matrix operator*(const Matrix& A, const Matrix& B) {
    Matrix C = multiply(A, B);
    return C;
}

std::ostream& operator<<(std::ostream& os, const Matrix& A) {
    os << A.rows << "x" << A.cols << std::endl;
    for (unsigned int i = 0; i < A.rows; i++) {
        os << "< ";
        for (unsigned int j = 0; j < A.cols; j++) {
            os << A.matrix[i][j] << " ";
        }
        os << ">" << std::endl;
    }
    return os;
}

std::vector<double> operator-(const std::vector<double>& A, const std::vector<double>& B) {
    std::vector<double> solution;
    if (A.size() != B.size()) return solution;
    for (unsigned int i = 0; i < A.size(); i++) {
        solution.push_back(A[i] - B[i]);
    }
    return solution;
}