#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <cmath>
#include <iomanip>
#include <iostream>

using namespace boost::numeric::ublas;

// elementary row-operations

// swaps row1 and row2
void rowInterchange(matrix<double> &m, int row1, int row2) {
  double *oldRow = nullptr;
  oldRow = new double[m.size2()];
  for (unsigned i = 0; i < m.size2(); i++) {
    oldRow[i] = m(row2, i);
    m.insert_element(row2, i, m(row1, i));
  }

  for (unsigned i = 0; i < m.size2(); i++) {
    m.insert_element(row1, i, oldRow[i]);
  }

  delete[] oldRow;
  oldRow = nullptr;
}

// adds "factor" multiple of row 1 to row 2
void addMultiple(matrix<double> &m, int row1, double factor, int row2) {
  for (unsigned i = 0; i < m.size2(); i++) {
    m.insert_element(row2, i, (m(row1, i) * factor) + m(row2, i));
  }
}

// scales an entire row by "factor"
void scaleRow(matrix<double> &m, int row, double factor) {
  for (unsigned i = 0; i < m.size2(); i++) {
    m.insert_element(row, i, m(row, i) * factor);
  }
}

// finds left-most nonzero column in matrix
int getLeftMostCol(const matrix<double> &m, bool &found, int col) {
  for (unsigned i = 0; i < m.size1(); i++) {
    if (m(i, col) != 0) {
      found = true;
    }
  }

  if (!found) {
    col = getLeftMostCol(m, found, col + 1);
  }

  return col;
}

// gets value for row containing entry with largest absolute value in the passed
// column
int getRowLarEnInCol(matrix<double> &m, int col) {
  int row = 0;

  for (unsigned i = 0; i < m.size1(); i++) {
    if (abs(m(i, col)) > row) {
      row = i;
    }
  }

  return row;
}

void printMatrix(const matrix<double> &m) {
  std::cout << std::endl;
  for (unsigned i = 0; i < m.size1(); i++) {
    for (unsigned j = 0; j < m.size2(); j++) {
      std::cout << "[" << m(i, j) << "] ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

int main(int argc, const char *argv[]) {
  int nRows(0), nCols(0);

  std::cout << "This is a simple row-reduction program." << std::endl;
  std::cout << "To enter a matrix, first enter the number of rows: ";
  std::cin >> nRows;
  std::cout << "Now enter the number of columns: ";
  std::cin >> nCols;

  matrix<double> m(nRows, nCols);

  std::cout << "Now enter numbers, one at a time, in row-major order: "
            << std::endl;
  double input = 0.0;

  for (unsigned i = 0; i < m.size1(); i++) {
    for (unsigned j = 0; j < m.size2(); j++) {
      std::cin >> input;
      m.insert_element(i, j, input);
    }
  }

  std::cout << std::endl
            << "The matrix you entered is: " << std::endl
            << std::endl;

  printMatrix(m);

  // std::cout << std::endl << "Now to row-reduce..." << std::endl;

  // bool found = false;
  // int mLNoZ = getLeftMostCol(m, found, 0);

  int row = 0;
  std::cout << "What row would you like to scale: ";
  std::cin >> row;
  std::cout << "By what factor?" << std::endl;
  double factor = 0;
  std::cin >> factor;
  std::cout << std::endl;

  scaleRow(m, row, factor);

  printMatrix(m);

  return 0;
}