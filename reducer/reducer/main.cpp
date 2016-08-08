#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>

using namespace boost::numeric::ublas;

// elementary row-operations
//
// swaps row1 and row2
void rowInterchange(matrix<double> &m, int row1, int row2) {
  double *oldRow(nullptr);
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

// adds "factor" multiple of row1 to row2
void addMultiple(matrix<double> &m, int row1, int row2, double factor) {
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

// other operations
//
// finds column to be used as pivot column
int getPivCol(const matrix<double> &m, bool &found, int startRow, int col) {
  if (col == m.size2()) { // right-most col is 0's
    col -= 1;
    found = true;
  } else {
    for (unsigned i = startRow; i < m.size1(); i++) {
      if (m(i, col) != 0) {
        found = true;
      }
    }
  }

  if (!found) { // some col is 0's, keep looking
    col = getPivCol(m, found, startRow, col + 1);
  }

  return col;
}

// gets value for row containing entry with largest absolute value in the passed
// column
int getPivRow(matrix<double> &m, int startRow, int col) {
  double max(0);
  int row(0);
  for (unsigned i = startRow; i < m.size1(); i++) {
    if (abs(m(i, col)) > max) {
      max = abs(m(i, col));
      row = i;
    }
  }

  if (row == 0) { // all remaining entries in col are 0, return current row
    row = startRow;
  }

  return row;
}

// finds pivot entry in row
int getPivEntry(matrix<double> &m, int row) {
  int pivot(0), i(0);
  bool found(false);
  do {
    if (m(row, i) != 0) {
      pivot = i;
      found = true;
    }
    i++;
  } while (!found && i < m.size2());

  if (!found && i == m.size2()) {
    pivot = -1;
  }

  return pivot;
}

// creates zeros in all entries under pivot in given column
void zeroCol(matrix<double> &m, int row, int col) {
  double pivot = m(row, col);
  for (unsigned i = row + 1; i < m.size1(); i++) {
    if (m(i, col) != 0) {
      addMultiple(m, row, i, -(m(i, col) / pivot));
    }
  }
}

// creates zeros in all entries above pivot in given column
void zeroColUp(matrix<double> &m, int row, int col) {
  for (unsigned i = row - 1; i >= 0 && i < m.size1(); i--) {
    if (m(i, col) != 0) {
      addMultiple(m, row, i, -m(i, col));
    }
  }
}

void printMatrix(const matrix<double> &m) {
  std::cout << std::endl;
  for (unsigned i = 0; i < m.size1(); i++) {
    for (unsigned j = 0; j < m.size2(); j++) {
      std::cout << "[" << m(i, j) << "] ";
    }
    std::cout << std::endl;
  }
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
  double input(0.0);

  for (unsigned i = 0; i < m.size1(); i++) {
    for (unsigned j = 0; j < m.size2(); j++) {
      std::cin >> input;
      m.insert_element(i, j, input);
    }
  }
  std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  std::cout << std::endl << "The matrix you entered is: " << std::endl;
  printMatrix(m);

  // create variables to store current and previous pivots
  bool found(false);
  size_t prevPivRow(0), prevPivCol(0), curPivRow(0), curPivCol(0);

  // find first pivot point to start the iterations correctly
  curPivCol = getPivCol(m, found, 0, 0);
  curPivRow = getPivRow(m, 0, curPivCol);

  // loop to fully reduce matrix to reduced echelon form
  for (unsigned i = 0; i < m.size1() && i < m.size2(); i++) {
    rowInterchange(m, curPivRow, i);
    curPivRow = i;
    std::cout << "after swap: " << std::endl;
    printMatrix(m);

    // zero out entries underneath pivot
    zeroCol(m, curPivRow, curPivCol);
    std::cout << "after zeroing down: " << std::endl;
    printMatrix(m);

    // scale entire row so that entry at pivot position is 1
    if (m(curPivRow, curPivCol) != 1 && m(curPivRow, curPivCol) != 0) {
      scaleRow(m, i, (1 / m(curPivRow, curPivCol)));
    }
    std::cout << "after scaling: " << std::endl;
    printMatrix(m);

    // zero out entries above pivot for reduced echelon form
    if (i > 0) {
      zeroColUp(m, curPivRow, curPivCol);
    }
    std::cout << "after zeroing up: " << std::endl;
    printMatrix(m);

    prevPivCol = curPivCol;
    prevPivRow = curPivRow;

    // get next pivot position
    found = false;
    curPivCol = getPivCol(m, found, i + 1, prevPivCol + 1);
    curPivRow = getPivRow(m, i + 1, curPivCol);
  }

  std::cout << std::endl;
  std::cout << "Your reduced matrix is: " << std::endl;
  printMatrix(m);
  std::cout << std::endl;

  std::cout << "Press ENTER or close the window to exit...";
  std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  return 0;
}
