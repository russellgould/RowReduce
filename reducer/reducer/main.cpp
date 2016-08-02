#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <iomanip>
#include <iostream>

using namespace boost::numeric::ublas;

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

  for (unsigned i = 0; i < m.size1(); i++) {
    for (unsigned j = 0; j < m.size2(); j++) {
      std::cout << "[" << m(i, j) << "] ";
    }
    std::cout << std::endl;
  }

  return 0;
}