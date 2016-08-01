#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <cmath>

class matrix
{
private:
  size_t mRows;
  size_t mCols;
  std::vector<double> mData;

public:
  // constructors and destructor
  matrix(size_t rows, size_t cols);
  ~matrix();
  
  // accessor functions


  // mutator functions

};

#endif