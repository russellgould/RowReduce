#include "matrix.h"

matrix::matrix(size_t rows, size_t cols)
    : mRows(rows), mCols(cols), mData(rows * cols) {}

matrix::~matrix() {}
