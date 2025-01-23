#include <pybind11/pybind11.h>  //includes also Python.h -> it has to be the first include!
                                // Indeed, Python.h defines some preprocessor variables
                                // that may affect the behavior of the standard headers.
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <pybind11/embed.h>  //python interpreter
#include <pybind11/stl.h>  // type conversion
#include <Eigen/Sparse>
#include <complex>
#include<array>
#include<vector>


/* Pybind11 is a library that exposes C++ type in Python and vice versa in order to make bindings.
It is lighter than the Boost.Python project (that was designed for the same purpose),
because it only matches Cpp11 compilers and above.*/


namespace py = pybind11;

const std::size_t k = 1; 
std::array<std::complex<double>, k> diagonalize(const Eigen::SparseMatrix<std::complex<double>>& H) {
    py::object diagonalize_py = py::module::import("diag").attr("diagonalize");
    auto result = diagonalize_py(H);
    return result.cast<std::array<std::complex<double>, k>>();
}
 ;


