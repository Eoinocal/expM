# expM

### Dubious Ways to Implement the Exponential of a Matrix

This repository is (will be) a small C++ library for computing the exponential of a matrix. It is inspired by the paper [Nineteen Dubious Ways to Compute the Exponential of a Matrix, Twenty-Five Years Later]() by Cleve Moler and Charles Van Loan. Note that I have no affiliation with either author. The abstract of their paper neatly sums up why there are so many *dubious* approaches to compute the exponential:

>In principle, the exponential of a matrix could be computed in many ways. Methods involving
approximation theory, differential equations, the matrix eigenvalues, and the matrix
characteristic polynomial have been proposed. In practice, consideration of computational
stability and efficiency indicates that some of the methods are preferable to others, but
that none are completely satisfactory.

## Dependencies

The library depends on [Boost](<http://www.boost.org>), mainly the [uBlas](<http://www.boost.org/doc/libs/release/libs/numeric/ublas/>) library, but also leverages other Boost libraries where appropriate. I aim the limit the dependencies to Boost for the CPU implementations but I also hope to provide GPGPU implementations too. It is likely the GPGPU implementations will rely on [VexCL](<https://github.com/ddemidov/vexcl>), [ViennaCL](<http://viennacl.sourceforge.net/index.html>) and possibly others as the library grows. There is a nice set of linear algebra and GPGPU libraries which all inter-operate, not availing of them would be silly.

## Goals

### Backends
The main GPGPU backend I plan to support is OpenCL, but seeing as many of the libraries I am building support multiple backends, for example CUDA, I would like to maintain that feature if at all possible.

### Code Style
The initial implementations will concentrate on clean and clear code, ideally aiming almost one-for-one rewrite of the mathematically presented algorithms in C++ code. However this *clean* implementations will break that restriction to avoid obvious inefficiencies. For example excessive copying, temporary variables and redundant recalculations will be avoided. This will inevitably mean that some algorithms will look, *superficially*, different from the theoretically presented ones, hoever they will perform the exact same calculation.

### Portability
Initial development is taking place on Windows with Visual Studio 2015 Community. But care is being taken to not use any platform specific code. I will test on additional platforms and compilers when for of the implementation is more fleshed out.

## License

The library is distributed under the [Boost Software License](<http://www.boost.org/users/license.html>).
