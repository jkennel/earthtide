# earthtide 0.1.2

* default to single threaded
* fix valgrind issue to get back on CRAN

# earthtide 0.1.1

* removed BH dependency

# earthtide 0.1.0

* switched from RcppParallel to RcppEigen. Eigen performs well regardless 
of the BLAS installed.
* switched to RcppThread from RcppParallel to simplify code base.
* Reordered C++ loops to improve performance.
* astro_update parameter is no longer used - they are updated for each time

# earthtide 0.0.13
* updated the documentation for latitude, longitude and elevation (thanks kohlerjl)

# earthtide 0.0.12

* updated the get_iers download locations (thanks kohlerjl)

# earthtide 0.0.9

* remove get_iers test as it can cause failure if no internet accesss, or the server is down.

# earthtide 0.0.8

* Added horizontal displacement calculation
* Added horizontal strain calculation
* Add scaling option for analyze function

# earthtide 0.0.7*

* Added a `NEWS.md` file to track changes to the package.
* Matrix output method
