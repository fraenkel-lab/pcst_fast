#! /bin/bash
git subtree pull --prefix external/pybind11 https://github.com/pybind/pybind11.git stable --squash
git subtree pull --prefix external/googletest https://github.com/google/googletest.git master --squash
