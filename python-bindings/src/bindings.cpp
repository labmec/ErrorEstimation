//
// Created by gustavo on 16/11/2019.
//

#include <pybind11/pybind11.h>
#include "TPZHDivErrorEstimatorH1.h"

namespace py = pybind11;

PYBIND11_MODULE(errorestimation, m) {
    m.doc() = R"pbdoc(
        -------------------------
        Python bindings for NeoPZ
        -------------------------
    )pbdoc"; // TODO: CHANGE THIS

    py::class_<TPZMultiphysicsCompMesh>(m, "TPZMultiphysicsCompMesh")
            .def(py::init<>())
            .def("AutoBuild", &TPZMultiphysicsCompMesh::AutoBuild)
            .def("__repr__",
                 [](const TPZMultiphysicsCompMesh &obj) {
                     return "It works, go home.";
                 }
            );
}