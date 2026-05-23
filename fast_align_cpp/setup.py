from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import find_packages, setup


ext_modules = [
    Pybind11Extension(
        "scan_splash_align._core",
        ["src/scan_splash_align/_align.cpp"],
        cxx_std=17,
    ),
]


setup(
    name="scan-splash-align",
    version="0.1.0",
    description="Fast C++ alignment helpers for SCAN_SPLASH",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
)
