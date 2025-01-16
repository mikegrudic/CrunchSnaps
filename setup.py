import setuptools, os
from setuptools import setup

# with open("README.md", "r", encoding="utf-8") as fh:
#    long_description = fh.read()

thelibFolder = os.path.dirname(os.path.realpath(__file__))
requirementPath = thelibFolder + "/requirements.txt"
install_requires = []
if os.path.isfile(requirementPath):
    with open(requirementPath) as f:
        install_requires = f.read().splitlines()

setup(
    name="CrunchSnaps",
    version="0.1",
    description="Package for running efficient passes on GIZMO simulation data",
    url="http://github.com/mikegrudic/CrunchSnaps",
    author="Mike Grudić",
    author_email="mike.grudich@gmail.com",
    license="MIT",
    project_urls={
        "Bug Tracker": "https://github.com/mikegrudic/CrunchSnaps",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.6",
    zip_safe=False,
    install_requires=install_requires,
    include_package_data=True
)
