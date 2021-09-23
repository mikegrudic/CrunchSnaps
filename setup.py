import setuptools
from setuptools import setup

#with open("README.md", "r", encoding="utf-8") as fh:
#    long_description = fh.read()

setup(name='CrunchSnaps',
      version='0.1',
      description='Package for running efficient passes on GIZMO simulation data',
#      long_description=long_description,
      url='http://github.com/mikegrudic/CrunchSnaps',
      author='Mike Grudić',
      author_email='mike.grudich@gmail.com',
      license='MIT',
#      packages=['meshoid'],
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
      zip_safe=False)
