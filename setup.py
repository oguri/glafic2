import os
import sys
import subprocess

from setuptools import setup, find_packages
from setuptools.command.build_py import build_py

# here = os.path.abspath(os.path.dirname(__file__))

def build_extensions():
    subprocess.run(
        "make clean; make python",
        shell = True,
        check = True,
    )


def make_cmd_class(base):
    class CmdClass(base):
        def run(self):
            build_extensions()
            super().run()

    return CmdClass

setup(
    name = 'glafic',
    version = '2.1.9',
    author = 'Masamune Oguri',
    description = 'Library for gravitational lensing analyses',
    url = 'https://github.com/oguri/glafic2',
    packages = find_packages(where = 'python'),
    package_dir = {'': 'python'},
    python_requires = '>=3',
    package_data = {'glafic': ['*.so']},
    cmdclass = {
       'build_py': make_cmd_class(build_py)
    },
)
