import os
import re
import sys
import platform
import subprocess
import versioneer

from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext
from distutils.version import LooseVersion


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError("CMake must be installed to build the following extensions: " + ", ".join(
                e.name for e in self.extensions))

        if platform.system() == "Windows":
            cmake_version = LooseVersion(re.search(r'version\s*([\d.]+)', out.decode()).group(1))
            if cmake_version < '3.1.0':
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        global cmake_args
        bypass_install = cmake_args.pop('-DBYPASS_INSTALL')

        internal_cmake_args = ['-DPYTHON_EXECUTABLE=' + sys.executable]
        internal_cmake_args += [k + "=" + v for k, v in cmake_args.items() if v]

        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]

        if platform.system() == "Windows":
            if sys.maxsize > 2**32:
                cmake_args += ['-A', 'x64']
            build_args += ['--', '/m']
        else:
            internal_cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
            build_args += ['--', '-j2']

        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(env.get('CXXFLAGS', ''), self.distribution.get_version())
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(['cmake', ext.sourcedir] + internal_cmake_args, cwd=self.build_temp, env=env)
        subprocess.check_call(['cmake', '--build', '.'] + build_args, cwd=self.build_temp)
        if not bypass_install:
            subprocess.check_call(['cmake', '--build', '.', '--target', 'install'], cwd=self.build_temp)


if __name__ == "__main__":

    # Valid CMake args
    valid_args = {
        '-DCMAKE_BUILD_TYPE': 'Release',
        '-DENABLE_XHOST': 'ON',
        '-DMAX_AM': '6',
        '-DCMAKE_C_FLAGS': False,
        '-DCMAKE_C_COMPILER': False,
        '-DCMAKE_PREFIX_PATH': False,
        '-DNATIVE_PYTHON_INSTALL_WITH_LIB': 'OFF',
        '-DBYPASS_INSTALL': False,
    }
    invalid_args = {
        '-DBUILD_SHARED_LIBS': 'ON',
        '-DENABLE_GENERIC': 'OFF',
        '-DBUILD_FPIC': 'ON',
        '-DINSTALL_PYMOD': 'ON',
        '-DNATIVE_PYTHON_INSTALL': 'ON'
    }
    cmake_args = valid_args.copy()
    cmake_args.update(invalid_args)

    # Parse out CMake args
    setup_args = []
    for arg in sys.argv:
        if "-D" not in arg:
            setup_args.append(arg)
            continue

        split_arg = [x.strip() for x in arg.split('=')]
        if len(split_arg) != 2:
            raise KeyError("CMake argument %s not understood." % arg)
        key, value = split_arg

        if key not in cmake_args:
            raise KeyError("CMake argument %s not understood." % arg)

        if key in invalid_args:
            raise KeyError("CMake argument %s cannot be changed with Python builds." % key)

        cmake_args[key] = value

    sys.argv = setup_args

    # Build full cmdclass
    cmdclass = versioneer.get_cmdclass()
    cmdclass["build_ext"] = CMakeBuild

    setup(
        name='gau2grid',
        version=versioneer.get_version(),
        description='Fast computation of a gaussian and its derivative on a grid.',
        author='Daniel G. A. Smith',
        author_email='dgasmith@icloud.com',
        url="https://github.com/dgasmith/gau2grid",
        license='BSD-3C',
        packages=find_packages(),
        include_package_data=True,
        ext_modules=[CMakeExtension('gau2grid.gg')],
        cmdclass=cmdclass,
        install_requires=[
            'numpy>=1.7',
        ],
        extras_require={
            'docs': [
                'sphinx==1.2.3',  # autodoc was broken in 1.3.1
                'sphinxcontrib-napoleon',
                'sphinx_rtd_theme',
                'numpydoc',
            ],
            'tests': [
                'pytest',
                'pytest-cov',
            ],
        },
        tests_require=[
            'pytest',
            'pytest-cov',
        ],
        classifiers=[
            'Development Status :: 4 - Beta',
            'Intended Audience :: Science/Research',
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3',
        ],
        zip_safe=False,
    )
