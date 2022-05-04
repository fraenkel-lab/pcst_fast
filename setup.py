from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import sys
import setuptools

# include external/pybind11/include/pybind11/*.h


class get_pybind_include(object):
    """Helper class to determine the pybind11 include path

    The purpose of this class is to postpone importing pybind11
    until it is actually installed, so that the ``get_include()``
    method can be invoked. """

    def __init__(self, user=False):
        self.user = user

    def __str__(self):
        import pybind11
        return pybind11.get_include(self.user)


ext_modules = [
    Extension(
        'pcst_fast',
        sources=['src/pcst_fast_pybind.cc', 'src/pcst_fast.cc'],
        include_dirs=[
            # Path to pybind11 headers
            get_pybind_include(),
            get_pybind_include(user=True)
        ],
        language='c++'
    ),
]

# As of Python 3.6, CCompiler has a `has_flag` method.
# cf http://bugs.python.org/issue26689
def has_flag(compiler, flagname):
    """Return a boolean indicating whether a flag name is supported on
    the specified compiler.
    """
    import tempfile
    with tempfile.NamedTemporaryFile('w', suffix='.cpp') as f:
        f.write('int main (int argc, char **argv) { return 0; }')
        try:
            compiler.compile([f.name], extra_postargs=[flagname])
        except setuptools.distutils.errors.CompileError:
            return False
    return True


def cpp_flag(compiler):

    if has_flag(compiler, '-std=c++11'): return '-std=c++11'
    else: raise RuntimeError('Unsupported compiler -- at least C++11 support is needed!')


class BuildExt(build_ext):
    """A custom build extension for adding compiler-specific options."""
    c_opts = {
        'msvc': ['/EHsc'],
        'unix': [],
    }

    if sys.platform == 'darwin':
        c_opts['unix'] += ['-stdlib=libc++', '-mmacosx-version-min=10.7']

    def build_extensions(self):
        ct = self.compiler.compiler_type
        opts = self.c_opts.get(ct, [])
        if ct == 'unix':
            opts.append('-DVERSION_INFO="%s"' % self.distribution.get_version())
            opts.append(cpp_flag(self.compiler))
            if has_flag(self.compiler, '-fvisibility=hidden'):
                opts.append('-fvisibility=hidden')
        elif ct == 'msvc':
            opts.append('/DVERSION_INFO=\\"%s\\"' % self.distribution.get_version())

        for ext in self.extensions:
            ext.extra_compile_args = opts

        build_ext.build_extensions(self)


# Read the contents of the README file
with open("README.md", encoding="utf-8") as f:
    README = f.read()

setup(
    name='pcst_fast',
    packages=['pcst_fast'],
    package_dir={'pcst_fast': 'src'},
    use_scm_version=True,
    url='https://github.com/fraenkel-lab/pcst_fast',
    license='GNU General Public License',
    author='ludwigschmidt',
    author_email='alex@lenail.org',
    description='A fast implementation of the Goemans-Williamson scheme for the prize-collecting Steiner tree / forest problem.',
    long_description=README,
    long_description_content_type="text/markdown",
    ext_modules=ext_modules,
    install_requires=['pybind11>=2.1.0'],
    setup_requires=['pybind11>=2.1.0', 'setuptools_scm'],
    cmdclass={'build_ext': BuildExt},
    zip_safe=False,
)
