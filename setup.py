from setuptools import setup, Extension

try:
    from Cython.Distutils import build_ext
    cmdclass = { 'build_ext': build_ext }
except ImportError:
    cmdclass = {}

import numpy
np_include = numpy.get_include()

setup(
    name = 'basekit',
    version = '0.0.2',
    cmdclass = cmdclass,
    ext_modules = [ 
        Extension(
            "basekit/utils/cgeom", 
            sources = [ "basekit/utils/cgeom.c" ],
            include_dirs = [ np_include ]
        ),
        Extension(
            "basekit/utils/align/calign", 
            sources = [ "basekit/utils/align/calign.c" ],
            include_dirs = [ np_include ]
        )
    ],
    url = 'www.weirdbyte.de',
    author = 'Alexander Rose',
    author_email = 'alexander.rose@weirdbyte.de',
    packages = [ 
        'basekit', 
        'basekit.utils',
        'basekit.utils.align'
    ],
    install_requires = [ 'numpy', 'matplotlib', 'poster' ]
)