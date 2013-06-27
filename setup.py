from setuptools import setup, Extension
#from distutils.core import setup
#from distutils.extension import Extension

# try:
#     from Cython.Distutils import build_ext
#     cmdclass = { 'build_ext': build_ext }
# except ImportError:
#     cmdclass = {}
# print cmdclass
# import numpy
# np_include = numpy.get_include()

setup(
    name = 'basekit',
    version = '0.0.2',
    #cmdclass = cmdclass,
    ext_modules = [ 
        Extension(
            "basekit/utils/cgeom", 
            sources = [ "basekit/utils/cgeom.c" ]
        )
    ],
    url = 'www.weirdbyte.de',
    author = 'Alexander Rose',
    author_email = 'alexander.rose@weirdbyte.de',
    packages = [ 'basekit', 'basekit.utils' ],
    install_requires = [ 'numpy', 'matplotlib' ]
)