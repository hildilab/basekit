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
    install_requires = [ 'numpy', 'matplotlib', 'poster', 'fastcluster' ],
    scripts=[
        'scripts/apbs.py',
        'scripts/capture.py',
        'scripts/dowser.py',
        'scripts/dssp.py',
        'scripts/hbexplore.py',
        'scripts/jmol.py',
        'scripts/linker.py',
        'scripts/mapman.py',
        'scripts/moderna_tools.py',
        'scripts/motif.py',
        'scripts/mppd.py',
        'scripts/mpstruc.py',
        'scripts/msa.py',
        'scripts/msms.py',
        'scripts/opm.py',
        'scripts/pdb.py',
        'scripts/project.py',
        'scripts/solvate.py',
        'scripts/spider.py',
        'scripts/sstruc.py',
        'scripts/tmdet.py',
        'scripts/voronoia.py',
        'scripts/voronoia_pipeline.py',
        'scripts/wine_watcher.py',
    ]
)