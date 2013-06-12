from setuptools import setup

setup(
    name = 'basekit',
    version = '0.0.2',
    url = 'www.weirdbyte.de',
    author = 'Alexander Rose',
    author_email = 'alexander.rose@weirdbyte.de',
    packages = ['basekit', 'basekit.utils'],
    install_requires = ['numpy', 'matplotlib']
)