from distutils.core import setup

setup(
    name='Brot',
    version='0.1',
    packages=['brot'],
    install_requires=['GitPython', 'numpy<2', 'scipy'],
    license='GNU GENERAL PUBLIC LICENSE, Version 3',
    long_description=open('README.md').read(),
)