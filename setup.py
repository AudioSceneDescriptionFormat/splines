from setuptools import setup

# "import" __version__
__version__ = 'unknown'
for line in open('src/splines.py'):
    if line.startswith('__version__'):
        exec(line)
        break

setup(
    name='splines',
    py_modules=['splines'],
    version=__version__,
    author='Matthias Geier',
    author_email='Matthias.Geier@gmail.com',
    description='Piecewise Polynomial Curves',
    long_description=open('README.rst').read(),
    package_dir={'': 'src'},
    install_requires=['NumPy'],
    license='MIT',
    keywords=''.split(),
    project_urls={
        'Documentation': 'http://splines.readthedocs.io/',
        'Source Code': 'https://github.com/AudioSceneDescriptionFormat/splines/',
        'Bug Tracker': 'https://github.com/AudioSceneDescriptionFormat/splines/issues/',
    },
    platforms='any',
    classifiers=[
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
    ],
    zip_safe=True,
)
