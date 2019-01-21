from setuptools import setup, find_packages
from codecs import open
import os

if sys.version < "3.6":
    sys.exit(
        "Error: Python 3.6 or greater required for Eniric (using {})".format(
            sys.version
        )
    )

base_dir = os.path.join(os.path.abspath(os.path.dirname(__file__)))

with open('requirements/requirements.txt') as f:
    requirements = f.read().splitlines()

with open(os.path.join(base_dir, "README.md")) as f:
    long_description = f.read()

setup(
    name='baraffe_tables',
    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version="0.2",

    description='Access the Baraffe evolutionary tables.',
    long_description=long_description,
    long_description_content_type="text/markdown",
    # The project's main homepage.
    url='https://github.com/jason-neal/baraffe_tables',
    author='Jason Neal',
    author_email='jason.neal@astro.up.pt',
    license='MIT Licence',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Topic :: Scientific/Engineering :: Physics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6+',
        'Natural Language :: English',
    ],

    # What does your project relate to?
    keywords=['astronomy'],

    packages=find_packages(exclude=['contrib', 'docs', 'test']),

    install_requires=requirements,
    # install_requires=[],
    setup_requires=['pytest-runner'],
    tests_require=['pytest', "hypothesis"],
    # $ pip install -e .[dev,test]
    extras_require={
        'dev': ['check-manifest'],
        'test': ['coverage', 'pytest', 'pytest-cov', 'python-coveralls', 'hypothesis'],
        'docs': ['sphinx >= 1.4', 'sphinx_rtd_theme', 'pyastronomy']
    },

    # If there are data files included in your packages that need to be
    # installed, specify them here.  If using Python 2.6 or less, then these
    # have to be included in MANIFEST.in as well.
    # package_data={"spectrum_overload": ["data/*.fits"]},
    package_data={"baraffe_tables": ["data/Baraffe20*/*.dat"]},

    data_files=[],

    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # pip to create the appropriate form of executable for the target platform.
    entry_points={
        #    'console_scripts': [
        #        'sample=sample:main',
        # ],
    },
    scripts=["baraffe_tables/mass_to_flux_ratio.py",
             "baraffe_tables/flux_ratio_to_mass.py",
             "baraffe_tables/query_baraffe.py",
             "baraffe_tables/teff2mass.py"],
)
