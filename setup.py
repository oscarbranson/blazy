from setuptools import setup, find_packages

try:
  from blazy import __version__
except:
  __version__ = "version_missing"

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='blazy',
      version=__version__,
      description='Tools for solution speciation calculations with PHREEQC',
      long_description=long_description,
      long_description_content_type="text/markdown",
      url='https://github.com/oscarbranson/blazy',
      author='Oscar Branson',
      author_email='ob266@cam.ac.uk',
      license='MIT',
      packages=find_packages(),
      classifiers=['Development Status :: 4 - Beta',
                   'Intended Audience :: Science/Research',
                   'Topic :: Scientific/Engineering',
                   'Programming Language :: Python :: 3',
                   ],
      python_requires='>3.6',
      install_requires=['numpy',
                        'pandas',
                        'matplotlib',
                        'uncertainties',
                        'scikit-learn',
                        'scipy',
                        'Ipython',
                        'configparser',
                        'tqdm',
                        'phreeqpy'
                        ],
      package_data={
        'latools': ['resources/*',
                    'resources/database/*'],
      },
      zip_safe=True)
