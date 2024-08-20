from setuptools import setup, find_packages
import codecs
import os

here = os.path.abspath(os.path.dirname(__file__))

with codecs.open(os.path.join(here, "README.md"), encoding="utf-8") as fh:
    long_description = fh.read()


setup(
    name='nuctools',
    version='{{VERSION_PLACEHOLDER}}', # updated in CI/CD
    description="nuclear data reduction tools",
    author="Jesse Brown",
    author_email="brownjm1968@gmail.com",
    packages=['nuctools'],
    include_package_data=True,
    license='BSD 3-clause',
    url='https://github.com/brownjm1968/nuctools',
    long_description_content_type="text/markdown",
    long_description=long_description,
    install_requires=['numpy','pandas','scipy','h5py','pyyaml',
                      'matplotlib','pytest'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: Unix",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows"
    ]
)

