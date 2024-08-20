from setuptools import setup, find_packages

# set the __version__
exec(open('nuctools/_version.py').read())

setup(
    name='nuctools',
    version='{{VERSION_PLACEHOLDER}}',
    description="nuclear data reduction tools",
    author="Jesse Brown",
    author_email="brownjm1968@gmail.com",
    packages=['nuctools'],
    include_package_data=True,
    license='BSD 3-clause',
    url='https://github.com/brownjm1968/nuctools',
    long_description_content_type="text/markdown",
    long_description=open('README.md').read(),
    install_requires=['numpy','pandas','scipy','h5py','pyyaml',
                      'matplotlib','pytest']
)

