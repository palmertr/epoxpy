from __future__ import print_function


from setuptools import setup, find_packages

setup(name='epoxpy',
        version='0.1',
        description='Simulation wrappers for epoxy simulations using HOOMD',
        url='https://bitbucket.org/cmelab/epoxy_sim',
        author='cmelab',
        author_email='stephenthomas1@boisestate.edu',
        license='MIT',
        packages=find_packages(),
        package_dir={'epoxpy':'epoxpy'},
        install_requires=[
            'numpy',
            'matplotlib', 'mbuild',
        ],
        zip_safe=False)
