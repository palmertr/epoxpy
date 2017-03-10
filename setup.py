from setuptools import setup

setup(name='epoxpy',
        version='0.1',
        description='Simulation wrappers for epoxy simulations using HOOMD',
        url='https://bitbucket.org/cmelab/epoxy_sim',
        author='cmelab',
        author_email='stephenthomas1@boisestate.edu',
        license='MIT',
        packages=['epoxpy'],
        install_requires=[
            'numpy',
            'matplotlib',
        ],
        zip_safe=False)
