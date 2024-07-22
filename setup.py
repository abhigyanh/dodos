from setuptools import setup, find_packages

setup(
    name='dodos',
    version='1.21',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'scipy',
        'matplotlib',
    ],
    entry_points={
        'console_scripts': [
            'dodos=src.launch:main',
        ],
    },
)