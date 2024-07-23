from setuptools import setup, find_packages
from src.version import project_version

setup(
    name='dodos',
    version=str(project_version),
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
