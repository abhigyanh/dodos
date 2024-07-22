from setuptools import setup, find_packages

setup(
    name='dodos',
    version='1.21',
    packages=find_packages(),
    install_requires=[
        # Add your dependencies here
    ],
    entry_points={
        'console_scripts': [
            'dodos=src.launch:main',
        ],
    },
)