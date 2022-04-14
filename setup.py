from setuptools import setup, find_packages


with open('requirements.txt', 'r') as f:
    required = f.read().splitlines()

setup(
    name = 'morph_utils',
    version = '0.1.0',
    description = """Functions for common and not so common morphology operations""",
    author = "Matthew Mallory",
    author_email = "matt.mallory@alleninstitute.org",
    url = '',
    packages = find_packages(),
    install_requires = required,
    include_package_data=True,
    setup_requires=['pytest-runner'],
)