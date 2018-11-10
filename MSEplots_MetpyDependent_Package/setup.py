import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="MSEplots_pkg",
    version="1.1.7",
    author="Wei-Ming Tsai",
    author_email="wxt108@rsmas.miami.edu",
    description="A package for static energy plots based on MetPy",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/weiming9115/Working-Space/tree/master/MSEplots_metpy",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=['metpy']
)
