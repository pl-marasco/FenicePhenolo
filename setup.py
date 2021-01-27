import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="Phenolo",
    version="2.5.2",
    author="Pier Lorenzo Marasco",
    author_email="pl.marasco@gmail.com",
    description="Phenology caclulator according to the shift methodology",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/pypa/Phenolo",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GPL-3.0-or-later",
        "Operating System :: OS Independent",
    ],
)