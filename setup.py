import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="xraytools-jorgefz",
    version="0.0.1",
    author="Jorge Fernandez",
    author_email="Jorge.Fernandez-Fernandez@warwick.ac.uk",
    description="A small example package",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/jorgefz/xraytools",
    project_urls={
        "Bug Tracker": "https://github.com/jorgefz/xraytools/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GPL License",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "xraytools"},
    packages=setuptools.find_packages(where="xraytools"),
    python_requires=">=3.6",
)