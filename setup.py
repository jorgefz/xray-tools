import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="xraytools",
    version="0.2",
    author="Jorge Fernandez",
    author_email="Jorge.Fernandez-Fernandez@warwick.ac.uk",
    
    package_dir={"": "."},
    packages=setuptools.find_packages(where="."),

    #packages=['xraytools'],
    description="Set of Python tools to work with X-ray data.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/jorgefz/xraytools",
    license="LICENSE",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GPL License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
    install_requires=[
       "astropy",
       "numpy",
       "scipy",
       "matplotlib"
   ],
)