import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="ERDC-MOOSE", # Replace with your own username
    version="1.0.0",
    author="Theodore Letcher",
    author_email="",
    description="A simple package designed to handle file I/O and geographic plotting for gridded meteorological datasets.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/wxted/ERDC_MOOSE",
    packages=['MOOSEplot','MOOSEpost','MOOSEfunctions'],
    package_dir={'MOOSEplot':'ERDCMOOSE/plotting','MOOSEpost':'ERDCMOOSE/post_processing','MOOSEfunctions':'ERDCMOOSE/functions/'},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.4',
)
