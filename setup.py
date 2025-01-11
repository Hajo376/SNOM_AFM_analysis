from setuptools import setup

setup(
    name="snom_analysis",
    version="0.1.0",
    description="Package for displaying and manipulating SNOM and AFM data.",
    packages=["snom_analysis", "snom_analysis/lib"],
    url="#",
    author="Hans-Joachim Schill",
    author_email="hajo.schill@acetovit.de",
    license="LICENSE",
    python_requires=">=3.12",
    zip_safe=False,
    include_package_data=True
)

# to create package just run <<python .\setup.py bdist_wheel>> in the command line, or add <<sdist>> for additional source distribution