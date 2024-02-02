from setuptools import setup

setup(
    name="SNOM_AFM_analysis",
    version="0.0.52",
    description="Package for displaying and manipulating SNOM and AFM data.",
    # package_dir={},
    packages=["SNOM_AFM_analysis", "SNOM_AFM_analysis/lib"],
    url="#",
    author="Hajo Schill",
    author_email="hajo.schill@acetovit.de",
    license="",
    install_requires=["numpy", "matplotlib", "scipy", "pandas"], # these ones cause errors? :, "pathlib", "tkinter", "enum", "datetime", "struct", "os" 
    python_requires=">=3.7",
    zip_safe=False,
    include_package_data=True
)

# to create package just run <<python .\setup.py bdist_wheel>> in the command line, or add <<sdist>> for additional source distribution