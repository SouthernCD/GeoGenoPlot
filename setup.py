# coding utf8
import setuptools
from geogenoplot.versions import get_versions

with open('README.md') as f:
    LONG_DESCRIPTION = f.read()

setuptools.setup(
    name="GeoGenoPlot",
    version=get_versions(),
    author="Yuxing Xu",
    author_email="xuyuxing@mail.kib.ac.cn",
    description="A python package for plotting correlations between variants and environmental factors.",
    long_description=LONG_DESCRIPTION,
    long_description_content_type='text/markdown',
    url="https://github.com/SouthernCD/GeoGenoPlot",
    include_package_data = True,

    entry_points={
        "console_scripts": ["geogenoplot = geogenoplot.cli:main"]
    },    

    packages=setuptools.find_packages(),

    install_requires=[
        "yxmap",
        "yxutil",
        "cyvcf2>=0.30.28"
    ],

    python_requires='>=3.5',
)