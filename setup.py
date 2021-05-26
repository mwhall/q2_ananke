from setuptools import setup, find_packages
import versioneer

setup(
    name="q2-ananke",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    packages=find_packages(),
    author="Michael Hall",
    author_email="mike.hall@dal.ca",
    description="Ananke time series clustering",
    license='BSD-3-Clause',
    url="https://qiime2.org",
    entry_points={
        'qiime2.plugins':
        ['q2-ananke=q2_ananke.plugin_setup:plugin']
    },
    package_data={'q2_ananke': ['citations.bib']},
    install_requires=['plotly','xxhash','bitarray','pandas','numpy','sklearn','biom-format','pybloom'],
    zip_safe=False,
)
