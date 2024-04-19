import io

from setuptools import setup, find_packages


def readme():
  with open('README.md') as f:
    return f.read()


setup(name='pypgatk',
      version='0.0.24',
      description='Python tools for proteogenomics',
      url='http://github.com/bigbio/py-pgatk',
      long_description=readme(),
      long_description_content_type='text/markdown',
      author='PgAtk Team',
      author_email='ypriverol@gmail.com',
      license='LICENSE.txt',
      include_package_data=True,
      install_requires=[
        'biopython',
        'Click',
        'gffutils',
        'numpy',
        'pandas',
        'PyYAML',
        'requests',
        'simplejson',
        'ratelimit',
        'pyteomics',
        'pathos',
        'pybedtools',
        'pyopenms',
        'matplotlib',
        'tqdm',
        'pyahocorasick'
      ],
      python_requires=">=3.6",
      scripts=['pypgatk/pypgatk_cli.py'],
      packages=find_packages(),
      entry_points={
        'console_scripts': [
          'pypgatk = pypgatk.pypgatk_cli:main'
        ]},
      package_data={'pypgatk': ['config/*.yaml', 'config/*.json']}, zip_safe=False)
