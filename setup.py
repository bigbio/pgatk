import io

from setuptools import setup, find_packages


def readme():
    with open('README.md') as f:
        return f.read()


def read_requirements_txt(default=None):
    if default is None:
        default = []
    with io.open('requirements.txt', 'r') as f:
        r = f.read().split()
    return r if r else default


install_requires = []
requirements_txt = read_requirements_txt()
install_requires.extend(requirements_txt)

setup(name='pypgatk',
      version='0.0.3',
      description='Python tools for proteogenomics ',
      url='http://github.com/bigbio/py-pgatk',
      long_description=readme(),
      long_description_content_type='text/markdown',
      author='PgAtk Team',
      author_email='ypriverol@gmail.com',
      license='Apache 2',
      install_requires=install_requires,
      scripts=['pypgatk/pypgatk_cli.py'],
      packages=find_packages(),
      zip_safe=False)
