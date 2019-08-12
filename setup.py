from setuptools import setup, find_packages


def readme():
    with open('README.md') as f:
        return f.read()

setup(name='pypgatk',
      version='0.0.2',
      description='Python tools for proteogenomics ',
      url='http://github.com/bigbio/py-pgatk',
      long_description=readme(),
      long_description_content_type='text/markdown',
      author='PgAtk Team',
      author_email='ypriverol@gmail.com',
      license='Apache 2',
      install_requires=[
            'Click==7.0',
            'requests==2.21.0',
            'biopython==1.73'],
      scripts=['pypgatk/pypgatk_cli.py'],
      packages=find_packages(),
      zip_safe=False)
