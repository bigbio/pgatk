import io

from setuptools import setup, find_packages


def readme():
    with open('README.md') as f:
        return f.read()

setup(name='pypgatk',
      version='0.0.9',
      description='Python tools for proteogenomics',
      url='http://github.com/bigbio/py-pgatk',
      long_description=readme(),
      long_description_content_type='text/markdown',
      author='PgAtk Team',
      author_email='ypriverol@gmail.com',
      license='LICENSE.txt',
      include_package_data=True,
      install_requires=[
          'argcomplete==1.9.5',
          'argh==0.26.2',
          'asn1crypto==0.24.0',
          'astroid==2.2.5',
          'bcrypt==3.1.6',
          'biopython==1.73',
          'certifi==2019.3.9',
          'cffi==1.12.3',
          'chardet==3.0.4',
          'Click==7.0',
          'cryptography==3.2',
          'gffutils==0.9',
          'idna==2.8',
          'isort==4.3.18',
          'lazy-object-proxy==1.4.0',
          'mccabe==0.6.1',
          'numpy==1.16.3',
          'paramiko==2.4.2',
          'pyasn1==0.4.5',
          'pycparser==2.19',
          'pyfaidx==0.5.5.2',
          'pylint==2.3.1',
          'PyNaCl==1.3.0',
          'pysftp==0.2.9',
          'PyVCF==0.6.8',
          'PyYAML==5.1',
          'requests==2.21.0',
          'simplejson==3.16.0',
          'six==1.12.0',
          'typed-ast==1.3.5',
          'urllib3==1.24.2',
          'wrapt==1.11.1',
          'ratelimit'
      ],
      scripts=['pypgatk/pypgatk_cli.py'],
      packages=find_packages(),
      entry_points={
           'console_scripts': [
             'pypgatk_cli = pypgatk.pypgatk_cli:main'
      ]},
      package_data={'pypgatk':['config/*.yaml', 'config/*.json']},
      zip_safe=False)
