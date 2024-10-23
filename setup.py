
from setuptools import setup, find_packages

setup(
    name='KmerDecon',
    version='0.1.0',
    description='Fast, memory-efficient decontamination of sequencing reads using Bloom filters',
    author='Yujia Feng, Xiaoyi Chen, Yuxiang Li',
    author_email='yli694@jh.edu',
    url='https://github.com/skysky2333/KmerDecon',
    license='MIT',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    install_requires=[
        'bitarray>=2.1.0',
        'biopython>=1.78',
        'mmh3>=2.5.1',
        'hyperloglog>=0.0.12',
    ],
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
    ],
    python_requires='>=3.6',
)
