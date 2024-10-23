from setuptools import setup, find_packages

setup(
    name="KmerDecon",
    version="0.1.0",
    author="Yuxiang Li, Yujia Feng, Xiaoyi Chen",
    description="A fast, memory-efficient tool for decontaminating sequencing reads using Bloom filters.",
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    url="https://github.com/skysky2333/KmerDecon",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    install_requires=[
        "bitarray>=2.1.0",
        "biopython>=1.78",
        "mmh3>=2.5.1",
        "hyperloglog>=0.0.12"
    ],
    entry_points={
        'console_scripts': [
            'kbuild=KmerDecon.build_bloom_filter:main',
            'kdecon=KmerDecon.decontaminate_reads:main',
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
)