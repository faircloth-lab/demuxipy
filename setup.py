import distribute_setup
distribute_setup.use_setuptools()
from setuptools import setup
from setuptools import find_packages
from setuptools.command.build_ext import build_ext
from setuptools.extension import Extension

if __name__ == '__main__':
    setup(
        name='demuxipy',
        version="1.0",
        description="Demultiplex hierarchically sequence-tagged massively"
            +"parallel sequencing reads",
        author="Brant Faircloth",
        author_email="brant.faircloth+demuxipy@gmail.com ",
        url="http://baddna.github.com/demuxipy/",
        license="BSD",
        classifiers=[
            'Development Status :: 4 - Beta',
            'Environment :: Console',
            'Operating System :: OS Independent',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: BSD License',
            'Programming Language :: Python',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
             ],
        requires=['NumPy (>=1.3)',],
        long_description=open('README.txt').read(),
        package_data = {
                # If any package contains *.txt or *.rst files, include them:
                '': ['*.txt', '*.rst'],
                'demuxipy': ['tests/test-data/*.txt'],
            },
        packages=['demuxipy',
                'demuxipy.tests',
                ],
        scripts=['bin/demuxi.py',
                ],
        ext_modules=[
                Extension(
                    'demuxipy/cpairwise2',
                    ['demuxipy/cpairwise2module.c']
                    ),
                ]
        )
