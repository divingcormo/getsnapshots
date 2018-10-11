#!/usr/bin/env python
# -*- encoding: utf-8 -*-
from datetime import date
from glob import glob ## for py_modules arg
import os

from setuptools import setup, find_packages

# This is only necessary when not using setuptools/distribute.
from sphinx.setup_command import BuildDoc
cmdclass = {'build_sphinx': BuildDoc}

package_name = 'getsnapshots'
version_from_module = 'snap2csv'
author = "Marius Wanko"
email = "marius.wanko@gmail.com"
description = "This is a test package wrapping the snap2csv script."
keyword_str = "CHARMM molecular dynamics trajectory data analysis"
project_url = "https://github.com/divingcormoran/getsnapshots/"
script_entry_points = [
    'get-snapshots = getsnapshots.snap2csv:main',
    'pyrmsd = getsnapshots.pyrmsd:main',
]

## how obtain this list automatically? tox? requirements.txt?
install_requires = [
    'pip>=9,<=18',
    'setuptools',
    'tqdm>=4.26',
    'numpy>=1.15',
]

classifiers = [
    # Comply with https://pypi.org/pypi?%3Aaction=list_classifiers
    # How mature is this project? Common values are
    #   1 - Planning
    #   2 - Pre-Alpha
    #   3 - Alpha
    #   4 - Beta
    #   5 - Production/Stable
    'Development Status :: 3 - Alpha',

    # Indicate who your project is intended for
    'Intended Audience :: Science/Research',
    'Topic :: Scientific/Engineering :: Chemistry'
    'Natural Language :: English',

    # Pick your license as you wish (should match "license" above)
    'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)'

    # Specify the Python versions you support here. In particular, ensure
    # that you indicate whether you support Python 2, Python 3 or both.
    'Programming Language :: Python :: 2',
    'Programming Language :: Python :: 2.7',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.4',
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: 3.6',
]


# Keep the code below generic -------------------------------------------------

package_pth = f'src/{package_name}'

# Get version.
""" 
In the src layout, the package modules cannot be imported like this:

    from mypackage.module import __version__

They could be found under "src.mypackage.module", but this would break the import
at a later point, when the module itself imports other modules as "mypackage.module".
In blog.ionelmc.ro/2014/06/25/python-packaging-pitfalls/ ionelmc warns (without 
mentioning this issue) to import package modules in setup.py and recommends to obtain 
the version number by read/parse that module.
"""
with open(f'{package_pth}/{version_from_module}.py') as fp:
    lines = fp.readlines()
for line in lines:
    if line.partition('=')[0].strip() in ('__version__', 'version'):
        __version__ = line.partition('=')[2].strip().split(' ')[0].strip('"').strip("'")
        # safer than: exec(line.strip())
        break
else:
    msg = "{} has no __version__ or version attribute"
    raise ValueError(msg.format(version_from_module))

for filename in 'README.rst README.md README.txt README'.split():
    if os.path.exists(f'{package_pth}/{filename}'):
        with open(f'{package_pth}/{filename}') as fp:
            README = fp.read()
        break
else:
    print("WARNING: no README file found")
    README = ''

# Get license (and version if applicable).
with open(f'{package_pth}/LICENSE') as fp:
    lines = fp.readlines()
for line in lines:
    if 'GNU GENERAL PUBLIC LICENSE' in line:
        license = 'GNU GPL'
        break
else:
    raise NotImplementedError("LICENSE file not recognized")
if license == 'GNU GPL':
    for line in lines:
        if line.strip().startswith('Version '):
            license += line.partition(',')[0].replace('Version ', 'v').strip()
            break

setup(
    name=package_name,
    version=__version__,
    package_dir={'': 'src'},
    packages=find_packages('src'),

    ## what are these good for?
    scripts=['src/getsnapshots/pyrmsd.py'],  # pyrmsd is included in built/scripts-3.6/
    # probably collects stand-alone script.py files. pyrmsd and get-snapshots could be
    # provided by either entry points (main function) or by using the module as script.
    #py_modules=[os.path.splitext(os.path.basename(path))[0] for path in glob('src/*.py')],
    # "packages" probably implies files that would be covered by "py_modules" and "scripts".
    # Also test/test*.py under src, README, README.txt, README.rst (3.7+), setup.py and setup.cfg
    # are included by default.

    install_requires=install_requires,

    # This may fail, better use MANIFEST.in and "include_package_data=True"
    #package_data={
    #    '': ['*.txt', '*.rst', '*.bash'],
    #    'tests': ['*'],
    #},
    include_package_data=True,

    # metadata to display on PyPI
    author=author,
    author_email=email,
    description=description,
    long_description=README,
    license=license,
    keywords=keyword_str,
    url=project_url,
    #project_urls={
    #    "Bug Tracker": "https://bugs.example.com/HelloWorld/",
    #    "Documentation": "https://docs.example.com/HelloWorld/",
    #    "Source Code": "https://code.example.com/HelloWorld/",
    #}
    classifiers=classifiers,

    # Sphinx (replace hyphens in conf.py variable names with underscores).
    cmdclass=cmdclass,
    command_options={
        'build_sphinx': {
            'project': ('setup.py', package_name),
            'version': ('setup.py', __version__),
            'release': ('setup.py', __version__),
            'source_dir': ('setup.py', 'docs/source'),
            'build_dir': ('setup.py', 'docs/build'),
            'copyright': ('setup.py', f'{date.today().year}, {author}'),
        },
    },

    # Entry points are grouped by custom strings (dict keys).
    entry_points={
        'console_scripts': script_entry_points,
    },
)
