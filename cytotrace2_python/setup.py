#!/usr/bin/env python

from setuptools import setup, find_packages 

setup(
    name='cytotrace2_py',
    version='1.1.0',
    packages=['cytotrace2_py',
             'cytotrace2_py.common'],
    package_dir={"cytotrace2_py": "cytotrace2_py"},
    package_data={"cytotrace2_py":["resources/*","resources/models/*","resources/background.pt"]},
    entry_points={'console_scripts':
                      ['cytotrace2=cytotrace2_py.cytotrace2_py:run_cytotrace2']
                  },
    include_package_data=True
)
