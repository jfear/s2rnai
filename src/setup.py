#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup

requirements = [i.strip() for i in open("requirements.txt").readlines()]

setup(
    name="s2rnai",
    version="0.0.1",
    description="Local library for the s2rnai project",
    author="Justin M Fear",
    author_email="justin.m.fear@gmail.com",
    url="https://github.com/jfear/s2rnai",
    packages=["s2rnai"],
    install_requires=requirements,
    license="MIT license",
)
