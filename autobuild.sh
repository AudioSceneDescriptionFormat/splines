#!/bin/sh

# Trying to avoid invalidating Sphinx's cache between commits and on errors

PYTHON=python3

$PYTHON -m sphinx_autobuild doc _build \
    -Drelease=dummy \
    -Dversion=dummy \
    -Dtoday=dummy \
    -Dhtml_title=splines-dummy \
    -Dnbsphinx_allow_errors=1 \
    $@
