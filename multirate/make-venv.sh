#!/usr/bin/env bash

python3 -m venv .venv
. .venv/bin/activate
pip install -e ../tooling
# still requires publishing to PyPI. Currently need to clone from https://github.com/BenjaminRodenberg/doConvergenceStudy
pip install -e ~/Programming/doConvergenceStudy/.
pip install -r requirements.txt
