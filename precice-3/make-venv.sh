#!/usr/bin/env bash

python3 -m venv .venv
. .venv/bin/activate
pip install -e ../tooling
pip install -r requirements.txt