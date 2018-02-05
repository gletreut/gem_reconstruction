#!/bin/bash
apt-get update
apt-get -y dist-upgrade
apt-get -y install build-essential checkinstall

# build python from source
apt-get install -y python python-pip python-tk

# upgrade pip
pip install wheel
pip install --upgrade pip

# install numpy and scipy
pip install numpy
pip install scipy

# install matplotlib
pip install matplotlib

# install PIL
pip install Pillow

# install YAML
pip install pyyaml

