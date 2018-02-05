FROM ubuntu:16.04

# update repositories
RUN apt-get update && \
    apt-get -y dist-upgrade

# install sudo
RUN apt-get install -y sudo

# development tools
RUN apt-get install -y vim less locales
RUN locale-gen en_US.UTF-8 en_GB.UTF-8 eu_FR.UTF-8

# for linking
ENV LD_LIBRARY_PATH /usr/local/lib
RUN export LD_LIBRARY_PATH

# run install script for cpp libraries
COPY install_cpplib.sh /root/
COPY silent_online.cfg /root/
RUN /bin/bash /root/install_cpplib.sh

# run install script for python
COPY install_python.sh /root/
COPY silent_online.cfg /root/
RUN /bin/bash /root/install_python.sh
