FROM python:2.7.9

MAINTAINER Haruka Ozaki <harukao.cb@gmail.com>

RUN apt-get update

RUN    pip install cython \
    && pip install numpy \
    && pip install Pillow \
    && pip install matplotlib \
    && pip install pandas \
    && pip install docopt

RUN pip install biopython

RUN curl -L http://sourceforge.net/projects/libpng/files/libpng15/1.5.28/libpng-1.5.28.tar.xz/download -o libpng.tar.xz \
	&& tar -Jxf libpng.tar.xz \
	&& cd libpng-1.5.28/ \
	&& ./configure --prefix=/usr \
	&& make \
	&& make install

RUN apt-get install -y libblas-dev liblapack-dev gfortran \
	&& pip install scipy \
	&& easy_install matplotlib-venn
