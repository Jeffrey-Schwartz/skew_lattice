#! /bin/sh

aclocal \
&& automake --add-missing \
&& autoconf \
&& ./configure \
&& dos2unix libtool \
&& dos2unix -n depcomp tmpFile \
&& mv tmpFile depcomp
