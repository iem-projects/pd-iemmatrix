#!/bin/sh

AUTORECONF=$(which autoreconf)

if [ -x "${AUTORECONF}" ]; then
  ${AUTORECONF} -fiv || exit 1
else
  aclocal && autoconf || exit 1
fi
echo "now run './configure'" 1>&2
echo "for help on args run './configure --help'" 1>&2
