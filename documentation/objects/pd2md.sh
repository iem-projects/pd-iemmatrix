#!/bin/sh

helppatch="$1"

if ! test -e "${helppatch}"; then
    echo "'${helppatch}' does not exist!" 1>&2
    exit 1
fi

title="${helppatch%.pd}"
title="${title%-help}"
title="${title##*/}"

description="$(grep -h DESCRIPTION "${helppatch}" | sed -e 's|^#X text [0-9]* [0-9]* DESCRIPTION \(.*\)|\1|' -e 's|\(, f [1-9][0-9]*\)*;||'  -e 's| \\, |, |g')"

cat <<EOF
---
title: ${title}
description: ${description}
categories:
- object
pdcategory: Generic
see_also:
EOF

echo "inlets:"
grep "^#X text [0-9]* [0-9]* INLET_" "${helppatch}" | sed -e 's|#X text [0-9]* [0-9]* INLET_\(.*\);|\1|' | sort -u | while read num type; do
    if [ "${num}" = "0" ]; then
        num="1st"
    fi
    cat <<EOF
  ${num}:
  - type: ${type}
    description: ...
EOF
    done
echo "outlets:"
grep "^#X text [0-9]* [0-9]* OUTLET_" "${helppatch}" | sed -e 's|#X text [0-9]* [0-9]* OUTLET_\(.*\);|\1|' | sort -u | while read num type; do
    if [ "${num}" = "0" ]; then
        num="1st"
    fi
    cat <<EOF
  ${num}:
  - type: ${type}
    description: ...
EOF
done

cat <<EOF
draft: true
---
EOF
