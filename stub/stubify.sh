#!/bin/sh

infile="$1"
module="${infile%.c}"
module="${module##*/}"

stubfuns="$(cat "${infile}" | egrep "^STUB" | sed -e 's|^STUB(\(.*\));|\1|')"


echo "/* ${module} */"
for fun in ${stubfuns}; do
	echo "DECLARE_STUB($fun);"
done
cat <<EOF

static void setup_${module}() {
  module_t mod = getmodule("${module}");

EOF
for fun in ${stubfuns}; do
	echo "  MAKE_STUB($fun, mod, $fun);"
done
cat <<EOF
}



// this goes into iemmatrix_get_stub():

EOF

echo "  /* ${module} */"
for fun in ${stubfuns}; do
	echo "  GET_STUB($module, $fun, name);"
done

