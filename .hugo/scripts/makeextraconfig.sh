#!/bin/sh

# creates additional configuration for hugo
# (to be used in a CI)

docfolder="documentation/objects"

scriptdir="$(cd $(dirname "$0"); pwd)"
hugodir="${scriptdir}/.."
objectsdir="${hugodir}/../${docfolder}"

total_objects="$(find "${objectsdir}" -type f -name "*.md" -exec echo .  {} ";" | grep -c .)"
draft_objects="$(cd "${hugodir}"; hugo list drafts 2>/dev/null | grep -c "${docfolder}")"


# generate output
cat <<EOF
[Params.objectdocs]
total = ${total_objects}
draft = ${draft_objects}

EOF
