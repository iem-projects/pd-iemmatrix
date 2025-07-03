#!/bin/sh


: "${PD:=pd}"
: "${PDARGS:=-noprefs -nosound -nrt}"

: "${libdir:=..}"
: "${absdir:=../abs}"

: "${IEMMATRIX:=-lib "${libdir}/iemmatrix" -path "${absdir}"}"

RUNTESTS_TXT=runtests.txt
runtests_log=runtests.log
XITCODE=0

for f in */*.pd; do
  test -e "${f}" && echo "${f%.pd};"
done | LC_ALL=C sort > "${RUNTESTS_TXT}"

run_nogui() {
 echo ${PD} ${PDARGS} ${IEMMATRIX} -nogui runtests_nogui.pd > "${RUNTESTS_LOG}" 2>&1
 ${PD} ${PDARGS} ${IEMMATRIX} -nogui runtests_nogui.pd >> "${RUNTESTS_LOG}" 2>&1
 NUMTESTS=$(grep -c . "${RUNTESTS_TXT}")
 echo "regression-test: ${NUMTESTS} tests total" >>  "${RUNTESTS_LOG}"

 cat "${RUNTESTS_LOG}" \
     | grep -E "^regression-test: " \
     | sed -e 's/^regression-test: //'
 FAILEDTESTS=$(cat "${RUNTESTS_LOG}" \
                   | grep -E "^regression-test: .*: failed$" \
                   | sed -e 's|^regression-test: ||' -e 's|: failed$||')
 echo "failed tests: ${FAILEDTESTS}"
 if [ "x${FAILEDTESTS}" != "x" ]; then
   XITCODE=1
 fi
}

run_withgui() {
 ${PD} ${PDARGS} ${IEMMATRIX} -stderr runtests.pd 2>&1 | tee "${RUNTESTS_LOG}"
}

if test "$1" = "-gui"; then
 : "${RUNTESTS_LOG:=${runtests_log}}"
 run_withgui
else
 : "${RUNTESTS_LOG:=${runtests_log}.$(date +%Y%m%d-%H%M).$$}"
 run_nogui
fi

exit ${XITCODE}
