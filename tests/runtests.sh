#!/bin/sh

if [ "x${PD}" = "x" ]
then
 PD=pd
fi

if [ "x${PDARGS}" = "x" ]
then
 PDARGS="-noprefs -nosound -nrt"
fi


RUNTESTS_TXT=runtests.txt
RUNTESTS_LOG=runtests.log
XITCODE=0

ls -1 */*.pd | sed 's/\.pd/;/' > $RUNTESTS_TXT

if [ "x${IEMMATRIX}" = "x" ]; then
 IEMMATRIX="-lib ../iemmatrix -path ../abs/"
fi

run_nogui() {
 ${PD} ${PDARGS} $IEMMATRIX -nogui runtests_nogui.pd > ${RUNTESTS_LOG} 2>&1
 NUMTESTS=`grep -c . $RUNTESTS_TXT`
 echo "regression-test: ${NUMTESTS} tests total" >>  ${RUNTESTS_LOG}

 cat ${RUNTESTS_LOG} | egrep "^regression-test: " | sed -e 's/^regression-test: //'
 FAILEDTESTS=$(cat ${RUNTESTS_LOG} | egrep "^regression-test: .*: failed$" | sed -e 's|^regression-test: ||' -e 's|: failed$||')
 echo -n "failed tests: "
 for ft in ${FAILEDTESTS}; do
   echo -n "${ft} "
 done
 echo
 if [ "x${FAILEDTESTS}" != "x" ]; then
   XITCODE=1
 fi
}

run_withgui() {
 ${PD} ${PDARGS} $IEMMATRIX -stderr runtests.pd 2>&1 | tee ${RUNTESTS_LOG}
}

if test "x$1" = "x-gui"; then
 run_withgui
else
 RUNTESTS_LOG=${RUNTESTS_LOG}.$(date +%Y%m%d-%H%M).$$
 run_nogui
fi

exit ${XITCODE}
