#!/bin/bash -e
REPORT="${1}"
SEARCHROC="${2}"
PROFILEROC="${3}"
DBPROFROC="${4}"
CLUNUM="${5}"
LINCLUNUM="${6}"

if [[ ! -s "${REPORT}" ]]; then
    exit 1
fi

AWKLEQ='{ print ($1 >= target) ? "GOOD" : "BAD" }'
AWKNEQ='{ p1 = target * 0.005; print ($1 <= (target + p1) && $1 >= (target - p1)) ? "GOOD" : "BAD" }'

ERROR=0
while IFS='' read -r line || [[ -n "$line" ]]; do
    IFS=$'\t' read -r -a values <<< "$line"
    TARGET=""
    TEST=""
    if   [[ "${values[2]}" -eq "0" ]]; then
        TARGET=$SEARCHROC
        TEST=AWKLEQ
    elif [[ "${values[2]}" -eq "1" ]]; then
        TARGET=$PROFILEROC
        TEST=AWKLEQ
    elif [[ "${values[2]}" -eq "9" ]]; then
        TARGET=$DBPROFROC
        TEST=AWKLEQ
    elif [[ "${values[2]}" -eq "3" ]]; then
        TARGET=$CLUNUM
        TEST=AWKNEQ
    elif [[ "${values[2]}" -eq "6" ]]; then
        TARGET=$LINCLUNUM
        TEST=AWKNEQ
    else
        continue
    fi

    GOOD=$(echo "${values[3]}" | awk -v target=$TARGET "${!TEST}")
    if [[ "$GOOD" != "GOOD" ]]; then
        ERROR=$((ERROR+1))
         >&2 echo "Failed check! Input: ${values[3]} Expected: $TARGET Comparison: ${TEST}"
         continue
    fi
done < "$REPORT"

if [ $ERROR -ne 0 ]; then
    exit 1
fi

exit 0
