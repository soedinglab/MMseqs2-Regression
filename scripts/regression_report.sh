#!/bin/bash -e
REPORT="${1}"
SEARCHROC="${2}"
PROFILEROC="${3}"

if [[ ! -s "${REPORT}" ]]; then
    exit 1
fi

ERROR=0
while IFS='' read -r line || [[ -n "$line" ]]; do
    IFS=$'\t' read -r -a values <<< "$line"
    if [[ "${values[2]}" -eq "0" ]]; then
        TARGET=$SEARCHROC
    elif [[ "${values[2]}" -eq "1" ]]; then
        TARGET=$PROFILEROC
    else
        continue
    fi
    GOOD=$(echo "${values[3]}" | awk -v target=$TARGET '{ print ($1 >= target) ? "GOOD" : "BAD" }')
    if [[ "$GOOD" != "GOOD" ]]; then
        ERROR=$((ERROR+1))
    fi
done < "$REPORT"

if [ $ERROR -ne 0 ]; then
    exit 1
fi

exit 0
