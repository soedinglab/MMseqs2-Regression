#!/bin/bash -e
REPORT="${1}"
SEARCHROC="${2}"
PROFILEROC="${3}"

ERROR=0
while IFS='' read -r line || [[ -n "$line" ]]; do
    IFS=$"\t" read -r -a values <<< "$string"
    if [[ "${values[2]}" -eq "0" ]]; then
        TARGET=$SEARCHROC
    else
        TARGET=$PROFILEROC
    fi
    GOOD=$(echo "${values[3]}" | awk -v target=$TARGET '{ print ($1 > target) ? "GOOD" : "BAD" }')
    if [[ "$GOOD" == "GOOD" ]]; then
        ERROR=$((ERROR+1))
    fi
done < "$REPORT"

if [ $ERROR -ne 0 ]; then
    exit 1
fi

exit 0