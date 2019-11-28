#!/bin/bash -e
SELF="$( cd "$(dirname "$0")" ; pwd -P )"
case "$(echo "$OSTYPE" | tr '[:upper:]' '[:lower:]')" in
  linux*) "$SELF/samtools-linux" "$@" ;;
  darwin*) "$SELF/samtools-darwin" "$@" ;;
  msys*|cygwin*) "$SELF/samtools-windows" "$@" ;;
  *) exit 1 ;;
esac
