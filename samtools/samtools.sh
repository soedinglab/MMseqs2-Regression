#!/bin/bash -e
SELF="$( cd "$(dirname "$0")" ; pwd -P )"
SUFFIX=""
case "$(uname -m)" in
  arm*|aarch*) SUFFIX="-aarch64" ;;
  ppc*) SUFFIX="-ppc64le" ;;
esac
case "$(echo "$OSTYPE" | tr '[:upper:]' '[:lower:]')" in
  linux*) "$SELF/samtools-linux$SUFFIX" "$@" ;;
  darwin*) "$SELF/samtools-darwin" "$@" ;;
  msys*|cygwin*) "$SELF/samtools-windows" "$@" ;;
  *) exit 1 ;;
esac
