#!/bin/bash -e
SELF="$( cd "$(dirname "$0")" ; pwd -P )"
SUFFIX=""
case "$(uname -m)" in
  arm*|aarch*)
   SUFFIX="-aarch64"
   # upx is not working with samtools on aarch64
   gunzip "$SELF/samtools-linux${SUFFIX}.gz"
   ;;
esac
case "$(echo "$OSTYPE" | tr '[:upper:]' '[:lower:]')" in
  linux*) "$SELF/samtools-linux$SUFFIX" "$@" ;;
  darwin*) "$SELF/samtools-darwin" "$@" ;;
  msys*|cygwin*) "$SELF/samtools-windows" "$@" ;;
  *) exit 1 ;;
esac
