#!/usr/bin/env bash

sed 's/\\ensuremath{\\backslash}/\\/g; s/\\\$/$/g; s/\\_/_/g; s/\\\^{}/^/g; s/\\{/{/g; s/\\}/}/g' "$@"
