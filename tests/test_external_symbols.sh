#!/bin/sh

# Check for external symbols without uvg_ prefix (or _uvg_ on macOS).

set -eu${BASH+o pipefail}

if nm -go --defined-only ../lib/libuvg266.a | grep -v ' uvg_' | grep -v ' _uvg_'; then
    printf '%s\n' 'Only symbols prefixed with "uvg_" should be exported from libuvg266.'
    false
fi
