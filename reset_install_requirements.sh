#!/usr/bin/env bash
# Build the checked-in Basilisk tree without replacing its local safety policy.
set -euo pipefail

if [[ $# -ne 0 ]]; then
    printf '%s\n' 'Usage: bash reset_install_requirements.sh (no destructive reset)' >&2
    exit 2
fi
PROJECT_ROOT="$(cd -- "$(dirname -- "$0")" && pwd -P)"
BASILISK="$PROJECT_ROOT/basilisk/src"
if [[ ! -f "$BASILISK/Makefile" ]] ||
   ! grep -Fq 'Basilisk network display is disabled in this repository' "$BASILISK/display.h"; then
    printf '%s\n' 'The hardened vendored Basilisk tree is missing. Restore it from a trusted checkout.' >&2
    exit 1
fi

(
    cd -- "$BASILISK"
    if [[ ! -e config && ! -L config ]]; then
        case "$(uname -s)" in
            Darwin) ln -s config.osx config ;;
            *) ln -s config.gcc config ;;
        esac
    fi
    make all
)
{
    printf 'export PROJECT_ROOT=%q\n' "$PROJECT_ROOT"
    printf 'export BASILISK=%q\n' "$BASILISK"
    printf 'export PATH="$BASILISK:$PATH"\n'
} > "$PROJECT_ROOT/.project_config"
printf '%s\n' 'Built the vendored Basilisk toolchain without network display targets.'
