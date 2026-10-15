#!/bin/sh
# Install an official prebuilt package; never invoke a compiler or sudo.
set -eu
version=0.2.0
bin_dir=${OVASM_INSTALL_DIR:-"$HOME/.local/bin"}
archive=
checksums=
die() { printf '%s\n' "$*" >&2; exit 1; }
while [ "$#" -gt 0 ]; do
    case "$1" in
        --version|--bin-dir|--archive|--checksums)
            [ "$#" -ge 2 ] || die "Missing value for $1"
            case "$1" in
                --version) version=${2#v} ;;
                --bin-dir) bin_dir=$2 ;;
                --archive) archive=$2 ;;
                --checksums) checksums=$2 ;;
            esac
            shift 2 ;;
        --help)
            echo 'Usage: install.sh [--bin-dir DIR] [--version 0.2.0] [--archive FILE --checksums SHA256SUMS]'
            exit 0 ;;
        *) die "Unknown option: $1" ;;
    esac
done
printf '%s' "$version" | grep -Eq '^[0-9]+\.[0-9]+\.[0-9]+$' || die 'Invalid release version'
[ "$(uname -s)" = Linux ] || die 'This installer provides Linux x86_64 only; see GitHub Releases.'
case "$(uname -m)" in x86_64|amd64) ;; *) die 'No prebuilt package for this CPU architecture.' ;; esac
if command -v getconf >/dev/null 2>&1; then
    libc=$(getconf GNU_LIBC_VERSION 2>/dev/null || true)
    [ -z "$libc" ] || printf '%s\n' "$libc" | awk '{split($2,v,".");exit !(v[1]>2 || (v[1]==2 && v[2]>=34))}' || die 'glibc 2.34 or newer is required.'
fi
package="ovasm-v$version-linux-x86_64"
asset="$package.tar.gz"
tempdir=$(mktemp -d "${TMPDIR:-/tmp}/ovasm-install.XXXXXXXX")
trap 'rm -rf "$tempdir"' EXIT HUP INT TERM
if [ -n "$archive" ]; then
    [ -n "$checksums" ] || die '--checksums is required with --archive'
    [ -f "$archive" ] && [ -f "$checksums" ] || die 'Local archive/checksum file missing'
else
    [ -z "$checksums" ] || die '--archive is required with --checksums'
    command -v curl >/dev/null 2>&1 || die 'curl is required; alternatively download and extract the release archive.'
    base="https://github.com/forageseed/ovasm/releases/download/v$version"
    archive="$tempdir/$asset"
    checksums="$tempdir/SHA256SUMS"
    curl --proto '=https' --tlsv1.2 -fsSL "$base/$asset" -o "$archive"
    curl --proto '=https' --tlsv1.2 -fsSL "$base/SHA256SUMS" -o "$checksums"
fi
expected=$(awk -v n="$asset" 'NF==2 && $2==n {print $1}' "$checksums")
printf '%s\n' "$expected" | grep -Eq '^[a-fA-F0-9]{64}$' || die 'Missing or invalid SHA256SUMS entry'
if command -v sha256sum >/dev/null 2>&1; then
    actual=$(sha256sum "$archive" | awk '{print $1}')
elif command -v shasum >/dev/null 2>&1; then
    actual=$(shasum -a 256 "$archive" | awk '{print $1}')
else
    die 'sha256sum or shasum is required for package verification.'
fi
[ "$actual" = "$expected" ] || die 'SHA-256 mismatch; nothing installed.'
tar -tzf "$archive" > "$tempdir/entries"
awk '/^\// || /(^|\/)\.\.(\/|$)/ {bad=1} END {exit bad}' "$tempdir/entries" || die 'Unsafe archive path'
tar -xzf "$archive" -C "$tempdir"
binary="$tempdir/$package/bin/ovasm"
[ -f "$binary" ] || die 'Release package has no bin/ovasm'
[ "$("$binary" --version)" = "ovasm $version" ] || die 'Unexpected executable version'
mkdir -p "$bin_dir"
install -m 755 "$binary" "$bin_dir/ovasm"
printf 'Installed %s to %s\n' "$("$bin_dir/ovasm" --version)" "$bin_dir/ovasm"
case ":$PATH:" in
    *":$bin_dir:"*) ;;
    *) printf 'For this terminal: export PATH="%s:$PATH"\n' "$bin_dir" ;;
esac
