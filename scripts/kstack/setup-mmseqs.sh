#!/bin/zsh
# Explicitly provision the pinned MMseqs2 release that Kstack bundles: verify the archive,
# thin the universal binary to arm64 and record release, archive and thinned-binary SHA-256 in
# .tools/mmseqs/provenance.json (the builder copies it and records it in the install manifest).
# A previously downloaded archive is reused when its checksum still matches.
set -euo pipefail
cd "${0:A:h:h:h}"
release=18-8cc5c
url=https://github.com/soedinglab/MMseqs2/releases/download/$release/mmseqs-osx-universal.tar.gz
expected=158102d4c3e3ce721424f455106fac348cd01d424ed5e66cb03eac3de5064eb6
archive=.tools/mmseqs-osx-universal.tar.gz
digest() { shasum -a 256 "$1" | cut -d ' ' -f 1 }

mkdir -p .tools
if [[ ! -f "$archive" || "$(digest "$archive")" != "$expected" ]]; then
  rm -f "$archive.partial"
  curl -fL --retry 3 "$url" -o "$archive.partial"
  mv -f "$archive.partial" "$archive"
fi
[[ "$(digest "$archive")" == "$expected" ]] || { print -u2 'MMseqs2 checksum mismatch'; exit 1; }

staging=$(mktemp -d .tools/.mmseqs-setup.XXXXXX)
trap 'rm -rf "$staging"' EXIT
tar -xzf "$archive" -C "$staging"
binary=$staging/mmseqs/bin/mmseqs
lipo "$binary" -thin arm64 -output "$binary.arm64"
mv -f "$binary.arm64" "$binary"
chmod 755 "$binary"
[[ "$(lipo -archs "$binary")" == arm64 ]] || { print -u2 'MMseqs2 is not thin arm64'; exit 1; }
reported=$("$binary" version)
[[ "$reported" =~ '^[0-9a-f]+$' ]] || { print -u2 "Unexpected MMseqs2 version: $reported"; exit 1; }
cat > "$staging/mmseqs/provenance.json" <<EOF
{
  "name": "MMseqs2",
  "version": "$release",
  "reported_version": "$reported",
  "url": "$url",
  "archive_sha256": "$expected",
  "architecture": "arm64",
  "binary": "bin/mmseqs",
  "binary_sha256": "$(digest "$binary")"
}
EOF
rm -rf .tools/mmseqs.previous
[[ -e .tools/mmseqs ]] && mv .tools/mmseqs .tools/mmseqs.previous
mv "$staging/mmseqs" .tools/mmseqs
rm -rf .tools/mmseqs.previous
print "MMseqs2 $release ($reported), thin arm64: .tools/mmseqs/bin/mmseqs"
