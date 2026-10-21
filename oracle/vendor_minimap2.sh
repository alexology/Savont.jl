#!/bin/sh
# Copy the minimap2 Rust wrapper (0.1.31, minimap2 2.30) from the cargo registry
# into rust_upstream/vendor/minimap2 and make it build on Windows: the wrapper
# converts index paths with a Unix-only API. The change affects path handling
# only, not alignment. Run once before building the oracle (needs the crate in
# the cargo cache: `cargo fetch` in rust_upstream downloads it).
set -e
here=$(cd "$(dirname "$0")" && pwd)
source=$(ls -d "${CARGO_HOME:-$HOME/.cargo}"/registry/src/*/minimap2-0.1.31+minimap2.2.30 | head -1)
target="$here/rust_upstream/vendor/minimap2"
rm -rf "$target"
mkdir -p "$(dirname "$target")"
cp -r "$source" "$target"
sed -i 's/^use std::os::unix::ffi::OsStrExt;/#[cfg(unix)]\nuse std::os::unix::ffi::OsStrExt;/' "$target/src/lib.rs"
sed -i 's/std::ffi::CString::new(path.as_ref().as_os_str().as_bytes())/std::ffi::CString::new(path.as_ref().to_str().expect("UTF-8 path").as_bytes())/' "$target/src/lib.rs"
grep -q 'UTF-8 path' "$target/src/lib.rs"
echo "patched $target"
