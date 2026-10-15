# Prebuilt installation (0.2.0)

[Download the release](https://github.com/forageseed/ovasm/releases/tag/v0.2.0).
These packages contain the executable and built-in seed references. Building from
source, installing Rust/C compilers and administrator permission are unnecessary.

## Windows x86_64

1. Download `ovasm-v0.2.0-windows-x86_64.zip` and extract it completely.
2. Double-click `install.cmd` in the extracted folder.
3. Open a new terminal and run `ovasm --version` (expected: `ovasm 0.2.0`).

The installer checks the executable against the packaged SHA-256 manifest,
copies it to `%LOCALAPPDATA%\Programs\ovasm` and adds that directory to the
current user's PATH. Existing PATH entries are preserved. No machine-wide PATH
change or administrator access is needed. The CRT is statically linked.

For an offline or custom-directory installation:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File .\install.ps1 -InstallDir 'D:\tools\ovasm'
```

Alternatively run `bin\ovasm.exe` in the extracted package directly. `-NoPathUpdate`
skips all persistent user PATH changes. Uninstall by removing the installed
`ovasm.exe` and its own PATH entry; the installer does not own your data/results.

## Linux x86_64

Requires glibc >=2.34 and normal system runtime libraries, including `libz.so.1`.
The package is built on Ubuntu 24.04; the minimum glibc symbol requirement is
checked on the actual executable. C/zlib development headers are build-time
requirements only and are not needed to run this binary.

```sh
curl -fsSL https://raw.githubusercontent.com/forageseed/ovasm/v0.2.0/scripts/install.sh | sh
export PATH="$HOME/.local/bin:$PATH"
ovasm --version
```

The script checks the platform, downloads the release archive plus `SHA256SUMS`,
verifies SHA-256, checks the extracted executable version and installs only the
binary in `~/.local/bin`. It never calls sudo or invokes a compiler. To use a
different location or an already-downloaded package:

```sh
sh scripts/install.sh --bin-dir /path/to/my/bin \
  --archive /path/to/ovasm-v0.2.0-linux-x86_64.tar.gz \
  --checksums /path/to/SHA256SUMS
```

You can also extract the tar.gz and run its `bin/ovasm` directly. No installer
modifies shell startup files; the exported PATH applies to the current terminal.
For future terminals add that directory to your shell PATH if not already present.

## Platform scope

Linux includes the native multichrom kernel. Windows includes the portable
stages, not the Linux-only kernel. This release provides no macOS, ARM64, or
musl/Alpine executable. Platform availability is separate from biological
acceptance: inspect evidence, unresolved branches and sequence scores for each run.
