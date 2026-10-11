use std::{env, fs, path::PathBuf, process::Command};

fn main() {
    println!("cargo:rustc-check-cfg=cfg(syncasm_native)");
    println!("cargo:rerun-if-changed=vendor/syncasm");
    println!("cargo:rerun-if-env-changed=CC");
    println!("cargo:rerun-if-env-changed=AR");
    // The scientific provider runs on Linux. Other targets retain the portable Rust stages
    // and return an explicit unsupported-provider error for the published graph kernel.
    if env::var("CARGO_CFG_TARGET_OS").as_deref() != Ok("linux") {
        return;
    }
    let compiler = env::var("CC").unwrap_or_else(|_| "cc".into());
    let archiver = env::var("AR").unwrap_or_else(|_| "ar".into());
    let root = PathBuf::from("vendor/syncasm");
    let out = PathBuf::from(env::var_os("OUT_DIR").unwrap());
    let sources = [
        "run_syncasm",
        "syncasm",
        "syncmer",
        "syncerr",
        "levdist",
        "graph",
        "alignment",
        "sstream",
        "misc",
        "kthread",
        "kalloc",
        "kopen",
        "ovasm_bridge",
    ];
    let mut objects = Vec::new();
    for source in sources {
        let object = out.join(format!("{source}.o"));
        let status = Command::new(&compiler)
            .args(["-O3", "-fPIC", "-std=gnu11", "-Wno-unused-function", "-I"])
            .arg(&root)
            .arg("-c")
            .arg(root.join(format!("{source}.c")))
            .arg("-o")
            .arg(&object)
            .status()
            .expect("cannot start C compiler for the native syncmer kernel");
        assert!(
            status.success(),
            "native syncmer kernel compilation failed: {source}"
        );
        objects.push(object);
    }
    let library = out.join("libovasm_syncasm.a");
    if library.exists() {
        fs::remove_file(&library).unwrap();
    }
    let status = Command::new(archiver)
        .arg("crs")
        .arg(&library)
        .args(objects)
        .status()
        .expect("cannot archive native syncmer kernel");
    assert!(status.success(), "native syncmer kernel archive failed");
    println!("cargo:rustc-link-search=native={}", out.display());
    println!("cargo:rustc-link-lib=static=ovasm_syncasm");
    for library in ["z", "m", "pthread"] {
        println!("cargo:rustc-link-lib={library}");
    }
    println!("cargo:rustc-cfg=syncasm_native");
}
