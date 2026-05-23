//! Compatibility tests: byte-identical output vs bedtools summary v2.31.1.
//!
//! Golden fixtures generated with:
//!   bedtools summary -i tests/golden/input.bed -g tests/golden/genome.txt
//!   bedtools summary -i tests/golden/small.bed  -g tests/golden/small_genome.txt
//!
//! The `bedtools_version_ok` guard skips when bedtools is absent or not v2.31
//! (CI runners don't have bedtools installed).

use std::path::Path;
use std::process::Command;

use rsomics_bed_summary::summary;

fn golden(name: &str) -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

fn run_summary(bed: &str, genome: &str) -> String {
    let mut buf = Vec::new();
    summary(&golden(bed), &golden(genome), &mut buf).expect("summary failed");
    String::from_utf8(buf).expect("valid utf8")
}

fn expected(name: &str) -> String {
    std::fs::read_to_string(golden(name)).expect("golden file missing")
}

#[test]
fn main_fixture_golden() {
    assert_eq!(
        run_summary("input.bed", "genome.txt"),
        expected("summary.tsv"),
        "main fixture mismatch vs bedtools summary"
    );
}

#[test]
fn small_fixture_golden() {
    assert_eq!(
        run_summary("small.bed", "small_genome.txt"),
        expected("small_summary.tsv"),
        "small fixture mismatch vs bedtools summary"
    );
}

/// Byte-identical live check against the installed bedtools binary.
///
/// Skipped gracefully when bedtools is absent or not v2.31.x.
#[test]
fn live_bedtools_compat() {
    let bt = match which_bedtools() {
        Some(p) => p,
        None => {
            eprintln!("SKIP live_bedtools_compat: bedtools not found");
            return;
        }
    };
    if !bedtools_version_ok(&bt) {
        eprintln!("SKIP live_bedtools_compat: bedtools version mismatch (want v2.31.x)");
        return;
    }

    let bed = golden("input.bed");
    let genome = golden("genome.txt");

    let upstream = Command::new(&bt)
        .args(["summary", "-i"])
        .arg(&bed)
        .arg("-g")
        .arg(&genome)
        .output()
        .expect("bedtools summary failed");
    assert!(upstream.status.success(), "bedtools exited non-zero");

    let mut ours = Vec::new();
    summary(&bed, &genome, &mut ours).expect("our summary failed");

    assert_eq!(
        String::from_utf8(ours).expect("utf8"),
        String::from_utf8(upstream.stdout).expect("utf8"),
        "live byte-identical check failed"
    );
}

fn which_bedtools() -> Option<String> {
    let out = Command::new("which").arg("bedtools").output().ok()?;
    if out.status.success() {
        let p = String::from_utf8(out.stdout).ok()?.trim().to_string();
        if p.is_empty() { None } else { Some(p) }
    } else {
        None
    }
}

fn bedtools_version_ok(bt: &str) -> bool {
    let out = Command::new(bt)
        .arg("--version")
        .output()
        .expect("bedtools --version failed");
    let v = String::from_utf8_lossy(&out.stdout);
    v.contains("v2.31")
}
