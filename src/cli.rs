use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, Tool, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};

use rsomics_bed_summary::summary;

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

/// bedtools summary uses multi-char single-dash flags (`-i`, `-g`);
/// clap supports single-char shorts so `-i`/`-g` map directly.
#[derive(Parser, Debug)]
#[command(
    name = "rsomics-bed-summary",
    version,
    about = "Statistical summary of BED intervals per chromosome — Rust port of bedtools summary",
    long_about = None,
    disable_help_flag = true
)]
pub struct Cli {
    /// Input BED file (intervals to summarise).
    ///
    /// Maps to bedtools `-i`.
    #[arg(short = 'i', long = "input", value_name = "FILE")]
    pub input: PathBuf,

    /// Genome file (`chrom<TAB>size` lines, one per chromosome).
    ///
    /// Maps to bedtools `-g`.
    #[arg(short = 'g', long = "genome", value_name = "FILE")]
    pub genome: PathBuf,

    /// Output file (default: stdout).
    #[arg(short = 'o', long, default_value = "-")]
    pub output: String,

    #[command(flatten)]
    pub common: CommonFlags,
}

impl Tool for Cli {
    fn meta() -> ToolMeta {
        META
    }

    fn common(&self) -> &CommonFlags {
        &self.common
    }

    fn execute(self) -> Result<()> {
        let mut out: Box<dyn std::io::Write> = if self.output == "-" {
            Box::new(std::io::stdout().lock())
        } else {
            Box::new(std::fs::File::create(&self.output).map_err(RsomicsError::Io)?)
        };
        summary(&self.input, &self.genome, &mut out)
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
    tagline: "Statistical summary of BED intervals per chromosome — Rust port of bedtools summary.",
    origin: Some(Origin {
        upstream: "bedtools summary",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: None,
    }),
    usage_lines: &["-i <in.bed> -g <genome> [-o out.tsv]"],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[
            FlagSpec {
                short: Some('i'),
                long: "input",
                aliases: &["-i"],
                value: Some("<path>"),
                type_hint: Some("PathBuf"),
                required: true,
                default: None,
                description: "Input BED file.",
                why_default: None,
            },
            FlagSpec {
                short: Some('g'),
                long: "genome",
                aliases: &["-g"],
                value: Some("<path>"),
                type_hint: Some("PathBuf"),
                required: true,
                default: None,
                description: "Genome file: tab-separated chrom and size per line.",
                why_default: None,
            },
            FlagSpec {
                short: Some('o'),
                long: "output",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("String"),
                required: false,
                default: Some("-"),
                description: "Output file (default: stdout).",
                why_default: None,
            },
        ],
    }],
    examples: &[
        Example {
            description: "Summarise ChIP-seq peaks per chromosome",
            command: "rsomics-bed-summary -i peaks.bed -g hg38.genome",
        },
        Example {
            description: "Write summary to file",
            command: "rsomics-bed-summary -i intervals.bed -g chrom.sizes -o summary.tsv",
        },
    ],
    json_result_schema_doc: None,
};

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn cli_debug_assert() {
        Cli::command().debug_assert();
    }
}
