#![allow(clippy::cast_precision_loss)] // u64→f64 matches bedtools' float-division semantics

//! Statistical summary of BED intervals per chromosome — bedtools summary equivalent.
//!
//! Chromosomes with no intervals emit `-1` for min/max/mean (exact upstream behaviour).
//! Rows with intervals get a trailing `\t` after the mean column (upstream quirk replicated).

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

struct ChromStats {
    len: u64,
    count: u64,
    total_bp: u64,
    min: u64,
    max: u64,
}

impl ChromStats {
    fn new(len: u64) -> Self {
        Self {
            len,
            count: 0,
            total_bp: 0,
            min: u64::MAX,
            max: 0,
        }
    }

    fn add(&mut self, interval_len: u64) {
        self.count += 1;
        self.total_bp += interval_len;
        if interval_len < self.min {
            self.min = interval_len;
        }
        if interval_len > self.max {
            self.max = interval_len;
        }
    }
}

type GenomeTable = (Vec<(String, ChromStats)>, HashMap<String, usize>);

fn load_genome(path: &Path) -> Result<GenomeTable> {
    let file = File::open(path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?;
    let reader = BufReader::new(file);
    let mut order: Vec<(String, ChromStats)> = Vec::new();
    let mut idx: HashMap<String, usize> = HashMap::new();
    for line in reader.lines() {
        let line = line.map_err(RsomicsError::Io)?;
        let line = line.trim_end_matches('\r');
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut cols = line.splitn(2, '\t');
        let chrom = cols
            .next()
            .ok_or_else(|| RsomicsError::InvalidInput(format!("genome line: {line}")))?
            .to_string();
        let size_s = cols
            .next()
            .ok_or_else(|| RsomicsError::InvalidInput(format!("genome line missing size: {line}")))?
            .trim();
        let size: u64 = size_s.parse().map_err(|e| {
            RsomicsError::InvalidInput(format!("genome size parse '{size_s}': {e}"))
        })?;
        let i = order.len();
        idx.insert(chrom.clone(), i);
        order.push((chrom, ChromStats::new(size)));
    }
    Ok((order, idx))
}

pub fn summary(bed: &Path, genome: &Path, out: &mut dyn Write) -> Result<()> {
    let (mut stats, idx) = load_genome(genome)?;

    let bed_file = File::open(bed)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", bed.display())))?;
    let reader = BufReader::new(bed_file);
    for line in reader.lines() {
        let line = line.map_err(RsomicsError::Io)?;
        let line = line.trim_end_matches('\r');
        if line.is_empty()
            || line.starts_with('#')
            || line.starts_with("track")
            || line.starts_with("browser")
        {
            continue;
        }
        let mut cols = line.splitn(4, '\t');
        let chrom = cols.next().unwrap_or("");
        let start_s = cols.next().unwrap_or("");
        let end_s = cols.next().unwrap_or("");
        let start: u64 = start_s
            .parse()
            .map_err(|e| RsomicsError::InvalidInput(format!("start parse '{start_s}': {e}")))?;
        let end: u64 = end_s
            .parse()
            .map_err(|e| RsomicsError::InvalidInput(format!("end parse '{end_s}': {e}")))?;
        let i = idx.get(chrom).ok_or_else(|| {
            RsomicsError::InvalidInput(format!(
                "requested chromosome {chrom} does not exist in the genome file {}. Exiting.",
                genome.display()
            ))
        })?;
        stats[*i].1.add(end.saturating_sub(start));
    }

    let total_genome: u64 = stats.iter().map(|(_, s)| s.len).sum();
    let total_ivls: u64 = stats.iter().map(|(_, s)| s.count).sum();
    let total_bp: u64 = stats.iter().map(|(_, s)| s.total_bp).sum();
    let all_min: u64 = stats
        .iter()
        .filter(|(_, s)| s.count > 0)
        .map(|(_, s)| s.min)
        .min()
        .unwrap_or(0);
    let all_max: u64 = stats
        .iter()
        .filter(|(_, s)| s.count > 0)
        .map(|(_, s)| s.max)
        .max()
        .unwrap_or(0);
    let all_mean = if total_ivls == 0 {
        0.0f64
    } else {
        total_bp as f64 / total_ivls as f64
    };

    let mut out = BufWriter::with_capacity(256 * 1024, out);

    writeln!(
        out,
        "chrom\tchrom_length\tnum_ivls\ttotal_ivl_bp\tchrom_frac_genome\tfrac_all_ivls\tfrac_all_bp\tmin\tmax\tmean"
    )
    .map_err(RsomicsError::Io)?;

    for (chrom, s) in &stats {
        let chrom_frac = s.len as f64 / total_genome as f64;
        let frac_ivls = if total_ivls == 0 {
            0.0f64
        } else {
            s.count as f64 / total_ivls as f64
        };
        let frac_bp = if total_bp == 0 {
            0.0f64
        } else {
            s.total_bp as f64 / total_bp as f64
        };

        if s.count == 0 {
            writeln!(
                out,
                "{chrom}\t{}\t0\t0\t{chrom_frac:.9}\t{frac_ivls:.9}\t{frac_bp:.9}\t-1\t-1\t-1",
                s.len,
            )
            .map_err(RsomicsError::Io)?;
        } else {
            let mean = s.total_bp as f64 / s.count as f64;
            // Trailing tab after mean: upstream quirk, replicated byte-for-byte.
            writeln!(
                out,
                "{chrom}\t{}\t{}\t{}\t{chrom_frac:.9}\t{frac_ivls:.9}\t{frac_bp:.9}\t{}\t{}\t{mean:.9}\t",
                s.len, s.count, s.total_bp, s.min, s.max,
            )
            .map_err(RsomicsError::Io)?;
        }
    }

    writeln!(
        out,
        "all\t{total_genome}\t{total_ivls}\t{total_bp}\t1.0\t1.0\t1.0\t{all_min}\t{all_max}\t{all_mean:.9}",
    )
    .map_err(RsomicsError::Io)?;

    out.flush().map_err(RsomicsError::Io)
}
