//! Dependency-free extraction of the primitive routines from savont 0.8.0.
//! Source commit: 5a5a1ce9053a867abc09e991f61b8ec7fb4fc6ef.
//! Kept structurally identical so the Windows minimap2 linker is not involved.

use std::collections::HashSet;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Blockmer { pub kmer: u64, pub is_forward: bool }

#[inline]
fn enc(base: u8) -> u64 {
    match base { b'C' | b'c' => 1, b'G' | b'g' => 2, b'T' | b't' | b'U' | b'u' => 3, 1..=3 => base as u64, _ => 0 }
}

#[inline]
pub fn mm_hash64(kmer: u64) -> u64 {
    let mut key = kmer;
    key = (!key).wrapping_add(key << 21);
    key ^= key >> 24;
    key = (key.wrapping_add(key << 3)).wrapping_add(key << 8);
    key ^= key >> 14;
    key = (key.wrapping_add(key << 2)).wrapping_add(key << 4);
    key ^= key >> 28;
    key.wrapping_add(key << 31)
}

#[inline]
pub fn rev_hash_64(hashed_key: u64) -> u64 {
    let mut key = hashed_key;
    let mut tmp = key.wrapping_sub(key << 31); key = key.wrapping_sub(tmp << 31);
    tmp = key ^ key >> 28; key ^= tmp >> 28;
    key = key.wrapping_mul(14933078535860113213u64);
    tmp = key ^ key >> 14; tmp = key ^ tmp >> 14; tmp = key ^ tmp >> 14; key ^= tmp >> 14;
    key = key.wrapping_mul(15244667743933553977u64);
    tmp = key ^ key >> 24; key ^= tmp >> 24;
    tmp = !key; tmp = !(key.wrapping_sub(tmp << 21)); tmp = !(key.wrapping_sub(tmp << 21));
    !(key.wrapping_sub(tmp << 21))
}

fn position_min(slice: &[u64]) -> usize {
    slice.iter().enumerate().max_by(|(_, a), (_, b)| b.cmp(a)).unwrap().0
}

pub fn fmh_seeds(string: &[u8], output: &mut Vec<u64>, c: usize, k: usize) {
    if string.len() < k { return; }
    let mut forward = 0u64; let mut reverse = 0u64;
    let shift = 2 * (k - 1); let mask = u64::MAX >> (64 - 2 * k); let rev_mask = !(3 << (2 * k - 2));
    let threshold = u64::MAX / c as u64;
    for &base in string.iter().take(k - 1) {
        let f = enc(base); forward = (forward << 2) | f; reverse = (reverse >> 2) | ((3 - f) << shift);
    }
    for &base in string.iter().skip(k - 1) {
        let f = enc(base); forward = ((forward << 2) | f) & mask; reverse = ((reverse >> 2) & rev_mask) | ((3 - f) << shift);
        let hash = mm_hash64(forward.min(reverse)); if hash < threshold { output.push(hash); }
    }
}

pub fn fmh_seeds_positions(string: &[u8], output: &mut Vec<u64>, positions: &mut Vec<u64>, c: usize, k: usize) {
    if string.len() < k { return; }
    let mut forward = 0u64; let mut reverse = 0u64;
    let shift = 2 * (k - 1); let mask = u64::MAX >> (64 - 2 * k); let rev_mask = !(3 << (2 * k - 2));
    let threshold = u64::MAX / c as u64;
    for &base in string.iter().take(k - 1) { let f = enc(base); forward = (forward << 2) | f; reverse = (reverse >> 2) | ((3 - f) << shift); }
    for (i, &base) in string.iter().enumerate().skip(k - 1) {
        let f = enc(base); forward = ((forward << 2) | f) & mask; reverse = ((reverse >> 2) & rev_mask) | ((3 - f) << shift);
        let canonical = forward.min(reverse); if mm_hash64(canonical) < threshold { output.push(canonical); positions.push((i + 1 - k) as u64); }
    }
}

pub fn minimizer_seeds_positions(string: &[u8], output: &mut Vec<u64>, positions: &mut Vec<u64>, w: usize, k: usize) {
    if string.len() < k + w - 1 { return; }
    let mut forward = 0u64; let mut reverse = 0u64; let mut canonical = 0u64;
    let shift = 2 * (k - 1); let mask = u64::MAX >> (64 - 2 * k); let rev_mask = !(3 << (2 * k - 2));
    let mut window = vec![u64::MAX; w];
    for i in 0..k + w - 1 {
        let f = enc(string[i]); forward = (forward << 2) | f; reverse = (reverse >> 2) | ((3 - f) << shift);
        if i >= k - 1 { canonical = forward.min(reverse); window[i + 1 - k] = mm_hash64(canonical); }
    }
    let mut min_pos = position_min(&window); let mut min_val = window[min_pos];
    output.push(canonical); positions.push(min_pos as u64);
    for i in k + w - 1..string.len() {
        let f = enc(string[i]); forward = ((forward << 2) | f) & mask; reverse = ((reverse >> 2) & rev_mask) | ((3 - f) << shift);
        let hash = mm_hash64(forward.min(reverse)); let global = i + 1 - k; window[global % w] = hash;
        if hash < min_val { min_val = hash; min_pos = global % w; output.push(hash); positions.push(global as u64); }
        else if min_pos == global % w {
            min_pos = position_min(&window); min_val = window[min_pos];
            let offset = (global as i64 % w as i64 - min_pos as i64).rem_euclid(w as i64) as usize;
            output.push(min_val); positions.push((global - offset) as u64);
        }
    }
}

pub fn split_kmer_mid(string: Vec<u8>, qualities: Option<Vec<u8>>, k: usize, minimum_bq: u8) -> Vec<u64> {
    if string.len() < k { return vec![]; } assert!(k % 2 == 1 && k <= 31);
    let shift = 2 * (k - 1); let mask = u64::MAX >> (64 - 2 * k); let rev_mask = !(3 << (2 * k - 2)); let split_mask = !(3 << (k - 1));
    let mut skip = HashSet::new();
    if let Some(q) = qualities.as_ref() {
        if q.iter().skip(1).any(|x| x != &q[0]) { for i in k - 1..q.len() { let middle = i + 1 + k / 2 - k; if q[middle] - 33 < minimum_bq { skip.insert(i); } } }
    }
    let mut forward = 0u64; let mut reverse = 0u64; let mut result = Vec::new();
    for (i, &base) in string.iter().enumerate() {
        let f = enc(base); forward = ((forward << 2) | f) & mask; reverse = ((reverse >> 2) & rev_mask) | ((3 - f) << shift);
        if i < k - 1 { continue; } let sf = forward & split_mask; let sr = reverse & split_mask;
        if sf == sr || skip.contains(&i) { continue; } let is_forward = sf < sr; let canonical = if is_forward { forward } else { reverse };
        result.push(canonical | ((is_forward as u64) << 63));
    }
    result
}

pub fn blockmer_kmers(string: Vec<u8>, qualities: Option<Vec<u8>>, k: usize, l: usize, minimum_bq: u8) -> Vec<Blockmer> {
    if string.len() < k + l { return vec![]; } assert!(k <= 31 && k + l <= 31);
    let shift = 2 * (k - 1); let mask = u64::MAX >> (64 - 2 * k); let mut forward = 0u64; let mut reverse = 0u64; let mut result = Vec::new();
    for (i, &base) in string.iter().enumerate() {
        let f = enc(base); forward = ((forward << 2) | f) & mask; reverse = (reverse >> 2) | ((3 - f) << shift);
        if i < k - 1 || forward == reverse { continue; }
        if forward < reverse {
            if i + l >= string.len() { continue; } let mut suffix = 0; let mut low = false;
            for pos in i + 1..=i + l { suffix = (suffix << 2) | enc(string[pos]); if qualities.as_ref().is_some_and(|q| q[pos] - 33 < minimum_bq) { low = true; break; } }
            if !low { result.push(Blockmer { kmer: (forward << (2 * l)) | suffix, is_forward: true }); }
        } else {
            let start = i - k + 1; if start < l { continue; } let mut suffix = 0; let mut low = false;
            for j in 1..=l { let pos = start - j; suffix = (suffix << 2) | (3 - enc(string[pos])); if qualities.as_ref().is_some_and(|q| q[pos] - 33 < minimum_bq) { low = true; break; } }
            if !low { result.push(Blockmer { kmer: (reverse << (2 * l)) | suffix, is_forward: false }); }
        }
    }
    result
}

pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|b| match b { b'A' | b'a' => b'T', b'T' | b't' => b'A', b'C' | b'c' => b'G', b'G' | b'g' => b'C', b'N' | b'n' => b'N', _ => b'N' }).collect()
}

pub fn homopolymer_compress(seq: &[u8], do_hpc: bool) -> (Vec<u8>, Vec<u8>) {
    if seq.is_empty() { return (vec![], vec![]); } if !do_hpc { return (seq.to_vec(), vec![1; seq.len()]); }
    let mut bases = vec![]; let mut lengths = vec![]; let mut current = seq[0]; let mut run = 1u8;
    for &base in &seq[1..] { if base == current && run < 255 { run += 1; } else { bases.push(current); lengths.push(run); current = base; run = 1; } }
    bases.push(current); lengths.push(run); (bases, lengths)
}

pub fn homopolymer_compress_with_quality(seq: &[u8], qualities: &[u8], do_hpc: bool) -> (Vec<u8>, Vec<u8>, Vec<u8>) {
    if seq.is_empty() || seq.len() != qualities.len() { return (vec![], vec![], vec![]); } if !do_hpc { return (seq.to_vec(), qualities.to_vec(), vec![1; seq.len()]); }
    let mut bases = vec![]; let mut qs = vec![]; let mut lengths = vec![]; let mut current = seq[0]; let mut run = 1u8; let mut min_q = qualities[0];
    for i in 1..seq.len() { if seq[i] == current && run < 255 { run += 1; min_q = min_q.min(qualities[i]); } else { bases.push(current); qs.push(min_q); lengths.push(run); current = seq[i]; run = 1; min_q = qualities[i]; } }
    bases.push(current); qs.push(min_q); lengths.push(run); (bases, qs, lengths)
}

pub fn quality_pool(qualities: Vec<u8>) -> Vec<u8> {
    let mut pool = vec![]; for i in 0..qualities.len() { if i > 2 && i < qualities.len() - 2 { pool.push(*qualities[i - 2..i + 2].iter().min().unwrap()); } else { pool.push(qualities[i]); } } pool
}

pub fn estimate_sequence_identity(qualities: &[u8]) -> Option<f64> {
    if qualities.is_empty() { return None; } let sum: f64 = qualities.iter().map(|q| 10f64.powf(-((*q - 33) as f64) / 10.)).sum(); Some(100. - sum / qualities.len() as f64 * 100.)
}
