mod factoring;

use factoring::{is_prime, FactoringStatus, Factorizer as FactorizerImpl};
use num_bigint::BigUint;
use num_traits::{One, Zero};
use rand::thread_rng;
use rand::RngCore;
use std::collections::HashMap;
use wasm_bindgen::prelude::*;

// ---------------------------------------------------------------------------
// Miller-Rabin (probabilistic, for interactive primality display)
// ---------------------------------------------------------------------------

fn miller_rabin_test(n: &BigUint, a: &BigUint) -> bool {
    let n_minus_1 = n - BigUint::one();
    let mut d = n_minus_1.clone();
    let mut r = 0u32;
    while &d % 2u32 == BigUint::zero() {
        d >>= 1;
        r += 1;
    }
    let mut x = a.modpow(&d, n);
    if x == BigUint::one() || x == n_minus_1 {
        return true;
    }
    for _ in 0..r - 1 {
        x = x.modpow(&BigUint::from(2u32), n);
        if x == n_minus_1 {
            return true;
        }
    }
    false
}

pub enum PrimeStatus {
    Composite,
    ProbablyPrime(usize),
}

pub struct PrimeCheckerInner {
    results: HashMap<BigUint, PrimeStatus>,
}

impl PrimeCheckerInner {
    pub fn new() -> Self {
        PrimeCheckerInner {
            results: HashMap::new(),
        }
    }

    pub fn is_probably_prime(&self, n: &BigUint) -> bool {
        if n == &BigUint::from(2u32) || n == &BigUint::from(3u32) {
            return true;
        }
        if n < &BigUint::from(2u32) || n % 2u32 == BigUint::zero() {
            return false;
        }
        match self.results.get(n) {
            Some(PrimeStatus::ProbablyPrime(tests)) if *tests >= 10 => true,
            _ => false,
        }
    }

    /// Do one primality step. Returns (primality_tag, probability):
    ///   ("prime", 1.0)          – certain prime (deterministic for n < ~3e23; beyond that,
    ///                             additional random MR witnesses are accumulated per call)
    ///   ("composite", 0.0)      – definite composite
    ///   ("probably_prime", p)   – high-confidence but not yet certain (very large n only)
    ///   ("checking", 0.0)       – first call, no result yet (returned once before caching)
    pub fn step(&mut self, n: &BigUint) -> (&'static str, f64) {
        // Check cache first
        match self.results.get(n) {
            Some(PrimeStatus::Composite) => return ("composite", 0.0),
            Some(PrimeStatus::ProbablyPrime(tests)) => {
                let tests = *tests;
                if tests >= 20 {
                    return ("prime", 1.0);
                }
                let prob = 1.0 - 0.25_f64.powi(tests as i32);
                return ("probably_prime", prob);
            }
            None => {}
        }

        if n == &BigUint::from(2u32) || n == &BigUint::from(3u32) {
            self.results.insert(n.clone(), PrimeStatus::ProbablyPrime(20));
            return ("prime", 1.0);
        }
        if n < &BigUint::from(2u32) || n % 2u32 == BigUint::zero() {
            self.results.insert(n.clone(), PrimeStatus::Composite);
            return ("composite", 0.0);
        }

        // Deterministic Miller-Rabin with 12 fixed witnesses.
        // Sufficient for n < 3.18×10²³; for larger n we accumulate random witnesses.
        if n.bits() <= 77 {
            // n < 2^77 ≈ 1.5×10²³ — deterministic result
            if is_prime(n) {
                self.results.insert(n.clone(), PrimeStatus::ProbablyPrime(20));
                ("prime", 1.0)
            } else {
                self.results.insert(n.clone(), PrimeStatus::Composite);
                ("composite", 0.0)
            }
        } else {
            // Large n: one random MR witness per call, building confidence
            let mut rng = thread_rng();
            let mut a_bytes = vec![0u8; (n.bits() as usize + 7) / 8];
            rng.fill_bytes(&mut a_bytes);
            let a_raw = BigUint::from_bytes_be(&a_bytes);
            let range = n - BigUint::from(3u32);
            let a = (a_raw % &range) + BigUint::from(2u32);

            if !miller_rabin_test(n, &a) {
                self.results.insert(n.clone(), PrimeStatus::Composite);
                return ("composite", 0.0);
            }
            // Also run deterministic witnesses to boost confidence quickly
            for &w in &[2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
                let wb = BigUint::from(w);
                if &wb < n && !miller_rabin_test(n, &wb) {
                    self.results.insert(n.clone(), PrimeStatus::Composite);
                    return ("composite", 0.0);
                }
            }
            let new_tests = 13; // 12 deterministic + 1 random
            self.results.insert(n.clone(), PrimeStatus::ProbablyPrime(new_tests));
            let prob = 1.0 - 0.25_f64.powi(new_tests as i32);
            ("probably_prime", prob)
        }
    }
}

impl Default for PrimeCheckerInner {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// JSON builder helpers (no external serde dependency)
// ---------------------------------------------------------------------------

fn json_escape(s: &str) -> String {
    s.replace('\\', "\\\\")
        .replace('"', "\\\"")
        .replace('\n', "\\n")
}

/// Build the combined JSON response for a given number.
///
/// Schema:
/// {
///   "primality": "prime" | "composite" | "probably_prime" | "checking",
///   "probability": <float 0..1>,       // 1.0 for certain, else MR confidence
///   "factors_complete": <bool>,
///   "factors": [
///     {"value": "<decimal>", "exp": <int>, "is_prime": true},
///     ...
///   ],
/// }
///
/// Or on bad input:
/// {"error": "<message>"}
fn build_json(
    primality: &str,
    probability: f64,
    factoring: &FactoringStatus,
) -> String {
    let mut all_factors: Vec<String> = factoring
        .prime_factors
        .iter()
        .map(|(val, exp)| {
            format!(
                r#"{{"value":"{}","exp":{},"is_prime":true}}"#,
                json_escape(&val.to_string()),
                exp
            )
        })
        .collect();

    for val in &factoring.pending_values {
        all_factors.push(format!(
            r#"{{"value":"{}","exp":1,"is_prime":false}}"#,
            json_escape(&val.to_string())
        ));
    }

    format!(
        r#"{{"primality":"{}","probability":{:.15},"factors_complete":{},"factors":[{}]}}"#,
        primality,
        probability,
        factoring.complete,
        all_factors.join(","),
    )
}

fn error_json(msg: &str) -> String {
    format!(r#"{{"error":"{}"}}"#, json_escape(msg))
}

// ---------------------------------------------------------------------------
// Combined Checker – the primary WASM export
// ---------------------------------------------------------------------------

pub struct CheckerInner {
    prime: PrimeCheckerInner,
    factorizer: FactorizerImpl,
}

impl CheckerInner {
    pub fn new() -> Self {
        CheckerInner {
            prime: PrimeCheckerInner::new(),
            factorizer: FactorizerImpl::new(),
        }
    }

    /// Do one step of both primality checking and factoring for `n_str`.
    /// Returns a JSON string with all info needed by Elm for display.
    pub fn query(&mut self, n_str: &str) -> String {
        let n = match BigUint::parse_bytes(n_str.trim().as_bytes(), 10) {
            Some(n) => n,
            None => return error_json("Please enter a natural number"),
        };

        if n.is_zero() {
            return error_json("Please enter a natural number");
        }
        if n < BigUint::from(2u32) {
            return error_json("Input should be at least 2");
        }

        // One step of primality
        let (primality, probability) = self.prime.step(&n);

        // One step of factoring
        let factoring = match self.factorizer.step(n_str) {
            Ok(s) => s,
            Err(msg) => return error_json(&msg),
        };

        build_json(primality, probability, &factoring)
    }

    /// For next/prev prime search: quickly check if n is probably prime.
    pub fn is_probably_prime(&mut self, n_str: &str) -> bool {
        let n = match BigUint::parse_bytes(n_str.trim().as_bytes(), 10) {
            Some(n) => n,
            None => return false,
        };
        if n < BigUint::from(2u32) {
            return false;
        }
        // Run several checks and use the deterministic result
        for _ in 0..5 {
            self.prime.step(&n);
        }
        self.prime.is_probably_prime(&n)
    }
}

impl Default for CheckerInner {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// WASM-bindgen export
// ---------------------------------------------------------------------------

#[wasm_bindgen]
pub struct Checker(CheckerInner);

#[wasm_bindgen]
impl Checker {
    #[wasm_bindgen(constructor)]
    pub fn new() -> Self {
        #[cfg(feature = "console_error_panic_hook")]
        console_error_panic_hook::set_once();
        Checker(CheckerInner::new())
    }

    /// Do one step of primality + factoring. Returns JSON string.
    pub fn query(&mut self, n_str: String) -> String {
        self.0.query(&n_str)
    }

    /// Returns true if n is probably prime (for next/prev prime search).
    pub fn is_probably_prime(&mut self, n_str: String) -> bool {
        self.0.is_probably_prime(&n_str)
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::factoring::is_prime;

    // ---------------------------------------------------------------------------
    // Test helpers
    // ---------------------------------------------------------------------------

    /// Parse a JSON query result. Returns (primality, probability, factors_complete, factors).
    /// factors: Vec<(value_str, exp)> — includes both prime and composite (pending) factors.
    /// Panics if the result has an "error" key.
    fn parse_result(json: &str) -> (String, f64, bool, Vec<(String, usize)>) {
        assert!(
            !json.contains("\"error\""),
            "Got error JSON: {}",
            json
        );
        let primality = extract_str(json, "primality");
        let probability = extract_f64(json, "probability");
        let factors_complete = json.contains("\"factors_complete\":true");
        let factors = parse_factors_array(json);
        (primality, probability, factors_complete, factors)
    }

    fn extract_str(json: &str, key: &str) -> String {
        let search = format!("\"{}\":\"", key);
        let start = json.find(&search).expect(&format!("key {} not found in {}", key, json))
            + search.len();
        let end = json[start..].find('"').expect("closing quote") + start;
        json[start..end].to_string()
    }

    fn extract_f64(json: &str, key: &str) -> f64 {
        let search = format!("\"{}\":", key);
        let start = json.find(&search).expect(&format!("key {} not found", key)) + search.len();
        let rest = &json[start..];
        let end = rest.find(|c: char| !c.is_ascii_digit() && c != '.' && c != 'e' && c != '-' && c != '+')
            .unwrap_or(rest.len());
        rest[..end].parse().expect("f64")
    }

    fn parse_factors_array(json: &str) -> Vec<(String, usize)> {
        // Very simple: find "value":"..." and "exp":N pairs
        let mut result = Vec::new();
        let mut pos = 0;
        while let Some(val_idx) = json[pos..].find("\"value\":\"") {
            let val_start = pos + val_idx + 9;
            let val_end = json[val_start..].find('"').unwrap() + val_start;
            let value = json[val_start..val_end].to_string();

            let exp_search = "\"exp\":";
            let exp_idx = json[val_end..].find(exp_search).unwrap() + val_end + exp_search.len();
            let exp_end = json[exp_idx..].find(|c: char| !c.is_ascii_digit()).unwrap_or(json.len() - exp_idx) + exp_idx;
            let exp: usize = json[exp_idx..exp_end].parse().expect("exp");

            result.push((value, exp));
            pos = exp_end;
        }
        result
    }

    /// Compute product of factors from a parsed result
    fn factor_product(factors: &[(String, usize)]) -> BigUint {
        let mut product = BigUint::one();
        for (val, exp) in factors {
            let v: BigUint = val.parse().unwrap();
            for _ in 0..*exp {
                product *= &v;
            }
        }
        product
    }

    /// Run query until factoring is complete, return final JSON.
    fn query_to_completion(checker: &mut CheckerInner, n: &str, max_steps: usize) -> String {
        let mut result = String::new();
        for _ in 0..max_steps {
            result = checker.query(n);
            if result.contains("\"factors_complete\":true") {
                return result;
            }
        }
        result
    }

    /// Simulate Elm tick loop: returns (last_json, steps_to_complete)
    fn simulate_elm_ticks(n: &str, max_steps: usize) -> (String, usize) {
        let mut checker = CheckerInner::new();
        let mut result = String::new();
        for i in 0..max_steps {
            result = checker.query(n);
            if result.contains("\"factors_complete\":true") {
                return (result, i + 1);
            }
        }
        (result, max_steps)
    }

    // ---------------------------------------------------------------------------
    // JSON structure tests
    // ---------------------------------------------------------------------------

    #[test]
    fn test_json_structure_prime() {
        let mut c = CheckerInner::new();
        // Warm up with many checks
        for _ in 0..25 {
            c.query("7");
        }
        let json = c.query("7");
        let (primality, prob, complete, factors) = parse_result(&json);
        assert_eq!(primality, "prime");
        assert!((prob - 1.0).abs() < 1e-9);
        assert!(complete);
        assert_eq!(factors, vec![("7".to_string(), 1)]);
    }

    #[test]
    fn test_json_structure_composite() {
        let mut c = CheckerInner::new();
        let json = c.query("12");
        let (primality, _prob, _complete, _factors) = parse_result(&json);
        assert_eq!(primality, "composite");
    }

    #[test]
    fn test_json_error_non_numeric() {
        let mut c = CheckerInner::new();
        let json = c.query("abc");
        assert!(json.contains("\"error\""), "Should have error: {}", json);
    }

    #[test]
    fn test_json_error_too_small() {
        let mut c = CheckerInner::new();
        let json = c.query("1");
        assert!(json.contains("\"error\""), "Should have error: {}", json);
        let json = c.query("0");
        assert!(json.contains("\"error\""), "Should have error: {}", json);
    }

    // ---------------------------------------------------------------------------
    // Primality tests
    // ---------------------------------------------------------------------------

    #[test]
    fn test_primality_small_primes() {
        let primes = [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47];
        for &p in &primes {
            let mut c = CheckerInner::new();
            // Warm up enough to pass deterministic check
            for _ in 0..3 {
                c.query(&p.to_string());
            }
            let json = c.query(&p.to_string());
            let (primality, _, _, _) = parse_result(&json);
            assert!(
                primality == "prime" || primality == "probably_prime",
                "{} should be prime, got {}",
                p,
                primality
            );
        }
    }

    #[test]
    fn test_primality_composites() {
        let composites = [4u64, 6, 8, 9, 10, 12, 15, 25, 35, 49, 77, 100];
        for &n in &composites {
            let mut c = CheckerInner::new();
            let json = c.query(&n.to_string());
            let (primality, _, _, _) = parse_result(&json);
            assert_eq!(primality, "composite", "{} should be composite", n);
        }
    }

    #[test]
    fn test_primality_large_prime() {
        let mut c = CheckerInner::new();
        for _ in 0..25 {
            c.query("1000000007");
        }
        let json = c.query("1000000007");
        let (primality, _, _, _) = parse_result(&json);
        assert!(primality == "prime", "10^9+7 should be prime, got {}", primality);
    }

    #[test]
    fn test_miller_rabin_1_to_1000() {
        let mut checker = PrimeCheckerInner::new();
        let naive_is_prime = |n: u64| -> bool {
            if n < 2 { return false; }
            if n == 2 { return true; }
            if n % 2 == 0 { return false; }
            for i in (3..=(n as f64).sqrt() as u64 + 1).step_by(2) {
                if n % i == 0 { return false; }
            }
            true
        };
        for n in 2u64..=1000 {
            let nb = BigUint::from(n);
            // Run enough times to get definite answer
            for _ in 0..25 {
                checker.step(&nb);
            }
            let (tag, _) = checker.step(&nb);
            let expected = naive_is_prime(n);
            if expected {
                assert!(tag == "prime" || tag == "probably_prime", "n={} expected prime, got {}", n, tag);
            } else {
                assert_eq!(tag, "composite", "n={} expected composite, got {}", n, tag);
            }
        }
    }

    // ---------------------------------------------------------------------------
    // Factoring correctness tests
    // ---------------------------------------------------------------------------

    fn verify_factoring_json(n: u64, json: &str) {
        assert!(
            json.contains("\"factors_complete\":true"),
            "Factoring not complete for {}: {}",
            n,
            json
        );
        let (_, _, _, factors) = parse_result(json);
        let product = factor_product(&factors);
        assert_eq!(
            product,
            BigUint::from(n),
            "Factor product mismatch for {}: {:?}",
            n,
            factors
        );
        for (val, _) in &factors {
            let v: BigUint = val.parse().unwrap();
            assert!(is_prime(&v), "Factor {} of {} should be prime", val, n);
        }
    }

    #[test]
    fn test_factor_small_composites() {
        let cases: &[(u64, &[(u64, usize)])] = &[
            (4, &[(2, 2)]),
            (6, &[(2, 1), (3, 1)]),
            (8, &[(2, 3)]),
            (12, &[(2, 2), (3, 1)]),
            (30, &[(2, 1), (3, 1), (5, 1)]),
            (360, &[(2, 3), (3, 2), (5, 1)]),
        ];
        let mut checker = CheckerInner::new();
        for &(n, expected) in cases {
            let json = query_to_completion(&mut checker, &n.to_string(), 500);
            verify_factoring_json(n, &json);
            let (_, _, _, factors) = parse_result(&json);
            let expected_str: Vec<(String, usize)> = expected
                .iter()
                .map(|(v, e)| (v.to_string(), *e))
                .collect();
            assert_eq!(factors, expected_str, "n={}", n);
        }
    }

    #[test]
    fn test_factor_prime_inputs() {
        let mut checker = CheckerInner::new();
        for &p in &[2u64, 3, 5, 7, 97, 997, 9973, 99991] {
            let json = query_to_completion(&mut checker, &p.to_string(), 200);
            verify_factoring_json(p, &json);
            let (_, _, _, factors) = parse_result(&json);
            assert_eq!(factors, vec![(p.to_string(), 1)], "Prime {} factors to itself", p);
        }
    }

    #[test]
    fn test_factor_squares_medium_primes() {
        // Primes near sqrt(1e9) ≈ 31623 (the "flickering squares" scenario)
        let primes_near_sqrt_1e9: &[u64] = &[31607, 31627, 31643, 31657, 31663, 31667];
        let mut checker = CheckerInner::new();
        for &p in primes_near_sqrt_1e9 {
            let n = p * p;
            let json = query_to_completion(&mut checker, &n.to_string(), 300);
            verify_factoring_json(n, &json);
        }
    }

    #[test]
    fn test_factor_large_semiprime() {
        let p = 9_999_999_967u64;
        let q = 9_999_999_929u64;
        let n = BigUint::from(p) * BigUint::from(q);
        let mut checker = CheckerInner::new();
        let json = query_to_completion(&mut checker, &n.to_string(), 3000);
        assert!(
            json.contains("\"factors_complete\":true"),
            "Should factor large semiprime: {}",
            json
        );
        let (_, _, _, factors) = parse_result(&json);
        let product = factor_product(&factors);
        assert_eq!(product, n, "Factor product should equal original");
    }

    // -----------------------------------------------------------------------
    // Elm interaction simulation tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_elm_simulation_prime() {
        let (json, _steps) = simulate_elm_ticks("9999999967", 500);
        assert!(json.contains("\"factors_complete\":true"), "Should complete: {}", json);
        let (primality, _, _, factors) = parse_result(&json);
        assert!(primality == "prime" || primality == "probably_prime");
        assert_eq!(factors, vec![("9999999967".to_string(), 1)]);
    }

    #[test]
    fn test_elm_simulation_composite() {
        let (json, _steps) = simulate_elm_ticks("123456789", 500);
        assert!(json.contains("\"factors_complete\":true"));
        let (primality, _, _, factors) = parse_result(&json);
        assert_eq!(primality, "composite");
        let product = factor_product(&factors);
        assert_eq!(product, BigUint::from(123456789u64));
    }

    #[test]
    fn test_elm_no_flickering_for_squares() {
        let p = 31627u64;
        let n = p * p;
        let n_str = n.to_string();
        let mut checker = CheckerInner::new();

        // Run to completion
        let mut done_step = None;
        let mut results = Vec::new();
        for i in 0..500 {
            let r = checker.query(&n_str);
            results.push(r.clone());
            if r.contains("\"factors_complete\":true") && done_step.is_none() {
                done_step = Some(i);
            }
        }
        let ds = done_step.expect("Should complete factoring");
        // All results after completion should be identical
        for i in ds..results.len() {
            assert_eq!(results[i], results[ds], "Result changed at step {}", i);
        }
    }

    #[test]
    fn test_elm_caching_stability() {
        let mut checker = CheckerInner::new();
        let n = "12345678";
        let json = query_to_completion(&mut checker, n, 1000);
        for _ in 0..10 {
            let r = checker.query(n);
            assert_eq!(r, json, "Cached result should be stable");
        }
    }

    #[test]
    fn test_elm_number_switch() {
        let mut checker = CheckerInner::new();
        // Work on 35 briefly
        for _ in 0..5 {
            checker.query("35");
        }
        // Switch to 77
        let json_77 = query_to_completion(&mut checker, "77", 200);
        verify_factoring_json(77, &json_77);
        // 35 still factorable (independent cache entries)
        let json_35 = query_to_completion(&mut checker, "35", 200);
        verify_factoring_json(35, &json_35);
    }

    #[test]
    fn test_elm_switch_back_resumes() {
        let mut checker = CheckerInner::new();
        // Partial work on large number
        for _ in 0..10 {
            checker.query("999999937");
        }
        // Complete a small number
        let _ = query_to_completion(&mut checker, "30", 100);
        // Switch back — should resume
        let json = query_to_completion(&mut checker, "999999937", 500);
        assert!(json.contains("\"factors_complete\":true"), "Should complete: {}", json);
    }

    // -----------------------------------------------------------------------
    // Factoring edge cases and squares
    // -----------------------------------------------------------------------

    #[test]
    fn test_factor_squares_large_primes() {
        let primes: &[u64] = &[1_000_003, 1_000_033, 9_999_991];
        let mut checker = CheckerInner::new();
        for &p in primes {
            let n = p * p;
            let json = query_to_completion(&mut checker, &n.to_string(), 500);
            verify_factoring_json(n, &json);
        }
    }

    #[test]
    fn test_factor_many_small_factors() {
        // 2^10 * 3^5 * 5^3 * 7
        let n: u64 = 2u64.pow(10) * 3u64.pow(5) * 5u64.pow(3) * 7;
        let mut checker = CheckerInner::new();
        let json = query_to_completion(&mut checker, &n.to_string(), 500);
        verify_factoring_json(n, &json);
        let (_, _, _, factors) = parse_result(&json);
        let two_exp = factors.iter().find(|(v, _)| v == "2").map(|(_, e)| *e).unwrap_or(0);
        assert_eq!(two_exp, 10);
    }

    #[test]
    fn test_factor_step_bounded_time() {
        let mut checker = CheckerInner::new();
        let start = std::time::Instant::now();
        for _ in 0..100 {
            checker.query("123456789");
        }
        let elapsed = start.elapsed();
        assert!(elapsed.as_secs() < 30, "100 steps took too long: {:?}", elapsed);
    }

    // -----------------------------------------------------------------------
    // ECM-specific tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_ecm_finds_factor_of_semiprimes() {
        use crate::factoring::{sieve_primes, try_ecm_curve};

        let primes = sieve_primes(200_000);
        // Known semiprime: 1299709 * 1299827 (both prime)
        let p = 1_299_709u64;
        let q = 1_299_827u64;
        if is_prime(&BigUint::from(p)) && is_prime(&BigUint::from(q)) {
            let n = BigUint::from(p) * BigUint::from(q);
            let mut found = false;
            for sigma in 6..200 {
                if let Some(f) = try_ecm_curve(&n, sigma, 5000, 100_000, &primes) {
                    assert!(
                        f == BigUint::from(p) || f == BigUint::from(q),
                        "Wrong factor: {}",
                        f
                    );
                    found = true;
                    break;
                }
            }
            assert!(found, "ECM should find a factor");
        }
    }

    #[test]
    fn test_is_probably_prime_method() {
        let mut checker = CheckerInner::new();
        assert!(checker.is_probably_prime("7"));
        assert!(checker.is_probably_prime("997"));
        assert!(!checker.is_probably_prime("9"));
        assert!(!checker.is_probably_prime("100"));
    }

    // ---------------------------------------------------------------------------
    // Responsiveness tests
    //
    // Elm calls query() at most once every 100 ms.  A single query() call that
    // blocks for longer than ~100 ms will freeze the browser.  These tests
    // simulate the Elm interaction model correctly: they measure the wall-clock
    // time of INDIVIDUAL query() calls (not the total time to complete factoring).
    // ---------------------------------------------------------------------------

    /// Time a single query() call.  Panics if it exceeds `budget_ms`.
    fn assert_single_call_fast(n: &str, budget_ms: u64) {
        use std::time::Instant;
        let mut checker = CheckerInner::new();
        // The very first call may do trial division (fast) or start Pollard/ECM.
        // Check a few calls to catch slow ones mid-factoring as well.
        for _ in 0..5 {
            let t = Instant::now();
            checker.query(n);
            let ms = t.elapsed().as_millis() as u64;
            assert!(
                ms <= budget_ms,
                "query({}) took {}ms, budget {}ms — UI would freeze",
                n, ms, budget_ms
            );
        }
    }

    #[test]
    fn test_responsiveness_large_with_small_factors() {
        // = 7 × 43 × 25301 × 2964807793319969 × 546782985124176271
        // Previously froze browser for 10+ s due to unbounded Pollard-Brent.
        assert_single_call_fast("12345699943999999909990999999999999999999", 200);
    }

    #[test]
    fn test_responsiveness_repunit_39() {
        // 39-digit repunit; previously caused "page unresponsive" warning.
        assert_single_call_fast("111111111111111111111111111111111111111", 200);
    }

    #[test]
    #[cfg(not(debug_assertions))]
    fn test_responsiveness_release_tighter_budget() {
        // In release mode the budget is tighter — should be well under 50ms.
        use std::time::Instant;
        let cases = [
            "12345699943999999909990999999999999999999",
            "111111111111111111111111111111111111111",
        ];
        let mut checker = CheckerInner::new();
        for n in &cases {
            for _ in 0..3 {
                let t = Instant::now();
                checker.query(n);
                let ms = t.elapsed().as_millis() as u64;
                assert!(ms <= 50, "query({}) took {}ms in release mode", n, ms);
            }
        }
    }
}
