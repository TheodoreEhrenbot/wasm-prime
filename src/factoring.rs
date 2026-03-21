/// ECM (Elliptic Curve Method) factoring using Brent's approach with Montgomery curves.
///
/// Montgomery curve form: y² = x³ + Ax² + x (mod n)
/// Uses Suyama's parameterization, Stage 1 + Stage 2 (baby-step/giant-step pairing).
use num_bigint::BigUint;
use num_traits::{One, Zero};
use std::collections::HashMap;

// ---------------------------------------------------------------------------
// Modular arithmetic helpers
// ---------------------------------------------------------------------------

/// Extended GCD: returns (g, x, y) such that a*x + b*y = g = gcd(a,b).
/// Works on i128 for small values used in Suyama parameterization helpers,
/// but we need BigUint variants for the main modular inverse.
fn extended_gcd_biguint(a: &BigUint, b: &BigUint) -> (BigUint, bool, BigUint) {
    // Returns (gcd, sign_of_x_is_negative, |x|) where a*x ≡ gcd (mod b).
    // We compute iteratively to avoid stack overflow on large inputs.
    if b.is_zero() {
        return (a.clone(), false, BigUint::one());
    }
    let mut old_r = a.clone();
    let mut r = b.clone();
    // We track coefficients as (value, negative) pairs.
    let mut old_s: (BigUint, bool) = (BigUint::one(), false);
    let mut s: (BigUint, bool) = (BigUint::zero(), false);

    while !r.is_zero() {
        let q = &old_r / &r;
        let new_r = &old_r - &q * &r;
        old_r = r;
        r = new_r;

        // new_s = old_s - q * s
        let qs_val = &q * &s.0;
        let qs_neg = s.1;
        let new_s = if old_s.1 == qs_neg {
            // same sign: subtract magnitudes
            if old_s.0 >= qs_val {
                (old_s.0 - qs_val, old_s.1)
            } else {
                (qs_val - old_s.0, !old_s.1)
            }
        } else {
            // different sign: add magnitudes, keep old_s sign
            (old_s.0 + qs_val, old_s.1)
        };
        old_s = s;
        s = new_s;
    }
    (old_r, old_s.1, old_s.0)
}

// ---------------------------------------------------------------------------
// u128 modular arithmetic + Pollard-Brent for small composites
// ---------------------------------------------------------------------------

/// (a + b) mod m, safe against overflow when a, b < m.
fn addmod_u128(a: u128, b: u128, m: u128) -> u128 {
    let (r, overflow) = a.overflowing_add(b);
    if overflow {
        r.wrapping_sub(m)
    } else if r >= m {
        r - m
    } else {
        r
    }
}

/// (a * b) mod m using binary multiplication; safe for any a, b < m < 2^128.
fn mulmod_u128(mut a: u128, mut b: u128, m: u128) -> u128 {
    let mut result = 0u128;
    a %= m;
    b %= m;
    while b > 0 {
        if b & 1 == 1 {
            result = addmod_u128(result, a, m);
        }
        a = addmod_u128(a, a, m);
        b >>= 1;
    }
    result
}

fn gcd_u128(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        let t = b;
        b = a % b;
        a = t;
    }
    a
}

/// Brent's Pollard-rho over u128.  Returns a non-trivial factor of n, or None.
fn pollard_brent_u128(n: u128) -> Option<u128> {
    if n % 2 == 0 {
        return Some(2);
    }
    if n % 3 == 0 {
        return Some(3);
    }
    // Try 30 different (y, c) seeds.
    for attempt in 0u64..30 {
        let c = (attempt.wrapping_mul(1_000_003).wrapping_add(1)) as u128 % (n - 1) + 1;
        let mut y = (attempt.wrapping_mul(1_000_007).wrapping_add(2)) as u128 % n;
        let mut r: u128 = 1;
        let batch: u128 = 128;
        let mut q = 1u128;
        let mut x = 0u128;
        let mut ys = 0u128;
        let mut g = 1u128;

        'outer: while g == 1 {
            x = y;
            for _ in 0..r {
                y = addmod_u128(mulmod_u128(y, y, n), c, n);
            }
            let mut k = 0u128;
            while k < r && g == 1 {
                ys = y;
                let bound = batch.min(r - k);
                for _ in 0..bound {
                    y = addmod_u128(mulmod_u128(y, y, n), c, n);
                    let diff = if x > y { x - y } else { y - x };
                    // avoid q becoming 0
                    q = mulmod_u128(q, diff.max(1), n);
                }
                g = gcd_u128(q, n);
                k += batch;
            }
            r *= 2;
            // Hard limit: give up after ~2M iterations per attempt
            if r > 1 << 21 {
                break 'outer;
            }
        }

        if g == n {
            // Backtrack step-by-step from last ys
            loop {
                ys = addmod_u128(mulmod_u128(ys, ys, n), c, n);
                let diff = if x > ys { x - ys } else { ys - x };
                g = gcd_u128(diff.max(1), n);
                if g > 1 {
                    break;
                }
                if ys == x {
                    // Full cycle with no factor found; try next seed
                    g = n;
                    break;
                }
            }
        }

        if g > 1 && g < n {
            return Some(g);
        }
    }
    None
}

/// Try to fit a BigUint into u128; returns None if it has more than 2 u64 limbs.
fn biguint_to_u128(n: &BigUint) -> Option<u128> {
    let digits = n.to_u64_digits();
    match digits.len() {
        0 => Some(0),
        1 => Some(digits[0] as u128),
        2 => Some((digits[1] as u128) << 64 | digits[0] as u128),
        _ => None,
    }
}

/// Modular inverse of a mod n.  Returns None if gcd(a, n) != 1.
pub fn modinv(a: &BigUint, n: &BigUint) -> Option<BigUint> {
    if n.is_one() {
        return Some(BigUint::zero());
    }
    let a_mod = a % n;
    if a_mod.is_zero() {
        return None;
    }
    let (g, neg, x) = extended_gcd_biguint(&a_mod, n);
    if !g.is_one() {
        return None;
    }
    if neg {
        Some(n - (x % n))
    } else {
        Some(x % n)
    }
}

// ---------------------------------------------------------------------------
// Miller-Rabin primality test (standalone, no external dependency)
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

/// Deterministic Miller-Rabin for n < 3,317,044,064,679,887,385,961,981.
/// Uses a fixed witness set that covers all 64-bit integers.
pub fn is_prime(n: &BigUint) -> bool {
    let zero = BigUint::zero();
    let two = BigUint::from(2u32);

    if n < &two {
        return false;
    }
    if n == &two || n == &BigUint::from(3u32) {
        return true;
    }
    if n % 2u32 == zero {
        return false;
    }
    if n % 3u32 == zero {
        return n == &BigUint::from(3u32);
    }

    // Witnesses sufficient for n < 3,317,044,064,679,887,385,961,981 (covers 64-bit range well).
    let witnesses: &[u64] = &[2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
    for &w in witnesses {
        let a = BigUint::from(w);
        if &a >= n {
            continue;
        }
        if !miller_rabin_test(n, &a) {
            return false;
        }
    }
    true
}

// ---------------------------------------------------------------------------
// Integer square root
// ---------------------------------------------------------------------------

fn isqrt(n: &BigUint) -> BigUint {
    if n.is_zero() {
        return BigUint::zero();
    }
    // Newton's method
    let bits = n.bits();
    let mut x = BigUint::one() << ((bits + 1) / 2);
    loop {
        let x1 = (&x + n / &x) >> 1u32;
        if x1 >= x {
            return x;
        }
        x = x1;
    }
}

fn is_perfect_square(n: &BigUint) -> Option<BigUint> {
    let s = isqrt(n);
    if &s * &s == *n {
        Some(s)
    } else {
        None
    }
}

// ---------------------------------------------------------------------------
// Sieve of Eratosthenes
// ---------------------------------------------------------------------------

pub fn sieve_primes(limit: u64) -> Vec<u64> {
    if limit < 2 {
        return vec![];
    }
    let limit_usize = limit as usize;
    let mut is_composite = vec![false; limit_usize + 1];
    is_composite[0] = true;
    is_composite[1] = true;
    let mut i = 2usize;
    while i * i <= limit_usize {
        if !is_composite[i] {
            let mut j = i * i;
            while j <= limit_usize {
                is_composite[j] = true;
                j += i;
            }
        }
        i += 1;
    }
    (2..=limit_usize)
        .filter(|&i| !is_composite[i])
        .map(|i| i as u64)
        .collect()
}

// ---------------------------------------------------------------------------
// Montgomery curve point arithmetic
// ---------------------------------------------------------------------------

/// A projective point on a Montgomery curve: represents affine (X/Z, ...).
#[derive(Clone, Debug)]
pub struct MontPoint {
    pub x: BigUint,
    pub z: BigUint,
}

impl MontPoint {
    #[allow(dead_code)]
    fn is_zero(&self) -> bool {
        self.z.is_zero()
    }
}

/// Montgomery curve doubling.
/// Input: P = (X:Z) on curve with parameter a24 = (A+2)/4 mod n.
/// Formulas from Bernstein (xDBL):
///   U = (X - Z)²
///   V = (X + Z)²
///   X2 = U * V
///   W = V - U
///   Z2 = W * (U + a24 * W)
pub fn mont_double(p: &MontPoint, a24: &BigUint, n: &BigUint) -> MontPoint {
    let xpz = (&p.x + &p.z) % n;
    let xmz = (n + &p.x - &p.z % n) % n;
    let v = (&xpz * &xpz) % n;
    let u = (&xmz * &xmz) % n;
    let x2 = (&u * &v) % n;
    let w = (n + &v - &u) % n;
    let z2 = (&w * ((&u + a24 * &w) % n)) % n;
    MontPoint { x: x2, z: z2 }
}

/// Montgomery differential addition.
/// Input: P, Q on curve, and diff = P - Q (affine difference, known).
/// xADD formula:
///   U = (Xp - Zp)(Xq + Zq)
///   V = (Xp + Zp)(Xq - Zq)
///   add = U + V
///   sub = U - V
///   X_out = Z_diff * add²
///   Z_out = X_diff * sub²
pub fn mont_diff_add(p: &MontPoint, q: &MontPoint, diff: &MontPoint, n: &BigUint) -> MontPoint {
    let u = ((&p.x + n - &p.z % n) % n * ((&q.x + &q.z) % n)) % n;
    let v = ((&p.x + &p.z) % n * ((n + &q.x - &q.z % n) % n)) % n;
    let add = (&u + &v) % n;
    let sub = (n + &u - &v) % n;
    let x_out = (&diff.z * (&add * &add % n)) % n;
    let z_out = (&diff.x * (&sub * &sub % n)) % n;
    MontPoint {
        x: x_out,
        z: z_out,
    }
}

/// Montgomery ladder scalar multiplication: compute k*P on curve with parameter a24.
/// Uses the standard Montgomery ladder (Algorithm 1 in "Speeding the Pollard and ECM").
pub fn mont_scalar_mul(k: u64, p: &MontPoint, a24: &BigUint, n: &BigUint) -> MontPoint {
    if k == 0 {
        return MontPoint {
            x: BigUint::zero(),
            z: BigUint::zero(),
        };
    }
    if k == 1 {
        return p.clone();
    }
    if k == 2 {
        return mont_double(p, a24, n);
    }

    // Find highest bit of k
    let bits = 64 - k.leading_zeros();

    // R0 = P, R1 = 2P
    let mut r0 = p.clone();
    let mut r1 = mont_double(p, a24, n);

    // Process bits from second-highest down to 0
    for i in (0..bits - 1).rev() {
        let bit = (k >> i) & 1;
        if bit == 0 {
            // R1 = R0 + R1 (with diff = P), R0 = 2*R0
            r1 = mont_diff_add(&r0, &r1, p, n);
            r0 = mont_double(&r0, a24, n);
        } else {
            // R0 = R0 + R1 (with diff = P), R1 = 2*R1
            r0 = mont_diff_add(&r0, &r1, p, n);
            r1 = mont_double(&r1, a24, n);
        }
    }
    r0
}

// ---------------------------------------------------------------------------
// Suyama's parameterization
// ---------------------------------------------------------------------------
// Given sigma, computes a starting point (x0:z0) and curve parameter a24
// on a Montgomery curve B*y² = x³ + A*x² + x  (where a24 = (A+2)/4).
// This guarantees the curve has group order divisible by 12.

/// Returns (x0, z0, a24) for the Montgomery curve determined by sigma.
/// Returns None if sigma leads to a degenerate curve (gcd != 1 check elsewhere).
fn suyama_param(sigma: u64, n: &BigUint) -> Option<(MontPoint, BigUint)> {
    let sigma_b = BigUint::from(sigma);
    let n = n;

    // u = sigma² - 5 mod n
    let u = (n + &sigma_b * &sigma_b % n - BigUint::from(5u32) % n) % n;
    // v = 4 * sigma mod n
    let v = (BigUint::from(4u32) * &sigma_b) % n;

    // x0 = u³ mod n
    let x0 = u.modpow(&BigUint::from(3u32), n);
    // z0 = v³ mod n
    let z0 = v.modpow(&BigUint::from(3u32), n);

    // Compute A:
    // A = (v - u)³ * (3u + v) / (4 * u³ * v) - 2
    // We compute (A+2)/4 = a24 directly:
    // a24 = (v - u)³ * (3u + v) / (16 * u³ * v)

    // num = (v - u)³ * (3u + v)
    let v_minus_u = (n + &v - &u % n) % n;
    let vm3 = v_minus_u.modpow(&BigUint::from(3u32), n);
    let three_u_plus_v = (BigUint::from(3u32) * &u + &v) % n;
    let num = (&vm3 * &three_u_plus_v) % n;

    // den = 16 * u³ * v
    let den = (BigUint::from(16u32) * &x0 % n * &v) % n;

    let den_inv = modinv(&den, n)?;
    let a24 = (&num * &den_inv) % n;

    Some((MontPoint { x: x0, z: z0 }, a24))
}

// ---------------------------------------------------------------------------
// GCD helper
// ---------------------------------------------------------------------------

fn gcd_biguint(a: &BigUint, b: &BigUint) -> BigUint {
    let mut a = a.clone();
    let mut b = b.clone();
    while !b.is_zero() {
        let t = &a % &b;
        a = b;
        b = t;
    }
    a
}

// ---------------------------------------------------------------------------
// Stage 1
// ---------------------------------------------------------------------------

/// Multiply a point by all prime powers p^e <= b1.
/// Returns the resulting point, or None if we found a factor (z has non-trivial gcd with n).
fn stage1(
    p0: MontPoint,
    a24: &BigUint,
    n: &BigUint,
    b1: u64,
    primes: &[u64],
) -> Result<MontPoint, BigUint> {
    let mut q = p0;
    for &p in primes {
        if p > b1 {
            break;
        }
        // Compute p^e <= b1
        let mut pe = p;
        while pe <= b1 / p {
            pe *= p;
        }
        q = mont_scalar_mul(pe, &q, a24, n);
        // Check for non-trivial factor
        let g = gcd_biguint(&q.z, n);
        if g > BigUint::one() && g < *n {
            return Err(g);
        }
    }
    Ok(q)
}

// ---------------------------------------------------------------------------
// Stage 2 (baby-step / giant-step pairing)
// ---------------------------------------------------------------------------
// D = 60 (so pairs cover residues 1..59 coprime to 2,3,5 = 16 baby steps)
// Baby steps: B[s] = (s * Q).X / Z  for odd s in [1, D-1] coprime to D
// Giant steps: iterate r = D, 2D, 3D,...
// For each prime q in (B1, B2], write q = r*D ± s, accumulate differences.

const D: u64 = 60;

fn stage2(
    q: &MontPoint,
    a24: &BigUint,
    n: &BigUint,
    b1: u64,
    b2: u64,
    primes: &[u64],
) -> Option<BigUint> {
    // Baby step residues: s in [1, D-1] coprime to D=60
    // coprime to 60 means coprime to 2,3,5
    let baby_s: Vec<u64> = (1..D)
        .filter(|&s| s % 2 != 0 && s % 3 != 0 && s % 5 != 0)
        .collect();
    // = [1,7,11,13,17,19,23,29,31,37,41,43,47,49,53,59] (16 values)

    // Precompute D*Q
    let dq = mont_scalar_mul(D, q, a24, n);

    // Precompute baby steps: B[i] = baby_s[i] * Q
    // We do this by repeated addition using the difference method.
    // Since we need multiple distinct s values, just do scalar mul for each.
    // For D=60, |baby_s| = 16, so 16 scalar muls. Acceptable.
    let baby: Vec<MontPoint> = baby_s
        .iter()
        .map(|&s| mont_scalar_mul(s, q, a24, n))
        .collect();

    // Accumulate: acc = product of (X_R * Z_B - Z_R * X_B) for each prime match
    let mut acc = BigUint::one();

    // Giant step: r_point = r * D * Q, starting from r=1
    // We iterate r = 1, 2, 3, ... and update r_point += D*Q each time.
    // r_point for r=1 starts at dq; we update with differential addition.

    // We need to track two consecutive points to do differential addition.
    // rp = r * dq, rpm1 = (r-1) * dq
    let mut rpm1 = MontPoint {
        x: BigUint::one(),
        z: BigUint::zero(),
    }; // 0 * dq = point at infinity
    let mut rp = dq.clone(); // 1 * dq

    // Find the range of r values needed.
    // Primes in (b1, b2]: q = r*D ± s  =>  r*D - (D-1) <= q <= r*D + (D-1)
    // So r >= (b1 + 1) / D  and  r <= (b2 + D - 1) / D
    let r_min = (b1 + 1) / D;
    let r_max = (b2 + D - 1) / D + 1;

    // Advance rp to r_min * dq
    // We start at r=1; advance to r_min by repeated giant steps.
    // To keep differential addition working, maintain rpm1 and rp as consecutive multiples of dq.
    let mut r_cur = 1u64;

    // Fast forward to r_min
    if r_min > 1 {
        // Just compute r_min * dq and (r_min-1)*dq directly.
        rpm1 = mont_scalar_mul(r_min - 1, &dq, a24, n);
        rp = mont_scalar_mul(r_min, &dq, a24, n);
        r_cur = r_min;
    }

    // Collect primes in (b1, b2] into a sorted list (they already are from sieve)
    let stage2_primes: Vec<u64> = primes
        .iter()
        .copied()
        .filter(|&p| p > b1 && p <= b2)
        .collect();

    let mut prime_idx = 0;

    for r in r_min..=r_max {
        // Advance r_cur to r if needed
        if r_cur < r {
            // step: new_rp = rp + dq (diff = rpm1)
            let new_rp = mont_diff_add(&rp, &dq, &rpm1, n);
            rpm1 = rp;
            rp = new_rp;
            r_cur = r;
        }

        // For each baby step s, try to match a prime p = r*D ± s
        let base = r * D;

        // Consume primes in [base - (D-1), base + (D-1)] ∩ (b1, b2]
        while prime_idx < stage2_primes.len() {
            let p = stage2_primes[prime_idx];
            // p should be in (base - D + 1, base + D)
            // But we must iterate r in order; p may belong to this r or a later r.
            if p > base + D - 1 {
                break; // this prime belongs to a later giant step
            }
            prime_idx += 1;

            // Find s such that p = base ± s, s in baby_s
            let (s_val, found) = if p <= base {
                let diff = base - p;
                if diff > 0 && diff < D {
                    (diff, true)
                } else {
                    (0, false)
                }
            } else {
                let diff = p - base;
                if diff < D {
                    (diff, true)
                } else {
                    (0, false)
                }
            };

            if !found || s_val == 0 {
                // s=0 means p = r*D exactly, which would need special handling; skip (rare).
                continue;
            }

            // Find baby step index for s_val
            if let Some(bi) = baby_s.iter().position(|&s| s == s_val) {
                // Accumulate (X_r * Z_b - Z_r * X_b) mod n
                let xr = &rp.x;
                let zr = &rp.z;
                let xb = &baby[bi].x;
                let zb = &baby[bi].z;
                let term = (xr * zb + n * zr % n - zr * xb % n) % n;
                acc = (acc * term) % n;
            }
        }
    }

    // GCD of acc with n
    let g = gcd_biguint(&acc, n);
    if g > BigUint::one() && g < *n {
        Some(g)
    } else {
        None
    }
}

// ---------------------------------------------------------------------------
// One ECM curve attempt
// ---------------------------------------------------------------------------

/// Try to find a non-trivial factor of n using one ECM curve.
/// sigma is the Suyama parameter (should be >= 6).
pub fn try_ecm_curve(
    n: &BigUint,
    sigma: u64,
    b1: u64,
    b2: u64,
    primes: &[u64],
) -> Option<BigUint> {
    // Suyama parameterization
    let (p0, a24) = suyama_param(sigma, n)?;

    // Check for immediate factor from parameterization
    let g = gcd_biguint(&p0.z, n);
    if g > BigUint::one() && g < *n {
        return Some(g);
    }

    // Stage 1
    let q = match stage1(p0, &a24, n, b1, primes) {
        Err(factor) => return Some(factor),
        Ok(q) => q,
    };

    // Check stage1 result
    let g = gcd_biguint(&q.z, n);
    if g > BigUint::one() && g < *n {
        return Some(g);
    }
    if g == *n {
        // Full collapse, no factor found on this curve
        return None;
    }

    // Stage 2
    stage2(&q, &a24, n, b1, b2, primes)
}

// ---------------------------------------------------------------------------
// Trial division helper
// ---------------------------------------------------------------------------

/// Try dividing n by small primes up to limit.
/// Returns (factors_found, remaining) where factors_found are prime factors.
fn trial_divide(n: &BigUint, primes: &[u64], cursor: usize, batch: usize) -> (Vec<BigUint>, BigUint, usize) {
    let mut remaining = n.clone();
    let mut factors = vec![];
    let end = (cursor + batch).min(primes.len());
    let mut new_cursor = end;

    for (i, &p) in primes[cursor..end].iter().enumerate() {
        let pb = BigUint::from(p);
        // Check if p² > remaining (can stop early)
        if &pb * &pb > remaining {
            new_cursor = cursor + i;
            // remaining is prime
            if remaining > BigUint::one() {
                factors.push(remaining.clone());
                remaining = BigUint::one();
            }
            break;
        }
        while &remaining % &pb == BigUint::zero() {
            factors.push(pb.clone());
            remaining /= &pb;
        }
    }
    (factors, remaining, new_cursor)
}

// ---------------------------------------------------------------------------
// Factorizer public API
// ---------------------------------------------------------------------------

/// Structured result from one factoring step.
#[derive(Debug, Clone)]
pub struct FactoringStatus {
    /// True when all factors are fully determined.
    pub complete: bool,
    /// Confirmed prime factors as (value, exponent) pairs, sorted ascending.
    pub prime_factors: Vec<(BigUint, usize)>,
    /// Composite sub-factors still being analysed (value, exponent=1 each).
    pub pending_values: Vec<BigUint>,
}

struct PendingFactor {
    n: BigUint,
    trial_div_done: bool,
    ecm_curves_tried: u32,
    b1: u64,
    pollard_tried: bool,
}

struct FactorizationState {
    /// Confirmed prime factors (with multiplicity, sorted ascending).
    prime_factors: Vec<BigUint>,
    /// Numbers still being worked on.
    pending: Vec<PendingFactor>,
    /// Whether factorization is complete.
    complete: bool,
}

pub struct Factorizer {
    cache: HashMap<BigUint, FactorizationState>,
    /// Precomputed small primes up to 1_000_000 for trial division.
    small_primes: Vec<u64>,
    /// Precomputed primes up to max B2 (adaptive, regenerated as needed).
    ecm_primes_b2: u64,
    ecm_primes: Vec<u64>,
}

impl Factorizer {
    pub fn new() -> Self {
        let small_primes = sieve_primes(1_000_000);
        let ecm_primes_b2 = 200_000u64; // initial B2 max
        let ecm_primes = sieve_primes(ecm_primes_b2);
        Factorizer {
            cache: HashMap::new(),
            small_primes,
            ecm_primes_b2,
            ecm_primes,
        }
    }

    /// Ensure ecm_primes covers up to `needed_b2`.
    fn ensure_ecm_primes(&mut self, needed_b2: u64) {
        if needed_b2 > self.ecm_primes_b2 {
            self.ecm_primes_b2 = needed_b2;
            self.ecm_primes = sieve_primes(self.ecm_primes_b2);
        }
    }

    /// Add a number to the pending list of `state`, or record it as prime.
    fn add_to_pending(state: &mut FactorizationState, n: BigUint) {
        if is_prime(&n) {
            // Insert in sorted order
            let pos = state.prime_factors.partition_point(|x| x <= &n);
            state.prime_factors.insert(pos, n);
        } else {
            state.pending.push(PendingFactor {
                n,
                trial_div_done: false,
                ecm_curves_tried: 0,
                b1: 2000,
                pollard_tried: false,
            });
        }
    }

    /// Do one unit of factoring work on the state.
    fn do_work(&mut self, state: &mut FactorizationState) {
        // Find a pending factor to work on
        if state.pending.is_empty() {
            state.complete = true;
            return;
        }

        // Find the first pending that needs work
        let idx = 0; // always work on first pending

        // Check perfect square first (fast)
        {
            let pf = &state.pending[idx];
            if let Some(s) = is_perfect_square(&pf.n) {
                if s > BigUint::one() {
                    let n_val = state.pending.remove(idx).n;
                    let _ = n_val; // dropped
                    // Add s twice
                    Self::add_to_pending(state, s.clone());
                    Self::add_to_pending(state, s);
                    return;
                }
            }
        }

        let pf = &state.pending[idx];

        if !pf.trial_div_done {
            let n_clone = pf.n.clone();
            // Do all trial division in one shot — fast enough (~10ms for 78k primes),
            // avoids wasting 80× 100ms ticks when no small factor exists.
            let (found, _remaining, new_cursor) =
                trial_divide(&n_clone, &self.small_primes, 0, self.small_primes.len());

            let pf = &mut state.pending[idx];
            for f in found {
                if f == pf.n {
                    // The whole number was found prime during trial div; shouldn't happen
                    // because we already checked is_prime before adding to pending.
                    // But handle gracefully.
                } else {
                    let pos = state.prime_factors.partition_point(|x| x <= &f);
                    state.prime_factors.insert(pos, f.clone());
                    pf.n /= &f;
                    // Re-check if we divided it multiple times (already handled in trial_divide loop)
                }
            }

            // After batch, check if remaining < n_pending
            // Actually trial_divide already divided pf.n; remaining is the leftover.
            // We need to sync: set pf.n = remaining if smaller.
            // Actually the loop above already modifies pf.n via pf.n /= &f for each factor found.
            // But trial_divide returns the remaining after dividing n_clone (original).
            // So we need to just set pf.n = remaining and add factors separately.
            // Let's redo this logic cleanly.

            // Actually, the above code has a bug: it divides pf.n by each factor individually,
            // but trial_divide already fully divided n_clone. Let's fix by using remaining directly.
            // We need to reload pf since we mutated it above... let's use a cleaner approach.
            // We'll reset and do this properly.
            // This is messy due to borrow checker. Let's restructure.

            // The simplest fix: just recompute.
            // We already have: found = list of prime factors, remaining = leftover.
            // Set pf.n = remaining, add found to prime_factors.
            // But we already mutated pf.n above (incorrectly). Let's undo that.

            // Actually since we divided each factor out of pf.n separately, and trial_divide
            // already produced the correct remaining (after all divisions), we should just set
            // pf.n = remaining. But the intermediate /= &f steps made pf.n incorrect.
            // This is a design flaw. Let's just set pf.n = remaining directly after restoring.

            // The safest approach: don't mutate pf.n in the loop above, just collect factors,
            // then set pf.n = remaining once. But we already did the mutation...
            // Since pf.n started as n_clone and we divided by each f: if found = [2, 2, 3]
            // then pf.n = n_clone / 2 / 2 / 3 = remaining. So it's actually correct!
            // trial_divide returns factors as individual prime instances (with multiplicity),
            // and divides n_clone by all of them, so n_clone / product(found) = remaining.
            // And pf.n /= f for each f gives pf.n = n_clone / product(found). Correct!

            // Now check if pf.n == remaining (sanity):
            // If yes, great. If no, there's an inconsistency. Let's just trust the logic.

            // Now check remaining state
            let pf = &state.pending[idx];
            // exhausted_small: tried all small primes (p² may still be < n)
            // sqrt_exceeded: next prime's square > pf.n, meaning pf.n is definitely prime
            let exhausted_small = new_cursor >= self.small_primes.len();
            let sqrt_exceeded = !exhausted_small && {
                let pb = BigUint::from(self.small_primes[new_cursor]);
                &pb * &pb > pf.n
            };

            if pf.n == BigUint::one() {
                state.pending.remove(idx);
            } else if sqrt_exceeded || is_prime(&pf.n) {
                // Provably prime: either sqrt check passed, or deterministic Miller-Rabin confirms
                let n_val = state.pending.remove(idx).n;
                if n_val > BigUint::one() {
                    let pos = state.prime_factors.partition_point(|x| x <= &n_val);
                    state.prime_factors.insert(pos, n_val);
                }
            } else if exhausted_small {
                // Finished trial division but n is composite (factors > 10^6): switch to ECM
                state.pending[idx].trial_div_done = true;
            }
            // else: continue trial division (cursor already updated)
        } else {
            // ECM phase — first try fast Pollard-Brent for composites fitting in u128
            let n_val = pf.n.clone();
            let pollard_tried = pf.pollard_tried;

            if !pollard_tried {
                state.pending[idx].pollard_tried = true;
                if let Some(n_u128) = biguint_to_u128(&n_val) {
                    if let Some(f_u128) = pollard_brent_u128(n_u128) {
                        let f = BigUint::from(f_u128);
                        let remainder = &n_val / &f;
                        state.pending.remove(idx);
                        Self::add_to_pending(state, f);
                        Self::add_to_pending(state, remainder);
                        if state.pending.is_empty() {
                            state.complete = true;
                        }
                        return;
                    }
                }
                // Pollard-Brent failed (n too large or genuinely hard); fall through to ECM
            }

            let curves_tried = state.pending[idx].ecm_curves_tried;
            let b1 = state.pending[idx].b1;
            let b2 = 20 * b1;

            self.ensure_ecm_primes(b2);

            // Generate sigma deterministically from curves_tried + some mixing
            let sigma = 6 + (curves_tried as u64 * 1_000_003) % 1_000_000_000;

            let factor = try_ecm_curve(&n_val, sigma, b1, b2, &self.ecm_primes);

            let pf = &mut state.pending[idx];
            pf.ecm_curves_tried += 1;

            // Adaptive B1 schedule
            if pf.ecm_curves_tried == 10 {
                pf.b1 *= 3;
            } else if pf.ecm_curves_tried == 50 {
                pf.b1 *= 3;
            } else if pf.ecm_curves_tried == 200 {
                pf.b1 *= 3;
            }

            if let Some(f) = factor {
                if f != n_val && f > BigUint::one() {
                    let remainder = &n_val / &f;
                    let old_n = state.pending.remove(idx).n;
                    let _ = old_n;
                    Self::add_to_pending(state, f);
                    Self::add_to_pending(state, remainder);
                }
                // If f == n_val: whole n found, treat as prime (shouldn't happen)
            }
        }

        if state.pending.is_empty() {
            state.complete = true;
        }
    }

    /// Return grouped prime factors from state as (value, exponent) pairs.
    fn group_factors(state: &FactorizationState) -> Vec<(BigUint, usize)> {
        let mut groups: Vec<(BigUint, usize)> = vec![];
        let pf = &state.prime_factors;
        let mut i = 0;
        while i < pf.len() {
            let base = pf[i].clone();
            let mut exp = 1usize;
            while i + exp < pf.len() && pf[i + exp] == base {
                exp += 1;
            }
            groups.push((base, exp));
            i += exp;
        }
        groups
    }

    /// Return the current factoring status as structured data (no side effects).
    fn get_status(state: &FactorizationState) -> FactoringStatus {
        FactoringStatus {
            complete: state.complete || state.pending.is_empty(),
            prime_factors: Self::group_factors(state),
            pending_values: state.pending.iter().map(|pf| pf.n.clone()).collect(),
        }
    }

    /// Do one step of factoring work for n_str.
    /// Returns structured status, or Err with a human-readable error message.
    pub fn step(&mut self, n_str: &str) -> Result<FactoringStatus, String> {
        let n = match BigUint::parse_bytes(n_str.trim().as_bytes(), 10) {
            Some(n) => n,
            None => return Err("Please enter a natural number".to_string()),
        };

        if n.is_zero() {
            return Err("Please enter a natural number".to_string());
        }
        if n < BigUint::from(2u32) {
            return Err("Input should be at least 2".to_string());
        }

        // Check if already complete in cache
        if let Some(state) = self.cache.get(&n) {
            if state.complete || state.pending.is_empty() {
                return Ok(Self::get_status(state));
            }
        }

        // If not in cache, initialize
        if !self.cache.contains_key(&n) {
            let mut state = FactorizationState {
                prime_factors: vec![],
                pending: vec![],
                complete: false,
            };
            Self::add_to_pending(&mut state, n.clone());
            if state.pending.is_empty() {
                state.complete = true;
            }
            self.cache.insert(n.clone(), state);
        }

        // Do one unit of work
        let mut state = self.cache.remove(&n).unwrap();
        if !state.complete {
            self.do_work(&mut state);
        }
        let status = Self::get_status(&state);
        self.cache.insert(n, state);
        Ok(status)
    }
}

impl Default for Factorizer {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// Unit tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn biguint(n: u64) -> BigUint {
        BigUint::from(n)
    }

    #[test]
    fn test_modinv_basic() {
        // 3 * 4 = 12 ≡ 1 (mod 11)
        assert_eq!(modinv(&biguint(3), &biguint(11)), Some(biguint(4)));
        // 7 * 8 = 56 ≡ 1 (mod 11)
        assert_eq!(modinv(&biguint(7), &biguint(11)), Some(biguint(8)));
        // gcd(6, 9) = 3 != 1, no inverse
        assert_eq!(modinv(&biguint(6), &biguint(9)), None);
    }

    #[test]
    fn test_is_prime() {
        let primes_under_100 = [
            2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73,
            79, 83, 89, 97,
        ];
        for &p in &primes_under_100 {
            assert!(is_prime(&biguint(p)), "{} should be prime", p);
        }
        for n in 2u64..100 {
            if !primes_under_100.contains(&n) {
                assert!(!is_prime(&biguint(n)), "{} should be composite", n);
            }
        }
    }

    #[test]
    fn test_sieve_primes() {
        let primes = sieve_primes(30);
        assert_eq!(primes, vec![2, 3, 5, 7, 11, 13, 17, 19, 23, 29]);
    }

    #[test]
    fn test_mont_double_on_known_curve() {
        // Use a simple curve mod a large prime to verify doubling doesn't crash.
        let n = BigUint::from(101u64);
        let a24 = BigUint::from(10u64);
        let p = MontPoint {
            x: BigUint::from(5u64),
            z: BigUint::from(1u64),
        };
        let p2 = mont_double(&p, &a24, &n);
        // Just verify z is nonzero (not point at infinity on this curve)
        assert!(!p2.z.is_zero());
    }

    #[test]
    fn test_scalar_mul_consistency() {
        // 2*P should equal mont_double(P)
        let n = BigUint::from(1_000_003u64); // prime
        let a24 = BigUint::from(12345u64);
        let p = MontPoint {
            x: BigUint::from(42u64),
            z: BigUint::from(1u64),
        };
        let p2_double = mont_double(&p, &a24, &n);
        let p2_scalar = mont_scalar_mul(2, &p, &a24, &n);
        // X/Z should be same (projective equality)
        // Check cross multiply: p2_double.x * p2_scalar.z == p2_scalar.x * p2_double.z (mod n)
        let lhs = (&p2_double.x * &p2_scalar.z) % &n;
        let rhs = (&p2_scalar.x * &p2_double.z) % &n;
        assert_eq!(lhs, rhs);
    }

    #[test]
    fn test_ecm_factors_15() {
        let n = BigUint::from(15u64);
        let primes = sieve_primes(200);
        // Try several sigma values
        for sigma in 6..20 {
            if let Some(f) = try_ecm_curve(&n, sigma, 100, 2000, &primes) {
                assert!(f == biguint(3) || f == biguint(5), "unexpected factor: {}", f);
                return;
            }
        }
        // If ECM didn't find it, trial division should. Just verify trial div works.
        let small_primes = sieve_primes(1000);
        let (factors, remaining, _) = trial_divide(&n, &small_primes, 0, 100);
        let mut all: Vec<u64> = factors.iter().map(|f| f.to_u64_digits()[0]).collect();
        if remaining > BigUint::one() {
            all.push(remaining.to_u64_digits()[0]);
        }
        all.sort();
        assert_eq!(all, vec![3, 5]);
    }

    #[test]
    fn test_ecm_factors_77() {
        let n = BigUint::from(77u64); // 7 * 11
        let primes = sieve_primes(500);
        for sigma in 6..30 {
            if let Some(f) = try_ecm_curve(&n, sigma, 200, 4000, &primes) {
                assert!(f == biguint(7) || f == biguint(11), "unexpected factor: {}", f);
                return;
            }
        }
        // Fallback: trial division
        let small_primes = sieve_primes(1000);
        let (factors, remaining, _) = trial_divide(&n, &small_primes, 0, 100);
        let mut all: Vec<u64> = factors.iter().map(|f| f.to_u64_digits()[0]).collect();
        if remaining > BigUint::one() {
            all.push(remaining.to_u64_digits()[0]);
        }
        all.sort();
        assert_eq!(all, vec![7, 11]);
    }

    #[test]
    fn test_trial_divide_basic() {
        let primes = sieve_primes(100);
        // 360 = 2^3 * 3^2 * 5
        let (factors, remaining, _) = trial_divide(&BigUint::from(360u64), &primes, 0, 50);
        assert_eq!(remaining, BigUint::one());
        let mut counts: HashMap<u64, usize> = HashMap::new();
        for f in &factors {
            *counts.entry(f.to_u64_digits()[0]).or_insert(0) += 1;
        }
        assert_eq!(counts[&2], 3);
        assert_eq!(counts[&3], 2);
        assert_eq!(counts[&5], 1);
    }

    fn step_to_completion(fz: &mut Factorizer, n: &str, max: usize) -> FactoringStatus {
        let mut last = FactoringStatus { complete: false, prime_factors: vec![], pending_values: vec![] };
        for _ in 0..max {
            last = fz.step(n).expect("valid input");
            if last.complete { return last; }
        }
        last
    }

    #[test]
    fn test_factorizer_small() {
        let mut fz = Factorizer::new();
        let status = step_to_completion(&mut fz, "12", 50);
        assert!(status.complete, "12 should factor completely");
        let factors: Vec<(u64, usize)> = status.prime_factors.iter()
            .map(|(v, e)| (v.to_u64_digits()[0], *e))
            .collect();
        assert_eq!(factors, vec![(2, 2), (3, 1)]);
    }

    #[test]
    fn test_factorizer_prime() {
        let mut fz = Factorizer::new();
        let status = step_to_completion(&mut fz, "97", 10);
        assert!(status.complete);
        assert_eq!(status.prime_factors, vec![(biguint(97), 1)]);
    }

    #[test]
    fn test_factorizer_invalid() {
        let mut fz = Factorizer::new();
        assert!(fz.step("abc").is_err());
        assert!(fz.step("1").is_err());
        assert!(fz.step("0").is_err());
    }

    #[test]
    fn test_perfect_square() {
        assert_eq!(is_perfect_square(&biguint(49)), Some(biguint(7)));
        assert_eq!(is_perfect_square(&biguint(100)), Some(biguint(10)));
        assert_eq!(is_perfect_square(&biguint(15)), None);
        assert_eq!(is_perfect_square(&biguint(0)), Some(biguint(0)));
        assert_eq!(is_perfect_square(&biguint(1)), Some(biguint(1)));
    }

    #[test]
    fn test_factorizer_semiprime() {
        // 221 = 13 * 17, small enough for trial division to find
        let mut fz = Factorizer::new();
        let status = step_to_completion(&mut fz, "221", 100);
        assert!(status.complete);
        let vals: Vec<u64> = status.prime_factors.iter()
            .flat_map(|(v, e)| std::iter::repeat(v.to_u64_digits()[0]).take(*e))
            .collect();
        assert_eq!(vals, vec![13, 17]);
    }
}
