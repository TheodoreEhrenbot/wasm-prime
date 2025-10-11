use gloo_timers::callback::Interval;
use num_bigint::{BigUint, RandBigInt};
use num_traits::{One, Zero};
use rand::thread_rng;
use std::collections::HashMap;
use web_sys::HtmlInputElement;
use yew::prelude::*;

// Miller-Rabin primality test
fn miller_rabin_test(n: &BigUint, a: &BigUint) -> bool {
    // Write n-1 as 2^r * d
    let n_minus_1 = n - BigUint::one();
    let mut d = n_minus_1.clone();
    let mut r = 0u32;

    while &d % 2u32 == BigUint::zero() {
        d /= 2u32;
        r += 1;
    }

    // Compute a^d mod n
    let mut x = a.modpow(&d, n);

    if x == BigUint::one() || x == n_minus_1 {
        return true;
    }

    for _ in 0..r-1 {
        x = x.modpow(&BigUint::from(2u32), n);
        if x == n_minus_1 {
            return true;
        }
    }

    false
}

pub enum PrimeStatus {
    Composite,
    ProbablyPrime(usize), // Number of Miller-Rabin tests passed
}

#[derive(Clone, Debug)]
pub enum Factor {
    Prime(BigUint),
    Composite(BigUint),
    Unknown(BigUint),
}

pub struct PrimeChecker {
    results: HashMap<BigUint, PrimeStatus>,
}

impl PrimeChecker {
    pub fn new() -> Self {
        PrimeChecker {
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

    pub fn is_composite(&self, n: &BigUint) -> bool {
        if n < &BigUint::from(2u32) {
            return true;
        }
        if n == &BigUint::from(2u32) || n == &BigUint::from(3u32) {
            return false;
        }
        matches!(self.results.get(n), Some(PrimeStatus::Composite))
    }

    pub fn factorize(&mut self, n: &BigUint) -> Vec<Factor> {
        let mut factors = Vec::new();
        let mut remaining = n.clone();
        let mut budget = 1_000_000u64;

        // Handle small cases
        if remaining < BigUint::from(2u32) {
            return vec![Factor::Unknown(remaining)];
        }

        // Factor out 2s
        while &remaining % 2u32 == BigUint::zero() {
            factors.push(Factor::Prime(BigUint::from(2u32)));
            remaining /= 2u32;
        }

        // Try odd divisors up to budget
        let mut divisor = BigUint::from(3u32);
        while &divisor * &divisor <= remaining && budget > 0 {
            while &remaining % &divisor == BigUint::zero() {
                // Check if this divisor is prime
                for _ in 0..10 {
                    self.check(&divisor.to_string());
                }

                if self.is_probably_prime(&divisor) {
                    factors.push(Factor::Prime(divisor.clone()));
                } else if self.is_composite(&divisor) {
                    factors.push(Factor::Composite(divisor.clone()));
                } else {
                    factors.push(Factor::Unknown(divisor.clone()));
                }

                remaining /= &divisor;
            }
            divisor += 2u32;
            budget -= 1;
        }

        // Check what's left
        if remaining > BigUint::one() {
            for _ in 0..10 {
                self.check(&remaining.to_string());
            }

            if self.is_probably_prime(&remaining) {
                factors.push(Factor::Prime(remaining));
            } else if self.is_composite(&remaining) {
                factors.push(Factor::Composite(remaining));
            } else {
                factors.push(Factor::Unknown(remaining));
            }
        }

        factors
    }

    pub fn check(&mut self, n_str: &str) -> String {
        let n = match BigUint::parse_bytes(n_str.as_bytes(), 10) {
            Some(n) => n,
            None => return "Please enter a natural number".to_string(),
        };

        if n == BigUint::zero() || n == BigUint::one() {
            self.results.insert(n, PrimeStatus::Composite);
            return "Input should be at least 2".to_string();
        }

        if n == BigUint::from(2u32) || n == BigUint::from(3u32) {
            return "prime (probability: 1.0)".to_string();
        }

        // Check if even
        if &n % 2u32 == BigUint::zero() {
            self.results.insert(n, PrimeStatus::Composite);
            return "composite".to_string();
        }

        // Check existing status
        let current_tests = match self.results.get(&n) {
            Some(PrimeStatus::Composite) => return "composite".to_string(),
            Some(PrimeStatus::ProbablyPrime(tests)) => *tests,
            None => 0,
        };

        // Run one more Miller-Rabin test
        let mut rng = thread_rng();
        let upper_bound = &n - BigUint::one();
        let lower_bound = BigUint::from(2u32);

        // For n=3, upper_bound=2, so range [2,2) is invalid. Already handled above.
        // For n>=5, this should be fine
        let a = rng.gen_biguint_range(&lower_bound, &upper_bound);

        if !miller_rabin_test(&n, &a) {
            self.results.insert(n, PrimeStatus::Composite);
            return "composite".to_string();
        }

        // Passed another test
        let new_tests = current_tests + 1;
        self.results.insert(n.clone(), PrimeStatus::ProbablyPrime(new_tests));

        // Probability of being composite after k tests is at most (1/4)^k
        let prob_composite = 0.25_f64.powi(new_tests as i32);
        let prob_prime = 1.0 - prob_composite;

        format!("prime (probability: {:.10})", prob_prime)
    }
}

#[function_component(PrimeCheckerApp)]
fn prime_checker_app() -> Html {
    let input = use_state(|| String::new());
    let result = use_state(|| String::new());
    let factorization = use_state(|| Vec::<Factor>::new());
    let checker = use_mut_ref(|| PrimeChecker::new());

    // Set up interval to recheck every 100ms
    {
        let result = result.clone();
        let checker = checker.clone();
        let input_clone = input.clone();

        use_effect(move || {
            let interval_handle = {
                let result = result.clone();
                let checker = checker.clone();
                let input_clone = input_clone.clone();

                Interval::new(100, move || {
                    let input_val = (*input_clone).clone();
                    if !input_val.is_empty() {
                        let mut checker_mut = checker.borrow_mut();
                        let new_result = checker_mut.check(&input_val);
                        result.set(new_result);
                    }
                })
            };

            // Return cleanup function
            move || drop(interval_handle)
        });
    }

    let on_input = {
        let input = input.clone();
        let result = result.clone();
        let factorization = factorization.clone();
        let checker = checker.clone();

        Callback::from(move |e: InputEvent| {
            let target: HtmlInputElement = e.target_unchecked_into::<HtmlInputElement>();
            let value = target.value();
            input.set(value.clone());

            // Check immediately on input
            let mut checker_mut = checker.borrow_mut();
            let new_result = checker_mut.check(&value);
            result.set(new_result);

            // Compute factorization
            if let Some(n) = BigUint::parse_bytes(value.as_bytes(), 10) {
                if n >= BigUint::from(2u32) {
                    let factors = checker_mut.factorize(&n);
                    factorization.set(factors);
                } else {
                    factorization.set(Vec::new());
                }
            } else {
                factorization.set(Vec::new());
            }
        })
    };

    // Determine if we should show prev/next buttons
    let current_n = BigUint::parse_bytes((*input).as_bytes(), 10);

    let show_prev = if let Some(ref n) = current_n {
        n > &BigUint::from(2u32)
    } else {
        false
    };

    let show_next = current_n.is_some();

    let on_prev = {
        let input = input.clone();
        let checker = checker.clone();

        Callback::from(move |_| {
            if let Some(mut n) = BigUint::parse_bytes((*input).as_bytes(), 10) {
                // Search backwards for previous prime
                loop {
                    if n <= BigUint::from(2u32) {
                        break;
                    }
                    n -= BigUint::one();

                    // Check this candidate
                    let mut checker_mut = checker.borrow_mut();
                    for _ in 0..10 {
                        checker_mut.check(&n.to_string());
                    }
                    if checker_mut.is_probably_prime(&n) {
                        input.set(n.to_string());
                        break;
                    }
                }
            }
        })
    };

    let on_next = {
        let input = input.clone();
        let checker = checker.clone();

        Callback::from(move |_| {
            if let Some(mut n) = BigUint::parse_bytes((*input).as_bytes(), 10) {
                // Search forwards for next prime
                loop {
                    n += BigUint::one();

                    // Check this candidate
                    let mut checker_mut = checker.borrow_mut();
                    for _ in 0..10 {
                        checker_mut.check(&n.to_string());
                    }
                    if checker_mut.is_probably_prime(&n) {
                        input.set(n.to_string());
                        break;
                    }
                }
            }
        })
    };

    html! {
        <div>
            <div style="display: flex; gap: 10px; align-items: center;">
                if show_prev {
                    <button onclick={on_prev}>{ "← Previous Prime" }</button>
                }
                <input
                    type="text"
                    placeholder="Enter a number"
                    value={(*input).clone()}
                    oninput={on_input}
                    style="flex: 1;"
                />
                if show_next {
                    <button onclick={on_next}>{ "Next Prime →" }</button>
                }
            </div>
            <div>{ (*result).clone() }</div>
            if !factorization.is_empty() {
                <div style="margin-top: 10px;">
                    <strong>{ "Factorization: " }</strong>
                    {
                        factorization.iter().enumerate().map(|(i, factor)| {
                            let separator = if i > 0 { " × " } else { "" };
                            match factor {
                                Factor::Prime(n) => html! {
                                    <span>
                                        { separator }
                                        <span style="color: green;">{ n.to_string() }</span>
                                    </span>
                                },
                                Factor::Composite(n) => html! {
                                    <span>
                                        { separator }
                                        <span style="color: red; font-weight: bold;">{ n.to_string() }</span>
                                    </span>
                                },
                                Factor::Unknown(n) => html! {
                                    <span>
                                        { separator }
                                        <span style="color: orange;">{ n.to_string() }</span>
                                    </span>
                                },
                            }
                        }).collect::<Html>()
                    }
                </div>
            }
        </div>
    }
}

fn main() {
    yew::Renderer::<PrimeCheckerApp>::new().render();
}

#[cfg(test)]
mod tests {
    use super::*;

    // Naive trial division to check primality (for testing)
    fn is_prime_naive(n: u64) -> bool {
        if n < 2 {
            return false;
        }
        if n == 2 {
            return true;
        }
        if n % 2 == 0 {
            return false;
        }
        let limit = (n as f64).sqrt() as u64 + 1;
        for i in (3..=limit).step_by(2) {
            if n % i == 0 {
                return false;
            }
        }
        true
    }

    #[test]
    fn test_miller_rabin_1_to_10000() {
        let mut checker = PrimeChecker::new();

        for n in 1..=10000 {
            let n_str = n.to_string();
            let expected_prime = is_prime_naive(n);

            // Run multiple checks to get high confidence
            for _ in 0..20 {
                let result = checker.check(&n_str);

                // If we get composite, check it matches expected
                if result == "composite" || result == "Input should be at least 2" {
                    assert!(!expected_prime, "Number {} is prime but Miller-Rabin said composite", n);
                    break;
                }
            }

            // After 20 tests, if it's still saying probably prime, it should actually be prime
            let final_result = checker.check(&n_str);
            if final_result.starts_with("prime") {
                assert!(expected_prime, "Number {} is composite but Miller-Rabin said prime", n);
            }
        }
    }
}
