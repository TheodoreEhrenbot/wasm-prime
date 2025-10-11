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

pub struct PrimeChecker {
    results: HashMap<BigUint, PrimeStatus>,
}

impl PrimeChecker {
    pub fn new() -> Self {
        PrimeChecker {
            results: HashMap::new(),
        }
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
        let checker = checker.clone();

        Callback::from(move |e: InputEvent| {
            let target: HtmlInputElement = e.target_unchecked_into::<HtmlInputElement>();
            let value = target.value();
            input.set(value.clone());

            // Check immediately on input
            let mut checker_mut = checker.borrow_mut();
            let new_result = checker_mut.check(&value);
            result.set(new_result);
        })
    };

    html! {
        <div>
            <input
                type="text"
                placeholder="Enter a number"
                value={(*input).clone()}
                oninput={on_input}
            />
            <div>{ (*result).clone() }</div>
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
