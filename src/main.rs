use gloo_timers::callback::Interval;
use num_bigint::BigUint;
use num_iter::range;
use std::collections::HashMap;
use web_sys::HtmlInputElement;
use yew::prelude::*;

pub enum PrimeStatus {
    Prime,
    Composite,
    CheckedUntil(BigUint),
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

        if n == BigUint::from(0_u32) || n == BigUint::from(1_u32) {
            self.results.insert(n, PrimeStatus::Composite);
            return "Input should be at least 2".to_string();
        }

        // Check if we've already computed this or have partial results
        let start = match self.results.get(&n) {
            Some(PrimeStatus::Prime) => return "prime".to_string(),
            Some(PrimeStatus::Composite) => return "composite".to_string(),
            Some(PrimeStatus::CheckedUntil(k)) => k + BigUint::from(1_u32),
            None => BigUint::from(2_u32),
        };

        let budget = BigUint::from(1_000_000_u64);

        let end = std::cmp::min(&start + budget, n.clone());

        for i in range(start, end.clone()) {
            if &n % i == BigUint::ZERO {
                self.results.insert(n, PrimeStatus::Composite);
                return "composite".to_string();
            }
        }

        // Update status based on how far we checked
        if end == n {
            // We've checked all possible divisors, it's prime
            self.results.insert(n, PrimeStatus::Prime);
            "prime".to_string()
        } else {
            // We've checked up to 'end', but not all divisors
            self.results.insert(n, PrimeStatus::CheckedUntil(end));
            "computing...".to_string()
        }
    }
}

#[function_component(PrimeCheckerApp)]
fn prime_checker_app() -> Html {
    let input = use_state(|| String::new());
    let result = use_state(|| String::new());
    let checker = use_mut_ref(|| PrimeChecker::new());

    // Set up interval to recheck every 100ms
    {
        let input = input.clone();
        let result = result.clone();
        let checker = checker.clone();

        use_effect_with((), move |_| {
            let interval = Interval::new(100, move || {
                let input_val = (*input).clone();
                if !input_val.is_empty() {
                    let mut checker_mut = checker.borrow_mut();
                    let new_result = checker_mut.check(&input_val);
                    result.set(new_result);
                }
            });

            // Return cleanup function
            move || drop(interval)
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
