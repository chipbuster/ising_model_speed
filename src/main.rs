use rand::prelude::*;
use rand_xoshiro::rand_core::SeedableRng;
use rand_xoshiro::Xoshiro256Plus;
use std::f64::consts::E;
use std::time::Instant;

/*
This code is largely based off of an implementation I fist saw written in
julia's Ising2D.jl written by genkuroki: https://github.com/genkuroki/Ising2D.jl/blob/master/src/Ising2D.jl
I wanted to see how a fairly naive port to rust would do as most the work
is done is "low level" ways using arrays and for loops. The performance of
this rust code is lackluster, taking twice the time of the julia implementation.
My rewritten stripped down julia implementation is included along side this one
for anyone iterested in testing for the comparison.
*/

fn main() {
    let t = 2. / (1. + 2f64.sqrt()).log(E); //The critical temperature of a 2D Ising model.

    run(true, t);
}
/*
List of constants that will apply for the rest of the program
J is the coupling constant of the interation, it can be a float
in principle, but I will use 1 for the classic ferromagnetic case.
*/
const J: i8 = 1;
const STEPS: usize = 1000;
const SIDE: usize = 500; // Making a default for square arrays
const NROWS: usize = SIDE;
const NCOLUMNS: usize = SIDE;
const LEN: usize = NROWS * NCOLUMNS;

/*
Defining a function to run the simulation inside of main. Order
determines if it is a "hot" or "cold" initial state.
*/
fn run(order: bool, t: f64) {
    /*
    initialize the array we will be using through the rest of the program.
    The rng is needed for determining if the state flips or not and for
    randomizing the initial state when order is passed as true.
    */
    let mut arr = [0i8; LEN];
    let mut rng = Xoshiro256Plus::from_entropy();

    let beta = 1. / t; // beta is a convenience variable for inverse Temp.

    /*
    When order is true the array is randomized between up (1) and down (-1)
    spin states. Otherwise the spins all point up.
    */
    if order == true {
        for x in &mut arr {
            if *&rng.gen::<f64>() > 0.5 {
                *x = 1i8
            } else {
                *x = -1i8
            }
        }
    } else {
        arr = [1i8; LEN];
    }

    let now = Instant::now();
    /*
    Create a static array of probabilites based on the local energy of a given
    site. This can be done because the local energy is one of a discrete number
    of possible outcomes and enables one to check the energy of the site as
    opposed to the whole lattice when determining whether or not to flip the spin.
    */
    let mut probs = [0f64; 9];
    let mut inc: f64 = -4.0;
    for prob in &mut probs {
        *prob = (-2. * beta * inc).exp();
        inc += 1.;
    }

    /*
    This is the major loop of the model. STEPS is the number of times we iterate through
    the lattice trying to flip each spin. j is the iterating over the "columns"
    of the array, and i is iterating over the "rows", though the array is technically
    1D. First the energy of each site is calculated based on the Ising model's
    hamiltonian, h = site * (nn + ss + ee + ww), where h is the energy of the site,
    site is the spin of the value of the site in question, and nn, ss, ee, ww are
    the site above, below, right and left of the current position respectively.
    I'm using wrapping boundary conditions to mitigate edge effects, meaning I
    am treating the right edge of the 2D array as touching the left edge, and
    the top wraps to the bottom. If the energy is lowered by fliping the site's spin we
    do so, otherwise we flip the spin with probability exp(-2 * energy at the site * beta).
    */

    for _ in 0..STEPS {
        for i in 0..NROWS {
            for j in 0..NCOLUMNS {
                let inorth = (i + 1) % NROWS;
                let isouth = if i == 0 { NROWS - 1 } else { i - 1 };
                let jeast = (j + 1) % NCOLUMNS;
                let jwest = if j == 0 { NCOLUMNS - 1 } else { j - 1 };

                let nn = &arr[inorth as usize * NROWS + j];
                let ss = &arr[isouth as usize * NROWS + j];
                let ee = &arr[i * NROWS + (jeast as usize)];
                let ww = &arr[i * NROWS + (jwest as usize)];
                let site = &arr[i * NROWS + j];

                let en = J * site * (nn + ss + ww + ee);
                let pcomp = &rng.gen::<f64>();

                let k = 4 + en;
                if *pcomp < probs[k as usize] {
                    arr[i * NROWS + j] = -1 * arr[i * NROWS + j];
                }
            }
        }
    }

    // Check how long the program took to run.
    let elapsed = Instant::now() - now;
    println!("the whole program took {:#?} seconds to run.", elapsed);
}
