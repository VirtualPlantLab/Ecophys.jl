# Ecophys release notes

We started keeping track of changes in the `NEWS.md` file after version 0.1.1.

# Ecophys 0.1.1 (2026-06-17)

Bug fixes and improve test suite (output will change gs when using both C3 adn C4 models, no changes to API).

- Switch to solving for Ci and then calculate gs inside the analytical functions of C3 and C4 photosynthesis instead of solving for gs directly (this was causing wrong values of gs).

- Add a bunch of tests to better verify several properties of response curves of CO2 assimilation and stomatal conductance.