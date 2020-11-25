# ForCE - A library for Forward Correction of Errors

ForCE provides a number type that is more accurate than hardware-supported
floating-point types and provides error estimates. For example, a system that
supports only single-precision IEEE numbers can use ForCE to obtain better
precision results, at the same time as obtaining a roundoff error estimate.

Test it by navigating to the `/tests` folder and typing

    make

## How it works

ForCE uses a variant of CENA, a method for the correction of floating-point
roundoff errors. The original CENA method uses previously-published formulas
to compute the local error introduced by operators (e.g. `+-*/`), and
takes inspiration from the reverse-mode of automatic differentiation to combine
local errors into a global error estimate. This error estimate is often so
precise that it can be subtracted from the result to obtain a better result.

Instead of the reverse-propagation that the original CENA method is built on,
ForCE uses forward-mode AD. This allows us to save memory, reduce complexity,
and offer an operator-overloading-based implementation that can be easily used
in an existing C++ program. ForCE naturally works in parallel OpenMP programs.