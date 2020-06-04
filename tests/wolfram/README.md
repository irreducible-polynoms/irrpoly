# Check methods' correctness verification using Wolfram
To verify that implemented methods are correct Wolfram could be used.
You need to register Wolfram ID and obtain a free developer license
[here](https://www.wolfram.com/engine/free-license/). Then download
WolframScript from [here](https://www.wolfram.com/wolframscript/).
Alternatively [Wolfram Engine](https://www.wolfram.com/engine/)
could be used.

Wolfram has `IrreduciblePolynomialQ` and `PrimitivePolynomialQ` methods
performing irreducibility and primitivity tests respectively. This methods
were tested using known polynomials for correctness. The result is that
`IrreduciblePolynomialQ` works correctly for all polynomials and
`PrimitivePolynomialQ` works correct for normalized polynomials
(polynomials with leading term equal to `1`). So in checks only normalized
polynomials could be used. Remember that multiplying the polynomial on
number won't change it's irreducibility and primitivity.

This test program generates random polynomials over some defined field,
checks their irreducibility and primitivity and outputs everything.
For performing verification save this output to file `output`
(be free to choose the other name, but remember that file
MUST HAVE NO EXTENSION - `output.txt`  or any other file type 
won't be read by Wolfram). Then open [verify.wsl](verify.wsl)
and replace `output` with the full path to file you created earlier.
Run the [verify.wsl](verify.wsl) with command `wolframscript -file verify.wsl`
providing full path to `verify.wsl` or just by clicking twice on file.
As result you must see the column of `True` - this means that Wolfram check
results equal to this library checks' results. If you see some `False` the
results of checks didn't match.
