# Contributing to ComplexRegions.jl

Thanks for your interest in contributing! Bug reports, feature requests,
documentation improvements, and pull requests are all welcome.

## Reporting issues

Please open an issue on the [issue tracker](https://github.com/complexvariables/ComplexRegions.jl/issues).
A good report includes:

- the version of ComplexRegions.jl and Julia you are using (`] status` and `versioninfo()`),
- a short, self-contained example that reproduces the problem, and
- what you expected to happen versus what actually happened.

## Development setup

1. Fork and clone the repository.
2. From the repository root, activate the project and instantiate its dependencies:

   ```julia
   julia --project=. -e 'using Pkg; Pkg.instantiate()'
   ```

3. Make your changes on a feature branch (not `master`).

## Running the tests

The test suite runs against both `Float64` and `BigFloat`, and includes
[Aqua](https://github.com/JuliaTesting/Aqua.jl) quality checks and
[JET](https://github.com/aviatesk/JET.jl) static analysis.

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

or from the Julia REPL package mode:

```julia
] test
```

Please make sure the full suite passes before opening a pull request, and add
tests covering any new behavior or bug fix.

## Coding conventions

- Follow the style of the surrounding code.
- Types are capitalized (e.g. `Segment`, `Circle`); methods that merely produce
  values of existing types are lowercase (e.g. `n_gon`).
- Use `tolerance(T)` for floating-point comparisons rather than hardcoded values,
  so that code works correctly across real types including `BigFloat`.
- See [`CLAUDE.md`](CLAUDE.md) for a tour of the type hierarchy, design
  conventions, and notes on extending the package with new curve or region types.

## Pull requests

- Keep each pull request focused on a single change.
- Include a clear description of what the change does and why.
- Update documentation and docstrings when behavior changes.
- The CI workflow runs the test suite and reports coverage; please check that it
  passes.

## AI-assisted contributions

Contributions developed with the help of large language models are welcome, but
you must understand and review the code you submit; please do not open pull
requests containing unreviewed generated code. If a contribution includes
substantial AI-generated content, note that in the pull request. This mirrors
the [Julia General registry guidance](https://github.com/JuliaRegistries/General).

## License

By contributing, you agree that your contributions will be licensed under the
[MIT License](LICENSE) that covers this project.
