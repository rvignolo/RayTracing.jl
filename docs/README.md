# RayTracing.jl Documentation

Build the documentation from the repository root:

```bash
julia --project=docs docs/make.jl
```

Or from this directory:

```bash
julia --project=. make.jl
```

The wrapper script is equivalent:

```bash
julia --project=docs docs/build_docs.jl
```

The generated site is written to `docs/build/`.

For local development, instantiate the docs environment first:

```bash
julia --project=docs -e 'using Pkg; Pkg.instantiate()'
```

Then rebuild whenever `docs/src/*.md` or source docstrings change.
