# Spectra to Queries

Repository to translate spectra to queries.

## Installation

```bash
git clone https://github.com/spectra-to-knowledge/spectra-to-queries.git
cd spectra-to-queries
```

### Python dependencies

💠 This is just an example for your own scripts using Python (based on the `pyprojet.toml` file).

#### If you do not have Poetry

```bash
if command -v poetry &> /dev/null; then
    echo "Poetry is already installed."
else
    curl -sSL https://install.python-poetry.org | python3 -
fi
```

#### Then

```
poetry install --no-root
```


### R dependencies

```bash
Rscript inst/scripts/install.R
```

## Use 🚀

```bash
Rscript inst/scripts/spectra_to_queries.R
```

Any feedback is welcome!
