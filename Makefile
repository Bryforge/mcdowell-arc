.PHONY: install dev test lint doctor rust-core clean

install:
	python -m pip install -e .

dev:
	python -m pip install --upgrade pip setuptools wheel
	python -m pip install -e ".[dev]"

test:
	python -m pytest -q

lint:
	python -m ruff check src tests

doctor:
	mcdowell-arc doctor

rust-core:
	python -m pip install maturin
	maturin develop --manifest-path crates/mcdowell-arc-core/Cargo.toml --release

clean:
	rm -rf .pytest_cache .ruff_cache build dist src/*.egg-info target profile.out
	find . -type d -name __pycache__ -prune -exec rm -rf {} +
