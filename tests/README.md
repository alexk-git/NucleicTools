# Tests for NucleicTools

This directory contains pytest tests for the NucleicTools package.

## Test Structure

tests/
 conftest.py # Shared fixtures (temporary files)
 test_core.py # Tests for DNASequence, RNASequence classes (4 tests)
 test_dna_rna_tools.py # Tests for run_dna_rna_tools function (2 tests)
 test_fastq_filter.py # Tests for FASTQ filtering (1 test)
 test_bio_processor.py # Tests for file I/O operations (1 test)

## Test Coverage

| Test File | Tests | Description |
|-----------|-------|-------------|
| `test_core.py` | 4 | Class methods + error handling |
| `test_dna_rna_tools.py` | 2 | Single/multiple sequence operations |
| `test_fastq_filter.py` | 1 | Quality-based filtering |
| `test_bio_processor.py` | 1 | File read/write operations |
| **Total** | **8** | **Includes error + file I/O tests** |

## Running Tests

```bash
# Install dev dependencies
pip install pytest pytest-cov

# Run all tests
pytest tests/ -v

# Run with coverage report
pytest tests/ --cov=nucleic_tools --cov-report=term-missing

# Run specific test file
pytest tests/test_core.py -v
```
