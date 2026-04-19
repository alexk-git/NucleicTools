# Tests for NucleicTools

This directory contains pytest tests for the NucleicTools package.

## Test Structure

```
tests/
 conftest.py          - Shared fixtures (temporary files)
 __init__.py
 test_core.py         - Tests for DNASequence, RNASequence, AminoAcidSequence classes
 test_fastq_filter.py - Tests for FASTQ filtering
 test_bio_processor.py - Tests for file I/O operations
```

## Test Coverage

| Test File | Tests | Description |
|-----------|-------|-------------|
| `test_core.py` | 6 | DNA complement, transcription, GC-count, validation; amino acid weight, validation |
| `test_fastq_filter.py` | 1 | Quality-based filtering |
| `test_bio_processor.py` | 1 | File read/write operations |

## Running Tests

```bash
# Install dev dependencies
pip install pytest pytest-cov

# Run all tests
pytest tests/ -v

# Run specific test file
pytest tests/test_core.py -v

# Run with coverage report
pytest tests/ --cov --cov-report=term-missing
```
