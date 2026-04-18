import pytest
from nucleic_tools.dna_rna_tools import run_dna_rna_tools

class TestRunDnaRnaTools:
    """Тесты для функции run_dna_rna_tools"""
    
    def test_single_sequence_transcribe(self):
        """run_dna_rna_tools с одной последовательностью"""
        result = run_dna_rna_tools('ATG', 'transcribe')
        assert result == 'AUG'
    
    def test_multiple_sequences_reverse(self):
        """run_dna_rna_tools с несколькими последовательностями"""
        result = run_dna_rna_tools('ATG', 'CGT', 'reverse')
        assert result == ['GTA', 'TGC']
    
    def test_complement_with_mixed_case(self):
        """проверка на комплиментарность с учётом регистра"""
        result = run_dna_rna_tools('AtGcCg', 'complement')
        assert result.upper() == 'TACGGC'

