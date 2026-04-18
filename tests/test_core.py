import pytest
from nucleic_tools.core import DNASequence, RNASequence, AminoAcidSequence

class TestDNASequence:
    """
    Тесты для ДНК последовательностей
    """
    
    def test_complement(self):
        """complement() возвращает правильную комплементарную цепь"""
        dna = DNASequence("ATGC")
        assert dna.complement() == "TACG"
    
    def test_transcribe(self):
        """transcribe() правильно конвертирует ДНК в РНК"""
        dna = DNASequence("ATGC")
        rna = dna.transcribe()
        assert str(rna) == "AUGC"
        assert isinstance(rna, RNASequence)
    
    def test_gc_content_rna(self):
        """gc_content() правильно считает GC-состав для РНК"""
        rna = RNASequence("AUGCCGCAU")
        # G+C = G(2) + C(2) = 4 из 9 = 44.44%
        assert abs(rna.gc_content() - 44.44) < 0.01
    
    def test_invalid_dna_sequence_raises_error(self):
        """Ошибка - при создании ДНК с неверными символами"""
        with pytest.raises(ValueError, match="Invalid DNA sequence"):
            DNASequence("ATGCX")  # X - неверный нуклеотид

