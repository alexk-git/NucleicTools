import pytest
from main import DNASequence, RNASequence, AminoAcidSequence

class TestDNASequence:
    """
    Тесты для ДНК последовательностей
    """
    
    def test_complement(self):
        """complement() возвращает правильную комплементарную цепь"""
        dna = DNASequence("ATGC")
        assert str(dna.complement()) == "TACG"
    
    def test_transcribe(self):
        """transcribe() правильно конвертирует ДНК в РНК"""
        dna = DNASequence("ATGC")
        rna = dna.transcribe()
        assert str(rna) == "AUGC"
        assert isinstance(rna, RNASequence)
    
    def test_gc_count(self):
        """gc_count() возвращает долю GC-состава от 0 до 1"""
        rna = RNASequence("AUGCCGCAU")
        gc_fraction = rna.gc_count()
    
        # Проверяем тип
        assert isinstance(gc_fraction, float)
    
        # Проверяем диапазон
        assert 0 <= gc_fraction <= 1
    
        # Проверяем, что для последовательности без GC возвращается 0
        rna_no_gc = RNASequence("AUAUAUAU")
        assert rna_no_gc.gc_count() == 0.0
    
        # Проверяем, что для последовательности только из GC возвращается 1
        rna_only_gc = RNASequence("GCGCGCGC")
        assert rna_only_gc.gc_count() == 1.0

    def test_invalid_dna_sequence_raises_error(self):
        """Ошибка при создании ДНК с неверными символами"""
        with pytest.raises(ValueError, match="is not DNA"):
            DNASequence("ATGCX")
