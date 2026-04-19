import pytest
import tempfile
import os
from bio_files_processor import convert_multiline_fasta_to_oneline


class TestConvertMultilineFasta:
    """Тесты для функции convert_multiline_fasta_to_oneline"""
    
    def test_convert_multiline_to_oneline(self):
        """Тест на чтение/запись: многострочный FASTA превращается в однострочный"""
        # Создаём временный входной файл
        input_filename = "test_multiline.fasta"
        output_filename = f"output_{input_filename}"
        
        try:
            # Записываем тестовый FASTA с многострочной последовательностью
            with open(input_filename, "w") as f:
                f.write(">seq1\n")
                f.write("ATGC\n")
                f.write("CGTA\n")
                f.write(">seq2\n")
                f.write("AAAA\n")
                f.write("CCCC\n")
            
            # Вызываем тестируемую функцию
            convert_multiline_fasta_to_oneline(input_filename)
            
            # Проверяем, что выходной файл создался
            assert os.path.exists(output_filename)
            
            # Читаем выходной файл и проверяем результат
            with open(output_filename, "r") as f:
                content = f.read()
            
            # Ожидаемый результат: последовательности на одной строке
            expected = ">seq1\nATGCCGTA\n>seq2\nAAAACCCC\n"
            assert content == expected
            
        finally:
            # Удаляем временные файлы
            if os.path.exists(input_filename):
                os.remove(input_filename)
            if os.path.exists(output_filename):
                os.remove(output_filename)
    
    def test_convert_single_line_fasta(self):
        """Тест: если последовательность уже в одну строку, ничего не ломается"""
        input_filename = "test_oneline.fasta"
        output_filename = f"output_{input_filename}"
        
        try:
            with open(input_filename, "w") as f:
                f.write(">seq1\n")
                f.write("ATGCGTA\n")
                f.write(">seq2\n")
                f.write("AAAACCCC\n")
            
            convert_multiline_fasta_to_oneline(input_filename)
            
            with open(output_filename, "r") as f:
                content = f.read()
            
            expected = ">seq1\nATGCGTA\n>seq2\nAAAACCCC\n"
            assert content == expected
            
        finally:
            if os.path.exists(input_filename):
                os.remove(input_filename)
            if os.path.exists(output_filename):
                os.remove(output_filename)
    
    def test_empty_fasta_file(self):
        """Тест на ошибку: пустой файл не вызывает исключений (или вызывает, если так задумано)"""
        input_filename = "test_empty.fasta"
        output_filename = f"output_{input_filename}"
        
        try:
            # Создаём пустой файл
            with open(input_filename, "w") as f:
                pass
            
            # Функция должна отработать без ошибок (или с ошибкой  проверяем)
            # По логике текущей реализации, она создаст пустой выходной файл
            convert_multiline_fasta_to_oneline(input_filename)
            
            assert os.path.exists(output_filename)
            
            with open(output_filename, "r") as f:
                content = f.read()
            
            # Пустой входной файл даёт пустой выходной
            assert content == ""
            
        finally:
            if os.path.exists(input_filename):
                os.remove(input_filename)
            if os.path.exists(output_filename):
                os.remove(output_filename)
    
    def test_file_not_found_error(self):
        """Тест на проверку ошибки: файл не существует"""
        with pytest.raises(FileNotFoundError):
            convert_multiline_fasta_to_oneline("nonexistent_file.fasta")
