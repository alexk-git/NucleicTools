import logging
from pathlib import Path

# Создаём папку для логов, если её нет
Path("logs").mkdir(exist_ok=True)

def setup_logger(name: str = None) -> logging.Logger:
    """Настраивает и возвращает логгер с указанным именем"""
    
    logger = logging.getLogger(name)
    
    # Чтобы не добавлять обработчики несколько раз
    if logger.handlers:
        return logger
    
    logger.setLevel(logging.INFO)
    
    # Формат логов
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    
    # Обработчик для файла
    file_handler = logging.FileHandler('logs/app.log')
    file_handler.setFormatter(formatter)
    
    # Обработчик для консоли
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)
    
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)
    
    return logger
