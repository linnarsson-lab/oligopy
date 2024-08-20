import logging

class SimpleLogger:
    def __init__(self, log_file='app.log', log_level=logging.INFO):
        # Create a logger
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(log_level)

        # Check if the logger has handlers already to avoid adding duplicates
        if not self.logger.handlers:
            # Create a file handler for logging to a file
            file_handler = logging.FileHandler(log_file)
            file_handler.setLevel(log_level)

            # Create a console handler for logging to the console
            console_handler = logging.StreamHandler()
            console_handler.setLevel(log_level)

            # Create a formatter for the log file (with more details)
            file_formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
            file_handler.setFormatter(file_formatter)

            # Create a simple formatter for the console (just the message)
            console_formatter = logging.Formatter('%(message)s')
            console_handler.setFormatter(console_formatter)

            # Add the handlers to the logger
            self.logger.addHandler(file_handler)
            self.logger.addHandler(console_handler)

    def get_logger(self):
        return self.logger