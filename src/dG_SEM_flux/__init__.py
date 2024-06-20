import sys
import logging
import os

# string for logging messages
logging_str = "[%(asctime)s: %(levelname)s: %(module)s: %(message)s]"

# setup a log directory to store log files
log_dir = "logs"
log_filepath = os.path.join(log_dir, "running_logs.log")
os.makedirs(log_dir, exist_ok=True)

logging.basicConfig(
    level= logging.INFO,
    format= logging_str,

    handlers=[
        logging.FileHandler(log_filepath),
        logging.StreamHandler(sys.stdout) # prints log messages in terminal
    ]
)

logger = logging.getLogger("dG_SEM_flux_Logger")
