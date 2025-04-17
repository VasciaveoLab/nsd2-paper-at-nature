# Import Libraries
from datetime import datetime

import IPython
from IPython.display import clear_output
ip = IPython.get_ipython()

# This function does timestamp printing
def print_timestamp(cell_info):
    clear_output(wait=False)
    print(f"🕒 Execution at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

# (Re)register the clean version
ip.events.register('pre_run_cell', print_timestamp)