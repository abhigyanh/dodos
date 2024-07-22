import subprocess
import sys
import os

def main():
    # Path to the script in the 'src' directory
    script_path = os.path.join(os.path.dirname(__file__), 'main.py')
    subprocess.run(['python', script_path] + sys.argv[1:])

if __name__ == '__main__':
    main()
