import signal
import subprocess
import os

PATH_TO_PROGRAM = os.path.join("..", "cmake-build-release", "src", "main")
PATH_TO_DATA = os.path.join("..", "data")
PATH_TO_EXPERIMENTS = os.path.join(PATH_TO_DATA, "samples", "experiments")
PATH_TO_REAL = os.path.join(PATH_TO_EXPERIMENTS, "real")


def run_with_timeout(filepath, max_area, max_diameter, timeout):
    process = subprocess.Popen(
        (PATH_TO_PROGRAM, '-g', '1',
         "-f", filepath, "-a", str(max_area), "-d", str(max_diameter), "-r"),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True,
        text=True,
        preexec_fn=lambda: signal.alarm(timeout))
    output = ""
    try:
        for line in iter(process.stdout.readline, ''):
            output += line.rstrip('\n')
        process.stdout.close()
    except subprocess.TimeoutExpired as te:
        process.kill()
    return output


for max_diameter in [2]:
    print("Testing with maximum diameter ", max_diameter, "mm")
    for i in range(1, 10):
        filepath = os.path.join(PATH_TO_REAL, f"{i}", "points_0.json")
        max_area = 4
        timeout = 60 * 60 * 2
        output = run_with_timeout(filepath, max_area, max_diameter, timeout)
        if output == "":
            print("[TIMEOUT]", filepath)
        else:
            print(output)
