import signal
import subprocess
import os
import constants


def run_with_timeout(filepath, max_area, max_diameter, timeout):
    path_to_program = os.path.join(constants.PATH_TO_REL_BUILD, "main")
    process = subprocess.Popen(
        (path_to_program, '-g', '1',
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


solutions = dict()

for max_diameter in [2]:
    print("Testing with maximum diameter ", max_diameter, "mm")
    solutions[str(max_diameter)] = dict()
    for i in range(10):
        filepath = os.path.join(constants.PATH_TO_EXPERIMENTS, "real", str(i), "points_0.json")
        max_area = 4
        timeout = 60 * 60 * 2
        output = run_with_timeout(filepath, max_area, max_diameter, timeout)
        if output == "":
            solutions[str(max_diameter)][filepath] = "{}"
            print("[TIMEOUT]", filepath)
        else:
            solutions[str(max_diameter)][filepath] = output
            print(output)

print(solutions)