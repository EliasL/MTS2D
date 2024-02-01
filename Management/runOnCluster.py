from connectToCluster import connectToCluster, uploadProject
import time
import sys

uploadProject()
ssh = connectToCluster()

# Connect to the cluster
ssh = connectToCluster()

# Assuming SSH connection is already established and command is executed
stdin, stdout, stderr = ssh.exec_command('python3 /home/elundheim/simulation/Management/runSimulation.py')

# Stream both stdout and stderr
while True:
    # Handle stdout
    if stdout.channel.recv_ready():
        output = stdout.channel.recv(4096).decode('utf-8')
        print(output, end='')

    # Handle stderr
    if stderr.channel.recv_ready():
        error = stderr.channel.recv(4096).decode('utf-8')
        print(error, end='', file=sys.stderr)

    if stdout.channel.exit_status_ready():
        # Check if there are any remaining data in the buffers
        if stdout.channel.recv_ready() or stderr.channel.recv_ready():
            continue  # Handle any remaining data
        else:
            break  # Exit the loop if no more data to read

    # Prevent tight loop
    time.sleep(0.1)

# No need to read from stdout or stderr outside the loop since all data is read within the loop
# Close the SSH connection
ssh.close()