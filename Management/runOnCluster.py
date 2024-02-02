from connectToCluster import connectToCluster, uploadProject, Servers
import time
import sys
# Use select for non-blocking I/O
import select

server = Servers.servers[0]
uploadProject(server)
ssh = connectToCluster(server)

# Assuming SSH connection is already established and command is executed
script = "benchmarking.py"
stdin, stdout, stderr = ssh.exec_command(f'python3 -u /home/elundheim/simulation/Management/{script}')

# Stream both stdout and stderr similar to run_command
def stream_ssh_output(stdout, stderr):
    data_read = False

    while True:
        # Use select to wait for output
        ready_to_read, _, _ = select.select([stdout.channel, stderr.channel], [], [], 0.1)
        if not ready_to_read and data_read:
            if stdout.channel.exit_status_ready():
                break  # Exit the loop if command execution is finished

        for channel in ready_to_read:
            line = channel.recv(4096).decode('utf-8')
            if line:
                data_read = True
                # Print to the appropriate stream based on the channel
                output_stream = sys.stdout if channel == stdout.channel else sys.stderr
                print(line, end='', file=output_stream)

        if not data_read:
            print("Waiting for data...", file=sys.stderr)  # Debugging line


stream_ssh_output(stdout, stderr)

# Get exit status
exit_status = stdout.channel.recv_exit_status()

# No need to read from stdout or stderr outside the loop since all data is read within the loop
# Close the SSH connection
ssh.close()