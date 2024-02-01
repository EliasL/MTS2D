import os
import subprocess
from connectToCluster import Servers
from getpass import getpass
import subprocess


def generate_ssh_key(key_path):
    """Generate an SSH key pair if it doesn't already exist."""
    if not os.path.exists(key_path) and not os.path.exists(key_path + ".pub"):
        subprocess.run(["ssh-keygen", "-t", "rsa", "-b", "4096", "-f", key_path, "-N", ""], check=True)
        print(f"SSH key generated at {key_path}")
    else:
        print("SSH key already exists.")

def copy_ssh_key_to_server(server_name, username, key_path, password):
    """Copy the public SSH key to the server's authorized keys using sshpass."""
    command = f"sshpass -p {password} ssh-copy-id -i {key_path}.pub {username}@{server_name}"
    subprocess.run(command, shell=True, check=True)
    print(f"SSH key copied to {server_name}")

def change_password_on_server(server_name, username, old_password, new_password):
    """Attempt to change the user's password on the server."""
    # Construct the command sequence for changing the password
    # This sequence is: old password, new password, confirm new password
    passwd_command = f'echo -e "{old_password}\n{new_password}\n{new_password}" | ssh {username}@{server_name} passwd'
    
    try:
        # Execute the command
        result = subprocess.run(passwd_command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print(f"Password changed on {server_name}")
    except subprocess.CalledProcessError as e:
        print(f"Failed to change password on {server_name}: {e.stderr.decode()}")

def main():
    username = "elundheim"  # Change this to your actual username on the servers
    key_path = "/home/elias/Work/ssh/eliasPmmhClusterKey.pub"  # SSH key path
    password = getpass("Enter your SSH password (will not be echoed): ")  # Securely enter password

    # Generate SSH key pair if it doesn't exist
    generate_ssh_key(key_path)

    # Loop through servers and copy the SSH key
    for server in Servers.servers:
        copy_ssh_key_to_server(server, username, key_path, password)

        # Prompt for the new password
    new_password = getpass("Enter the new password (will not be echoed): ")

    # Change password on each server
    for server in Servers.servers:
        change_password_on_server(server, username, new_password)

if __name__ == "__main__":
    main()
