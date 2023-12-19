import os

def parse_cpp_header(header_path):
    settings = {}
    
    with open(header_path, 'r') as file:
        for line in file:
            if line.startswith('#define'):
                parts = line.split()
                if len(parts) == 3:
                    settings[parts[1]] = parts[2].strip('\"')
    return settings

settings = parse_cpp_header('src/settings.h')