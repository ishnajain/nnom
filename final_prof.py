import csv




stack=[]

# Function to read the function table CSV file
def read_function_table(filename):
    functions = {}
    with open(filename, 'r') as file:
        reader = csv.reader(file)
        next(reader)  # Skip header
        for row in reader:
            
            function_name, start_address, return_addresses = row
            
            functions[function_name.strip()] = {
                'start_address': start_address.strip(),
                'return_addresses': [addr.strip() for addr in return_addresses.split(',') if addr.strip()]
            }
        # print(functions)
    return functions

# Function to read the log file
def read_log_file(filename):
    log_data = []
    with open(filename, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        next(reader)  # Skip header
        for row in reader:
            log_data.append({
                'Time': int(row[0]),
                'Cycle': int(row[1]),
                'PC': row[2][2:],
                'Insn': row[3],
                'Decoded_instruction': row[4],
                'Register_and_memory_contents': row[5]
            })
            # print(row)
    return log_data

# Function to process log data for each function
def process_log_data(function, log_data):
    start_address = function['start_address']
    # print(start_address)
    return_addresses = function['return_addresses']
    start_counts = 0
    return_counts = 0
    start_times = []
    return_times = []
    time_cycles = []
    # infunc=0

    for entry in log_data:
        # print(entry['PC'])
        if entry['PC'] == start_address:
            start_counts += 1
            # infunc=1
            start_times.append((entry['Time'], entry['Cycle']))

        elif (entry['PC'] in return_addresses):
            return_counts += 1
            # infunc=0
            return_times.append((entry['Time'], entry['Cycle']))

    for start_time, start_cycle in start_times:
        for return_time, return_cycle in return_times:
            if start_time <= return_time:
                time_diff = return_time - start_time
                cycle_diff = return_cycle - start_cycle
                time_cycles.append((time_diff, cycle_diff))
                break

    return start_counts, time_cycles

# Function to write results to a new CSV file
def write_results(filename, functions, log_data):
    with open(filename, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Function Name', 'Start Count', 'Time Difference', 'Cycle Difference'])
        # print(functions.items())
        for function_name, function_data in functions.items():
            # print(function_data)
            # print(function_name)
            if not function_data['return_addresses']:
                continue
            # print(function_name)
            start_count, time_cycles = process_log_data(function_data, log_data)
            # print(start_count)
            for time_diff, cycle_diff in time_cycles:
                writer.writerow([function_name, start_count, time_diff, cycle_diff])

# Main function
def main():
    function_table = read_function_table('function_table.csv')
    filename='tracecore_00000000.log'
    log_file_path = '../../../../../../../../media/ishna/Windows-SSD/FYP_windows/' + filename
    log_data = read_log_file(log_file_path)
    out_filename ='prof_log.csv'
    out_file_path='../'+ out_filename
    write_results('prof_log.csv', function_table, log_data)

if __name__ == "__main__":
    main()
