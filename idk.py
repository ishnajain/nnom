from itertools import islice
import re
import csv


stack=[]
def parse_assembly(assembly_code):
    print("in func")
    functions = {}  
    current_function = None


    function_regex = r'^\s*(\S+) <(.*?)>:'
    address_regex = r'^\s*(\S+):'

  
    for line in assembly_code.split('\n'):

        
        match = re.match(function_regex, line)
    
        if match:
            current_function = match.group(2)
            # functions[current_function] = {'start_address': None, 'ret_addresses': []}
            functions[current_function] = {'start_address': None, 'ret_addresses': [], 'execution': {'time':0,'cycle':0,'Ttime': 0, 'Tcycle': 0, 'count': 0}}
            continue

      
        match = re.match(address_regex, line)
        if match and current_function:
            address = match.group(1)
            # print(address)
            if not functions[current_function]['start_address']:
                functions[current_function]['start_address'] = address
            if 'ret' in line or 'mret' in line:
                # print(line)
                # print("found ret or mret")
                dec = re.search(r'\b(ret|mret)\b', line)
                # print(dec)

                # print(address_match)
                address = match.group(1)
                # print(match.group(2),match.group(3))
                # if address_match:
                #     address = address_match.group(1)
                if(dec):
                    functions[current_function]['ret_addresses'].append(address)
            continue
        
    return functions

def write_to_file(functions):
    with open('function_table.csv', 'w', newline='') as csvfile:
        fieldnames = ['Function Name', 'Start Address', 'Return Addresses']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        
        writer.writeheader()
        
        for function, details in functions.items():
            start_address = details['start_address']
            ret_addresses = ', '.join(details['ret_addresses'])
            # Adjust field widths to ensure proper alignment
            writer.writerow({'Function Name': function.ljust(20), 'Start Address': start_address.ljust(10), 'Return Addresses': ret_addresses})




assembly_code = """
<putchar>:
  100084:	0ff57793          	andi	a5,a0,255
  100088:	00020737          	lui	a4,0x20
  10008c:	c31c                	sw	a5,0(a4)
  10008e:	8082                	ret

<puts>:
  100090:	00020737          	lui	a4,0x20
  100094:	00054783          	lbu	a5,0(a0)
  100098:	e399                	bnez	a5,10009e <puts+0xe>
  10009a:	4501                	li	a0,0
  10009c:	8082                	ret
  10009e:	0505                	addi	a0,a0,1
  1000a0:	c31c                	sw	a5,0(a4)
  1000a2:	bfcd                	j	100094 <puts+0x4>
# """
with open('nnom_dump.txt', 'r') as file:
    assembly_code = file.read()
# print(assembly_code)


functions_dict = parse_assembly(assembly_code)
#print(functions_dict)
write_to_file(functions_dict)

def print_function_details(functions_dict):
  """Prints the details of each function in the dictionary on a separate line,
  excluding functions with a 'count' of zero."""
#   sorted_by_time = sorted(updated_func_dict.items(), key=lambda x: x[1]['execution']['Ttime'])
  for function_name, details in functions_dict.items():
    execution_data = details.get('execution', {})  # Get execution data (or empty dict)
    count = execution_data.get('count', 0)

    # Check if count is zero before printing
    if count != 0:
      print(f"Function Name: {function_name}")
      print(f"Start Address: {details.get('start_address')}")
      print(f"Return Addresses: {', '.join(details.get('ret_addresses', []))}")

      print(f"Execution:")
      print(f"\tTime: {execution_data.get('time', 0)}")
      print(f"\tCycle: {execution_data.get('cycle', 0)}")
      print(f"\tTtime: {execution_data.get('Ttime', 0)}")
      print(f"\tTcycle: {execution_data.get('Tcycle', 0)}")
      print(f"\tCount: {count}")
      print()  # Add an 
# Example usage


def read_and_compare_pc(filename,functions_dict):
  with open(filename, 'r') as file:
    next(file)  # Skip the header line
    for line in islice(file, None):  # Iterate through remaining lines
            fields = line.strip().split('\t')
            pc_value = fields[2][2:]  # Assuming PC value is in the third column (index 2)
            time = int(fields[0])  # Assuming time is in the first column (index 0)
            cycle = int(fields[1])  # Assuming cycle is in the second column (index 1)

            for function_name, details in functions_dict.items():
                start_address = details["start_address"]
                ret_addresses = details["ret_addresses"]

                # Retrieve execution dictionary (avoiding unnecessary variable)
                execution = details.get('execution', {})

                # # Ensure execution dictionary exists, initialize if needed
                # if not execution:
                #     execution = {'time': 0, 'cycle': 0, 'Ttime': 0, 'Tcycle': 0, 'count': 0}
                #     details['executions'] = execution  # Update function details

                if pc_value == start_address:
                    # Update execution data for start address
                    execution['time'] = time
                    execution['cycle'] = cycle
                    execution['count'] += 1
                    details['executions'] = execution
                    # print(function_name,execution)
                    stack.append(execution)  
                    #print(function_name, stack)

                elif pc_value in ret_addresses:
                    # Update execution data for return address
                    execution['Ttime'] += (time - execution['time'])
                    execution['Tcycle'] += (cycle - execution['cycle'])
                    # execution['count'] += 1
                    details['executions'] = execution
                    cur_ex=stack.pop()
                    #print(function_name,time-cur_ex['time'])

                    for temp in stack:
                         temp['time'] += time-cur_ex['time']
                         temp['cycle'] += cycle - cur_ex['cycle']
                    
                        #  print(temp['time'])
                

                        


              


    return functions_dict

              # Stop checking functions for this PC value (optional)
    # print_function_details(functions_dict)
   

      
    #   if pc_value == pc_to_compare:
    #     # Found a matching line, process it (e.g., print entire line)
    #     print(f"Found matching PC: {line.strip()}")
    #     break  # Stop iterating if only one match is needed



# Example usage
pc_to_compare = "00106c68"
# func_dict=read_and_compare_pc("../../../../trace_core_00000000.log")
# file_name= "../../../../../../../../media/ishna/Windows-SSD/FYP_windows/tracecore_00000000.log"

updated_func_dict = read_and_compare_pc("../../../../../../../../media/ishna/Windows-SSD/FYP_windows/tracecore_00000000.log" ,functions_dict)
print_function_details(updated_func_dict)
sorted_by_time = sorted(updated_func_dict.items(), key=lambda x: x[1]['execution']['Ttime'],reverse = True)
print(sorted_by_time)
with open('function_data.csv', 'w', newline='') as csvfile:
  writer = csv.writer(csvfile)
  # Write header row
  writer.writerow(['Function Name','Ttime', 'Tcycle', 'Count'])
  for function_name, details in sorted_by_time:
    # Extract details for each row
    # start_address = details.get('start_address')
    # ret_addresses = ','.join(details.get('ret_addresses', []))  # Join list into comma-separated string
    execution_data = details.get('execution', {})
    # execution_time = execution_data.get('time', 0)
    # execution_cycle = execution_data.get('cycle', 0)
    Ttime = execution_data.get('Ttime', 0)
    Tcycle = execution_data.get('Tcycle', 0)
    count = execution_data.get('count', 0)
    # Write each function data as a row
    writer.writerow([function_name , Ttime, Tcycle, count])
# print_function_details(sorted_by_time)

