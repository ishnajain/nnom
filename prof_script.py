import re
import csv


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
            functions[current_function] = {'start_address': None, 'ret_addresses': []}
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

    
        # match = re.match(address_regex, line)
        # if match and current_function:
        #     address = match.group(1)
        #     print(address)
        #     if 'ret' in line or 'mret' in line:
        #         functions[current_function]['ret_addresses'].append(address)
        # if current_function and ('ret' in line or 'mret' in line):
        #     if address:
        #         functions[current_function]['ret_addresses'].append(address)
      
        # if current_function:
        #     # print(line)
        #     if ' ret ' in line or 'mret' in line:
                
        #         print("found ret or mret")
        #         address_match = re.search(r'(\S+):\s*\b(ret|mret)\b', line)
        #         if address_match:
        #             address = address_match.group(1)
        #             functions[current_function]['ret_addresses'].append(address)
    return functions

# Function to write function details to a file
# def write_to_file(functions):
#     with open('function_table.txt', 'w') as f:
#         f.write("Function Name\tStart Address\tReturn Addresses\n")
#         for function, details in functions.items():
#             start_address = details['start_address']
#             ret_addresses = ', '.join(details['ret_addresses'])
#             f.write(f"{function}\t{start_address}\t{ret_addresses}\n")
# def write_to_file(functions):
#     with open('function_table.txt', 'w') as f:
#         f.write("Function Name\tStart Address\tReturn Addresses\n")
#         for function, details in functions.items():
#             start_address = details['start_address']
#             ret_addresses = ', '.join(details['ret_addresses'])
#             f.write(f"{function}\t\t\t{start_address}\t\t\t{ret_addresses}\n")
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


functions = parse_assembly(assembly_code)
print(functions)
write_to_file(functions)
