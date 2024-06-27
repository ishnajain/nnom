for function_name, details in functions_dict.items():
      start_address = details["start_address"]
      ret_addresses = details["ret_addresses"]
      # Check for match with start address
      if pc_value == start_address:
        print(f"Found matching PC for function: {function_name} (start_address: {start_address})")
      # Check for match with return addresses (if any)
      if pc_value in ret_addresses:
        print(f"Found matching PC as return address for function: {function_name}")
