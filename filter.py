import argparse

def read_file_lines(file_path):
    with open(file_path, 'r', encoding='utf-8') as file:
        return file.readlines()

def filter_objects_by_parameter(lines, a_min, a_max):
    header = lines[:6]
    filtered_lines = []
    object_names = []
    
    for line in lines[6:]:
        parts = line.split()
        if len(parts) >= 3:
            try:
                a_value = float(parts[2])
                if a_min <= a_value <= a_max:
                    filtered_lines.append(line)
                    object_names.append(parts[0])
            except ValueError:
                continue  # Skip lines where the third column is not a valid float

    return header + filtered_lines, object_names

def save_to_file(lines, file_path):
    with open(file_path, 'w', encoding='utf-8') as file:
        file.writelines(lines)

def save_names_to_file(names, file_path):
    with open(file_path, 'w', encoding='utf-8') as file:
        file.write(f"# Total object names: {len(names)}\n")
        for name in names:
            file.write(f"{name}\n")

def main(input_file, output_cat_file, output_names_file, a_min, a_max):
    lines = read_file_lines(input_file)
    filtered_lines, object_names = filter_objects_by_parameter(lines, a_min, a_max)
    
    save_to_file(filtered_lines, output_cat_file)
    save_names_to_file(object_names, output_names_file)
    
    print(f"Filtered objects have been saved to {output_cat_file}")
    print(f"Object names have been saved to {output_names_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter objects by parameter 'a' in a .cat file and save results to new files.")
    parser.add_argument("input_file", help="Path to the input .cat file")
    parser.add_argument("output_cat_file", help="Path to the output .cat file")
    parser.add_argument("output_names_file", help="Path to the output names .txt file")
    parser.add_argument("--a_min", type=float, default=41.0, help="Minimum value of parameter 'a'")
    parser.add_argument("--a_max", type=float, default=46.0, help="Maximum value of parameter 'a'")

    args = parser.parse_args()
    main(args.input_file, args.output_cat_file, args.output_names_file, args.a_min, args.a_max)

