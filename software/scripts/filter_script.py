"""
A command-line tool to filter tab-separated files based on a specific column value.

Usage:
  ./filter_script.py -a T|F -b <col_index> -c <value> -d 1|2 -i <input_file> -o <output_file>

Arguments:
  -a: Header presence. 'T' for header exists, 'F' otherwise.
  -b: Column index (1-based) to apply filter on.
  -c: Value to compare against in the specified column.
  -d: Filter rule. 1 = equal to -c, 2 = not equal to -c.
  -i: Input tab-separated file path.
  -o: Output file path.
"""

import argparse
import sys


def main():
    parser = argparse.ArgumentParser(description="Filter tab-separated file by column value.")
    parser.add_argument('-a', required=True, choices=['T', 'F'], help="Header presence: T or F")
    parser.add_argument('-b', type=int, required=True, help="Column index (1-based) to filter on")
    parser.add_argument('-c', required=True, help="Value to match in column -b")
    parser.add_argument('-d', type=int, required=True, choices=[1, 2], help="1: equal, 2: not equal")
    parser.add_argument('-i', required=True, help="Input tab-separated file")
    parser.add_argument('-o', required=True, help="Output file")

    args = parser.parse_args()

    has_header = (args.a == 'T')
    col_index = args.b - 1  # Convert to 0-based index
    match_value = args.c
    equal_match = (args.d == 1)

    try:
        with open(args.i, 'r', encoding='utf-8') as infile, \
             open(args.o, 'w', encoding='utf-8') as outfile:

            first_line = True
            for line in infile:
                line = line.rstrip('\n\r')
                fields = line.split('\t')

                # Check if current line is header
                if has_header and first_line:
                    outfile.write(line + '\n')
                    first_line = False
                    continue

                # Skip lines with insufficient columns
                if col_index >= len(fields):
                    continue

                field_value = fields[col_index]
                match = (field_value == match_value)
                if equal_match == match:
                    outfile.write(line + '\n')

    except FileNotFoundError:
        sys.stderr.write(f"Error: Input file '{args.i}' not found.\n")
        sys.exit(1)
    except Exception as e:
        sys.stderr.write(f"Error: {e}\n")
        sys.exit(1)


if __name__ == '__main__':
    main()

