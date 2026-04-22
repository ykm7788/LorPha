
import argparse
import gzip
import sys
import re
from typing import Dict, Any

def main():
    parser = argparse.ArgumentParser(
        description="Find identical/different entries between two files (table/fasta)",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="""Example:
python find_identical.py xx.fa xx.txt -format 2 -site2 1
python find_identical.py xx.txt xx.txt -format 1
python find_identical.py xx.fa xx.fa -format 3"""
    )
    parser.add_argument("file1", help="First input file (support .gz)")
    parser.add_argument("file2", help="Second input file (support .gz)")
    parser.add_argument("-format", type=int, default=1, 
                        help="Input file type: 1=both table, 2=first fa/second table, 3=both fasta (default:1)")
    parser.add_argument("-site1", type=int, default=1, 
                        help="Matching column/site for file1 (1-based, default:1)")
    parser.add_argument("-site2", type=int, default=1, 
                        help="Matching column/site for file2 (1-based, default:1)")
    parser.add_argument("-type", type=int, default=1, 
                        help="Output type:\n1=identical in file1\n2=identical in file2\n3=identical in both (tab-sep)\n4=different in file1\n5=different in file2 (default:1)")
    parser.add_argument("-help", action="store_true", help="Show help information")
    
    args = parser.parse_args()
    
    if args.help or len(sys.argv) < 3:
        parser.print_help()
        sys.exit(0)
    
    type_ = args.type
    format_ = args.format
    site1 = args.site1 - 1
    site2 = args.site2 - 1
    hash1: Dict[str, str] = {}
    hash2: Dict[str, str] = {}
    rec1: Dict[int, str] = {}
    rec2: Dict[int, str] = {}
    
    def open_file(file_path):
        if file_path.endswith(".gz"):
            return gzip.open(file_path, "rt", encoding="utf-8")
        else:
            return open(file_path, "r", encoding="utf-8")
    
    with open_file(args.file1) as f1:
        if format_ == 1:
            line_num = 0
            for line in f1:
                line_num += 1
                line = line.rstrip("\n")
                if not line:
                    continue
                parts = re.split(r"\s+", line)
                if len(parts) > site1:
                    key = parts[site1]
                    hash1[key] = line
                    rec1[line_num] = key
        elif format_ == 2 or format_ == 3:
            fasta_content = f1.read()
            fasta_records = fasta_content.split(">")[1:]
            line_num = 0
            for record in fasta_records:
                line_num += 1
                record = record.rstrip("\n")
                if not record:
                    continue
                record_parts = re.split(r"\n+", record, maxsplit=1)
                id_line = record_parts[0]
                id_parts = re.split(r"\s+", id_line)
                key = id_parts[0] if id_parts else ""
                hash1[key] = record
                rec1[line_num] = key
    
    with open_file(args.file2) as f2:
        if format_ == 1 or format_ == 2:
            line_num = 0
            for line in f2:
                line_num += 1
                line = line.rstrip("\n")
                if not line:
                    continue
                parts = re.split(r"\s+", line)
                if len(parts) > site2:
                    key = parts[site2]
                    hash2[key] = line
                    rec2[line_num] = key
        elif format_ == 3:
            fasta_content = f2.read()
            fasta_records = fasta_content.split(">")[1:]
            line_num = 0
            for record in fasta_records:
                line_num += 1
                record = record.rstrip("\n")
                if not record:
                    continue
                record_parts = re.split(r"\n+", record, maxsplit=1)
                id_line = record_parts[0]
                id_parts = re.split(r"\s+", id_line)
                key = id_parts[0] if id_parts else ""
                hash2[key] = record
                rec2[line_num] = key
    
    if type_ == 4:
        for key in hash2.keys():
            if key in hash1:
                del hash1[key]
        for line_num in sorted(rec1.keys(), key=int):
            key = rec1[line_num]
            if key in hash1:
                content = hash1[key]
                if format_ == 2 or format_ == 3:
                    print(f">{content}", end="")
                else:
                    print(content)
    
    elif type_ == 5:
        for key in hash1.keys():
            if key in hash2:
                del hash2[key]
        for line_num in sorted(rec2.keys(), key=int):
            key = rec2[line_num]
            if key in hash2:
                content = hash2[key]
                if format_ == 3:
                    print(f">{content}", end="")
                else:
                    print(content)
    
    else:
        if type_ == 1 or type_ == 3:
            sorted_line_nums = sorted(rec1.keys(), key=int)
        else:
            sorted_line_nums = sorted(rec2.keys(), key=int)
        
        for line_num in sorted_line_nums:
            if type_ == 1 or type_ == 3:
                key = rec1[line_num]
                if key in hash2:
                    if type_ == 1:
                        content = hash1[key]
                        if format_ == 2 or format_ == 3:
                            print(f">{content}", end="")
                        else:
                            print(content)
                    elif type_ == 3:
                        content1 = hash1[key]
                        content2 = hash2[key]
                        print(f"{content1}\t{content2}")
            elif type_ == 2:
                key = rec2[line_num]
                if key in hash1:
                    content = hash2[key]
                    if format_ == 3:
                        print(f">{content}", end="")
                    else:
                        print(content)

if __name__ == "__main__":
    main()

