import os
import sys

def merge_filtered_vcfs(input_directory, output_directory, target_sample=None):
    # List all files in the input directory that end with '_filtered.vcf'
    all_files = [file for file in os.listdir(input_directory) if file.endswith('_filtered.vcf')]

    input_files = []
    
    # Logic to handle specific sample target
    if target_sample:
        # Construct expected filename
        target_file = f"{target_sample}_filtered.vcf"
        if target_file in all_files:
            input_files = [target_file]
            print(f"Processing specific sample: {target_file}")
        else:
            print(f"Warning: Target file {target_file} not found in {input_directory}")
            return # Exit if target not found
    else:
        # Process all files (Batch mode)
        input_files = all_files
        print(f"Processing all {len(input_files)} files in directory.")

    for input_file in input_files:
        input_path = os.path.join(input_directory, input_file)
        output_file = input_file.replace('_filtered.vcf', '_filtered_sum.vcf')
        output_path = os.path.join(output_directory, output_file)

        if os.path.isfile(input_path):
            try:
                with open(input_path, 'rt') as f, open(output_path, 'w') as fsum:
                    datack = False
                    olditem = []
                    pos = -1

                    for rows in f:
                        newitem = rows.split('\t')

                        if not datack:
                            if rows[:6] == "#CHROM":
                                datack = True
                            fsum.write(rows)
                        elif datack and (len(newitem) >= 5): # Safety check for line length
                            # Logic to merge sequential SNPs (pos, pos+1)
                            if (pos + 1 == int(newitem[1])) and (len(newitem[3]) == 1) and (len(newitem[4]) == 1):
                                pos = pos + 1
                                olditem[3] = olditem[3] + newitem[3] # Merge REF
                                olditem[4] = olditem[4] + newitem[4] # Merge ALT
                            elif len(olditem) > 0:
                                stem = "\t".join(olditem)
                                fsum.write(stem)
                                pos = int(newitem[1])
                                olditem = newitem[:]
                            else:
                                pos = int(newitem[1])
                                olditem = newitem[:]
                        
                    # Write the last item
                    if len(olditem) > 0:
                        stem = "\t".join(olditem)
                        fsum.write(stem)
                        
            except Exception as e:
                print(f"Error processing {input_file}: {str(e)}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python merge_vcf.py <input_dir> <output_dir> [sample_name]")
        sys.exit(1)

    input_directory = sys.argv[1]
    output_directory = sys.argv[2]
    
    # Optional 3rd argument for sample name
    target_sample = None
    if len(sys.argv) > 3:
        target_sample = sys.argv[3]

    merge_filtered_vcfs(input_directory, output_directory, target_sample)